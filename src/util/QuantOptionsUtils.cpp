#include "salmon/internal/util/QuantOptionsUtils.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "salmon/internal/config/SalmonDefaults.hpp"
#include "salmon/internal/util/IOUtils.hpp"
#include "salmon/internal/util/SalmonMath.hpp"
#include "salmon/internal/util/SalmonUtils.hpp"

namespace salmon::utils {

namespace {

bool createDirectoryVerbose_(boost::filesystem::path& dirPath) {
  namespace bfs = boost::filesystem;
  if (bfs::exists(dirPath)) {
    // If it already exists and isn't a directory, then complain
    if (!bfs::is_directory(dirPath)) {
      fmt::print(stderr,
                 "{}ERROR{}: Path [{}] already exists "
                 "and is not a directory.\n"
                 "Please either remove this file or choose another "
                 "auxiliary directory.\n",
                 ioutils::SET_RED, ioutils::RESET_COLOR, dirPath.string());
      return false;
    }
  } else { // If the path doesn't exist, then create it
    if (!bfs::create_directories(dirPath)) { // creation failed for some reason
      fmt::print(stderr,
                 "{}ERROR{}: Could not create the directory [{}]. "
                 "Please check that doing so is valid.",
                 ioutils::SET_RED, ioutils::RESET_COLOR, dirPath.string());
      return false;
    }
  }
  return true;
}

// taken from :
// https://www.boost.org/doc/libs/1_67_0/libs/program_options/example/real.cpp
void conflicting_options(const boost::program_options::variables_map& vm,
                         const char* opt1, const char* opt2) {
  if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) &&
      !vm[opt2].defaulted()) {
    throw std::logic_error(std::string("Conflicting options '") + opt1 +
                           "' and '" + opt2 + "'.");
  }
}

// taken from :
// https://www.boost.org/doc/libs/1_67_0/libs/program_options/example/real.cpp
void option_dependency(const boost::program_options::variables_map& vm,
                       const char* for_what, const char* required_option) {
  if (vm.count(for_what) && !vm[for_what].defaulted()) {
    if (vm.count(required_option) == 0 || vm[required_option].defaulted()) {
      throw std::logic_error(std::string("Option '") + for_what +
                             "' requires option '" + required_option + "'.");
    }
  }
}

} // namespace

bool validateOptionsAlignment_(SalmonOpts& sopt,
                               boost::program_options::variables_map& vm) {
  (void)vm;
  if (!sopt.sampleOutput and sopt.sampleUnaligned) {
    sopt.jointLog->warn(
        "You passed in the (-u/--sampleUnaligned) flag, but did not request a "
        "sampled "
        "output file (-s/--sampleOut).  This flag will be ignored!");
  }

  if (sopt.useErrorModel and sopt.rangeFactorizationBins < 4) {
    uint32_t nbins{4};
    sopt.jointLog->info(
        "Usage of --useErrorModel implies use of range factorization. "
        "rangeFactorization bins is being set to {}",
        nbins);
    sopt.rangeFactorizationBins = nbins;
    sopt.useRangeFactorization = true;
  }
  return true;
}

bool validateOptionsMapping_(SalmonOpts& sopt,
                             boost::program_options::variables_map& vm) {
  auto numUnpaired = sopt.unmatedReadFiles.size();
  auto numLeft = sopt.mate1ReadFiles.size();
  auto numRight = sopt.mate2ReadFiles.size();

  if (sopt.mimicBT2 or sopt.mimicStrictBT2 or sopt.hardFilter) {
    sopt.jointLog->info("The --mimicBT2, --mimicStrictBT2 and --hardFilter flags imply mapping validation (--validateMappings). "
                        "Enabling mapping validation.");
    sopt.validateMappings = true;
  }

  if (!sopt.validateMappings) {
    sopt.jointLog->warn("\n\n"
                        "NOTE: It appears you are running salmon without the `--validateMappings` option.\n"
                        "Mapping validation can generally improve both the sensitivity and specificity of mapping,\n"
                        "with only a moderate increase in use of computational resources. \n"
                        "Mapping validation is planned to become a default option (i.e. turned on by default) in\n"
                        "the next release of salmon.\n"
                        "Unless there is a specific reason to do this (e.g. testing on clean simulated data),\n"
                        "`--validateMappings` is generally recommended.\n");
  }

  bool is_pe_library = (numLeft + numRight > 0);
  bool is_se_library = (numUnpaired > 0);

  if (is_pe_library and is_se_library) {
      sopt.jointLog->warn("You seem to have passed in both un-paired reads and paired-end reads. "
                          "It is not currently possible to quantify hybrid library types in salmon.");
  }

  if (is_pe_library) {
    if (numLeft != numRight) {
      sopt.jointLog->error("You passed paired-end files to salmon, but you passed {} files to --mates1 "
                           "and {} files to --mates2.  You must pass the same number of files to both flags",
                           numLeft, numRight);
      return false;
    }
  }

  auto checkScoreValue = [&sopt](int16_t score, std::string sname) -> bool {
                           using score_t = int8_t;
                           auto minval = static_cast<int16_t>(std::numeric_limits<score_t>::min());
                           auto maxval = static_cast<int16_t>(std::numeric_limits<score_t>::max());
                           if (score <  minval or score > maxval) {
                             sopt.jointLog->error("You set the {} as {}, but it must be in "
                                                  "the range [{}, {}].", sname, score, minval, maxval);
                             return false;
                           }
                           return true;
                         };

  if(!checkScoreValue(sopt.matchScore, "match score")) { return false; }
  if(!checkScoreValue(sopt.mismatchPenalty, "mismatch penalty")) { return false; }
  if(!checkScoreValue(sopt.gapOpenPenalty, "gap open penalty")) { return false; }
  if(!checkScoreValue(sopt.gapExtendPenalty, "gap extend penalty")) { return false; }

  if (sopt.mismatchPenalty > 0) {
    sopt.jointLog->warn(
                        "You set the mismatch penalty as {}, but it should be negative.  It is being negated to {}.",
                        sopt.mismatchPenalty, -sopt.mismatchPenalty);
    sopt.mismatchPenalty = -sopt.mismatchPenalty;
  }

  if (sopt.consensusSlack < 0 or sopt.consensusSlack >= 1.0) {
    sopt.jointLog->error("You set consensusSlack as {}, but it must in [0,1).", sopt.consensusSlack);
    return false;
  }

  if (sopt.mismatchSeedSkip < 1) {
    sopt.jointLog->warn("The mismatchSeedSkip was set to {}, but it cannot be < 1.  Setting mismatchSeedSkip to 1");
    sopt.mismatchSeedSkip = 1;
  }

  if (sopt.mismatchSeedSkip > 31) {
    sopt.jointLog->warn("Setting the mismatchSeedSkip too high can hurt the sensitivity of mapping.  Consider "
    "setting a lower mismatchSeedSkip.");
  }

  if (sopt.validateMappings) {
    if (!vm.count("minScoreFraction")) {
      sopt.minScoreFraction = salmon::defaults::minScoreFraction;
      sopt.jointLog->info(
                          "Usage of --validateMappings implies use of minScoreFraction. "
                          "Since not explicitly specified, it is being set to {}", sopt.minScoreFraction
                          );
    }

    if (sopt.hardFilter) {
      if (sopt.rangeFactorizationBins > 0) {
        sopt.jointLog->info("The use of range-factorized equivalence classes does not make sense "
                            "in conjunction with --hardFilter.  Disabling range-factorized equivalence classes. ");
        sopt.rangeFactorizationBins = 0;
        sopt.useRangeFactorization = false;
      }
    } else if (sopt.rangeFactorizationBins < 4) {
      uint32_t nbins{4};
      sopt.jointLog->info(
                          "Usage of --validateMappings, without --hardFilter implies use of range factorization. "
                          "rangeFactorizationBins is being set to {}", nbins
                          );
      sopt.rangeFactorizationBins = nbins;
      sopt.useRangeFactorization = true;
    }

    bool consensusSlackExplicit = !vm["consensusSlack"].defaulted();
    if (!consensusSlackExplicit) {
      sopt.jointLog->info(
                          "Setting consensusSlack to selective-alignment default of {}.",
                          sopt.consensusSlack);
    }

    bool pre_merge_chain_sub_thresh_explicit = !vm["preMergeChainSubThresh"].defaulted();
    bool post_merge_chain_sub_thresh_explicit = !vm["postMergeChainSubThresh"].defaulted();
    bool orphan_chain_sub_thresh_explicit = !vm["orphanChainSubThresh"].defaulted();

    if (is_se_library) {
      if (!pre_merge_chain_sub_thresh_explicit) {
        sopt.pre_merge_chain_sub_thresh = 1.0;
      }

      if (post_merge_chain_sub_thresh_explicit) {
        sopt.jointLog->warn("The postMergeChainSubThresh is not meaningful for single-end "
        "libraries.  Setting this value "
        "to 1.0 and ignoring");
      }
      if (orphan_chain_sub_thresh_explicit) {
        sopt.jointLog->warn("The orphanChainSubThresh is not meaningful for single-end "
        "libraries.  Setting this value "
        "to 1.0 and ignoring");
      }
      sopt.post_merge_chain_sub_thresh = 1.0;
      sopt.orphan_chain_sub_thresh = 1.0;
    }

    if (sopt.pre_merge_chain_sub_thresh < 0 or sopt.pre_merge_chain_sub_thresh > 1.0) {
      sopt.jointLog->error("You set preMergeChainSubThresh as {}, but it must in [0,1].",
        sopt.pre_merge_chain_sub_thresh);
      return false;
    }
    if (sopt.post_merge_chain_sub_thresh < 0 or sopt.post_merge_chain_sub_thresh > 1.0) {
      sopt.jointLog->error("You set postMergeChainSubThresh as {}, but it must in [0,1].",
        sopt.post_merge_chain_sub_thresh);
      return false;
    }
    if (sopt.orphan_chain_sub_thresh < 0 or sopt.orphan_chain_sub_thresh > 1.0) {
      sopt.jointLog->error("You set orphanChainSubThresh as {}, but it must in [0,1].",
        sopt.orphan_chain_sub_thresh);
      return false;
    }

    if (sopt.mimicBT2 and sopt.mimicStrictBT2) {
      sopt.jointLog->error("You passed both the --mimicBT2 and --mimicStrictBT2 parameters.  These are mutually exclusive. "
                            "Please select only one of these flags.");
      return false;
    }

    if (sopt.mimicBT2 or sopt.mimicStrictBT2) {
      sopt.maxReadOccs = 1000;
      sopt.jointLog->info("The --mimicBT2 and --mimicStrictBT2 flags increases maxReadOccs to {}.", sopt.maxReadOccs);

      sopt.consensusSlack = 0.5;
      sopt.jointLog->info("The --mimicBT2 and --mimicStrictBT2 flags increases consensusSlack to {}.", sopt.consensusSlack);

      if (sopt.mimicBT2) {
        sopt.jointLog->info(
                            "Usage of --mimicBT2 overrides other settings for mapping validation. Setting "
                            "Bowtie2-like parameters now.");
        sopt.discardOrphansQuasi = true;
        if (sopt.softclipOverhangs) {
          sopt.jointLog->info("Softclipping of overhangs is not allowed in mimicBT2 mode; setting to false.");
          sopt.softclipOverhangs = false;
        }

        sopt.matchScore = 2;
        sopt.mismatchPenalty = -4;
        sopt.gapOpenPenalty = 5;
        sopt.gapExtendPenalty = 3;
      }

      if (sopt.mimicStrictBT2) {
        sopt.jointLog->info(
                            "Usage of --mimicStrictBT2 overrides other settings for mapping validation. Setting "
                            "strict RSEM+Bowtie2-like parameters now.");
        if (sopt.softclipOverhangs) {
          sopt.jointLog->info("Softclipping of overhangs is not allowed in mimicStrictBT2 mode; setting to false.");
          sopt.softclipOverhangs = false;
        }
        sopt.discardOrphansQuasi = true;
        sopt.minScoreFraction = 0.8;
        sopt.matchScore = 1;
        sopt.mismatchPenalty = 0;
        sopt.gapOpenPenalty = 25;
        sopt.gapExtendPenalty = 25;
      }
    }
  }
  return true;
}

bool createAuxMapLoggers_(SalmonOpts& sopt,
                          boost::program_options::variables_map& vm) {
  (void)vm;
  namespace bfs = boost::filesystem;

  auto jointLog = sopt.jointLog;

  if (sopt.writeUnmappedNames) {
    boost::filesystem::path auxDir = sopt.outputDirectory / sopt.auxDir;
    bool auxSuccess = bfs::exists(auxDir) and bfs::is_directory(auxDir);
    if (!auxSuccess) {
      return false;
    }
    bfs::path unmappedNameFile = auxDir / "unmapped_names.txt";
    std::ofstream* outFile = new std::ofstream(unmappedNameFile.string());
    if (!outFile->is_open()) {
      jointLog->error("Could not create file for unmapped read names [{}]",
                      unmappedNameFile.string());
      delete outFile;
      return false;
    }
    auto outputSink =
        std::make_shared<spdlog::sinks::ostream_sink_mt>(*outFile);

    std::shared_ptr<spdlog::logger> outLog =
        std::make_shared<spdlog::logger>("unmappedLog", outputSink);
    spdlog::register_logger(outLog);
    outLog->set_pattern("%v");
    sopt.unmappedFile.reset(outFile);
    sopt.unmappedLog = outLog;
  }

  if (sopt.writeOrphanLinks) {
    boost::filesystem::path auxDir = sopt.outputDirectory / sopt.auxDir;
    bool auxSuccess = bfs::exists(auxDir) and bfs::is_directory(auxDir);
    if (!auxSuccess) {
      return false;
    }
    bfs::path orphanLinkFile = auxDir / "orphan_links.txt";
    std::ofstream* outFile = new std::ofstream(orphanLinkFile.string());
    if (!outFile->is_open()) {
      jointLog->error("Could not create file for orphan links [{}]",
                      orphanLinkFile.string());
      delete outFile;
      return false;
    }

    auto outputSink =
        std::make_shared<spdlog::sinks::ostream_sink_mt>(*outFile);

    std::shared_ptr<spdlog::logger> outLog =
        std::make_shared<spdlog::logger>("orphanLinkLog", outputSink);
    spdlog::register_logger(outLog);
    outLog->set_pattern("%v");
    sopt.orphanLinkFile.reset(outFile);
    sopt.orphanLinkLog = outLog;
  }

  bool writeQuasimappings = (sopt.qmFileName != "");

  if (writeQuasimappings) {
    std::streambuf* qmBuf{nullptr};
    if (sopt.qmFileName == "-") {
      qmBuf = std::cout.rdbuf();
    } else {
      sopt.qmFileName = boost::filesystem::absolute(sopt.qmFileName).string();
      bfs::path qmDir = boost::filesystem::path(sopt.qmFileName).parent_path();
      bool qmDirSuccess = boost::filesystem::is_directory(qmDir);
      if (!qmDirSuccess) {
        qmDirSuccess = boost::filesystem::create_directories(qmDir);
      }
      if (qmDirSuccess) {
        sopt.qmFile.open(sopt.qmFileName);
        if (!sopt.qmFile.is_open()) {
          jointLog->error(
              "Could not create file for writing selective-alignments [{}]",
              sopt.qmFileName);
          return false;
        }
        qmBuf = sopt.qmFile.rdbuf();
      } else {
        bfs::path qmFileName =
            boost::filesystem::path(sopt.qmFileName).filename();
        jointLog->error("Couldn't create requested directory {} in which "
                        "to place the mapping output {}",
                        qmDir.string(), qmFileName.string());
        return false;
      }
    }
    sopt.qmStream.reset(new std::ostream(qmBuf));

    auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(
        *(sopt.qmStream.get()));
    sopt.qmLog = std::make_shared<spdlog::logger>("qmStream", outputSink);
    sopt.qmLog->set_pattern("%v");
  }
  return true;
}

bool processQuantOptions(SalmonOpts& sopt,
                         boost::program_options::variables_map& vm,
                         int32_t numBiasSamples) {
  using std::string;
  namespace bfs = boost::filesystem;

  sopt.numBiasSamples.store(numBiasSamples);
  sopt.runStartTime = getCurrentTimeAsString();

  bfs::path geneMapPath;
  if (vm.count("geneMap")) {
    geneMapPath = vm["geneMap"].as<std::string>();
    if (!bfs::exists(geneMapPath)) {
      std::cerr << "ERROR: Could not find transcript <=> gene map file "
                << geneMapPath << "\n";
      std::cerr << "Exiting now: please either omit the \'geneMap\' option or "
                   "provide a valid file\n";
      return false;
    }
    sopt.geneMapPath = geneMapPath;
  }

  bfs::path outputDirectory(vm["output"].as<std::string>());
  bool outputDirOK = createDirectoryVerbose_(outputDirectory);
  if (!outputDirOK) {
    return false;
  }
  sopt.outputDirectory = outputDirectory;

  bfs::path logDirectory = outputDirectory / "logs";
  bool logDirOK = createDirectoryVerbose_(logDirectory);
  if (!logDirOK) {
    return false;
  }
  if (!sopt.quiet) {
    std::cerr << "Logs will be written to " << logDirectory.string() << "\n";
  }

  bfs::path paramsDir = outputDirectory / "libParams";
  bool paramDirOK = createDirectoryVerbose_(paramsDir);
  if (!paramDirOK) {
    return false;
  }
  sopt.paramsDirectory = paramsDir;

  bfs::path auxDir = sopt.outputDirectory / sopt.auxDir;
  bool auxDirOK = createDirectoryVerbose_(auxDir);
  if (!auxDirOK) {
    return false;
  }

  if (sopt.meta) {
    sopt.initUniform = true;
    sopt.noRichEqClasses = true;
    sopt.useEM = true;
    sopt.useVBOpt = false;
  }

  size_t max_q_size = 131072;
  bfs::path logPath = logDirectory / "salmon_quant.log";

  if (sopt.quantMode == SalmonQuantMode::MAP) {
    bfs::path indexDirectory(vm["index"].as<string>());
    sopt.indexDirectory = indexDirectory;
    bool writeQuasimappings = (sopt.qmFileName != "");
    if (writeQuasimappings or sopt.writeUnmappedNames or sopt.writeOrphanLinks) {
      max_q_size = 2097152;
    }
  }

  spdlog::init_thread_pool(max_q_size, 1);
  auto fileSink =
      std::make_shared<spdlog::sinks::basic_file_sink_mt>(logPath.string(), true);
  auto consoleSink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  consoleSink->set_color(spdlog::level::warn, consoleSink->magenta);
  auto consoleLog = std::make_shared<spdlog::async_logger>(
      "stderrLog", consoleSink, spdlog::thread_pool(),
      spdlog::async_overflow_policy::block);
  spdlog::register_or_replace(consoleLog);
  auto fileLog = std::make_shared<spdlog::async_logger>(
      "fileLog", fileSink, spdlog::thread_pool(),
      spdlog::async_overflow_policy::block);
  spdlog::register_or_replace(fileLog);
  std::vector<spdlog::sink_ptr> sinks{consoleSink, fileSink};
  auto jointLog = std::make_shared<spdlog::async_logger>(
      "jointLog", std::begin(sinks), std::end(sinks), spdlog::thread_pool(),
      spdlog::async_overflow_policy::block);
  spdlog::register_or_replace(jointLog);

  if (sopt.quiet) {
    jointLog->set_level(spdlog::level::err);
  }

  sopt.jointLog = jointLog;
  sopt.fileLog = fileLog;

  if (sopt.quantMode == SalmonQuantMode::MAP) {
    bool auxLoggersOK = createAuxMapLoggers_(sopt, vm);
    if (!auxLoggersOK) {
      return auxLoggersOK;
    }
  }

  if (sopt.rangeFactorizationBins > 0) {
    sopt.useRangeFactorization = true;
  }

  if (sopt.gcBiasCorrect and !sopt.biasCorrect) {
    sopt.numConditionalGCBins = 1;
  }

  std::transform(sopt.hitFilterPolicyStr.begin(), sopt.hitFilterPolicyStr.end(),
                 sopt.hitFilterPolicyStr.begin(), ::toupper);
  if ( sopt.hitFilterPolicyStr == "BEFORE" ) {
    sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::FILTER_BEFORE_CHAINING;
  } else if ( sopt.hitFilterPolicyStr == "AFTER" ) {
    sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::FILTER_AFTER_CHAINING;
  } else if ( sopt.hitFilterPolicyStr == "BOTH" ) {
    sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::FILTER_BEFORE_AND_AFTER_CHAINING;
  } else if ( sopt.hitFilterPolicyStr == "NONE" ) {
    sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::DO_NOT_FILTER;
  } else {
    jointLog->critical("The argument {} for --hitFilterPolicy is invalid. Valid options are "
                       "BEFORE, AFTER, BOTH and NONE.", sopt.hitFilterPolicyStr);
    jointLog->flush();
    return false;
  }

  try {
    conflicting_options(vm, "validateMappings", "noSA");
    conflicting_options(vm, "mimicBT2", "noSA");
    conflicting_options(vm, "mimicStrictBT2", "noSA");
    conflicting_options(vm, "hardFilter", "noSA");
  } catch (std::logic_error& e) {
    jointLog->critical("{}", e.what());
    jointLog->flush();
    return false;
  }

  if (sopt.disableSA) {
    jointLog->critical("Note: Alignment-free mapping (i.e. mapping without subsequent selective-alignment) "
                       "has not yet been throughly tested under the pufferfish-based index and using the "
                       "pufferfish-based mapping strategies.  Thus, disabling of selective-alignment "
                       "is not currently allowed.  We may, potentially explore re-enabling this option in future "
                       "versions of salmon.");
    jointLog->flush();
    return false;
  }

  try {
    conflicting_options(vm, "useVBOpt", "useEM");
  } catch (std::logic_error& e) {
    jointLog->critical("{}", e.what());
    jointLog->flush();
    return false;
  }
  if(sopt.useEM) {
    sopt.useVBOpt = false;
  }

  try {
    conflicting_options(vm, "perNucleotidePrior", "perTranscriptPrior");
  } catch (std::logic_error& e) {
    jointLog->critical("{}", e.what());
    jointLog->flush();
    return false;
  }
  if(sopt.perNucleotidePrior) {
    sopt.perTranscriptPrior = false;
  }

  if (sopt.useVBOpt and sopt.perNucleotidePrior and vm["vbPrior"].defaulted()) {
    sopt.vbPrior = 1e-5;
    jointLog->info("Using per-nucleotide prior with the default VB prior.  Setting the default prior to {}",sopt.vbPrior);
  }

  if (vm["maxHashResizeThreads"].defaulted()) {
    sopt.maxHashResizeThreads = sopt.numThreads;
    jointLog->info("setting maxHashResizeThreads to {}", sopt.maxHashResizeThreads);
  }

  try {
    conflicting_options(vm, "alignments", "eqclasses");
  } catch (std::logic_error& e) {
    jointLog->critical("{}", e.what());
    jointLog->flush();
    return false;
  }

  try {
    option_dependency(vm, "alignments", "targets");
  } catch (std::logic_error& e) {
    jointLog->critical("{}", e.what());
    jointLog->flush();
    return false;
  }

  if (sopt.numBurninFrags < sopt.numPreBurninFrags) {
    sopt.jointLog->warn("You set the number of burnin fragments "
                        "(--numAuxModelSamples) to be less than the number "
                        "of \n"
                        "pre-burnin fragments (--numPreAuxModelSamples), but "
                        "it must be at least as large.  The \n"
                        "number of pre-burnin fragments and burnin fragments "
                        "is being set to the same value "
                        "({})",
                        sopt.numBurninFrags);
    sopt.numPreBurninFrags = sopt.numBurninFrags;
  }

  if (sopt.incompatPrior < 1e-100 or sopt.incompatPrior == 0.0) {
    jointLog->info("Fragment incompatibility prior below threshold.  "
                   "Incompatible fragments will be ignored.");
    sopt.incompatPrior = salmon::math::LOG_0;
    sopt.ignoreIncompat = true;
  } else {
    sopt.incompatPrior = std::log(sopt.incompatPrior);
    sopt.ignoreIncompat = false;
  }

  if (sopt.dumpEqWeights and !sopt.dumpEq) {
    sopt.dumpEq = true;
    jointLog->info("You specified --dumpEqWeights, which implies --dumpEq; "
                   "that option has been enabled.");
  }

  if (sopt.numGibbsSamples > 0 and sopt.numBootstraps > 0) {
    jointLog->critical(
        "You cannot perform both Gibbs sampling and bootstrapping. "
        "Please choose one.");
    jointLog->flush();
    return false;
  }
  if (sopt.numGibbsSamples > 0) {
    if (!(sopt.thinningFactor >= 1)) {
      jointLog->critical(
          "The Gibbs sampling thinning factor (--thinningFactor) "
          "cannot be smaller than 1.");
      jointLog->flush();
      return false;
    }
  }

  if (sopt.noFragLengthDist and !sopt.noEffectiveLengthCorrection) {
    jointLog->critical(
        "You cannot enable --noFragLengthDist without "
        "also enabling --noEffectiveLengthCorrection; exiting.");
    jointLog->flush();
    return false;
  }

  if (sopt.noLengthCorrection) {
    bool anyBiasCorrect =
        sopt.gcBiasCorrect or sopt.biasCorrect or sopt.posBiasCorrect;
    if (anyBiasCorrect) {
      jointLog->critical(
          "Since bias correction relies on modifying "
          "effective lengths, you cannot enable bias "
          "correction simultaneously with the --noLengthCorrection "
          "option.");
      jointLog->flush();
      return false;
    }
  }

  if (sopt.dontExtrapolateCounts) {
    if (sopt.numGibbsSamples == 0) {
      sopt.jointLog->critical("You passed the --noExtrapolateCounts flag, "
                              "but are not using Gibbs sampling. "
                              "The fomer implies the latter.  Please enable "
                              "Gibbs sampling to use this flag.");
      return false;
    }
  }

  bool perModeValidate{true};
  if (sopt.quantMode == SalmonQuantMode::ALIGN) {
    perModeValidate = validateOptionsAlignment_(sopt, vm);
  } else if (sopt.quantMode == SalmonQuantMode::MAP) {
    perModeValidate = validateOptionsMapping_(sopt, vm);
  }

  return perModeValidate;
}

} // namespace salmon::utils
