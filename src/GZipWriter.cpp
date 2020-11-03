#include <ctime>
#include <fstream>
#include <numeric>
#include <memory>

#include "parallel_hashmap/phmap.h"
#include "cereal/archives/json.hpp"

#include "AlignmentLibrary.hpp"
#include "AlevinTypes.hpp"
#include "DistributionUtils.hpp"
#include "GZipWriter.hpp"
#include "ReadExperiment.hpp"
#include "ReadPair.hpp"
#include "SalmonOpts.hpp"
#include "SalmonUtils.hpp"
#include "UnpairedRead.hpp"
#include "TranscriptGroup.hpp"
#include "SingleCellProtocols.hpp"
#include "zstr.hpp"

GZipWriter::GZipWriter(const boost::filesystem::path path,
                       std::shared_ptr<spdlog::logger> logger)
    : path_(path), logger_(logger) {}

GZipWriter::~GZipWriter() {
  if (bsStream_) {
    bsStream_->reset();
  }
  if (cellEQStream_){
    cellEQStream_->reset();
  }
  if (cellDedupEQStream_){
    cellDedupEQStream_->reset();
  }
  if (umiGraphStream_){
    umiGraphStream_->reset();
  }
  if (countMatrixStream_) {
    countMatrixStream_->reset();
  }
  if (tierMatrixStream_) {
    tierMatrixStream_->reset();
  }
  if (meanMatrixStream_) {
    meanMatrixStream_->reset();
  }
  if (varMatrixStream_) {
    varMatrixStream_->reset();
  }
  if (arboMatrixStream_) {
    arboMatrixStream_->reset();
  }
  if (fullBootstrapMatrixStream_) {
    fullBootstrapMatrixStream_->reset();
  }
  if (bcBootNameStream_) {
    bcBootNameStream_.reset();
  }
  if (bcFeaturesStream_) {
    bcFeaturesStream_.reset();
  }
  if (bcNameStream_) {
    bcNameStream_.reset();
  }
}

void GZipWriter::close_all_streams(){
  if (bsStream_) {
    bsStream_->reset();
  }
  if (cellEQStream_){
    cellEQStream_->reset();
  }
  if (cellDedupEQStream_){
    cellDedupEQStream_->reset();
  }
  if (umiGraphStream_){
    umiGraphStream_->reset();
  }
  if (countMatrixStream_) {
    countMatrixStream_->reset();
  }
  if (tierMatrixStream_) {
    tierMatrixStream_->reset();
  }
  if (meanMatrixStream_) {
    meanMatrixStream_->reset();
  }
  if (varMatrixStream_) {
    varMatrixStream_->reset();
  }
  if (arboMatrixStream_) {
    arboMatrixStream_->reset();
  }
  if (fullBootstrapMatrixStream_) {
    fullBootstrapMatrixStream_->reset();
  }
  if (bcBootNameStream_) {
    bcBootNameStream_->close();
  }
  if (bcFeaturesStream_) {
    bcFeaturesStream_->close();
  }
  if (bcNameStream_) {
    bcNameStream_->close();
  }
}

/**
 * Creates a new gzipped file (path) and writes the contents
 * of the vector (vec) to the file in binary.
 */
template <typename T>
bool writeVectorToFile(boost::filesystem::path path,
                       const std::vector<T>& vec) {

  {
    bool binary = std::is_same<T, std::string>::value;
    auto flags = std::ios_base::out | std::ios_base::binary;

    boost::iostreams::filtering_ostream out;
    out.push(boost::iostreams::gzip_compressor(6));
    out.push(boost::iostreams::file_sink(path.string(), flags));

    size_t num = vec.size();
    size_t elemSize = sizeof(typename std::vector<T>::value_type);
    // We have to get rid of constness below, but this should be OK
    out.write(reinterpret_cast<char*>(const_cast<T*>(vec.data())),
              num * elemSize);
    out.reset();
  }
  return true;
}

/**
 * Write the equivalence class information to file.
 * The header will contain the transcript / target ids in
 * a fixed order, then each equivalence class will consist
 * of a line / row.
 */
template <typename ExpT>
bool GZipWriter::writeEquivCounts(const SalmonOpts& opts, ExpT& experiment) {
  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / opts.auxDir;
  bool auxSuccess = boost::filesystem::create_directories(auxDir);

  // write gzipped equivalence classes
  bfs::path eqFilePath = auxDir / "eq_classes.txt.gz";
  std::unique_ptr<std::ostream> equivFilePtr(new zstr::ofstream(eqFilePath.string()));

  // write uncompressed equivalence classes
  //bfs::path eqFilePath = auxDir / "eq_classes.txt";
  //std::unique_ptr<std::ostream> equivFilePtr(new std::ofstream(eqFilePath.string()));
  //std::ofstream equivFile(eqFilePath.string());

  auto& transcripts = experiment.transcripts();
  auto& eqBuilder = experiment.equivalenceClassBuilder();
  auto& eqVec = eqBuilder.eqVec();
  size_t numEqClasses = eqVec.size();
  bool dumpRichWeights = opts.dumpEqWeights;

  // we need this scope, but will fill it in only if we need to
  phmap::flat_hash_map<TranscriptGroup, uint64_t, TranscriptGroupHasher> collapsedMap;

  if (!dumpRichWeights) {
    logger_->info("Collapsing factorization information into simplified equivalence classes.");

    // if we are using range-factorization, but don't want weights,
    // collapse the equivalence classes into naive ones.
    collapsedMap.reserve(eqVec.size());
    for (size_t eqIdx = 0; eqIdx < eqVec.size(); ++eqIdx) {
      auto& eq = eqVec[eqIdx];
      const TranscriptGroup& tgroup = eq.first;
      const std::vector<uint32_t>& txps = tgroup.txps;
      uint64_t count = eq.second.count;
      const uint32_t groupSize = eqBuilder.getNumTranscriptsForClass(eqIdx);

      // make a key from the IDs, and copy over to
      std::vector<uint32_t> txpIDs(groupSize);
      auto bit = txps.begin();
      auto eit = txps.begin()+groupSize;
      std::copy(bit, eit, txpIDs.begin());

      TranscriptGroup key(txpIDs);
      collapsedMap[key] += count;
    }

    numEqClasses = collapsedMap.size();
    logger_->info("done.");
  }

  // Number of transcripts
  (*equivFilePtr) << transcripts.size() << '\n';

  // Number of equivalence classes
  (*equivFilePtr) << numEqClasses << '\n';

  // Transcript names
  for (auto& t : transcripts) {
    (*equivFilePtr) << t.RefName << '\n';
  }

  // If we are dumping eq weights, then just go over what
  // we have
  if (dumpRichWeights) {

    for (size_t eqIdx = 0; eqIdx < eqVec.size(); ++eqIdx) {
      auto& eq = eqVec[eqIdx];
      uint64_t count = eq.second.count;
      const TranscriptGroup& tgroup = eq.first;
      const std::vector<uint32_t>& txps = tgroup.txps;
      const auto& auxs = eq.second.combinedWeights;
      const uint32_t groupSize = eqBuilder.getNumTranscriptsForClass(eqIdx);

      (*equivFilePtr) << groupSize << '\t';
      for (uint32_t i = 0; i < groupSize; ++i) {
        (*equivFilePtr) << txps[i] << '\t';
      }
      for (uint32_t i = 0; i < groupSize; ++i) {
        (*equivFilePtr) << auxs[i] << '\t';
      }
      (*equivFilePtr) << count << '\n';
    }
  } else {
    // Otherwise, go over the collapsed map with the simplified
    // classes.

    // dump the size, txp ids and count
    for (auto&& kv : collapsedMap) {
      auto& txps = kv.first.txps;
      auto& count = kv.second;
      const uint32_t groupSize = txps.size();
      (*equivFilePtr) << groupSize << '\t';
      for (uint32_t i = 0; i < groupSize; ++i) {
        (*equivFilePtr) << txps[i] << '\t';
      }
      (*equivFilePtr) << count << '\n';
    }
  }

  // NOTE: no close via the ostream interface, will happen upon destruction.
  // equivFilePtr->close();
  return true;
}

// FIXME(@k3yavi): We should probably be writing down the umiLength in the BFH files
// somewhere.
template <typename ExpT>
bool GZipWriter::writeBFH(boost::filesystem::path& outDir,
                          ExpT& experiment, size_t /*umiLength*/,
                          std::vector<std::string>& bcSeqVec) {
  namespace bfs = boost::filesystem;

  bfs::path eqFilePath = outDir / "bfh.txt";
  std::ofstream equivFile(eqFilePath.string());

  auto& transcripts = experiment.transcripts();
  const auto& eqVec =
    experiment.equivalenceClassBuilder().eqMap().lock_table();

  // Number of transcripts
  equivFile << transcripts.size() << '\n';

  // Number of barcodes
  equivFile << bcSeqVec.size() << '\n';

  // Number of equivalence classes
  equivFile << eqVec.size() << '\n';

  for (auto& t : transcripts) {
    equivFile << t.RefName << '\n';
  }

  for (auto& b : bcSeqVec) {
    equivFile << b << '\n';
  }

  alevin::types::AlevinUMIKmer umiObj;

  for (auto& eq : eqVec) {
    uint64_t count = eq.second.count;
    // for each transcript in this class
    const TranscriptGroup& tgroup = eq.first;
    const std::vector<uint32_t>& txps = tgroup.txps;

    // group size
    equivFile << txps.size() << '\t';
    // each group member
    for (auto tid : txps) { equivFile << tid << '\t'; }
    const auto& bgroup = eq.second.barcodeGroup;
    equivFile << count << "\t" << bgroup.size();
    for (auto  bcIt : bgroup){
      auto bc = bcIt.first;
      auto ugroup = bcIt.second;
      equivFile << "\t" << bc << "\t" << ugroup.size();
      for (auto umiIt : ugroup){
        auto umi = umiIt.first;
        umiObj.word__(0) = umi;
        auto count = umiIt.second;

        std::string s = umiObj.toStr();
        equivFile << "\t" << s << "\t" << count;
      }
    }
    equivFile << "\n";
  }
  equivFile.close();
  return true;
}


/**
 * Write the equivalence class information to file.
 * The header will contain the transcript / target ids in
 * a fixed order, then each equivalence class will consist
 * of a line / row.
 */
template <typename ExpT, typename ProtocolT>
bool GZipWriter::writeEquivCounts(
    const AlevinOpts<ProtocolT>& aopts,
    ExpT& experiment) {

  namespace bfs = boost::filesystem;

  bfs::path eqFilePath = aopts.outputDirectory / "cell_eq_info.txt";
  std::ofstream equivFile(eqFilePath.string());

  auto& transcripts = experiment.transcripts();
  const auto& eqVec =
    experiment.equivalenceClassBuilder().eqMap().lock_table();

  // Number of transcripts
  equivFile << transcripts.size() << '\n';

  // Number of equivalence classes
  equivFile << eqVec.size() << '\n';

  for (auto& t : transcripts) {
    equivFile << t.RefName << '\n';
  }
  for (auto& eq : eqVec) {
    uint64_t count = eq.second.count;
    // for each transcript in this class
    const TranscriptGroup& tgroup = eq.first;
    const std::vector<uint32_t>& txps = tgroup.txps;

    // group size
    equivFile << txps.size() ;
    // each group member
    for (auto tid : txps) { equivFile << '\t' << tid; }
    equivFile << '\n';
  }

  equivFile.close();
  return true;
}

template <typename ExpT>
std::vector<std::string> getLibTypeStrings(const ExpT& experiment) {
  auto& libs = experiment.readLibraries();
  std::vector<std::string> libStrings;
  for (auto& rl : libs) {
    libStrings.push_back(rl.getFormat().toString());
  }
  return libStrings;
}

template <typename AlnT, typename EQBuilderT>
std::vector<std::string>
getLibTypeStrings(const AlignmentLibrary<AlnT, EQBuilderT>& experiment) {
  std::vector<std::string> libStrings;
  libStrings.push_back(experiment.format().toString());
  return libStrings;
}

/**
 * Write the ``main'' metadata to file when the quantifications will be empty.
 */
template <typename ExpT>
bool GZipWriter::writeEmptyMeta(const SalmonOpts& opts, const ExpT& experiment,
                                std::vector<std::string>& errors) {
  namespace bfs = boost::filesystem;
  using salmon::utils::DuplicateTargetStatus;

  bfs::path auxDir = path_ / opts.auxDir;
  bool auxSuccess = boost::filesystem::create_directories(auxDir);

  bfs::path info = auxDir / "meta_info.json";

  {
    std::ofstream os(info.string());
    cereal::JSONOutputArchive oa(os);

    std::string sampType = "none";
    std::string optType = "none";

    auto& transcripts = experiment.transcripts();
    oa(cereal::make_nvp("salmon_version", std::string(salmon::version)));
    oa(cereal::make_nvp("samp_type", sampType));
    oa(cereal::make_nvp("opt_type", optType));

    oa(cereal::make_nvp("quant_errors", errors));
    auto libStrings = getLibTypeStrings(experiment);
    oa(cereal::make_nvp("num_libraries", libStrings.size()));
    oa(cereal::make_nvp("library_types", libStrings));

    oa(cereal::make_nvp("frag_dist_length", 0));
    oa(cereal::make_nvp("frag_length_mean", 0.0));
    oa(cereal::make_nvp("frag_length_sd", 0.0));
    oa(cereal::make_nvp("seq_bias_correct", false));
    oa(cereal::make_nvp("gc_bias_correct", false));
    oa(cereal::make_nvp("num_bias_bins", 0));

    std::string mapTypeStr = opts.alnMode ? "alignment" : "mapping";
    oa(cereal::make_nvp("mapping_type", mapTypeStr));

    auto has_dups = experiment.index_retains_duplicates();
    switch(has_dups) {
      case DuplicateTargetStatus::RETAINED_DUPLICATES:
        oa(cereal::make_nvp("keep_duplicates", true));
        break;
      case DuplicateTargetStatus::REMOVED_DUPLICATES:
        oa(cereal::make_nvp("keep_duplicates", false));
        break;
      case DuplicateTargetStatus::UNKNOWN:
      default:
        break;
    }

    auto numValidTargets = transcripts.size();
    auto numDecoys = experiment.getNumDecoys();
    oa(cereal::make_nvp("num_targets", numValidTargets));
    oa(cereal::make_nvp("num_decoy_targets", numDecoys));

    // True if we dumped the equivalence classes, false otherwise
    oa(cereal::make_nvp("serialized_eq_classes", false));
    oa(cereal::make_nvp("num_eq_classes", 0));

    // For now, this vector is empty unless we dumped the equivalence classes
    // with weights.  In which case it contains the string "scalar_weights".
    std::vector<std::string> props;

    bool isRangeFactorizationOn = (opts.rangeFactorizationBins > 0);
    bool dumpRichWeights = opts.dumpEqWeights;
    if(isRangeFactorizationOn){
      props.push_back("range_factorized") ;
    }
    if (dumpRichWeights) {
      props.push_back("scalar_weights");
    }
    // write down that we are gzipping the eq classes
    props.push_back("gzipped");

    oa(cereal::make_nvp("eq_class_properties", props));

    oa(cereal::make_nvp("length_classes", experiment.getLengthQuantiles()));
    oa(cereal::make_nvp("index_seq_hash", experiment.getIndexSeqHash256()));
    oa(cereal::make_nvp("index_name_hash", experiment.getIndexNameHash256()));
    oa(cereal::make_nvp("index_seq_hash512", experiment.getIndexSeqHash512()));
    oa(cereal::make_nvp("index_name_hash512", experiment.getIndexNameHash512()));
    oa(cereal::make_nvp("index_decoy_seq_hash", experiment.getIndexDecoySeqHash256()));
    oa(cereal::make_nvp("index_decoy_name_hash", experiment.getIndexDecoyNameHash256()));
    oa(cereal::make_nvp("num_bootstraps", 0));
    oa(cereal::make_nvp("num_processed", experiment.numObservedFragments()));
    oa(cereal::make_nvp("num_mapped", experiment.numMappedFragments()));
    oa(cereal::make_nvp("percent_mapped",
                        experiment.effectiveMappingRate() * 100.0));
    oa(cereal::make_nvp("call", std::string("quant")));
    oa(cereal::make_nvp("start_time", opts.runStartTime));
    oa(cereal::make_nvp("end_time", opts.runStopTime));
  }
  return true;
}

template <typename ProtocolT>
bool GZipWriter::writeMetaAlevin(const AlevinOpts<ProtocolT>& opts,
                                 boost::filesystem::path aux_dir) {
  namespace bfs = boost::filesystem;
  bfs::path auxDir = path_ / aux_dir;
  bfs::path info = auxDir / "alevin_meta_info.json";

  {
    std::ofstream os(info.string());
    cereal::JSONOutputArchive oa(os);

    oa(cereal::make_nvp("total_reads", opts.totalReads));
    oa(cereal::make_nvp("reads_with_N",
                        opts.totalReads - opts.totalUsedReads));
    oa(cereal::make_nvp("noisy_cb_reads", opts.readsThrown));
    oa(cereal::make_nvp("noisy_umi_reads", opts.noisyUmis));
    oa(cereal::make_nvp("used_reads", opts.totalUsedReads
                        - opts.readsThrown
                        - opts.noisyUmis));
    oa(cereal::make_nvp("mapping_rate", opts.mappingRate));
    oa(cereal::make_nvp("reads_in_eqclasses", opts.eqReads));

    oa(cereal::make_nvp("total_cbs", opts.totalCBs));
    oa(cereal::make_nvp("used_cbs", opts.totalUsedCBs));

    if(not boost::filesystem::exists(opts.whitelistFile)){
      oa(cereal::make_nvp("initial_whitelist", opts.kneeCutoff));
      oa(cereal::make_nvp("low_conf_cbs", opts.totalLowConfidenceCBs));
      oa(cereal::make_nvp("num_features", opts.numFeatures));
    }
    if(opts.numNoMapCB > 0) {
      oa(cereal::make_nvp("no_read_mapping_cbs", opts.numNoMapCB));
    }
    oa(cereal::make_nvp("final_num_cbs", opts.intelligentCutoff));

    if (opts.numBootstraps > 0) {
      oa(cereal::make_nvp("num_bootstraps",
                          opts.numBootstraps));
    }

    oa(cereal::make_nvp("deduplicated_umis", opts.totalDedupUMIs));
    oa(cereal::make_nvp("mean_umis_per_cell",
                        opts.totalDedupUMIs / opts.intelligentCutoff));
    oa(cereal::make_nvp("mean_genes_per_cell",
                        opts.totalExpGenes / opts.intelligentCutoff));
  }

  return true;
}

/**
 * Write the ``main'' metadata to file.  Currently this includes:
 *   -- Names of the target id's if bootstrapping / gibbs is performed
 *   -- The fragment length distribution
 *   -- The expected and observed bias values
 *   -- A json file with information about the run
 */
template <typename ExpT>
bool GZipWriter::writeMeta(const SalmonOpts& opts, const ExpT& experiment, const MappingStatistics& mstats) {

  namespace bfs = boost::filesystem;
  using salmon::utils::DuplicateTargetStatus;

  bfs::path auxDir = path_ / opts.auxDir;
  bool auxSuccess = boost::filesystem::create_directories(auxDir);

  auto numBootstraps = opts.numBootstraps;
  auto numSamples = (numBootstraps > 0) ? numBootstraps : opts.numGibbsSamples;
  if (numSamples > 0) {
    bsPath_ = auxDir / "bootstrap";
    bool bsSuccess = boost::filesystem::create_directories(bsPath_);
    {

      boost::iostreams::filtering_ostream nameOut;
      nameOut.push(boost::iostreams::gzip_compressor(6));
      auto bsFilename = bsPath_ / "names.tsv.gz";
      nameOut.push(
          boost::iostreams::file_sink(bsFilename.string(), std::ios_base::out));

      auto& transcripts = experiment.transcripts();
      size_t numTxps = transcripts.size();
      if (numTxps == 0) {
        return false;
      }
      for (size_t tn = 0; tn < numTxps; ++tn) {
        auto& t = transcripts[tn];
        nameOut << t.RefName;
        if (tn < numTxps - 1) {
          nameOut << '\t';
        }
      }
      nameOut << '\n';
      nameOut.reset();
    }
  }

  bfs::path fldPath = auxDir / "fld.gz";
  int32_t numFLDSamples{10000};
  auto fldSummary = distribution_utils::samplesFromLogPMF(
      experiment.fragmentLengthDistribution(), numFLDSamples);
  writeVectorToFile(fldPath, fldSummary.samples);

  bfs::path normBiasPath = auxDir / "expected_bias.gz";
  writeVectorToFile(normBiasPath, experiment.expectedSeqBias());

  bfs::path obsBiasPath = auxDir / "observed_bias.gz";
  // TODO: dump both sense and anti-sense models
  const auto& bcounts =
      experiment.readBias(salmon::utils::Direction::FORWARD).counts;
  std::vector<int32_t> observedBias(bcounts.size(), 0);
  std::copy(bcounts.begin(), bcounts.end(), observedBias.begin());
  writeVectorToFile(obsBiasPath, observedBias);

  bfs::path obsBiasPath3p = auxDir / "observed_bias_3p.gz";
  const auto& bcounts3p =
      experiment.readBias(salmon::utils::Direction::REVERSE_COMPLEMENT).counts;
  std::vector<int32_t> observedBias3p(bcounts3p.size(), 0);
  std::copy(bcounts3p.begin(), bcounts3p.end(), observedBias3p.begin());
  writeVectorToFile(obsBiasPath3p, observedBias3p);

  if (opts.biasCorrect) {
    // 5' observed
    {
      bfs::path obs5Path = auxDir / "obs5_seq.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(obs5Path.string(), flags));
      auto& obs5 =
          experiment.readBiasModelObserved(salmon::utils::Direction::FORWARD);
      obs5.writeBinary(out);
    }
    // 3' observed
    {
      bfs::path obs3Path = auxDir / "obs3_seq.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(obs3Path.string(), flags));
      auto& obs3 = experiment.readBiasModelObserved(
          salmon::utils::Direction::REVERSE_COMPLEMENT);
      obs3.writeBinary(out);
    }

    // 5' expected
    {
      bfs::path exp5Path = auxDir / "exp5_seq.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(exp5Path.string(), flags));
      auto& exp5 =
          experiment.readBiasModelExpected(salmon::utils::Direction::FORWARD);
      exp5.writeBinary(out);
    }
    // 3' expected
    {
      bfs::path exp3Path = auxDir / "exp3_seq.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(exp3Path.string(), flags));
      auto& exp3 = experiment.readBiasModelExpected(
          salmon::utils::Direction::REVERSE_COMPLEMENT);
      exp3.writeBinary(out);
    }
  }

  if (opts.gcBiasCorrect) {
    // GC observed
    {
      bfs::path obsGCPath = auxDir / "obs_gc.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(obsGCPath.string(), flags));
      auto& obsgc = experiment.observedGC();
      obsgc.writeBinary(out);
    }
    // GC expected
    {
      bfs::path expGCPath = auxDir / "exp_gc.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(expGCPath.string(), flags));
      auto& expgc = experiment.expectedGCBias();
      expgc.writeBinary(out);
    }
  }

  if (opts.posBiasCorrect) {
    // the length classes
    const auto& lenBounds = experiment.getLengthQuantiles();

    // lambda to write out a vector of SimplePosBias models (along with the
    // length bounds) to file.
    auto writePosModel =
        [&lenBounds, this](bfs::path fpath,
                           const std::vector<SimplePosBias>& model) -> bool {
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(fpath.string(), flags));
      // Write out the number of different models
      uint32_t numModels = static_cast<uint32_t>(lenBounds.size());
      out.write(reinterpret_cast<char*>(&numModels), sizeof(numModels));
      // Write out the length class for each model
      for (const auto& b : lenBounds) {
        out.write(reinterpret_cast<char*>(const_cast<uint32_t*>(&b)),
                  sizeof(b));
      }
      // write out each
      for (auto& pb : model) {
        bool success = pb.writeBinary(out);
        if (!success) {
          this->logger_->error(
              "Could not write out positional bias model to {}!",
              fpath.string());
        }
      }
      return true;
    };

    // 5' observed
    {
      bfs::path obsPosPath = auxDir / "obs5_pos.gz";
      // Get the pos bias vector
      auto& posBiases = experiment.posBias(salmon::utils::Direction::FORWARD);
      writePosModel(obsPosPath, posBiases);
    }
    // 3' observed
    {
      bfs::path obsPosPath = auxDir / "obs3_pos.gz";
      // Get the pos bias vector
      auto& posBiases =
          experiment.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
      writePosModel(obsPosPath, posBiases);
    }
    // 5' expected
    {
      bfs::path expPosPath = auxDir / "exp5_pos.gz";
      // Get the pos bias vector
      auto& posBiases =
          experiment.posBiasExpected(salmon::utils::Direction::FORWARD);
      writePosModel(expPosPath, posBiases);
    }
    // 3' expected
    {
      bfs::path expPosPath = auxDir / "exp3_pos.gz";
      // Get the pos bias vector
      auto& posBiases = experiment.posBiasExpected(
          salmon::utils::Direction::REVERSE_COMPLEMENT);
      writePosModel(expPosPath, posBiases);
    }
  }

  /*
  bfs::path normGCPath = auxDir / "expected_gc.gz";
  writeVectorToFile(normGCPath, experiment.expectedGCBias());

  bfs::path obsGCPath = auxDir / "observed_gc.gz";
  const auto& gcCounts = experiment.observedGC();
  std::vector<double> observedGC(gcCounts.size(), 0.0);
  std::copy(gcCounts.begin(), gcCounts.end(), observedGC.begin());
  writeVectorToFile(obsGCPath, observedGC);
  */

  bfs::path info = auxDir / "meta_info.json";

  {
    std::ofstream os(info.string());
    cereal::JSONOutputArchive oa(os);

    std::string sampType = "none";
    if (numBootstraps == 0 and numSamples > 0) {
      sampType = "gibbs";
    }
    if (numBootstraps > 0) {
      sampType = "bootstrap";
    }

    auto& transcripts = experiment.transcripts();
    oa(cereal::make_nvp("salmon_version", std::string(salmon::version)));
    oa(cereal::make_nvp("samp_type", sampType));

    std::string optType = "none";
    if (opts.useEM) {
      optType = "em";
    } else {
      optType = "vb";
    }
    oa(cereal::make_nvp("opt_type", optType));

    std::vector<std::string> errors;
    oa(cereal::make_nvp("quant_errors", errors));

    auto libStrings = getLibTypeStrings(experiment);
    oa(cereal::make_nvp("num_libraries", libStrings.size()));
    oa(cereal::make_nvp("library_types", libStrings));

    oa(cereal::make_nvp("frag_dist_length", fldSummary.samples.size()));
    oa(cereal::make_nvp("frag_length_mean", fldSummary.mean));
    oa(cereal::make_nvp("frag_length_sd", fldSummary.sd));
    oa(cereal::make_nvp("seq_bias_correct", opts.biasCorrect));
    oa(cereal::make_nvp("gc_bias_correct", opts.gcBiasCorrect));
    oa(cereal::make_nvp("num_bias_bins", bcounts.size()));

    std::string mapTypeStr = opts.alnMode ? "alignment" : "mapping";
    oa(cereal::make_nvp("mapping_type", mapTypeStr));

    auto has_dups = experiment.index_retains_duplicates();
    switch(has_dups) {
      case DuplicateTargetStatus::RETAINED_DUPLICATES:
        oa(cereal::make_nvp("keep_duplicates", true));
        break;
      case DuplicateTargetStatus::REMOVED_DUPLICATES:
        oa(cereal::make_nvp("keep_duplicates", false));
        break;
      case DuplicateTargetStatus::UNKNOWN:
      default:
        break;
    }

    auto numValidTargets = transcripts.size();
    auto numDecoys = experiment.getNumDecoys();
    oa(cereal::make_nvp("num_valid_targets", numValidTargets));
    oa(cereal::make_nvp("num_decoy_targets", numDecoys));

    auto& eqBuilder = const_cast<ExpT&>(experiment).equivalenceClassBuilder();
    auto num_eq_classes_opt = eqBuilder.numEqClasses();
    size_t num_eq_classes = num_eq_classes_opt ? (*num_eq_classes_opt) : 0;
    oa(cereal::make_nvp("num_eq_classes", num_eq_classes));

    // True if we dumped the equivalence classes, false otherwise
    oa(cereal::make_nvp("serialized_eq_classes", opts.dumpEq));

    // For now, this vector is empty unless we dumped the equivalence classes
    // with weights.  In which case it contains the string "scalar_weights".
    std::vector<std::string> props;
    bool isRangeFactorizationOn = opts.rangeFactorizationBins ;
    bool dumpRichWeights = opts.dumpEqWeights;
    if(isRangeFactorizationOn){
      props.push_back("range_factorized") ;
    }
    if (opts.dumpEqWeights) {
      props.push_back("scalar_weights");
    }
    // write down that we are gzipping the eq classes
    props.push_back("gzipped");

    oa(cereal::make_nvp("eq_class_properties", props));

    oa(cereal::make_nvp("length_classes", experiment.getLengthQuantiles()));
    oa(cereal::make_nvp("index_seq_hash", experiment.getIndexSeqHash256()));
    oa(cereal::make_nvp("index_name_hash", experiment.getIndexNameHash256()));
    oa(cereal::make_nvp("index_seq_hash512", experiment.getIndexSeqHash512()));
    oa(cereal::make_nvp("index_name_hash512", experiment.getIndexNameHash512()));
    oa(cereal::make_nvp("index_decoy_seq_hash", experiment.getIndexDecoySeqHash256()));
    oa(cereal::make_nvp("index_decoy_name_hash", experiment.getIndexDecoyNameHash256()));
    oa(cereal::make_nvp("num_bootstraps", numSamples));
    oa(cereal::make_nvp("num_processed", experiment.numObservedFragments()));
    oa(cereal::make_nvp("num_mapped", experiment.numMappedFragments()));
    oa(cereal::make_nvp("num_decoy_fragments", mstats.numDecoyFragments.load()));
    oa(cereal::make_nvp("num_dovetail_fragments", mstats.numDovetails.load()));
    oa(cereal::make_nvp("num_fragments_filtered_vm", mstats.numFragmentsFiltered.load()));
    oa(cereal::make_nvp("num_alignments_below_threshold_for_mapped_fragments_vm",
                        mstats.numMappingsFiltered.load()));
    oa(cereal::make_nvp("percent_mapped",
                        experiment.effectiveMappingRate() * 100.0));
    oa(cereal::make_nvp("call", std::string("quant")));
    oa(cereal::make_nvp("start_time", opts.runStartTime));
    oa(cereal::make_nvp("end_time", opts.runStopTime));
  }

  {
    bfs::path ambigInfo = auxDir / "ambig_info.tsv";
    std::ofstream os(ambigInfo.string());
    os << "UniqueCount\tAmbigCount\n";

    auto& transcripts = experiment.transcripts();
    auto& eqBuilder = const_cast<ExpT&>(experiment).equivalenceClassBuilder();
    auto& eqVec = eqBuilder.eqVec();

    class CountPair {
    public:
      uint32_t unique = 0;
      uint32_t potential = 0;
    };

    std::vector<CountPair> counts(transcripts.size());
    for (size_t eqIdx = 0; eqIdx < eqVec.size(); ++eqIdx) {
      auto& eq = eqVec[eqIdx];

      // NOTE: this function is safe only if we've populated eqVec_ member.
      auto groupSize = eqBuilder.getNumTranscriptsForClass(eqIdx);
      uint64_t count = eq.second.count;
      const TranscriptGroup& tgroup = eq.first;
      const std::vector<uint32_t>& txps = tgroup.txps;

      if (groupSize > 1) {
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          counts[tid].potential += count;
        }
      } else {
        counts[txps.front()].unique += count;
      }
    }
    for (size_t i = 0; i < transcripts.size(); ++i) {
      os << counts[i].unique << '\t' << counts[i].potential << '\n';
    }
  }

  return true;
}

bool GZipWriter::writeAbundances(
                                 std::vector<double>& alphas,
                                 std::vector<Transcript>& transcripts) {
  namespace bfs = boost::filesystem;

  bfs::path fname = path_ / "quant.sf";
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(fname.c_str(), "w"), std::fclose);

  fmt::print(output.get(), "Name\tLength\tNumMolecules\n");

  // Now posterior has the transcript fraction
  for (size_t i=0; i < transcripts.size(); i++) {
    fmt::print(output.get(), "{}\t{}\t{}\n",
               transcripts[i].RefName,
               transcripts[i].CompleteLength,
               alphas[i]);
  }
  return true;
}

bool GZipWriter::writeBootstraps(std::string& bcName,
                                 std::vector<double>& alphas,
                                 std::vector<double>& variance,
                                 bool useAllBootstraps,
                                 std::vector<std::vector<double>>& sampleEstimates){
#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  namespace bfs = boost::filesystem;
  if (!meanMatrixStream_) {
    meanMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    meanMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto meanMatFilename = path_ / "alevin" / "quants_mean_mat.gz";
    meanMatrixStream_->push(boost::iostreams::file_sink(meanMatFilename.string(),
                                                         std::ios_base::out | std::ios_base::binary));

    varMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    varMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto varMatFilename = path_ / "alevin" / "quants_var_mat.gz";
    varMatrixStream_->push(boost::iostreams::file_sink(varMatFilename.string(),
                                                       std::ios_base::out | std::ios_base::binary));

    auto bcBootNameFilename = path_ / "alevin" / "quants_boot_rows.txt";
    bcBootNameStream_.reset(new std::ofstream);
    bcBootNameStream_->open(bcBootNameFilename.string());
  }

  if (useAllBootstraps and !fullBootstrapMatrixStream_) {
    auto fullBootstrapsFilename = path_ / "alevin" / "quants_boot_mat.gz";
    fullBootstrapMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    fullBootstrapMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    fullBootstrapMatrixStream_->push(boost::iostreams::file_sink(fullBootstrapsFilename.string(),
                                                      std::ios_base::out | std::ios_base::binary));
  }

  boost::iostreams::filtering_ostream& countfile = *meanMatrixStream_;
  boost::iostreams::filtering_ostream& varfile = *varMatrixStream_;
  std::ofstream& namefile = *bcBootNameStream_;

  size_t num = alphas.size();
  if (alphas.size() != variance.size()){
    std::cerr<<"ERROR: Quants matrix and varicance matrix size differs"<<std::flush;
    exit(74);
  }
  size_t elSize = sizeof(typename std::vector<double>::value_type);
  countfile.write(reinterpret_cast<char*>(alphas.data()),
                  elSize * num);
  varfile.write(reinterpret_cast<char*>(variance.data()),
                  elSize * num);

  if (useAllBootstraps){
    boost::iostreams::filtering_ostream& fullBootsfile = *fullBootstrapMatrixStream_;
    for(auto& sample: sampleEstimates){
      fullBootsfile.write(reinterpret_cast<char*>(sample.data()),
                          elSize * num);
    }
  }

  namefile << std::endl;
  namefile.write(bcName.c_str(), bcName.size());
  return true;
}

bool GZipWriter::writeSparseBootstraps(std::string& bcName,
                                       std::vector<double>& alphas,
                                       std::vector<double>& variance,
                                       bool useAllBootstraps,
                                       std::vector<std::vector<double>>& sampleEstimates){
  size_t num = alphas.size();
  std::vector<float> meanSparse, varSparse;
  std::vector<std::vector<float>> bootSparse;
  meanSparse.reserve(num/2);
  varSparse.reserve(num/2);
  if (useAllBootstraps){ bootSparse.reserve(num/2); }

  std::vector<uint8_t> meanFlag, varFlag;
  std::vector<std::vector<uint8_t>> bootFlag;
  meanFlag.reserve(static_cast<size_t>(std::ceil(num/8.0)));
  varFlag.reserve(static_cast<size_t>(std::ceil(num/8.0)));
  if (useAllBootstraps) { bootFlag.reserve(num); }

  size_t elSize = sizeof(decltype(meanSparse)::value_type);
  size_t flagSize = sizeof(decltype(meanFlag)::value_type);

  for (size_t i=0; i<num; i+=8) {
    uint8_t flag {0};
    for (size_t j=0; j<8; j++) {
      size_t vectorIndex = i+j;
      if (vectorIndex >= num) { break; }

      if (alphas[vectorIndex] > std::numeric_limits<float>::min()) {
        meanSparse.emplace_back(alphas[vectorIndex]);
        flag |= 128 >> j;
      }
    }
    meanFlag.emplace_back(flag);
  }

  for (size_t i=0; i<num; i+=8) {
    uint8_t flag {0};
    for (size_t j=0; j<8; j++) {
      size_t vectorIndex = i+j;
      if (vectorIndex >= num) { break; }

      if (variance[vectorIndex] > std::numeric_limits<uint8_t>::min()) {
        varSparse.emplace_back(variance[vectorIndex]);
        flag |= 128 >> j;
      }
    }
    varFlag.emplace_back(flag);
  }

  if (useAllBootstraps) {
    for(auto& sample: sampleEstimates){
      std::vector<float> sampleSparse;
      std::vector<uint8_t> sampleFlag;
      for (size_t i=0; i<num; i+=8) {
        uint8_t flag {0};
        for (size_t j=0; j<8; j++) {
          size_t vectorIndex = i+j;
          if (vectorIndex >= num) { break; }

          if (sample[vectorIndex] > std::numeric_limits<uint8_t>::min()) {
            sampleSparse.emplace_back(sample[vectorIndex]);
            flag |= 128 >> j;
          }
        }
        sampleFlag.emplace_back(flag);
      }

      bootFlag.emplace_back(sampleFlag);
      bootSparse.emplace_back(sampleSparse);
    }
  }

#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  namespace bfs = boost::filesystem;
  if (!meanMatrixStream_) {
    meanMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    meanMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto meanMatFilename = path_ / "alevin" / "quants_mean_mat.gz";
    meanMatrixStream_->push(boost::iostreams::file_sink(meanMatFilename.string(),
                                                         std::ios_base::out | std::ios_base::binary));

    varMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    varMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto varMatFilename = path_ / "alevin" / "quants_var_mat.gz";
    varMatrixStream_->push(boost::iostreams::file_sink(varMatFilename.string(),
                                                       std::ios_base::out | std::ios_base::binary));

    auto bcBootNameFilename = path_ / "alevin" / "quants_boot_rows.txt";
    bcBootNameStream_.reset(new std::ofstream);
    bcBootNameStream_->open(bcBootNameFilename.string());
  }

  if (useAllBootstraps and !fullBootstrapMatrixStream_) {
    auto fullBootstrapsFilename = path_ / "alevin" / "quants_boot_mat.gz";
    fullBootstrapMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    fullBootstrapMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    fullBootstrapMatrixStream_->push(boost::iostreams::file_sink(fullBootstrapsFilename.string(),
                                                      std::ios_base::out | std::ios_base::binary));
  }

  boost::iostreams::filtering_ostream& countfile = *meanMatrixStream_;
  boost::iostreams::filtering_ostream& varfile = *varMatrixStream_;
  std::ofstream& namefile = *bcBootNameStream_;

  if (alphas.size() != variance.size()){
    std::cerr<<"ERROR: Quants matrix and varicance matrix size differs"<<std::flush;
    exit(74);
  }

  countfile.write(reinterpret_cast<char*>(meanFlag.data()),
                  flagSize * meanFlag.size());
  countfile.write(reinterpret_cast<char*>(meanSparse.data()),
                  elSize * meanSparse.size());

  varfile.write(reinterpret_cast<char*>(varFlag.data()),
                flagSize * varFlag.size());
  varfile.write(reinterpret_cast<char*>(varSparse.data()),
                elSize * varSparse.size());

  if (useAllBootstraps) {
    boost::iostreams::filtering_ostream& fullBootsfile = *fullBootstrapMatrixStream_;
    if (bootSparse.size() != bootSparse.size()) {
      std::cerr<< "ERROR: all bootstrap dumping error" << std::flush;
      exit(74);
    }

    for (size_t i=0; i<bootSparse.size(); i++) {
      fullBootsfile.write(reinterpret_cast<char*>(bootFlag[i].data()),
                          flagSize * bootFlag[i].size());
      fullBootsfile.write(reinterpret_cast<char*>(bootSparse[i].data()),
                          elSize * bootSparse[i].size());
    }
  }

  namefile.write(bcName.c_str(), bcName.size());
  namefile << std::endl;

  return true;
}

bool GZipWriter::writeSparseAbundances(std::string& bcName,
                                       std::string& features,
                                       std::string& arboData,
                                       uint8_t featureCode,
                                       std::vector<double>& alphas,
                                       std::vector<uint8_t>& tiers,
                                       bool dumpArborescences,
                                       bool dumpUmiGraph){

  // construct the output vectors outside of the critical section
  // since e.g. this is more non-trivial work than in the dense case.
  size_t num = alphas.size();
  std::vector<float> alphasSparse;
  alphasSparse.reserve(num/2);
  std::vector<uint8_t> alphasFlag;
  alphasFlag.reserve(static_cast<size_t>(std::ceil(num/8.0)));

  std::vector<uint8_t> tiersSparse;
  tiersSparse.reserve(num/2);
  std::vector<uint8_t> tiersFlag;
  tiersFlag.reserve(static_cast<size_t>(std::ceil(num/8.0)));

  size_t elSize = sizeof(decltype(alphasSparse)::value_type);
  size_t flagSize = sizeof(decltype(alphasFlag)::value_type);
  size_t trSize = sizeof(decltype(tiersSparse)::value_type);

  for (size_t i=0; i<num; i+=8) {
    uint8_t flag {0};
    for (size_t j=0; j<8; j++) {
      size_t vectorIndex = i+j;
      if (vectorIndex >= num) { break; }

      if (alphas[vectorIndex] > std::numeric_limits<float>::min()) {
        alphasSparse.emplace_back(alphas[vectorIndex]);
        flag |= 128 >> j;
      }
    }
    alphasFlag.emplace_back(flag);
  }

  for (size_t i=0; i<num; i+=8) {
    uint8_t flag {0};
    for (size_t j=0; j<8; j++) {
      size_t vectorIndex = i+j;
      if (vectorIndex >= num) { break; }

      if (tiers[vectorIndex] > std::numeric_limits<uint8_t>::min()) {
        tiersSparse.emplace_back(tiers[vectorIndex]);
        flag |= 128 >> j;
      }
    }
    tiersFlag.emplace_back(flag);
  }

#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  namespace bfs = boost::filesystem;
  if (!countMatrixStream_) {
    countMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    countMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto countMatFilename = path_ / "alevin" / "quants_mat.gz";
    countMatrixStream_->push(boost::iostreams::file_sink(countMatFilename.string(),
                                                         std::ios_base::out | std::ios_base::binary));
    tierMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    tierMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto tierMatFilename = path_ / "alevin" / "quants_tier_mat.gz";
    tierMatrixStream_->push(boost::iostreams::file_sink(tierMatFilename.string(),
                                                        std::ios_base::out | std::ios_base::binary));

    if (dumpArborescences) {
      arboMatrixStream_.reset(new boost::iostreams::filtering_ostream);
      arboMatrixStream_->push(boost::iostreams::gzip_compressor(6));
      auto arboMatFilename = path_ / "alevin" / "arborescence_dump.txt.gz";
      arboMatrixStream_->push(boost::iostreams::file_sink(arboMatFilename.string(),
                                                          std::ios_base::out | std::ios_base::binary));
    }
  }

  if (!bcNameStream_) {
    auto bcNameFilename = path_ / "alevin" / "quants_mat_rows.txt";
    bcNameStream_.reset(new std::ofstream);
    bcNameStream_->open(bcNameFilename.string());

    auto bcFeaturesFilename = path_ / "alevin" / "featureDump.txt";
    bcFeaturesStream_.reset(new std::ofstream);
    bcFeaturesStream_->open(bcFeaturesFilename.string());

    std::string header = "CB\tCorrectedReads\tMappedReads\tDeduplicatedReads"
      "\tMappingRate\tDedupRate\tMeanByMax\tNumGenesExpressed\tNumGenesOverMean";
    if (featureCode == 3) {
      header += "\tmRnaFraction\trRnaFraction";
    } else if (featureCode == 1) {
      header += "\tmRnaFraction";
    } else if (featureCode == 2) {
      header += "\trRnaFraction";
    } else if (featureCode != 0) {
      std::cerr<<"Error: Wrong feature code: " << featureCode << std::flush;
      exit(74);
    } else if (dumpUmiGraph) {
      header += "\tArborescenceCount";
    } header += "\n";
    bcFeaturesStream_->write(header.c_str(), header.size());
  }

  boost::iostreams::filtering_ostream& countfile = *countMatrixStream_;
  boost::iostreams::filtering_ostream& tierfile = *tierMatrixStream_;
  boost::iostreams::filtering_ostream& arbofile = *arboMatrixStream_;

  std::ofstream& namefile = *bcNameStream_;
  std::ofstream& featuresfile = *bcFeaturesStream_;

  countfile.write(reinterpret_cast<char*>(alphasFlag.data()),
                  flagSize * alphasFlag.size());
  countfile.write(reinterpret_cast<char*>(alphasSparse.data()),
                  elSize * alphasSparse.size());

  if (dumpArborescences) {
    arbofile.write(arboData.c_str(), arboData.size());
  }

  tierfile.write(reinterpret_cast<char*>(tiersFlag.data()),
                 flagSize * tiersFlag.size());
  tierfile.write(reinterpret_cast<char*>(tiersSparse.data()),
                 trSize * tiersSparse.size());

  namefile.write(bcName.c_str(), bcName.size());
  namefile << std::endl;

  featuresfile.write(features.c_str(), features.size());
  featuresfile << std::endl;

  return true;
}

bool GZipWriter::writeAbundances(std::string& bcName,
                                 std::string& features,
                                 uint8_t featureCode,
                                 std::vector<double>& alphas,
                                 std::vector<uint8_t>& tiers,
                                 bool dumpUmiGraph){
#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  namespace bfs = boost::filesystem;
  if (!countMatrixStream_) {
    countMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    countMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto countMatFilename = path_ / "alevin" / "quants_mat.gz";
    countMatrixStream_->push(boost::iostreams::file_sink(countMatFilename.string(),
                                                         std::ios_base::out | std::ios_base::binary));

    tierMatrixStream_.reset(new boost::iostreams::filtering_ostream);
    tierMatrixStream_->push(boost::iostreams::gzip_compressor(6));
    auto tierMatFilename = path_ / "alevin" / "quants_tier_mat.gz";
    tierMatrixStream_->push(boost::iostreams::file_sink(tierMatFilename.string(),
                                                        std::ios_base::out | std::ios_base::binary));
  }

  if (!bcNameStream_) {
    auto bcNameFilename = path_ / "alevin" / "quants_mat_rows.txt";
    bcNameStream_.reset(new std::ofstream);
    bcNameStream_->open(bcNameFilename.string());

    auto bcFeaturesFilename = path_ / "alevin" / "featureDump.txt";
    bcFeaturesStream_.reset(new std::ofstream);
    bcFeaturesStream_->open(bcFeaturesFilename.string());

    std::string header = "CB\tCorrectedReads\tMappedReads\tDeduplicatedReads"
      "\tMappingRate\tDedupRate\tMeanByMax\tNumGenesExpressed\tNumGenesOverMean";
    if (featureCode == 3) {
      header += "\tmRnaFraction\trRnaFraction";
    } else if (featureCode == 1) {
      header += "\tmRnaFraction";
    } else if (featureCode == 2) {
      header += "\trRnaFraction";
    } else if (dumpUmiGraph) {
      header += "\tArborescenceCount";
    } header += "\n";
    bcFeaturesStream_->write(header.c_str(), header.size());
  }

  boost::iostreams::filtering_ostream& countfile = *countMatrixStream_;
  boost::iostreams::filtering_ostream& tierfile = *tierMatrixStream_;
  std::ofstream& namefile = *bcNameStream_;
  std::ofstream& featuresfile = *bcFeaturesStream_;

  size_t num = alphas.size();
  size_t elSize = sizeof(typename std::vector<double>::value_type);
  size_t trSize = sizeof(typename std::vector<uint8_t>::value_type);
  countfile.write(reinterpret_cast<char*>(alphas.data()),
                  elSize * num);
  tierfile.write(reinterpret_cast<char*>(tiers.data()),
                  trSize * num);

  namefile.write(bcName.c_str(), bcName.size());
  namefile << std::endl;

  featuresfile.write(features.c_str(), features.size());
  featuresfile << std::endl;

  return true;
}


template <typename ExpT>
bool GZipWriter::writeEmptyAbundances(const SalmonOpts& sopt, ExpT& readExp) {

  namespace bfs = boost::filesystem;

  bfs::path fname = path_ / "quant.sf";
  std::unique_ptr<std::FILE, int (*)(std::FILE*)> output(
      std::fopen(fname.c_str(), "w"), std::fclose);
  auto* outputRaw = output.get();
  fmt::print(outputRaw, "Name\tLength\tEffectiveLength\tTPM\tNumReads\n");
  // Now posterior has the transcript fraction
  std::vector<Transcript>& transcripts_ = readExp.transcripts();
  for (auto& transcript : transcripts_) {
    fmt::print(outputRaw, "{}\t{}\t", transcript.RefName, transcript.CompleteLength);
    fmt::print(outputRaw, "{:.{}f}\t", static_cast<float>(transcript.CompleteLength), sopt.sigDigits);
    fmt::print(outputRaw, "{:f}\t", 0.0);
    fmt::print(outputRaw, "{:.{}f}\n", 0.0, sopt.sigDigits);
  }
  return true;
}

template <typename ExpT>
bool GZipWriter::writeAbundances(const SalmonOpts& sopt, ExpT& readExp, bool explicitSum) {

  namespace bfs = boost::filesystem;

  using salmon::math::LOG_0;
  using salmon::math::LOG_1;

  // If we're using lightweight-alignment (FMD)
  // and not allowing orphans.
  bool useScaledCounts = !(sopt.useQuasi or sopt.allowOrphans or sopt.alnMode);
  bfs::path fname = path_ / "quant.sf";

  std::unique_ptr<std::FILE, int (*)(std::FILE*)> output(
      std::fopen(fname.c_str(), "w"), std::fclose);
  auto* outputRaw = output.get();

  fmt::print(outputRaw, "Name\tLength\tEffectiveLength\tTPM\tNumReads\n");
  std::vector<Transcript>& transcripts_ = readExp.transcripts();

  double numMappedFrags {0.0};
  if ( explicitSum ) {
    for (auto& transcript : transcripts_) {
        numMappedFrags += transcript.sharedCount();
    }
  } else {
    numMappedFrags = readExp.upperBoundHits();
  }

  for (auto& transcript : transcripts_) {
    transcript.projectedCounts = useScaledCounts
                                     ? (transcript.mass(false) * numMappedFrags)
                                     : transcript.sharedCount();
  }

  double tfracDenom{0.0};
  for (auto& transcript : transcripts_) {
    double refLength = transcript.EffectiveLength;
    tfracDenom += (transcript.projectedCounts / numMappedFrags) / refLength;
  }

  double million = 1000000.0;
  // Now posterior has the transcript fraction
  for (auto& transcript : transcripts_) {
    double count = transcript.projectedCounts;
    double npm = (transcript.projectedCounts / numMappedFrags);
    double effLength = transcript.EffectiveLength;
    double tfrac = (npm / effLength) / tfracDenom;
    double tpm = tfrac * million;
    fmt::print(outputRaw, "{}\t{}\t", transcript.RefName, transcript.CompleteLength);
    fmt::print(outputRaw, "{:.{}f}\t", effLength, sopt.sigDigits);
    fmt::print(outputRaw, "{:f}\t", tpm);
    fmt::print(outputRaw, "{:.{}f}\n", count, sopt.sigDigits);
  }
  return true;
}

bool GZipWriter::setSamplingPath(const SalmonOpts& sopt) {
  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / sopt.auxDir;
  if (!bfs::exists(auxDir)) {
    bool auxSuccess = boost::filesystem::create_directories(auxDir);
    if (!auxSuccess) {
      sopt.jointLog->critical("Could not create auxiliary directory {}",
                              auxDir.string());
      return false;
    }
  }
  bsPath_ = auxDir / "bootstrap";
  if (!bfs::exists(bsPath_)) {
    bool bsSuccess = boost::filesystem::create_directories(bsPath_);
    if (!bsSuccess) {
      sopt.jointLog->critical("Could not create sampling directory {}",
                              bsPath_.string());
      return false;
    }
  }
  return true;
}

template <typename T>
bool GZipWriter::writeBootstrap(const std::vector<T>& abund, bool quiet) {
#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  if (!bsStream_) {
    bsStream_.reset(new boost::iostreams::filtering_ostream);
    bsStream_->push(boost::iostreams::gzip_compressor(6));
    auto bsFilename = bsPath_ / "bootstraps.gz";
    bsStream_->push(boost::iostreams::file_sink(
        bsFilename.string(), std::ios_base::out | std::ios_base::binary));
  }

  boost::iostreams::filtering_ostream& ofile = *bsStream_;
  size_t num = abund.size();
  size_t elSize = sizeof(typename std::vector<T>::value_type);
  ofile.write(reinterpret_cast<char*>(const_cast<T*>(abund.data())),
              elSize * num);
  if (!quiet) {
    logger_->info("wrote {} bootstraps", numBootstrapsWritten_.load() + 1);
  }
  ++numBootstrapsWritten_;
  return true;
}

bool GZipWriter::writeCellEQVec(size_t barcode, const std::vector<uint32_t>& offsets,
                                const std::vector<uint32_t>& counts, bool quiet) {
#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  if (!cellEQStream_) {
    cellEQStream_.reset(new boost::iostreams::filtering_ostream);
    cellEQStream_->push(boost::iostreams::gzip_compressor(6));
    auto ceqFilename = path_ / "alevin" / "cell_eq_mat.gz";
    cellEQStream_->push(boost::iostreams::file_sink(ceqFilename.string(),
                                                    std::ios_base::out | std::ios_base::binary));
  }

  boost::iostreams::filtering_ostream& ofile = *cellEQStream_;
  size_t num = offsets.size();
  size_t elSize = sizeof(typename std::vector<uint32_t>::value_type);
  // write the barcode
  ofile.write(reinterpret_cast<char*>(&barcode), sizeof(barcode));
  // write the number of elements in the list
  ofile.write(reinterpret_cast<char*>(&num), sizeof(barcode));
  // write the offsets and counts
  ofile.write(reinterpret_cast<char*>(const_cast<uint32_t*>(offsets.data())), elSize * num);
  ofile.write(reinterpret_cast<char*>(const_cast<uint32_t*>(counts.data())), elSize * num);

  if (!quiet) {
    logger_->info("wrote EQ vector for barcode ID {}", barcode);
  }
  return true;
}

// TODO: complete this method
bool GZipWriter::writeDedupCellEQVec(
  size_t barcode, 
  const std::vector<std::vector<uint32_t>>& labels,
  const std::vector<uint32_t>& counts,
  bool quiet
  ) {
#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  if (!cellDedupEQStream_) {
    cellDedupEQStream_.reset(new boost::iostreams::filtering_ostream);
    cellDedupEQStream_->push(boost::iostreams::gzip_compressor(6));
    auto ceqFilename = path_ / "alevin" / "cell_dedup_eq_mat.gz";
    cellDedupEQStream_->push(boost::iostreams::file_sink(ceqFilename.string(),
                                                    std::ios_base::out | std::ios_base::binary));
  }
  
  boost::iostreams::filtering_ostream& ofile = *cellDedupEQStream_;
  size_t num = labels.size();
  size_t elSize = sizeof(typename std::vector<uint32_t>::value_type);
  // write the barcode
  ofile.write(reinterpret_cast<char*>(&barcode), sizeof(barcode));
  // write the number of elements in the list
  ofile.write(reinterpret_cast<char*>(&num), sizeof(barcode));
  // write the individual labels
  for (size_t i = 0 ; i < num ; ++i) {
    size_t labelLength = labels[i].size();
    // write the length of the equivalence class vector
    ofile.write(reinterpret_cast<char*>(&labelLength), sizeof(barcode));
    // write the label vector
    ofile.write(reinterpret_cast<char*>(const_cast<uint32_t*>(labels[i].data())), elSize * labelLength);
  }
  ofile.write(reinterpret_cast<char*>(const_cast<uint32_t*>(counts.data())), elSize * num);
  // write the offsets and counts
  if (!quiet) {
    logger_->info("wrote EQ vector for barcode ID {}", barcode);
  }
  return true;
}

bool GZipWriter::writeUmiGraph(alevin::graph::Graph& g,
                               std::string& cbString) {
#if defined __APPLE__
  spin_lock::scoped_lock sl(writeMutex_);
#else
  std::lock_guard<std::mutex> lock(writeMutex_);
#endif
  if (!umiGraphStream_) {
    umiGraphStream_.reset(new boost::iostreams::filtering_ostream);
    umiGraphStream_->push(boost::iostreams::gzip_compressor(6));
    auto graphFilename = path_ / "alevin" / "cell_umi_graphs.gz";
    umiGraphStream_->push(boost::iostreams::file_sink(graphFilename.string(),
                                                      std::ios_base::out));
  }

  boost::iostreams::filtering_ostream& ofile = *umiGraphStream_;
  ofile << cbString << "\t" << g.num_vertices() << "\t";
  for(auto& edge: g.edges){
    uint32_t source = edge.first;
    ofile<< source;
    for(uint32_t sink: edge.second) {
      ofile<< "," << sink;
    }
    ofile<< "\t";
  }

  ofile << std::endl;
  return true;
}

void GZipWriter::writeMtx(std::shared_ptr<spdlog::logger>& jointLog, 
              boost::filesystem::path& outputDirectory, 
              size_t numGenes, size_t numCells, size_t totalExpGeneCounts) {
    jointLog->info("Starting dumping cell v gene counts in mtx format");
    boost::filesystem::path qFilePath = outputDirectory / "quants_mat.mtx.gz";

    boost::iostreams::filtering_ostream qFile;
    qFile.push(boost::iostreams::gzip_compressor(6));
    qFile.push(boost::iostreams::file_sink(qFilePath.string(),
                                           std::ios_base::out | std::ios_base::binary));

    // mtx header
    qFile << "%%MatrixMarket\tmatrix\tcoordinate\treal\tgeneral" << std::endl
          << numCells << "\t" << numGenes << "\t" << totalExpGeneCounts << std::endl;

    {

      auto popcount = [](uint8_t n) {
        size_t count {0};
        while (n) {
          n &= n-1;
          ++count;
        }
        return count;
      };

      uint32_t zerod_cells {0};
      size_t numFlags = std::ceil(numGenes/8.0);
      std::vector<uint8_t> alphasFlag (numFlags, 0);
      size_t flagSize = sizeof(decltype(alphasFlag)::value_type);

      std::vector<float> alphasSparse;
      alphasSparse.reserve(numFlags/2);
      size_t elSize = sizeof(decltype(alphasSparse)::value_type);

      auto countMatFilename = outputDirectory / "quants_mat.gz";
      if(not boost::filesystem::exists(countMatFilename)){
        jointLog->error("Can't import Binary file quants.mat.gz, it doesn't exist");
        jointLog->flush();
        exit(84);
      }

      boost::iostreams::filtering_istream countMatrixStream;
      countMatrixStream.push(boost::iostreams::gzip_decompressor());
      countMatrixStream.push(boost::iostreams::file_source(countMatFilename.string(),
                                                           std::ios_base::in | std::ios_base::binary));

      for (size_t cellCount=0; cellCount<numCells; cellCount++){
        countMatrixStream.read(reinterpret_cast<char*>(alphasFlag.data()), flagSize * numFlags);

        size_t numExpGenes {0};
        std::vector<size_t> indices;
        for (size_t j=0; j<alphasFlag.size(); j++) {
          uint8_t flag = alphasFlag[j];
          size_t numNonZeros = popcount(flag);
          numExpGenes += numNonZeros;

          for (size_t i=0; i<8; i++){
            if (flag & (128 >> i)) {
              indices.emplace_back( i+(8*j) );
            }
          }
        }

        if (indices.size() != numExpGenes) {
          jointLog->error("binary format reading error {}: {}: {}",
                           indices.size(), numExpGenes);
          jointLog->flush();
          exit(84);
        }


        alphasSparse.clear();
        alphasSparse.resize(numExpGenes);
        countMatrixStream.read(reinterpret_cast<char*>(alphasSparse.data()), elSize * numExpGenes);

        float readCount {0.0};
        readCount += std::accumulate(alphasSparse.begin(), alphasSparse.end(), 0.0);
        if (readCount > 1000000) {
          jointLog->warn("A cell has more 1M count, Possible error");
          jointLog->flush();
        }

        for(size_t i=0; i<numExpGenes; i++) {
          qFile << std::fixed
                << cellCount + 1 << "\t"
                << indices[i] + 1 << "\t"
                << alphasSparse[i] <<  std::endl;
        }

        if (readCount == 0.0){
          zerod_cells += 1;
        }
      } // end-for each cell

      if (zerod_cells > 0) {
        jointLog->warn("Found {} cells with 0 counts", zerod_cells);
      }
    }

    boost::iostreams::close(qFile);
    jointLog->info("Finished dumping counts into mtx");
}

using SCExpT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;
using BulkExpT = ReadExperiment<EquivalenceClassBuilder<TGValue>>;
template <typename FragT>
using BulkAlignLibT = AlignmentLibrary<FragT, EquivalenceClassBuilder<TGValue>>;

template bool
GZipWriter::writeBootstrap<double>(const std::vector<double>& abund,
                                   bool quiet);

template bool GZipWriter::writeBootstrap<int>(const std::vector<int>& abund,
                                              bool quiet);

template bool GZipWriter::writeEquivCounts<BulkExpT>(const SalmonOpts& sopt,
                                             BulkExpT& readExp);

template bool GZipWriter::writeEquivCounts<BulkAlignLibT<UnpairedRead>>(
    const SalmonOpts& sopt, BulkAlignLibT<UnpairedRead>& readExp);

template bool GZipWriter::writeEquivCounts<BulkAlignLibT<ReadPair>>(
    const SalmonOpts& sopt, BulkAlignLibT<ReadPair>& readExp);

template bool GZipWriter::writeBFH<SCExpT>(boost::filesystem::path& outDir,
                                           SCExpT& experiment, size_t umiLength,
                                           std::vector<std::string>& bcSeqVec);

template bool GZipWriter::writeAbundances<BulkExpT>(const SalmonOpts& sopt,
                                                    BulkExpT& readExp,
                                                    bool explicitSum);

template bool GZipWriter::writeAbundances<BulkAlignLibT<UnpairedRead>>(const SalmonOpts& sopt,
                                                                       BulkAlignLibT<UnpairedRead>& readExp,
                                                                       bool explicitSum);
template bool GZipWriter::writeAbundances<BulkAlignLibT<ReadPair>>(const SalmonOpts& sopt,
                                                                   BulkAlignLibT<ReadPair>& readExp,
                                                                   bool explicitSum);

template bool GZipWriter::writeEmptyAbundances<BulkExpT>(const SalmonOpts& sopt,
                                                 BulkExpT& readExp);
template bool GZipWriter::writeEmptyAbundances<SCExpT>(const SalmonOpts& sopt,
                                                               SCExpT& readExp);

template bool GZipWriter::writeEmptyAbundances<BulkAlignLibT<UnpairedRead>>(
    const SalmonOpts& sopt, BulkAlignLibT<UnpairedRead>& readExp);
template bool GZipWriter::writeEmptyAbundances<BulkAlignLibT<ReadPair>>(
    const SalmonOpts& sopt, BulkAlignLibT<ReadPair>& readExp);

template bool
GZipWriter::writeMeta<BulkExpT>(const SalmonOpts& opts,
                                const BulkExpT& experiment,
                                const MappingStatistics& mstats);
template bool
GZipWriter::writeMeta<SCExpT>(const SalmonOpts& opts,
                              const SCExpT& experiment,
                              const MappingStatistics& mstats);

template bool
GZipWriter::writeMeta<BulkAlignLibT<UnpairedRead>>(const SalmonOpts& opts, const BulkAlignLibT<UnpairedRead>& experiment,
                                                   const MappingStatistics& mstats);

template bool
GZipWriter::writeMeta<BulkAlignLibT<ReadPair>>(const SalmonOpts& opts, const BulkAlignLibT<ReadPair>& experiment,
                                               const MappingStatistics& mstats);

template bool
GZipWriter::writeEmptyMeta<BulkExpT>(const SalmonOpts& opts,
                                           const BulkExpT& experiment,
                                           std::vector<std::string>& errors);
template bool
GZipWriter::writeEmptyMeta<SCExpT>(const SalmonOpts& opts,
                                           const SCExpT& experiment,
                                           std::vector<std::string>& errors);

template bool GZipWriter::writeEmptyMeta<BulkAlignLibT<UnpairedRead>>(
    const SalmonOpts& opts, const BulkAlignLibT<UnpairedRead>& experiment,
    std::vector<std::string>& errors);

template bool GZipWriter::writeEmptyMeta<BulkAlignLibT<ReadPair>>(
    const SalmonOpts& opts, const BulkAlignLibT<ReadPair>& experiment,
    std::vector<std::string>& errors);

template std::vector<std::string>
getLibTypeStrings(const BulkAlignLibT<UnpairedRead>& experiment);

template std::vector<std::string>
getLibTypeStrings(const BulkAlignLibT<ReadPair>& experiment);

namespace apt = alevin::protocols;

template
bool GZipWriter::writeEquivCounts<SCExpT, apt::DropSeq>(
                                                        const AlevinOpts<apt::DropSeq>& aopts,
                                                        SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::CITESeq>(
                                                        const AlevinOpts<apt::CITESeq>& aopts,
                                                        SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::InDrop>(
                                                       const AlevinOpts<apt::InDrop>& aopts,
                                                       SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::ChromiumV3>(
                                                         const AlevinOpts<apt::ChromiumV3>& aopts,
                                                         SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::Chromium>(
                                                         const AlevinOpts<apt::Chromium>& aopts,
                                                         SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::Gemcode>(
                                                        const AlevinOpts<apt::Gemcode>& aopts,
                                                        SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::CELSeq>(
                                                        const AlevinOpts<apt::CELSeq>& aopts,
                                                        SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::CELSeq2>(
                                                        const AlevinOpts<apt::CELSeq2>& aopts,
                                                        SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::QuartzSeq2>(
                                                        const AlevinOpts<apt::QuartzSeq2>& aopts,
                                                        SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::Custom>(
                                                       const AlevinOpts<apt::Custom>& aopts,
                                                       SCExpT& readExp);
template
bool GZipWriter::writeEquivCounts<SCExpT, apt::CustomGeometry>(
                                                       const AlevinOpts<apt::CustomGeometry>& aopts,
                                                       SCExpT& readExp);
template bool
GZipWriter::writeMetaAlevin<apt::DropSeq>(const AlevinOpts<apt::DropSeq>& opts,
                                          boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::CITESeq>(const AlevinOpts<apt::CITESeq>& opts,
                                          boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::InDrop>(const AlevinOpts<apt::InDrop>& opts,
                                         boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::Chromium>(const AlevinOpts<apt::Chromium>& opts,
                                           boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::ChromiumV3>(const AlevinOpts<apt::ChromiumV3>& opts,
                                             boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::CELSeq>(const AlevinOpts<apt::CELSeq>& opts,
                                         boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::CELSeq2>(const AlevinOpts<apt::CELSeq2>& opts,
                                          boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::QuartzSeq2>(const AlevinOpts<apt::QuartzSeq2>& opts,
                                             boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::Custom>(const AlevinOpts<apt::Custom>& opts,
                                         boost::filesystem::path aux_dir);
template bool
GZipWriter::writeMetaAlevin<apt::Gemcode>(const AlevinOpts<apt::Gemcode>& opts,
                                          boost::filesystem::path aux_dir);

template bool
GZipWriter::writeMetaAlevin<apt::CustomGeometry>(const AlevinOpts<apt::CustomGeometry>& opts,
                                          boost::filesystem::path aux_dir);


