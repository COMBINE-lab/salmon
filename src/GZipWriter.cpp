#include <ctime>
#include <fstream>

#include "cereal/archives/json.hpp"

#include "DistributionUtils.hpp"
#include "GZipWriter.hpp"
#include "SalmonOpts.hpp"
#include "ReadExperiment.hpp"
#include "AlignmentLibrary.hpp"
#include "ReadPair.hpp"
#include "UnpairedRead.hpp"

GZipWriter::GZipWriter(const boost::filesystem::path path, std::shared_ptr<spdlog::logger> logger) :
  path_(path), logger_(logger) {
}

GZipWriter::~GZipWriter() {
  if (bsStream_) {
    bsStream_->reset();
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
bool GZipWriter::writeEquivCounts(
    const SalmonOpts& opts,
    ExpT& experiment) {

  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / opts.auxDir;
  bool auxSuccess = boost::filesystem::create_directories(auxDir);
  bfs::path eqFilePath = auxDir / "eq_classes.txt";

  std::ofstream equivFile(eqFilePath.string());

  auto& transcripts = experiment.transcripts();
  std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        experiment.equivalenceClassBuilder().eqVec();
  bool dumpRichWeights = opts.dumpEqWeights;

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
    equivFile << txps.size() << '\t';
    // each group member
    for (auto tid : txps) { equivFile << tid << '\t'; }
    if (dumpRichWeights) {
      const auto& auxs = eq.second.combinedWeights;
      for (auto aux : auxs) { equivFile << aux << '\t'; }
    }
    // count for this class
    equivFile << count << '\n';
  }

  equivFile.close();
  return true;
}

std::vector<std::string> getLibTypeStrings(const ReadExperiment& experiment) {
  auto& libs = experiment.readLibraries();
  std::vector<std::string> libStrings;
  for (auto& rl : libs) {
    libStrings.push_back(rl.getFormat().toString());
  }
  return libStrings;
}

template <typename AlnT>
std::vector<std::string> getLibTypeStrings(const AlignmentLibrary<AlnT>& experiment) {
  std::vector<std::string> libStrings;
  libStrings.push_back(experiment.format().toString());
  return libStrings;
}

/**
 * Write the ``main'' metadata to file.  Currently this includes:
 *   -- Names of the target id's if bootstrapping / gibbs is performed
 *   -- The fragment length distribution
 *   -- The expected and observed bias values
 *   -- A json file with information about the run
 */
template <typename ExpT>
bool GZipWriter::writeMeta(
    const SalmonOpts& opts,
    const ExpT& experiment) {

  namespace bfs = boost::filesystem;

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
          nameOut.push(boost::iostreams::file_sink(bsFilename.string(), std::ios_base::out));

          auto& transcripts = experiment.transcripts();
          size_t numTxps = transcripts.size();
          if (numTxps == 0) { return false; }
          for (size_t tn = 0; tn < numTxps; ++tn) {
              auto& t  = transcripts[tn];
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
  auto fldSamples = distribution_utils::samplesFromLogPMF(
                        experiment.fragmentLengthDistribution(), numFLDSamples);
  writeVectorToFile(fldPath, fldSamples);

  bfs::path normBiasPath = auxDir / "expected_bias.gz";
  writeVectorToFile(normBiasPath, experiment.expectedSeqBias());

  bfs::path obsBiasPath = auxDir / "observed_bias.gz";
  // TODO: dump both sense and anti-sense models
  const auto& bcounts = experiment.readBias(salmon::utils::Direction::FORWARD).counts;
  std::vector<int32_t> observedBias(bcounts.size(), 0);
  std::copy(bcounts.begin(), bcounts.end(), observedBias.begin());
  writeVectorToFile(obsBiasPath, observedBias);

  bfs::path obsBiasPath3p = auxDir / "observed_bias_3p.gz";
  const auto& bcounts3p = experiment.readBias(salmon::utils::Direction::REVERSE_COMPLEMENT).counts;
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
      auto& obs5 = experiment.readBiasModelObserved(salmon::utils::Direction::FORWARD);
      obs5.writeBinary(out);
    }
    // 3' observed
    {
      bfs::path obs3Path = auxDir / "obs3_seq.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(obs3Path.string(), flags));
      auto& obs3 = experiment.readBiasModelObserved(salmon::utils::Direction::REVERSE_COMPLEMENT);
      obs3.writeBinary(out);
    }

    // 5' expected
    {
      bfs::path exp5Path = auxDir / "exp5_seq.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(exp5Path.string(), flags));
      auto& exp5 = experiment.readBiasModelExpected(salmon::utils::Direction::FORWARD);
      exp5.writeBinary(out);
    }
    // 3' expected
    {
      bfs::path exp3Path = auxDir / "exp3_seq.gz";
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(exp3Path.string(), flags));
      auto& exp3 = experiment.readBiasModelExpected(salmon::utils::Direction::REVERSE_COMPLEMENT);
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
    
    // lambda to write out a vector of SimplePosBias models (along with the length bounds) to file.
    auto writePosModel= [&lenBounds, this](bfs::path fpath, const std::vector<SimplePosBias>& model) -> bool {
      auto flags = std::ios_base::out | std::ios_base::binary;
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor(6));
      out.push(boost::iostreams::file_sink(fpath.string(), flags));
      // Write out the number of different models
      uint32_t numModels = static_cast<uint32_t>(lenBounds.size());
      out.write(reinterpret_cast<char*>(&numModels), sizeof(numModels));
      // Write out the length class for each model
      for (const auto& b : lenBounds) {
        out.write(reinterpret_cast<char*>(const_cast<uint32_t*>(&b)), sizeof(b));
      }
      // write out each
      for (auto& pb : model) { 
	bool success = pb.writeBinary(out);
	if (!success) {
	  this->logger_->error("Could not write out positional bias model to {}!", fpath.string());
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
    //3' observed
    {
      bfs::path obsPosPath = auxDir / "obs3_pos.gz";
      // Get the pos bias vector
      auto& posBiases = experiment.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
      writePosModel(obsPosPath, posBiases);
    }
    // 5' expected
    {
      bfs::path expPosPath = auxDir / "exp5_pos.gz";
      // Get the pos bias vector
      auto& posBiases = experiment.posBiasExpected(salmon::utils::Direction::FORWARD);
      writePosModel(expPosPath, posBiases);
    }
    // 3' expected
    {
      bfs::path expPosPath = auxDir / "exp3_pos.gz";
      // Get the pos bias vector
      auto& posBiases = experiment.posBiasExpected(salmon::utils::Direction::REVERSE_COMPLEMENT);
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

      auto libStrings = getLibTypeStrings(experiment);
      oa(cereal::make_nvp("num_libraries", libStrings.size()));
      oa(cereal::make_nvp("library_types", libStrings));

      oa(cereal::make_nvp("frag_dist_length", fldSamples.size()));
      oa(cereal::make_nvp("seq_bias_correct", opts.biasCorrect));
      oa(cereal::make_nvp("gc_bias_correct", opts.gcBiasCorrect));
      oa(cereal::make_nvp("num_bias_bins", bcounts.size()));

      std::string mapTypeStr = opts.alnMode ? "alignment" : "mapping";
      oa(cereal::make_nvp("mapping_type", mapTypeStr));

      oa(cereal::make_nvp("num_targets", transcripts.size()));

      // True if we dumped the equivalence classes, false otherwise
      oa(cereal::make_nvp("serialized_eq_classes", opts.dumpEq));
      
      // For now, this vector is empty unless we dumped the equivalence classes
      // with weights.  In which case it contains the string "scalar_weights".
      std::vector<std::string> props;
      if (opts.dumpEqWeights) { props.push_back("scalar_weights"); }
      oa(cereal::make_nvp("eq_class_properties", props));

      oa(cereal::make_nvp("length_classes", experiment.getLengthQuantiles()));
      oa(cereal::make_nvp("index_seq_hash", experiment.getIndexSeqHash()));
      oa(cereal::make_nvp("index_name_hash", experiment.getIndexNameHash()));
      oa(cereal::make_nvp("num_bootstraps", numSamples));
      oa(cereal::make_nvp("num_processed", experiment.numObservedFragments()));
      oa(cereal::make_nvp("num_mapped", experiment.numMappedFragments()));
      oa(cereal::make_nvp("percent_mapped", experiment.effectiveMappingRate() * 100.0));
      oa(cereal::make_nvp("call", std::string("quant")));
      oa(cereal::make_nvp("start_time", opts.runStartTime));
      oa(cereal::make_nvp("end_time", opts.runStopTime));
  }

  {
    bfs::path ambigInfo = auxDir / "ambig_info.tsv";
    std::ofstream os(ambigInfo.string());
    os << "UniqueCount\tAmbigCount\n";

    auto& transcripts = experiment.transcripts();
    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
      const_cast<ExpT&>(experiment).equivalenceClassBuilder().eqVec();

    class CountPair {
    public:
      uint32_t unique=0;
      uint32_t potential=0;
    };

    std::vector<CountPair> counts(transcripts.size());
    for (auto& eq : eqVec) {
      uint64_t count = eq.second.count;
      const TranscriptGroup& tgroup = eq.first;
      const std::vector<uint32_t>& txps = tgroup.txps;
      if (txps.size() > 1) {
        for (auto tid : txps) {
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

template <typename ExpT>
bool GZipWriter::writeAbundances(
    const SalmonOpts& sopt,
    ExpT& readExp) {

  namespace bfs = boost::filesystem;

  using salmon::math::LOG_0;
  using salmon::math::LOG_1;

  // If we're using lightweight-alignment (FMD)
  // and not allowing orphans.
  bool useScaledCounts = (!sopt.useQuasi and sopt.allowOrphans == false);
  bfs::path fname = path_ / "quant.sf";

  std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(fname.c_str(), "w"), std::fclose);

  fmt::print(output.get(), "Name\tLength\tEffectiveLength\tTPM\tNumReads\n");

  double numMappedFrags = readExp.upperBoundHits();

  std::vector<Transcript>& transcripts_ = readExp.transcripts();
  for (auto& transcript : transcripts_) {
      transcript.projectedCounts = useScaledCounts ?
          (transcript.mass(false) * numMappedFrags) : transcript.sharedCount();
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
      fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n",
              transcript.RefName, transcript.CompleteLength, effLength,
              tpm, count);
  }
  return true;
}

bool GZipWriter::setSamplingPath(const SalmonOpts& sopt) {
  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / sopt.auxDir;
  if (!bfs::exists(auxDir)) {
    bool auxSuccess = boost::filesystem::create_directories(auxDir);
    if (!auxSuccess) {
      sopt.jointLog->critical("Could not create auxiliary directory {}", auxDir.string());
      return false;
    }
  }
  bsPath_ = auxDir / "bootstrap";
  if (!bfs::exists(bsPath_)) {
    bool bsSuccess = boost::filesystem::create_directories(bsPath_);
    if (!bsSuccess) {
      sopt.jointLog->critical("Could not create sampling directory {}", bsPath_.string());
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
	      bsStream_->push(
                  boost::iostreams::file_sink(bsFilename.string(),
                                              std::ios_base::out | std::ios_base::binary));
	    }

	    boost::iostreams::filtering_ostream& ofile = *bsStream_;
	    size_t num = abund.size();
        size_t elSize = sizeof(typename std::vector<T>::value_type);
        ofile.write(reinterpret_cast<char*>(const_cast<T*>(abund.data())),
                    elSize * num);
        if (!quiet){
          logger_->info("wrote {} bootstraps", numBootstrapsWritten_.load()+1);
        }
        ++numBootstrapsWritten_;
        return true;
}

template
bool GZipWriter::writeBootstrap<double>(const std::vector<double>& abund, bool quiet);

template
bool GZipWriter::writeBootstrap<int>(const std::vector<int>& abund, bool quiet);

template
bool GZipWriter::writeEquivCounts<ReadExperiment>(const SalmonOpts& sopt,
                                                 ReadExperiment& readExp);
template
bool GZipWriter::writeEquivCounts<AlignmentLibrary<UnpairedRead>>(const SalmonOpts& sopt,
                                                 AlignmentLibrary<UnpairedRead>& readExp);
template
bool GZipWriter::writeEquivCounts<AlignmentLibrary<ReadPair>>(const SalmonOpts& sopt,
                                                 AlignmentLibrary<ReadPair>& readExp);
template
bool GZipWriter::writeAbundances<ReadExperiment>(const SalmonOpts& sopt,
                                                 ReadExperiment& readExp);
template
bool GZipWriter::writeAbundances<AlignmentLibrary<UnpairedRead>>(const SalmonOpts& sopt,
                                                 AlignmentLibrary<UnpairedRead>& readExp);
template
bool GZipWriter::writeAbundances<AlignmentLibrary<ReadPair>>(const SalmonOpts& sopt,
                                                 AlignmentLibrary<ReadPair>& readExp);

template
bool GZipWriter::writeMeta<ReadExperiment>(
    const SalmonOpts& opts,
    const ReadExperiment& experiment);

template
bool GZipWriter::writeMeta<AlignmentLibrary<UnpairedRead>>(
    const SalmonOpts& opts,
    const AlignmentLibrary<UnpairedRead>& experiment);

template
bool GZipWriter::writeMeta<AlignmentLibrary<ReadPair>>(
    const SalmonOpts& opts,
    const AlignmentLibrary<ReadPair>& experiment);

template
std::vector<std::string> getLibTypeStrings(const AlignmentLibrary<UnpairedRead>& experiment);

template
std::vector<std::string> getLibTypeStrings(const AlignmentLibrary<ReadPair>& experiment);
