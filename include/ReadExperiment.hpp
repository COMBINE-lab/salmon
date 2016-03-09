#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
}

// Our includes
#include "ClusterForest.hpp"
#include "Transcript.hpp"
#include "ReadLibrary.hpp"
#include "FragmentLengthDistribution.hpp"
#include "FragmentStartPositionDistribution.hpp"
#include "SequenceBiasModel.hpp"
#include "SalmonOpts.hpp"
#include "SalmonIndex.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "SpinLock.hpp" // RapMap's with try_lock
#include "UtilityFunctions.hpp"
#include "ReadKmerDist.hpp"

// Logger includes
#include "spdlog/spdlog.h"

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

// Standard includes
#include <vector>
#include <memory>
#include <fstream>


/**
  *  This class represents a library of alignments used to quantify
  *  a set of target transcripts.  The AlignmentLibrary contains info
  *  about both the alignment file and the target sequence (transcripts).
  *  It is used to group them together and track information about them
  *  during the quantification procedure.
  */
class ReadExperiment {

    public:

    ReadExperiment(std::vector<ReadLibrary>& readLibraries,
                   //const boost::filesystem::path& transcriptFile,
                   const boost::filesystem::path& indexDirectory,
		           SalmonOpts& sopt) :
        readLibraries_(readLibraries),
        //transcriptFile_(transcriptFile),
        transcripts_(std::vector<Transcript>()),
        totalAssignedFragments_(0),
        fragStartDists_(5),
        seqBiasModel_(1.0),
	eqBuilder_(sopt.jointLog),
        expectedBias_(constExprPow(4, readBias_.getK()), 1.0),
        expectedGC_(101, 0.0),
        observedGC_(101, 1e-5) {
            namespace bfs = boost::filesystem;

            // Make sure the read libraries are valid.
            for (auto& rl : readLibraries_) { rl.checkValid(); }

            size_t maxFragLen = sopt.fragLenDistMax;
            size_t meanFragLen = sopt.fragLenDistPriorMean;
            size_t fragLenStd = sopt.fragLenDistPriorSD;
            size_t fragLenKernelN = 4;
            double fragLenKernelP = 0.5;
            fragLengthDist_.reset(new FragmentLengthDistribution(1.0, maxFragLen,
                    meanFragLen, fragLenStd,
                    fragLenKernelN,
                    fragLenKernelP, 1));

            // Make sure the transcript file exists.
            /*
            if (!bfs::exists(transcriptFile_)) {
                std::stringstream ss;
                ss << "The provided transcript file: " << transcriptFile_ <<
                    " does not exist!\n";
                throw std::invalid_argument(ss.str());
            }
            */

            // ==== Figure out the index type
            boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
            SalmonIndexVersionInfo versionInfo;
            versionInfo.load(versionPath);
            if (versionInfo.indexVersion() == 0) {
                fmt::MemoryWriter infostr;
                infostr << "Error: The index version file " << versionPath.string()
                    << " doesn't seem to exist.  Please try re-building the salmon "
                    "index.";
                throw std::invalid_argument(infostr.str());
            }
            // Check index version compatibility here
            auto indexType = versionInfo.indexType();
            // ==== Figure out the index type

            salmonIndex_.reset(new SalmonIndex(sopt.jointLog, indexType));
            salmonIndex_->load(indexDirectory);

	    // Now we'll have either an FMD-based index or a QUASI index
	    // dispatch on the correct type.

	    switch (salmonIndex_->indexType()) {
            case SalmonIndexType::QUASI:
                if (salmonIndex_->is64BitQuasi()) {
                    if (salmonIndex_->isPerfectHashQuasi()) {
                        loadTranscriptsFromQuasi(salmonIndex_->quasiIndexPerfectHash64(), sopt);
                    } else {
                        loadTranscriptsFromQuasi(salmonIndex_->quasiIndex64(), sopt);
                    }
                } else {
                    if (salmonIndex_->isPerfectHashQuasi()) {
                        loadTranscriptsFromQuasi(salmonIndex_->quasiIndexPerfectHash32(), sopt);
                    } else {
                        loadTranscriptsFromQuasi(salmonIndex_->quasiIndex32(), sopt);
                    }
                }
                break;
            case SalmonIndexType::FMD:
                loadTranscriptsFromFMD();
                break;
	    }


            // Create the cluster forest for this set of transcripts
            clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));
        }

    EquivalenceClassBuilder& equivalenceClassBuilder() {
        return eqBuilder_;
    }

    std::vector<Transcript>& transcripts() { return transcripts_; }
    const std::vector<Transcript>& transcripts() const { return transcripts_; }

    void updateTranscriptLengthsAtomic(std::atomic<bool>& done) {
        if (sl_.try_lock()) {
            if (!done) {
                auto& fld = *(fragLengthDist_.get());

                std::vector<double> logPMF;
                size_t minVal;
                size_t maxVal;
                double logFLDMean = fld.mean();
                fld.dumpPMF(logPMF, minVal, maxVal);
                double sum = salmon::math::LOG_0;
                for (auto v : logPMF) {
                    sum = salmon::math::logAdd(sum, v);
                }
                for (auto& v : logPMF) {
                    v -= sum;
                }
                // Update the effective length of *every* transcript
                for( auto& t : transcripts_ ) {
                    t.updateEffectiveLength(logPMF, logFLDMean, minVal, maxVal);
                }
                // then declare that we are done
                done = true;
                sl_.unlock();
            }
        }
    }

    uint64_t numAssignedFragments() { return numAssignedFragments_; }
    uint64_t numMappedFragments() const { return numAssignedFragments_; }

    uint64_t upperBoundHits() { return upperBoundHits_; }
    void setUpperBoundHits(uint64_t ubh) { upperBoundHits_ = ubh; }

    std::atomic<uint64_t>& numAssignedFragmentsAtomic() { return numAssignedFragments_; }

    void setNumObservedFragments(uint64_t numObserved) { numObservedFragments_ = numObserved; }

    uint64_t numObservedFragments() const {
        return numObservedFragments_;
    }

    double mappingRate() {
        if (quantificationPasses_ > 0) {
            return static_cast<double>(numAssignedFragsInFirstPass_) / numObservedFragsInFirstPass_;
        } else {
            return static_cast<double>(numAssignedFragments_) / numObservedFragments_;
        }
    }

    SalmonIndex* getIndex() { return salmonIndex_.get(); }

    template <typename QuasiIndexT>
    void loadTranscriptsFromQuasi(QuasiIndexT* idx_, const SalmonOpts& sopt) {
	    size_t numRecords = idx_->txpNames.size();

	    fmt::print(stderr, "Index contained {} targets\n", numRecords);
	    //transcripts_.resize(numRecords);
	    double alpha = 0.005;
	    for (auto i : boost::irange(size_t(0), numRecords)) {
		    uint32_t id = i;
		    const char* name = idx_->txpNames[i].c_str();
		    uint32_t len = idx_->txpLens[i];
		    // copy over the length, then we're done.
		    transcripts_.emplace_back(id, name, len, alpha);
		    auto& txp = transcripts_.back();
		    // The transcript sequence
		    //auto txpSeq = idx_->seq.substr(idx_->txpOffsets[i], len);

		    // Set the transcript sequence
		    txp.setSequenceBorrowed(idx_->seq.c_str() + idx_->txpOffsets[i],
                                    sopt.gcBiasCorrect, sopt.gcSampFactor);
		    // Length classes taken from
		    // ======
		    // Roberts, Adam, et al.
		    // "Improving RNA-Seq expression estimates by correcting for fragment bias."
		    // Genome Biol 12.3 (2011): R22.
		    // ======
		    // perhaps, define these in a more data-driven way
        if (txp.RefLength <= 1334) {
          txp.lengthClassIndex(0);
        } else if (txp.RefLength <= 2104) {
          txp.lengthClassIndex(0);
        } else if (txp.RefLength <= 2988) {
          txp.lengthClassIndex(0);
        } else if (txp.RefLength <= 4389) {
          txp.lengthClassIndex(0);
        } else {
          txp.lengthClassIndex(0);
        }
      }
	    // ====== Done loading the transcripts from file
    }

    void loadTranscriptsFromFMD() {
	    bwaidx_t* idx_ = salmonIndex_->bwaIndex();
	    size_t numRecords = idx_->bns->n_seqs;
	    std::vector<Transcript> transcripts_tmp;
        //transcripts_tmp.reserve(numRecords);
        //transcripts_.reserve(numRecords);

	    fmt::print(stderr, "Index contained {} targets\n", numRecords);
	    //transcripts_.resize(numRecords);
	    for (auto i : boost::irange(size_t(0), numRecords)) {
		    uint32_t id = i;
		    char* name = idx_->bns->anns[i].name;
		    uint32_t len = idx_->bns->anns[i].len;
		    // copy over the length, then we're done.
		    transcripts_tmp.emplace_back(id, name, len);
	    }

	    std::sort(transcripts_tmp.begin(), transcripts_tmp.end(),
			    [](const Transcript& t1, const Transcript& t2) -> bool {
			    return t1.id < t2.id;
			    });


	    double alpha = 0.005;
	    char nucTab[256];
	    nucTab[0] = 'A'; nucTab[1] = 'C'; nucTab[2] = 'G'; nucTab[3] = 'T';
	    for (size_t i = 4; i < 256; ++i) { nucTab[i] = 'N'; }

        size_t tnum = 0;
	    // Load the transcript sequence from file
	    for (auto& t : transcripts_tmp) {
		    transcripts_.emplace_back(t.id, t.RefName.c_str(), t.RefLength, alpha);
		    /* from BWA */
		    uint8_t* rseq = nullptr;
		    int64_t tstart, tend, compLen, l_pac = idx_->bns->l_pac;
		    tstart  = idx_->bns->anns[t.id].offset;
		    tend = tstart + t.RefLength;
		    rseq = bns_get_seq(l_pac, idx_->pac, tstart, tend, &compLen);
		    if (compLen != t.RefLength) {
			    fmt::print(stderr,
					    "For transcript {}, stored length ({}) != computed length ({}) --- index may be corrupt. exiting\n",
					    t.RefName, compLen, t.RefLength);
			    std::exit(1);
		    }
		    std::string seq(t.RefLength, ' ');
		    if (rseq != 0) {
			    for (int64_t i = 0; i < compLen; ++i) { seq[i] = nucTab[rseq[i]]; }
		    }

            auto& txp = transcripts_.back();

            // allocate space for the new copy
            char* seqCopy = new char[seq.length()+1];
            std::strcpy(seqCopy, seq.c_str());
            txp.setSequenceOwned(seqCopy);
		    txp.setSAMSequenceOwned(salmon::stringtools::encodeSequenceInSAM(seq.c_str(), t.RefLength));

		    // Length classes taken from
		    // ======
		    // Roberts, Adam, et al.
		    // "Improving RNA-Seq expression estimates by correcting for fragment bias."
		    // Genome Biol 12.3 (2011): R22.
		    // ======
		    // perhaps, define these in a more data-driven way
		    if (t.RefLength <= 1334) {
			    txp.lengthClassIndex(0);
		    } else if (t.RefLength <= 2104) {
			    txp.lengthClassIndex(0);
		    } else if (t.RefLength <= 2988) {
			    txp.lengthClassIndex(0);
		    } else if (t.RefLength <= 4389) {
			    txp.lengthClassIndex(0);
		    } else {
			    txp.lengthClassIndex(0);
		    }
		    /*
		       std::cerr << "TS = " << t.RefName << " : \n";
		       std::cerr << seq << "\n VS \n";
		       for (size_t i = 0; i < t.RefLength; ++i) {
		       std::cerr << transcripts_.back().charBaseAt(i);
		       }
		       std::cerr << "\n\n";
		       */
		    free(rseq);
		    /* end BWA code */
            ++tnum;
	    }

	    // Since we have the de-coded reference sequences, we no longer need
	    // the encoded sequences, so free them.
	    /** TEST OPT **/
	    // free(idx_->pac); idx_->pac = nullptr;
	    /** END TEST OPT **/
	    transcripts_tmp.clear();
	    // ====== Done loading the transcripts from file
    }


    template <typename CallbackT>
    bool processReads(const uint32_t& numThreads, const SalmonOpts& sopt, CallbackT& processReadLibrary) {
        std::atomic<bool> burnedIn{totalAssignedFragments_ + numAssignedFragments_ > sopt.numBurninFrags};
        for (auto& rl : readLibraries_) {
            processReadLibrary(rl, salmonIndex_.get(), transcripts_, clusterForest(),
                               *(fragLengthDist_.get()), numAssignedFragments_,
                               numThreads, burnedIn);
        }
        return true;
    }

    ~ReadExperiment() {
        // ---- Get rid of things we no longer need --------
        // bwa_idx_destroy(idx_);
    }

    ClusterForest& clusterForest() { return *clusters_.get(); }

    std::string readFilesAsString() {
        std::stringstream sstr;
        size_t ln{0};
        size_t numReadLibraries{readLibraries_.size()};

        for (auto &rl : readLibraries_) {
            sstr << rl.readFilesAsString();
            if (ln++ < numReadLibraries) { sstr << "; "; }
        }
        return sstr.str();
    }

    uint64_t numAssignedFragsInFirstPass() {
        return numAssignedFragsInFirstPass_;
    }

    uint64_t numObservedFragsInFirstPass() {
        return numObservedFragsInFirstPass_;
    }

    double effectiveMappingRate() const {
        return effectiveMappingRate_;
    }

    void setEffectiveMappingRate(double emr) {
        effectiveMappingRate_ = emr;
    }

    std::vector<FragmentStartPositionDistribution>& fragmentStartPositionDistributions() {
        return fragStartDists_;
    }

    SequenceBiasModel& sequenceBiasModel() {
        return seqBiasModel_;
    }

    bool softReset() {
        if (quantificationPasses_ == 0) {
            numAssignedFragsInFirstPass_ = numAssignedFragments_;
            numObservedFragsInFirstPass_ = numObservedFragments_;
        }
        numObservedFragments_ = 0;
        totalAssignedFragments_ += numAssignedFragments_;
        numAssignedFragments_ = 0;
        quantificationPasses_++;
        return true;
    }

    bool reset() {
        namespace bfs = boost::filesystem;
        for (auto& rl : readLibraries_) {
            if (!rl.isRegularFile()) { return false; }
        }

        if (quantificationPasses_ == 0) {
            numAssignedFragsInFirstPass_ = numAssignedFragments_;
            numObservedFragsInFirstPass_ = numObservedFragments_;
        }

        numObservedFragments_ = 0;
        totalAssignedFragments_ += numAssignedFragments_;
        numAssignedFragments_ = 0;
        quantificationPasses_++;
        return true;
    }

    void summarizeLibraryTypeCounts(boost::filesystem::path& opath){
        LibraryFormat fmt1(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U);
        LibraryFormat fmt2(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U);

        std::ofstream ofile(opath.string());

        fmt::MemoryWriter errstr;

        auto log = spdlog::get("jointLog");

        uint64_t numFmt1{0};
        uint64_t numFmt2{0};
        uint64_t numAgree{0};
        uint64_t numDisagree{0};

        for (auto& rl : readLibraries_) {
            auto fmt = rl.format();
            auto& counts = rl.libTypeCounts();

            // If the format is un-stranded, check that
            // we have a similar number of mappings in both
            // directions and then aggregate the forward and
            // reverse counts.
            if (fmt.strandedness == ReadStrandedness::U) {
                std::vector<ReadStrandedness> strands;
                switch (fmt.orientation) {
                    case ReadOrientation::SAME:
                    case ReadOrientation::NONE:
                        strands.push_back(ReadStrandedness::S);
                        strands.push_back(ReadStrandedness::A);
                        break;
                    case ReadOrientation::AWAY:
                    case ReadOrientation::TOWARD:
                        strands.push_back(ReadStrandedness::AS);
                        strands.push_back(ReadStrandedness::SA);
                        break;
                }

                fmt1.type = fmt.type; fmt1.orientation = fmt.orientation;
                fmt1.strandedness = strands[0];
                fmt2.type = fmt.type; fmt2.orientation = fmt.orientation;
                fmt2.strandedness = strands[1];

                numFmt1 = 0;
                numFmt2 = 0;
                numAgree = 0;
                numDisagree = 0;

                for (size_t i = 0; i < counts.size(); ++i) {
                    if (i == fmt1.formatID()) {
                        numFmt1 = counts[i];
                    } else if (i == fmt2.formatID()) {
                        numFmt2 = counts[i];
                    } else {
                        numDisagree += counts[i];
                    }
                }
                numAgree = numFmt1 + numFmt2;
                double ratio = static_cast<double>(numFmt1) / (numFmt1 + numFmt2);

                if ( std::abs(ratio - 0.5) > 0.01) {
                    errstr << "NOTE: Read Lib [" << rl.readFilesAsString() << "] :\n";
                    errstr << "\nDetected a strand bias > 1\% in an unstranded protocol "
                           << "check the file: " << opath.string() << " for details\n";

                    log->warn() << errstr.str();
                    errstr.clear();
                }

                ofile << "========\n"
                      << "Read library consisting of files: "
                      << rl.readFilesAsString()
                      << "\n\n"
                      << "Expected format: " << rl.format()
                      << "\n\n"
                      << "# of consistent alignments: " << numAgree << "\n"
                      << "# of inconsistent alignments: " << numDisagree << "\n"
                      << "strand bias = " << ratio << " (0.5 is unbiased)\n"
                      << "# alignments with format " << fmt1 << ": " << numFmt1 << "\n"
                      << "# alignments with format " << fmt2 << ": " << numFmt2 << "\n"
                      << "\n========\n";
            } else {
                numAgree = 0;
                numDisagree = 0;

                for (size_t i = 0; i < counts.size(); ++i) {
                    if (i == fmt.formatID()) {
                        numAgree = counts[i];
                    } else {
                        numDisagree += counts[i];
                    }
                } // end for

                ofile << "========\n"
                      << "Read library consisting of files: "
                      << rl.readFilesAsString()
                      << "\n\n"
                      << "Expected format: " << rl.format()
                      << "\n\n"
                      << "# of consistent alignments: " << numAgree << "\n"
                      << "# of inconsistent alignments: " << numDisagree << "\n"
                      << "\n========\n";

            } //end else

            double disagreeRatio = static_cast<double>(numDisagree) / (numAgree + numDisagree);
            if (disagreeRatio > 0.05) {
                errstr << "NOTE: Read Lib [" << rl.readFilesAsString() << "] :\n";
                errstr << "\nGreater than 5\% of the alignments (but not, necessarily reads) "
                       << "disagreed with the provided library type; "
                       << "check the file: " << opath.string() << " for details\n";

                log->warn() << errstr.str();
                errstr.clear();
            }

            ofile << "---- counts for each format type ---\n";
            for (size_t i = 0; i < counts.size(); ++i) {
                ofile << LibraryFormat::formatFromID(i) << " : " << counts[i] << "\n";
            }
            ofile << "------------------------------------\n\n";
        }
        ofile.close();
    }

    std::vector<ReadLibrary>& readLibraries() { return readLibraries_; }
    FragmentLengthDistribution* fragmentLengthDistribution() const { return fragLengthDist_.get(); }

    void setGCFracForward(double fracForward) { gcFracFwd_ = fracForward; }

    double gcFracFwd() const { return gcFracFwd_; }
    double gcFracRC() const { return 1.0 - gcFracFwd_; }

    void setExpectedSeqBias(const std::vector<double>& expectedBiasIn) {
        expectedBias_ = expectedBiasIn;
    }

    std::vector<double>& expectedSeqBias() {
        return expectedBias_;
    }

    const std::vector<double>& expectedSeqBias() const {
        return expectedBias_;
    }

    void setExpectedGCBias(const std::vector<double>& expectedBiasIn) {
        expectedGC_ = expectedBiasIn;
    }

    std::vector<double>& expectedGCBias() {
        return expectedGC_;
    }

    const std::vector<double>& expectedGCBias() const {
        return expectedGC_;
    }

    const std::vector<double>& observedGC() const {
        return observedGC_;
    }

    std::vector<double>& observedGC() {
        return observedGC_;
    }

    ReadKmerDist<6, std::atomic<uint32_t>>& readBias() { return readBias_; }
    const ReadKmerDist<6, std::atomic<uint32_t>>& readBias() const { return readBias_; }

    private:
    /**
     * The file from which the alignments will be read.
     * This can be a SAM or BAM file, and can be a regular
     * file or a fifo.
     */
    std::vector<ReadLibrary> readLibraries_;
    /**
     * The file from which the transcripts are read.
     * This is expected to be a FASTA format file.
     */
    //boost::filesystem::path transcriptFile_;
    /**
     * The targets (transcripts) to be quantified.
     */
    std::vector<Transcript> transcripts_;
    /**
     * The index we've built on the set of transcripts.
     */
    std::unique_ptr<SalmonIndex> salmonIndex_{nullptr};
    //bwaidx_t *idx_{nullptr};
    /**
     * The cluster forest maintains the dynamic relationship
     * defined by transcripts and reads --- if two transcripts
     * share an ambiguously mapped read, then they are placed
     * in the same cluster.
     */
    std::unique_ptr<ClusterForest> clusters_;
    /**
      *
      *
      */
    std::vector<FragmentStartPositionDistribution> fragStartDists_;

    SequenceBiasModel seqBiasModel_;

    /** Keeps track of the number of passes that have been
     *  made through the alignment file.
     */
    std::atomic<uint64_t> numObservedFragments_{0};
    std::atomic<uint64_t> numAssignedFragments_{0};
    uint64_t totalAssignedFragments_{0};
    size_t quantificationPasses_{0};
    uint64_t numAssignedFragsInFirstPass_{0};
    uint64_t numObservedFragsInFirstPass_{0};
    uint64_t upperBoundHits_{0};
    double effectiveMappingRate_{0.0};
    SpinLock sl_;
    std::unique_ptr<FragmentLengthDistribution> fragLengthDist_;
    EquivalenceClassBuilder eqBuilder_;

    /** GC-fragment bias things **/
    // One bin for each percentage GC content
    double gcFracFwd_{-1.0};
    std::vector<double> observedGC_;
    std::vector<double> expectedGC_;

    /** Sequence specific bias things **/
    // Since multiple threads can touch this dist, we
    // need atomic counters.
    ReadKmerDist<6, std::atomic<uint32_t>> readBias_;
    std::vector<double> expectedBias_;
};

#endif // EXPERIMENT_HPP
