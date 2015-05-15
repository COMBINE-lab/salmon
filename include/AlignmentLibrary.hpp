#ifndef ALIGNMENT_LIBRARY_HPP
#define ALIGNMENT_LIBRARY_HPP

extern "C" {
#include "io_lib/scram.h"
#include "io_lib/os.h"
#undef max
#undef min
}

// Our includes
#include "ClusterForest.hpp"
#include "Transcript.hpp"
#include "BAMQueue.hpp"
#include "SalmonUtils.hpp"
#include "LibraryFormat.hpp"
#include "SalmonOpts.hpp"
#include "FragmentLengthDistribution.hpp"
#include "FragmentStartPositionDistribution.hpp"
#include "AlignmentGroup.hpp"
#include "ErrorModel.hpp"
#include "AlignmentModel.hpp"
#include "FASTAParser.hpp"
#include "concurrentqueue.h"

// Boost includes
#include <boost/filesystem.hpp>

// Standard includes
#include <vector>
#include <memory>
#include <functional>

template <typename T>
class NullFragmentFilter;

/**
  *  This class represents a library of alignments used to quantify
  *  a set of target transcripts.  The AlignmentLibrary contains info
  *  about both the alignment file and the target sequence (transcripts).
  *  It is used to group them together and track information about them
  *  during the quantification procedure.
  */
template <typename FragT>
class AlignmentLibrary {

    public:

    AlignmentLibrary(std::vector<boost::filesystem::path>& alnFiles,
                     const boost::filesystem::path& transcriptFile,
                     LibraryFormat libFmt,
                     SalmonOpts& salmonOpts) :
        alignmentFiles_(alnFiles),
        transcriptFile_(transcriptFile),
        libFmt_(libFmt),
        transcripts_(std::vector<Transcript>()),
	fragStartDists_(5),
        quantificationPasses_(0) {
            namespace bfs = boost::filesystem;

            // Make sure the alignment file exists.
            for (auto& alignmentFile : alignmentFiles_) {
                if (!bfs::exists(alignmentFile)) {
                    std::stringstream ss;
                    ss << "The provided alignment file: " << alignmentFile <<
                        " does not exist!\n";
                    throw std::invalid_argument(ss.str());
                }
            }


            // Make sure the transcript file exists.
            if (!bfs::exists(transcriptFile_)) {
                std::stringstream ss;
                ss << "The provided transcript file: " << transcriptFile_ <<
                    " does not exist!\n";
                throw std::invalid_argument(ss.str());
            }

            // The alignment file existed, so create the alignment queue
            size_t numParseThreads = salmonOpts.numParseThreads;
            std::cerr << "parseThreads = " << numParseThreads << "\n";
            bq = std::unique_ptr<BAMQueue<FragT>>(new BAMQueue<FragT>(alnFiles, libFmt_, numParseThreads,
                                                                      salmonOpts.mappingCacheMemoryLimit));

            std::cerr << "Checking that provided alignment files have consistent headers . . . ";
            if (! salmon::utils::headersAreConsistent(bq->headers()) ) {
                std::stringstream ss;
                ss << "\nThe multiple alignment files provided had inconsistent headers.\n";
                ss << "Currently, we require that if multiple SAM/BAM files are provided,\n";
                ss << "they must have identical @SQ records.\n";
                throw std::invalid_argument(ss.str());
            }
            std::cerr << "done\n";

            SAM_hdr* header = bq->header();

            // The transcript file existed, so load up the transcripts
            double alpha = 0.005;
            for (size_t i = 0; i < header->nref; ++i) {
                transcripts_.emplace_back(i, header->ref[i].name, header->ref[i].len, alpha);
            }

            FASTAParser fp(transcriptFile.string());

            fmt::print(stderr, "Populating targets from aln = {}, fasta = {} . . .",
                       alnFiles.front(), transcriptFile_);
            fp.populateTargets(transcripts_);
	    for (auto& txp : transcripts_) {
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
            fmt::print(stderr, "done\n");

            // Create the cluster forest for this set of transcripts
            clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));

            // Initialize the fragment length distribution
            size_t maxFragLen = 800;
            size_t meanFragLen = 200;
            size_t fragLenStd = 80;
            size_t fragLenKernelN = 4;
            double fragLenKernelP = 0.5;
            flDist_.reset(new
                    FragmentLengthDistribution(
                    1.0, maxFragLen,
                    meanFragLen, fragLenStd,
                    fragLenKernelN,
                    fragLenKernelP, 1)
                    );

            alnMod_.reset(new AlignmentModel(1.0, salmonOpts.numErrorBins));
            alnMod_->setLogger(salmonOpts.jointLog);

            // Start parsing the alignments
           NullFragmentFilter<FragT>* nff = nullptr;
           bool onlyProcessAmbiguousAlignments = false;
           bq->start(nff, onlyProcessAmbiguousAlignments);
        }

    std::vector<Transcript>& transcripts() { return transcripts_; }

    inline bool getAlignmentGroup(AlignmentGroup<FragT>*& ag) { return bq->getAlignmentGroup(ag); }

    //inline t_pool* threadPool() { return threadPool_.get(); }

    inline SAM_hdr* header() { return bq->header(); }

    inline std::vector<FragmentStartPositionDistribution>& fragmentStartPositionDistributions() {
	    return fragStartDists_;
    }

    inline FragmentLengthDistribution& fragmentLengthDistribution() {
        return *flDist_.get();
    }

    inline AlignmentModel& alignmentModel() {
        return *alnMod_.get();
    }

//    inline tbb::concurrent_queue<FragT*>& fragmentQueue() {
    inline tbb::concurrent_queue<FragT*>& fragmentQueue() {
        return bq->getFragmentQueue();
    }

//    inline tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*>& alignmentGroupQueue() {
    inline moodycamel::ConcurrentQueue<AlignmentGroup<FragT*>*>& alignmentGroupQueue() {
        return bq->getAlignmentGroupQueue();
    }

    inline BAMQueue<FragT>& getAlignmentGroupQueue() { return *bq.get(); }

    inline size_t upperBoundHits() { return bq->numMappedReads(); }
    inline size_t numMappedReads() { return bq->numMappedReads(); }
    inline size_t numUniquelyMappedReads() { return bq->numUniquelyMappedReads(); }

    //const boost::filesystem::path& alignmentFile() { return alignmentFile_; }

    ClusterForest& clusterForest() { return *clusters_.get(); }

    template <typename FilterT>
    bool reset(bool incPasses=true, FilterT filter=nullptr, bool onlyProcessAmbiguousAlignments=false) {
        namespace bfs = boost::filesystem;

        for (auto& alignmentFile : alignmentFiles_) {
            if (!bfs::is_regular_file(alignmentFile)) {
                return false;
            }
        }

        bq->reset();
        bq->start(filter, onlyProcessAmbiguousAlignments);
        if (incPasses) {
            quantificationPasses_++;
            fmt::print(stderr, "Current iteration = {}\n", quantificationPasses_);
        }
        return true;
    }

    inline LibraryFormat format() { return libFmt_; }

    private:
    /**
     * The file from which the alignments will be read.
     * This can be a SAM or BAM file, and can be a regular
     * file or a fifo.
     */
    std::vector<boost::filesystem::path> alignmentFiles_;
    /**
     * The file from which the transcripts are read.
     * This is expected to be a FASTA format file.
     */
    boost::filesystem::path transcriptFile_;
    /**
     * Describes the expected format of the sequencing
     * fragment library.
     */
    LibraryFormat libFmt_;
    /**
     * The targets (transcripts) to be quantified.
     */
    std::vector<Transcript> transcripts_;
    /**
     * A pointer to the queue from which the fragments
     * will be read.
     */
    //std::unique_ptr<t_pool, std::function<void(t_pool*)>> threadPool_;
    std::unique_ptr<BAMQueue<FragT>> bq;

    /**
     * The cluster forest maintains the dynamic relationship
     * defined by transcripts and reads --- if two transcripts
     * share an ambiguously mapped read, then they are placed
     * in the same cluster.
     */
    std::unique_ptr<ClusterForest> clusters_;

    /**
      * The emperical fragment start position distribution
      */
    std::vector<FragmentStartPositionDistribution> fragStartDists_;

    /**
     * The emperical fragment length distribution.
     *
     */
    std::unique_ptr<FragmentLengthDistribution> flDist_;
    /**
      * The emperical error model
      */
    std::unique_ptr<AlignmentModel> alnMod_;

    /** Keeps track of the number of passes that have been
     *  made through the alignment file.
     */
    size_t quantificationPasses_;
};

#endif // ALIGNMENT_LIBRARY_HPP
