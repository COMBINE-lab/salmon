#include <vector>
#include <set> // std::set
#include <deque> // std::deque
//#include <queue> // std::priority_queue
#include <memory>
#include "spdlog/spdlog.h"
#include "sparsepp/spp.h"
#include <thread>
#include <mutex>
//#include "tbb/tbb.h"


#include "cedar/Taxa.h"
struct ReadInfo {
    std::string rid;
    uint32_t cnt = 0;
    uint32_t len = 0;
    std::vector<TaxaNode> mappings;
};
#include "cedar/EquivalenceClassBuilder.hpp"
#include "cedar/PuffMappingReader.hpp"
#include "cedar/SAMReader.hpp"



constexpr uint32_t ALIGNMENTS_PER_BATCH{20};
namespace util {
/*
 * Use atomic compare-and-swap to update val to
 * val + inc (*in log-space*).  Update occurs in a loop in case other
 * threads update in the meantime.
 */
    /*inline void incLoopLog(tbb::atomic<double> &val, double inc) {
        double oldMass = val.load();
        double returnedMass = oldMass;
        double newMass{salmon::math::LOG_0};
        do {
            oldMass = returnedMass;
            newMass = salmon::math::logAdd(oldMass, inc);
            returnedMass = val.compare_and_swap(newMass, oldMass);
        } while (returnedMass != oldMass);
    }*/

/*
 * Same as above, but overloaded for "plain" doubles
 */
    inline void incLoop(double &val, double inc) { val += inc; }

/*
 * Use atomic compare-and-swap to update val to
 * val + inc.  Update occurs in a loop in case other
 * threads update in the meantime.
 */
    inline void incLoop(std::atomic<double> &val, double inc) {
        double oldMass = val.load();
        //double returnedMass = oldMass;
        double newMass{oldMass + inc};
        bool successful = false;
        do {
            //oldMass = val.load();
            newMass = oldMass + inc;
            //returnedMass = val.compare_and_swap(newMass, oldMass);
            successful = val.compare_exchange_weak(oldMass, newMass);
        } while (!successful);//while (returnedMass != oldMass);
//        std::cerr << inc << "\n";
    }

    inline void update(std::atomic<double> &val, double newval) {      // Update x and return old value of x.
        double oldval = val;
        if (oldval == newval) return;
        do {
            // Read globalX
            oldval = val;
            // Store new value if another thread has not changed globalX.
        } while (val.compare_exchange_strong/*compare_and_swap*/(newval, oldval) != oldval);
    }

}

struct Stats {
    uint64_t readCnt{0};
    double globalprobsum{0};
    uint64_t totalReadCnt{0};
    uint64_t seqNotFound{0};
    uint64_t totalMultiMappedReads{0};
    uint64_t totalUnmappedReads{0};
    uint64_t totalReadsNotPassingCond{0};
    uint64_t conflicting{0};
    uint64_t discordantMappings{0};
    void update(Stats& s) {
        readCnt+=s.readCnt;
        globalprobsum+=s.globalprobsum;
        totalReadCnt+=s.totalReadCnt;
        seqNotFound+=s.seqNotFound;
        totalMultiMappedReads+=s.totalMultiMappedReads;
        totalUnmappedReads+=s.totalUnmappedReads;
        totalReadsNotPassingCond+=s.totalReadsNotPassingCond;
        conflicting+=s.conflicting;
        discordantMappings+=s.discordantMappings;
    }
};

template<class ReaderType, class FileReaderType>
class Cedar {
public:
    Cedar(std::string &taxonomyTree_filename, std::string &
    refId2TaxId_filename, std::string pruneLevelIn, double filteringThresholdIn,
          bool flatAbund,
          std::shared_ptr<spdlog::logger> loggerIn);

    void run(std::string mapperOutput_filename,
             bool requireConcordance,
             size_t maxIter,
             double eps,
             double minCnt,
             std::string &output_filename,
             bool onlyUniq,
             bool onlyPerf,
             uint32_t segmentSize,
             uint32_t rangeFactorizationBins,
             uint32_t numThreads);

private:

    void loadMappingInfo(std::string mapperOutput_filename, bool requireConcordance, bool onlyUniq, bool onlyPerfect,
                         uint32_t segmentSize, uint32_t rangeFactorizationBins);

    void loadMappingInfo(std::string mapperOutput_filename,
                         bool requireConcordance,
                         bool onlyUniq,
                         bool onlyPerfect,
                         uint32_t segmentSize,
                         uint32_t rangeFactorizationBins,
                         uint32_t nThreads);

    void processAlignmentBatch(uint32_t threadID,
                               std::vector<ReadInfo> &alignmentGrp,
                               Stats &stats,
                               EquivalenceClassBuilder& eqb,
                               std::mutex &iomutex,
                               bool requireConcordance,
                               bool onlyUniq,
                               bool onlyPerfect,
                               uint32_t segmentSize,
                               uint32_t rangeFactorizationBins,
                               bool getReadName);

    bool applySetCover(std::vector<double> &strainCnt, std::vector<bool> &strainValid,
                       std::vector<bool> &strainPotentiallyRemovable, double minCnt, bool canHelp,
                       bool verbose = false);

    bool basicEM(size_t maxIter, double eps, double minCnt, uint32_t numThreads, bool verbose = false);

    void serialize(std::string &output_filename);

    void serializeFlat(std::string &output_filename);

    void calculateCoverage();

    spp::sparse_hash_map<uint32_t, TaxaNode> taxaNodeMap;
    spp::sparse_hash_map<std::string, uint32_t> refId2taxId;
    spp::sparse_hash_map<uint32_t, uint32_t> seqToTaxMap;
    Rank pruningLevel = Rank::SPECIES;
    std::set<uint64_t> activeTaxa;
    uint64_t rootId = 1;
    double filteringThreshold = 0;
    spp::sparse_hash_map<uint64_t, double> strain;
    uint64_t readCnt = 0;
    bool flatAbund = false;
    std::vector<std::vector<std::pair<uint64_t, double>>> readPerStrainProb;
    EquivalenceClassBuilder eqb;
    ReaderType fileReader;
    std::shared_ptr<spdlog::logger> logger;

    spp::sparse_hash_map<uint32_t, double> cov;

    spp::sparse_hash_map<uint64_t, double> strain_coverage;
    spp::sparse_hash_map<uint64_t, std::vector<uint32_t> > strain_coverage_bins;
};

