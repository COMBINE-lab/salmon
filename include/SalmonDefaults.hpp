#ifndef SALMON_DEFAULTS_HPP
#define SALMON_DEFAULTS_HPP

#include <thread>

namespace salmon {
namespace defaults {
  // general
  constexpr const bool discardOrphansAln{false};
  constexpr const bool discardOrphansQuasi{false};
  constexpr const bool allowOrphansFMD{false};
  constexpr const bool seqBiasCorrect{false};
  constexpr const bool gcBiasCorrect{false};
  constexpr const bool posBiasCorrect{false};
  constexpr const uint32_t numThreads{8};
  constexpr const double incompatPrior{0.0};
  constexpr const char quasiMappingDefaultFile[] = "";
  constexpr const char quasiMappingImplicitFile[] = "-";
  constexpr const bool metaMode{false};
  constexpr const bool disableMappingCache{true};

  // advanced
  constexpr const bool validateMappings{true};
  constexpr const bool disableSA{false};
  constexpr const float consensusSlack{0.35};
  constexpr const double minScoreFraction{0.65};
  constexpr const double pre_merge_chain_sub_thresh{0.75};
  constexpr const double post_merge_chain_sub_thresh{0.9};
  constexpr const double orphan_chain_sub_thresh{0.95};
  constexpr const double scoreExp{1.0};
  constexpr const int16_t matchScore{2};
  constexpr const int16_t mismatchPenalty{-4};
  constexpr const int16_t gapOpenPenalty{6};
  constexpr const int16_t gapExtendPenalty{2};
  constexpr const int32_t dpBandwidth{15};
  constexpr const uint32_t mismatchSeedSkip{3};
  constexpr const bool disableChainingHeuristic{false};
  constexpr const bool eqClassMode{false};
  constexpr const bool hardFilter{false};
  constexpr const bool mimicStrictBT2{false};
  constexpr const bool mimicBT2{false};
  constexpr const bool softclip{false};
  constexpr const bool softclipOverhangs{false};
  constexpr const bool fullLengthAlignment{false};
  constexpr const bool allowDovetail{false};
  constexpr const bool recoverOrphans{false};
  constexpr const bool alternativeInitMode{false};
  constexpr const char auxDir[] = "aux_info";
  constexpr const bool consistentHits{false};
  constexpr const bool skipQuant{false};
  constexpr const bool dumpEq{false};
  constexpr const bool dumpEqWeights{false};
  constexpr const bool fasterMapping{false};
  constexpr const uint32_t minAssignedFrags{10};
  constexpr const bool reduceGCMemory{false};
  constexpr const uint32_t biasSpeedSamp{5};
  constexpr const uint32_t maxFragLength{1000};
  constexpr const double fragLenPriorMean{250.0};
  constexpr const double fragLenPriorSD{25.0};
  constexpr const double ffactor{0.65};
  constexpr const uint32_t maxSMEMOccs{200};
  constexpr const bool initUniform{false};
  constexpr const uint32_t maxReadOccs{200};
  constexpr const uint32_t maxOccsPerHit{1000};
  constexpr const bool noLengthCorrection{false};
  constexpr const bool noEffectiveLengthCorrection{false};
  constexpr const bool noFragLengthDist{false};
  constexpr const bool noSingleFragProb{false};
  constexpr const bool noBiasLengthThreshold{false};
  constexpr const uint32_t numBiasSamples{2000000};
  constexpr const uint32_t numBurninFrags{5000000};
  constexpr const uint32_t numPreBurninFrags{5000};
  constexpr const bool useEM{false};
  constexpr const bool useVBOpt{true};
  constexpr const uint32_t sigDigits{3};
  constexpr const uint32_t rangeFactorizationBins{4};
  constexpr const uint32_t numGibbsSamples{0};
  constexpr const bool noGammaDraw{false};
  constexpr const bool bootstrapReproject{false};
  constexpr const uint32_t thinningFactor{16};
  constexpr const uint32_t numBootstraps{0};
  constexpr const bool quiet{false};
  constexpr const bool perTranscriptPrior{true};
  constexpr const bool perNucleotidePrior{false};
  constexpr const double vbPrior{1e-2};
  constexpr const bool writeOrphanLinks{false};
  constexpr const bool writeUnmappedNames{false};
  constexpr const double quasiCoverage{0.0};
  constexpr const double decoyThreshold{1.0};
  const std::string hitFilterPolicyStr{"AFTER"};
  constexpr const double minAlnProb{1e-5};
  const std::string auxTargetFile{""};

   // FMD-specific options
  constexpr const int fmdMinSeedLength{19};
  constexpr const bool fmdSensitive{false};
  constexpr const bool fmdExtraSeedPass{false};
  constexpr const double fmdCoverageThresh{0.7};
  constexpr const int fmdSplitWidth{0};
  constexpr const bool fmdSplitSpanningSeeds{false};

  // options not shown by default
  constexpr const size_t numFragGCBins{25};
  constexpr const size_t numConditionalGCBins{3};
  constexpr const size_t numRequiredFrags{50000000}; // deprecated
  const uint32_t maxHashResizeThreads{std::thread::hardware_concurrency()};

  // experimental / testing
  constexpr const bool noRichEqClasses{false};
  constexpr const bool noFragLengthFactor{false};
  constexpr const bool rankEqClasses{false};
  constexpr const bool dontExtrapolateCounts{false};
  constexpr const bool disableAlignmentCache{false};

  // alignment-based mode
  //constexpr const bool useErrorModel{true};
  constexpr const bool noErrorModel{false};
  constexpr const bool useMassBanking{false};
  constexpr const bool gencodeRef{false};
  constexpr const uint32_t mappingCacheMemoryLimit{2000000};
  constexpr const uint32_t numErrorBins{6};
  constexpr const bool sampleOutput{false};
  constexpr const bool sampleUnaligned{false};
}
}

namespace alevin {
namespace defaults {
  constexpr const bool naiveEqclass{false};
  constexpr const bool noDedup{false};
  constexpr const bool noWhitelist{false};
  constexpr const bool txpLevel{false};
  constexpr const bool eqClassLevel{false};
  constexpr const bool isDropseq{false};
  constexpr const bool isChromium{false};
  constexpr const bool isChromiumV3{false};
  constexpr const bool isInDrop{false};
  constexpr const bool isGemcode{false};
  constexpr const bool isCITESeq{false};
  constexpr const bool isCELSeq{false};
  constexpr const bool isCELSeq2{false};
  constexpr const bool isQuartzSeq2{false};
  constexpr const bool noQuant{false};
  constexpr const bool dumpFQ{false};
  constexpr const bool dumpArborescences{false};
  constexpr const bool dumpBarcodeEq{false};
  constexpr const bool dumpFeatures{false};
  constexpr const bool dumpBFH{false};
  constexpr const bool dumpUmiGraph{false};
  constexpr const bool dumpCellEq{false};
  constexpr const bool dumpMtx{false};
  constexpr const bool noEM{false};
  constexpr const bool debug{true};
  constexpr const bool just_align{false};
  constexpr const bool sketch_mode{false};
  constexpr const uint32_t trimRight{0};
  constexpr const uint32_t numBootstraps{0};
  constexpr const uint32_t numGibbsSamples{0};
  constexpr const uint32_t lowRegionMinNumBarcodes{200};
  constexpr const uint32_t maxNumBarcodes{100000};
  constexpr const uint32_t expectCells{0};
  constexpr const uint32_t forceCells{0};
  constexpr const uint32_t freqThreshold{10};
  constexpr const uint32_t umiEditDistance{1};
  constexpr const float consensusSlack{0.6};
  constexpr const double keepCBFraction{0.0};
  constexpr const double minScoreFraction{0.87};
  constexpr const double vbemNorm{0.0};
}
}

#endif // SALMON_DEFAULTS_HPP
