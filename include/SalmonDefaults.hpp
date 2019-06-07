#ifndef SALMON_DEFAULTS_HPP
#define SALMON_DEFAULTS_HPP

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
  constexpr const bool validateMappings{false};
  constexpr const float consensusSlack{0.0};
  constexpr const double minScoreFraction{0.65};
  constexpr const int16_t matchScore{2};
  constexpr const int16_t mismatchPenalty{-4};
  constexpr const int16_t gapOpenPenalty{4};
  constexpr const int16_t gapExtendPenalty{2};
  constexpr const int32_t dpBandwidth{15};
  constexpr const bool hardFilter{false};
  constexpr const bool mimicStrictBT2{false};
  constexpr const bool mimicBT2{false};
  constexpr const bool allowDovetail{false};
  constexpr const bool recoverOrphans{false};
  constexpr const int32_t maxMMPExtension{7};
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
  constexpr const uint32_t fragLenPriorMean{250};
  constexpr const uint32_t fragLenPriorSD{25};
  constexpr const double ffactor{0.65};
  constexpr const uint32_t maxSMEMOccs{200};
  constexpr const bool initUniform{false};
  constexpr const uint32_t maxReadOccs{200};
  constexpr const bool noLengthCorrection{false};
  constexpr const bool noEffectiveLengthCorrection{false};
  constexpr const bool noFragLengthDist{false};
  constexpr const bool noBiasLengthThreshold{false};
  constexpr const uint32_t numBiasSamples{2000000};
  constexpr const uint32_t numBurninFrags{5000000};
  constexpr const uint32_t numPreBurninFrags{5000};
  constexpr const bool useEM{false};
  constexpr const bool useVBOpt{true};
  constexpr const uint32_t sigDigits{3};
  constexpr const uint32_t rangeFactorizationBins{0};
  constexpr const uint32_t numGibbsSamples{0};
  constexpr const bool noGammaDraw{false};
  constexpr const bool bootstrapReproject{false};
  constexpr const uint32_t thinningFactor{16};
  constexpr const uint32_t numBootstraps{0};
  constexpr const bool quiet{false};
  constexpr const bool perTranscriptPrior{false};
  constexpr const double vbPrior{1e-5};
  constexpr const bool writeOrphanLinks{false};
  constexpr const bool writeUnmappedNames{false};
  constexpr const double quasiCoverage{0.0};

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
  constexpr const uint32_t maxHashResizeThreads{std::numeric_limits<uint32_t>::max()};

  // experimental / testing
  constexpr const bool noRichEqClasses{false};
  constexpr const bool noFragLengthFactor{false};
  constexpr const bool rankEqClasses{false};
  constexpr const bool dontExtrapolateCounts{false};

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
  constexpr const bool txpLevel{false};
  constexpr const bool eqClassLevel{false};
  constexpr const bool isDropseq{false};
  constexpr const bool isChromium{false};
  constexpr const bool isChromiumV3{false};
  constexpr const bool isInDrop{false};
  constexpr const bool isGemcode{false};
  constexpr const bool isCELSeq{false};
  constexpr const bool isCELSeq2{false};
  constexpr const bool noQuant{false};
  constexpr const bool dumpFQ{false};
  constexpr const bool dumpBarcodeEq{false};
  constexpr const bool dumpFeatures{false};
  constexpr const bool dumpBFH{false};
  constexpr const bool dumpUmiGraph{false};
  constexpr const bool dumpMtx{false};
  constexpr const bool noEM{false};
  constexpr const bool debug{true};
  constexpr const uint32_t trimRight{0};
  constexpr const uint32_t numBootstraps{0};
  constexpr const uint32_t lowRegionMinNumBarcodes{200};
  constexpr const uint32_t maxNumBarcodes{100000};
  constexpr const double minScoreFraction{0.87};
  constexpr const float consensusSlack{0.6};
  constexpr const uint32_t expectCells{0};
  constexpr const uint32_t forceCells{0};
  constexpr const double keepCBFraction{0.0};
  constexpr const uint32_t freqThreshold{10};
}
}

#endif // SALMON_DEFAULTS_HPP
