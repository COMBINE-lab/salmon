#ifndef ALEVIN_SALMON_OPT_HPP
#define ALEVIN_SALMON_OPT_HPP

#include "SalmonOpts.hpp"
#include "AlevinOpts.hpp"
#include "SalmonUtils.hpp"

template <typename ProtocolT>
void setSalmonOpt(SalmonOpts& sopt, AlevinOpts<ProtocolT>& aopt){

  //IMPORTANT FLAG changes have to check this up with Rob
  if (aopt.numThreads > 3){
    sopt.numThreads = aopt.numThreads + 3;
  }
  else{
    sopt.numThreads = aopt.numThreads;
  }
  sopt.allowOrphans = true; //impt
  sopt.biasCorrect = false; //impt
  sopt.posBiasCorrect = false; //impt
  sopt.gcBiasCorrect = false; //impt
  sopt.noRichEqClasses = true; //impt
  sopt.initUniform = true; //impt
  sopt.noLengthCorrection = true; //impt
  sopt.noEffectiveLengthCorrection = true; //impt
  sopt.noFragLengthDist = false; //impt
  sopt.fragLenDistPriorMean = 250; //impt
  sopt.fragLenDistPriorSD = 25; //impt
  sopt.noBiasLengthThreshold = true; //impt

  sopt.jointLog = aopt.jointLog;
  sopt.incompatPrior = std::log(1e-20);
  sopt.ignoreIncompat = false;
  sopt.meta = false;
  sopt.disableMappingCache = false;
  sopt.alternativeInitMode = false;
  sopt.auxDir = "aux_info";
  sopt.consistentHits = false;
  sopt.dumpEq = false;
  sopt.dumpEqWeights = false;
  sopt.fasterMapping = false;
  sopt.gcSampFactor = 1;
  sopt.pdfSampFactor = 1;
  sopt.strictIntersect = false;
  sopt.fragLenDistMax = 1000;
  sopt.forgettingFactor = 0.65;
  sopt.maxReadOccs = 100;
  sopt.numBurninFrags = 5000000;
  sopt.numPreBurninFrags = 1000000;
  sopt.useVBOpt = false;
  sopt.numGibbsSamples = 0;
  sopt.numBootstraps = 0;
  sopt.thinningFactor = 16;
  sopt.quiet = false;
  sopt.perTranscriptPrior = false;
  sopt.vbPrior = 1e-3;
  sopt.writeOrphanLinks = false;
  sopt.writeUnmappedNames = false;
  sopt.quasiCoverage = 0.0;
  sopt.sensitive = false;
  sopt.extraSeedPass = false;
  sopt.splitSpanningSeeds = false;
  sopt.numFragGCBins = 25;
  sopt.numConditionalGCBins = 3;
  sopt.numRequiredFragments = 50000000;
  sopt.noFragLenFactor = false;
  sopt.rankEqClasses = false;
  sopt.dontExtrapolateCounts = false;
  sopt.useFSPD = false;
  sopt.quantMode = SalmonQuantMode::MAP;
  sopt.numBiasSamples = 2000000;
  sopt.runStartTime = salmon::utils::getCurrentTimeAsString();
  sopt.outputDirectory = aopt.outputDirectory;
  //sopt.indexDirectory = aopt.indexDirectory;

  if (sopt.gcBiasCorrect and !sopt.biasCorrect) {
    sopt.numConditionalGCBins = 1;
  }

}

#endif // ALEVIN_SALMON_OPT
