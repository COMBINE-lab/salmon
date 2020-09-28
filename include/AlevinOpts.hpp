#ifndef ALEVIN_OPTS_HPP
#define ALEVIN_OPTS_HPP

#include <string>
#include "spdlog/spdlog.h"
#include <boost/filesystem.hpp>

/**
  * A structure to hold some common options used
  * by Salmon so that we don't have to pass them
  * all around as separate arguments.
  */

enum class BarcodeEnd { FIVE = 5, THREE = 3 };
enum class Sequence { BARCODE, UMI };

template <class protocolT>
struct AlevinOpts {
  AlevinOpts(): numParsingThreads(1),
                numConsumerThreads(2),
                useVBEM{false},
                numNoMapCB(0),
                initUniform{false}{}

  //IUPAC code for the cell-barcodes
  std::string iupac;
  //protocol to use for single cell
  protocolT protocol;
  //dump barcodes fq files
  bool dumpfq;
  // dump Arborescence Fragment Counts
  bool dumpArborescences;
  //dump CB features for whitelisting
  bool dumpfeatures;
  //eqclass level barcode count
  bool dumpBarcodeEq;
  //dump cellvtxp matrix in csv
  bool dumpMtx;
  // dump big fishing hash
  bool dumpBFH;
  // dump per cell level umi-graph
  bool dumpUmiGraph;
  // dump per cell level de-duplicated 
  // equivalence class
  bool dumpCellEq;
  //Stop progress sumps
  bool quiet;
  // just perform alignment and produce 
  // an output directory with a RAD file. 
  bool just_align;
  bool sketch_mode;
  //flag for deduplication
  bool noDedup;
  //flag for not performing external whitelisting
  bool noWhitelist;
  //Number of generator threads
  uint32_t numParsingThreads;
  //Number of consumer threads
  uint32_t numConsumerThreads;
  //total num of threads
  uint32_t numThreads;
  //threshold for the frequency of the barcodes
  uint32_t freqThreshold;
  // maximum allowable edit distance for UMI collapsing
  uint32_t umiEditDistance;
  // sequences to trim from right in the read sequences
  uint32_t trimRight;
  //no downstream salmon quant
  bool noQuant;
  // don't run EM flag
  bool noEM;
  // use vbem
  bool useVBEM;
  // Avoid segfaults based on no whitelist mapping
  bool debug;
  // perform naive deduplication
  bool naiveEqclass;
  // perform eqclass level analysis instead of gene or txp level minsets
  bool eqClassLevel;
  // perform txp level analysis instead of gene level
  bool txpLevel;
  // initialize EM with uniform prior
  bool initUniform;
  //number of cells
  uint32_t numCells;
  // minimum number of CB to use for low confidence region
  uint32_t lowRegionMinNumBarcodes;
  // maximum number of barcodes to use
  uint32_t maxNumBarcodes;
  // number of bootstraps to perform
  uint32_t numBootstraps;
  // number of gibbs samples to perform
  uint32_t numGibbsSamples;
  // force the number of cells
  uint32_t forceCells;
  // define a close upper bound on expected number of cells
  uint32_t expectCells;

  // Related to the logger
  std::shared_ptr<spdlog::logger> jointLog{nullptr};

  // barcode output directory
  boost::filesystem::path outputDirectory;
  // barcode white-list File path
  boost::filesystem::path whitelistFile;
  // alevin matrix eds file with vbem Priors
  boost::filesystem::path vbemPriorFile;
  // barcode mitochondrial genes File path
  boost::filesystem::path mRnaFile;
  // barcode ribosomal gene File path
  boost::filesystem::path rRnaFile;
  // Txp to gene map tsv file
  boost::filesystem::path geneMapFile;
  // Bfh file
  boost::filesystem::path bfhFile;

  //meta-info related tags
  uint32_t totalReads;
  uint32_t totalUsedReads;
  uint32_t readsThrown;

  uint32_t totalCBs;
  uint32_t totalUsedCBs;

  uint32_t kneeCutoff;
  uint32_t intelligentCutoff;
  uint32_t totalLowConfidenceCBs;
  uint32_t numFeatures;
  uint32_t numNoMapCB;

  uint32_t eqReads;
  uint32_t noisyUmis;
  double mappingRate;
  double keepCBFraction;
  double vbemNorm;

  uint32_t totalDedupUMIs;
  uint32_t totalExpGenes;
};

#endif // ALEVIN_OPTS_HPP
