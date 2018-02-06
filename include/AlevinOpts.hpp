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
                freqThreshold(10){}

  //IUPAC code for the cell-barcodes
  std::string iupac;
  //protocol to use for single cell
  protocolT protocol;
  //Dump soft-mapping of corrected barcodes
  bool dumpBarcodeMap;
  //umi_tools supported whitelist
  bool dumpUmiToolsMap;
  //dump barcodes fq files
  bool dumpfq;
  //dump CB features for whitelisting
  bool dumpfeatures;
  //flag for using barcoding or just single cell dedup
  bool nobarcode;
  //eqclass level barcode count
  bool dumpBarcodeEq;
  //Stop progress sumps
  bool quiet;
  //flag for deduplication
  bool dedup;
  //Number of generator threads
  uint32_t numParsingThreads;
  //Number of consumer threads
  uint32_t numConsumerThreads;
  //total num of threads
  uint32_t numThreads;
  //threshold for the frequency of the barcodes
  uint32_t freqThreshold;
  //no downstream salmon quant
  bool noQuant;
  // don't run EM flag
  bool noEM;
  //do hard-assignment of error bcs
  bool noSoftMap;
  //number of cells
  uint32_t numCells;

  // Related to the logger
  std::shared_ptr<spdlog::logger> jointLog{nullptr};

  // barcode output directory
  boost::filesystem::path outputDirectory;
  // barcode white-list File path
  boost::filesystem::path whitelistFile;
  // Txp to gene map tsv file
  boost::filesystem::path geneMapFile;
};

#endif // ALEVIN_OPTS_HPP
