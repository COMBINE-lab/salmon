/**
>HEADER
    Copyright (c) 2013, 2014, 2015, 2016 Rob Patro rob.patro@cs.stonybrook.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#include <random>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <iterator>
#include <map>
#include <mutex>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

// Boost Includes
#include <boost/container/flat_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// logger includes
#include "spdlog/spdlog.h"

#include "cereal/types/vector.hpp"

// utility includes
#include "nonstd/string_view.hpp"
#include "nonstd/optional.hpp"

//alevin include
#include "AlevinOpts.hpp"
#include "AlevinHash.hpp"
#include "AlevinUtils.hpp"
#include "BarcodeModel.hpp"
#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"
#include "tsl/array_map.h"

// salmon includes
#include "FastxParser.hpp"
#include "SalmonConfig.hpp"
#include "SalmonDefaults.hpp"
#include "SalmonOpts.hpp"
#include "SalmonUtils.hpp"
#include "ProgramOptionsGenerator.hpp"

using paired_parser_qual = fastx_parser::FastxParser<fastx_parser::ReadQualPair>;
using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

/* ALEVIN DECLERATIONS*/
using bcEnd = BarcodeEnd;
namespace apt = alevin::protocols;
namespace aut = alevin::utils;

template <typename ProtocolT>
int alevin_sc_align(AlevinOpts<ProtocolT>& aopt,
                    SalmonOpts& sopt,
                    boost::program_options::parsed_options& orderedOptions);

template <typename ProtocolT>
int alevinQuant(AlevinOpts<ProtocolT>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);

//colors for progress monitoring
const char RESET_COLOR[] = "\x1b[0m";
char green[] = "\x1b[30m";
char red[] = "\x1b[30m";

/*
  Parse through mate1 file Rapidly and counts the density of each barcode
 */
template <typename ProtocolT>
void densityCalculator(single_parser* parser,
                       AlevinOpts<ProtocolT>& aopt,
                       CFreqMapT& freqCounter,
                       std::atomic<uint64_t>& usedNumBarcodes,
                       std::atomic<uint64_t>& totNumBarcodes){
  size_t rangeSize{0};
  uint32_t index;

  auto rg = parser->getReadGroup();
  auto log = aopt.jointLog;

  uint64_t usedNumBarcodesLocal{0};
  uint64_t totNumBarcodesLocal{0};
  std::string extractedBarcode( aopt.protocol.barcodeLength, 'N' );
  while (parser->refill(rg)) {
    rangeSize = rg.size();
    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      //Sexy Progress monitor
      ++totNumBarcodesLocal;
      if (not aopt.quiet and totNumBarcodesLocal % 500000 == 0) {
        aopt.jointLog->flush();
        fmt::print(stderr, "\r\r{}processed{} {} Million {}barcodes{}",
                   green, red, totNumBarcodesLocal/1000000, green, RESET_COLOR);
      }

      auto& rp = rg[i];
      std::string& seq = rp.seq;
      //if (aopt.protocol.end == bcEnd::THREE) {
      //  std::reverse(seq.begin(), seq.end());
      //}
      //extractedBarcode.clear();

      // NOTE: This means custom barcode extraction can't be on
      // on more than one read in this mode!  Come back and fix this later
      bool extraction_ok = aut::extractBarcode(seq, seq, aopt.protocol, extractedBarcode);
      bool seqOk = (extraction_ok) ? aut::sequenceCheck(extractedBarcode) : false;
      if (not seqOk){
        bool recovered = aut::recoverBarcode(extractedBarcode);
        if (not recovered) { continue; }
      }
      freqCounter[extractedBarcode] += 1;
      ++usedNumBarcodesLocal;
    }//end-for
  }//end-while
  totNumBarcodes += totNumBarcodesLocal;
  usedNumBarcodes += usedNumBarcodesLocal;
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

// reference from https://github.com/scipy/scipy/blob/master/scipy/stats/kde.py
// and https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/umi_methods.py#L193
bool gaussianKDE(const std::vector<uint32_t>& freqCounter,
                 const std::vector<size_t>& sortedIdx,
                 double& invCovariance, double& normFactor,
                 uint32_t expect_cells, uint32_t& predicted_cells){
  double covariance{0.0}, mean{0.0}, bw_method {0.01}, threshold = 0.001*freqCounter[sortedIdx[0]];
  std::vector<double> logDataset;
  size_t numElem, xSpace {10000};

  // extract counts above threshold
  for(size_t i=0; i<freqCounter.size(); i++){
    double count = static_cast<double>(freqCounter[ sortedIdx[ i ] ]);
    if (count <= threshold){
      break;
    }
    count = std::log10(count);
    mean += count;
    logDataset.emplace_back(count);
  }

  // size of the reference
  numElem = logDataset.size();

  // get mean of the data
  mean /= numElem;

  // generate the covariance
  for (auto count: logDataset){
    covariance += std::pow((count - mean), 2);
  }
  covariance = (covariance * bw_method) / (numElem - 1.0);

  if (covariance == 0){
    std::cout << "0 Covariance error for Gaussian kde"<<std::flush;
    exit(64);
  }

  const double PI = 3.1415926535897;
  invCovariance = 1.0 / covariance;
  normFactor = std::sqrt( 2.0*PI*covariance ) * numElem;

  // Evaluate Step starting from here
  double decrement = (logDataset[0] - logDataset[logDataset.size()-1]) / static_cast<double>(xSpace);
  std::vector<double> density(xSpace, 0.0);

  for (size_t i=0; i<numElem; i++){
    double pred = logDataset[0];
    for(size_t j=0; j<xSpace; j++,pred -= decrement){
      double diff = logDataset[ i ] - pred;
      double energy = ( std::pow(diff,2)*invCovariance ) / 2.0;
      density[j] += std::exp(-energy);
    }
  }

  //calculating the argrelextrema
  std::vector<uint32_t> local_mins;
  for (size_t i=1; i<xSpace-1; i++){
    if (density[i-1] > density[i] and density[i] < density[i+1]){
      local_mins.emplace_back(i);
    }
  }

  for (auto minIdx: local_mins){
    double freqThreshold = std::pow(10, logDataset[0]-(minIdx*decrement));
    size_t boundary {0};
    while( freqThreshold <= freqCounter[sortedIdx[boundary]] ){
      boundary ++;
    }
    if (boundary > expect_cells){
      predicted_cells = boundary;
      return false;
    }
    else if (expect_cells*0.1 > boundary){
      continue;
    }
    else{
      predicted_cells = boundary;
      return true;
    }
  }
  predicted_cells = 0;
  return false;
}

uint32_t getLeftBoundary(std::vector<size_t>& sortedIdx,
                         uint32_t topxBarcodes,
                         const std::vector<uint32_t>& freqCounter){
  // iterate in reverse order since sortedIdx is sorted in decreasing
  // order
  double cumCount{0.0};
  std::vector<double> freqs(topxBarcodes);
  for(uint32_t i = 0; i < topxBarcodes; i++){
    size_t ind = sortedIdx[topxBarcodes-i-1];
    cumCount += freqCounter[ind];
    freqs[i] = std::log(cumCount);
  }

  std::vector<uint32_t> X(topxBarcodes);
  std::iota(X.begin(), X.end(), 0);

  bool isUp;
  uint32_t x;
  double y, slope, leftExtreme{freqs[0]};
  for(uint32_t j=0; j<topxBarcodes; j++){
    x = X[j];
    y = freqs[j];

    if(y == leftExtreme){
      continue;
    }

    size_t nextBcIdx(j+1);
    std::vector<double> Y(topxBarcodes-nextBcIdx);
    slope = y/x;
    // fill in the values for fitted line
    std::transform(X.begin()+nextBcIdx, X.end(), Y.begin(),
                   [slope](uint32_t i) {return i*slope;} );

    isUp = false;
    double curveY, lineY;
    for(auto i=nextBcIdx; i<topxBarcodes; i++){
      curveY = freqs[i];
      lineY = Y[i-nextBcIdx];
      if (lineY > curveY){
        isUp = true;
        break;
      }
    }

    if(isUp == false){
      // ignoring all the frequencies having same frequency as cutoff
      uint32_t cutoff = topxBarcodes-j;
      uint32_t cutoffFrequency = freqCounter[sortedIdx[cutoff]];
      uint32_t nearestLeftFrequency = cutoffFrequency;
      while(nearestLeftFrequency == cutoffFrequency){
        nearestLeftFrequency = freqCounter[sortedIdx[--cutoff]];
      }
      return cutoff;
    }
  }

  return 0;
}


void dumpFeatures(std::string& frequencyFileName,
                  std::vector<size_t>& sortedIdx,
                  const std::vector<uint32_t>& freqCounter,
                  std::unordered_map<uint32_t, nonstd::basic_string_view<char>>& colMap
                  ){
  std::ofstream freqFile;
  freqFile.open(frequencyFileName);
  for (auto i:sortedIdx){
    uint32_t count = freqCounter[i];
    if (count == 0){
      break;
    }
    freqFile << colMap[i] << "\t"
             << count << "\n";
  }
  freqFile.close();
}

/*
  Knee calculation and sample true set of barcodes
 */
template <typename ProtocolT>
void sampleTrueBarcodes(const std::vector<uint32_t>& freqCounter,
                        TrueBcsT& trueBarcodes, size_t& lowRegionNumBarcodes,
                        std::unordered_map<uint32_t, nonstd::basic_string_view<char>> colMap,
                        AlevinOpts<ProtocolT>& aopt){
  std::vector<size_t> sortedIdx = sort_indexes(freqCounter);
  size_t maxNumBarcodes { aopt.maxNumBarcodes }, lowRegionMaxNumBarcodes { 1000 };
  size_t lowRegionMinNumBarcodes { aopt.lowRegionMinNumBarcodes };
  double lowConfidenceFraction { 0.5 };
  uint32_t topxBarcodes = std::min(maxNumBarcodes, freqCounter.size());
  uint64_t history { 0 };

  if (aopt.forceCells > 0) {
    topxBarcodes = aopt.forceCells - 1;
    size_t belowThresholdCb {0};
    while( freqCounter[sortedIdx[topxBarcodes]] < aopt.freqThreshold ) {
      --topxBarcodes;
      ++belowThresholdCb;
    }

    // converting 0 index to 1 index
    ++topxBarcodes;
    aopt.jointLog->info("Throwing {} barcodes with < {} reads", belowThresholdCb, aopt.freqThreshold);
    aopt.jointLog->flush();
  }
  else if (aopt.expectCells > 0){
    // Expect Cells algorithm is taken from
    // https://github.com/10XGenomics/cellranger/blob/e5396c6c444acec6af84caa7d3655dd33a162852/lib/python/cellranger/stats.py#L138

    constexpr uint32_t MAX_FILTERED_CELLS = 2;
    constexpr double UPPER_CELL_QUANTILE = 0.01;
    constexpr double FRACTION_MAX_FREQUENCY = 0.1;

    uint32_t baselineBcs = std::max(1, static_cast<int32_t>( aopt.expectCells*UPPER_CELL_QUANTILE ));
    double cutoffFrequency = std::max(1.0, freqCounter[sortedIdx[baselineBcs]] * FRACTION_MAX_FREQUENCY);
    uint32_t maxNumCells = std::min(static_cast<uint32_t>(freqCounter.size()), aopt.expectCells * MAX_FILTERED_CELLS);

    topxBarcodes = maxNumCells;
    for(size_t i=baselineBcs; i<maxNumCells; i++){
      if (freqCounter[sortedIdx[i]] < cutoffFrequency) {
        topxBarcodes = i + 1;
        break;
      }
    }

    if (topxBarcodes == maxNumCells){
      aopt.jointLog->error("Can't find right Boundary.\n"
                           "Please Report this issue on github.");
      aopt.jointLog->flush();
      exit(64);
    }
  } else{
    topxBarcodes = getLeftBoundary(sortedIdx,
                                   topxBarcodes,
                                   freqCounter);
    if (topxBarcodes == 0){
      aopt.jointLog->error("Can't find left Boundary.\n"
                           "Please Report this issue on github.");
      aopt.jointLog->flush();
      exit(64);
    }
    else{
      aopt.jointLog->info("Knee found left boundary at {} {} {}",
                          green, topxBarcodes, RESET_COLOR);

      double invCovariance {0.0}, normFactor{0.0};
      uint32_t gaussThreshold;
      bool isGaussOk = gaussianKDE(freqCounter, sortedIdx,
                                   invCovariance, normFactor,
                                   topxBarcodes, gaussThreshold);

      if ( isGaussOk ){
        topxBarcodes = gaussThreshold;
        // consider only if within 10% of current prediction
        aopt.jointLog->info("Gauss Corrected Boundary at {} {} {}",
                            green, gaussThreshold, RESET_COLOR);
      }
      else{
        aopt.jointLog->warn("Gauss Prediction {} Too far from knee prediction skipping it",
                            gaussThreshold);
      }

      aopt.jointLog->info("Learned InvCov: {} normfactor: {}",
                          invCovariance, normFactor);
      if (invCovariance == 0.0 or normFactor == 0.0){
        aopt.jointLog->error("Wrong invCovariance/Normfactor");
        aopt.jointLog->flush();
        exit(64);
      }
    }
  }//end-left-knee finding case

  uint32_t fractionTrueBarcodes = static_cast<int>(lowConfidenceFraction * topxBarcodes);

  if (fractionTrueBarcodes < lowRegionMinNumBarcodes){
    lowRegionNumBarcodes = lowRegionMinNumBarcodes;
  } else if (fractionTrueBarcodes > lowRegionMaxNumBarcodes){
    lowRegionNumBarcodes = lowRegionMaxNumBarcodes;
  } else {
    lowRegionNumBarcodes = fractionTrueBarcodes;
  }


  // ignoring all the frequencies having same frequency as cutoff
  // to imitate stable sort
  aopt.kneeCutoff = topxBarcodes;
  size_t totalUsableBarcodes = lowRegionNumBarcodes + topxBarcodes + 1;
  if (totalUsableBarcodes > freqCounter.size()) {
    size_t offset = freqCounter.size() - topxBarcodes - 1;
    lowRegionNumBarcodes = offset;
  }
  topxBarcodes += lowRegionNumBarcodes;

  uint32_t cutoffFrequency = freqCounter[sortedIdx[ topxBarcodes ]];
  uint32_t nearestLeftFrequency = cutoffFrequency;
  while(nearestLeftFrequency == cutoffFrequency && lowRegionNumBarcodes != 0){
    nearestLeftFrequency = freqCounter[sortedIdx[--topxBarcodes]];
    if (lowRegionNumBarcodes > lowRegionMinNumBarcodes)  { lowRegionNumBarcodes--; }
  }

  if ( lowRegionNumBarcodes != 0 ) {
    lowRegionNumBarcodes++;
    topxBarcodes++;
  }

  aopt.totalLowConfidenceCBs = lowRegionNumBarcodes;
  // keeping some cells left of the left boundary for learning
  aopt.jointLog->info("Total {}{}{}(has {}{}{} low confidence)"
                      " barcodes",
                      green, topxBarcodes+1, RESET_COLOR,
                      green, lowRegionNumBarcodes, RESET_COLOR);
  aopt.jointLog->flush();

  if (aopt.dumpfeatures) {
    auto frequencyFileName = (aopt.outputDirectory / "raw_cb_frequency.txt").string();
    dumpFeatures(frequencyFileName, sortedIdx, freqCounter, colMap);
  }

  for (size_t i=0; i<topxBarcodes; i++){
    trueBarcodes.insert(nonstd::to_string(colMap[sortedIdx[i]]));
  }
  aopt.numCells = trueBarcodes.size();
}

/*
  Index barcodes: map each possible 16M barcodes to one
  element of the set of true Barcode.
 */
template <typename ProtocolT>
void indexBarcodes(AlevinOpts<ProtocolT>& aopt,
                   CFreqMapT& freqCounter,
                   TrueBcsT& trueBarcodes,
                   SoftMapT& barcodeSoftMap){
  std::unordered_set<std::string> neighbors;
  std::unordered_map<std::string, std::vector<std::string>> ZMatrix;
  uint32_t wrngWhiteListCount{0};
  for (const auto trueBarcode: trueBarcodes){
    neighbors.clear();
    //find all neighbors to this true barcode
    aut::findNeighbors(aopt.protocol.barcodeLength,
                       trueBarcode,
                       neighbors);

    for(const auto& neighbor : neighbors){
      uint32_t freq{0};
      auto indexIt = freqCounter.find(neighbor);
      bool indexOk = indexIt != freqCounter.end();
      //bool indexOk = freqCounter.find(neighbor, freq);
      bool inTrueBc = trueBarcodes.find(neighbor) != trueBarcodes.end();
      if(not inTrueBc and indexOk and *indexIt > aopt.freqThreshold){
        ZMatrix[neighbor].push_back(trueBarcode);
      }
    }

    //bool inTrueBc = freqCounter.contains(trueBarcode);
    bool inTrueBc = freqCounter.find(trueBarcode) != freqCounter.end();
    if (not inTrueBc){
      wrngWhiteListCount += 1;
    }
  }//end-for
  //Done filling ZMatrix

  aopt.jointLog->info("Done populating Z matrix");
  if(trueBarcodes.size() - wrngWhiteListCount < 50){
    aopt.jointLog->warn("{} Whitelisted Barcodes with 0 frequency", wrngWhiteListCount);
  }

  std::string barcode;
  size_t numBarcodesCorrected {0};
  std::vector<std::pair<std::string, double>> dumpPair;

  for(auto& ZRow:ZMatrix){ //loop over every row of sparse matrix
    barcode = ZRow.first;
    //Avi -> have to work more on model
    dumpPair.clear();
    alevin::model::coinTossBarcodeModel(barcode,
                                        ZRow.second,
                                        dumpPair);
    numBarcodesCorrected += dumpPair.size();
    for(auto& updateVal: dumpPair){
      barcodeSoftMap[barcode].emplace_back(updateVal);
    }
  }//end-for

  aopt.jointLog->info("Total {} CB got sequence corrected", numBarcodesCorrected);
}

template <typename ProtocolT>
bool writeFastq(AlevinOpts<ProtocolT>& aopt,
                paired_parser_qual* parser,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes){
  size_t rangeSize{0};
  uint32_t totNumBarcodes{0};
  std::string barcode( aopt.protocol.barcodeLength, 'N');
  std::string umi( aopt.protocol.umiLength, 'N');

  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    std::random_device rd("/dev/urandom");
  #else
    std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 2020012

  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  try{
    auto rg = parser->getReadGroup();
    auto log = aopt.jointLog;

    fmt::print(stderr, "\n\n");
    while (parser->refill(rg)) {
      rangeSize = rg.size();
      for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
        auto& rp = rg[i];
        // barcode.clear();
        if(aopt.protocol.end == bcEnd::FIVE){
          {
            bool extracted_barcode = aut::extractBarcode(rp.first.seq, rp.second.seq, aopt.protocol, barcode);
            bool seqOk = (extracted_barcode) ? aut::sequenceCheck(barcode) : false;
            if (not seqOk){
              bool recovered = aut::recoverBarcode(barcode);
              if (not recovered) { continue; }
            }
          }

          {
            bool seqOk = aut::extractUMI(rp.first.seq, rp.second.seq, aopt.protocol, umi);
            if (not seqOk){
              std::cerr<<"Can't extract UMI";
              return false;
            }
          }
        }
        //else if (aopt.protocol.end == bcEnd::THREE) {
        //  std::string seq = rp.first.seq;
        //  std::reverse(seq.begin(), seq.end());
        //  barcode = rp.first.seq.substr(0, aopt.protocol.barcodeLength);
        //  umi = rp.first.seq.substr(aopt.protocol.barcodeLength,
        //                            aopt.protocol.umiLength);
        //}

        bool inTrueBc = trueBarcodes.find(barcode) != trueBarcodes.end();
        auto it = barcodeMap.find(barcode);
        bool inBarMap = it != barcodeMap.end();
        std::string corrBarcode;

        if(inTrueBc){
          corrBarcode = barcode;
        }
        else if(inBarMap and it->second.size() > 1){
          //toss a [0,1) real valued coin
          double rn = dist(mt);

          for (auto bcPair: it->second){
            if (rn < bcPair.second){
              corrBarcode = bcPair.first;
              break;
            }
          }
        }
        else if(inBarMap){
          corrBarcode = it->second.front().first;
        }
        else{
          continue;
        }

        std::cout << "@" << rp.second.name << "_" << corrBarcode << "_" << umi << "\n"
                  << rp.second.seq << "\n"
                  << "+" << "\n"
                  << rp.second.qual << "\n";

        totNumBarcodes += 1;
        if (totNumBarcodes % 500000 == 0) {
          char lred[] = "\x1b[30m";
          lred[3] = '0' + static_cast<char>(fmt::RED);
          fmt::print(stderr, "\r\r{}Dumped{} {} {}reads{}", green, lred,
                     totNumBarcodes, green, RESET_COLOR);
        }
      }//end-for
    }//end-while
    fmt::print(stderr, "\n");
    return true;
  }
  catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n" << std::flush;
    return false;
  }
}

/*
  function to rapidly parse through the barcode file, generate the density
  of each Unique barcode, use the knee method to select true barcodes and
  use our model to generate mapping of each 16M barcodes to true/null
  barcode.
 */
template <typename ProtocolT>
void processBarcodes(std::vector<std::string>& barcodeFiles,
                     std::vector<std::string>& readFiles,
                     AlevinOpts<ProtocolT>& aopt,
                     SoftMapT& barcodeSoftMap,
                     TrueBcsT& trueBarcodes,
                     CFreqMapT& freqCounter,
                     size_t& numLowConfidentBarcode){
  // Barcode processing always uses 2 threads now,
  // one parsing thread and one consumer thread.
  std::unique_ptr<single_parser> singleParserPtr{nullptr};
  constexpr uint32_t miniBatchSize{5000};
  uint32_t numParsingThreads{aopt.numParsingThreads}, numThreads{aopt.numConsumerThreads};
  std::vector<std::thread> threads;
  std::atomic<uint64_t> totNumBarcodes{0}, usedNumBarcodes{0};

  if (aopt.numThreads <= 3) {
    numThreads = 1;
  }

  //Populating Barcode Density Vector
  singleParserPtr.reset(new single_parser(barcodeFiles, numThreads,
                                          numParsingThreads, miniBatchSize));

  singleParserPtr->start();
  densityCalculator(singleParserPtr.get(), aopt,
                    freqCounter, usedNumBarcodes, totNumBarcodes);
  singleParserPtr->stop();

  fmt::print(stderr, "\n\n");
  aopt.jointLog->info("Done barcode density calculation.");
  aopt.jointLog->info("# Barcodes Used: {}{}{} / {}{}{}.",
                      green, usedNumBarcodes, RESET_COLOR,
                      red,totNumBarcodes, RESET_COLOR);

  aopt.totalReads = totNumBarcodes;
  aopt.totalUsedReads = usedNumBarcodes;

  //import whitelist barcodes if present
  if(boost::filesystem::exists(aopt.whitelistFile)){
    aut::readWhitelist(aopt.whitelistFile,
                       trueBarcodes);
    aopt.jointLog->info("Done importing white-list Barcodes");
    if (trueBarcodes.size() == 737280) {
      aopt.jointLog->error("Wrong whitelist provided\n"
                           "Please check https://salmon.readthedocs.io/en/develop/alevin.html#whitelist");
      aopt.jointLog->flush();
      exit(64);
    }

    std::vector<std::string> skippedTrueBarcodes ;
    for ( auto trueBarcode: trueBarcodes ) {
      auto it = freqCounter.find(trueBarcode);
      if (it == freqCounter.end() ) {
        skippedTrueBarcodes.emplace_back(trueBarcode);
      }
    }

    if ( skippedTrueBarcodes.size() > 0 ) {
      aopt.jointLog->warn("Skipping {} Barcodes as no read was mapped",
                          skippedTrueBarcodes.size());
      for (auto trueBarcode: skippedTrueBarcodes) {
        trueBarcodes.erase(trueBarcode);
      }
    }

    aopt.jointLog->info("Total {} white-listed Barcodes", trueBarcodes.size());

    if (aopt.dumpfeatures) {
      aopt.jointLog->info("Sorting and dumping raw barcodes");
      auto frequencyFileName = (aopt.outputDirectory / "raw_cb_frequency.txt").string();
      std::vector<uint32_t> collapsedfrequency;
      std::unordered_map<uint32_t, nonstd::basic_string_view<char>> collapMap;
      size_t ind{0};
      for(auto fqIt = freqCounter.begin(); fqIt != freqCounter.end(); ++fqIt){
        collapsedfrequency.push_back(fqIt.value());
        collapMap[ind] = fqIt.key_sv();
        ind += 1;
      }

      std::vector<size_t> sortedIdx = sort_indexes(collapsedfrequency);
      dumpFeatures(frequencyFileName, sortedIdx, collapsedfrequency, collapMap);
    }
  }
  else {
    std::vector<uint32_t> collapsedfrequency;
    std::unordered_map<uint32_t, nonstd::basic_string_view<char>> collapMap;
    size_t ind{0};
    for(auto fqIt = freqCounter.begin(); fqIt != freqCounter.end(); ++fqIt){
      collapsedfrequency.push_back(fqIt.value());
      collapMap[ind] = fqIt.key_sv();
      ind += 1;
    }

    if (aopt.keepCBFraction != 0.0) {
      aopt.forceCells = std::min(static_cast<uint32_t>(aopt.keepCBFraction * freqCounter.size()),
                                 aopt.maxNumBarcodes);
      aopt.jointLog->info("Forcing to use {} cells", aopt.forceCells);
      aopt.jointLog->flush();
    }

    //Calculate the knee using the frequency distribution
    //and get the true set of barcodes
    sampleTrueBarcodes(collapsedfrequency, trueBarcodes,
                       numLowConfidentBarcode, collapMap, aopt);

    aopt.jointLog->info("Done True Barcode Sampling");
  }

  uint64_t readsThrown{0};
  for(auto fqIt = freqCounter.begin(); fqIt != freqCounter.end(); ++fqIt){
    std::string cb_seq = std::string(fqIt.key_sv());
    if (trueBarcodes.find(cb_seq) == trueBarcodes.end() ) {
      readsThrown += fqIt.value();
    }
  }

  double percentThrown = readsThrown*100.0 / totNumBarcodes;
  if (percentThrown > 50.0) {
    aopt.jointLog->warn("Total {}% reads will be thrown away"
                        " because of noisy Cellular barcodes.",
                        percentThrown);
  }
  else{
    aopt.jointLog->info("Total {}% reads will be thrown away"
                        " because of noisy Cellular barcodes.",
                        percentThrown);
  }
  aopt.readsThrown = readsThrown;
  aopt.intelligentCutoff = trueBarcodes.size();

  indexBarcodes(aopt, freqCounter, trueBarcodes, barcodeSoftMap);
  aopt.jointLog->info("Done indexing Barcodes");

  aopt.jointLog->info("Total Unique barcodes found: {}", freqCounter.size());
  aopt.jointLog->info("Used Barcodes except Whitelist: {}", barcodeSoftMap.size());

  aopt.totalCBs = freqCounter.size();
  aopt.totalUsedCBs = barcodeSoftMap.size()+trueBarcodes.size();

  for(auto& softMapIt: barcodeSoftMap ){
    auto indexIt = freqCounter.find(softMapIt.first);
    bool indexOk = indexIt != freqCounter.end();
    if ( not indexOk){
      aopt.jointLog->error("Error: index not find in freq Counter\n"
                           "Please Report the issue on github");
      aopt.jointLog->flush();
      exit(64);
    }

    // NOTE: Need more clever way for ambiguity resolution.
    std::vector<std::pair<std::string, double>>& trBcVec = softMapIt.second;
    while(trBcVec.size() != 1){
      trBcVec.pop_back();
    }
    trBcVec.front().second = 1.0;

    std::string trBc = trBcVec.front().first;
    freqCounter[trBc] += *indexIt;
    freqCounter.erase(softMapIt.first);
  }

  if (aopt.dumpfq){
    std::unique_ptr<paired_parser_qual> pairedParserQualPtr{nullptr};
    pairedParserQualPtr.reset(new paired_parser_qual(barcodeFiles, readFiles,
                                                     1, 1, miniBatchSize));
    pairedParserQualPtr->start();
    bool isDumpok = writeFastq(aopt, pairedParserQualPtr.get(),
                               barcodeSoftMap, trueBarcodes);
    pairedParserQualPtr->stop();
    if(!isDumpok){
      aopt.jointLog->error("Not able to dump fastq."
                           "Something went wrong.\n"
                           "Please report this issue to github");
      aopt.jointLog->flush();
      std::exit(64);
    }
    aopt.jointLog->info("Done dumping fastq File");
  }
}

template <typename ProtocolT, typename OrderedOptionsT>
void initiatePipeline(AlevinOpts<ProtocolT>& aopt,
                      SalmonOpts& sopt,
                      OrderedOptionsT& orderedOptions,
                      boost::program_options::variables_map& vm,
                      std::string commentString, bool noTgMap,
                      std::vector<std::string> barcodeFiles,
                      std::vector<std::string> readFiles){
  bool isOptionsOk = aut::processAlevinOpts(aopt, sopt, noTgMap, vm);
  if (!isOptionsOk){
    aopt.jointLog->flush();
    exit(1);
  }

  spp::sparse_hash_map<uint32_t, uint32_t> txpToGeneMap;
  spp::sparse_hash_map<std::string, uint32_t> geneIdxMap;

  auto txpInfoFile = sopt.indexDirectory / "ctable.bin";
  if(!boost::filesystem::exists(txpInfoFile)){
    aopt.jointLog->error("Index Directory or the ctable.bin: {} doesn't exist.",
                         txpInfoFile.string());
    aopt.jointLog->flush();
    exit(1);
  }

  auto txpLengthFile = sopt.indexDirectory / "reflengths.bin";
  if(!boost::filesystem::exists(txpLengthFile)){
    aopt.jointLog->error("Index Directory or the reflengths.bin: {} doesn't exist.",
                         txpLengthFile.string());
    aopt.jointLog->flush();
    exit(1);
  }

  auto headerFile = sopt.indexDirectory / "info.json";
  if(!boost::filesystem::exists(headerFile)){
    aopt.jointLog->error("Index Directory or the info.json: {} doesn't exist.",
                         headerFile.string());
    aopt.jointLog->flush();
    exit(1);
  }

  aut::getTxpToGeneMap(txpToGeneMap, geneIdxMap,
                       aopt.geneMapFile.string(),
                       txpInfoFile.string(),
                       txpLengthFile.string(),
                       headerFile.string(),
                       aopt.jointLog, noTgMap);

  // If we're supposed to be quiet, set the global logger level to >= warn
  if (aopt.quiet) {
    spdlog::set_level(spdlog::level::warn); //Set global log level to warn
  }
  else {
    fmt::print(stderr, "{}\n\n", commentString);
  }

  if (aopt.just_align) {
    // if we are just aligning 
    auto rc = alevin_sc_align(aopt, sopt, orderedOptions);
    if (rc == 0) {
      aopt.jointLog->info("sc-align successful.");
      std::exit(0);
    } else {
      aopt.jointLog->error("sc-align exited with return code {}", rc);
      std::exit(rc);
    }
  }

  /*
    Barcode Knee generation
  */
  SoftMapT barcodeSoftMap;
  TrueBcsT trueBarcodes;
  //frequency counter
  CFreqMapT freqCounter;
  //freqCounter.set_max_resize_threads(sopt.maxHashResizeThreads);
  freqCounter.reserve(2097152);

  size_t numLowConfidentBarcode {0};

  if(boost::filesystem::exists(aopt.bfhFile)) {
    if (aopt.noQuant) {
      aopt.jointLog->error("Can't use --noQuant with hashMode (--hash)");
      aopt.jointLog->flush();
      exit(1);
    }

    salmonHashQuantify(aopt, sopt.outputDirectory, freqCounter);
    aopt.noQuant = true;

    aopt.jointLog->info("Done Processing");
  } else {
    aopt.jointLog->info("Processing barcodes files (if Present) \n\n ");
    processBarcodes(barcodeFiles,
                    readFiles,
                    aopt,
                    barcodeSoftMap,
                    trueBarcodes,
                    freqCounter,
                    numLowConfidentBarcode);
    aopt.jointLog->flush();
  }

  if(!aopt.noQuant){
    aopt.jointLog->info("Done with Barcode Processing; Moving to Quantify\n");
    alevinQuant(aopt, sopt, barcodeSoftMap, trueBarcodes,
                txpToGeneMap, geneIdxMap, orderedOptions,
                freqCounter, numLowConfidentBarcode);
  }
  else{
    boost::filesystem::path cmdInfoPath = vm["output"].as<std::string>();
    // Write out information about the command / run
    bool isWriteOk = aut::writeCmdInfo(cmdInfoPath / "cmd_info.json", orderedOptions);
    if(!isWriteOk){
      fmt::print(stderr, "writing in output directory failed\n Exiting Now");
      exit(1);
    }
  }
}

int salmonBarcoding(int argc, const char* argv[]) {
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  std::vector<std::string> barcodeFiles, readFiles, unmateFiles;
  bool optChain{false};
  int32_t numBiasSamples{0};

  double coverageThresh;
  //vector<string> unmatedReadFiles;
  //vector<string> mate1ReadFiles;
  //vector<string> mate2ReadFiles;

  SalmonOpts sopt;
  auto tot_cores = std::thread::hardware_concurrency();
  sopt.numThreads = std::max(1, static_cast<int>(tot_cores/4.0));

  salmon::ProgramOptionsGenerator pogen;

  auto inputOpt = pogen.getMappingInputOptions(sopt);
  auto mapSpecOpt = pogen.getMappingSpecificOptions(sopt);
  auto advancedOpt = pogen.getAdvancedOptions(numBiasSamples, sopt);
  auto hiddenOpt = pogen.getHiddenOptions(sopt);
  auto testingOpt = pogen.getTestingOptions(sopt);
  auto deprecatedOpt = pogen.getDeprecatedOptions(sopt);
  auto alevinBasicOpt = pogen.getAlevinBasicOptions(sopt);
  auto alevinDevsOpt = pogen.getAlevinDevsOptions();

  po::options_description all("alevin options");
  all.add(inputOpt).add(alevinBasicOpt).add(alevinDevsOpt).add(mapSpecOpt).add(advancedOpt).add(testingOpt).add(hiddenOpt).add(deprecatedOpt);

  po::options_description visible("alevin options");
  visible.add(inputOpt).add(alevinBasicOpt);

  po::variables_map vm;
  try {
    auto orderedOptions =
        po::command_line_parser(argc, argv).options(all).run();

    po::store(orderedOptions, vm);

    if (vm.count("help")) {
      auto hstring = R"(
alevin
==========
salmon-based processing of single-cell RNA-seq data.
)";

      std::cout << hstring << std::endl;
      std::cout << visible << std::endl;
      std::exit(0);
    }

    po::notify(vm);

    green[3] = '0' + static_cast<char>(fmt::GREEN);
    red[3] = '0' + static_cast<char>(fmt::RED);


    bool noTgMap {false};
    bool dropseq = vm["dropseq"].as<bool>();
    bool indrop = vm["indrop"].as<bool>();
    bool citeseq = vm["citeseq"].as<bool>();
    bool chromV3 = vm["chromiumV3"].as<bool>();
    bool chrom = vm["chromium"].as<bool>();
    bool gemcode = vm["gemcode"].as<bool>();
    bool celseq = vm["celseq"].as<bool>();
    bool celseq2 = vm["celseq2"].as<bool>();
    bool quartzseq2 = vm["quartzseq2"].as<bool>();
    bool custom_old =  vm.count("barcodeLength") and
                   vm.count("umiLength") and
                   vm.count("end");
    bool custom_new = (vm.count("bc-geometry") and 
                   vm.count("umi-geometry"));
    bool custom = custom_old or custom_new;

    uint8_t validate_num_protocols {0};
    if (dropseq) validate_num_protocols += 1;
    if (indrop) validate_num_protocols += 1;
    if (citeseq) { validate_num_protocols += 1; noTgMap = true;}
    if (chromV3) validate_num_protocols += 1;
    if (chrom) validate_num_protocols += 1;
    if (gemcode) validate_num_protocols += 1;
    if (celseq) validate_num_protocols += 1;
    if (celseq2) validate_num_protocols += 1;
    if (quartzseq2) validate_num_protocols += 1;
    if (custom and !noTgMap) validate_num_protocols += 1;

    if ( validate_num_protocols != 1 ) {
      fmt::print(stderr, "ERROR: Please specify one and only one scRNA protocol;");
      exit(1);
    }

    std::stringstream commentStream;
    commentStream << "### alevin (dscRNA-seq quantification) v" << salmon::version << "\n";
    commentStream << "### [ program ] => salmon \n";
    commentStream << "### [ command ] => alevin \n";
    for (auto& opt : orderedOptions.options) {
      commentStream << "### [ " << opt.string_key << " ] => {";
      for (auto& val : opt.value) {
        commentStream << " " << val;
      }
      commentStream << " }\n";
    }
    std::string commentString = commentStream.str();

    // Until we can figure out a better way to generify our parsing
    barcodeFiles = sopt.mate1ReadFiles;
    readFiles = sopt.mate2ReadFiles;
    unmateFiles = sopt.unmatedReadFiles;

    if (dropseq){
      AlevinOpts<apt::DropSeq> aopt;
      //aopt.jointLog->warn("Using DropSeq Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    }
    else if(indrop){
      std::cout<<"Indrop get neighbors removed, please use other protocols";
      exit(1);
      if(vm.count("w1") != 0){
        std::string w1 = vm["w1"].as<std::string>();
        AlevinOpts<apt::InDrop> aopt;
        aopt.protocol.setW1(w1);
        //aopt.jointLog->warn("Using InDrop Setting for Alevin");
        initiatePipeline(aopt, sopt, orderedOptions,
                         vm, commentString, noTgMap,
                         barcodeFiles, readFiles);
      }
      else{
        fmt::print(stderr, "ERROR: indrop needs w1 flag too.\n Exiting Now");
        exit(1);
      }
    }
    else if(citeseq){
      if(vm.count("featureStart") != 0 and vm.count("featureLength") != 0){
        AlevinOpts<apt::CITESeq> aopt;
        aopt.protocol.setFeatureLength(vm["featureLength"].as<size_t>());
        aopt.protocol.setFeatureStart(vm["featureStart"].as<size_t>());

        //aopt.jointLog->warn("Using InDrop Setting for Alevin");
        initiatePipeline(aopt, sopt, orderedOptions,
                         vm, commentString, noTgMap,
                         barcodeFiles, readFiles);
      }
      else{
        fmt::print(stderr, "ERROR: citeseq needs featureStart and featureLength flag too.\n Exiting Now");
        exit(1);
      }
    }
    else if(chromV3){
      AlevinOpts<apt::ChromiumV3> aopt;
      //aopt.jointLog->warn("Using 10x v3 Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    }
    else if(chrom){
      AlevinOpts<apt::Chromium> aopt;
      //aopt.jointLog->warn("Using 10x v2 Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    }
    else if(gemcode){
      AlevinOpts<apt::Gemcode> aopt;
      //aopt.jointLog->warn("Using 10x v1 Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       unmateFiles, readFiles);
    }
    else if(celseq){
      AlevinOpts<apt::CELSeq> aopt;
      //aopt.jointLog->warn("Using CEL-Seq Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    }
    else if(celseq2){
      AlevinOpts<apt::CELSeq2> aopt;
      //aopt.jointLog->warn("Using CEL-Seq2 Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    }
    else if(quartzseq2){
      AlevinOpts<apt::QuartzSeq2> aopt;
      //aopt.jointLog->warn("Using Quartz-Seq2 Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    } else if (custom_old) {
      AlevinOpts<apt::Custom> aopt;
      //aopt.jointLog->warn("Using Custom Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    } else if (custom_new) {
      AlevinOpts<apt::CustomGeometry> aopt;
      //aopt.jointLog->warn("Using Custom Setting for Alevin");
      initiatePipeline(aopt, sopt, orderedOptions,
                       vm, commentString, noTgMap,
                       barcodeFiles, readFiles);
    }

  } catch (po::error& e) {
    std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (const spdlog::spdlog_ex& ex) {
    std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n";
    std::cerr << argv[0] << " alevin was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0]
              << " alevin --help\nExiting.\n";
    std::exit(1);
  }

  return 0;
}
