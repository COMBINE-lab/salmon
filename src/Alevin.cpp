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

//alevin include
#include "Filter.hpp"
#include "AlevinOpts.hpp"
#include "AlevinUtils.hpp"
#include "BarcodeModel.hpp"
#include "SingleCellProtocols.hpp"
#include "BarcodeGroup.hpp"

// salmon includes
#include "FastxParser.hpp"
#include "SalmonConfig.hpp"

using paired_parser_qual = fastx_parser::FastxParser<fastx_parser::ReadQualPair>;
using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

constexpr uint32_t miniBatchSize{5000};

/* ALEVIN DECLERATIONS*/
using bcEnd = BarcodeEnd;
namespace apt = alevin::protocols;
namespace aut = alevin::utils;

template <typename ProtocolT>
int alevinQuant(AlevinOpts<ProtocolT>& aopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                int argc, char* argv[]);

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
                       std::mutex& ioMutex,
                       CFreqMapT& freqCounter,
                       std::atomic<uint64_t>& usedNumBarcodes,
                       std::atomic<uint64_t>& totNumBarcodes){
  size_t rangeSize{0};
  uint32_t index;
  std::string barcode;

  auto rg = parser->getReadGroup();
  auto log = aopt.jointLog;

  auto updatefn = [](uint32_t &num) { ++num; };

  while (parser->refill(rg)) {
    rangeSize = rg.size();
    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      //Sexy Progress monitor
      totNumBarcodes += 1;
      if (not aopt.quiet and totNumBarcodes % 500000 == 0) {
        ioMutex.lock();
        fmt::print(stderr, "\r\r{}processed{} {} Million {}barcodes{}",
                   green, red, totNumBarcodes/1000000, green, RESET_COLOR);
        ioMutex.unlock();
      }

      auto& rp = rg[i];
      std::string seq = rp.seq;
      if (aopt.protocol.end == bcEnd::THREE) {
        std::reverse(seq.begin(), seq.end());
      }
      bool isExtractOk = aut::extractBarcode(seq, aopt.protocol, barcode);
      if(!isExtractOk){
        continue;
      }

      bool seqOk = aut::sequenceCheck<ProtocolT>(barcode,
                                                 aopt,
                                                 ioMutex);
      if (not seqOk){
        continue;
      }
      freqCounter.upsert(barcode, updatefn, 1);
      usedNumBarcodes += 1;
    }//end-for
  }//end-while
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

//Avi -> Currently using Sircel's windowing method for knee
//calculation, lot's of predefined parameters.
//Should switch to a better way.
//[UPDATE:] Modified convergence criteria
//[UPDATE:] Not using this at all
//Returns number of top x barcodes in pdf to use as true Barcodes.
//uint32_t calculateKnee(std::vector<uint64_t> cdfDist){
//  //some pre-defined parameters
//  uint32_t PDF_THRESHOLD{100000}, LOCAL_WINDOW_LEN {25};
//
//  std::vector<uint32_t> xWindow(LOCAL_WINDOW_LEN);
//  std::vector<double> yWindow(LOCAL_WINDOW_LEN);
//  uint32_t noWindows = PDF_THRESHOLD - LOCAL_WINDOW_LEN + 1;
//  std::vector<double> slopes, filter;
//  double threshold;
//
//  //std::for_each(cdfDist.begin(),
//  //              cdfDist.begin() + PDF_THRESHOLD,
//  //              [](double &n){ n = log(n); } );
//
//  for (size_t i=0; i<noWindows; i++){
//    std::iota (std::begin(xWindow), std::end(xWindow), i);
//    yWindow.assign ( cdfDist.begin() + i,
//                     cdfDist.begin() + i + LOCAL_WINDOW_LEN);
//    auto slope = alevin::filter::fitLine( xWindow, yWindow );
//    slopes.push_back(-slope);
//  }
//
//  uint32_t SAVGOL_FRAME_SIZE{251}, SAVGOL_DIM{4};
//  size_t CONVERGENCE_WINDOW_RANGE{10}, CONVERGENCE_DIFF_THRESHOLD{5};
//  std::pair<uint32_t, uint32_t> SAVGOL_WINDOW{100, 5000};
//
//  alevin::filter::sgSmooth(slopes, SAVGOL_FRAME_SIZE, SAVGOL_DIM, filter);
//
//  //threshold = std::distance(filter.begin(),
//  //                          std::max_element(filter.begin() + SAVGOL_WINDOW.first,
//  //                                           filter.begin() + SAVGOL_WINDOW.second))
//  threshold = (SAVGOL_FRAME_SIZE-1)/2 + (LOCAL_WINDOW_LEN/2.0);
//
//  int32_t histVal{0}, histInd{0}, ind{SAVGOL_WINDOW.first}, sCount{0},diff;
//  for (auto it = filter.begin()+SAVGOL_WINDOW.first;
//       it < filter.begin()+SAVGOL_WINDOW.second;
//       it++){
//    ind += 1;
//    diff = std::abs(histVal - *it);
//    histVal = *it;
//    if (diff < CONVERGENCE_DIFF_THRESHOLD){
//      if (histInd == ind-1 or sCount == 0){
//        sCount += 1;
//        histInd = ind;
//      }
//      else{
//        sCount = 0;
//      }
//      if (sCount == CONVERGENCE_WINDOW_RANGE){
//        break;
//      }
//    }//end-diff-if
//  }//end-for
//  threshold += ind;
//
//  return static_cast<uint32_t>(threshold);
//}

uint32_t getLeftBoundary(std::vector<size_t>& sortedIdx,
                         uint32_t topxBarcodes,
                         const std::vector<uint32_t>& freqCounter){
  // iterate in reverse order since sortedIdx is sorted in decreasing
  // order
  double cumCount{0.0};
  std::vector<double> freqs(topxBarcodes);
  for(auto i=0; i<topxBarcodes; i++){
    size_t ind = sortedIdx[topxBarcodes-i];
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
    isUp = false;
    slope = y/x;
    // fill in the values for fitted line
    std::transform(X.begin()+nextBcIdx, X.end(), Y.begin(),
                   [slope](uint32_t i) {return i*slope;} );

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
      return topxBarcodes-j;
    }
  }

  return 0;
}

/*
  Knee calculation and sample true set of barcodes
 */
template <typename ProtocolT>
void sampleTrueBarcodes(const std::vector<uint32_t>& freqCounter,
                        TrueBcsT& trueBarcodes,
                        std::unordered_map<uint32_t, std::string> colMap,
                        AlevinOpts<ProtocolT>& aopt){
  std::vector<size_t> sortedIdx = sort_indexes(freqCounter);
  std::vector<uint64_t> cdfDist;
  size_t maxNumBarcodes = 100000;
  size_t lowRegionNumBarcodes = 1000;
  uint32_t topxBarcodes = std::min(maxNumBarcodes, freqCounter.size());
  uint64_t history{0};
  uint32_t threshold;

  topxBarcodes = getLeftBoundary(sortedIdx,
                                 topxBarcodes,
                                 freqCounter);
  if (topxBarcodes == 0){
    aopt.jointLog->error("Can't find left Boundary.\n"
                         "Please Report this issue on github.");
    exit(1);
  }
  else{
    aopt.jointLog->info("Found left boundary at {}{}{} / {}{}{}.",
                        green, maxNumBarcodes-topxBarcodes,
                        RESET_COLOR, red, maxNumBarcodes,
                        RESET_COLOR);

    // keeping 1000 cells left of the left boundary for learning
    aopt.jointLog->info("Total {}{}{}(+{}{}{} low confidence)"
                        " barcodes",
                        green, topxBarcodes, RESET_COLOR,
                        green, lowRegionNumBarcodes, RESET_COLOR);
    topxBarcodes += lowRegionNumBarcodes;
  }

  threshold = topxBarcodes;

  if(aopt.dumpfeatures){
    auto frequencyFileName = aopt.outputDirectory / "frequency.txt";
    std::ofstream freqFile;
    freqFile.open(frequencyFileName.string());
    for (auto i:sortedIdx){
      uint32_t count = freqCounter[i];
      if (topxBarcodes == 0 or count == 0){
        break;
      }
      freqFile << colMap[i] << "\t"
               << count << "\n";
      topxBarcodes--;
    }
    freqFile.close();
  }

  for (size_t i=0; i<threshold; i++){
    trueBarcodes.insert(colMap[sortedIdx[i]]);
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
      bool indexOk = freqCounter.find(neighbor, freq);
      bool inTrueBc = trueBarcodes.find(neighbor) != trueBarcodes.end();
      if(not inTrueBc and indexOk and freq > aopt.freqThreshold){
        ZMatrix[neighbor].push_back(trueBarcode);
      }
    }

    bool inTrueBc = freqCounter.contains(trueBarcode);
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
  std::vector<std::pair<std::string, double>> dumpPair;

  for(auto& ZRow:ZMatrix){ //loop over every row of sparse matrix
    barcode = ZRow.first;
    //Avi -> have to work more on model
    dumpPair.clear();
    alevin::model::coinTossBarcodeModel(barcode,
                                        aopt,
                                        ZRow.second,
                                        freqCounter,
                                        dumpPair);
    for(auto updateVal: dumpPair){
      barcodeSoftMap[barcode].push_back(updateVal);
    }
  }//end-for

  if(aopt.dumpBarcodeMap){
    auto dumpMapFile = aopt.outputDirectory / "barcodeSoftMaps.txt";
    std::ofstream mapFile;
    mapFile.open(dumpMapFile.string());

    for(const auto& softMapIt: barcodeSoftMap){
      auto trBcVec = softMapIt.second;
      mapFile << softMapIt.first << "\t" << trBcVec.size();
      for (auto trBc: trBcVec){
        mapFile << "\t" << trBc.first << "\t" << trBc.second;
      }
      mapFile << "\n";
    }
    mapFile.close();
  }

  if(aopt.dumpUmiToolsMap){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::unordered_map<std::string, std::vector<std::string>> umitoolsMap;

    for(auto trBc: trueBarcodes){
      umitoolsMap[trBc] = std::vector<std::string>();
    }

    auto umitoolsMapFile = aopt.outputDirectory / "umitoolsMap.txt";
    std::ofstream utFile;
    utFile.open(umitoolsMapFile.string());
    for(const auto& softMapIt: barcodeSoftMap){
      std::string trBc, bc{softMapIt.first};
      auto trBcVec = softMapIt.second;

      if(trBcVec.size() == 1){
        trBc = trBcVec.front().first;
      }
      else{
        double rn = dist(mt);
        for(auto dp: trBcVec){
          if(rn < dp.second){
            trBc = dp.first;
            break;
          }
        }
      }

      umitoolsMap[trBc].push_back(bc);
    }

    for (auto elem: umitoolsMap){
      auto trBc = elem.first;
      utFile << trBc << "\t";
      for (auto bc : elem.second){
        utFile << bc << ",";
      }
      utFile << "\b\n";
    }
    utFile.close();
  }
}

template <typename ProtocolT>
bool writeFastq(AlevinOpts<ProtocolT>& aopt,
                paired_parser_qual* parser,
                SoftMapT& barcodeMap,
                std::mutex& ioMutex,
                TrueBcsT& trueBarcodes){
  size_t rangeSize{0};
  uint32_t totNumBarcodes{0};
  std::string barcode;
  std::string umi;

  std::random_device rd;
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
        if(aopt.protocol.end == bcEnd::FIVE){
          barcode = rp.first.seq.substr(0, aopt.protocol.barcodeLength);
          umi = rp.first.seq.substr(aopt.protocol.barcodeLength,
                                    aopt.protocol.umiLength);
        }
        else if (aopt.protocol.end == bcEnd::THREE) {
          std::string seq = rp.first.seq;
          std::reverse(seq.begin(), seq.end());
          barcode = rp.first.seq.substr(0, aopt.protocol.barcodeLength);
          umi = rp.first.seq.substr(aopt.protocol.barcodeLength,
                                    aopt.protocol.umiLength);
        }

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
          char red[] = "\x1b[30m";
          red[3] = '0' + static_cast<char>(fmt::RED);
          fmt::print(stderr, "\r\r{}Dumped{} {} {}reads{}", green, red,
                     totNumBarcodes, green, RESET_COLOR);
        }
      }//end-for
    }//end-while
    fmt::print(stderr, "\n");
    return true;
  }
  catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n";
    return false;
  }
}

/*
  function to Rapidly parse through the barcode file, generate density
  of each Unique barcode, use knee method to select true barcodes and
  use our model to generate mapping of each 16M barcodes to true/null
  barcode.
 */
template <typename ProtocolT>
void processBarcodes(std::vector<std::string>& barcodeFiles,
                     std::vector<std::string>& readFiles,
                     AlevinOpts<ProtocolT>& aopt,
                     SoftMapT& barcodeSoftMap,
                     TrueBcsT& trueBarcodes){
  if (not aopt.nobarcode){
    //Avi -> HardCoding threads for Barcode Parsing
    //2 for consuming 1 for generating since
    //consumer thread is almost as fast as generator.
    std::unique_ptr<single_parser> singleParserPtr{nullptr};
    constexpr uint32_t miniBatchSize{5000};
    uint32_t numParsingThreads{aopt.numParsingThreads},
      numThreads{aopt.numConsumerThreads};
    std::vector<std::thread> threads;
    std::mutex ioMutex;
    std::atomic<uint64_t> totNumBarcodes{0}, usedNumBarcodes{0};

    //frequency counter
    CFreqMapT freqCounter;

    if (aopt.numThreads <= 3) {
      numThreads = 1;
    }

    //Populating Barcode Density Vector
    singleParserPtr.reset(new single_parser(barcodeFiles, numThreads,
                                            numParsingThreads, miniBatchSize));

    singleParserPtr->start();

    for (int i = 0; i < numThreads; ++i) {
      // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
      // change value before the lambda below is evaluated --- crazy!
      auto threadFun = [&, i]() -> void {
        densityCalculator(singleParserPtr.get(), aopt, ioMutex,
                          freqCounter, usedNumBarcodes, totNumBarcodes);
      };
      threads.emplace_back(threadFun);
    }

    for (int i = 0; i < numThreads; ++i) {
      threads[i].join();
    }

    singleParserPtr->stop();

    fmt::print(stderr, "\n\n");
    aopt.jointLog->info("Done barcode density calculation.");
    aopt.jointLog->info("# Barcodes Used: {}{}{} / {}{}{}.",
                        green, usedNumBarcodes, RESET_COLOR,
                        red,totNumBarcodes, RESET_COLOR);

    //import whitelist barcodes if present
    if(boost::filesystem::exists(aopt.whitelistFile)){
      std::ifstream whiteFile(aopt.whitelistFile.string());
      std::string whtBc;
      if(whiteFile.is_open()) {
        while(getline(whiteFile, whtBc)) {
          trueBarcodes.insert(whtBc);
        }
        whiteFile.close();
      }
      aopt.jointLog->info("Done importing white-list Barcodes");
      aopt.jointLog->info("Total {} white-listed Barcodes", trueBarcodes.size());
    }
    else {
      std::vector<uint32_t> collapsedfrequency;
      std::unordered_map<uint32_t, std::string> collapMap;
      size_t ind{0};
      for(const auto& fqIt : freqCounter.lock_table()){
        collapsedfrequency.push_back(fqIt.second);
        collapMap[ind] = fqIt.first;
        ind += 1;
      }

      //Calculate the knee using the frequency distribution
      //and get the true set of barcodes
      sampleTrueBarcodes(collapsedfrequency, trueBarcodes, collapMap, aopt);
      aopt.jointLog->info("Done True Barcode Sampling");
    }

    indexBarcodes(aopt, freqCounter, trueBarcodes, barcodeSoftMap);
    aopt.jointLog->info("Done indexing Barcodes");

    aopt.jointLog->info("Total Unique barcodes found: {}", freqCounter.size());
    aopt.jointLog->info("Used Barcodes except Whitelist: {}", barcodeSoftMap.size());

    uint32_t mmBcCounts{0}, mmBcReadCount{0};
    std::unordered_set<std::string> softMapWhiteBcSet;
    for(auto& softMapIt: barcodeSoftMap){
      std::vector<std::pair<std::string, double>>& trBcVec = softMapIt.second;
      if (trBcVec.size() > 1){
        mmBcCounts += 1;
        uint32_t numReads;
        bool indexOk = freqCounter.find(softMapIt.first, numReads);
        if ( not indexOk){
          aopt.jointLog->error("Error: index not find in freq Counter\n"
                               "Please Report the issue on github");
          exit(1);
        }
        for(std::pair<std::string, double> whtBc : trBcVec){
          softMapWhiteBcSet.insert(whtBc.first);
        }
        mmBcReadCount += numReads;
      }

      if (aopt.noSoftMap){
        while(trBcVec.size() != 1){
          trBcVec.pop_back();
        }
        trBcVec.front().second = 1.0;
      }
    }

    if (not aopt.noSoftMap){
      aopt.jointLog->info("Total Ambiguous Barcodes(soft-assigned):  {}", mmBcCounts);
      aopt.jointLog->info("Total CB-level Soft-Assignable Reads:  {}", mmBcReadCount);
      aopt.jointLog->info("Total whitelist-cells ambiguous reads can be assigned to: {}",
                          softMapWhiteBcSet.size());
      aopt.jointLog->info("Expected gain/cell using Alevin: {}", mmBcReadCount/softMapWhiteBcSet.size());
    }

    if (aopt.dumpfq){
      std::unique_ptr<paired_parser_qual> pairedParserQualPtr{nullptr};
      pairedParserQualPtr.reset(new paired_parser_qual(barcodeFiles, readFiles,
                                                       1, 1, miniBatchSize));
      pairedParserQualPtr->start();
      bool isDumpok = writeFastq(aopt, pairedParserQualPtr.get(),
                                 barcodeSoftMap, ioMutex, trueBarcodes);
      pairedParserQualPtr->stop();
      if(!isDumpok){
        aopt.jointLog->error("Not able to dump fastq."
                             "Something went wrong.\n"
                             "Please report this issue to github");
        aopt.jointLog->flush();
        std::exit(1);
      }
      aopt.jointLog->info("Done dumping fastq File");
    }
  }
  else{
    trueBarcodes.insert("AAA");
  }
}

template <typename ProtocolT, typename OrderedOptionsT>
void initiatePipeline(AlevinOpts<ProtocolT>& aopt,
                      OrderedOptionsT& orderedOptions,
                      boost::program_options::variables_map& vm,
                      std::string commentString,
                      std::vector<std::string> barcodeFiles,
                      std::vector<std::string> readFiles,
                      int argc, char* argv[]){
  bool isOptionsOk = aut::processAlevinOpts(aopt, vm);
  if (!isOptionsOk){
    exit(1);
  }

  // If we're supposed to be quiet, set the global logger level to >= warn
  if (aopt.quiet) {
    spdlog::set_level(spdlog::level::warn); //Set global log level to warn
  }
  else {
    fmt::print(stderr, "{}\n\n", commentString);
  }

  /*
    Barcode Knee generation
  */
  SoftMapT barcodeSoftMap;
  TrueBcsT trueBarcodes;
  aopt.jointLog->info("Processing barcodes files (if Present) \n\n ");

  processBarcodes(barcodeFiles,
                  readFiles,
                  aopt,
                  barcodeSoftMap,
                  trueBarcodes);

  aopt.jointLog->flush();
  spdlog::drop_all();

  if(!aopt.noQuant){
    aopt.jointLog->info("Done with Barcode Processing; Moving to Quantify\n");
    alevinQuant(aopt, barcodeSoftMap, trueBarcodes, argc, argv);
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

int salmonBarcoding(int argc, char* argv[]) {
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  std::vector<std::string> barcodeFiles, readFiles;

  green[3] = '0' + static_cast<char>(fmt::GREEN);
  red[3] = '0' + static_cast<char>(fmt::RED);

  po::options_description alevin_generic("\nAlevin Options");
  alevin_generic.add_options()
    ("version,v", "print version string")
    (
     "help,h", "produce help message")
    (
     "barcodeFiles,1", po::value<std::vector<std::string>>(&barcodeFiles)->multitoken()->required(),
     "List of files containing the barcodes+umi from drop seq dataset, usually the first file")
    (
     "readFiles,2", po::value<std::vector<std::string>>(&readFiles)->multitoken(),
     "List of files containing the barcodes+umi from drop seq dataset, usually the first file")
    (
     "output,o", po::value<std::string>()->required(),
     "Parent directory for dumping barcode mapping")
    (
     "dedup", po::bool_switch()->default_value(false),
     "Perform Directional per-cell deduplication")
    (
     "dropseq", po::bool_switch()->default_value(false),
     "Use DropSeq Single Cell protocol for the library")
    (
     "chromium", po::bool_switch()->default_value(false),
     "Use 10x chromium (under-development) Single Cell protocol for the library.")
    (
     "indrop", po::bool_switch()->default_value(false),
     "Use inDrop Single Cell protocol for the library. must specify w1 too.")
    (
     "w1", po::value<std::string>(),
     "Must be used in conjunction with inDrop;")
    (
     "index,i", po::value<std::string>(),
     "Salmon index")
    (
     "libType,l", po::value<std::string>(),
     "Format string describing the library type")
    (
     "whitelist", po::value<std::string>(),
     "File containing white-list barcodes")
    (
     "threads,p",
     po::value<uint32_t>(),
     "The number of threads to use concurrently.")
    (
     "dumpbarcodeeq", po::bool_switch()->default_value(false),
     "Dump JointEqClas with umi-barcode count.(Only DropSeq)")
    (
     "noquant", po::bool_switch()->default_value(false),
     "Don't run downstream barcode-Salmon model.")
    (
     "nosoftmap", po::bool_switch()->default_value(false),
     "Don't use soft-assignment for quant instead do hard-assignment.")
    (
     "dumpfq", po::bool_switch()->default_value(false),
     "Dump barcode modified fastq file for downstream analysis by"
     "using coin toss for multi-mapping.")
    (
     "dumpfeatures", po::bool_switch()->default_value(false),
     "Dump features for whitelist and downstream analysis.")
    (
     "dumpumitoolsmap", po::bool_switch()->default_value(false),
     "Dump umi_tools readable whitelist map for downstream analysis.")
    (
     "dumpbarcodemap", po::bool_switch()->default_value(false),
     "Dump BarcodeMap for downstream analysis.")
    (
     "quiet,q", po::bool_switch()->default_value(false),
     "Be quiet while doing quantification (don't write informative "
     "output to the console unless something goes wrong).")
    (
     "iupac,u",po::value<std::string>(),
     "<Deprecated>iupac code for cell-level barcodes.")
    (
     "end",po::value<uint32_t>(),
     "Cell-Barcodes end (5 or 3) location in the read sequence from where barcode has to"
     "be extracted. (end, umilength, barcodelength)"
     " should all be provided if using this option")
    (
     "umilength",po::value<uint32_t>(),
     "umi length Parameter for unknown protocol. (end, umilength, barcodelength)"
     " should all be provided if using this option")
    (
     "barcodelength",po::value<uint32_t>(),
     "umi length Parameter for unknown protocol. (end, umilength, barcodelength)"
     " should all be provided if using this option")
    (
     "noem",po::bool_switch()->default_value(false),
     "do not run em")
    (
     "nobarcode",po::bool_switch()->default_value(false),
     "this flag should be used when there is no barcode i.e. only one cell deduplication.")
    (
     "geneMap,g", po::value<std::string>(), "transcript to gene map tsv file")
    (
    "freqthreshold",po::value<uint32_t>(),
    "threshold for the frequency of the barcodes");


  po::variables_map vm;
  try {
    auto orderedOptions =
        po::command_line_parser(argc, argv).options(alevin_generic).run();

    po::store(orderedOptions, vm);

    if (vm.count("help")) {
      auto hstring = R"(
Alevin
==========
Salmon-based processing of single-cell RNA-seq data.
)";

      std::cerr << hstring << std::endl;
      std::cerr << alevin_generic << std::endl;
      std::exit(0);
    }

    po::notify(vm);

    bool dropseq = vm["dropseq"].as<bool>();
    bool indrop = vm["indrop"].as<bool>();
    bool chrom = vm["chromium"].as<bool>();

    if((dropseq and indrop) or
       (dropseq and chrom) or
       (chrom and indrop)){
      fmt::print(stderr, "ERROR: Please specify only one scRNA protocol;");
      exit(1);
    }

    std::stringstream commentStream;
    commentStream << "### salmon (single-cell-based) v" << salmon::version << "\n";
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

    if (dropseq){
      AlevinOpts<apt::DropSeq> aopt;
      initiatePipeline(aopt, orderedOptions,
                       vm, commentString,
                       barcodeFiles, readFiles,
                       argc, argv);
    }
    else if(indrop){
      std::cout<<"Indrop get neighbors removed, please use other protocols";
      exit(1);
      if(vm.count("w1") != 0){
        std::string w1 = vm["w1"].as<std::string>();
        AlevinOpts<apt::InDrop> aopt;
        aopt.protocol.setW1(w1);
        initiatePipeline(aopt, orderedOptions,
                         vm, commentString,
                         barcodeFiles, readFiles,
                         argc, argv);
      }
      else{
        fmt::print(stderr, "ERROR: indrop needs w1 flag too.\n Exiting Now");
        exit(1);
      }
    }
    else if(chrom){
      AlevinOpts<apt::Chromium> aopt;
      initiatePipeline(aopt, orderedOptions,
                       vm, commentString,
                       barcodeFiles, readFiles,
                       argc, argv);
    }
    else{
      AlevinOpts<apt::Custom> aopt;
      initiatePipeline(aopt, orderedOptions,
                       vm, commentString,
                       barcodeFiles, readFiles,
                       argc, argv);
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
