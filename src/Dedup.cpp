#include "Dedup.hpp"

// idea of deduplication of read in *each and every* eqclass is from
//Ntranos, Vasilis, et al. "Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts." Genome biology 17.1 (2016): 112
// we build upon the idea proposed on the above paper.

uint32_t edLibCollapse(const std::unordered_set<uint64_t>& umiList,
                       std::vector<uint64_t>& vList,
                       const size_t& length,
                       const UGroupT& ugroup,
                       AlignerEngine& ae,
                       std::vector<std::string>& umiSeqs){
  uint32_t numCollapsed{0};

  uint64_t qUmiJellyIdx = vList.back();
  std::string qUmiStr = umiSeqs.back();

  uint64_t rUmiJellyIdx = vList.front();
  std::string rUmiStr = umiSeqs.front();

  //convert back to char* pointer
  ae(qUmiStr.c_str(), length, rUmiStr.c_str(), length, edlibNewAlignConfig(1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
  auto distance = ae.result().editDistance;

  if (distance == 1){
    uint32_t rFreq  = ugroup.find(rUmiJellyIdx)->second;
    uint32_t qFreq = ugroup.find(qUmiJellyIdx)->second;
    if ( (qFreq/2.0)+1 > rFreq ){
      numCollapsed += 1;
    }
  }

  return numCollapsed;
}

uint32_t neighborCollapse(const std::unordered_set<uint64_t>& umiList,
                          std::vector<uint64_t>& vList,
                          const size_t& length,
                          const UGroupT& ugroup,
                          std::vector<std::string>& umiSeqs){
  //get index/string of current most frequent UMI
  //size_t qUmiIdx {vList[0]};
  //vList.erase(vList.begin());

  uint64_t qUmiJellyIdx = vList.back();
  std::string qUmiStr = umiSeqs.back();
  vList.pop_back();
  umiSeqs.pop_back();

  uint32_t numCollapsed {0};
  std::vector<uint32_t> collapseList;

  std::unordered_set<std::string> neighbors;
  alevin::utils::findNeighbors(length,
                               qUmiStr,
                               neighbors);

  for (int32_t i=vList.size()-1; i>=0; i--){
    //size_t rUmiIdx {vList[i]};
    uint64_t rUmiJellyIdx = vList[ i ];
    std::string rUmiStr = umiSeqs[ i ];

    auto got = neighbors.find(rUmiStr);

    if (got != neighbors.end()){
      uint32_t qFreq{1}, rFreq{1};
      auto rIt = ugroup.find(rUmiJellyIdx);
      auto qIt = ugroup.find(qUmiJellyIdx);

      if(rIt != ugroup.end() and qIt != ugroup.end()){
        rFreq = rIt->second;
        qFreq = qIt->second;
      }
      else{
        std::cerr<<"Wrong jelly inx found in the collapse stage."
                 <<"Please Report this on github \n";
        exit(1);
      }
      if ( (qFreq/2.0)+1 > rFreq ){
        numCollapsed += 1;
        collapseList.emplace_back(i);
      }
    }
  }

  for (auto rIt = collapseList.begin();
       rIt!= collapseList.end(); ++rIt){
    vList.erase(vList.begin() + *rIt);
    umiSeqs.erase(umiSeqs.begin() + *rIt);
  }

  return numCollapsed;
}


uint32_t dedupReads(
                    const size_t umiLength,
                    std::shared_ptr<spdlog::logger>& jointLog,
                    const UGroupT& ugroup,
                    std::unordered_map<uint32_t,
                    std::unordered_set<uint64_t>>& umiBiasList,
                    const TranscriptGroup& tgroup){

  //Aligner engine for edlib
  AlignerEngine ae;

  uint32_t bias{0};

  //lambda for extracting keys (umi)
  //auto key_selector = [](std::pair<uint64_t, uint32_t> pair){return pair.first;};
  std::unordered_set<uint64_t> umiList;

  //make a vector of umis
  for(auto& kpair: ugroup){
    umiList.insert(kpair.first);
  }
  //transform(ugroup.begin(), ugroup.end(), umiList.begin(), key_selector);

  //keep the umi for unique eqclasses otherwise remove from the list
  auto& txps  = tgroup.txps;
  if (txps.size()==1){
    for(auto umi: umiList){
      umiBiasList[txps[0]].insert(umi);
    }
  }
  else{
    for (auto txp : txps){
      auto got = umiBiasList.find(txp);
      if ( got != umiBiasList.end() ){
        auto biasUmis = got->second;
        for (auto umi : biasUmis){
          //auto it = std::find(umiList.begin(),
          //                    umiList.end(),
          //                    umi);
          auto it = umiList.find(umi);
          if(it != umiList.end()){
            umiList.erase(it);
            bias += 1;
          }
        }
      }
    }
  }

  //check for single UMI class
  if (umiList.size() == 1){
    return 1;
  }
  else if(umiList.size() == 0){
    return 0;
  }

  // A list maintaining visited nodes
  std::vector<uint64_t> visitList (umiList.begin(), umiList.end());

  //sort the order of umis by increasing frequency
  sort(visitList.begin(), visitList.end(),
       [&ugroup](uint64_t i1, uint64_t i2)
       {return ugroup.find(i1)->second < ugroup.find(i2)->second;});

  //  iota(visitList.begin(), visitList.end(), 0);
  uint32_t numMolecules {static_cast<uint32_t>(umiList.size())};

  // making a vector umi sequences
  std::vector<std::string> umiSeqs;
  alevin::kmer::AlvKmer jellyObj(umiLength);
  for(auto& umi: visitList){
    jellyObj.fromNum(umi);
    umiSeqs.emplace_back(jellyObj.to_str());
  }

  if (umiList.size() == 2){
    //keep counting until collapse everything
    auto numCollapsed = edLibCollapse(umiList,
                                      visitList,
                                      umiLength,
                                      ugroup,
                                      ae,
                                      umiSeqs);
    numMolecules -= numCollapsed;
    if (numMolecules < 1){
      jointLog->error("Collapse Procedure Error");
      exit(1);
    }
  }
  else{
    while(not visitList.empty()){
      auto numCollapsed = neighborCollapse(umiList,
                                           visitList,
                                           umiLength,
                                           ugroup,
                                           umiSeqs);
      numMolecules -= numCollapsed;
      if (numMolecules < 1){
        jointLog->error("Collapse Procedure Error");
        exit(1);
      }
    }
  }

  //if(!quiet){
  //  jointLog->info("# biased UMI count {}", bias);
  //  jointLog->info("# unique eqclasses {}", umiBiasList.size());
  //}

  // Write the main results
  //GZipWriter gzw(outDir, aopt.jointLog);
  //gzw.writeAbundances(sopt, experiment);


  //returns the number of molecules
  return numMolecules;
}
