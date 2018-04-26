#include "Dedup.hpp"

// idea of deduplication of read in *each and every* eqclass is from
//Ntranos, Vasilis, et al. "Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts." Genome biology 17.1 (2016): 112
// we build upon the idea proposed on the above paper.

uint32_t edLibCollapse(const size_t& length,
                       const UGroupT& ugroup,
                       AlignerEngine& ae,
                       std::vector<std::pair<uint64_t, std::string>>& umiSeqs){
  uint32_t numCollapsed{0};

  uint64_t qUmiJellyIdx = umiSeqs.back().first;
  std::string qUmiStr = umiSeqs.back().second;

  uint64_t rUmiJellyIdx = umiSeqs.front().first;
  std::string rUmiStr = umiSeqs.front().second;

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

uint32_t neighborCollapse(const size_t& length,
                          const UGroupT& ugroup,
                          std::vector<std::pair<uint64_t, std::string>>& vList){
  uint32_t numCollapsed {0};
  std::deque<std::pair<uint64_t, std::string>> flagged_umis {vList.back()};
  vList.pop_back();

  while(flagged_umis.size() != 0){
    std::vector<uint32_t> collapseList;
    uint64_t qUmiJellyIdx = flagged_umis.front().first;
    std::string qUmiStr = flagged_umis.front().second;
    flagged_umis.pop_front();

    std::unordered_set<std::string> neighbors;
    alevin::utils::findNeighbors(length,
                                 qUmiStr,
                                 neighbors);

    // DON't change int to size_t: stupid wrapper error
    for (int32_t i = vList.size()-1; i>=0; i--){
      uint64_t rUmiJellyIdx = vList[ i ].first;
      std::string rUmiStr = vList[ i ].second;
      auto got = neighbors.find(rUmiStr);

      if (got != neighbors.end()){
        uint32_t qFreq, rFreq;
        auto rIt = ugroup.find(rUmiJellyIdx);
        auto qIt = ugroup.find(qUmiJellyIdx);

        if(rIt != ugroup.end() and qIt != ugroup.end()){
          rFreq = rIt->second;
          qFreq = qIt->second;
        }
        else{
          std::cerr<<"Wrong jelly inx found in the collapse stage."
                   <<"Please Report this on github \n" << std::flush;
          exit(1);
        }
        if ( (qFreq/2.0)+1 > rFreq ){
          numCollapsed += 1;
          collapseList.emplace_back(i);
          flagged_umis.push_back(vList[i]);
        }
      }
    }

    for (auto rIt = collapseList.begin();
         rIt!= collapseList.end(); ++rIt){
      vList.erase(vList.begin() + *rIt);
    }
  }

  return numCollapsed;
}


uint32_t dedupReads(
                    const size_t umiLength,
                    std::shared_ptr<spdlog::logger>& jointLog,
                    const UGroupT& ugroup,
                    spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint64_t>>& umiBiasList,
                    const TranscriptGroup& tgroup){

  //Aligner engine for edlib
  AlignerEngine ae;

  uint32_t bias{0};

  //lambda for extracting keys (umi)
  //auto key_selector = [](std::pair<uint64_t, uint32_t> pair){return pair.first;};
  spp::sparse_hash_set<uint64_t> umiList;

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
    // make a list of already deduplicated UMI
    spp::sparse_hash_set<uint64_t> seen_umis;
    for (auto txp : txps){
      if(umiBiasList.contains(txp)){
        for(auto umi: umiBiasList[txp]){
          seen_umis.insert(umi);
        }
      }
    }
    // extract indices to remove
    std::vector<uint64_t> removeUmis;
    for(auto umi: umiList){
      if (seen_umis.contains(umi)){
        removeUmis.emplace_back(umi);
      }
    }

    // remove the UMI from the list
    for (auto umi: removeUmis){
      umiList.erase(umi);
      bias += 1;
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
  std::vector<std::pair<uint64_t, std::string>> umis;
  for(auto& umi: visitList){
    alevin::kmer::AlvKmer jellyObj(umiLength);
    jellyObj.fromNum(umi);
    std::string umiseq = jellyObj.to_str();
    umis.emplace_back(std::make_pair(umi, umiseq));
    if (umiseq.size() > umiLength or umiseq.size() < umiLength){
        std::cout << "Size mismatch from Jelly Object\n"
            << "Expected :" << umiLength << "\tGot: " << umiseq.size()
            << std::flush;
        exit(1);
    }
  }

  if (umis.size() == 2){
    //keep counting until collapse everything
    auto numCollapsed = edLibCollapse(umiLength,
                                      ugroup,
                                      ae,
                                      umis);
    numMolecules -= numCollapsed;
    if (numMolecules < 1){
      jointLog->error("Collapse Procedure Error");
      exit(1);
    }
  }
  else{
    while(not umis.empty()){
      auto numCollapsed = neighborCollapse(umiLength,
                                           ugroup,
                                           umis);
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
