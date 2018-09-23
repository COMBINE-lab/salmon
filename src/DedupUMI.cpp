#include "DedupUMI.hpp"

enum EdgeType {
  NoEdge,
  BiDirected,
  XToY,
  YToX,
};

EdgeType hasEdge(std::pair<std::string, uint32_t> &x,
                 std::pair<std::string, uint32_t> &y,
                 AlignerEngine& ae) {
  if ( x.first.compare(y.first) == 0 ) {
    return EdgeType::BiDirected;
  }
  if ( x.second > (2*y.second-1) ) {
    ae(x.first.c_str(), x.first.size(),
       y.first.c_str(), y.first.size(),
       edlibNewAlignConfig(1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));

    if ( ae.result().editDistance < 2 ){
      return EdgeType::XToY;
    }
    else {
      return EdgeType::NoEdge;
    }
  }
  else if (y.second > (2*x.second-1) ) {
    ae(x.first.c_str(), x.first.size(),
       y.first.c_str(), y.first.size(),
       edlibNewAlignConfig(1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));

    if ( ae.result().editDistance < 2 ){
      return EdgeType::YToX;
    }
    else {
      return EdgeType::NoEdge;
    }
  }
  else{
    return EdgeType::NoEdge;
  }
}


uint32_t getGeneId(spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap,
                   uint32_t tId ) {
  if( txpToGeneMap.contains(tId) ){
    return txpToGeneMap.at(tId);
  }
  else{
    std::cerr << "Out of Range error for txp to gene Map: "
              << '\n' << std::flush;
    std::cerr << tId << "\t not found" << std::flush;
    exit(1);
  }
}


typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::directedS> DirectedGraph;

void graphFromCell(std::vector<TGroupT> txpGroups,
                   std::vector<UGroupT> umiGroups,
                   spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap) {
  spp::sparse_hash_map<uint32_t, std::vector<uint32_t>> tidMap;
  std::vector<bool> multiGeneVec;

  // Get a map from each transcript to it's list of eq class
  size_t numClasses = txpGroups.size();
  for (size_t eqId=0; eqId<numClasses; eqId++) {
    spp::sparse_hash_set<uint32_t> gSet;
    for ( uint32_t txp : txpGroups[eqId] ) {
      tidMap[txp].emplace_back(eqId);
      uint32_t gId = getGeneId( txpToGeneMap, txp );
      gSet.insert(gId);
    }

    if ( gSet.size() > 1 ) {
      multiGeneVec.emplace_back(true);
    }
    else {
      multiGeneVec.emplace_back(false);
    }
  }

  //DirectedGraph g;
  //AlignerEngine ae;
  //for (size_t eqId=0; eqId<numClasses; eqId++) {
  //  size_t numUmis = umiGroups[eqId].size();
  //  for ( size_t uId=0; uId<numUmis; uId++ ){
  //    
  //  }
  //}//end-for
}

bool dedupClasses(std::vector<TGroupT> txpGroups,
                  std::vector<UGroupT> umiGroups,
                  spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap){
  graphFromCell(txpGroups, umiGroups, txpToGeneMap);
  return true;
}
