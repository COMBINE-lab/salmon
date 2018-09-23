#include "DedupUMI.hpp"


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


// choosing list for edges and vector for adjacency container
DirectedGraph graphFromCell(std::vector<TGroupT> txpGroups,
                            std::vector<UGroupT> umiGroups) {
  spp::sparse_hash_map<uint32_t, std::vector<uint32_t>> tidMap;

  // Get a map from each transcript to it's list of eq class
  size_t numClasses = txpGroups.size();
  for (size_t eqId=0; eqId<numClasses; eqId++) {
    spp::sparse_hash_set<uint32_t> gSet;
    for ( uint32_t txp : txpGroups[eqId] ) {
      tidMap[txp].emplace_back(eqId);
    }
  }

  DirectedGraph g;
  AlignerEngine ae;
  // alevin kmer object
  alevin::types::AlevinUMIKmer umiObj;

  //iterating over all eqclasses
  for (size_t eqId=0; eqId<numClasses; eqId++) {
    size_t numUmis = umiGroups[eqId].size();

    // extracting umi sequences
    std::vector<std::pair<std::string, uint32_t>> umiSeqCounts;

    for(auto& it: umiGroups[eqId]) {
      umiObj.word__(0) = it.first;
      umiSeqCounts.emplace_back(std::make_pair(umiObj.toStr(), it.second));
    }

    for ( size_t uId=0; uId<numUmis; uId++ ){
      std::pair<uint32_t, uint32_t> node (static_cast<uint32_t>(eqId),
                                          static_cast<uint32_t>(uId));
      VertexT v1 = add_vertex(node, g);

      for ( size_t uId_second=uId+1; uId_second<numUmis; uId_second++ ){
        std::pair<uint32_t, uint32_t> node_second (static_cast<uint32_t>(eqId),
                                                   static_cast<uint32_t>(uId_second));
        VertexT v2 = add_vertex(node_second, g);

        //check if two UMI can be connected
        EdgeType edge = hasEdge( umiSeqCounts[uId], umiSeqCounts[uId_second], ae );

        switch ( edge ) {
        case EdgeType::BiDirected:
          add_edge(v1, v2, EdgeType::BiDirected, g);
          add_edge(v2, v1, EdgeType::BiDirected, g);
          break;
        case EdgeType::XToY:
          add_edge(v1, v2, EdgeType::XToY, g);
          break;
        case EdgeType::YToX:
          add_edge(v2, v1, EdgeType::YToX, g);
          break;
        case EdgeType::NoEdge:
          break;
        };
      }
    }//end-inner-for

    spp::sparse_hash_set<uint32_t> hSet;
    // iterate over all the transcripts
    for ( auto& txpEqPair: tidMap ) {
      uint32_t txp = txpEqPair.first;
      for (uint32_t eq2Id: txpEqPair.second) {
        if (eq2Id < eqId) {
          continue;
        }

        if ( hSet.contains(eq2Id) ) {
          continue;
        }
        hSet.insert(eq2Id);

        size_t num2Umis = umiGroups[eq2Id].size();

        // extracting umi sequences
        std::vector<std::pair<std::string, uint32_t>> umi2SeqCounts;

        for(auto& it: umiGroups[eq2Id]) {
          umiObj.word__(0) = it.first;
          umi2SeqCounts.emplace_back(std::make_pair(umiObj.toStr(), it.second));
        }

        for ( size_t uId=0; uId<numUmis; uId++ ){
          std::pair<uint32_t, uint32_t> node (static_cast<uint32_t>(eqId),
                                              static_cast<uint32_t>(uId));
          VertexT v1 = add_vertex(node, g);

          for ( size_t uId_second=0; uId_second<num2Umis; uId_second++ ){
            std::pair<uint32_t, uint32_t> node_second (static_cast<uint32_t>(eq2Id),
                                                       static_cast<uint32_t>(uId_second));
            VertexT v2 = add_vertex(node_second, g);

            //check if two UMI can be connected
            EdgeType edge = hasEdge( umiSeqCounts[uId], umi2SeqCounts[uId_second], ae );

            switch ( edge ) {
            case EdgeType::BiDirected:
              add_edge(v1, v2, EdgeType::BiDirected, g);
              add_edge(v2, v1, EdgeType::BiDirected, g);
              break;
            case EdgeType::XToY:
              add_edge(v1, v2, EdgeType::XToY, g);
              break;
            case EdgeType::YToX:
              add_edge(v2, v1, EdgeType::YToX, g);
              break;
            case EdgeType::NoEdge:
              break;
            };
          } //end-for inner UMI
        }//end-for outerUMI
      }//end-for eq2Id
    }//end-inner for for txp
  }//end-outer-for

  return g;
}

bool dedupClasses(std::vector<TGroupT> txpGroups,
                  std::vector<UGroupT> umiGroups,
                  spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap){
  // make directed graph from eqclasses
  DirectedGraph g = graphFromCell(txpGroups, umiGroups, txpToGeneMap);
  return true;
}
