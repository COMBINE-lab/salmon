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
Graph graphFromCell(std::vector<TGroupT> txpGroups,
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

  Graph g;
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
      VertexType node = {
        static_cast<uint32_t>(eqId),
        static_cast<uint32_t>(uId)
      };
      VertexT v1 = add_vertex( Graph::vertex_property_type(node), g);

      for ( size_t uId_second=uId+1; uId_second<numUmis; uId_second++ ){
        VertexType node_second = {
          static_cast<uint32_t>(eqId),
          static_cast<uint32_t>(uId_second)
        };
        VertexT v2 = add_vertex(node_second, g);

        //check if two UMI can be connected
        EdgeType edge = hasEdge( umiSeqCounts[uId], umiSeqCounts[uId_second], ae );

        switch ( edge ) {
        case EdgeType::BiDirected:
          add_edge(v1, v2, g);
          add_edge(v2, v1, g);
          break;
        case EdgeType::XToY:
          add_edge(v1, v2, g);
          break;
        case EdgeType::YToX:
          add_edge(v2, v1, g);
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
          VertexType node = {
            static_cast<uint32_t>(eqId),
            static_cast<uint32_t>(uId)
          };
          VertexT v1 = add_vertex( Graph::vertex_property_type(node), g);

          for ( size_t uId_second=0; uId_second<num2Umis; uId_second++ ){
            VertexType node_second = {
              static_cast<uint32_t>(eq2Id),
              static_cast<uint32_t>(uId_second)
            };
            VertexT v2 = add_vertex( Graph::vertex_property_type(node_second), g);

            //check if two UMI can be connected
            EdgeType edge = hasEdge( umiSeqCounts[uId], umi2SeqCounts[uId_second], ae );

            switch ( edge ) {
            case EdgeType::BiDirected:
              add_edge(v1, v2, g);
              add_edge(v2, v1, g);
              break;
            case EdgeType::XToY:
              add_edge(v1, v2, g);
              break;
            case EdgeType::YToX:
              add_edge(v2, v1, g);
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

void collapseVertices(uint32_t vertex,
                      Graph& g,
                      std::vector<TGroupT>& txpGroups,
                      boost::property_map<Graph,boost::vertex_name_t>::type& vertName,
                      uint32_t& chosenTxp,
                      std::vector<uint32_t>& largestMcc) {
  VertexType vertexName = vertName[boost::vertex(vertex, g)];

  for (uint32_t txp: txpGroups[vertexName.eqclassId]){
    std::deque<uint32_t> bfsList;
    bfsList.push_back(vertex);

    spp::sparse_hash_set<uint32_t> visitedSet;
    visitedSet.insert(vertex);

    std::vector<uint32_t> currentMcc;
    while ( bfsList.size() != 0 ){
      uint32_t cv = bfsList.front();
      bfsList.pop_front();
      currentMcc.emplace_back(cv);

      typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = out_edges(vertex, g); ei != ei_end; ++ei) {
        auto source = boost::source ( *ei, g );
        auto nextVertex = boost::target ( *ei, g );

        if (visitedSet.contains(nextVertex)) {
          continue;
        }
        else{
          visitedSet.insert(nextVertex);
        }

        // extract transcripts from new vertex
        VertexType nvertexName = vertName[boost::vertex(nextVertex, g)];
        for (uint32_t ntxp: txpGroups[nvertexName.eqclassId]) {
          if (ntxp == txp){
            bfsList.emplace_back(nextVertex);
            break;
          }
        }//end-txp group for
      }//end-neighbors-for
    }//end-while

    if (largestMcc.size() < currentMcc.size()) {
      largestMcc = currentMcc;
      chosenTxp = txp;
    }
  } //end-for
}

void getNumMolecules(Graph& g,
                     std::vector<TGroupT>& txpGroups,
                     spp::sparse_hash_map<uint32_t, uint32_t>& t2gMap,
                     std::vector<SalmonEqClass>& salmonEqclasses){
  // get connected components
  std::vector<uint32_t> component( num_vertices(g) );
  uint32_t numComps = connected_components(g, component.data());
  spp::sparse_hash_map<std::vector<uint32_t>,
                       uint32_t,
                       container_hash<std::vector<uint32_t>>> eqclassHash;

  // making sets of relevant connected vertices
  std::vector<std::vector<uint32_t>> comps (numComps);
  for (size_t i=0; i<component.size(); i++) {
    comps[component[i]].emplace_back(static_cast<uint32_t>(i));
  }

  //property accessors
  boost::property_map<Graph, boost::vertex_name_t>::type vertName
    = boost::get(boost::vertex_name, g);

  // iterating over connected components
  for (auto& comp: comps) {
    // more than one vertex in the component
    if ( comp.size() > 1 ) {
      spp::sparse_hash_set<uint32_t> vset(comp.begin(), comp.end());

      while ( vset.size() != 0 ){
        std::vector<uint32_t> bestMcc;
        uint32_t bestCoveringTxp = std::numeric_limits<uint32_t>::max();
        for (uint32_t vertex: vset) {
          uint32_t coveringTxp;
          std::vector<uint32_t> newMcc;

          collapseVertices(vertex, g, txpGroups, vertName,
                           coveringTxp, newMcc);
          // choose the longer collapse: Greedy
          if (bestMcc.size() < newMcc.size()) {
            bestMcc = newMcc;
            bestCoveringTxp = coveringTxp;
          }
        }// end-vset for

        assert( bestCoveringTxp != std::numeric_limits<uint32_t>::max() );

        // get the gene id
        uint32_t bestCoveringGene = getGeneId(t2gMap, bestCoveringTxp);

        spp::sparse_hash_set<uint32_t> globalGenes ;
        for (size_t vId=0; vId<bestMcc.size(); vId++){
          uint32_t vertex = bestMcc[vId];
          spp::sparse_hash_set<uint32_t> localGenes;
          VertexType vertexName = vertName[boost::vertex(vertex, g)];

          for (uint32_t txp: txpGroups[vertexName.eqclassId]){
            uint32_t gId = getGeneId(t2gMap, txp);
            localGenes.insert(gId);
          }

          if (vId == 0) {
            globalGenes = localGenes;
          }
          else {
            spp::sparse_hash_set<uint32_t> intersect;
            std::set_intersection (globalGenes.begin(),
                                   globalGenes.end(),
                                   localGenes.begin(),
                                   localGenes.end(),
                                   std::inserter(intersect,
                                                 intersect.begin()));
            globalGenes = intersect;
          }
        }//end-mcc for

        assert(globalGenes.size() > 0);
        assert(globalGenes.contains(bestCoveringGene));

        for (auto rv: bestMcc){
          vset.erase(rv);
        }

        std::vector<uint32_t> genesVec (globalGenes.begin(),
                                        globalGenes.end());
        std::sort (genesVec.begin(), genesVec.end());
        eqclassHash[genesVec] += 1;
      }//end-while
    } // end-if comp.size()>1
    else{
      assert(comp.size() == 1);
      uint32_t vertex = comp[0];
      VertexType vertexName = vertName[boost::vertex(vertex, g)];
      TGroupT txps = txpGroups[vertexName.eqclassId];

      spp::sparse_hash_set<uint32_t> genes;
      for (auto txp: txps) {
        uint32_t gId = getGeneId(t2gMap, txp);
        genes.insert(gId);
      }

      assert(genes.size() > 0);

      std::vector<uint32_t> genesVec (genes.begin(), genes.end());
      std::sort (genesVec.begin(), genesVec.end());
      eqclassHash[genesVec] += 1;
    }//end-else comp.size()==1
  } //end-outer for comps iterator

  for (auto& it: eqclassHash) {
    SalmonEqClass eqclass = {
      it.first,
      it.second,
    };
    salmonEqclasses.emplace_back(eqclass);
  }
}

bool dedupClasses(std::vector<TGroupT> txpGroups,
                  std::vector<UGroupT> umiGroups,
                  spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap){
  // make directed graph from eqclasses
  Graph g = graphFromCell(txpGroups, umiGroups);
  std::vector<SalmonEqClass> salmonEqclasses;

  getNumMolecules(g, txpGroups, txpToGeneMap, salmonEqclasses);
  return true;
}
