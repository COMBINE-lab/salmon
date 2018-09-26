#include "DedupUMI.hpp"

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
void graphFromCell(std::vector<TGroupT>& txpGroups,
                   std::vector<UGroupT>& umiGroups,
                   alevin::graph::Graph& g) {
  using namespace alevin::graph;
  spp::sparse_hash_map<uint32_t, std::vector<uint32_t>> tidMap;

  // Get a map from each transcript to it's list of eq class
  size_t numClasses = txpGroups.size();
  for (size_t eqId=0; eqId<numClasses; eqId++) {
    spp::sparse_hash_set<uint32_t> gSet;
    for ( uint32_t txp : txpGroups[eqId] ) {
      tidMap[txp].emplace_back(eqId);
    }
  }

  AlignerEngine ae;
  // alevin kmer object
  alevin::types::AlevinUMIKmer umiObj;
  spp::sparse_hash_map<VertexT, uint32_t, boost::hash<VertexT>> vertexIndexMap;

  //iterating over all eqclasses
  for (size_t eqId=0; eqId<numClasses; eqId++) {
    size_t numUmis = umiGroups[eqId].size();

    //// extracting umi sequences
    std::vector<std::pair<std::string, uint32_t>> umiSeqCounts;

    for(auto& it: umiGroups[eqId]) {
      umiObj.word__(0) = it.first;
      umiSeqCounts.emplace_back(std::make_pair(umiObj.toStr(), it.second));
    }

    for ( size_t uId=0; uId<numUmis; uId++ ){
      VertexT node (static_cast<uint32_t>(eqId), static_cast<uint32_t>(uId));
      uint32_t v1 = alevin::graph::getVertexIndex(vertexIndexMap, node);

      for ( size_t uId_second=uId+1; uId_second<numUmis; uId_second++ ){
        VertexT node_second (static_cast<uint32_t>(eqId), static_cast<uint32_t>(uId_second));
        uint32_t v2 = alevin::graph::getVertexIndex(vertexIndexMap, node_second);

        //check if two UMI can be connected
        EdgeType edge = alevin::graph::hasEdge( umiSeqCounts[uId], umiSeqCounts[uId_second], ae );

        switch ( edge ) {
        case EdgeType::BiDirected:
          g.add_edge(v1, v2);
          g.add_edge(v2, v1);
          break;
        case EdgeType::XToY:
          g.add_edge(v1, v2);
          break;
        case EdgeType::YToX:
          g.add_edge(v2, v1);
          break;
        case EdgeType::NoEdge:
          break;
        };
      }
    }//end-inner-for

    spp::sparse_hash_set<uint32_t> hSet;
    TGroupT& tgroup = txpGroups[eqId];
    size_t numTxps = tgroup.size();

    // iterate over all the transcripts
    for ( auto& txp: tgroup ) {
      for (uint32_t eq2Id: tidMap[txp]) {
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
          VertexT node (static_cast<uint32_t>(eqId), static_cast<uint32_t>(uId));
          uint32_t v1 = alevin::graph::getVertexIndex(vertexIndexMap, node);

          for ( size_t uId_second=0; uId_second<num2Umis; uId_second++ ){
            VertexT node_second (static_cast<uint32_t>(eq2Id), static_cast<uint32_t>(uId_second));
            if ( node == node_second ) {
              continue;
            }
            uint32_t v2 = alevin::graph::getVertexIndex(vertexIndexMap, node_second);

            //check if two UMI can be connected
            EdgeType edge = alevin::graph::hasEdge( umiSeqCounts[uId], umi2SeqCounts[uId_second], ae );

            switch ( edge ) {
            case EdgeType::BiDirected:
              g.add_edge(v1, v2);
              g.add_edge(v2, v1);
              break;
            case EdgeType::XToY:
              g.add_edge(v1, v2);
              break;
            case EdgeType::YToX:
              g.add_edge(v2, v1);
              break;
            case EdgeType::NoEdge:
              break;
            };
          } //end-for inner UMI
        }//end-for outerUMI
      }//end-for eq2Id
    }//end-inner for for txp
  }//end-outer-for

  size_t num_vertices = vertexIndexMap.size();
  g.vertexNames.resize(num_vertices);
  for (auto& it: vertexIndexMap) {
    g.vertexNames[it.second] = it.first;
  } // Done populating graph object
}

void collapseVertices(uint32_t vertex,
                      alevin::graph::Graph& g,
                      std::vector<TGroupT>& txpGroups,
                      uint32_t& chosenTxp,
                      std::vector<uint32_t>& largestMcc) {
  uint32_t eqclassId = g.getEqclassId(vertex);
  for (uint32_t txp: txpGroups[eqclassId]){
    std::deque<uint32_t> bfsList;
    bfsList.push_back(vertex);

    spp::sparse_hash_set<uint32_t> visitedSet;
    visitedSet.insert(vertex);

    std::vector<uint32_t> currentMcc;
    while ( bfsList.size() != 0 ){
      uint32_t cv = bfsList.front();
      bfsList.pop_front();
      currentMcc.emplace_back(cv);

      for (auto nextVertex: g.getNeighbors(vertex)) {
        if (visitedSet.contains(nextVertex)) {
          continue;
        }
        else{
          visitedSet.insert(nextVertex);
        }

        // extract transcripts from new vertex
        eqclassId = g.getEqclassId( nextVertex );
        for (uint32_t ntxp: txpGroups[eqclassId]) {
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

void getNumMolecules(alevin::graph::Graph& g,
                     std::vector<TGroupT>& txpGroups,
                     spp::sparse_hash_map<uint32_t, uint32_t>& t2gMap,
                     std::vector<SalmonEqClass>& salmonEqclasses){
  // get connected components
  std::vector<uint32_t> component;
  uint32_t numComps = g.connected_components(component);
  spp::sparse_hash_map<std::vector<uint32_t>,
                       uint32_t,
                       boost::hash<std::vector<uint32_t>>> eqclassHash;

  // making sets of relevant connected vertices
  std::vector<std::vector<uint32_t>> comps (numComps);
  for (size_t i=0; i<component.size(); i++) {
    comps[component[i]].emplace_back(static_cast<uint32_t>(i));
  }

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

          collapseVertices(vertex, g, txpGroups,
                           coveringTxp, newMcc);
          //choose the longer collapse: Greedy
          if (bestMcc.size() < newMcc.size()) {
            bestMcc = newMcc;
            bestCoveringTxp = coveringTxp;
          }
        }// end-vset for

        assert( bestCoveringTxp != std::numeric_limits<uint32_t>::max() );

        // get the gene id
        uint32_t bestCoveringGene = getGeneId(t2gMap, bestCoveringTxp);

        std::unordered_set<uint32_t> globalGenes ;
        for (size_t vId=0; vId<bestMcc.size(); vId++){
          uint32_t vertex = bestMcc[vId];
          std::unordered_set<uint32_t> localGenes;
          uint32_t eqclassId = g.getEqclassId(vertex);

          for (uint32_t txp: txpGroups[eqclassId]){
            uint32_t gId = getGeneId(t2gMap, txp);
            localGenes.insert(gId);
          }

          if (vId == 0) {
            globalGenes = localGenes;
          }
          else {
            std::unordered_set<uint32_t> intersect;
            unordered_set_intersection (globalGenes.begin(),
                                        globalGenes.end(),
                                        localGenes.begin(),
                                        localGenes.end(),
                                        std::inserter(intersect,
                                                      intersect.begin()));

            globalGenes = intersect;
          }
        }//end-mcc for

        if( globalGenes.size() == 0 ) {
          std::cerr << "can't find a representative gene for a molecule\n"
                    << "Please report this on gothub";
          exit(1);
        }

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
      uint32_t eqclassId = g.getEqclassId(vertex);
      TGroupT txps = txpGroups[eqclassId];

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

bool dedupClasses(std::vector<double>& geneAlphas,
                  uint64_t& totalUMICount,
                  std::vector<TGroupT>& txpGroups,
                  std::vector<UGroupT>& umiGroups,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap){
  // make directed graph from eqclasses
  alevin::graph::Graph g;
  graphFromCell(txpGroups, umiGroups, g);

  // make gene based eqclasses
  getNumMolecules(g, txpGroups, txpToGeneMap, salmonEqclasses);

  for( auto& eqclass: salmonEqclasses ) {
    if ( eqclass.labels.size() == 1 ) {
      totalUMICount += eqclass.count;
      geneAlphas[eqclass.labels.front()] += eqclass.count;
    }
    else if (eqclass.labels.size() == 0){
      std::cerr<<"Eqclasses with No gene labels\n";
      exit(1);
    }
  }

  return true;
}
