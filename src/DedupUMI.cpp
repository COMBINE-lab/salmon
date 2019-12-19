#include "DedupUMI.hpp"
#include "tsl/hopscotch_map.h"

uint32_t getGeneId(spp::sparse_hash_map<uint32_t, uint32_t> &txpToGeneMap,
                   uint32_t tId ) {
  if( txpToGeneMap.contains(tId) ){
    return txpToGeneMap.at(tId);
  }
  else{
    std::cerr << "\n\n\nOut of Range error for txp to gene Map: "
              << tId << "\t not found" << std::flush;
    exit(74);
  }
}


// choosing list for edges and vector for adjacency container
void graphFromCell(std::vector<TGroupT>& txpGroups,
                   std::vector<UGroupT>& umiGroups,
                   alevin::graph::Graph& g,
                   uint32_t umiEditDistance, uint32_t umiLength,
                   std::atomic<uint64_t>& totalUniEdgesCounts,
                   std::atomic<uint64_t>& totalBiEdgesCounts) {
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

  // alevin kmer object
  alevin::types::AlevinUMIKmer umiObj;
  spp::sparse_hash_map<VertexT, uint32_t, boost::hash<VertexT>> vertexIndexMap;
  vertexIndexMap.reserve(numClasses);

  //iterating over all eqclasses
  for (size_t eqId=0; eqId<numClasses; eqId++) {
    size_t numUmis = umiGroups[eqId].size();

    //// extracting umi sequences
    std::vector<std::pair<uint64_t, uint32_t>> umiSeqCounts;

    for(auto& it: umiGroups[eqId]) {
      umiSeqCounts.emplace_back(std::make_pair(it.first, it.second));
    }

    for ( size_t uId=0; uId<numUmis; uId++ ){
      VertexT node (static_cast<uint32_t>(eqId), static_cast<uint32_t>(uId));
      uint32_t v1 = alevin::graph::getVertexIndex(vertexIndexMap, node);

      for ( size_t uId_second=uId+1; uId_second<numUmis; uId_second++ ){

        //check if two UMI can be connected
        EdgeType edge = alevin::graph::hasEdge( umiSeqCounts[uId], umiSeqCounts[uId_second], umiEditDistance, umiLength );
        if ( edge == EdgeType::NoEdge ) { continue; }

        VertexT node_second (static_cast<uint32_t>(eqId), static_cast<uint32_t>(uId_second));
        uint32_t v2 = alevin::graph::getVertexIndex(vertexIndexMap, node_second);

        switch ( edge ) {
        case EdgeType::BiDirected:
          totalBiEdgesCounts += 1;
          g.add_edge(v1, v2);
          g.add_edge(v2, v1);
          break;
        case EdgeType::XToY:
          totalUniEdgesCounts += 1;
          g.add_edge(v1, v2);
          break;
        case EdgeType::YToX:
          totalUniEdgesCounts += 1;
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
        if (eq2Id <= eqId) {
          continue;
        }

        if ( hSet.contains(eq2Id) ) {
          continue;
        }
        hSet.insert(eq2Id);

        size_t num2Umis = umiGroups[eq2Id].size();

        // extracting umi sequences
        std::vector<std::pair<uint64_t, uint32_t>> umi2SeqCounts;

        for(auto& it: umiGroups[eq2Id]) {
          umi2SeqCounts.emplace_back(std::make_pair(it.first, it.second));
        }

        for ( size_t uId=0; uId<numUmis; uId++ ){
          VertexT node (static_cast<uint32_t>(eqId), static_cast<uint32_t>(uId));
          uint32_t v1 = alevin::graph::getVertexIndex(vertexIndexMap, node);

          for ( size_t uId_second=0; uId_second<num2Umis; uId_second++ ){
            //check if two UMI can be connected
            EdgeType edge = alevin::graph::hasEdge( umiSeqCounts[uId], umi2SeqCounts[uId_second], umiEditDistance, umiLength );
            VertexT node_second (static_cast<uint32_t>(eq2Id), static_cast<uint32_t>(uId_second));
            if ( node == node_second or edge == EdgeType::NoEdge ) {
              continue;
            }

            uint32_t v2 = alevin::graph::getVertexIndex(vertexIndexMap, node_second);

            switch ( edge ) {
            case EdgeType::BiDirected:
              totalBiEdgesCounts += 1;
              g.add_edge(v1, v2);
              g.add_edge(v2, v1);
              break;
            case EdgeType::XToY:
              totalUniEdgesCounts += 1;
              g.add_edge(v1, v2);
              break;
            case EdgeType::YToX:
              totalUniEdgesCounts += 1;
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
                      std::vector<uint32_t>& largestMcc,
                      spp::sparse_hash_set<uint32_t>& processedSet) {
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
        if (visitedSet.contains(nextVertex) || !processedSet.contains(nextVertex)) {
          continue;
        }
        else{
          visitedSet.insert(nextVertex);
        }

        // extract transcripts from new vertex
        uint32_t neqclassId = g.getEqclassId( nextVertex );
        for (uint32_t ntxp: txpGroups[neqclassId]) {
          if (ntxp == txp){
            bfsList.emplace_back(nextVertex);
            break;
          }
        }//end-txp group for
      }//end-neighbors-for
    }//end-while

    if (largestMcc.size() < currentMcc.size()) {
      largestMcc = currentMcc;
    }
  } //end-for
}

void getNumMoleculesWithArborescence(alevin::graph::Graph& g,
                                     std::vector<TGroupT>& txpGroups,
                                     std::vector<UGroupT>& umiGroups,
                                     spp::sparse_hash_map<uint32_t, uint32_t>& t2gMap,
                                     std::vector<spp::sparse_hash_map<uint16_t, uint16_t>>& arboEqClassCount,
                                     std::vector<SalmonEqClass>& salmonEqclasses){
  // get connected components
  std::vector<uint32_t> component;
  uint32_t numComps = g.connected_components(component);
  spp::sparse_hash_map<std::vector<uint32_t>,
                       std::pair<uint32_t, spp::sparse_hash_map<uint16_t, uint16_t>>,
                       boost::hash<std::vector<uint32_t>>> eqclassHash;


  // making sets of relevant connected vertices
  std::vector<std::vector<uint32_t>> comps (numComps);
  for (size_t i=0; i<component.size(); i++) {
    comps[component[i]].emplace_back(static_cast<uint32_t>(i));
  }

  spp::sparse_hash_map<size_t, uint32_t> compCounter;
  // iterating over connected components
  for (auto& comp: comps) {
    compCounter[comp.size()] += 1;

    // more than one vertex in the component
    if ( comp.size() > 1 ) {
      spp::sparse_hash_set<uint32_t> vset(comp.begin(), comp.end());

      while ( vset.size() != 0 ){
        std::vector<uint32_t> bestMcc;
        for (uint32_t vertex: vset) {
          std::vector<uint32_t> newMcc;

          collapseVertices(vertex, g, txpGroups,
                           newMcc,
                           vset);
          //choose the longer collapse: Greedy
          if (bestMcc.size() < newMcc.size()) {
            bestMcc = newMcc;
          }
        }// end-vset for

        tsl::hopscotch_map<uint32_t, uint32_t> globalTxpCounts;
        for (size_t vId=0; vId<bestMcc.size(); vId++){
          uint32_t vertex = bestMcc[vId];
          uint32_t eqclassId = g.getEqclassId(vertex);

          for (uint32_t txp: txpGroups[eqclassId]){
            globalTxpCounts[txp] += 1;
          }
        }//end-mcc for

        // only transcripts that occur in *every* vertex, and hence
        // have a count of bestMcc.size(), are in the proper intersection.
        uint32_t requiredCount = bestMcc.size();
        std::vector<uint32_t> globalTxps;
        globalTxps.reserve(globalTxpCounts.size());
        for (auto kv : globalTxpCounts) {
          if (kv.second == requiredCount) { globalTxps.push_back(kv.first); }
        }

        if( globalTxps.size() == 0 ) {
          std::cerr << "can't find a representative transcript for a molecule\n"
                    << "Please report this on github ";
          exit(74);
        }

        uint16_t readspmol {0};
        for (auto rv: bestMcc){
          uint32_t eId = g.getEqclassId(rv);
          uint32_t uId = g.getUId(rv);
          auto it = umiGroups[eId].begin();
          std::advance(it, uId);
          readspmol += it->second;
          vset.erase(rv);
        }

        spp::sparse_hash_set<uint32_t> globalGenes;
        for(auto txp: globalTxps){
          uint32_t geneId = getGeneId(t2gMap, txp);
          globalGenes.insert(geneId);
        }

        std::vector<uint32_t> genesVec (globalGenes.begin(),
                                        globalGenes.end());
        std::sort (genesVec.begin(), genesVec.end());
        eqclassHash[genesVec].first += 1;
        eqclassHash[genesVec].second[readspmol] += 1;
      }//end-while
    } // end-if comp.size()>1
    else{
      assert(comp.size() == 1);
      uint32_t vertex = comp[0];
      uint32_t eId = g.getEqclassId(vertex);
      uint32_t uId = g.getUId(vertex);
      auto it = umiGroups[eId].begin();
      std::advance(it, uId);
      uint16_t readspmol = it->second;

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
      eqclassHash[genesVec].first += 1;
      eqclassHash[genesVec].second[readspmol] += 1;
    }//end-else comp.size()==1
  } //end-outer for comps iterator

  for (auto& it: eqclassHash) {
    SalmonEqClass eqclass = {
      it.first,
      it.second.first,
    };

    salmonEqclasses.emplace_back(eqclass);
    arboEqClassCount.emplace_back(it.second.second);
  }
}

void getNumMolecules(alevin::graph::Graph& g,
                     std::vector<TGroupT>& txpGroups,
                     std::vector<UGroupT>& umiGroups,
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

  spp::sparse_hash_map<size_t, uint32_t> compCounter;
  // iterating over connected components
  for (auto& comp: comps) {
    compCounter[comp.size()] += 1;

    // more than one vertex in the component
    if ( comp.size() > 1 ) {
      spp::sparse_hash_set<uint32_t> vset(comp.begin(), comp.end());

      while ( vset.size() != 0 ){
        std::vector<uint32_t> bestMcc;
        for (uint32_t vertex: vset) {
          std::vector<uint32_t> newMcc;

          collapseVertices(vertex, g, txpGroups,
                           newMcc,
                           vset);
          //choose the longer collapse: Greedy
          if (bestMcc.size() < newMcc.size()) {
            bestMcc = newMcc;
          }
        }// end-vset for

        tsl::hopscotch_map<uint32_t, uint32_t> globalTxpCounts;
        for (size_t vId=0; vId<bestMcc.size(); vId++){
          uint32_t vertex = bestMcc[vId];
          uint32_t eqclassId = g.getEqclassId(vertex);

          for (uint32_t txp: txpGroups[eqclassId]){
            globalTxpCounts[txp] += 1;
          }
        }//end-mcc for

        // only transcripts that occur in *every* vertex, and hence
        // have a count of bestMcc.size(), are in the proper intersection.
        uint32_t requiredCount = bestMcc.size();
        std::vector<uint32_t> globalTxps;
        globalTxps.reserve(globalTxpCounts.size());
        for (auto kv : globalTxpCounts) {
          if (kv.second == requiredCount) { globalTxps.push_back(kv.first); }
        }

        if( globalTxps.size() == 0 ) {
          std::cerr << "can't find a representative transcript for a molecule\n"
                    << "Please report this on github ";
          exit(74);
        }

        uint16_t readspmol {0};
        for (auto rv: bestMcc){
          uint32_t eId = g.getEqclassId(rv);
          uint32_t uId = g.getUId(rv);
          auto it = umiGroups[eId].begin();
          std::advance(it, uId);
          readspmol += it->second;
          vset.erase(rv);
        }

        spp::sparse_hash_set<uint32_t> globalGenes;
        for(auto txp: globalTxps){
          uint32_t geneId = getGeneId(t2gMap, txp);
          globalGenes.insert(geneId);
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
      uint32_t eId = g.getEqclassId(vertex);
      uint32_t uId = g.getUId(vertex);
      auto it = umiGroups[eId].begin();
      std::advance(it, uId);
      uint16_t readspmol = it->second;

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

void assignTiers(std::vector<TGroupT>& txpGroups,
                 spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                 std::vector<uint8_t>& tiers) {
  // adding tiers to the genes
  std::vector<std::vector<uint32_t>> geneClasses;
  spp::sparse_hash_map<uint32_t, uint32_t> vertexIndices;
  for (auto& eclass: txpGroups) {
    spp::sparse_hash_set<uint32_t> genes;
    for(auto txp: eclass){
      uint32_t gene = getGeneId(txpToGeneMap, txp);
      genes.insert(gene);
    }

    // first tier
    if (genes.size() == 1){
      tiers[*genes.begin()] = 1;
    }
    else {
      // have to re parse for second and third tier
      std::vector<uint32_t> geneClass (genes.begin(), genes.end());
      geneClasses.emplace_back(geneClass);

      // populate gene indices
      for(auto gene: geneClass){
        if (!vertexIndices.contains(gene)){
          auto gid = vertexIndices.size();
          vertexIndices[gene] = gid;
        }
      }
    }//end-else
  }//end-for

  //generating edges for the graph
  spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t>> edges;
  for (auto& geneClass: geneClasses) {
    for (size_t i=0; i<geneClass.size()-1; i++) {
      for(size_t j=i+1; j<geneClass.size(); j++) {
        uint32_t gene_from = geneClass[i];
        uint32_t gene_to = geneClass[j];

        if (gene_from == gene_to) {
          continue;
        }

        if (!vertexIndices.contains(gene_from) or
            !vertexIndices.contains(gene_to)){
          std::cerr<<"ERROR: Tier creation can't match indexToGene"<<std::flush;
          std::exit(74);
        }
        uint32_t gfromIndex = vertexIndices[gene_from];
        uint32_t gtoIndex = vertexIndices[gene_to];
        edges[ gfromIndex ].insert( gtoIndex );
      }//end-j
    }//end-i
  }//end-geneclass

  // iterating over edges and filling the graph
  typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > AdjList;
  AdjList adjList;
  for (auto& it: edges) {
    uint32_t source = it.first;
    for(uint32_t target: it.second) {
      boost::add_edge(source, target, adjList);
    }
  }

  //find the connected component
  std::vector<uint32_t> component(num_vertices(adjList));
  uint32_t numComps = boost::connected_components(adjList, component.data());

  // making sets of relevant connected vertices
  std::vector<std::vector<uint32_t>> comps (numComps);
  for (size_t i=0; i<component.size(); i++) {
    comps[component[i]].emplace_back(static_cast<uint32_t>(i));
  }

  std::vector<uint32_t> indexToGene(vertexIndices.size());
  for(auto& it: vertexIndices){
    indexToGene[it.second] = it.first;
  }

  if(component.size() != vertexIndices.size()){
    std::cerr<<"ERROR: tiers size doesn't match";
    std::exit(74);
  }

  // iterating over connected components and assigning tiers
  for (auto& comp: comps) {
    bool tier2flag = false;
    for(auto geneIndex: comp) {
      if (geneIndex >= indexToGene.size()){
        std::cerr<<"ERROR:" << geneIndex
                 <<" gene Index > indexToGene size"
                 << indexToGene.size()
                 << std::flush;
        std::exit(74);
      }
      uint32_t gene = indexToGene[geneIndex];
      if (tiers[gene] == 1){
        tier2flag = true;
        break;
      }
    }//end gene for

    uint8_t tierCategory = tier2flag ? 2:3;
    for(auto geneIndex: comp){
      uint32_t gene = indexToGene[geneIndex];
      tiers[gene] = tierCategory;
    } //end-for
  }//end eclass for
}

bool dedupClasses(std::vector<double>& geneAlphas,
                  double& totalUMICount,
                  std::vector<TGroupT>& txpGroups,
                  std::vector<UGroupT>& umiGroups,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  uint32_t umiLength,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  std::vector<uint8_t>& tiers,
                  GZipWriter& gzw, uint32_t umiEditDistance,
                  bool dumpUmiGraph, std::string& trueBarcodeStr,
                  std::vector<spp::sparse_hash_map<uint16_t, uint16_t>>& arboEqClassCount,
                  bool dumpArborescences,
                  std::atomic<uint64_t>& totalUniEdgesCounts,
                  std::atomic<uint64_t>& totalBiEdgesCounts){
  // make directed graph from eqclasses
  alevin::graph::Graph g;
  graphFromCell(txpGroups, umiGroups, g,
                umiEditDistance,
                umiLength,
                totalUniEdgesCounts,
                totalBiEdgesCounts);

  if (dumpUmiGraph){
    gzw.writeUmiGraph(g, trueBarcodeStr);
  }

  // assign tiers to the genes
  assignTiers(txpGroups, txpToGeneMap, tiers);

  if (dumpArborescences) {
    // make gene based eqclasses
    getNumMoleculesWithArborescence(g, txpGroups, umiGroups,
                                    txpToGeneMap, arboEqClassCount,
                                    salmonEqclasses);
  } else {
    // make gene based eqclasses
    getNumMolecules(g, txpGroups, umiGroups,
                    txpToGeneMap,
                    salmonEqclasses);
  }

  for( auto& eqclass: salmonEqclasses ) {
    if ( eqclass.labels.size() == 1 ) {
      totalUMICount += eqclass.count;
      geneAlphas[eqclass.labels.front()] += eqclass.count;
    }
    else if (eqclass.labels.size() == 0){
      std::cerr<<"Eqclasses with No gene labels\n";
      return false;
    }
  }

  if (dumpArborescences and arboEqClassCount.size() != salmonEqclasses.size()) {
    std::cerr<<"Incorrect Arborescence features\n";
    return false;
  }

  return true;
}
