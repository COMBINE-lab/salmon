//#include "FatPufferGraph.hpp"
#include "zstr/zstr.hpp"
#include <fstream>
#include <iostream>
//#include "OurGFAReader.hpp"
#include "CLI/CLI.hpp"
#include "Util.hpp"
#include "sparsepp/spp.h"

#define _verbose(fmt, args...) fprintf(stderr, fmt, ##args)

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}



class Node{

public:


  uint64_t id ;
  std::set<uint64_t> incoming ;
  std::set<uint64_t> outgoing ;

  Node() {}
  Node(uint64_t id_){
    id = id_ ;
  }

  bool addIncoming(uint64_t neighbor){
      incoming.insert(neighbor) ;
      return true;
  }

  bool addOutgoing(uint64_t neighbor){
      outgoing.insert(neighbor) ;
      return true;
  }

  int getDegree(){
    return (incoming.size() + outgoing.size() ) ;
  }


};



void parseGFA(std::string& gfaFile, spp::sparse_hash_map<uint64_t, Node>& g){
  std::string ln;
  std::string tag, id, value;
  //size_t contig_cnt{0};

  std::unique_ptr<zstr::ifstream> file;
  file.reset(new zstr::ifstream(gfaFile));

  std::cerr << "Start loading segments... \n";

  //we just need ids and not sequences
  //spp::sparse_hash_map<std::string, uint64_t> seq2id;

  //uint64_t idCntr = 0;
  uint64_t maxnid{0};
  uint64_t contig_cnt{0} ;
  std::vector<uint64_t> singleNodes ;
  uint64_t pathCount{0} ;

  //int countDebugNodes{0} ;

  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    //skipping everything other than path and contig

    pathCount++ ;
    _verbose("\rNumber of processed paths: %lu", pathCount);

    if (firstC != 'S' and firstC != 'P')
      continue;

    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = pufferfish::util::split(lnview, '\t');
    tag = splited[0].to_string();
    id = splited[1].to_string();

    // count vertex
    if (tag == "S") {
      continue ;
      try {
        uint64_t nid = std::stoll(id);
        if(nid > maxnid)
          maxnid = nid ;
      } catch (std::exception& e) {
        // not a numeric contig id
      }
      contig_cnt++;
    }

		if (firstC != 'P')
			continue;
    if (splited.size() != 3) {
      std::cerr << "the path line should have 3 fields (tab separated) --- skipping!\n";
      continue;
    }

    if (id == "") {
      std::cerr << "interesting; the ID is empty.  Skipping\n";
      continue;
    }
    // A path line
    auto pvalue = splited[2];
    // parse value and add all conitgs to contigVec
    std::vector<std::pair<uint64_t, bool>> contigVec =
      pufferfish::util::explode(pvalue, ',');

		// update graph and add the new set of segments to it
    if(contigVec.size() == 1){
      //singleNodes.push_back(contigVec[0].first) ;
      Node n(contigVec[0].first) ;
      if(g.find(contigVec[0].first) == g.end()){
        g[contigVec[0].first] = n ;
      }
    }
    for (size_t i = 0; i < contigVec.size()-1; i++){
      auto fromNode = contigVec[i].first ;
      auto toNode = contigVec[i+1].first ;

      bool fromNodeSign = contigVec[i].second ;
      bool toNodeSign = contigVec[i+1].second ;

      if(g.find(fromNode) == g.end()){
        Node n(fromNode);
        g[fromNode] = n ;
      }
      if(g.find(toNode) == g.end()){
        Node n(toNode);
        g[toNode] = n ;
      }

      if(toNodeSign){
        g[toNode].addIncoming(fromNode) ;
      }else{
        g[toNode].addOutgoing(fromNode) ;
      }
      if(fromNodeSign){
        g[fromNode].addOutgoing(toNode) ;
      }else{
        g[fromNode].addIncoming(toNode) ;
      }
      if(g[toNode].getDegree() > 8 or g[fromNode].getDegree() > 8){
        std::cerr << "After parsing line "<<lnview<<" bad\n" ;
        std::cerr << toNode << "\t" << g[toNode].getDegree() << "\n"
                  << fromNode << "\t" << g[fromNode].getDegree() << "\n" ;

        std::exit(1) ;
      }

      //g.addEdge(contigVec[i].first, contigVec[i].second, contigVec[i+1].first, contigVec[i+1].second);
    }

    /*
    std::cerr << "After parsing line "<<lnview<<"\n" ;
    if(contigVec.size() > 1){
      for(auto it : g){
        std::cerr << it.first << "\t" << it.second.getDegree() << "\n";
      }
      //std::cerr << "\n" ;
      std::exit(1) ;
      }*/
    //<< "processed paths  "<< pathCount << " %\r"; }
	}
  std::cerr << "\nNumber of vertices " << contig_cnt << "\n" ;
  std::cerr << "Done Reading gfa file\n";
 

}

int main(int argc, char* argv[]) {
  (void)argc;

  CLI::App app{"Calculate Edge density"};

  std::string gfaFile ;// = "/mnt/scratch7/pufferfish_data/gencode.v25.pc_transcripts_fixed.k27.pufferized.gfa";
  std::string edgeFile ;

  app.add_option("-g,--gfa", gfaFile, "path to the input GFA file")
    ->required();
  app.add_option("-o,--out", edgeFile, "path to the out file")
    ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }



  spp::sparse_hash_map<uint64_t, Node> g ;
  parseGFA(gfaFile, g) ;

  //Go over all the vertices and 
  //calculate the statistics
  
  //spp::sparse_hash_map<uint64_t, pufg::Node>& nodes = g.getVertices();
  double avgDensity{0} ;

  /*
  for(const auto& x : nodes){
    auto node = x.second ;
    auto inEdges = node.getPredecessors() ;
    auto outEdges = node.getSuccessors() ;

    size_t totalEdges = inEdges.size() + outEdges.size() ;
    avgDensity += double(totalEdges)/8.0 ; 
    }*/

  std::ofstream outfile(edgeFile) ;
  for(auto n : g){
    avgDensity += double(n.second.getDegree())/8.0 ;
    outfile << n.first << "\t" << n.second.getDegree() << "\n" ;
  }
  outfile.close() ;

  avgDensity = avgDensity / double(g.size()) ;

  std::cout << "Edge density for "<<gfaFile<<" "<<avgDensity << "\n" ;


}
