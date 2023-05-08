#include "FatPufferGraph.hpp"
#include "zstr/zstr.hpp"
#include <fstream>
#include <iostream>
//#include "OurGFAReader.hpp"
#include "CLI/CLI.hpp"
#include "Util.hpp"





void parseGFA(std::string& gfaFile, pufg::Graph& g){
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

  //int countDebugNodes{0} ;

  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    //skipping everything other than path and contig
    if (firstC != 'S' and firstC != 'P')
      continue;

    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = pufferfish::util::split(lnview, '\t');


    tag = splited[0].to_string();
    id = splited[1].to_string();

    // add vertex
    if (tag == "S") {
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

		//contigVec = convertOld2newPath(contigVec);

		// update graph and add the new set of segments to it
    if(contigVec.size() == 1){
      singleNodes.push_back(contigVec[0].first) ;
    }
    for (size_t i = 0; i < contigVec.size()-1; i++) 
      g.addEdge(contigVec[i].first, contigVec[i].second, contigVec[i+1].first, contigVec[i+1].second);


    //debugging line 
    //sanityCheck(g) ;
    if(contigVec.size() > 1){
      auto nodes = g.getVertices() ;
      //std::cerr << "num of nodes " << nodes.size() << "\n" ;
      for(const auto& x : nodes){
        auto node = x.second ;
        auto inEdges = node.getPredecessors() ;
        auto outEdges = node.getSuccessors() ;

        size_t totalEdges = inEdges.size() + outEdges.size() ;
        if(totalEdges > 8){
          std::cerr << "After adding line "<<lnview << "\n" ;
          std::cerr << "num of nodes " << nodes.size() << "\n" ;
          std::cerr << node.getId() << ": " << totalEdges << "\n" ;
          std::exit(1) ;
        }
        //avgDensity += double(totalEdges)/8.0 ; 
      }
    }

		// set pathStart and pathEnd (do we need this)
		// setPathStart(id, contigVec[0]); (do we need this)
		// setDonTouchSide(contigVec); (do we need this)
	}
  std::cerr << "Number of vertices " << contig_cnt << "\n" ;
  std::cerr << "Done Reading gfa file\n";
 

}

int main(int argc, char* argv[]) {
  (void)argc;

  CLI::App app{"Calculate Edge density"};

  std::string gfaFile ;// = "/mnt/scratch7/pufferfish_data/gencode.v25.pc_transcripts_fixed.k27.pufferized.gfa";

  app.add_option("-g,--gfa", gfaFile, "path to the input GFA file")
    ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }



  pufg::Graph g ;
  parseGFA(gfaFile, g) ;

  //Go over all the vertices and 
  //calculate the statistics

  spp::sparse_hash_map<uint64_t, pufg::Node>& nodes = g.getVertices();
  double avgDensity{0} ;

  for(const auto& x : nodes){
    auto node = x.second ;
    auto inEdges = node.getPredecessors() ;
    auto outEdges = node.getSuccessors() ;

    size_t totalEdges = inEdges.size() + outEdges.size() ;
    avgDensity += double(totalEdges)/8.0 ; 
  }

  avgDensity = avgDensity / double(nodes.size()) ;

  std::cout << "Edge density for "<<gfaFile<<" "<<avgDensity << "\n" ;


}
