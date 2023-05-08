#include "FastxParser.hpp"
#include "jellyfish/mer_dna.hpp"
#include "string_view.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <type_traits>
#include <vector>
//#include "PufferFS.hpp"
#include "BooPHF.hpp"
//#include "OurGFAReader.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "CLI/CLI.hpp"


int main(int argc, char* argv[]) {
  (void)argc;
  std::vector<std::string> read_file;//= {argv[1]};
  std::string rfname;//= {argv[1]};
  std::string gfa_file;//= argv[2];

  //std::cerr << "\n fasta file " << argv[1] << "\n";
  //std::cerr << "\n gfa file " << argv[2] << "\n";

  int k = 31;
  CLI::App app{"Pufferfish : An efficient dBG index."};
  app.add_option("-g,--gfa", gfa_file, "GFA file")->required();
  app.add_option("-r,--ref", rfname, "reference file")->required();
  app.add_option("-k,--klen", k, "kmer length", 31);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }
  std::cerr << "k = " << k << "\n";

  size_t overlap = k - 1;
  read_file.push_back(rfname);
  // std::string outfile = "contigs.fa" ;

  spp::sparse_hash_map<std::string, std::string> fastaMap;

  std::cout << "\n starting .. dodo\n";
  {
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    uint32_t rn = 0;
    while (parser.refill(rg)) {
      // Here, rg will contain a chunk of read pairs
      // we can process.
      for (auto& rp : rg) {
        if (rn % 10000 == 0) {
          std::cerr << "rn : " << rn << "\n";
        }
        ++rn;

        auto& r1 = rp.seq;
        fastaMap[rp.name] = r1;
      }
    }
  }

  std::cerr << "\n fasta file contains " << fastaMap.size() << " \n";

  std::ifstream file(gfa_file);

  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  // size_t ref_cnt{0};
  spp::sparse_hash_map<std::string, bool> touchedSegment;
  spp::sparse_hash_map<uint64_t, std::string> contigid2seq;
  spp::sparse_hash_map<std::string, std::string> reconstructedTr;

  //int k = 27;
  //size_t overlap = k - 1;

  std::stringstream reconStream;

  while (std::getline(file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S' and firstC != 'P')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = pufferfish::util::split(lnview, '\t');
    tag = splited[0].to_string();
    id = splited[1].to_string();
    value = splited[2].to_string();
    // A segment line
    if (tag == "S") {
      if (pufferfish::util::is_number(id)) {
		uint64_t contigId = std::stoll(id);
        contigid2seq[contigId] = value;
        // std::cerr << value << "\n" ;
        // std::exit(2) ;
        // touchedSegment[id] = false ;
      }
      contig_cnt++;
    }

    // A path line
    if (tag == "P") {
      auto pvalue = splited[2];
      std::vector<std::pair<uint64_t, bool>> contigVec =
          pufferfish::util::explode(pvalue, ',');
      // parse value and add all conitgs to contigVec
      // if(reconstructedTr[id] != "") continue ;
      // reconstructedTr[id] = "";
      auto& refStr = fastaMap[splited[1].to_string()];
      reconStream.clear();
      reconStream.str(std::string());
      size_t i = 0;
      for (auto core : contigVec) {
        auto contig_id = core.first;
        auto ore = core.second;
        std::string added;
        if (i != contigVec.size() - 1) {
          if (!ore) {
            added = pufferfish::util::revcomp(contigid2seq[contig_id]);
          } else {
            added = contigid2seq[contig_id];
          }
          // contigid2seq[contig_id].erase(contigid2seq[contig_id].size()-31+1,31)
          // ;
          if (added.size() > overlap - 1) {
            added.erase(added.size() - overlap, overlap);
            reconStream << added;
            // reconstructedTr[id] += added ;
          } else {
            std::cerr << "\nBAAAAAAAD\n" << contig_id << " : " << added.size();
          }
        } else {
          if (!ore) {
            added = pufferfish::util::revcomp(contigid2seq[contig_id]);
          } else {
            added = contigid2seq[contig_id];
          }
          // reconstructedTr[id] += added ;
          reconStream << added;
        }
        i++;
      }

      if (reconStream.str() != refStr) {
        std::cerr
            << "Reference string and reconstructed string do not match for " << splited[1].to_string() << "\n";
	std::cerr << "reconstructed = " << reconStream.str() << "\n";
	std::cerr << "ref str       = " << refStr << "\n";
      }

      // std::cerr << id << "\n" ;
      /*
      if(fastaMap[id] != reconstructedTr[id]){
          std::cerr << id << "\n" ;
          std::cerr << "true\n " << fastaMap[id] << " " << fastaMap[id].size()
      << "\n" ;
          std::cerr << "reconstructed\n " << reconstructedTr[id] << " " <<
      reconstructedTr[id].size() << "\n" ;
          std::cerr << " number of contigs " << contigVec.size() << "\n" ;
          size_t j = 0;
          for(auto core : contigVec){
              auto contig_id = core.first ;
              auto ore = core.second ;
              std::string added ;
              if (j != contigVec.size()-1){
                  if(!ore){
                      added =  pufferfish::util::revcomp(contigid2seq[contig_id]) ;
                  }else{
                      added = contigid2seq[contig_id] ;
                  }
                  //contigid2seq[contig_id].erase(contigid2seq[contig_id].size()-31+1,31)
      ;
                  if(added.size() > overlap){
                      //added.erase(added.size()-overlap,overlap) ;
                      //reconstructedTr[id] += added ;
                      std::cerr << contig_id << " " << added.substr(0,
      added.size()-overlap) << "\t" << added.substr(added.size()-overlap) <<
      "\n" ;
                      //std::cerr << added << "\n" ;
                  }
                      //std::cerr << contig_id << " " << added << "\n" ;
              }else{
                  if(!ore){
                      added =  pufferfish::util::revcomp(contigid2seq[contig_id]) ;
                  }else{
                      added = contigid2seq[contig_id] ;
                  }
                  //reconstructedTr[id] += added ;
                  std::cerr << contig_id << " " << added << "\n" ;
              }
              j++ ;
          }


          std::exit(1) ;
      }
*/
    }
  }

  return 0;

  /*
  int found = 0;
  int notFound = 0;
  int notSub = 0;
  for (auto& kv : fastaMap) {
    if (kv.second == reconstructedTr[kv.first]) {
      found++;
    } else if (kv.second.find(reconstructedTr[kv.first]) != -1 and
               (kv.second.size() - 1 == reconstructedTr[kv.first].size())) {
      // std::cerr << "\n" << notSub << "\n\t" << kv.second << "\n\t" <<
      // reconstructedTr[kv.first];
      notSub++;
    } else {
      // std::cerr << "tid " << kv.first << "\n" ;
      // std::cerr << "true seq " << kv.second << "\n" ;
      // std::cerr << "our seq " << reconstructedTr[kv.first] << "\n" ;
      //				std::cout << kv.second.size() << "\t" <<
      //reconstructedTr[kv.first].size() << "\n" ;
      notFound++;
      // std::exit(1) ;
    }
  }

  std::cerr << "\n Constructed " << found << " Not Constructed " << notFound
            << " Sub-strings " << notSub << "\n";

  return 0;
  */
}
