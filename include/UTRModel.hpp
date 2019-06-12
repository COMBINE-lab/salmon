#ifdef __UTR_MODEL_HPP__
#define __UTR_MODEL_HPP__

/** A class that contains the UTR structure
    reads the UTR tsv file with the following
    header that tells where the three prime
    UTR is.
    <tid> <pos> <tid list>
**/

#include <vector>
#include <string>
#include <unordered_map>

#include "SalmonUtils.hpp"
#include "Transcripts.hpp"
#include "iostream"

using std::string

inline string stripVersion(string& trName){
  string::size_type pos = trName.find('.') ;

  if(pos != string::npos){
    return trName.substr(0, pos) ;
  }else{
    return trName ;
  }
}


inline void split(
   const std::string& str,
   std::vector<std::string>& tokens,
   const std::string& delim
){
  const std::string whiteSpace = " " ;
  size_t prev = 0, pos = 0;
  do
    {
      pos = str.find(delim, prev);
      if (pos == std::string::npos) pos = str.length();
      std::string token = str.substr(prev, pos-prev);
      if (!token.empty() and token != whiteSpace) tokens.push_back(token);
      prev = pos + delim.length();
    }
  while (pos < str.length() && prev < str.length());
}

class UTR{
public:
  UTR(uint32_t trnascriptIdIn, size_t startPosIn,
      std::vector<uint32_t> compatibleTranscriptIdsIn):
    UTRId(UTRIdIn),
    transcriptId(transcriptIdIn),
    startPos(startPosIn),
    compatibleTranscriptIds(compatibleTranscriptIdsIn){}

  uint32_t transcriptId ;
  size_t startPos ;
  std::vector<uint32_t> compatibleTranscriptIds ;
};

class TranscriptUTRMap {
public:
  TranscriptUTRMap(const std::string& fname,
                   std::vector<Transcript>& transcripts,
                   size_t thresholdLength
                   )
  {
   namespace bfs = boost::filesystem;
   if(bfs::exists(fname)){
     std::ifstream UTRMapStream(fname.c_str()) ;
     std::string line ;

     std::unordered_map<std::string, uint32_t> nameToIdMap ;
     for(auto& it = transcripts.begin(); it != transcripts.end(); ++it){
       nameToIdMap[stripVersion(it->RefName)] = it->id ;
     }

     uint32_t utrId{0} ;
     while(std::getline(UTRMapStream, line)){
       bool addUTR{false} ;

       line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
       std::vector<std::string> tokens ;
       split(line, tokens, "\t") ;
       if(tokens.size() > 3){
         std::string tidName = tokens[0] ;
         uint32_t tid = nameToIdMap[tidName] ;
         uint32_t startPos = std::stoul(tokens[1]) ;
         size_t width = std::stoul(tokens[2]) ;
         std::vector<uint32_t> compatibleTranscriptIds ;

         if(width > thresholdLength){
           addUTR = true ;
           for(size_t i = 3; i < tokens.size() ; ++i){
             compatibleTranscriptsIds.push_back(nameToIdMap[tokens[i]]) ;
           }
         }

       }
       if(addUTR){
         UTRs.emplace_back(tid, startPos, compatibleTranscriptIds) ;
         if(transcriptUTRMap.find(tid) != transcriptUTRMap.end()){
           transcriptUTRMap[tid].emplace_back(utrId) ;
         }else{
           transcriptUTRMap[tid] = {utrId} ;
         }
         utrId += 1 ;
       }
     }


     std::cerr << "# number of UTRs " << UTRs.size() << "\n" ;
   }else{
     std::cerr << "file " << fname << "does not exist \n" ;
   }

  }

  std::vector<UTR> UTRs ; // vector of UTRs
  std::unordered_map<uint32_t, uint32_t> transcriptToUTRMap ; 

} ;

#endif // UTR_MODEL_HPP
