/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#ifndef TRANSCRIPT_ALLELE_MAP_HPP
#define TRANSCRIPT_ALLELE_MAP_HPP

// Allows for the serialization of this class
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <algorithm>
#include <unordered_map>
#include <vector>

class TranscriptAlleleMap {
  using Index = size_t;
  using Size = size_t;
  using NameVector = std::vector<std::string>;
  using IndexVector = std::vector<size_t>;
  using IndexVectorList = std::vector<std::vector<size_t>>;
  using NameMap = std::unordered_map<std::string, size_t>;

private:
  NameMap _transcriptNameMap;
  IndexVectorList _transcriptToAlleleMap;


  friend class cereal::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar(_transcriptNameMap, _transcriptToAlleleMap);
  }

public:
  TranscriptAlleleMap()
      : _transcriptNameMap(NameMap()),
        _transcriptToAlleleMap(IndexVectorList()){}

  TranscriptAlleleMap(const NameMap& transcriptNameMap,
                    const IndexVectorList& transcriptToAlleleMap)
      : _transcriptNameMap(transcriptNameMap),
        _transcriptToAlleleMap(transcriptToAlleleMap){}

  TranscriptAlleleMap(const TranscriptAlleleMap& other) = default;
  TranscriptAlleleMap& operator=(const TranscriptAlleleMap& other) = default;

  Index INVALID{std::numeric_limits<Index>::max()};

  Index findTranscriptID(const std::string& tname) {
    auto trIt = _transcriptNameMap.find(tname) ;
    if(trIt != _transcriptNameMap.end()){
      return trIt->second ; 
    }else{
      return INVALID ;
    }
  }

  IndexVector getAlleleVector(const size_t tid){
    return _transcriptToAlleleMap[tid] ;
  }

  bool isPosIntersectingAllele(const std::string& tname,
                               const size_t pos){
    auto tid = findTranscriptID(tname) ;
    auto allelePosVector = _transcriptToAlleleMap[tid] ;
    size_t range{10} ;

    // Search for position in the allelePosVector
    // std::lower_bound gives the a value x such that
    // n >= a
    auto nextBigger = std::lower_bound(allelePosVector.begin(),
                                        allelePosVector.end(), pos) ;
    auto indexOfNextBigger = std::distance(nextBigger,
                                              allelePosVector.begin()) ;

    if (indexOfNextBigger == allelePosVector.size()){
      return false ;
    }else{
      size_t lowerIndex{0};
      size_t upperIndex{indexOfNextBigger};
      if(*nextBigger == pos){
        return true ;
      }
      if(indexOfNextBigger > 0){
        lowerIndex = indexOfNextBigger - 1 ;
      }
      // Check if lower bound is within the range
      auto distanceFromAllele = std::abs(static_cast<int>(pos-lowerIndex)) ;
      if(distanceFromAllele <= range){
        return true ;
      }
    }

    return false ;

  }

  Size numTranscripts() { return _transcriptNameMap.size(); }

};

#endif // TRANSCRIPT_ALLELE_MAP_HPP
