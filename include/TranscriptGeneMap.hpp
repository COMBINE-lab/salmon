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

#ifndef TRANSCRIPT_GENE_MAP_HPP
#define TRANSCRIPT_GENE_MAP_HPP

// Allows for the serialization of this class
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <algorithm>
#include <unordered_map>
#include <vector>

class TranscriptGeneMap {
  using Index = size_t;
  using Size = size_t;
  using NameVector = std::vector<std::string>;
  using IndexVector = std::vector<size_t>;
  using IndexVectorList = std::vector<std::vector<size_t>>;

private:
  NameVector _transcriptNames;
  NameVector _geneNames;
  IndexVector _transcriptsToGenes;
  IndexVectorList _genesToTranscripts;
  bool _haveReverseMap;

  void _computeReverseMap() {

    _genesToTranscripts.resize(_geneNames.size(), {});

    Index transcriptID = 0;
    size_t maxNumTrans = 0;
    Index maxGene;
    for (transcriptID = 0; transcriptID < _transcriptsToGenes.size();
         ++transcriptID) {
      _genesToTranscripts[_transcriptsToGenes[transcriptID]].push_back(
          transcriptID);
      if (maxNumTrans <
          _genesToTranscripts[_transcriptsToGenes[transcriptID]].size()) {
        maxNumTrans =
            _genesToTranscripts[_transcriptsToGenes[transcriptID]].size();
        maxGene = _transcriptsToGenes[transcriptID];
      }
    }
    std::cerr << "max # of transcripts in a gene was " << maxNumTrans
              << " in gene " << _geneNames[maxGene] << "\n";
  }

  friend class cereal::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar(_transcriptNames, _geneNames, _transcriptsToGenes, _genesToTranscripts,
       _haveReverseMap);
  }

public:
  TranscriptGeneMap()
      : _transcriptNames(NameVector()), _geneNames(NameVector()),
        _transcriptsToGenes(IndexVector()), _haveReverseMap(false) {}

  TranscriptGeneMap(const NameVector& transcriptNames,
                    const NameVector& geneNames,
                    const IndexVector& transcriptsToGenes)
      : _transcriptNames(transcriptNames), _geneNames(geneNames),
        _transcriptsToGenes(transcriptsToGenes), _haveReverseMap(false) {}

  TranscriptGeneMap(const TranscriptGeneMap& other) = default;
  TranscriptGeneMap& operator=(const TranscriptGeneMap& other) = default;

  Index INVALID{std::numeric_limits<Index>::max()};

  Index findTranscriptID(const std::string& tname) {
    using std::distance;
    using std::lower_bound;
    auto it =
        lower_bound(_transcriptNames.begin(), _transcriptNames.end(), tname);
    if (it == _transcriptNames.end() or *it != tname) {
      return INVALID;
    } else {
      return distance(_transcriptNames.begin(), it);
    }
  }

  Size numTranscripts() { return _transcriptNames.size(); }
  Size numGenes() { return _geneNames.size(); }

  bool needReverse() {
    if (_haveReverseMap) {
      return false;
    } else {
      _computeReverseMap();
      return true;
    }
  }

  const IndexVector& transcriptsForGene(Index geneID) {
    return _genesToTranscripts[geneID];
  }

  inline std::string nameFromGeneID(Index geneID) { return _geneNames[geneID]; }
  inline Index gene(Index transcriptID) {
    return _transcriptsToGenes[transcriptID];
  }
  inline std::string geneName(Index transcriptID) {
    return _geneNames[_transcriptsToGenes[transcriptID]];
  }
  inline std::string geneName(const std::string& transcriptName, bool& found) {
    found = false;
    auto tid = findTranscriptID(transcriptName);
    if (tid != INVALID) {
      found = true;
      return geneName(tid);
    } else {
      found = false;
      return transcriptName;
    }
  }

  inline std::string transcriptName(Index transcriptID) {
    return _transcriptNames[transcriptID];
  }
};

#endif // TRANSCRIPT_GENE_MAP_HPP
