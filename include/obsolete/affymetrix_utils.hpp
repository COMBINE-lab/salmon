
/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* vim: set ts=2 et sw=2 tw=80: */
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>

#include <seqan/sequence.h>

namespace affymetrix {

  namespace utils {

    using std::string;
    using std::vector;


class ProbeTarget{
public:
friend std::ostream& operator<< (std::ostream&, const ProbeTarget&);
ProbeTarget( const string& _geneName, const seqan::DnaString& _probeSeq );
ProbeTarget();
  std::string geneName;
  seqan::DnaString probeSeq;
};

    typedef std::unordered_map<string, ProbeTarget> StringMapT;
    typedef std::unordered_set<string> StringSetT;

    /**
     * Given a file containing the probe information, build a reverse mapping from
     * probe name to ENSEMBL gene name.  Any probe mapping to > 1 gene will be discarded.
     */
    StringMapT ComputeReverseMap( const string& fname );

    /**
     * Given a file containing a list of strings, return a vector of corresponding
     * paths.  This is the list of probe files to be read.
     */
    vector< boost::filesystem::path > GetProbeFiles( const string& fname );

    class AffyEntry {
    public:
      AffyEntry( const string& _probeName, float _expression, float _background );

      std::string probeName;
      float expression;
      float background;
    };

  }

}
