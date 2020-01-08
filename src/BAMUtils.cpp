#include "BAMUtils.hpp"
#include <unordered_map>
#include <string>
#include <functional>
#include <algorithm>
#include <cctype>

namespace salmon {
  namespace bam_utils {
    std::string to_string(AlignerDetails a) {
      switch (a) {
      case AlignerDetails::BOWTIE2:
        return "Bowtie2";
      case AlignerDetails::BOWTIE2_LOCAL:
        return "Bowtie2-local";
      case AlignerDetails::BWA_MEM:
        return "BWA";
      case AlignerDetails::STAR:
        return "STAR";
      case AlignerDetails::MINIMAP2:
        return "minimap2";
      case AlignerDetails::RAPMAP:
        return "RapMap";
      case AlignerDetails::PUFFERFISH:
        return "pufferfish";
      case AlignerDetails::UNKNOWN:
        return "unknown";
      default:
        return "unknown";
      }
    }

    AlignerDetails inferAlignerFromHeader(SAM_hdr* header) {
      std::unordered_map<std::string, AlignerDetails> alignerNameToType = 
        {
         {"bowtie2", AlignerDetails::BOWTIE2},
         {"star", AlignerDetails::STAR},
         {"bwa", AlignerDetails::BWA_MEM},
         {"minimap2", AlignerDetails::MINIMAP2},
         {"rapmap", AlignerDetails::RAPMAP},
         {"pufferfish", AlignerDetails::PUFFERFISH}
        };

      AlignerDetails aligner = AlignerDetails::UNKNOWN;

      char pgtag[] = {'P', 'G'};
      char pntag[] = {'P', 'N'};
      char idtag[] = {'I', 'D'};
      char vntag[] = {'V', 'N'};
      char cltag[] = {'C', 'L'};

      auto pginfo = sam_hdr_find(header, pgtag, NULL, NULL);
      if (pginfo) {
        auto pntagval = sam_hdr_find_key(header, pginfo, pntag, NULL);
        if (pntagval) {
          std::string progname(pntagval->str+3, pntagval->len-3);
          std::transform(progname.begin(), progname.end(), progname.begin(), [](unsigned char c) { return std::tolower(c); });
          if (alignerNameToType.find(progname) != alignerNameToType.end()) {
            aligner = alignerNameToType[progname];
          }
        }
        auto cltagval = sam_hdr_find_key(header, pginfo, cltag, NULL);
        if (cltagval) {
          std::string clstring(cltagval->str+3, cltagval->len-3);
          // if the aligner is bowtie2, check there is a `--local` flag
          if (aligner == AlignerDetails::BOWTIE2) {
            auto lpos = clstring.find("--local");
            if (lpos != clstring.npos) { aligner = AlignerDetails::BOWTIE2_LOCAL; }
          }
        }
      }

      return aligner;
    }

  }
}

