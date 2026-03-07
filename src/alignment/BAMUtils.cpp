#include "salmon/internal/alignment/BAMUtils.hpp"
#include <unordered_map>
#include <string>
#include <functional>
#include <algorithm>
#include <cctype>
#include <sstream>

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

    AlignerDetails inferAlignerFromHeader(AlignmentHeader* header) {
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
      if (header == nullptr || header->raw == nullptr) {
        return aligner;
      }

      const char* headerText = sam_hdr_str(header->raw);
      if (headerText == nullptr) {
        return aligner;
      }

      std::istringstream lines(headerText);
      std::string line;
      while (std::getline(lines, line)) {
        if (line.rfind("@PG\t", 0) != 0) {
          continue;
        }

        std::string progname;
        std::string commandLine;
        std::istringstream fields(line);
        std::string field;
        while (std::getline(fields, field, '\t')) {
          if (field.rfind("PN:", 0) == 0) {
            progname = field.substr(3);
          } else if (field.rfind("CL:", 0) == 0) {
            commandLine = field.substr(3);
          }
        }

        std::transform(progname.begin(), progname.end(), progname.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        auto alignerIt = alignerNameToType.find(progname);
        if (alignerIt != alignerNameToType.end()) {
          aligner = alignerIt->second;
        }
        if (aligner == AlignerDetails::BOWTIE2 &&
            commandLine.find("--local") != commandLine.npos) {
          aligner = AlignerDetails::BOWTIE2_LOCAL;
        }
        if (aligner != AlignerDetails::UNKNOWN) {
          break;
        }
      }

      return aligner;
    }

  }
}
