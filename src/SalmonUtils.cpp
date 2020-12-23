#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/range/join.hpp>
#include <boost/thread/thread.hpp>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "json.hpp"
#include "tbb/combinable.h"
#include "tbb/parallel_for.h"

#include "AlignmentLibrary.hpp"
#include "DistributionUtils.hpp"
#include "GCFragModel.hpp"
#include "KmerContext.hpp"
#include "LibraryFormat.hpp"
#include "ReadExperiment.hpp"
#include "ReadPair.hpp"
#include "SBModel.hpp"
#include "SalmonMath.hpp"
#include "SalmonUtils.hpp"
#include "TryableSpinLock.hpp"
#include "UnpairedRead.hpp"
#include "TranscriptGroup.hpp"
#include "Transcript.hpp"

#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/spdlog.h"

#include "gff.h"

#include "FastxParser.hpp"
//#include "jellyfish/mer_dna.hpp"

#include "GenomicFeature.hpp"
#include "SGSmooth.hpp"
#include "TranscriptGeneMap.hpp"

#include "StadenUtils.hpp"
#include "SalmonDefaults.hpp"

#include "pufferfish/Util.hpp"

#include "zstr.hpp"

namespace salmon {
namespace utils {

  using MateStatus = pufferfish::util::MateStatus;

std::string str(const MappingType& mt) {
  switch (mt) {
  case MappingType::UNMAPPED:
    return "u";
  case MappingType::LEFT_ORPHAN:
    return "m1";
  case MappingType::RIGHT_ORPHAN:
    return "m2";
  case MappingType::BOTH_ORPHAN:
    return "m12";
  case MappingType::PAIRED_MAPPED:
    return "mp";
  case MappingType::SINGLE_MAPPED:
    return "ms";
  case MappingType::DECOY:
    return "d";
  }
  // should never get here!
  return "E";
}

bool headersAreConsistent(SAM_hdr* h1, SAM_hdr* h2) {

  bool consistent{true};
  // Both files must contain the same number of targets
  if (h1->nref != h2->nref) {
    consistent = false;
  }

  // Check each target to ensure that the name and length are the same.
  size_t i = 0;
  size_t n = h1->nref;
  while (consistent and i < n) {
    size_t l1 = h1->ref[i].len;
    size_t l2 = h2->ref[i].len;
    consistent = (l1 == l2) and (strcmp(h1->ref[i].name, h2->ref[i].name) == 0);
    ++i;
  }

  return consistent;
}

bool headersAreConsistent(std::vector<SAM_hdr*>&& headers) {
  if (headers.size() == 1) {
    return true;
  }

  // Ensure that all of the headers are consistent (i.e. the same), by
  // comparing each with the first.
  bool consistent{true};
  auto itFirst = headers.begin();
  auto it = itFirst;
  while (++it != headers.end()) {
    if (!headersAreConsistent(*itFirst, *it)) {
      consistent = false;
      break;
    }
  }
  return consistent;
}

std::ostream& operator<<(std::ostream& os, OrphanStatus s) {
  switch (s) {
  case OrphanStatus::LeftOrphan:
    os << "left orphan";
    break;
  case OrphanStatus::RightOrphan:
    os << "right orphan";
    break;
  case OrphanStatus::Paired:
    os << "paired";
    break;
  }
  return os;
}

bool isCompatible(const LibraryFormat observed, const LibraryFormat expected,
                  int32_t start, bool isForward, MateStatus ms) {
  // If we're dealing with a single end read.
  bool compat{false};
  if (ms != MateStatus::PAIRED_END_PAIRED) {
    compat = compatibleHit(expected, start, isForward, ms);
  } else {
    compat = compatibleHit(expected, observed);
  }
  return compat;
}

double logAlignFormatProb(const LibraryFormat observed,
                          const LibraryFormat expected, int32_t start,
                          bool isForward, MateStatus ms,
                          double incompatPrior) {
  // If we're dealing with a single end read.
  bool compat{false};
  if (ms != MateStatus::PAIRED_END_PAIRED) {
    compat = compatibleHit(expected, start, isForward, ms);
  } else {
    compat = compatibleHit(expected, observed);
  }
  return (compat) ? salmon::math::LOG_1 : incompatPrior;
  /** Old compat code
  if (expected.type == ReadType::PAIRED_END and
      observed.type == ReadType::SINGLE_END) {
      double logOrphanProb = salmon::math::LOG_ORPHAN_PROB;
      if (expected.strandedness == ReadStrandedness::U or
          expected.strandedness == ReadStrandedness::AS or
          expected.strandedness == ReadStrandedness::SA) {
          return salmon::math::LOG_1;
      } else {
          return (expected.strandedness == observed.strandedness) ?
  logOrphanProb : incompatPrior;
      }
  } else if (observed.type != expected.type or
      observed.orientation != expected.orientation ) {
      return incompatPrior;
  } else {
      if (expected.strandedness == ReadStrandedness::U) {
          return salmon::math::LOG_ONEHALF;
      } else {
          if (expected.strandedness == observed.strandedness) {
              return salmon::math::LOG_1;
          } else {
              return incompatPrior;
          }
      }
  }

  fmt::print(stderr, "WARNING: logAlignFormatProb --- should not get here");
  return salmon::math::LOG_0;
  */
}

// for single end reads or orphans
bool compatibleHit(const LibraryFormat expected, int32_t start, bool isForward,
                   MateStatus ms) {
  auto expectedStrand = expected.strandedness;
  auto expectedType = expected.type;
  switch (ms) {
  case MateStatus::SINGLE_END:
    if (isForward) { // U, SF
      return (expectedStrand == ReadStrandedness::U or
              expectedStrand == ReadStrandedness::S);
    } else { // U, SR
      return (expectedStrand == ReadStrandedness::U or
              expectedStrand == ReadStrandedness::A);
    }
    break;
    // The next two cases are for *orphaned* PE reads

    // This case is where the mapped read belongs to the "left" (i.e. 1) end
    // of the pair.
  case MateStatus::PAIRED_END_LEFT:
    // "M"atching or same orientation is a special case
    /*
    if (expectedType == ReadType::SINGLE_END) {
      return (expectedStrand == ReadStrandedness::U or
              (expectedStrand == ReadStrandedness::S and isForward) or
              (expectedStrand == ReadStrandedness::A and !isForward));
    } else
    */
    if (expected.orientation == ReadOrientation::SAME) {
      return (expectedStrand == ReadStrandedness::U or
              (expectedStrand == ReadStrandedness::S and isForward) or
              (expectedStrand == ReadStrandedness::A and !isForward));
    } else if (isForward) { // IU, ISF, OU, OSF, MU, MSF
      return (expectedStrand == ReadStrandedness::U or
              expectedStrand == ReadStrandedness::SA);
    } else { // IU, ISR, OU, OSR, MU, MSR
      return (expectedStrand == ReadStrandedness::U or
              expectedStrand == ReadStrandedness::AS);
    }
    break;

    // This case is where the mapped read belongs to the "right" (i.e. 2) end
    // of the pair.
  case MateStatus::PAIRED_END_RIGHT:
    // "M"atching or same orientation is a special case
    /*
    if (expectedType == ReadType::SINGLE_END) {
      return (expectedStrand == ReadStrandedness::U or
              (expectedStrand == ReadStrandedness::A and isForward) or
              (expectedStrand == ReadStrandedness::S and !isForward));
    } else
    */
    if (expected.orientation == ReadOrientation::SAME) {
      return (expectedStrand == ReadStrandedness::U or
              (expectedStrand == ReadStrandedness::S and isForward) or
              (expectedStrand == ReadStrandedness::A and !isForward));
    } else if (isForward) { // IU, ISR, OU, OSR, MU, MSR
      return (expectedStrand == ReadStrandedness::U or
              expectedStrand == ReadStrandedness::AS);
    } else { // IU, ISF, OU, OSF, MU, MSF
      return (expectedStrand == ReadStrandedness::U or
              expectedStrand == ReadStrandedness::SA);
    }
    break;
  default:
    // SHOULD NOT GET HERE
    fmt::print(stderr,
               "WARNING: Could not associate known library type with read!\n");
    return false;
    break;
  }
  // SHOULD NOT GET HERE
  fmt::print(stderr,
             "WARNING: Could not associate known library type with read!\n");
  return false;
}

// for paired-end reads
bool compatibleHit(const LibraryFormat expected, const LibraryFormat observed) {
  if (observed.type != ReadType::PAIRED_END) {
    // SHOULD NOT GET HERE
    fmt::print(stderr,
               "WARNING: PE compatibility function called with SE read!\n");
    fmt::print(stderr, "expected: {}, observed: {}\n", expected, observed);
    return false;
  }

  auto es = expected.strandedness;
  auto eo = expected.orientation;

  auto os = observed.strandedness;
  auto oo = observed.orientation;

  // If the orientations are different, they are incompatible
  if (eo != oo) {
    return false;
  } else { // In this branch, the orientations are always compatible
    return (es == ReadStrandedness::U or es == os);
  }
  // SHOULD NOT GET HERE
  fmt::print(stderr, "WARNING: Could not determine strand compatibility!");
  fmt::print(stderr, "please report this.\n");
  return false;
}

template <typename ExpLib>
void writeAbundancesFromCollapsed(const SalmonOpts& sopt, ExpLib& alnLib,
                                  boost::filesystem::path& fname,
                                  std::string headerComments) {
  using salmon::math::LOG_0;
  using salmon::math::LOG_1;

  // If we're using lightweight-alignment (FMD)
  // and not allowing orphans.
  bool useScaledCounts = !(sopt.useQuasi or sopt.allowOrphans);

  std::unique_ptr<std::FILE, int (*)(std::FILE*)> output(
      std::fopen(fname.c_str(), "w"), std::fclose);

  fmt::print(output.get(), "{}", headerComments);
  fmt::print(output.get(), "Name\tLength\tEffectiveLength\tTPM\tNumReads\n");

  double numMappedFrags = alnLib.upperBoundHits();

  std::vector<Transcript>& transcripts_ = alnLib.transcripts();
  for (auto& transcript : transcripts_) {
    transcript.projectedCounts = useScaledCounts
                                     ? (transcript.mass(false) * numMappedFrags)
                                     : transcript.sharedCount();
  }

  double tfracDenom{0.0};
  for (auto& transcript : transcripts_) {
    double refLength = sopt.noEffectiveLengthCorrection
                           ? transcript.RefLength
                           : std::exp(transcript.getCachedLogEffectiveLength());
    if (sopt.noLengthCorrection) {
      refLength = 100.0;
    }
    tfracDenom += (transcript.projectedCounts / numMappedFrags) / refLength;
  }

  double million = 1000000.0;
  // Now posterior has the transcript fraction
  for (auto& transcript : transcripts_) {
    double logLength = sopt.noEffectiveLengthCorrection
                           ? std::log(transcript.RefLength)
                           : transcript.getCachedLogEffectiveLength();
    double count = transcript.projectedCounts;
    double npm = (transcript.projectedCounts / numMappedFrags);
    double effLength = std::exp(logLength);
    if (sopt.noLengthCorrection) {
      effLength = 100.0;
    }
    double tfrac = (npm / effLength) / tfracDenom;
    double tpm = tfrac * million;
    fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n", transcript.RefName,
               transcript.CompleteLength, effLength, tpm, count);
  }
}

template <typename ExpLib>
void writeAbundances(const SalmonOpts& sopt, ExpLib& alnLib,
                     boost::filesystem::path& fname,
                     std::string headerComments) {
  using salmon::math::LOG_0;
  using salmon::math::LOG_1;

  std::unique_ptr<std::FILE, int (*)(std::FILE*)> output(
      std::fopen(fname.c_str(), "w"), std::fclose);

  fmt::print(output.get(), "{}", headerComments);
  fmt::print(output.get(), "# Name\tLength\tTPM\tFPKM\tNumReads\n");

  auto& refs = alnLib.transcripts();
  auto numMappedFragments = alnLib.numMappedFragments();
  const double logBillion = std::log(1000000000.0);
  const double million = 1000000.0;
  const double logNumFragments =
      std::log(static_cast<double>(numMappedFragments));
  const double upperBoundFactor =
      static_cast<double>(alnLib.upperBoundHits()) / numMappedFragments;

  auto clusters = alnLib.clusterForest().getClusters();
  size_t clusterID = 0;
  for (auto cptr : clusters) {

    // double logClusterMass = cptr->logMass();
    // EDIT
    double logClusterMass = salmon::math::LOG_0;
    double logClusterCount =
        std::log(upperBoundFactor * static_cast<double>(cptr->numHits()));

    bool requiresProjection{false};

    auto& members = cptr->members();
    size_t clusterSize{0};
    for (auto transcriptID : members) {
      Transcript& t = refs[transcriptID];
      t.uniqueCounts = t.uniqueCount();
      t.totalCounts = t.totalCount();
      logClusterMass = salmon::math::logAdd(logClusterMass, t.mass(false));
      ++clusterSize;
    }

    if (logClusterMass == LOG_0) {
      // std::cerr << "Warning: cluster " << clusterID << " has 0 mass!\n";
    }

    for (auto transcriptID : members) {
      Transcript& t = refs[transcriptID];
      double logTranscriptMass = t.mass(false);
      // Try bias
      /*
      double logBias = t.bias();
      logTranscriptMass += t.bias();
      */

      if (logTranscriptMass == LOG_0) {
        t.projectedCounts = 0;
      } else {
        double logClusterFraction = logTranscriptMass - logClusterMass;
        t.projectedCounts = std::exp(logClusterFraction + logClusterCount);
        requiresProjection |=
            t.projectedCounts > static_cast<double>(t.totalCounts) or
            t.projectedCounts < static_cast<double>(t.uniqueCounts);
      }
    }

    if (clusterSize > 1 and requiresProjection) {
      cptr->projectToPolytope(refs);
    }
    ++clusterID;
  }

  auto& transcripts_ = refs;
  double tfracDenom{0.0};
  for (auto& transcript : transcripts_) {
    double refLength = sopt.noEffectiveLengthCorrection
                           ? transcript.RefLength
                           : std::exp(transcript.getCachedLogEffectiveLength());
    // refLength = transcript.RefLength;
    tfracDenom += (transcript.projectedCounts / numMappedFragments) / refLength;
  }

  // Now posterior has the transcript fraction
  for (auto& transcript : transcripts_) {
    double logLength = sopt.noEffectiveLengthCorrection
                           ? std::log(transcript.RefLength)
                           : transcript.getCachedLogEffectiveLength();
    if (sopt.noLengthCorrection) {
      logLength = 1.0;
    }

    // logLength = std::log(transcript.RefLength);
    double fpkmFactor = std::exp(logBillion - logLength - logNumFragments);
    double count = transcript.projectedCounts;
    // double countTotal = transcripts_[transcriptID].totalCounts;
    // double countUnique = transcripts_[transcriptID].uniqueCounts;
    double fpkm = count > 0 ? fpkmFactor * count : 0.0;
    double npm = (transcript.projectedCounts / numMappedFragments);
    double refLength = std::exp(logLength);
    double tfrac = (npm / refLength) / tfracDenom;
    double tpm = tfrac * million;

    fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n", transcript.RefName,
               transcript.RefLength, tpm, fpkm, count);
  }
}

template <typename AlnLibT>
void normalizeAlphas(const SalmonOpts& sopt, AlnLibT& alnLib) {

  using salmon::math::LOG_0;
  using salmon::math::LOG_1;

  auto& refs = alnLib.transcripts();
  auto numMappedFragments = alnLib.numMappedFragments();
  const double logNumFragments =
      std::log(static_cast<double>(numMappedFragments));
  auto clusters = alnLib.clusterForest().getClusters();
  size_t clusterID = 0;
  for (auto cptr : clusters) {

    // double logClusterMass = cptr->logMass();
    // EDIT
    double logClusterMass = salmon::math::LOG_0;
    double logClusterCount = std::log(static_cast<double>(cptr->numHits()));

    bool requiresProjection{false};

    auto& members = cptr->members();
    size_t clusterSize{0};
    for (auto transcriptID : members) {
      Transcript& t = refs[transcriptID];
      t.uniqueCounts = t.uniqueCount();
      t.totalCounts = t.totalCount();
      logClusterMass = salmon::math::logAdd(logClusterMass,
                                            t.mass(false)); // + t.bias());
      ++clusterSize;
    }

    if (logClusterMass == LOG_0) {
      // std::cerr << "Warning: cluster " << clusterID << " has 0 mass!\n";
    }

    for (auto transcriptID : members) {
      Transcript& t = refs[transcriptID];
      double logTranscriptMass = t.mass(false);
      // Try bias
      // double logBias = t.bias();
      // logTranscriptMass += t.bias();

      if (logTranscriptMass == LOG_0) {
        t.projectedCounts = 0;
      } else {
        double logClusterFraction = logTranscriptMass - logClusterMass;
        t.projectedCounts = std::exp(logClusterFraction + logClusterCount);
        requiresProjection |=
            t.projectedCounts > static_cast<double>(t.totalCounts) or
            t.projectedCounts < static_cast<double>(t.uniqueCounts);
      }
    }

    if (clusterSize > 1 and requiresProjection) {
      cptr->projectToPolytope(refs);
    }
    ++clusterID;
  }

  auto& transcripts_ = refs;
  double nFracDenom{0.0};
  for (auto& transcript : transcripts_) {
    nFracDenom += (transcript.projectedCounts / numMappedFragments);
  }

  double invNFracTotal = 1.0 / nFracDenom;
  for (auto& transcript : transcripts_) {
    double v = transcript.projectedCounts / numMappedFragments;
    // transcript.setMass(v * invNFracTotal);
    transcript.setMass(transcript.projectedCounts);
  }
}

LibraryFormat hitType(int32_t end1Start, bool end1Fwd, int32_t end2Start,
                      bool end2Fwd) {

  // If the reads come from opposite strands
  if (end1Fwd != end2Fwd) {
    // and if read 1 comes from the forward strand
    if (end1Fwd) {
      // then if read 1 start < read 2 start ==> ISF
      if (end1Start <= end2Start) {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                             ReadStrandedness::SA);
      } // otherwise read 2 start < read 1 start ==> OSF
      else {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                             ReadStrandedness::SA);
      }
    }
    // and if read 2 comes from the forward strand
    if (end2Fwd) {
      // then if read 2 start <= read 1 start ==> ISR
      if (end2Start <= end1Start) {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                             ReadStrandedness::AS);
      } // otherwise, read 2 start > read 1 start ==> OSR
      else {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                             ReadStrandedness::AS);
      }
    }
  } else {         // Otherwise, the reads come from the same strand
    if (end1Fwd) { // if it's the forward strand ==> MSF
      return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                           ReadStrandedness::S);
    } else { // if it's the reverse strand ==> MSR
      return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                           ReadStrandedness::A);
    }
  }
  // SHOULD NOT GET HERE
  spdlog::get("jointLog")
      ->error("ERROR: Could not associate any known library type with read! "
              "Please report this bug!\n");
  std::exit(-1);
  return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE,
                       ReadStrandedness::U);
}

LibraryFormat hitType(int32_t end1Start, bool end1Fwd, uint32_t len1,
                      int32_t end2Start, bool end2Fwd, uint32_t len2,
                      bool canDovetail) {

  // If the reads come from opposite strands
  if (end1Fwd != end2Fwd) {
    // and if read 1 comes from the forward strand
    if (end1Fwd) {
      // then if read 1 start < read 2 start ==> ISF
      // NOTE: We can't really delineate between inward facing reads that
      // stretch
      // past each other and outward facing reads --- the purpose of stretch is
      // to help
      // make this determinateion.
      int32_t stretch = canDovetail ? len2 : 0;
      if (end1Start <= end2Start + stretch) {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                             ReadStrandedness::SA);
      } // otherwise read 2 start < read 1 start ==> OSF
      else {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                             ReadStrandedness::SA);
      }
    }
    // and if read 2 comes from the forward strand
    if (end2Fwd) {
      // then if read 2 start <= read 1 start ==> ISR
      // NOTE: We can't really delineate between inward facing reads that
      // stretch
      // past each other and outward facing reads --- the purpose of stretch is
      // to help
      // make this determinateion.
      int32_t stretch = canDovetail ? len1 : 0;
      if (end2Start <= end1Start + stretch) {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                             ReadStrandedness::AS);
      } // otherwise, read 2 start > read 1 start ==> OSR
      else {
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                             ReadStrandedness::AS);
      }
    }
  } else {         // Otherwise, the reads come from the same strand
    if (end1Fwd) { // if it's the forward strand ==> MSF
      return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                           ReadStrandedness::S);
    } else { // if it's the reverse strand ==> MSR
      return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                           ReadStrandedness::A);
    }
  }
  // SHOULD NOT GET HERE
  spdlog::get("jointLog")
      ->error("ERROR: Could not associate any known library type with read! "
              "Please report this bug!\n");
  std::exit(-1);
  return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE,
                       ReadStrandedness::U);
}

LibraryFormat hitType(int32_t start, bool isForward) {
  // If the read comes from the forward strand
  if (isForward) {
    return LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                         ReadStrandedness::S);
  } else {
    return LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                         ReadStrandedness::A);
  }
  // SHOULD NOT GET HERE
  fmt::print(stderr,
             "WARNING: Could not associate known library type with read!\n");
  return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE,
                       ReadStrandedness::U);
}

using std::string;
using NameVector = std::vector<string>;
using IndexVector = std::vector<size_t>;
using KmerVector = std::vector<uint64_t>;

/**
 * This function parses the library format string that specifies the format in
 * which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatStringNew(std::string& fmt) {
  using std::vector;
  using std::string;
  using std::map;
  using std::stringstream;

  map<string, LibraryFormat> formatMap = {
      {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                           ReadStrandedness::U)},
      {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                            ReadStrandedness::SA)},
      {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                            ReadStrandedness::AS)},
      {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                           ReadStrandedness::U)},
      {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                            ReadStrandedness::SA)},
      {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                            ReadStrandedness::AS)},
      {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                           ReadStrandedness::U)},
      {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                            ReadStrandedness::S)},
      {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                            ReadStrandedness::A)},
      {"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                          ReadStrandedness::U)},
      {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                           ReadStrandedness::S)},
      {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                           ReadStrandedness::A)}};

  // inspired by
  // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
  // first convert the string to upper-case
  for (auto& c : fmt) {
    c = std::toupper(c);
  }

  auto libFmtIt = formatMap.find(fmt);

  if (libFmtIt == formatMap.end()) {
    stringstream errstr;
    errstr << "unknown library format string : " << fmt;
    throw std::invalid_argument(errstr.str());
  }

  return libFmtIt->second;
}

/**
 * Parses a set of __ordered__ command line options and extracts the relevant
 * read libraries from them.
 */
std::vector<ReadLibrary>
extractReadLibraries(boost::program_options::parsed_options& orderedOptions) {
  // The current (default) format for paired end data
  LibraryFormat peFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                         ReadStrandedness::U);
  // The current (default) format for single end data
  LibraryFormat seFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                         ReadStrandedness::U);

  auto isAutoLibType = [](std::string& fmt) -> bool {
    return (fmt.length() == 1 and (fmt.front() == 'a' or fmt.front() == 'A'));
  };

  auto log = spdlog::get("jointLog");

  bool sawLibType{false};
  bool sawPairedLibrary{false};
  bool sawUnpairedLibrary{false};
  bool autoLibType{false};

  std::vector<ReadLibrary> peLibs{peFormat};
  std::vector<ReadLibrary> seLibs{seFormat};

  for (auto& opt : orderedOptions.options) {
    // Update the library type
    if (opt.string_key == "libType") {
      if (!isAutoLibType(opt.value[0])) {
        auto libFmt = parseLibraryFormatStringNew(opt.value[0]);
        if (libFmt.type == ReadType::PAIRED_END) {
          peFormat = libFmt;
          peLibs.emplace_back(libFmt);
        } else {
          seFormat = libFmt;
          seLibs.emplace_back(libFmt);
        }
      } else {
        autoLibType = true;
      }
      sawLibType = true;
    }

    if (opt.string_key == "mates1") {
      if (!sawLibType) {
        log->warn("Encountered a read file (--mates1/-1) before a "
                  "library type specification.  The (--libType/-l) "
                  "option must precede the input files.");
        peLibs.clear();
        return peLibs;
      }
      peLibs.back().addMates1(opt.value);
      if (autoLibType) {
        peLibs.back().enableAutodetect();
      }
      sawPairedLibrary = true;
    }
    if (opt.string_key == "mates2") {
      if (!sawLibType) {
        log->warn("Encountered a read file (--mates2/-2) before a "
                  "library type specification.  The (--libType/-l) "
                  "option must precede the input files.");
        peLibs.clear();
        return peLibs;
      }
      peLibs.back().addMates2(opt.value);
      if (autoLibType) {
        peLibs.back().enableAutodetect();
      }
      sawPairedLibrary = true;
    }
    if (opt.string_key == "unmatedReads") {
      if (!sawLibType) {
        log->warn("Encountered a read file (--unmatedReads/-r) before a "
                  "library type specification.  The (--libType/-l) "
                  "option must precede the input files.");
        seLibs.clear();
        return seLibs;
      }
      seLibs.back().addUnmated(opt.value);
      if (autoLibType) {
        seLibs.back().enableAutodetect();
      }
      sawUnpairedLibrary = true;
    }
  }

  std::vector<ReadLibrary> libs;

  // @Avi : Allow this temporarily for now, since there is some use to hijack
  // this behavior in Alevin.  However, we should figure out a proper parsing
  // strategy for that rather than abusing single & PE library types.  Once we
  // fix that, we should uncomment the below.
  /*
  if (sawPairedLibrary and sawUnpairedLibrary) {
    log->warn("It seems you have specified both paired-end and unpaired read "
              "libraries.  Salmon does not accepted mixed library types, and "
              "different library types should typically not be quantified together "
              "anyway.  Please quantifiy distinct library types separately.");
    return libs;
  }
  */
  (void)sawPairedLibrary;
  (void)sawUnpairedLibrary;

  libs.reserve(peLibs.size() + seLibs.size());
  for (auto& lib : boost::range::join(seLibs, peLibs)) {
    if (lib.format().type == ReadType::SINGLE_END) {
      if (lib.unmated().size() == 0) {
        // Didn't use default single end library type
        continue;
      }
    } else if (lib.format().type == ReadType::PAIRED_END) {
      if (lib.mates1().size() == 0 or lib.mates2().size() == 0) {
        // Didn't use default paired-end library type
        continue;
      }
    }
    libs.push_back(lib);
  }

  size_t numLibs = libs.size();
  if (numLibs == 1) {
    log->info("There is 1 library.");
  } else if (numLibs > 1) {
    log->info("There are {} libraries.", numLibs);
  }
  return libs;
}

/**
 * This function parses the library format string that specifies the format in
 * which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatString(std::string& fmt) {
  using std::vector;
  using std::string;
  using std::map;
  using std::stringstream;

  // inspired by
  // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c

  // first convert the string to upper-case
  for (auto& c : fmt) {
    c = std::toupper(c);
  }
  // split on the delimiter ':', and put the key, value (k=v) pairs into a map
  stringstream ss(fmt);
  string item;
  map<string, string> kvmap;
  while (std::getline(ss, item, ':')) {
    auto splitPos = item.find('=', 0);
    string key{item.substr(0, splitPos)};
    string value{item.substr(splitPos + 1)};
    kvmap[key] = value;
  }

  map<string, ReadType> readType = {{"SE", ReadType::SINGLE_END},
                                    {"PE", ReadType::PAIRED_END}};
  map<string, ReadOrientation> orientationType = {
      {">>", ReadOrientation::SAME},
      {"<>", ReadOrientation::AWAY},
      {"><", ReadOrientation::TOWARD},
      {"*", ReadOrientation::NONE}};
  map<string, ReadStrandedness> strandType = {{"SA", ReadStrandedness::SA},
                                              {"AS", ReadStrandedness::AS},
                                              {"A", ReadStrandedness::A},
                                              {"S", ReadStrandedness::S},
                                              {"U", ReadStrandedness::U}};
  auto it = kvmap.find("T");
  string typeStr = "";
  if (it != kvmap.end()) {
    typeStr = it->second;
  } else {
    it = kvmap.find("TYPE");
    if (it != kvmap.end()) {
      typeStr = it->second;
    }
  }

  if (typeStr != "SE" and typeStr != "PE") {
    string e = typeStr + " is not a valid read type; must be one of {SE, PE}";
    throw std::invalid_argument(e);
  }

  ReadType type =
      (typeStr == "SE") ? ReadType::SINGLE_END : ReadType::PAIRED_END;
  ReadOrientation orientation = (type == ReadType::SINGLE_END)
                                    ? ReadOrientation::NONE
                                    : ReadOrientation::TOWARD;
  ReadStrandedness strandedness{ReadStrandedness::U};
  // Construct the LibraryFormat class from the key, value map
  for (auto& kv : kvmap) {
    auto& k = kv.first;
    auto& v = kv.second;
    if (k == "O" or k == "ORIENTATION") {
      auto it = orientationType.find(v);
      if (it != orientationType.end()) {
        orientation = orientationType[it->first];
      } else {
        string e =
            v + " is not a valid orientation type; must be one of {>>, <>, ><}";
        throw std::invalid_argument(e);
      }
    }
    if (k == "S" or k == "STRAND") {
      auto it = strandType.find(v);
      if (it != strandType.end()) {
        strandedness = strandType[it->first];
      } else {
        string e =
            v + " is not a valid strand type; must be one of {SA, AS, S, A, U}";
        throw std::invalid_argument(e);
      }
    }
  }
  LibraryFormat lf(type, orientation, strandedness);
  return lf;
}

bool peekBAMIsPaired(const boost::filesystem::path& file) {
  namespace bfs = boost::filesystem;
  std::string readMode = "r";

  if (bfs::is_regular_file(file)) {
    if (bfs::is_empty(file)) {
      fmt::MemoryWriter errstr;
      errstr << "file [" << file.string()
             << "] appears to be empty "
                "(i.e. it has size 0).  This is likely an error. "
                "Please re-run salmon with a corrected input file.\n\n";
      throw std::invalid_argument(errstr.str());
      return false;
    }
  }
  if (file.extension() == ".bam") {
    readMode = "rb";
  }

  auto* fp = scram_open(file.c_str(), readMode.c_str());

  // If we couldn't open the file, then report this and exit.
  if (fp == NULL) {
    fmt::MemoryWriter errstr;
    errstr << "ERROR: Failed to open file " << file.string() << ", exiting!\n";
    throw std::invalid_argument(errstr.str());
    return false;
  }

  bam_seq_t* read = nullptr;
  read = staden::utils::bam_init();

  bool didRead = (scram_get_seq(fp, &read) >= 0);
  bool isPaired{false};

  if (didRead) {
    isPaired = bam_flag(read) & BAM_FPAIRED;
  } else {
    fmt::MemoryWriter errstr;
    errstr << "ERROR: Failed to read alignment from " << file.string()
           << ", exiting!\n";
    staden::utils::bam_destroy(read);
    throw std::invalid_argument(errstr.str());
    return false;
  }

  scram_close(fp);
  staden::utils::bam_destroy(read);
  return isPaired;
}

uint64_t encode(uint64_t tid, uint64_t offset) {
  uint64_t res = (((tid & 0xFFFFFFFF) << 32) | (offset & 0xFFFFFFFF));
  return res;
}

uint32_t transcript(uint64_t enc) {
  uint32_t t = (enc & 0xFFFFFFFF00000000) >> 32;
  return t;
}

uint32_t offset(uint64_t enc) {
  uint32_t o = enc & 0xFFFFFFFF;
  return o;
}

size_t numberOfReadsInFastaFile(const std::string& fname) {
  constexpr size_t bufferSize = 16184;
  char buffer[bufferSize];
  std::ifstream ifile(fname, std::ifstream::in);
  ifile.rdbuf()->pubsetbuf(buffer, bufferSize);

  size_t numReads = 0;
  std::string s;
  while (ifile >> s) {
    if (s.front() == '>') {
      ++numReads;
    }
  }

  ifile.close();

  return numReads;
}

bool readKmerOrder(const std::string& fname, std::vector<uint64_t>& kmers) {

  std::ifstream mlist(fname, std::ios::in | std::ios::binary);
  // Get the number of kmers from file
  size_t numKmers{0};
  mlist.read(reinterpret_cast<char*>(&numKmers), sizeof(size_t));

  // Resize the array that will hold the sorted kmers
  kmers.resize(numKmers, 0);
  mlist.read(reinterpret_cast<char*>(&kmers[0]),
             sizeof(uint64_t) * kmers.size());

  mlist.close();

  return true;
}

template <template <typename> class S, typename T>
bool overlap(const S<T>& a, const S<T>& b) {
  // Query from the smaller set to the larger set
  if (a.size() <= b.size()) {
    for (auto& ae : a) {
      if (b.find(ae) != b.end()) {
        return true;
      }
    }
  } else {
    for (auto& be : b) {
      if (a.find(be) != b.end()) {
        return true;
      }
    }
  }
  // If nothing from the smaller set is in the larger set, then they don't
  // overlap
  return false;
}

TranscriptGeneMap transcriptGeneMapFromGTF(const std::string& fname,
                                           std::string key) {

  using std::unordered_set;
  using std::unordered_map;
  using std::vector;
  using std::tuple;
  using std::string;
  using std::get;

  // Get the logger
  auto logger = spdlog::get("jointLog");

  // Use GffReader to read the file
  GffReader reader(const_cast<char*>(fname.c_str()), true, false);
  // Remember the optional attributes
  reader.readAll(true);

  struct TranscriptKeyPair {
    const char* transcript_id;
    const char* key;
    TranscriptKeyPair(const char* t, const char* k)
        : transcript_id(t), key(k) {}
  };

  // The user can group transcripts by gene_id, gene_name, or
  // an optional attribute that they provide as a string.
  enum class TranscriptKey { GENE_ID, GENE_NAME, DYNAMIC };

  // Select the proper attribute by which to group
  TranscriptKey tkey = TranscriptKey::GENE_ID;

  if (key == "gene_id") {
    // This is the default initalization above.
  } else if (key == "gene_name") {
    tkey = TranscriptKey::GENE_NAME;
  } else {
    tkey = TranscriptKey::DYNAMIC;
  }

  // Iterate over all transcript features and build the
  // transcript <-> key vector.
  auto nfeat = reader.gflst.Count();
  std::vector<TranscriptKeyPair> feats;
  for (int i = 0; i < nfeat; ++i) {
    auto f = reader.gflst[i];
    if (f->isTranscript()) {
      const char* keyStr = nullptr;
      switch (tkey) {
      case TranscriptKey::GENE_ID:
        keyStr = f->getGeneID();
        break;
      case TranscriptKey::GENE_NAME:
        keyStr = f->getGeneName();
        break;
      case TranscriptKey::DYNAMIC:
        keyStr = f->getAttr(key.c_str());
        break;
      }
      if (keyStr != nullptr and keyStr != NULL and f->hasGffID()) {
        feats.emplace_back(f->getID(), keyStr);
      } else {
        if (!f->hasGffID()) {
          logger->warn("Feature has no GFF ID");
        }
        if (keyStr == NULL) {
          const char* fid = f->hasGffID() ? f->getID() : "NO_GFF_ID";
          logger->warn("Could not find key for feature {}", fid);
        }
      }
    }
  }

  // Given the transcript <-> key vector, build the
  // TranscriptGeneMap.

  IndexVector t2g;
  NameVector transcriptNames;
  NameVector geneNames;

  // holds the mapping from transcript ID to gene ID
  IndexVector t2gUnordered;
  // holds the set of gene IDs
  unordered_map<string, size_t> geneNameToID;

  // To read the input and assign ids
  size_t transcriptCounter = 0;
  size_t geneCounter = 0;
  string transcript;
  string gene;

  std::sort(feats.begin(), feats.end(),
            [](const TranscriptKeyPair& a, const TranscriptKeyPair& b) -> bool {
              return std::strcmp(a.transcript_id, b.transcript_id) < 0;
            });

  std::string currentTranscript = "";
  for (auto& feat : feats) {

    std::string gene(feat.key);
    std::string transcript(feat.transcript_id);

    if (transcript != currentTranscript) {
      auto geneIt = geneNameToID.find(gene);
      size_t geneID = 0;

      if (geneIt == geneNameToID.end()) {
        // If we haven't seen this gene yet, give it a new ID
        geneNameToID[gene] = geneCounter;
        geneID = geneCounter;
        geneNames.push_back(gene);
        ++geneCounter;
      } else {
        // Otherwise lookup the ID
        geneID = geneIt->second;
      }

      transcriptNames.push_back(transcript);
      t2g.push_back(geneID);

      //++transcriptID;
      currentTranscript = transcript;
    }
  }

  return TranscriptGeneMap(transcriptNames, geneNames, t2g);
}

TranscriptGeneMap readTranscriptToGeneMap(std::ifstream& ifile) {

  using std::unordered_set;
  using std::unordered_map;
  using std::vector;
  using std::tuple;
  using std::string;
  using std::get;

  using NameID = tuple<string, size_t>;

  IndexVector t2g;
  NameVector transcriptNames;
  NameVector geneNames;

  // holds the transcript name ID mapping
  vector<NameID> transcripts;
  // holds the mapping from transcript ID to gene ID
  IndexVector t2gUnordered;
  // holds the set of gene IDs
  unordered_map<string, size_t> geneNameToID;

  // To read the input and assign ids
  size_t transcriptCounter = 0;
  size_t geneCounter = 0;
  string transcript;
  string gene;

  while (ifile >> transcript >> gene) {
    // The transcript and it's ID
    transcripts.push_back(make_tuple(transcript, transcriptCounter));

    auto geneIt = geneNameToID.find(gene);
    size_t geneID = 0;

    if (geneIt == geneNameToID.end()) {
      // If we haven't seen this gene yet, give it a new ID
      geneNameToID[gene] = geneCounter;
      geneID = geneCounter;
      geneNames.push_back(gene);
      ++geneCounter;
    } else {
      // Otherwise lookup the ID
      geneID = geneIt->second;
    }

    // Map the transcript to the gene in terms of their IDs
    t2gUnordered.push_back(geneID);

    ++transcriptCounter;
  }

  std::sort(transcripts.begin(), transcripts.end(),
            [](const NameID& a, const NameID& b) -> bool {
              return get<0>(a) < get<0>(b);
            });

  // Resize these vectors for fast access
  transcriptNames.resize(t2gUnordered.size());
  t2g.resize(t2gUnordered.size());

  for (size_t newID = 0; newID < transcripts.size(); ++newID) {
    // For each transcript, map it to the appropriate gene
    string oldName;
    size_t oldID;
    std::tie(oldName, oldID) = transcripts[newID];
    t2g[newID] = t2gUnordered[oldID];
    transcriptNames[newID] = oldName;
  }

  return TranscriptGeneMap(transcriptNames, geneNames, t2g);
}

TranscriptGeneMap
transcriptToGeneMapFromFasta(const std::string& transcriptsFile) {
  using std::vector;
  using sequence_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;
  namespace bfs = boost::filesystem;

  NameVector transcriptNames;
  NameVector geneNames{"gene"};

  vector<bfs::path> paths{transcriptsFile};
  size_t maxReadGroupSize{100};
  std::vector<std::string> readFiles{transcriptsFile};
  sequence_parser parser(readFiles, 1, 1, maxReadGroupSize);
  parser.start();

  // while there are transcripts left to process
  auto rg = parser.getReadGroup();
  while (parser.refill(rg)) {
    for (auto& read : rg) {
      // The transcript name
      std::string fullHeader(read.name);
      std::string header = fullHeader.substr(0, fullHeader.find(' '));
      transcriptNames.emplace_back(header);
    }
  }

  // Sort the transcript names
  std::sort(transcriptNames.begin(), transcriptNames.end());

  // Since we have no real gene groupings, the t2g vector is trivial,
  // everything maps to gene 0.
  IndexVector t2g(transcriptNames.size(), 0);

  return TranscriptGeneMap(transcriptNames, geneNames, t2g);
}

class ExpressionRecord {
public:
  ExpressionRecord(const std::string& targetIn, uint32_t lengthIn,
                   double effLengthIn, std::vector<double>& expValsIn)
      : target(targetIn), length(lengthIn), effLength(effLengthIn),
        expVals(expValsIn) {}

  ExpressionRecord(ExpressionRecord&& other) {
    std::swap(target, other.target);
    length = other.length;
    effLength = other.effLength;
    std::swap(expVals, other.expVals);
  }

  ExpressionRecord(std::vector<std::string>& inputLine) {
    if (inputLine.size() < 3) {
      std::string err("Any expression line must contain at least 3 tokens");
      throw std::invalid_argument(err);
    } else {
      auto it = inputLine.begin();
      target = *it;
      ++it;
      length = std::stoi(*it);
      ++it;
      effLength = std::stod(*it);
      ++it;
      for (; it != inputLine.end(); ++it) {
        expVals.push_back(std::stod(*it));
      }
    }
  }

  std::string target;
  uint32_t length;
  double effLength;
  std::vector<double> expVals;
};

// From : http://stackoverflow.com/questions/9435385/split-a-string-using-c11
std::vector<std::string> split(const std::string& str,
                               int delimiter(int) = ::isspace) {
  using namespace std;
  vector<string> result;
  auto e = str.end();
  auto i = str.begin();
  while (i != e) {
    i = find_if_not(i, e, delimiter);
    if (i == e)
      break;
    auto j = find_if(i, e, delimiter);
    result.push_back(string(i, j));
    i = j;
  }
  return result;
}

std::string getCurrentTimeAsString() {
  // Get the current time as a string
  std::time_t result = std::time(NULL);

  // old non-threadsafe version
  //auto time = std::string(std::asctime(std::localtime(&result)));
  //time.pop_back(); // remove the newline
  //return time;

  struct tm local_tm;
  // NOTE: localtime_r may not exist on windows systems.  This is OK
  // as salmon doesn't support windows 
  ::localtime_r(&result, &local_tm);
  char buffer[80] = {0};
  std::strftime(buffer, sizeof(buffer),"%a %b %d %H:%M:%S %Y", &local_tm);
  std::string str(buffer);
  return str;
}

  bool validateOptionsAlignment_(
                                 SalmonOpts& sopt,
                                 boost::program_options::variables_map& vm
                                 ) {
  if (!sopt.sampleOutput and sopt.sampleUnaligned) {
    sopt.jointLog->warn(
        "You passed in the (-u/--sampleUnaligned) flag, but did not request a "
        "sampled "
        "output file (-s/--sampleOut).  This flag will be ignored!");
  }

  if (sopt.useErrorModel and sopt.rangeFactorizationBins < 4) {
    uint32_t nbins{4};
    sopt.jointLog->info(
                        "Usage of --useErrorModel implies use of range factorization. "
                        "rangeFactorization bins is being set to {}", nbins
                        );
    sopt.rangeFactorizationBins = nbins;
    sopt.useRangeFactorization = true;
  }
  return true;
}

  bool validateOptionsMapping_(
                               SalmonOpts& sopt,
                               boost::program_options::variables_map& vm
                               ) {
  auto numUnpaired = sopt.unmatedReadFiles.size();
  auto numLeft = sopt.mate1ReadFiles.size();
  auto numRight = sopt.mate2ReadFiles.size();

  /** Info, not warnings or errors, but informative messages to the user **/
  if (sopt.mimicBT2 or sopt.mimicStrictBT2 or sopt.hardFilter) {
    sopt.jointLog->info("The --mimicBT2, --mimicStrictBT2 and --hardFilter flags imply mapping validation (--validateMappings). "
                        "Enabling mapping validation.");
    sopt.validateMappings = true;
  }

  // If not in alevin mode, inform the user about validateMappings
  if (!sopt.alevinMode and !sopt.validateMappings) {
    sopt.jointLog->warn("\n\n"
                        "NOTE: It appears you are running salmon without the `--validateMappings` option.\n"
                        "Mapping validation can generally improve both the sensitivity and specificity of mapping,\n"
                        "with only a moderate increase in use of computational resources. \n"
                        "Mapping validation is planned to become a default option (i.e. turned on by default) in\n"
                        "the next release of salmon.\n"
                        "Unless there is a specific reason to do this (e.g. testing on clean simulated data),\n"
                        "`--validateMappings` is generally recommended.\n");
  }

  bool is_pe_library = (numLeft + numRight > 0);
  bool is_se_library = (numUnpaired > 0);

  // currently there is some strange use for this in alevin, I think ...
  // check with avi.
  if (is_pe_library and is_se_library) {
      sopt.jointLog->warn("You seem to have passed in both un-paired reads and paired-end reads. "
                          "It is not currently possible to quantify hybrid library types in salmon.");
  }

  if (is_pe_library) {
    if (numLeft != numRight) {
      sopt.jointLog->error("You passed paired-end files to salmon, but you passed {} files to --mates1 "
                           "and {} files to --mates2.  You must pass the same number of files to both flags",
                           numLeft, numRight);
      return false;
    }
  } 

  auto checkScoreValue = [&sopt](int16_t score, std::string sname) -> bool {
                           using score_t = int8_t;
                           auto minval = static_cast<int16_t>(std::numeric_limits<score_t>::min());
                           auto maxval = static_cast<int16_t>(std::numeric_limits<score_t>::max());
                           if (score <  minval or score > maxval) {
                             sopt.jointLog->error("You set the {} as {}, but it must be in "
                                                  "the range [{}, {}].", sname, score, minval, maxval);
                             return false;
                           }
                           return true;
                         };

  if(!checkScoreValue(sopt.matchScore, "match score")) { return false; }
  if(!checkScoreValue(sopt.mismatchPenalty, "mismatch penalty")) { return false; }
  if(!checkScoreValue(sopt.gapOpenPenalty, "gap open penalty")) { return false; }
  if(!checkScoreValue(sopt.gapExtendPenalty, "gap extend penalty")) { return false; }

  if (sopt.mismatchPenalty > 0) {
    sopt.jointLog->warn(
                        "You set the mismatch penalty as {}, but it should be negative.  It is being negated to {}.",
                        sopt.mismatchPenalty, -sopt.mismatchPenalty);
    sopt.mismatchPenalty = -sopt.mismatchPenalty;
  }

  // Make sure that consensusSlack is not negative
  if (sopt.consensusSlack < 0 or sopt.consensusSlack >= 1.0) {
    sopt.jointLog->error("You set consensusSlack as {}, but it must in [0,1).", sopt.consensusSlack);
    return false;
  }

  if (sopt.mismatchSeedSkip < 1) {
    sopt.jointLog->warn("The mismatchSeedSkip was set to {}, but it cannot be < 1.  Setting mismatchSeedSkip to 1");
    sopt.mismatchSeedSkip = 1;
  }

  if (sopt.mismatchSeedSkip > 31) {
    sopt.jointLog->warn("Setting the mismatchSeedSkip too high can hurt the sensitivity of mapping.  Consider "
    "setting a lower mismatchSeedSkip.");
  }

  // If we have validate mappings, then make sure we automatically enable
  // range factorization
  if (sopt.validateMappings) {
    if (!vm.count("minScoreFraction")) {
      sopt.minScoreFraction = salmon::defaults::minScoreFraction;
      sopt.jointLog->info(
                          "Usage of --validateMappings implies use of minScoreFraction. "
                          "Since not explicitly specified, it is being set to {}", sopt.minScoreFraction
                          );
    }

    if (sopt.hardFilter) {
      // range factorization doesn't make sense with hard filtering
      if (sopt.rangeFactorizationBins > 0) {
        sopt.jointLog->info("The use of range-factorized equivalence classes does not make sense "
                            "in conjunction with --hardFilter.  Disabling range-factorized equivalence classes. ");
        sopt.rangeFactorizationBins = 0;
        sopt.useRangeFactorization = false;
      }
    } else if (sopt.rangeFactorizationBins < 4) {
      uint32_t nbins{4};
      sopt.jointLog->info(
                          "Usage of --validateMappings, without --hardFilter implies use of range factorization. "
                          "rangeFactorizationBins is being set to {}", nbins
                          );
      sopt.rangeFactorizationBins = nbins;
      sopt.useRangeFactorization = true;
    }

    // If the consensus slack was not set explicitly, then it defaults to 0.35 with
    // validateMappings
    bool consensusSlackExplicit = !vm["consensusSlack"].defaulted();
    if (!consensusSlackExplicit) {
      sopt.jointLog->info(
                          "Setting consensusSlack to selective-alignment default of {}.", 
                          sopt.consensusSlack);
    }

    bool pre_merge_chain_sub_thresh_explicit = !vm["preMergeChainSubThresh"].defaulted();
    bool post_merge_chain_sub_thresh_explicit = !vm["postMergeChainSubThresh"].defaulted();
    bool orphan_chain_sub_thresh_explicit = !vm["orphanChainSubThresh"].defaulted();

    // for a single-end library (or effectively so by being single-cell), we set 
    // pre_merge_chain_sub_thresh to 1.0 by default
    if ( is_se_library or sopt.alevinMode ) {

      // The default of preMergeChainSubThresh for single-end libraries is 1.0, so set that here
      if (!pre_merge_chain_sub_thresh_explicit) {
        sopt.pre_merge_chain_sub_thresh = 1.0;
      }

      // for single-end libraries, postMergeChainSubThresh and orphanChainSubThresh are meaningless 
      if (post_merge_chain_sub_thresh_explicit) {
        sopt.jointLog->warn("The postMergeChainSubThresh is not meaningful for single-end "
        "(or effectively single-end — e.g. tagged-end single-cell) libraries.  Setting this value "
        "to 1.0 and ignoring");
      }
      if (orphan_chain_sub_thresh_explicit) {
        sopt.jointLog->warn("The orphanChainSubThresh is not meaningful for single-end "
        "(or effectively single-end — e.g. tagged-end single-cell) libraries.  Setting this value "
        "to 1.0 and ignoring");
      }
      sopt.post_merge_chain_sub_thresh = 1.0;
      sopt.orphan_chain_sub_thresh = 1.0;
    }

    // value range check for filters
    // pre-merge
    if (sopt.pre_merge_chain_sub_thresh < 0 or sopt.pre_merge_chain_sub_thresh > 1.0) {
      sopt.jointLog->error("You set preMergeChainSubThresh as {}, but it must in [0,1].", 
        sopt.pre_merge_chain_sub_thresh);
      return false;
    }
    // post-merge
    if (sopt.post_merge_chain_sub_thresh < 0 or sopt.post_merge_chain_sub_thresh > 1.0) {
      sopt.jointLog->error("You set postMergeChainSubThresh as {}, but it must in [0,1].", 
        sopt.post_merge_chain_sub_thresh);
      return false;
    }
    // orphan
    if (sopt.orphan_chain_sub_thresh < 0 or sopt.orphan_chain_sub_thresh > 1.0) {
      sopt.jointLog->error("You set orphanChainSubThresh as {}, but it must in [0,1].", 
        sopt.orphan_chain_sub_thresh);
      return false;
    }


    if (sopt.mimicBT2 and sopt.mimicStrictBT2) {
      sopt.jointLog->error("You passed both the --mimicBT2 and --mimicStrictBT2 parameters.  These are mutually exclusive. "
                            "Please select only one of these flags.");
      return false;
    }

    if (sopt.mimicBT2 or sopt.mimicStrictBT2) {
      /*
      sopt.jointLog->info("The --mimicBT2 and --mimicStrictBT2 flags imply orphan recovery (--recoverOrphans). "
                          "Enabling orphan recovery.");
      sopt.recoverOrphans = true;
      */

      sopt.maxReadOccs = 1000;
      sopt.jointLog->info("The --mimicBT2 and --mimicStrictBT2 flags increases maxReadOccs to {}.", sopt.maxReadOccs);

      sopt.consensusSlack = 0.5;
      sopt.jointLog->info("The --mimicBT2 and --mimicStrictBT2 flags increases consensusSlack to {}.", sopt.consensusSlack);

      if (sopt.mimicBT2) {
        sopt.jointLog->info(
                            "Usage of --mimicBT2 overrides other settings for mapping validation. Setting "
                            "Bowtie2-like parameters now.");
        sopt.discardOrphansQuasi = true;
        if (sopt.softclipOverhangs) {
          sopt.jointLog->info("Softclipping of overhangs is not allowed in mimicBT2 mode; setting to false.");
          sopt.softclipOverhangs = false;
        }

        sopt.matchScore = 2;
        sopt.mismatchPenalty = -4;
        sopt.gapOpenPenalty = 5;
        sopt.gapExtendPenalty = 3;
        /*
        sopt.matchScore = 0;
        sopt.mismatchPenalty = -6;
        sopt.gapOpenPenalty = 5;
        sopt.gapExtendPenalty = 3;
        */
      }

      if (sopt.mimicStrictBT2) {
        sopt.jointLog->info(
                            "Usage of --mimicStrictBT2 overrides other settings for mapping validation. Setting "
                            "strict RSEM+Bowtie2-like parameters now.");
        if (sopt.softclipOverhangs) {
          sopt.jointLog->info("Softclipping of overhangs is not allowed in mimicStrictBT2 mode; setting to false.");
          sopt.softclipOverhangs = false;
        }
        sopt.discardOrphansQuasi = true;
        sopt.minScoreFraction = 0.8;
        sopt.matchScore = 1;
        sopt.mismatchPenalty = 0;
        // NOTE: as a limitation of ksw2, we can't have
        // (gapOpenPenalty + gapExtendPenalty) * 2 + matchScore < numeric_limits<int8_t>::max()
        // these parameters below are sufficiently large penalties to
        // prohibit gaps, while not overflowing the above condition
        sopt.gapOpenPenalty = 25;
        sopt.gapExtendPenalty = 25;
      }
    }
  }
  return true;
}

/**
 * In mapping mode, depending on what the user has requested, we may have to
 * write out some files.  Prepare loggers so we can do this asynchronously.
 **/
bool createAuxMapLoggers_(SalmonOpts& sopt,
                          boost::program_options::variables_map& vm) {
  using std::cerr;
  using std::vector;
  using std::string;
  namespace bfs = boost::filesystem;

  auto jointLog = sopt.jointLog;

  // Create the file (and logger) for outputting unmapped reads, if the user has
  // asked for it.
  if (sopt.writeUnmappedNames) {
    boost::filesystem::path auxDir = sopt.outputDirectory / sopt.auxDir;
    bool auxSuccess = bfs::exists(auxDir) and bfs::is_directory(auxDir);
    if (!auxSuccess) {
      return false;
    }
    bfs::path unmappedNameFile = auxDir / "unmapped_names.txt";
    std::ofstream* outFile = new std::ofstream(unmappedNameFile.string());
    // Make sure file opened successfully.
    if (!outFile->is_open()) {
      jointLog->error("Could not create file for unmapped read names [{}]",
                      unmappedNameFile.string());
      delete outFile;
      return false;
    }
    auto outputSink =
        std::make_shared<spdlog::sinks::ostream_sink_mt>(*outFile);

    std::shared_ptr<spdlog::logger> outLog =
        std::make_shared<spdlog::logger>("unmappedLog", outputSink);
    spdlog::register_logger(outLog);
    outLog->set_pattern("%v");
    sopt.unmappedFile.reset(outFile);
    sopt.unmappedLog = outLog;
  }

  // Create the file (and logger) for outputting unmapped reads, if the user has
  // asked for it.
  if (sopt.writeOrphanLinks) {
    boost::filesystem::path auxDir = sopt.outputDirectory / sopt.auxDir;
    bool auxSuccess = bfs::exists(auxDir) and bfs::is_directory(auxDir);
    if (!auxSuccess) {
      return false;
    }
    bfs::path orphanLinkFile = auxDir / "orphan_links.txt";
    std::ofstream* outFile = new std::ofstream(orphanLinkFile.string());
    // Make sure file opened successfully.
    if (!outFile->is_open()) {
      jointLog->error("Could not create file for orphan links [{}]",
                      orphanLinkFile.string());
      delete outFile;
      return false;
    }

    auto outputSink =
        std::make_shared<spdlog::sinks::ostream_sink_mt>(*outFile);

    std::shared_ptr<spdlog::logger> outLog =
        std::make_shared<spdlog::logger>("orphanLinkLog", outputSink);
    spdlog::register_logger(outLog);
    outLog->set_pattern("%v");
    sopt.orphanLinkFile.reset(outFile);
    sopt.orphanLinkLog = outLog;
  }

  // Determine what we'll do with selective-alignment results
  bool writeQuasimappings = (sopt.qmFileName != "");

  if (writeQuasimappings) {
    std::streambuf* qmBuf{nullptr};
    // output to stdout
    if (sopt.qmFileName == "-") {
      qmBuf = std::cout.rdbuf();
    } else { // output to the requested path, making the directory if it doesn't
             // exist
      // get the absolute file path
      sopt.qmFileName = boost::filesystem::absolute(sopt.qmFileName).string();
      // get the parent directory
      bfs::path qmDir = boost::filesystem::path(sopt.qmFileName).parent_path();
      // if it's not already a directory that exists
      bool qmDirSuccess = boost::filesystem::is_directory(qmDir);
      // try to create it
      if (!qmDirSuccess) {
        qmDirSuccess = boost::filesystem::create_directories(qmDir);
      }
      // if the directory already existed, or we created it successfully, open
      // the file
      if (qmDirSuccess) {
        sopt.qmFile.open(sopt.qmFileName);
        // Make sure file opened successfully.
        if (!sopt.qmFile.is_open()) {
          jointLog->error(
              "Could not create file for writing selective-alignments [{}]",
              sopt.qmFileName);
          return false;
        }
        qmBuf = sopt.qmFile.rdbuf();
      } else {
        bfs::path qmFileName =
            boost::filesystem::path(sopt.qmFileName).filename();
        jointLog->error("Couldn't create requested directory {} in which "
                        "to place the mapping output {}",
                        qmDir.string(), qmFileName.string());
        return false;
      }
    }
    // Now set the output stream to the buffer, which is
    // either std::cout, or a file.
    sopt.qmStream.reset(new std::ostream(qmBuf));

    auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(
        *(sopt.qmStream.get()));
    sopt.qmLog = std::make_shared<spdlog::logger>("qmStream", outputSink);
    sopt.qmLog->set_pattern("%v");
  }
  return true;
}

/**
 *  Try to create the directory named by `dirPath`.  If it already
 *  exists, make sure it is a directory.  If it already exists
 *  but is not a directory, or if we are unable to create it, then
 *  print an appropriate message to stderr.  This function returns
 *  true if the directory was created successfully or already existed, and
 *  false otherwise.
 **/
bool createDirectoryVerbose_(boost::filesystem::path& dirPath) {
  namespace bfs = boost::filesystem;
  if (bfs::exists(dirPath)) {
    // If it already exists and isn't a directory, then complain
    if (!bfs::is_directory(dirPath)) {
      fmt::print(stderr,
                 "{}ERROR{}: Path [{}] already exists "
                 "and is not a directory.\n"
                 "Please either remove this file or choose another "
                 "auxiliary directory.\n",
                 ioutils::SET_RED, ioutils::RESET_COLOR, dirPath);
      return false;
    }
  } else { // If the path doesn't exist, then create it
    if (!bfs::create_directories(dirPath)) { // creation failed for some reason
      fmt::print(stderr,
                 "{}ERROR{}: Could not create the directory [{}]. "
                 "Please check that doing so is valid.",
                 ioutils::SET_RED, ioutils::RESET_COLOR, dirPath);
      return false;
    }
  }
  return true;
}


/* Function used to check that 'opt1' and 'opt2' are not specified
    at the same time. */
// taken from : https://www.boost.org/doc/libs/1_67_0/libs/program_options/example/real.cpp
void conflicting_options(const boost::program_options::variables_map& vm,
                          const char* opt1, const char* opt2){
  if (vm.count(opt1) && !vm[opt1].defaulted()
      && vm.count(opt2) && !vm[opt2].defaulted()) {
    throw std::logic_error(std::string("Conflicting options '")
                           + opt1 + "' and '" + opt2 + "'.");
  }
}

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
// taken from : https://www.boost.org/doc/libs/1_67_0/libs/program_options/example/real.cpp
void option_dependency(const boost::program_options::variables_map& vm,
                       const char* for_what, const char* required_option){
  if (vm.count(for_what) && !vm[for_what].defaulted())
    if (vm.count(required_option) == 0 || vm[required_option].defaulted())
      throw std::logic_error(string("Option '") + for_what 
                             + "' requires option '" + required_option + "'.");
}

/**
 * Validate the options for salmon, and create the necessary
 * output directories and logging infrastructure.
 **/
bool processQuantOptions(SalmonOpts& sopt,
                         boost::program_options::variables_map& vm,
                         int32_t numBiasSamples) {
  using std::cerr;
  using std::vector;
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;
  // Set the atomic variable numBiasSamples from the local version
  sopt.numBiasSamples.store(numBiasSamples);

  // Get the time at the start of the run
  sopt.runStartTime = getCurrentTimeAsString();

  // Verify the geneMap before we start doing any real work.
  bfs::path geneMapPath;
  if (vm.count("geneMap")) {
    // Make sure the provided file exists
    geneMapPath = vm["geneMap"].as<std::string>();
    if (!bfs::exists(geneMapPath)) {
      std::cerr << "ERROR: Could not find transcript <=> gene map file "
                << geneMapPath << "\n";
      std::cerr << "Exiting now: please either omit the \'geneMap\' option or "
                   "provide a valid file\n";
      return false;
    }
    sopt.geneMapPath = geneMapPath;
  }

  /**
   * Create some necessary directories
   **/

  // output directory
  bfs::path outputDirectory(vm["output"].as<std::string>());
  bool outputDirOK = createDirectoryVerbose_(outputDirectory);
  // set the output directory
  if (!outputDirOK) {
    return false;
  }
  sopt.outputDirectory = outputDirectory;

  // log directory
  bfs::path logDirectory = outputDirectory / "logs";
  bool logDirOK = createDirectoryVerbose_(logDirectory);
  if (!logDirOK) {
    return false;
  }
  if (!sopt.quiet) {
    std::cerr << "Logs will be written to " << logDirectory.string() << "\n";
  }

  // parameter directory
  bfs::path paramsDir = outputDirectory / "libParams";
  bool paramDirOK = createDirectoryVerbose_(paramsDir);
  if (!paramDirOK) {
    return false;
  }
  sopt.paramsDirectory = paramsDir;

  // auxiliary directory
  bfs::path auxDir = sopt.outputDirectory / sopt.auxDir;
  bool auxDirOK = createDirectoryVerbose_(auxDir);
  if (!auxDirOK) {
    return false;
  }
  /** Done creating directories **/

  // Metagenomic option
  if (sopt.meta) {
    sopt.initUniform = true;
    sopt.noRichEqClasses = true;
    // for now, meta mode uses the EM.
    sopt.useEM = true;
    sopt.useVBOpt = false;
    // sopt.incompatPrior = salmon::math::LOG_0;
    // sopt.ignoreIncompat = true;
  }

  //  Size for the logger buffer
  size_t max_q_size = 131072; // 2097152;
  bfs::path logPath = logDirectory / "salmon_quant.log";

  if (sopt.quantMode == SalmonQuantMode::MAP) {
    bfs::path indexDirectory(vm["index"].as<string>());
    sopt.indexDirectory = indexDirectory;

    // Determine what we'll do with selective-alignment results
    bool writeQuasimappings = (sopt.qmFileName != "");

    // make it larger if we're writing mappings or
    // unmapped names.
    if (writeQuasimappings or sopt.writeUnmappedNames or
        sopt.writeOrphanLinks) {
      max_q_size = 2097152; // 4194304;//16777216;
    }
  }

  spdlog::set_async_mode(max_q_size);
  auto fileSink =
    std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
  // auto rawConsoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
  // auto consoleSink =
  //    std::make_shared<spdlog::sinks::ansicolor_sink>(rawConsoleSink);
  auto consoleSink =
      std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>();
  consoleSink->set_color(spdlog::level::warn, consoleSink->magenta);
  auto consoleLog = spdlog::create("stderrLog", {consoleSink});
  auto fileLog = spdlog::create("fileLog", {fileSink});
  std::vector<spdlog::sink_ptr> sinks{consoleSink, fileSink};
  auto jointLog =
      spdlog::create("jointLog", std::begin(sinks), std::end(sinks));

  // If we're being quiet, then only emit errors.
  if (sopt.quiet) {
    jointLog->set_level(spdlog::level::err);
  }

  sopt.jointLog = jointLog;
  sopt.fileLog = fileLog;

  if (sopt.quantMode == SalmonQuantMode::MAP) {
    bool auxLoggersOK = createAuxMapLoggers_(sopt, vm);
    if (!auxLoggersOK) {
      return auxLoggersOK;
    }
  }

  /**
   *  SETTING CONDITIONAL OPTIONS
   *
   *
   **/

  // If the user is using range factorization, set the appropriate parameters
  // here
  if (sopt.rangeFactorizationBins > 0) {
    sopt.useRangeFactorization = true;
  }

  // If the user is enabling *just* GC bias correction
  // i.e. without seq-specific bias correction, then disable
  // the conditional model.
  if (sopt.gcBiasCorrect and !sopt.biasCorrect) {
    sopt.numConditionalGCBins = 1;
  }

  {
    std::transform(sopt.hitFilterPolicyStr.begin(), sopt.hitFilterPolicyStr.end(),
                   sopt.hitFilterPolicyStr.begin(), ::toupper);
    if ( sopt.hitFilterPolicyStr == "BEFORE" ) {
      sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::FILTER_BEFORE_CHAINING;
    } else if ( sopt.hitFilterPolicyStr == "AFTER" ) {
      sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::FILTER_AFTER_CHAINING;
    } else if ( sopt.hitFilterPolicyStr == "BOTH" ) {
      sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::FILTER_BEFORE_AND_AFTER_CHAINING;
    } else if ( sopt.hitFilterPolicyStr == "NONE" ) {
      sopt.hitFilterPolicy = pufferfish::util::HitFilterPolicy::DO_NOT_FILTER;
    } else {
      jointLog->critical("The argument {} for --hitFilterPolicy is invalid. Valid options are "
                         "BEFORE, AFTER, BOTH and NONE.", sopt.hitFilterPolicyStr);
      jointLog->flush();
      return false;
    }
  }

  // The growing list of thou shalt nots
  {
    try {
      conflicting_options(vm, "validateMappings", "noSA");
      conflicting_options(vm, "mimicBT2", "noSA");
      conflicting_options(vm, "mimicStrictBT2", "noSA");
      conflicting_options(vm, "hardFilter", "noSA");
    } catch (std::logic_error& e) {
      jointLog->critical(e.what());
      jointLog->flush();
      return false;
    }

    // currently, with the pufferfish based index
    // we have not tested SA-free methods.  So let's
    // disable this ability and give a warning.
    if (sopt.disableSA) {
      jointLog->critical("Note: Alignment-free mapping (i.e. mapping without subsequent selective-alignment) "
                         "has not yet been throughly tested under the pufferfish-based index and using the "
                         "pufferfish-based mapping strategies.  Thus, disabling of selective-alignment "
                         "is not currently allowed.  We may, potentially explore re-enabling this option in future "
                         "versions of salmon.");
      jointLog->flush();
      return false;
      //sopt.validateMappings = false;
    }
  }

  {
    try {
      conflicting_options(vm, "useVBOpt", "useEM");
    } catch (std::logic_error& e) {
      jointLog->critical(e.what());
      jointLog->flush();
      return false;
    }
    // If the user passed useEM, but not useVBOpt, then
    // turn off VB.  The fact that there is not a better
    // way to handle this suggests a potential shortcoming
    // of boost::program_options.
    if(sopt.useEM) {
      sopt.useVBOpt = false;
    }
  }

  {
    try {
      conflicting_options(vm, "perNucleotidePrior", "perTranscriptPrior");
    } catch (std::logic_error& e) {
      jointLog->critical(e.what());
      jointLog->flush();
      return false;
    }
    // If the user passed perNucleotidePrior, then
    // turn off perTranscriptPrior.
    if(sopt.perNucleotidePrior) {
      sopt.perTranscriptPrior = false;
    }
  }

  // If we are in vbOpt mode, and using a perNucleotidePrior, and the user didn't change the default
  // vbPrior, then divide this value by 1000.
  if (sopt.useVBOpt and sopt.perNucleotidePrior and vm["vbPrior"].defaulted()) {
    sopt.vbPrior = 1e-5;
    jointLog->info("Using per-nucleotide prior with the default VB prior.  Setting the default prior to {}",sopt.vbPrior);
  }

  // If the maxHashResizeThreads was defaulted, then set it equal to the regular number 
  // of threads.
  if (vm["maxHashResizeThreads"].defaulted()) {
    sopt.maxHashResizeThreads = sopt.numThreads;
    jointLog->info("setting maxHashResizeThreads to {}", sopt.maxHashResizeThreads);
  }

  {
    try {
      conflicting_options(vm, "alignments", "eqclasses");
    } catch (std::logic_error& e) {
      jointLog->critical(e.what());
      jointLog->flush();
      return false;
    }
  }

  {
    try {
      option_dependency(vm, "alignments", "targets");
    } catch (std::logic_error& e) {
      jointLog->critical(e.what());
      jointLog->flush();
      return false;
    }
  }

  /** Warnings, not errors **/
  {
    if (sopt.numBurninFrags < sopt.numPreBurninFrags) {
      sopt.jointLog->warn("You set the number of burnin fragments "
                          "(--numAuxModelSamples) to be less than the number "
                          "of \n"
                          "pre-burnin fragments (--numPreAuxModelSamples), but "
                          "it must be at least as large.  The \n"
                          "number of pre-burnin fragments and burnin fragments "
                          "is being set to the same value "
                          "({})",
                          sopt.numBurninFrags);
      sopt.numPreBurninFrags = sopt.numBurninFrags;
    }

    // maybe arbitrary, but if it's smaller than this, consider it
    // equal to LOG_0.
    if (sopt.incompatPrior < 1e-100 or sopt.incompatPrior == 0.0) {
      jointLog->info("Fragment incompatibility prior below threshold.  "
                     "Incompatible fragments will be ignored.");
      sopt.incompatPrior = salmon::math::LOG_0;
      sopt.ignoreIncompat = true;
    } else {
      sopt.incompatPrior = std::log(sopt.incompatPrior);
      sopt.ignoreIncompat = false;
    }

    // Dumping equivalnce class weights implies dumping equivalence classes
    if (sopt.dumpEqWeights and !sopt.dumpEq) {
      sopt.dumpEq = true;
      jointLog->info("You specified --dumpEqWeights, which implies --dumpEq; "
                     "that option has been enabled.");
    }
  }

  /** Errors -- will prevent Salmon from running **/
  {
    if (sopt.numGibbsSamples > 0 and sopt.numBootstraps > 0) {
      jointLog->critical(
          "You cannot perform both Gibbs sampling and bootstrapping. "
          "Please choose one.");
      jointLog->flush();
      return false;
    }
    if (sopt.numGibbsSamples > 0) {
      if (!(sopt.thinningFactor >= 1)) {
        jointLog->critical(
            "The Gibbs sampling thinning factor (--thinningFactor) "
            "cannot be smaller than 1.");
        jointLog->flush();
        return false;
      }
    }

    if (sopt.noFragLengthDist and !sopt.noEffectiveLengthCorrection) {
      jointLog->critical(
          "You cannot enable --noFragLengthDist without "
          "also enabling --noEffectiveLengthCorrection; exiting.");
      jointLog->flush();
      return false;
    }

    if (sopt.noLengthCorrection) {
      bool anyBiasCorrect =
          sopt.gcBiasCorrect or sopt.biasCorrect or sopt.posBiasCorrect;
      if (anyBiasCorrect) {
        jointLog->critical(
            "Since bias correction relies on modifying "
            "effective lengths, you cannot enable bias "
            "correction simultaneously with the --noLengthCorrection "
            "option.");
        jointLog->flush();
        return false;
      }
    }

    if (sopt.dontExtrapolateCounts) { // If the user has provided this option,
                                      // (s)he must be using Gibbs sampling
      if (sopt.numGibbsSamples == 0) {
        sopt.jointLog->critical("You passed the --noExtrapolateCounts flag, "
                                "but are not using Gibbs sampling. "
                                "The fomer implies the latter.  Please enable "
                                "Gibbs sampling to use this flag.");
        return false;
      }
    }
  }
  /** End of generic Error validation **/

  // Validation that is different for alignment and mapping based modes.
  bool perModeValidate{true};
  if (sopt.quantMode == SalmonQuantMode::ALIGN) {
    perModeValidate = validateOptionsAlignment_(sopt, vm);
  } else if (sopt.quantMode == SalmonQuantMode::MAP) {
    perModeValidate = validateOptionsMapping_(sopt, vm);
  }

  return perModeValidate;
}

// TODO: Check the use-case of this.  If we still want to support it, then update
// it to read from a potentially gzipped equivalence class file.  Also, the eq file 
// is in a non-standard format (effective lengths at the end), so get this info some
// other way.
bool readEquivCounts(boost::filesystem::path& eqFilePathString,
                     std::vector<string>& tnames,
                     std::vector<double>& tefflens,
                     std::vector<std::vector<uint32_t>>& eqclasses,
                     std::vector<std::vector<double>>& auxs_vals,
                     std::vector<uint32_t>& eqclass_counts ) {

  auto l = spdlog::get("jointLog");

  namespace bfs = boost::filesystem;
  bfs::path eqFilePath {eqFilePathString};

  std::unique_ptr<std::istream> equivFilePtr = nullptr;
  if (eqFilePath.extension() == ".gz") {
    equivFilePtr.reset(new zstr::ifstream(eqFilePath.string()));
  } else {
    equivFilePtr.reset(new std::ifstream(eqFilePath.string()));
  }
  //std::ifstream equivFile(eqFilePath.string());

  size_t numTxps, numEqClasses;
  // Number of transcripts
  (*equivFilePtr) >> numTxps;

  // Number of equivalence classes
  (*equivFilePtr) >> numEqClasses;

  tnames.reserve(numTxps);
  eqclasses.reserve(numEqClasses);
  string tname;
  std::unordered_map<string, size_t> nameToIndex;
  for (size_t i = 0; i < numTxps; ++i) {
    (*equivFilePtr) >> tname;
    tnames.emplace_back(tname);
    nameToIndex[tname] = i;
  }

  for (size_t i= 0; i < numEqClasses; ++i) {
    size_t classLength;
    (*equivFilePtr) >> classLength;

    // each group member
    uint64_t tid;
    std::vector<uint32_t> tids; tids.reserve(classLength);
    for (size_t i = 0; i < classLength; i++) {
      (*equivFilePtr) >> tid;
      tids.emplace_back(tid);
    }

    double aux;
    std::vector<double> auxs; auxs.reserve(classLength);
    for (size_t i = 0; i < classLength; i++) {
      (*equivFilePtr) >> aux;
      auxs.emplace_back(aux);
    }

    // count for this class
    uint64_t count;
    (*equivFilePtr) >> count;

    eqclasses.emplace_back(tids);
    auxs_vals.emplace_back(auxs);
    eqclass_counts.emplace_back(count);
  }

  tefflens.resize(nameToIndex.size());
  std::vector<size_t> indexList(nameToIndex.size());
  std::iota(indexList.begin(), indexList.end(), 0);
  std::unordered_set<size_t> indexSet(indexList.begin(), indexList.end());

  double tlen;
  while ((*equivFilePtr) >> tname >> tlen) {
    size_t index {0};
    auto it = nameToIndex.find(tname);
    if ( it != nameToIndex.end() ) {
      index = it->second;
    } else {
      l->warn("Missing effective lens for {}", it->first);
      return false;
    }

    indexSet.erase(index);
    tefflens[index] = tlen;
  }

  if (indexSet.size() > 0) {
    l->warn("Missing effective lens for {} transcripts; setting to 100.0.", indexSet.size());
    l->warn("NOTE: Since effective lengths are not provided, please do not rely on the TPM field \n"
            "in the ouput quantifications.  Only the NumReads field will be reliable.");
    for (auto& idx : indexSet) {
      tefflens[idx] = 100.0;
    }
  }

  //equivFile.close();
  return true;
}

/**
 * @param sopt : The salmon options object that tells us the relevant files and contains the pointer 
 *              to the logger object
 * @param transcripts : The list of transcript objects
 * 
 * If the auxTargetFile is not empty (i.e. if the file exists), then read the auxiliary targets
 * and mark them as such in the transcripts vector.  This function will also write a file containing
 * the IDs of the targets marked as auxiliary to the file aux_target_ids.json in the `aux_info` directory.
 **/
void markAuxiliaryTargets(SalmonOpts& sopt, std::vector<Transcript>& transcripts) { 
  
  namespace bfs = boost::filesystem;
  auto& log = sopt.jointLog;
  const std::string& auxTargetFile = sopt.auxTargetFile;

  // If the aux file is empty or doesn't exist
  if (auxTargetFile == "") { 
    return; 
  } else if (!bfs::exists(auxTargetFile)) { 
    log->warn("The auxiliary target file {}, does not exist.  No targets will be treated as auxiliary.", 
               auxTargetFile);
    return; 
  }

  std::ifstream auxFile(auxTargetFile);
  if (!auxFile.good()) { 
    log->warn("Could not open the auxiliary target file {}. No targets will be treated as auxiliary.",
               auxTargetFile);
    return;
  }

  spp::sparse_hash_set<std::string> auxTargetNames;
  std::string tname;
  while (auxFile >> tname) { auxTargetNames.insert(tname); }
  auxFile.close();
  log->info("Parsed {:n} auxiliary targets from {}", auxTargetNames.size(), auxTargetFile);

  size_t tid = 0;
  std::vector<size_t> auxIDs;
  for (auto& txp : transcripts) {
    bool isAux = auxTargetNames.contains(txp.RefName);
    txp.setSkipBiasCorrection(isAux);
    if (isAux) { auxIDs.push_back(tid); }
    ++tid;
  }
  size_t numAuxFound = auxIDs.size();

  if (numAuxFound != auxTargetNames.size()) {
    log->warn("While {:n} auxiliary target names were found in {}, only {:n} were actually found "
              "among tanscripts in the index.  Please make sure that the names in {} match the "
              "transcript names in the index as expected.", auxTargetNames.size(), auxTargetFile, 
              numAuxFound, auxTargetFile);
  }

  // write down the aux target ids in the output directory
  bfs::path auxDir = sopt.outputDirectory / sopt.auxDir;
  if (!bfs::exists(auxDir)) {
    log->warn("The salmon aux directory {} did not exist.  Cannot write aux_target_ids.json!", auxDir);
  }
  bfs::path auxIDFilePath = sopt.outputDirectory / sopt.auxDir / "aux_target_ids.json";
  std::ofstream auxIDFile(auxIDFilePath.string());
  if (auxIDFile.is_open()) {
    nlohmann::json o;
    o["aux_target_ids"] = auxIDs;
    auxIDFile << o;
  } else {
    log->warn("Could not properly open the aux_target_ids file {}.", auxIDFilePath.string());
  }
  auxIDFile.close();
}

/**
 * Computes (and returns) new effective lengths for the transcripts
 * based on the current abundance estimates (alphas) and the current
 * effective lengths (effLensIn).  This approach to sequence-specifc bias is
 * based on the one taken in Roberts et al. (2011) [1].
 * Here, we also consider fragment-GC bias which uses a novel method extending
 * the idea of adjusting the effective lengths.
 *
 * [1] Roberts, Adam, et al. "Improving RNA-Seq expression estimates by
 * correcting for fragment bias."
 *     Genome Biol 12.3 (2011): R22.
 */
template <typename AbundanceVecT, typename ReadExpT>
Eigen::VectorXd
updateEffectiveLengths(SalmonOpts& sopt, ReadExpT& readExp,
                       Eigen::VectorXd& effLensIn, AbundanceVecT& alphas,
                       std::vector<bool>& available, bool writeBias) {

  using std::vector;
  using BlockedIndexRange = tbb::blocked_range<size_t>;
  using salmon::math::EPSILON;
  using salmon::math::LOG_EPSILON;

  // A read cutoff for a txp to be present, adopted from Bray et al. 2016
  double minAlpha = 1e-8;

  double minCDFMass = 1e-10;
  uint32_t gcSamp{sopt.pdfSampFactor};
  bool gcBiasCorrect{sopt.gcBiasCorrect};
  bool seqBiasCorrect{sopt.biasCorrect};
  bool posBiasCorrect{sopt.posBiasCorrect};

  double probFwd = readExp.gcFracFwd();
  double probRC = readExp.gcFracRC();

  if (gcBiasCorrect and probFwd < 0.0) {
    sopt.jointLog->warn("Had no fragments from which to estimate "
                        "fwd vs. rev-comp mapping rate.  Skipping "
                        "sequence-specific / fragment-gc bias correction");
    return effLensIn;
  }

  // calculate read bias normalization factor -- total count in read
  // distribution.
  auto& obs5 = readExp.readBiasModelObserved(salmon::utils::Direction::FORWARD);
  auto& obs3 = readExp.readBiasModelObserved(
      salmon::utils::Direction::REVERSE_COMPLEMENT);
  obs5.normalize();
  obs3.normalize();

  auto& pos5Obs = readExp.posBias(salmon::utils::Direction::FORWARD);
  auto& pos3Obs = readExp.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);

  int32_t K =
      seqBiasCorrect ? static_cast<int32_t>(obs5.getContextLength()) : 1;
  int32_t contextUpstream = seqBiasCorrect ? obs5.contextBefore(false) : 0;

  FragmentLengthDistribution& fld = *(readExp.fragmentLengthDistribution());

  // The *expected* biases from GC effects
  auto& transcriptGCDist = readExp.expectedGCBias();
  auto& gcCounts = readExp.observedGC();
  double readGCNormFactor = 0.0;
  int32_t fldLow{0};
  int32_t fldHigh{1};

  double quantileCutoffLow = 0.005;
  double quantileCutoffHigh = 1.0 - quantileCutoffLow;

  // The CDF and PDF of the fragment length distribution
  std::vector<double> cdf(fld.maxVal() + 1, 0.0);
  std::vector<double> pdf(fld.maxVal() + 1, 0.0);
  {
    transcriptGCDist.reset(distribution_utils::DistributionSpace::LINEAR);

    bool lb{false};
    bool ub{false};
    for (size_t i = 0; i <= fld.maxVal(); ++i) {
      pdf[i] = std::exp(fld.pmf(i));
      cdf[i] = (i > 0) ? cdf[i - 1] + pdf[i] : pdf[i];
      auto density = cdf[i];

      if (!lb and density >= quantileCutoffLow) {
        lb = true;
        fldLow = i;
      }
      if (!ub and density >= quantileCutoffHigh) {
        ub = true;
        fldHigh = i;
      }
    }

    /*
    if (gcBiasCorrect) {
      for (auto& c : gcCounts) {
        readGCNormFactor += c;
      }
    }
    */
  }

  // Make this const so there are no shenanigans
  const auto& transcripts = readExp.transcripts();

  // The effective lengths adjusted for bias
  Eigen::VectorXd effLensOut(effLensIn.size());

  // How much to cut off
  int32_t trunc = K;

  using GCBiasVecT = std::vector<double>;
  using SeqBiasVecT = std::vector<double>;

  /**
   * These will store "thread local" parameters
   * for the appropriate bias terms.
   */
  class CombineableBiasParams {
  public:
    CombineableBiasParams(uint32_t K, size_t numCondBins, size_t numGCBins)
        : expectGC(numCondBins, numGCBins,
                   distribution_utils::DistributionSpace::LINEAR),
          expectPos5(std::vector<SimplePosBias>(5)),
          expectPos3(std::vector<SimplePosBias>(5)) {
      //expectPos5 = std::vector<SimplePosBias>(5);
      //expectPos3 = std::vector<SimplePosBias>(5);
    }

    std::vector<SimplePosBias> expectPos5;
    std::vector<SimplePosBias> expectPos3;
    SBModel expectSeqFW;
    SBModel expectSeqRC;
    GCFragModel expectGC;
  };

  auto revComplement = [](const char* s, int32_t l, std::string& o) -> void {
    if (l > static_cast<int32_t>(o.size())) {
      o.resize(l, 'A');
    }
    int32_t j = 0;
    for (int32_t i = l - 1; i >= 0; --i, ++j) {
      switch (s[i]) {
      case 'A':
      case 'a':
        o[j] = 'T';
        break;
      case 'C':
      case 'c':
        o[j] = 'G';
        break;
      case 'T':
      case 't':
        o[j] = 'A';
        break;
      case 'G':
      case 'g':
        o[j] = 'C';
        break;
      default:
        o[j] = 'N';
        break;
      }
    }
  };

  int outsideContext{3};
  int insideContext{2};

  /**
   * New context counting
   */

  int contextSize = outsideContext + insideContext;
  // double cscale = 100.0 / (2 * contextSize);
  auto populateContextCounts =
      [outsideContext, insideContext, contextSize](
          const Transcript& txp, const char* tseq,
          Eigen::VectorXd& contextCountsFP, Eigen::VectorXd& contextCountsTP,
          Eigen::VectorXd& windowLensFP, Eigen::VectorXd& windowLensTP) {
        auto refLen = static_cast<int32_t>(txp.RefLength);
        auto lastPos = refLen - 1;
        if (refLen > contextSize) {
          // window starts like this
          // -3 === -2 === -1 === 0 === 1
          //         3'           5'
          // and then shifts to the right one base at a time.
          int windowEnd = insideContext - 1;
          int windowStart = -outsideContext;
          int fp = 0;
          int tp = windowStart + (insideContext - 1);
          double count = txp.gcAt(windowEnd - 1);
          for (; tp < refLen; ++fp, ++tp) {
            if (windowStart > 0) {
              switch (tseq[windowStart - 1]) {
              case 'G':
              case 'g':
              case 'C':
              case 'c':
                count -= 1;
              }
            }
            if (windowEnd < refLen) {
              switch (tseq[windowEnd]) {
              case 'G':
              case 'g':
              case 'C':
              case 'c':
                count += 1;
              }
            }
            double actualWindowLength = (windowEnd < contextSize)
                                            ? windowEnd + 1
                                            : (windowEnd - windowStart + 1);
            if (fp < refLen) {
              contextCountsFP[fp] = count;
              windowLensFP[fp] = actualWindowLength;
            }
            if (tp >= 0) {
              contextCountsTP[tp] = count;
              windowLensTP[tp] = actualWindowLength;
            }
            // Shift the end of the window right 1 base
            if (windowEnd < refLen - 1) {
              ++windowEnd;
            }
            ++windowStart;
          }
        }
      };

  /**
   * orig context counting
   **/
  /*
int contextSize = outsideContext + insideContext;
  double cscale = 100.0 / (2 * contextSize);
  auto populateContextCounts = [outsideContext, insideContext, contextSize](
      const Transcript& txp, const char* tseq, Eigen::VectorXd& contextCountsFP,
      Eigen::VectorXd& contextCountsTP) {
    auto refLen = static_cast<int32_t>(txp.RefLength);
    auto lastPos = refLen - 1;
    if (refLen > contextSize) {
      int windowStart = -1;
      int windowEnd = contextSize - 1;
      int fp = outsideContext;
      int tp = insideContext - 1;
      double count = txp.gcAt(windowEnd);
      contextCountsFP[fp] = count;
      contextCountsTP[tp] = count;
      ++windowStart;
      ++windowEnd;
      ++fp;
      ++tp;
      for (; tp < refLen; ++windowStart, ++windowEnd, ++fp, ++tp) {
        switch (tseq[windowStart]) {
        case 'G':
        case 'g':
        case 'C':
        case 'c':
          count -= 1;
        }
        if (windowEnd < refLen) {
          switch (tseq[windowEnd]) {
          case 'G':
          case 'g':
          case 'C':
          case 'c':
            count += 1;
          }
        }
        if (fp < refLen) {
          contextCountsFP[fp] = count;
        }
        contextCountsTP[tp] = count;
      }
    }
  };
  */

  /**
   * The local bias terms from each thread can be combined
   * via simple summation.
   */
  auto getBiasParams = [K, &sopt]() -> CombineableBiasParams {
    return CombineableBiasParams(K, sopt.numConditionalGCBins,
                                 sopt.numFragGCBins);
  };
  tbb::combinable<CombineableBiasParams> expectedDist(getBiasParams);
  std::atomic<size_t> numBackgroundTranscripts{0};

  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(transcripts.size())),
      [&](const BlockedIndexRange& range) -> void {

        auto& expectSeqFW = expectedDist.local().expectSeqFW;
        auto& expectSeqRC = expectedDist.local().expectSeqRC;
        auto& expectGC = expectedDist.local().expectGC;
        auto& expectPos5 = expectedDist.local().expectPos5;
        auto& expectPos3 = expectedDist.local().expectPos3;

        std::string rcSeq;
        // For each transcript
        for (auto it : boost::irange(range.begin(), range.end())) {

          // Get the transcript
          const auto& txp = transcripts[it];

          // If this transcript is in the list of transcripts for which the user 
          // has requested we skip bias correction, then do so
          if (txp.skipBiasCorrection()) { continue; }

          // Get the reference length and the
          // "initial" effective length (not considering any biases)
          int32_t refLen = static_cast<int32_t>(txp.RefLength);
          int32_t elen = static_cast<int32_t>(txp.EffectiveLength);

          // The difference between the actual and effective length
          int32_t unprocessedLen = std::max(0, refLen - elen);

          int32_t cdfMaxArg =
              std::min(static_cast<int32_t>(cdf.size() - 1), refLen);
          double cdfMaxVal = cdf[cdfMaxArg];
          // need a reliable CDF
          if (cdfMaxVal < minCDFMass) {
            continue;
          }
          auto conditionalCDF = [cdfMaxArg, cdfMaxVal,
                                 &cdf](double x) -> double {
            return (x > cdfMaxArg) ? 1.0 : (cdf[x] / cdfMaxVal);
          };

          // Skip transcripts with trivial expression or that are too
          // short
          if (alphas[it] < minAlpha or unprocessedLen <= 0) { // or txp.uniqueUpdateFraction() < 0.90) {
            continue;
          }
          ++numBackgroundTranscripts;

          // Otherwise, proceed giving this transcript the following weight
          double weight = (alphas[it] / effLensIn(it));

          Eigen::VectorXd contextCountsFP(refLen);
          Eigen::VectorXd contextCountsTP(refLen);
          Eigen::VectorXd windowLensFP(refLen);
          Eigen::VectorXd windowLensTP(refLen);
          contextCountsFP.setZero();
          contextCountsTP.setZero();
          windowLensFP.setZero();
          windowLensTP.setZero();

          // This transcript's sequence
          bool have_seq = txp.have_sequence();
          SBMer fwmer;
          SBMer rcmer;

          const char* tseq = have_seq ? txp.Sequence() : nullptr;
          // only do this if we have the sequence, we may not if just pos. bias.
          if (have_seq) { 
            revComplement(tseq, refLen, rcSeq);
            fwmer.fromChars(tseq);
            rcmer.fromChars(rcSeq);
          }
          // may be empty, but we shouldn't actually
          // use it unless it is meaningful.
          const char* rseq = rcSeq.c_str();

          int32_t contextLength{expectSeqFW.getContextLength()};

          if (gcBiasCorrect and seqBiasCorrect) {
            populateContextCounts(txp, tseq, contextCountsFP, contextCountsTP,
                                  windowLensFP, windowLensTP);
          }

          // The smallest and largest values of fragment
          // lengths we'll consider for this transcript.
          int32_t locFLDLow = (refLen < cdfMaxArg) ? 1 : fldLow;
          int32_t locFLDHigh = (refLen < cdfMaxArg) ? cdfMaxArg : fldHigh;

          // For each position along the transcript
          // Starting from the 5' end and moving toward the 3' end
          for (int32_t fragStartPos = 0; fragStartPos < refLen - K;
               ++fragStartPos) {
            // Seq-specific bias
            if (seqBiasCorrect) {
              int32_t contextEndPos =
                  fragStartPos + K - 1; // -1 because pos is *inclusive*

              if (contextEndPos >= 0 and contextEndPos < refLen) {
                int32_t maxFragLen =
                    refLen - (fragStartPos + expectSeqFW.contextBefore(false));
                if (maxFragLen >= 0 and maxFragLen < refLen) {
                  auto cdensity = conditionalCDF(maxFragLen);
                  expectSeqFW.addSequence(fwmer, weight * cdensity);
                  expectSeqRC.addSequence(rcmer, weight * cdensity);
                }
              }

              // shift the context one nucleotide to the right
              //fwmer.shift_left(tseq[fragStartPos + contextLength]);
              //rcmer.shift_left(rseq[fragStartPos + contextLength]);
              fwmer.append(tseq[fragStartPos + contextLength]);
              rcmer.append(rseq[fragStartPos + contextLength]);
            } // end: Seq-specific bias

            // fragment-GC bias
            if (gcBiasCorrect) {
              size_t sp =
                  static_cast<size_t>((locFLDLow > 0) ? locFLDLow - 1 : 0);
              double prevFLMass = conditionalCDF(sp);
              int32_t fragStart = fragStartPos;
              for (int32_t fl = locFLDLow; fl <= locFLDHigh; fl += gcSamp) {
                int32_t fragEnd = fragStart + fl - 1;
                if (fragEnd < refLen) {
                  // The GC fraction for this putative fragment
                  auto gcFrac = txp.gcFrac(fragStart, fragEnd);
                  /*
                  int32_t contextFrac = std::lrint(
                                                   (contextCountsFP[fragStart] +
                  contextCountsTP[fragEnd]) * cscale);
                  */
                  double contextLength =
                      (windowLensFP[fragStart] + windowLensTP[fragEnd]);
                  int32_t contextFrac =
                      (contextLength > 0)
                          ? (std::lrint(100.0 *
                                        (contextCountsFP[fragStart] +
                                         contextCountsTP[fragEnd]) /
                                        contextLength))
                          : 0;

                  GCDesc desc{gcFrac, contextFrac};
                  expectGC.inc(desc,
                               weight * (conditionalCDF(fl) - prevFLMass));
                  prevFLMass = conditionalCDF(fl);
                } else {
                  break;
                } // no more valid positions
              }   // end: for each fragment length
            }     // end: fragment GC bias

            // positional bias
            if (posBiasCorrect) {
              int32_t maxFragLenFW = refLen - fragStartPos + 1;
              int32_t maxFragLenRC = fragStartPos;
              auto densityFW = conditionalCDF(maxFragLenFW);
              auto densityRC = conditionalCDF(maxFragLenRC);
              if (weight * densityFW > EPSILON) {
                expectPos5[txp.lengthClassIndex()].addMass(
                    fragStartPos, txp.RefLength, std::log(weight * densityFW));
              }
              if (weight * densityRC > EPSILON) {
                expectPos3[txp.lengthClassIndex()].addMass(
                    fragStartPos, txp.RefLength, std::log(weight * densityRC));
              }
            }
          } // end: for every fragment start position
        }   // end for each transcript

      } // end tbb for function
  );

  size_t bgCutoff =
      std::min(static_cast<size_t>(150),
               static_cast<size_t>(numBackgroundTranscripts * 0.1));
  if (numBackgroundTranscripts < bgCutoff) {
    sopt.jointLog->warn("I found only {} transcripts meeting the necessary "
                        "conditions to contribute to "
                        "the bias background distribution.  This is likely too "
                        "small to safely do bias correction. "
                        "I'm skipping bias correction",
                        numBackgroundTranscripts.load());
    sopt.biasCorrect = false;
    sopt.gcBiasCorrect = false;
    sopt.posBiasCorrect = false;
    return effLensIn;
  }

  /**
   * The local bias terms from each thread can be combined
   * via simple summation.  Here, we combine the locally-computed
   * bias terms.
   */
  SBModel& exp5 = readExp.readBiasModelExpected(salmon::utils::Direction::FORWARD);
  SBModel& exp3 = readExp.readBiasModelExpected(salmon::utils::Direction::REVERSE_COMPLEMENT);

  auto& pos5Exp = readExp.posBiasExpected(salmon::utils::Direction::FORWARD);
  auto& pos3Exp =
      readExp.posBiasExpected(salmon::utils::Direction::REVERSE_COMPLEMENT);

  auto combineBiasParams =
      [seqBiasCorrect, gcBiasCorrect, posBiasCorrect, &pos5Exp, &pos3Exp, &exp5,
       &exp3, &transcriptGCDist](const CombineableBiasParams& p) -> void {
    if (seqBiasCorrect) {
      exp5.combineCounts(p.expectSeqFW);
      exp3.combineCounts(p.expectSeqRC);
    }
    if (gcBiasCorrect) {
      transcriptGCDist.combineCounts(p.expectGC);
    }
    if (posBiasCorrect) {
      for (size_t i = 0; i < p.expectPos5.size(); ++i) {
        pos5Exp[i].combine(p.expectPos5[i]);
        pos3Exp[i].combine(p.expectPos3[i]);
      }
    }
  };
  expectedDist.combine_each(combineBiasParams);

  // finalize expected positional biases
  if (posBiasCorrect) {
    for (size_t i = 0; i < pos5Exp.size(); ++i) {
      pos5Exp[i].finalize();
      pos3Exp[i].finalize();
    }
  }
  if (gcBiasCorrect) {
    transcriptGCDist.normalize();
  }

  sopt.jointLog->info("Computed expected counts (for bias correction)");

  auto gcBias = gcCounts.ratio(transcriptGCDist, 1000.0);

  exp5.normalize();
  exp3.normalize();

  bool noThreshold = sopt.noBiasLengthThreshold;
  std::atomic<size_t> numCorrected{0};
  std::atomic<size_t> numUncorrected{0};

  std::atomic<uint32_t> numProcessed{0};
  size_t numTranscripts = transcripts.size();
  size_t stepSize = static_cast<size_t>(transcripts.size() * 0.1);
  size_t nextUpdate{0};

  // std::mutex updateMutex;
  TryableSpinLock tsl;
  /**
   * Compute the effective lengths of each transcript (in parallel)
   */
  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(transcripts.size())),
      [&](const BlockedIndexRange& range) -> void {

        std::string rcSeq;
        // For each transcript
        for (auto it : boost::irange(range.begin(), range.end())) {

          auto& txp = transcripts[it];

          // eff. length starts out as 0
          double effLength = 0.0;

          // Reference length
          int32_t refLen = static_cast<int32_t>(txp.RefLength);
          // Effective length before any bias correction
          int32_t elen = static_cast<int32_t>(txp.EffectiveLength);

          // How much of this transcript (beginning and end) should
          // not be considered
          int32_t unprocessedLen = std::max(0, refLen - elen);
          int32_t cdfMaxArg =
              std::min(static_cast<int32_t>(cdf.size() - 1), refLen);
          double cdfMaxVal = cdf[cdfMaxArg];
          auto conditionalCDF = [cdfMaxArg, cdfMaxVal,
                                 &cdf](double x) -> double {
            return (x > cdfMaxArg) ? 1.0 : (cdf[x] / cdfMaxVal);
          };
          // The smallest and largest values of fragment
          // lengths we'll consider for this transcript.
          int32_t locFLDLow = (refLen < cdfMaxArg) ? 1 : fldLow;
          int32_t locFLDHigh = (refLen < cdfMaxArg) ? cdfMaxArg : fldHigh;
          bool wasProcessed{false};
          if (alphas[it] >= minAlpha
              // available[it]
              and unprocessedLen > 0 and cdfMaxVal > minCDFMass) {

            Eigen::VectorXd seqFactorsFW(refLen);
            Eigen::VectorXd seqFactorsRC(refLen);
            seqFactorsFW.setOnes();
            seqFactorsRC.setOnes();

            Eigen::VectorXd contextCountsFP(refLen);
            Eigen::VectorXd contextCountsTP(refLen);
            Eigen::VectorXd windowLensFP(refLen);
            Eigen::VectorXd windowLensTP(refLen);
            contextCountsFP.setZero();
            contextCountsTP.setZero();
            windowLensFP.setZero();
            windowLensTP.setZero();

            std::vector<double> posFactorsFW(refLen, 1.0);
            std::vector<double> posFactorsRC(refLen, 1.0);

            // This transcript's sequence
            bool have_seq = txp.have_sequence();
            const char* tseq = have_seq ? txp.Sequence() : nullptr;
            // only do this if we have the sequence, we may not if just pos.
            // bias.
            if (have_seq) {
              revComplement(tseq, refLen, rcSeq);
            }
            // may be empty, but we shouldn't actually
            // use it unless it is meaningful.
            const char* rseq = rcSeq.c_str();

            int32_t fl = locFLDLow;
            auto maxLen = std::min(refLen, locFLDHigh + 1);
            bool done{fl >= maxLen};

            if (gcBiasCorrect and seqBiasCorrect) {
              populateContextCounts(txp, tseq, contextCountsFP, contextCountsTP,
                                    windowLensFP, windowLensTP);
            }

            if (posBiasCorrect) {
              std::vector<double> posFactorsObs5(refLen, 1.0);
              std::vector<double> posFactorsObs3(refLen, 1.0);
              std::vector<double> posFactorsExp5(refLen, 1.0);
              std::vector<double> posFactorsExp3(refLen, 1.0);
              auto li = txp.lengthClassIndex();
              auto& p5O = pos5Obs[li];
              auto& p3O = pos3Obs[li];
              auto& p5E = pos5Exp[li];
              auto& p3E = pos3Exp[li];
              p5O.projectWeights(posFactorsObs5);
              p3O.projectWeights(posFactorsObs3);
              p5E.projectWeights(posFactorsExp5);
              p3E.projectWeights(posFactorsExp3);
              for (int32_t fragStart = 0; fragStart < refLen - K; ++fragStart) {
                posFactorsFW[fragStart] =
                    posFactorsObs5[fragStart] / posFactorsExp5[fragStart];
                posFactorsRC[fragStart] =
                    posFactorsObs3[fragStart] / posFactorsExp3[fragStart];
              }
            }

            // Evaluate the sequence specific bias (5' and 3') over the length
            // of the transcript.  After this loop,
            // seqFactorsFW will contain the sequence-specific bias for each
            // position on the 5' strand
            // and seqFactorsRC will contain the sequence-specific bias for each
            // position on the 3' strand.
            if (seqBiasCorrect) {
              //Mer mer;
              //Mer rcmer;
              //mer.from_chars(tseq);
              //rcmer.from_chars(rseq);
              SBMer mer;
              SBMer rcmer;
              mer.fromChars(tseq);
              rcmer.fromChars(rseq);
 
              int32_t contextLength{exp5.getContextLength()};

              for (int32_t fragStart = 0; fragStart < refLen - K; ++fragStart) {
                int32_t readStart = fragStart + obs5.contextBefore(false);
                int32_t kmerEndPos =
                    fragStart + K - 1; // -1 because pos is *inclusive*

                if (kmerEndPos >= 0 and kmerEndPos < refLen and
                    readStart < refLen) {
                  seqFactorsFW[readStart] =
                      std::exp(obs5.evaluateLog(mer) - exp5.evaluateLog(mer));
                  seqFactorsRC[readStart] = std::exp(obs3.evaluateLog(rcmer) -
                                                     exp3.evaluateLog(rcmer));
                }
                // shift the context one nucleotide to the right
                //mer.shift_left(tseq[fragStart + contextLength]);
                //rcmer.shift_left(rseq[fragStart + contextLength]);
                mer.append(tseq[fragStart + contextLength]);
                rcmer.append(rseq[fragStart + contextLength]);
              }
              // We need these in 5' -> 3' order, so reverse them
              seqFactorsRC.reverseInPlace();
            } // end sequence-specific factor calculation

            if (numProcessed > nextUpdate) {
              if (tsl.try_lock()) {
                if (numProcessed > nextUpdate) {
                  sopt.jointLog->info(
                      "processed bias for {:3.1f}% of the transcripts",
                      100.0 *
                          (numProcessed / static_cast<double>(numTranscripts)));
                  nextUpdate += stepSize;
                  if (nextUpdate > numTranscripts) {
                    nextUpdate = numTranscripts - 1;
                  }
                }
                tsl.unlock();
              }
            }

            size_t sp = static_cast<size_t>((fl > 0) ? fl - 1 : 0);
            double prevFLMass = conditionalCDF(sp);
            double unbiasedMass{0.0};

            // For every possible fragment length
            while (!done) {
              if (fl >= maxLen) {
                done = true;
                fl = maxLen - 1;
              }
              double flWeight = conditionalCDF(fl) - prevFLMass;
              prevFLMass = conditionalCDF(fl);

              double flMassTotal{0.0};
              // For every position a fragment of length fl could start
              for (int32_t kmerStartPos = 0; kmerStartPos < refLen - fl;
                   ++kmerStartPos) {
                int32_t fragStart = kmerStartPos;
                int32_t fragEnd = fragStart + fl - 1;

                // If the 3' end is within the transcript
                if (fragStart < refLen and fragEnd < refLen) {
                  double fragFactor =
                      seqFactorsFW[fragStart] * seqFactorsRC[fragEnd];
                  if (gcBiasCorrect) {
                    auto gcFrac = txp.gcFrac(fragStart, fragEnd);
                    /*
                    int32_t contextFrac =
                      std::lrint((contextCountsFP[fragStart] +
                                  contextCountsTP[fragEnd]) *
                                 cscale);
                    */
                    double contextLength =
                        (windowLensFP[fragStart] + windowLensTP[fragEnd]);
                    int32_t contextFrac =
                        (contextLength > 0)
                            ? (std::lrint(100.0 *
                                          (contextCountsFP[fragStart] +
                                           contextCountsTP[fragEnd]) /
                                          contextLength))
                            : 0;

                    GCDesc desc{gcFrac, contextFrac};
                    fragFactor *= gcBias.get(desc);
                    /*
                    fragFactor *= gcBias[gcFrac];
                    */
                  }
                  if (posBiasCorrect) {
                    fragFactor *=
                        posFactorsFW[fragStart] * posFactorsRC[fragEnd];
                  }
                  flMassTotal += fragFactor;
                } else {
                  break;
                }
              }

              effLength += (flWeight * flMassTotal);
              fl += gcSamp;
            }
            wasProcessed = true;
          } // for the processed transcript

          if (wasProcessed) {
            // throw caution to the wind
            double thresh = noThreshold ? 1.0 : unprocessedLen;
            if (noThreshold) {
              if (unprocessedLen > 0.0 and effLength > thresh) {
                effLensOut(it) = effLength;
              } else {
                effLensOut(it) = effLensIn(it);
              }
            } else {
              double offset = std::max(1.0, thresh);
              double effLengthNoBias = static_cast<double>(elen);
              auto barrierLength = [effLengthNoBias, offset](double x) -> double {
                                     return std::max(x, std::min(effLengthNoBias, offset));
                                   };
              effLensOut(it) = barrierLength(effLength);
            }
          } else {
            effLensOut(it) = static_cast<double>(elen);
          }
          ++numProcessed;
        }
      } // end parallel_for lambda
  );

  sopt.jointLog->info("processed bias for 100.0% of the transcripts");
  return effLensOut;
}

void aggregateEstimatesToGeneLevel(TranscriptGeneMap& tgm,
                                   boost::filesystem::path& inputPath) {
  using std::vector;
  using std::string;
  using std::ofstream;
  using std::unordered_map;
  using std::move;
  using std::cerr;
  using std::max;

  auto logger = spdlog::get("jointLog");

  logger->info("NOTE: We recommend using tximport (https://bioconductor.org/packages/release/bioc/html/tximport.html) "
               "for aggregating transcript-level salmon abundance estimates to the gene level.  It is more versatile, "
               "exposes more features, and allows considering multi-sample information during aggregation.");

  constexpr double minTPM = std::numeric_limits<double>::denorm_min();
  std::ifstream expFile(inputPath.string());

  if (!expFile.is_open()) {
    perror("Error reading file");
  }

  //====================== From GeneSum ====================
  vector<string> comments;
  unordered_map<string, vector<ExpressionRecord>> geneExps;
  string l;
  size_t ln{0};

  bool headerLine{true};
  while (getline(expFile, l)) {
    auto it =
        find_if(l.begin(), l.end(), [](char c) -> bool { return !isspace(c); });
    if (it != l.end()) {
      if (*it == '#') {
        comments.push_back(l);
      } else {
        // If this isn't the first non-comment line
        if (!headerLine) {
          vector<string> toks = split(l);
          ExpressionRecord er(toks);
          bool foundGene{false};
          auto gn = tgm.geneName(er.target, foundGene);
          if (!foundGene) {
            logger->warn("couldn't find transcript named [{}] in transcript "
                         "<-> gene map; "
                         "returning transcript as it's own gene",
                         er.target);
          }
          geneExps[gn].push_back(move(er));
        } else { // treat the header line as a comment
          comments.push_back(l);
          headerLine = false;
        }
      }
    }
  }
  expFile.close();

  logger->info("Aggregating expressions to gene level");
  boost::filesystem::path outputFilePath(inputPath);
  outputFilePath.replace_extension(".genes.sf");
  ofstream outFile(outputFilePath.string());

  // preserve any comments in the output
  for (auto& c : comments) {
    outFile << c << '\n';
  }

  for (auto& kv : geneExps) {
    auto& gn = kv.first;

    double geneLength = kv.second.front().length;
    double geneEffLength = kv.second.front().effLength;
    vector<double> expVals(kv.second.front().expVals.size(), 0);
    const size_t NE{expVals.size()};

    size_t tpmIdx{0};
    for (auto& tranExp : kv.second) {
      // expVals[0] = TPM
      // expVals[1] = count
      for (size_t i = 0; i < NE; ++i) {
        expVals[i] += tranExp.expVals[i];
      }
    }
    double totalTPM = expVals[tpmIdx];

    // If this gene was expressed
    if (totalTPM > minTPM) {
      geneLength = 0.0;
      geneEffLength = 0.0;
      for (auto& tranExp : kv.second) {
        double frac = tranExp.expVals[tpmIdx] / totalTPM;
        geneLength += tranExp.length * frac;
        geneEffLength += tranExp.effLength * frac;
      }
    } else {
      geneLength = 0.0;
      geneEffLength = 0.0;
      double frac = 1.0 / kv.second.size();
      for (auto& tranExp : kv.second) {
        geneLength += tranExp.length * frac;
        geneEffLength += tranExp.effLength * frac;
      }
    }

    // Otherwise, if the gene wasn't expressed, the length
    // is reported as the longest transcript length.

    outFile << gn << '\t' << geneLength << '\t' << geneEffLength;
    for (size_t i = 0; i < NE; ++i) {
      outFile << '\t' << expVals[i];
    }
    outFile << '\n';
  }

  outFile.close();
  logger->info("done");
  //====================== From GeneSum =====================
}

void generateGeneLevelEstimates(boost::filesystem::path& geneMapPath,
                                boost::filesystem::path& estDir) {
  namespace bfs = boost::filesystem;
  auto logger = spdlog::get("jointLog");
  logger->info("Computing gene-level abundance estimates");
  std::set<std::string> validGTFExtensions = {".gtf", ".gff", ".gff3",
                                              ".GTF", ".GFF", ".GFF3"};
  auto extension = geneMapPath.extension();

  TranscriptGeneMap tranGeneMap;
  // parse the map as a GTF file
  if (validGTFExtensions.find(extension.string()) != validGTFExtensions.end()) {
    // Using libgff
    tranGeneMap = salmon::utils::transcriptGeneMapFromGTF(geneMapPath.string(),
                                                          "gene_id");
  } else { // parse the map as a simple format files
    std::ifstream tgfile(geneMapPath.string());
    tranGeneMap = salmon::utils::readTranscriptToGeneMap(tgfile);
    tgfile.close();
  }

  logger->info("There were {} transcripts mapping to {} genes",
               tranGeneMap.numTranscripts(), tranGeneMap.numGenes());

  bfs::path estFilePath = estDir / "quant.sf";
  if (!bfs::exists(estFilePath)) {
    std::stringstream errstr;
    errstr << "Attempting to compute gene-level esimtates, but could not \n"
           << "find isoform-level file " << estFilePath;
    throw std::invalid_argument(errstr.str());
  } else {
    salmon::utils::aggregateEstimatesToGeneLevel(tranGeneMap, estFilePath);
  }

  /** Create a gene-level summary of the bias-corrected estimates as well if
   * these exist **/
  /*
  if (haveBiasCorrectedFile) {
      bfs::path biasCorrectEstFilePath = estDir / "quant_bias_corrected.sf";
      if (!bfs::exists(biasCorrectEstFilePath)) {
          std::stringstream errstr;
          errstr << "Attempting to compute gene-level esimtates, but could not
  \n"
              << "find bias-corrected isoform-level file " <<
  biasCorrectEstFilePath;
          throw std::invalid_argument(errstr.str());
      } else {
          salmon::utils::aggregateEstimatesToGeneLevel(tranGeneMap,
  biasCorrectEstFilePath);
      }
  }
  */
}

double compute_1_edit_threshold(int32_t l, const SalmonOpts& sopt) {
  int32_t match = sopt.matchScore;
  int32_t mismatch = (sopt.mismatchPenalty < 0) ? sopt.mismatchPenalty : -sopt.mismatchPenalty;
  int32_t go = (sopt.gapOpenPenalty < 0) ? sopt.gapOpenPenalty : -sopt.gapOpenPenalty;
  int32_t ge = (sopt.gapExtendPenalty < 0) ? sopt.gapExtendPenalty : -sopt.gapExtendPenalty;

  // cost of subst = mismatch - cost of losing a match
  // cost of deletion = gap open + gap extend - cost of losing a match
  int32_t edit_cost = std::min(mismatch - match, go + ge - match);
  int32_t max_score = l * match;
  return  (static_cast<double>(max_score + edit_cost) - 0.5) / max_score;
}

  /**
   * The fix below is from : https://github.com/coryan/google-cloud-cpp-common/blob/a6e7b6b362d72451d6dc1fec5bc7643693dbea96/google/cloud/internal/random.cc
   **/
  /*
// This function returns a random device.  Ideally, we would simply use 
// the default constructor.  However, apparently certain versions of 
// libstdc++ break this when you need mulitple random divices across 
// threads.  This function should work in a robust way across 
// multiple libstdc++ implementations.
std::random_device get_random_device() {

    // We use the default C++ random device to generate entropy.  The quality of
  // this source of entropy is implementation-defined:
  //     http://en.cppreference.com/w/cpp/numeric/random/random_device/random_device
  // However, in all the platforms we care about, this seems to be a reasonably
  // non-deterministic source of entropy.
  //
  // On Linux with libstdc++ (the most common library on Linux) is based on
  // either `/dev/urandom`, or (if available) the RDRAND [1], or the RDSEED [1]
  // CPU instructions.
  //
  // On Linux with libc++ the default seems to be using `/dev/urandom`, but the
  // library could have been compiled [2] to use `getentropy(3)` [3],
  // `arc4random()` [4], or even `nacl` [5]. It does not seem like the library
  // uses the RDRAND or RDSEED instructions directly. In any case, it seems that
  // all the choices are good entropy sources.
  //
  // With MSVC the documentation says that the numbers are not deterministic,
  // and cryptographically secure, but no further details are available:
  //     https://docs.microsoft.com/en-us/cpp/standard-library/random-device-class
  //
  // On macOS the library is libc++ implementation so the previous comments
  // apply.
  //
  // One would want to simply pass a `std::random_device` to the constructor for
  // the random bit generators, but the C++11 approach is annoying, see this
  // critique for the details:
  //     http://www.pcg-random.org/posts/simple-portable-cpp-seed-entropy.html
  //
  // [1]: https://en.wikipedia.org/wiki/RDRAND
  // [2]: https://github.com/llvm-mirror/libcxx/blob/master/src/random.cpp
  // [3]: http://man7.org/linux/man-pages/man3/getentropy.3.html
  // [4]: https://linux.die.net/man/3/arc4random
  // [5]: https://en.wikipedia.org/wiki/NaCl_(software)
  //
#if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
  // Workaround for a libstd++ bug:
  //     https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94087
  // I (@coryan) would like to use `rdrand` as the token, but:
  //   - `rdrand` is not supported by all versions of libstdc++
  //   - even for the versions that do support it, that is CPU specific
  //   - I know of no reliable way to detect if the library version supports
  //     `rdrand` (other than trying and getting an exception), for
  //     example:
  //     * __GLIBCXX__ is 2020318 on Fedora:31 with g++-9.3.1, but it does
  //       *not* support `rdrand`
  //     * __GLIBCXX__ is 20200306 on openSUSE/Tumbleweed, with g++-9.2.1,
  //       but it *does* support `rdrand`
  return std::random_device rd("/dev/urandom");
#else
  return std::random_device rd;
#endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
}
*/


} // namespace utils
} // namespace salmon

// === Explicit instantiations
using SCExpT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;
using BulkExpT = ReadExperiment<EquivalenceClassBuilder<TGValue>>;
template <typename FragT>
using BulkAlignLibT = AlignmentLibrary<FragT, EquivalenceClassBuilder<TGValue>>;

// explicit instantiations for writing abundances ---
template void salmon::utils::writeAbundances<BulkAlignLibT<ReadPair>>(
    const SalmonOpts& opts, BulkAlignLibT<ReadPair>& alnLib,
    boost::filesystem::path& fname, std::string headerComments);

template void salmon::utils::writeAbundances<BulkAlignLibT<UnpairedRead>>(
    const SalmonOpts& opts, BulkAlignLibT<UnpairedRead>& alnLib,
    boost::filesystem::path& fname, std::string headerComments);

template void salmon::utils::writeAbundances<BulkExpT>(
    const SalmonOpts& opts, BulkExpT& alnLib,
    boost::filesystem::path& fname, std::string headerComments);

template void salmon::utils::writeAbundances<SCExpT>(
                                                             const SalmonOpts& opts, SCExpT& alnLib,
                                                             boost::filesystem::path& fname, std::string headerComments);


template void
salmon::utils::writeAbundancesFromCollapsed<BulkAlignLibT<ReadPair>>(
    const SalmonOpts& opts, BulkAlignLibT<ReadPair>& alnLib,
    boost::filesystem::path& fname, std::string headerComments);

template void
salmon::utils::writeAbundancesFromCollapsed<BulkAlignLibT<UnpairedRead>>(
    const SalmonOpts& opts, BulkAlignLibT<UnpairedRead>& alnLib,
    boost::filesystem::path& fname, std::string headerComments);

template void salmon::utils::writeAbundancesFromCollapsed<BulkExpT>(
    const SalmonOpts& opts, BulkExpT& alnLib,
    boost::filesystem::path& fname, std::string headerComments);

template void salmon::utils::writeAbundancesFromCollapsed<SCExpT>(
                                                                          const SalmonOpts& opts, SCExpT& alnLib,
                                                                          boost::filesystem::path& fname, std::string headerComments);
// explicit instantiations for normalizing alpha vectors ---
template void
salmon::utils::normalizeAlphas<BulkExpT>(const SalmonOpts& sopt,
                                               BulkExpT& alnLib);
template void
salmon::utils::normalizeAlphas<SCExpT>(const SalmonOpts& sopt,
                                               SCExpT& alnLib);

template void salmon::utils::normalizeAlphas<BulkAlignLibT<UnpairedRead>>(
    const SalmonOpts& sopt, BulkAlignLibT<UnpairedRead>& alnLib);
template void salmon::utils::normalizeAlphas<BulkAlignLibT<ReadPair>>(
    const SalmonOpts& sopt, BulkAlignLibT<ReadPair>& alnLib);

// explicit instantiations for effective length updates ---
template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<std::atomic<double>>,
                                      BulkExpT>(
    SalmonOpts& sopt, BulkExpT& readExp, Eigen::VectorXd& effLensIn,
    std::vector<std::atomic<double>>& alphas, std::vector<bool>& available,
    bool finalRound);

template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<std::atomic<double>>,
                                      SCExpT>(
                                                      SalmonOpts& sopt, SCExpT& readExp, Eigen::VectorXd& effLensIn,
                                                      std::vector<std::atomic<double>>& alphas, std::vector<bool>& available,
                                                      bool finalRound);

template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<double>, BulkExpT>(
    SalmonOpts& sopt, BulkExpT& readExp, Eigen::VectorXd& effLensIn,
    std::vector<double>& alphas, std::vector<bool>& available, bool finalRound);

template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<double>, SCExpT>(
                                                                           SalmonOpts& sopt, SCExpT& readExp, Eigen::VectorXd& effLensIn,
                                                                           std::vector<double>& alphas, std::vector<bool>& available, bool finalRound);

template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<std::atomic<double>>,
                                      BulkAlignLibT<ReadPair>>(
    SalmonOpts& sopt, BulkAlignLibT<ReadPair>& readExp,
    Eigen::VectorXd& effLensIn, std::vector<std::atomic<double>>& alphas,
    std::vector<bool>& available, bool finalRound);

template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<double>,
                                      BulkAlignLibT<ReadPair>>(
    SalmonOpts& sopt, BulkAlignLibT<ReadPair>& readExp,
    Eigen::VectorXd& effLensIn, std::vector<double>& alphas,
    std::vector<bool>& available, bool finalRound);

template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<std::atomic<double>>,
                                      BulkAlignLibT<UnpairedRead>>(
    SalmonOpts& sopt, BulkAlignLibT<UnpairedRead>& readExp,
    Eigen::VectorXd& effLensIn, std::vector<std::atomic<double>>& alphas,
    std::vector<bool>& available, bool finalRound);

template Eigen::VectorXd
salmon::utils::updateEffectiveLengths<std::vector<double>,
                                      BulkAlignLibT<UnpairedRead>>(
    SalmonOpts& sopt, BulkAlignLibT<UnpairedRead>& readExp,
    Eigen::VectorXd& effLensIn, std::vector<double>& alphas,
    std::vector<bool>& available, bool finalRound);

//// 0th order model --- code for computing bias factors.

/*
double scoreExpected5 = 1.0;
double scoreExpected3 = 1.0;
double scoreObserved5 = 1.0;
double scoreObserved3 = 1.0;
for (size_t i = 0; i < 6; ++i) {
                    scoreExpected5 *= zeroOrderModel5e(idx(tseq[fragStart + i]),
i);
                    scoreObserved5 *= zeroOrderModel5o(idx(tseq[fragStart + i]),
i);
                    scoreExpected3 *= zeroOrderModel3e(idx(rseq[fragStart + i]),
i);
                    scoreObserved3 *= zeroOrderModel3o(idx(rseq[fragStart + i]),
i);
                }
                seqFactorsFW[fragStart + 2] = scoreObserved5 / scoreExpected5;
                seqFactorsRC[fragStart + 3] = scoreObserved3 / scoreExpected3;
                */

// seqFactorsRC[fragStart + 1] = seqBias3p[dis(gen)];
/*
              seqFactorsFW[fragStart] +=
                  (readBiasFW.counts[idxFW]/ (transcriptKmerDistFW[idxFW] +
   seqPriorFW));
              seqFactorsRC[fragStart] +=
                  (readBiasRC.counts[idxRC]/ (transcriptKmerDistRC[idxRC] +
   seqPriorRC));
*/

/**
 * Computes (and returns) new effective lengths for the transcripts
 * based on the current abundance estimates (alphas) and the current
 * effective lengths (effLensIn).  This approach to sequence-specifc bias is
 * based on the one taken in Roberts et al. (2011) [1].
 * Here, we also consider fragment-GC bias which uses a novel method extending
 * the idea of adjusting the effective lengths.
 *
 * [1] Roberts, Adam, et al. "Improving RNA-Seq expression estimates by
 * correcting for fragment bias."
 *     Genome Biol 12.3 (2011): R22.
 */
/*
template <typename AbundanceVecT, typename ReadExpT>
Eigen::VectorXd updateEffectiveLengths(SalmonOpts& sopt, ReadExpT& readExp,
                                       Eigen::VectorXd& effLensIn,
                                       AbundanceVecT& alphas, bool writeBias) {

  using std::vector;
  using BlockedIndexRange = tbb::blocked_range<size_t>;

  double minAlpha = 1e-8;
  uint32_t gcSamp{sopt.pdfSampFactor};
  bool gcBiasCorrect{sopt.gcBiasCorrect};
  bool seqBiasCorrect{sopt.biasCorrect};
  bool posBiasCorrect{sopt.posBiasCorrect};

  double probFwd = readExp.gcFracFwd();
  double probRC = readExp.gcFracRC();

  if (gcBiasCorrect and probFwd < 0.0) {
    sopt.jointLog->warn("Had no fragments from which to estimate "
                        "fwd vs. rev-comp mapping rate.  Skipping "
                        "sequence-specific / fragment-gc bias correction");
    return effLensIn;
  }

  // calculate read bias normalization factor -- total count in read
  // distribution.
  auto& obs5 = readExp.readBiasModel(salmon::utils::Direction::FORWARD);
  auto& obs3 =
      readExp.readBiasModel(salmon::utils::Direction::REVERSE_COMPLEMENT);
  obs5.normalize();
  obs3.normalize();

  auto& pos5Obs = readExp.posBias(salmon::utils::Direction::FORWARD);
  auto& pos3Obs = readExp.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);

  int32_t K = static_cast<int32_t>(obs5.getContextLength());

  FragmentLengthDistribution& fld = *(readExp.fragmentLengthDistribution());

  // The *expected* biases from GC effects
  auto& transcriptGCDist = readExp.expectedGCBias();
  auto& gcCounts = readExp.observedGC();
  double readGCNormFactor = 0.0;
  int32_t fldLow{0};
  int32_t fldHigh{1};

  // The CDF and PDF of the fragment length distribution
  std::vector<double> cdf(fld.maxVal() + 1, 0.0);
  std::vector<double> pdf(fld.maxVal() + 1, 0.0);
  {
    transcriptGCDist.clear();
    transcriptGCDist.resize(101, 0.0);

    bool lb{false};
    bool ub{false};
    for (size_t i = 0; i <= fld.maxVal(); ++i) {
      pdf[i] = std::exp(fld.pmf(i));
      cdf[i] = (i > 0) ? cdf[i - 1] + pdf[i] : pdf[i];
      auto density = cdf[i];

      if (!lb and density >= 0.005) {
        lb = true;
        fldLow = i;
      }
      if (!ub and density >= 0.995) {
        ub = true;
        fldHigh = i;
      }
    }

    if (gcBiasCorrect) {
      for (auto& c : gcCounts) {
        readGCNormFactor += c;
      }
    }
  }

  // Make this const so there are no shenanigans
  const auto& transcripts = readExp.transcripts();

  double minObservedLength = effLensIn.minCoeff();

  // The effective lengths adjusted for bias
  Eigen::VectorXd effLensOut(effLensIn.size());

  // How much to cut off
  int32_t trunc = K;

  using GCBiasVecT = std::vector<double>;
  using SeqBiasVecT = std::vector<double>;

  //
  // These will store "thread local" parameters
  // for the appropriate bias terms.
  //
  class CombineableBiasParams {
  public:
    CombineableBiasParams(uint32_t K) {
      expectGC = std::vector<double>(101, 0.0);
      expectPos5 = std::vector<SimplePosBias>(5);
      expectPos3 = std::vector<SimplePosBias>(5);
    }

    std::vector<SimplePosBias> expectPos5;
    std::vector<SimplePosBias> expectPos3;
    SBModel expectSeqFW;
    SBModel expectSeqRC;
    std::vector<double> expectGC;
  };

  auto revComplement = [](const char* s, int32_t l, std::string& o) -> void {
    if (l > o.size()) {
      o.resize(l, 'A');
    }
    int32_t j = 0;
    for (int32_t i = l - 1; i >= 0; --i, ++j) {
      switch (s[i]) {
      case 'A':
      case 'a':
        o[j] = 'T';
        break;
      case 'C':
      case 'c':
        o[j] = 'G';
        break;
      case 'T':
      case 't':
        o[j] = 'A';
        break;
      case 'G':
      case 'g':
        o[j] = 'C';
        break;
      default:
        o[j] = 'N';
        break;
      }
    }
  };

  //
  // The local bias terms from each thread can be combined
  // via simple summation.
  //
  auto getBiasParams = [K]() -> CombineableBiasParams {
    return CombineableBiasParams(K);
  };
  tbb::combinable<CombineableBiasParams> expectedDist(getBiasParams);
  std::atomic<size_t> numBackgroundTranscripts{0};
  std::atomic<size_t> numExpressedTranscripts{0};

  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(transcripts.size())),
      [&](const BlockedIndexRange& range) -> void {

        auto& expectSeqFW = expectedDist.local().expectSeqFW;
        auto& expectSeqRC = expectedDist.local().expectSeqRC;
        auto& expectGC = expectedDist.local().expectGC;
        auto& expectPos5 = expectedDist.local().expectPos5;
        auto& expectPos3 = expectedDist.local().expectPos3;

        std::string rcSeq;
        // For each transcript
        for (auto it : boost::irange(range.begin(), range.end())) {

          // Get the transcript
          const auto& txp = transcripts[it];

          // Get the reference length and the
          // "initial" effective length (not considering any biases)
          int32_t refLen = static_cast<int32_t>(txp.RefLength);
          int32_t elen = static_cast<int32_t>(txp.EffectiveLength);

          // The difference between the actual and effective length
          int32_t unprocessedLen = std::max(0, refLen - elen);

          // Skip transcripts with trivial expression or that are too
          // short
          if (alphas[it] < minAlpha or unprocessedLen <= 0) {  // or
txp.uniqueUpdateFraction() < 0.90) {
            if (alphas[it] >= minAlpha) {
              ++numExpressedTranscripts;
            }
            continue;
          }
          ++numBackgroundTranscripts;

          // Otherwise, proceed giving this transcript the following weight
          double weight = (alphas[it] / effLensIn(it));

          // This transcript's sequence
          const char* tseq = txp.Sequence();
          revComplement(tseq, refLen, rcSeq);
          const char* rseq = rcSeq.c_str();

          Mer fwmer;
          fwmer.from_chars(tseq);
          Mer rcmer;
          rcmer.from_chars(rseq);
          int32_t contextLength{expectSeqFW.getContextLength()};

          // For each position along the transcript
          // Starting from the 5' end and moving toward the 3' end
          for (int32_t fragStartPos = 0; fragStartPos < refLen - K;
               ++fragStartPos) {
            // Seq-specific bias
            if (seqBiasCorrect) {
              int32_t contextEndPos =
                  fragStartPos + K - 1; // -1 because pos is *inclusive*

              if (contextEndPos >= 0 and contextEndPos < refLen) {
                int32_t maxFragLen =
                    refLen - (fragStartPos + expectSeqFW.contextBefore(false));
                if (maxFragLen >= 0 and maxFragLen < refLen) {
                  auto cdensity =
                      (maxFragLen >= cdf.size()) ? 1.0 : cdf[maxFragLen];
                  expectSeqFW.addSequence(fwmer, weight * cdensity);
                  expectSeqRC.addSequence(rcmer, weight * cdensity);
                }
              }

              // shift the context one nucleotide to the right
              fwmer.shift_left(tseq[fragStartPos + contextLength]);
              rcmer.shift_left(rseq[fragStartPos + contextLength]);
            } // end: Seq-specific bias

            // fragment-GC bias
            if (gcBiasCorrect) {
              size_t sp = static_cast<size_t>((fldLow > 0) ? fldLow - 1 : 0);
              double prevFLMass = cdf[sp];
              int32_t fragStart = fragStartPos;
              for (int32_t fl = fldLow; fl <= fldHigh; fl += gcSamp) {
                int32_t fragEnd = fragStart + fl - 1;
                if (fragEnd < refLen) {
                  // The GC fraction for this putative fragment
                  auto gcFrac = txp.gcFrac(fragStart, fragEnd);
                  expectGC[gcFrac] += weight * (cdf[fl] - prevFLMass);
                  prevFLMass = cdf[fl];
                } else {
                  break;
                } // no more valid positions
              }   // end: for each fragment length
            }     // end: fragment GC bias

            // positional bias
            if (posBiasCorrect) {
              int32_t maxFragLenFW = refLen - fragStartPos + 1;
              int32_t maxFragLenRC = fragStartPos;
              auto densityFW =
                  (maxFragLenFW >= cdf.size()) ? 1.0 : cdf[maxFragLenFW];
              auto densityRC =
                  (maxFragLenRC >= cdf.size()) ? 1.0 : cdf[maxFragLenRC];
              if (weight * densityFW > 1e-8) {
                expectPos5[txp.lengthClassIndex()].addMass(
                    fragStartPos, txp.RefLength, std::log(weight * densityFW));
              }
              if (weight * densityRC > 1e-8) {
                expectPos3[txp.lengthClassIndex()].addMass(
                    fragStartPos, txp.RefLength, std::log(weight * densityRC));
              }
            }
          } // end: for every fragment start position
        }   // end for each transcript

      } // end tbb for function
      );

  size_t bgCutoff =
      std::min(static_cast<size_t>(150),
               static_cast<size_t>(numBackgroundTranscripts *  0.1));
  if (numBackgroundTranscripts < bgCutoff) {
    sopt.jointLog->warn("I found only {} transcripts meeting the necessary "
                        "conditions to contribute to "
                        "the bias background distribution.  This is likely too "
                        "small to safely do bias correction. "
                        "I'm skipping bias correction",
                        numBackgroundTranscripts.load());
    sopt.biasCorrect = false;
    sopt.gcBiasCorrect = false;
    sopt.posBiasCorrect = false;
    return effLensIn;
  }

   //
   // The local bias terms from each thread can be combined
   // via simple summation.  Here, we combine the locally-computed
   // bias terms.
   //
  SBModel exp5;
  SBModel exp3;
  std::vector<SimplePosBias> pos5Exp(5);
  std::vector<SimplePosBias> pos3Exp(5);
  auto combineBiasParams =
      [seqBiasCorrect, gcBiasCorrect, posBiasCorrect, &pos5Exp, &pos3Exp, &exp5,
       &exp3, &transcriptGCDist](const CombineableBiasParams& p) -> void {
    if (seqBiasCorrect) {
      exp5.combineCounts(p.expectSeqFW);
      exp3.combineCounts(p.expectSeqRC);
    }
    if (gcBiasCorrect) {
      for (size_t i = 0; i < p.expectGC.size(); ++i) {
        transcriptGCDist[i] += p.expectGC[i];
      }
    }
    if (posBiasCorrect) {
      for (size_t i = 0; i < p.expectPos5.size(); ++i) {
        pos5Exp[i].combine(p.expectPos5[i]);
        pos3Exp[i].combine(p.expectPos3[i]);
      }
    }
  };
  expectedDist.combine_each(combineBiasParams);

  // finalize expected positional biases
  if (posBiasCorrect) {
    for (size_t i = 0; i < pos5Exp.size(); ++i) {
      pos5Exp[i].finalize();
      pos3Exp[i].finalize();
    }
  }

  auto smoothDist = [](std::vector<double>& v, int w, int d) -> void {
    v = sg_smooth(v, w, d);
    double gcSum = 0.0;
    for (size_t i = 0; i < v.size(); ++i) {
      if (v[i] < 1e-5) { v[i] = 1e-5; }
      gcSum += v[i];
    }
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] /= gcSum;
    }
  };

  sopt.jointLog->info("Computed expected counts (for bias correction)");

  // Compute appropriate priors and normalization factors
  double txomeGCNormFactor = 0.0;
  double gcPrior = 0.0;
  if (gcBiasCorrect) {
    for (auto m : transcriptGCDist) {
      txomeGCNormFactor += m;
    }
    auto pmass = 1e-5 * 101.0;
    gcPrior =
        ((pmass / (readGCNormFactor - pmass)) * txomeGCNormFactor) / 101.0;
    txomeGCNormFactor += gcPrior * 101.0;
  }
  double scaleGCBias = (txomeGCNormFactor / readGCNormFactor);

  // Compute the bias weights for each fragment-GC bin
  Eigen::VectorXd gcBias(101);
  double gcBiasMax = 100.0;
  double gcBiasMin = 1.0 / gcBiasMax;
  if (gcBiasCorrect) {
    for (size_t i = 0; i < 101; ++i) {
      gcBias[i] = scaleGCBias * (gcCounts[i] / (gcPrior + transcriptGCDist[i]));
      gcBias[i] = (gcBias[i] > gcBiasMax) ? gcBiasMax :
          ((gcBias[i] < gcBiasMin) ? gcBiasMin : gcBias[i]);
    }
    //gcBias *= scaleGCBias;
  }


  exp5.normalize();
  exp3.normalize();

  bool noThreshold = sopt.noBiasLengthThreshold;
  std::atomic<size_t> numCorrected{0};
  std::atomic<size_t> numUncorrected{0};

  // Write out the bias model parameters we learned
  if (writeBias) {
    boost::filesystem::path auxDir = sopt.outputDirectory / sopt.auxDir;
    bool auxSuccess = boost::filesystem::is_directory(auxDir);
    if (!auxSuccess) {
      auxSuccess = boost::filesystem::create_directories(auxDir);
    }
    if (auxSuccess) {
      auto exp5fn = auxDir / "exp5_marginals.txt";
      auto& exp5m = exp5.marginals();
      std::ofstream exp5f(exp5fn.string());
      exp5f << exp5m.rows() << '\t' << exp5m.cols() << '\n';
      exp5f << exp5m;
      exp5f.close();

      auto exp3fn = auxDir / "exp3_marginals.txt";
      auto& exp3m = exp3.marginals();
      std::ofstream exp3f(exp3fn.string());
      exp3f << exp3m.rows() << '\t' << exp3m.cols() << '\n';
      exp3f << exp3m;
      exp3f.close();

      auto obs5fnc = auxDir / "obs5_conditionals.txt";
      std::ofstream obs5fc(obs5fnc.string());
      obs5.dumpConditionalProbabilities(obs5fc);
      obs5fc.close();

      auto obs3fnc = auxDir / "obs3_conditionals.txt";
      std::ofstream obs3fc(obs3fnc.string());
      obs3.dumpConditionalProbabilities(obs3fc);
      obs3fc.close();

      auto exp5fnc = auxDir / "exp5_conditionals.txt";
      std::ofstream exp5fc(exp5fnc.string());
      exp5.dumpConditionalProbabilities(exp5fc);
      exp5fc.close();

      auto exp3fnc = auxDir / "exp3_conditionals.txt";
      std::ofstream exp3fc(exp3fnc.string());
      exp3.dumpConditionalProbabilities(exp3fc);
      exp3fc.close();
    } else {
      sopt.jointLog->warn(
          "Couldn't create auxiliary directory {} to write bias parameters",
          auxDir);
    }
  }

  std::atomic<uint32_t> numProcessed{0};
  size_t numTranscripts = transcripts.size();
  size_t stepSize = static_cast<size_t>(transcripts.size() * 0.1);
  size_t nextUpdate{0};

  std::mutex updateMutex;
  // Compute the effective lengths of each transcript (in parallel)
  tbb::parallel_for(
      BlockedIndexRange(size_t(0), size_t(transcripts.size())),
      [&](const BlockedIndexRange& range) -> void {

        std::string rcSeq;
        // For each transcript
        for (auto it : boost::irange(range.begin(), range.end())) {

          auto& txp = transcripts[it];

          // eff. length starts out as 0
          double effLength = 0.0;

          // Reference length
          int32_t refLen = static_cast<int32_t>(txp.RefLength);
          // Effective length before any bias correction
          int32_t elen = static_cast<int32_t>(txp.EffectiveLength);

          // How much of this transcript (beginning and end) should
          // not be considered
          int32_t unprocessedLen = std::max(0, refLen - elen);

          if (alphas[it] >= minAlpha and unprocessedLen > 0) {

            Eigen::VectorXd seqFactorsFW(refLen);
            Eigen::VectorXd seqFactorsRC(refLen);
            seqFactorsFW.setOnes();
            seqFactorsRC.setOnes();

            std::vector<double> posFactorsFW(refLen, 1.0);
            std::vector<double> posFactorsRC(refLen, 1.0);

            // This transcript's sequence
            const char* tseq = txp.Sequence();
            revComplement(tseq, refLen, rcSeq);
            const char* rseq = rcSeq.c_str();

            int32_t fl = fldLow;
            auto maxLen = std::min(refLen, fldHigh + 1);
            bool done{fl >= maxLen};

            if (posBiasCorrect) {
              std::vector<double> posFactorsObs5(refLen, 1.0);
              std::vector<double> posFactorsObs3(refLen, 1.0);
              std::vector<double> posFactorsExp5(refLen, 1.0);
              std::vector<double> posFactorsExp3(refLen, 1.0);
              auto li = txp.lengthClassIndex();
              auto& p5O = pos5Obs[li];
              auto& p3O = pos3Obs[li];
              auto& p5E = pos5Exp[li];
              auto& p3E = pos3Exp[li];
              p5O.projectWeights(posFactorsObs5);
              p3O.projectWeights(posFactorsObs3);
              p5E.projectWeights(posFactorsExp5);
              p3E.projectWeights(posFactorsExp3);
              for (int32_t fragStart = 0; fragStart < refLen - K; ++fragStart) {
                posFactorsFW[fragStart] =
                    posFactorsObs5[fragStart] / posFactorsExp5[fragStart];
                posFactorsRC[fragStart] =
                    posFactorsObs3[fragStart] / posFactorsExp3[fragStart];
              }
            }

            // Evaluate the sequence specific bias (5' and 3') over the length
            // of the transcript.  After this loop,
            // seqFactorsFW will contain the sequence-specific bias for each
            // position on the 5' strand
            // and seqFactorsRC will contain the sequence-specific bias for each
            // position on the 3' strand.
            if (seqBiasCorrect) {
              Mer mer;
              Mer rcmer;
              mer.from_chars(tseq);
              rcmer.from_chars(rseq);
              int32_t contextLength{exp5.getContextLength()};

              for (int32_t fragStart = 0; fragStart < refLen - K; ++fragStart) {
                int32_t readStart = fragStart + obs5.contextBefore(false);
                int32_t kmerEndPos =
                    fragStart + K - 1; // -1 because pos is *inclusive*

                if (kmerEndPos >= 0 and kmerEndPos < refLen and
                    readStart < refLen) {
                  seqFactorsFW[readStart] =
                      std::exp(obs5.evaluateLog(mer) - exp5.evaluateLog(mer));
                  seqFactorsRC[readStart] = std::exp(obs3.evaluateLog(rcmer) -
                                                     exp3.evaluateLog(rcmer));
                }
                // shift the context one nucleotide to the right
                mer.shift_left(tseq[fragStart + contextLength]);
                rcmer.shift_left(rseq[fragStart + contextLength]);
              }
              // We need these in 5' -> 3' order, so reverse them
              seqFactorsRC.reverseInPlace();
            } // end sequence-specific factor calculation

            if (numProcessed > nextUpdate) {
              updateMutex.try_lock();
          if (numProcessed > nextUpdate) {
        sopt.jointLog->info(
                    "processed bias for {:3.1f}% of the transcripts",
                    100.0 * (numProcessed /
static_cast<double>(numTranscripts)));
        nextUpdate += stepSize;
        if (nextUpdate > numTranscripts) {
          nextUpdate = numTranscripts - 1;
        }
          }
          updateMutex.unlock();
        }
        ++numProcessed;

        size_t sp = static_cast<size_t>((fl > 0) ? fl - 1 : 0);
            double prevFLMass = cdf[sp];
            double unbiasedMass{0.0};

            // For every possible fragment length
            while (!done) {
              if (fl >= maxLen) {
                done = true;
                fl = maxLen - 1;
              }
              double flWeight = cdf[fl] - prevFLMass;
              prevFLMass = cdf[fl];

              double flMassTotal{0.0};
              // For every position a fragment of length fl could start
              for (int32_t kmerStartPos = 0; kmerStartPos < refLen - fl;
                   ++kmerStartPos) {
                int32_t fragStart = kmerStartPos;
                int32_t fragEnd = fragStart + fl - 1;

                // If the 3' end is within the transcript
                if (fragStart < refLen and fragEnd < refLen) {
                  double fragFactor =
                      seqFactorsFW[fragStart] * seqFactorsRC[fragEnd];
                  if (gcBiasCorrect) {
                    auto gcFrac = txp.gcFrac(fragStart, fragEnd);
                    fragFactor *= gcBias[gcFrac];
                  }
                  if (posBiasCorrect) {
                    fragFactor *=
                        posFactorsFW[fragStart] * posFactorsRC[fragEnd];
                  }
                  flMassTotal += fragFactor;
                } else {
                  break;
                }
              }

              effLength += (flWeight * flMassTotal);
              fl += gcSamp;
            }
            // effLength = flMassTotal;//flMasses.sum();

          } // for the processed transcript

          // throw caution to the wind
          double thresh = noThreshold ? 1.0 : unprocessedLen;

          // JUNE 17
          if (noThreshold) {
              if (unprocessedLen > 0.0 and effLength > thresh) {
                  effLensOut(it) = effLength;
              } else {
                  effLensOut(it) = effLensIn(it);
              }
          } else {
              double offset = std::max(1.0, thresh);
              double effLengthNoBias = static_cast<double>(elen);
              auto barrierLength = [effLengthNoBias, offset](double x) -> double
{
                  return std::max(x, std::min(effLengthNoBias, offset));
                  //return x + ((unprocessedLen * unprocessedLen) /
                  //            (x + unprocessedLen));
              };
              effLength = barrierLength(effLength);
              effLensOut(it) = effLength;
          }

          // END: JUNE 17
//
//          // commented out JUNE 20
//          // To correct the transcript length, we require it to be
//          // "sufficiently" long to begin with.
//          // WORKING MODEL: if (unprocessedLen > 0.0 and elen > thresh and
effLength > thresh) {
//	  if (unprocessedLen > 0.0 and effLength > thresh) {
//            ++numCorrected;
//            effLensOut(it) = effLength;
//          } else {
//            ++numUncorrected;
//            effLensOut(it) = effLensIn(it);
//          }



        }
      } // end parallel_for lambda
      );
  sopt.jointLog->info("processed bias for 100.0% of the transcripts");
  return effLensOut;
}
*/
