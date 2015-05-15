#include <boost/thread/thread.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/range/join.hpp>

#include "SalmonUtils.hpp"
#include "AlignmentLibrary.hpp"
#include "ReadPair.hpp"
#include "UnpairedRead.hpp"
#include "SalmonMath.hpp"
#include "LibraryFormat.hpp"
#include "ReadExperiment.hpp"

#include "spdlog/spdlog.h"

#include "gff.h"

#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "TranscriptGeneMap.hpp"
#include "GenomicFeature.hpp"

namespace salmon {
namespace utils {

    bool headersAreConsistent(SAM_hdr* h1, SAM_hdr* h2) {

        bool consistent{true};
        // Both files must contain the same number of targets
        if (h1->nref != h2->nref) { consistent = false; }

        // Check each target to ensure that the name and length are the same.
        size_t i = 0;
        size_t n = h1->nref;
        while (consistent and i < n) {
            size_t l1 = h1->ref[i].len;
            size_t l2 = h2->ref[i].len;
            consistent = (l1 == l2) and
                         (strcmp(h1->ref[i].name, h2->ref[i].name) == 0);
            ++i;
        }

        return consistent;
    }

    bool headersAreConsistent(std::vector<SAM_hdr*>&& headers) {
        if (headers.size() == 1) { return true; }

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

    double logAlignFormatProb(const LibraryFormat observed, const LibraryFormat expected, double incompatPrior) {
        // Allow orphaned reads in a paired-end library, but
        // decrease their a priori probability.
        if (expected.type == ReadType::PAIRED_END and
            observed.type == ReadType::SINGLE_END) {
            double logOrphanProb = salmon::math::LOG_ORPHAN_PROB;
            if (expected.strandedness == ReadStrandedness::U or
                expected.strandedness == ReadStrandedness::AS or
                expected.strandedness == ReadStrandedness::SA) {
                return salmon::math::LOG_1;
            } else {
                return (expected.strandedness == observed.strandedness) ? logOrphanProb : incompatPrior;
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
		    std::cerr << "expected = " << expected << "\n";
		    std::cerr << "observed = " << observed << "\n";
                    return incompatPrior;
                }
            }
        }

        fmt::print(stderr, "WARNING: logAlignFormatProb --- should not get here");
        return salmon::math::LOG_0;
    }

    template <typename ExpLib>
    void writeAbundances(const SalmonOpts& sopt,
                         ExpLib& alnLib,
                         boost::filesystem::path& fname,
                         std::string headerComments) {
        using salmon::math::LOG_0;
        using salmon::math::LOG_1;

        std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(fname.c_str(), "w"), std::fclose);

        fmt::print(output.get(), "{}", headerComments);
        fmt::print(output.get(), "# Name\tLength\tTPM\tFPKM\tNumReads\n");


        auto& refs = alnLib.transcripts();
        auto numMappedReads = alnLib.numMappedReads();
        const double logBillion = std::log(1000000000.0);
        const double million = 1000000.0;
        const double logNumFragments = std::log(static_cast<double>(numMappedReads));
        const double upperBoundFactor = static_cast<double>(alnLib.upperBoundHits()) /
                                        numMappedReads;
        auto clusters = alnLib.clusterForest().getClusters();
        size_t clusterID = 0;
        for(auto cptr : clusters) {

            double logClusterMass = cptr->logMass();
            double logClusterCount = std::log(upperBoundFactor * static_cast<double>(cptr->numHits()));

            if (logClusterMass == LOG_0) {
                std::cerr << "Warning: cluster " << clusterID << " has 0 mass!\n";
            }

            bool requiresProjection{false};

            auto& members = cptr->members();
            size_t clusterSize{0};
            for (auto transcriptID : members) {
                Transcript& t = refs[transcriptID];
                t.uniqueCounts = t.uniqueCount();
                t.totalCounts = t.totalCount();
                ++clusterSize;
            }

            for (auto transcriptID : members) {
                Transcript& t = refs[transcriptID];
                double logTranscriptMass = t.mass(false);
                if (logTranscriptMass == LOG_0) {
                    t.projectedCounts = 0;
                } else {
                    double logClusterFraction = logTranscriptMass - logClusterMass;
                    t.projectedCounts = std::exp(logClusterFraction + logClusterCount);
                    requiresProjection |= t.projectedCounts > static_cast<double>(t.totalCounts) or
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
            double refLength = sopt.noEffectiveLengthCorrection ?
                               transcript.RefLength :
                               std::exp(transcript.getCachedEffectiveLength());
            //refLength = transcript.RefLength;
            tfracDenom += (transcript.projectedCounts / numMappedReads) / refLength;
        }

        // Now posterior has the transcript fraction
        for (auto& transcript : transcripts_) {
            double logLength = sopt.noEffectiveLengthCorrection ?
                               std::log(transcript.RefLength) :
                               transcript.getCachedEffectiveLength();
            //logLength = std::log(transcript.RefLength);
            double fpkmFactor = std::exp(logBillion - logLength - logNumFragments);
            double count = transcript.projectedCounts;
            //double countTotal = transcripts_[transcriptID].totalCounts;
            //double countUnique = transcripts_[transcriptID].uniqueCounts;
            double fpkm = count > 0 ? fpkmFactor * count : 0.0;
            double npm = (transcript.projectedCounts / numMappedReads);
            double refLength = std::exp(logLength);
            double tfrac = (npm / refLength) / tfracDenom;
            double tpm = tfrac * million;

            fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n",
                    transcript.RefName, transcript.RefLength,
                    tpm, fpkm, count);
        }

    }

    LibraryFormat hitType(int32_t end1Start, bool end1Fwd,
                          int32_t end2Start, bool end2Fwd) {

        // If the reads come from opposite strands
        if (end1Fwd != end2Fwd) {
            // and if read 1 comes from the forward strand
            if (end1Fwd) {
                // then if read 1 start < read 2 start ==> ISF
                if (end1Start <= end2Start) {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);
                } // otherwise read 2 start < read 1 start ==> OSF
                else {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);
                }
            }
            // and if read 2 comes from the forward strand
            if (end2Fwd) {
                // then if read 2 start <= read 1 start ==> ISR
                if (end2Start <= end1Start) {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);
                } // otherwise, read 2 start > read 1 start ==> OSR
                else {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);
                }
            }
        } else { // Otherwise, the reads come from the same strand
            if (end1Fwd) { // if it's the forward strand ==> MSF
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);
            } else { // if it's the reverse strand ==> MSR
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);
            }
        }
        // SHOULD NOT GET HERE
        spdlog::get("jointLog")->error("ERROR: Could not associate any known library type with read! "
                                       "Please report this bug!\n");
        std::exit(-1);
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U);
    }



    LibraryFormat hitType(int32_t end1Start, bool end1Fwd, uint32_t len1,
                          int32_t end2Start, bool end2Fwd, uint32_t len2, bool canDovetail) {

        // If the reads come from opposite strands
        if (end1Fwd != end2Fwd) {
            // and if read 1 comes from the forward strand
            if (end1Fwd) {
                // then if read 1 start < read 2 start ==> ISF
                // NOTE: We can't really delineate between inward facing reads that stretch
                // past each other and outward facing reads --- the purpose of stretch is to help
                // make this determinateion.
                int32_t stretch = canDovetail ? len2 : 0;
                if (end1Start <= end2Start + stretch) {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);
                } // otherwise read 2 start < read 1 start ==> OSF
                else {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);
                }
            }
            // and if read 2 comes from the forward strand
            if (end2Fwd) {
                // then if read 2 start <= read 1 start ==> ISR
                // NOTE: We can't really delineate between inward facing reads that stretch
                // past each other and outward facing reads --- the purpose of stretch is to help
                // make this determinateion.
                int32_t stretch = canDovetail ? len1 : 0;
                if (end2Start <= end1Start + stretch) {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);
                } // otherwise, read 2 start > read 1 start ==> OSR
                else {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);
                }
            }
        } else { // Otherwise, the reads come from the same strand
            if (end1Fwd) { // if it's the forward strand ==> MSF
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);
            } else { // if it's the reverse strand ==> MSR
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);
            }
        }
        // SHOULD NOT GET HERE
        spdlog::get("jointLog")->error("ERROR: Could not associate any known library type with read! "
                                       "Please report this bug!\n");
        std::exit(-1);
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U);
    }


    LibraryFormat hitType(int32_t start, bool isForward) {
        // If the read comes from the forward strand
        if (isForward) {
            return LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S);
        } else {
            return LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A);
        }
        // SHOULD NOT GET HERE
        fmt::print(stderr, "WARNING: Could not associate known library type with read!\n");
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U);

    }

using std::string;
using NameVector = std::vector<string>;
using IndexVector = std::vector<size_t>;
using KmerVector = std::vector<uint64_t>;

/**
 * This function parses the library format string that specifies the format in which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatStringNew(std::string& fmt) {
	using std::vector;
	using std::string;
	using std::map;
	using std::stringstream;

    map<string, LibraryFormat> formatMap = {
        {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U)},
        {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA)},
        {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS)},
        {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)},
        {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA)},
        {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS)},
        {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)},
        {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)},
        {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)},
        {"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U)},
        {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S)},
        {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A)}};

	// inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
	// first convert the string to upper-case
	for (auto& c : fmt) { c = std::toupper(c); }


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
std::vector<ReadLibrary> extractReadLibraries(boost::program_options::parsed_options& orderedOptions) {
	// The current (default) format for paired end data
	LibraryFormat peFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U);
	// The current (default) format for single end data
	LibraryFormat seFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U);

	std::vector<ReadLibrary> peLibs{ReadLibrary(peFormat)};
	std::vector<ReadLibrary> seLibs{ReadLibrary(seFormat)};
	for (auto& opt : orderedOptions.options) {
		// Update the library type
		if (opt.string_key == "libType") {
			auto libFmt = parseLibraryFormatStringNew(opt.value[0]);
			if (libFmt.type == ReadType::PAIRED_END) {
				peFormat = libFmt;
				peLibs.emplace_back(libFmt);
			} else {
				seFormat = libFmt;
				seLibs.emplace_back(libFmt);
			}
		}
		if (opt.string_key == "mates1") {
			peLibs.back().addMates1(opt.value);
		}
		if (opt.string_key == "mates2") {
			peLibs.back().addMates2(opt.value);
		}
		if (opt.string_key == "unmatedReads") {
			seLibs.back().addUnmated(opt.value);
		}
	}

	std::vector<ReadLibrary> libs;
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
	std::cerr << "there " << ((numLibs > 1) ? "are " : "is ") << libs.size() << ((numLibs > 1) ? " libs\n" : " lib\n");
	return libs;
}



/**
 * This function parses the library format string that specifies the format in which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatString(std::string& fmt) {
    using std::vector;
    using std::string;
    using std::map;
    using std::stringstream;

    // inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c

    // first convert the string to upper-case
    for (auto& c : fmt) { c = std::toupper(c); }
    // split on the delimiter ':', and put the key, value (k=v) pairs into a map
    stringstream ss(fmt);
    string item;
    map<string, string> kvmap;
    while (std::getline(ss, item, ':')) {
        auto splitPos = item.find('=', 0);
        string key{item.substr(0, splitPos)};
        string value{item.substr(splitPos+1)};
        kvmap[key] = value;
    }

    map<string, ReadType> readType = {{"SE", ReadType::SINGLE_END}, {"PE", ReadType::PAIRED_END}};
    map<string, ReadOrientation> orientationType = {{">>", ReadOrientation::SAME},
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

    ReadType type = (typeStr == "SE") ? ReadType::SINGLE_END : ReadType::PAIRED_END;
    ReadOrientation orientation = (type == ReadType::SINGLE_END) ? ReadOrientation::NONE : ReadOrientation::TOWARD;
    ReadStrandedness strandedness{ReadStrandedness::U};
    // Construct the LibraryFormat class from the key, value map
    for (auto& kv : kvmap) {
        auto& k = kv.first; auto& v = kv.second;
        if (k == "O" or k == "ORIENTATION") {
            auto it = orientationType.find(v);
            if (it != orientationType.end()) { orientation = orientationType[it->first]; } else {
                string e = v + " is not a valid orientation type; must be one of {>>, <>, ><}";
                throw std::invalid_argument(e);
            }

        }
        if (k == "S" or k == "STRAND") {
            auto it = strandType.find(v);
            if (it != strandType.end()) { strandedness = strandType[it->first]; } else {
                string e = v + " is not a valid strand type; must be one of {SA, AS, S, A, U}";
                throw std::invalid_argument(e);
            }
        }

    }
    LibraryFormat lf(type, orientation, strandedness);
    return lf;
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
    while (ifile >> s) { if (s.front() == '>') { ++numReads; } }

    ifile.close();

    return numReads;
}

bool readKmerOrder( const std::string& fname, std::vector<uint64_t>& kmers ) {

  std::ifstream mlist(fname, std::ios::in | std::ios::binary);
  // Get the number of kmers from file
  size_t numKmers{0};
  mlist.read( reinterpret_cast<char*>( &numKmers ), sizeof( size_t ) );

  // Resize the array that will hold the sorted kmers
  kmers.resize(numKmers, 0);
  mlist.read( reinterpret_cast<char*>( &kmers[0] ), sizeof( uint64_t) * kmers.size() );

  mlist.close();

  return true;
}

template <template<typename> class S, typename T>
bool overlap( const S<T> &a, const S<T> &b ) {
    // Query from the smaller set to the larger set
    if ( a.size() <= b.size() ) {
        for ( auto & ae : a ) {
            if (b.find(ae) != b.end()) {
                return true;
            }
        }
    } else {
        for ( auto & be : b ) {
            if (a.find(be) != b.end()) {
                return true;
            }
        }
    }
    // If nothing from the smaller set is in the larger set, then they don't overlap
    return false;
}


TranscriptGeneMap transcriptGeneMapFromGTF(const std::string& fname, std::string key) {

    using std::unordered_set;
    using std::unordered_map;
    using std::vector;
    using std::tuple;
    using std::string;
    using std::get;

    // Use GffReader to read the file
    GffReader reader(const_cast<char*>(fname.c_str()));
    // Remember the optional attributes
    reader.readAll(true);

    struct TranscriptKeyPair {
        const char* transcript_id;
        const char* key;
        TranscriptKeyPair(const char* t, const char* k) :
            transcript_id(t), key(k) {}
    };

    // The user can group transcripts by gene_id, gene_name, or
    // an optinal attribute that they provide as a string.
    enum class TranscriptKey { GENE_ID, GENE_NAME, DYNAMIC };

    // Select the proper attribute by which to group
    TranscriptKey tkey = TranscriptKey::GENE_ID;

    if (key == "gene_id") {
    } else if (key == "gene_name") {
        tkey = TranscriptKey::GENE_NAME;
    } else {
        tkey = TranscriptKey::DYNAMIC;
    }

    // Iterate over all transcript features and build the
    // transcript <-> key vector.
    auto nfeat = reader.gflst.Count();
    std::vector<TranscriptKeyPair> feats;
    for (int i=0; i < nfeat; ++i) {
        auto f = reader.gflst[i];
        if (f->isTranscript()) {
            const char* keyStr;
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
            feats.emplace_back(f->getID(), keyStr);
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

    std::sort( feats.begin(), feats.end(),
    []( const TranscriptKeyPair & a, const TranscriptKeyPair & b) -> bool {
        return std::strcmp(a.transcript_id, b.transcript_id) < 0;
    } );

    std::string currentTranscript = "";
    for ( auto & feat : feats ) {

        std::string gene(feat.key);
        std::string transcript(feat.transcript_id);

        if ( transcript != currentTranscript ) {
            auto geneIt = geneNameToID.find(gene);
            size_t geneID = 0;

            if ( geneIt == geneNameToID.end() ) {
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


TranscriptGeneMap readTranscriptToGeneMap( std::ifstream &ifile ) {

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

    while ( ifile >> transcript >> gene ) {
        // The transcript and it's ID
        transcripts.push_back( make_tuple(transcript, transcriptCounter) );

        auto geneIt = geneNameToID.find(gene);
        size_t geneID = 0;

        if ( geneIt == geneNameToID.end() ) {
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

    std::sort( transcripts.begin(), transcripts.end(),
               []( const NameID & a, const NameID & b) -> bool { return get<0>(a) < get<0>(b); } );

    // Resize these vectors for fast access
    transcriptNames.resize(t2gUnordered.size());
    t2g.resize(t2gUnordered.size());

    for ( size_t newID = 0; newID < transcripts.size(); ++newID ) {
        // For each transcript, map it to the appropriate gene
        string oldName; size_t oldID;
        std::tie(oldName, oldID) = transcripts[newID];
        t2g[newID] = t2gUnordered[oldID];
        transcriptNames[newID] = oldName;
    }

    return TranscriptGeneMap(transcriptNames, geneNames, t2g);
}


TranscriptGeneMap transcriptToGeneMapFromFasta( const std::string& transcriptsFile ) {
    using std::vector;
    using stream_manager = jellyfish::stream_manager<char**>;
    using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;
    namespace bfs = boost::filesystem;

    NameVector transcriptNames;
    NameVector geneNames {"gene"};

    vector<bfs::path> paths{transcriptsFile};

    // Create a jellyfish parser
    const int concurrentFile{1};
    char** fnames = new char*[1];
    fnames[0] = const_cast<char*>(transcriptsFile.c_str());
    stream_manager streams(fnames, fnames + 1, concurrentFile);

    size_t maxReadGroupSize{100};
    sequence_parser parser(4, maxReadGroupSize, concurrentFile, streams);

    // while there are transcripts left to process
    while (true) {
        sequence_parser::job j(parser);
        // If this job is empty, then we're done
        if (j.is_empty()) { break; }

        for (size_t i=0; i < j->nb_filled; ++i) {
            // The transcript name
            std::string fullHeader(j->data[i].header);
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
                         std::vector<double>& expValsIn) :
            target(targetIn), length(lengthIn), expVals(expValsIn) {}

        ExpressionRecord( ExpressionRecord&& other ) {
            std::swap(target, other.target);
            length = other.length;
            std::swap(expVals, other.expVals);
        }

        ExpressionRecord(std::vector<std::string>& inputLine) {
            if (inputLine.size() < 3) {
                std::string err ("Any expression line must contain at least 3 tokens");
                throw std::invalid_argument(err);
            } else {
                auto it = inputLine.begin();
                target = *it; ++it;
                length = std::stoi(*it); ++it;
                for (; it != inputLine.end(); ++it) {
                    expVals.push_back(std::stod(*it));
                }
            }
        }

        std::string target;
        uint32_t length;
        std::vector<double> expVals;
};

// From : http://stackoverflow.com/questions/9435385/split-a-string-using-c11
std::vector<std::string> split(const std::string& str, int delimiter(int) = ::isspace){
    using namespace std;
    vector<string> result;
    auto e=str.end();
    auto i=str.begin();
    while (i != e) {
        i = find_if_not(i,e, delimiter);
        if (i == e) break;
        auto j = find_if(i,e, delimiter);
        result.push_back(string(i,j));
        i = j;
    }
    return result;
}

void aggregateEstimatesToGeneLevel(TranscriptGeneMap& tgm, boost::filesystem::path& inputPath) {

    using std::vector;
    using std::string;
    using std::ofstream;
    using std::unordered_map;
    using std::move;
    using std::cerr;
    using std::max;

    std::ifstream expFile(inputPath.string());

    if (!expFile.is_open()) {
        perror("Error reading file");
    }

    //====================== From GeneSum ====================
    vector<string> comments;
    unordered_map<string, vector<ExpressionRecord>> geneExps;
    string l;
    size_t ln{0};

    while (getline(expFile, l)) {
        if (++ln % 1000 == 0) {
            cerr << "\r\rParsed " << ln << " expression lines";
        }
        auto it = find_if(l.begin(), l.end(),
                    [](char c) -> bool {return !isspace(c);});
        if (it != l.end()) {
            if (*it == '#') {
                comments.push_back(l);
            } else {
                vector<string> toks = split(l);
                ExpressionRecord er(toks);
                auto gn = tgm.geneName(er.target);
                geneExps[gn].push_back(move(er));
            }
        }
    }
    cerr << "\ndone\n";
    expFile.close();

    cerr << "Aggregating expressions to gene level . . .";
    boost::filesystem::path outputFilePath(inputPath);
    outputFilePath.replace_extension(".genes.sf");
    ofstream outFile(outputFilePath.string());

    // preserve any comments in the output
    for (auto& c : comments) {
        outFile << c << '\n';
    }

    for (auto& kv : geneExps) {
        auto& gn = kv.first;

        uint32_t geneLength{kv.second.front().length};
        vector<double> expVals(kv.second.front().expVals.size(), 0);
        const size_t NE{expVals.size()};

        for (auto& tranExp : kv.second) {
            geneLength = max(geneLength, tranExp.length);
            for (size_t i = 0; i < NE; ++i) { expVals[i] += tranExp.expVals[i]; }
        }

        outFile << gn << '\t' << geneLength;
        for (size_t i = 0; i < NE; ++i) {
            outFile << '\t' << expVals[i];
        }
        outFile << '\n';
    }

    outFile.close();
    cerr << " done\n";
    //====================== From GeneSum =====================
}

void generateGeneLevelEstimates(boost::filesystem::path& geneMapPath,
                                boost::filesystem::path& estDir,
                                bool haveBiasCorrectedFile) {
    namespace bfs = boost::filesystem;
    std::cerr << "Computing gene-level abundance estimates\n";
    bfs::path gtfExtension(".gtf");
    auto extension = geneMapPath.extension();

    TranscriptGeneMap tranGeneMap;
    // parse the map as a GTF file
    if (extension == gtfExtension) {
        // Using libgff
        tranGeneMap = salmon::utils::transcriptGeneMapFromGTF(geneMapPath.string(), "gene_id");
    } else { // parse the map as a simple format files
        std::ifstream tgfile(geneMapPath.string());
        tranGeneMap = salmon::utils::readTranscriptToGeneMap(tgfile);
        tgfile.close();
    }

    std::cerr << "There were " << tranGeneMap.numTranscripts() << " transcripts mapping to "
        << tranGeneMap.numGenes() << " genes\n";

    bfs::path estFilePath = estDir / "quant.sf";
    if (!bfs::exists(estFilePath)) {
        std::stringstream errstr;
        errstr << "Attempting to compute gene-level esimtates, but could not \n"
            << "find isoform-level file " << estFilePath;
        throw std::invalid_argument(errstr.str());
    } else {
        salmon::utils::aggregateEstimatesToGeneLevel(tranGeneMap, estFilePath);
    }

    /** Create a gene-level summary of the bias-corrected estimates as well if these exist **/
    if (haveBiasCorrectedFile) {
        bfs::path biasCorrectEstFilePath = estDir / "quant_bias_corrected.sf";
        if (!bfs::exists(biasCorrectEstFilePath)) {
            std::stringstream errstr;
            errstr << "Attempting to compute gene-level esimtates, but could not \n"
                << "find bias-corrected isoform-level file " << biasCorrectEstFilePath;
            throw std::invalid_argument(errstr.str());
        } else {
            salmon::utils::aggregateEstimatesToGeneLevel(tranGeneMap, biasCorrectEstFilePath);
        }
    }
}

}
}

template
void salmon::utils::writeAbundances<AlignmentLibrary<ReadPair>>(
                                              const SalmonOpts& opts,
                                              AlignmentLibrary<ReadPair>& alnLib,
                                              boost::filesystem::path& fname,
                                              std::string headerComments);

template
void salmon::utils::writeAbundances<AlignmentLibrary<UnpairedRead>>(
                                                  const SalmonOpts& opts,
                                                  AlignmentLibrary<UnpairedRead>& alnLib,
                                                  boost::filesystem::path& fname,
                                                  std::string headerComments);
template
void salmon::utils::writeAbundances<ReadExperiment>(
                                                  const SalmonOpts& opts,
                                                  ReadExperiment& alnLib,
                                                  boost::filesystem::path& fname,
                                                  std::string headerComments);

