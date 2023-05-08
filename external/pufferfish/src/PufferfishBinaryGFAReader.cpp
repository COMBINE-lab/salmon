#include "PufferfishBinaryGFAReader.hpp"
#include "CanonicalKmer.hpp"
#include "cereal/archives/binary.hpp"
#include "xxhash.h"
#include "Kmer.hpp"
#include "string_view.hpp"
#include <chrono>
#include <algorithm>
#include <string>
#include <bitset>


namespace kmers = combinelib::kmers;

namespace pufferfish {

    enum class Direction : bool {
        PREPEND, APPEND
    };


    namespace gfa_reader {
        namespace detail {
#define R -1
#define I -2
#define O -3
#define A 3
#define C 0
#define G 1
#define T 2
            static constexpr int shift_table[256] = {
                    O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O, O, O, O, O, O, O, O, O,
                    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O,
                    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, A, R, C, R, O, O, G,
                    R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O,
                    O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, O, O, R, R, T, O, R, R,
                    R, R, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
                    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
                    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
                    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
                    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
                    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O};
#undef R
#undef I
#undef O
#undef A
#undef C
#undef G
#undef T
        }
    }

    uint8_t encodeEdge(char c, Direction dir) {
        //std::map<char,uint8_t> shift_table = {{'A',3}, {'T',2}, {'G',1}, {'C',0}};
        uint8_t val{1};
        if (dir == Direction::APPEND)
            return (val << gfa_reader::detail::shift_table[static_cast<uint8_t>(c)]);
        else
            return (val << (gfa_reader::detail::shift_table[static_cast<uint8_t>(c)] + 4));
    }

// Avoiding un-necessary stream creation + replacing strings with string view
// is a bit > than a 2x win!
// implementation from : https://marcoarena.wordpress.com/tag/string_view/
    std::vector<stx::string_view> BinaryGFAReader::split(stx::string_view str,
                                                         char delims) {
        std::vector<stx::string_view> ret;

        stx::string_view::size_type start = 0;
        auto pos = str.find_first_of(delims, start);
        while (pos != stx::string_view::npos) {
            if (pos != start) {
                ret.push_back(str.substr(start, pos - start));
            }
            start = pos + 1;
            pos = str.find_first_of(delims, start);
        }
        if (start < str.length()) {
            ret.push_back(str.substr(start, str.length() - start));
        }
        return ret;
    }

    BinaryGFAReader::BinaryGFAReader(const char *binDir, size_t input_k,
                                     bool buildEqClses, bool buildEdgeVec,
                                     std::shared_ptr<spdlog::logger> logger) {
        logger_ = logger;
        filename_ = std::string(binDir);
        logger_->info("Setting the index/BinaryGfa directory {}", binDir);
        k = input_k;
        buildEdgeVec_ = buildEdgeVec;
        buildEqClses_ = buildEqClses;
        {
            CLI::AutoTimer timer{"Loading contigs", CLI::Timer::Big};
            std::string bfile = filename_ + "/seq.bin";
            seqVec_.deserialize(bfile, false);
        }
        {
            CLI::AutoTimer timer{"Loading contig boundaries", CLI::Timer::Big};
            std::string bfile = filename_ + "/rank.bin";
            rankVec_.deserialize(bfile, false);
        }

    }

    void BinaryGFAReader::encodeSeq(compact::vector<uint64_t, 2> &seqVec, size_t offset, stx::string_view str) {
      for (size_t i = 0; i < str.length(); ++i) {
        auto c = kmers::codeForChar(str[i]);
        seqVec[offset + i] = c;
      }
    }

    compact::vector<uint64_t, 2> &BinaryGFAReader::getContigSeqVec() { return seqVec_; }

    compact::vector<uint64_t, 1> &BinaryGFAReader::getRankVec() { return rankVec_; }

    compact::vector<uint64_t, 8> &BinaryGFAReader::getEdgeVec() { return edgeVec_; }


    void BinaryGFAReader::parseFile() {
        std::string refId;
        uint64_t contigCntr{0}, prevPos{0}, nextPos{1}, ref_cnt{0};
        uint64_t refIdLen;

        k = k + 1;
        CanonicalKmer::k(k);

        // fill out contigId2Seq
        std::unique_ptr<rank9sel> rankSelDict{nullptr};
        rankSelDict.reset(new rank9sel(&rankVec_, rankVec_.size()));
        auto contig2seqSize = rankSelDict->rank(rankVec_.size()-1) + 1;
        std::cerr << contig2seqSize << "\n";
        logger_->info("Done wrapping the rank vector with a rank9sel structure.");
        contigid2seq = pufferfish::util::PackedContigInfoVec(
                                        rankVec_.size(),contig2seqSize);
        while (nextPos < rankVec_.size() and nextPos != 0) {
            nextPos = static_cast<uint64_t>(rankSelDict->select(contigCntr)) + 1;// select(0) is meaningful
            contigid2seq.add(prevPos);
            prevPos = nextPos;
            contigCntr++;
        }
        logger_->info("contig count for validation: {:n}", contigCntr);

        // start and end kmer-hash over the contigs
        // might get deprecated later
        std::ifstream file(filename_ + "/path.bin", std::ios::binary);
        uint64_t refIdSize{0};
        if (file.good()) {
            file.read(reinterpret_cast<char *>(&refIdSize), sizeof(refIdSize));
        }
        uint64_t contigCntPerPath;
        while (file.good()) {
            refIdLen = 0;
            file.read(reinterpret_cast<char *>(&refIdLen), refIdSize);
            if (!file.good()) break;
            refId.assign(refIdLen, '\0');
            file.read(&(refId[0]), refIdLen);
            file.read(reinterpret_cast<char *>(&contigCntPerPath), sizeof(contigCntPerPath));
//            std::cerr << "pathlen: " << contigCntPerPath << "\n";
//            std::cerr << refId << " " << contigCntPerPath << "\n";
            path[ref_cnt].resize(contigCntPerPath);
            for (uint64_t i = 0; i < contigCntPerPath; i++) {
                int64_t contigIdAndOri;
                file.read(reinterpret_cast<char *>(&contigIdAndOri), sizeof(contigIdAndOri));
                if (abs(contigIdAndOri) == 0) {
                    logger_->error("Should never ever happen. ContigId = 0 ");
                    std::exit(3);
                }
                uint64_t contigId = abs(contigIdAndOri)-1;
                if (contigId >= contigid2seq.size()) {
                    logger_->error("Should never ever happen. "
                                   "Found a contigId in a path that was greater than the max contigId: {}, {}",
                                   contigId, contigid2seq.size());
                    std::exit(3);
                }
                bool ori = contigIdAndOri > 0; // ori is set to 1 for case fw where the contigId is positive
                path[ref_cnt][i] = std::make_pair(contigId, ori);
            }
            uint32_t refLength{0};
            bool firstContig{true};
            for (auto &ctig : path[ref_cnt]) {
                uint64_t l = contigid2seq[ctig.first].length - (firstContig ? 0 : (k - 1));
                refLength += l;
                firstContig = false;
            }
            refLengths.push_back(refLength);
            refMap.push_back(refId);
            ref_cnt++;
        }

        //Initialize edgeVec_
        //bad way, have to re-think
        if (buildEdgeVec_) {
            edgeVec_.resize(contigCntr);
            for (uint64_t i = 0; i < contigCntr; i++) edgeVec_[i] = 0;

            for (auto const &ent: path) {
                const std::vector<std::pair<uint64_t, bool>> &contigs = ent.second;

                for (size_t i = 0; i < contigs.size() - 1; i++) {
                    auto cid = contigs[i].first;
                    bool ore = contigs[i].second;
                    size_t forder = cid;//contigid2seq[cid].fileOrder;
                    auto nextcid = contigs[i + 1].first;
                    bool nextore = contigs[i + 1].second;

                    bool nextForder = nextcid;//contigid2seq[nextcid].fileOrder;
                    // a+,b+ end kmer of a , start kmer of b
                    // a+,b- end kmer of a , rc(end kmer of b)
                    // a-,b+ rc(start kmer of a) , start kmer of b
                    // a-,b- rc(start kmer of a) , rc(end kmer of b)
                    //  1. `a+,(*)` we need to append a nucl to the `end-kmer`
                    //  2. `a-,(*)` we need to prepend rc(nucl) to the `start-kmer`
                    //  3. `(*),a+` we need to prepend a nucl to `start-kmer`
                    //  4. `(*),a-` we need to append a rc(nucl) to `end-kmer`

                    CanonicalKmer lastKmerInContig;
                    CanonicalKmer firstKmerInNextContig;
                    Direction contigDirection;
                    Direction nextContigDirection;
                    // If a is in the forward orientation, the last k-mer comes from the end, otherwise it is the reverse complement of the first k-mer
                    if (ore) {
                        auto cinfo = contigid2seq[cid];
                        lastKmerInContig.fromNum(
                                seqVec_.get_int(2 * (cinfo.offset + cinfo.length - k), 2 * k));
                        contigDirection = Direction::APPEND;
                    } else {
                        lastKmerInContig.fromNum(seqVec_.get_int(2 * contigid2seq[cid].offset, 2 * k));
                        lastKmerInContig.swap();
                        contigDirection = Direction::PREPEND;
                    }

                    // If a is in the forward orientation, the first k-mer comes from the beginning, otherwise it is the reverse complement of the last k-mer
                    if (nextore) {
                        firstKmerInNextContig.fromNum(seqVec_.get_int(2 * contigid2seq[nextcid].offset, 2 * k));
                        nextContigDirection = Direction::PREPEND;
                    } else {
                        auto cinfo = contigid2seq[nextcid];
                        firstKmerInNextContig.fromNum(
                                seqVec_.get_int(2 * (cinfo.offset + cinfo.length - k),
                                                2 * k));
                        firstKmerInNextContig.swap();
                        nextContigDirection = Direction::APPEND;
                    }

                    // The character to append / prepend to contig to get to next contig
                    const char contigChar = firstKmerInNextContig.to_str()[k - 1];
                    // The character to prepend / append to next contig to get to contig
                    const char nextContigChar = lastKmerInContig.to_str()[0];

                    edgeVec_[forder] = edgeVec_[forder] | encodeEdge(contigChar, contigDirection);
                    edgeVec_[nextForder] = edgeVec_[nextForder] | encodeEdge(nextContigChar, nextContigDirection);
                }
            }
        }
        k = k - 1;
        logger_->info("Total # of Contigs : {:n}", contigCntr);
        logger_->info("Total # of numerical Contigs : {:n}", contigid2seq.size());
    }


    pufferfish::util::PackedContigInfoVec&
    BinaryGFAReader::getContigNameMap() {
        return contigid2seq;
    }

    std::vector<std::string> &BinaryGFAReader::getRefIDs() { return refMap; }

    std::vector<uint32_t> &BinaryGFAReader::getRefLengths() { return refLengths; }

    void BinaryGFAReader::mapContig2Pos() {
        uint64_t pos = 0;
        uint64_t accumPos;
        uint64_t currContigLength = 0;
        std::vector<uint64_t> cposOffsetvec(contigid2seq.size() + 1, 0);
        uint64_t totalPosCnt = 0;
        for (auto const &ent : path) {
            const std::vector<std::pair<uint64_t, bool>> &contigs = ent.second;
            for (const auto &contig : contigs) {
                cposOffsetvec[contig.first + 1]++;
                totalPosCnt++;
            }
        }
        for (uint64_t i = 0; i < cposOffsetvec.size() - 1; i++) {
            cposOffsetvec[i+1] = cposOffsetvec[i]+cposOffsetvec[i+1];
        }
        logger_->info("Total # of contig vec entries: {:n}", totalPosCnt);
        auto w = static_cast<uint32_t >(std::ceil(std::log2(totalPosCnt+1)));
        logger_->info("bits per offset entry {:n}", w);

        contig2pos.resize(totalPosCnt);
        for (auto const &ent : path) {
            const uint64_t tr = ent.first;
            const std::vector<std::pair<uint64_t, bool>> &contigs = ent.second;
            accumPos = 0;
            for (const auto &contig : contigs) {
                pos = accumPos;
                contig2pos[cposOffsetvec[contig.first]].update(tr, pos, contig.second);
//                std::cerr << cposOffsetvec[contig.first] << ":" << tr << " " << contig2pos[cposOffsetvec[contig.first]].transcript_id() << "\n";
                cposOffsetvec[contig.first]++;
                currContigLength = contigid2seq[contig.first].length;
                //currContigLength = getContigLength(contig.first);
                accumPos += currContigLength - k;
            }
        }
        // We need the +1 here because we store the last entry that is 1 greater than the last offset
        // so we must be able to represent of number of size contigVecSize+1, not contigVecSize.
        cpos_offsets.reset(new compact::vector<uint64_t>(w, contigid2seq.size()+1));
        cpos_offsets->clear_mem();
        (*cpos_offsets)[0] = 0;
        for (uint64_t i = 0; i < cpos_offsets->size()-1; i++) {
            (*cpos_offsets)[i+1] = cposOffsetvec[i];
        }
        cposOffsetvec.clear();
        cposOffsetvec.shrink_to_fit();
        logger_->info("Done constructing the contig vector. {}", cpos_offsets->size());
    }

    void BinaryGFAReader::clearContigTable() {
        refMap.clear(); refMap.shrink_to_fit();
        refLengths.clear(); refLengths.shrink_to_fit();
        contig2pos.clear(); contig2pos.shrink_to_fit();
        cpos_offsets.reset(nullptr);
    }

// Note : We assume that odir is the name of a valid (i.e., existing) directory.
    void BinaryGFAReader::serializeContigTable(const std::string &odir,
                                               const std::vector<std::pair<std::string, uint16_t>>& shortRefsNameLen,
                                               const std::vector<uint32_t>& refIdExtensions) {
        std::string ofile = odir + "/ctable.bin";
        std::string rlfile = odir + "/reflengths.bin";
        std::ofstream ct(ofile);
        std::ofstream rl(rlfile);
        cereal::BinaryOutputArchive ar(ct);
        cereal::BinaryOutputArchive rlAr(rl);
        decltype(refLengths) tmpRefLen;
        tmpRefLen.reserve(refLengths.size() + shortRefsNameLen.size());

        uint64_t shortIdx{0}, longIdx{0}, prevExt{0};
//        for (uint64_t i = 0; i < refLengths.size(); i++) {
        while (longIdx < refLengths.size()) {
            if (refIdExtensions[longIdx] == prevExt) { // If no new short reference
                tmpRefLen.push_back(refLengths[longIdx]); // add the next long reference
                longIdx++;
            } else {
                tmpRefLen.push_back(shortRefsNameLen[shortIdx].second); // add the next short reference
                prevExt++;
                shortIdx++;
            }
        }
        // add all remaining short references
        for (uint64_t i = shortIdx; i < shortRefsNameLen.size(); i++) {
            tmpRefLen.push_back(shortRefsNameLen[i].second);
        }

        {
            // Write out the reference lengths
            rlAr(tmpRefLen);

            // We want to iterate over the contigs in precisely the
            // order they appear in the contig array (i.e., the iterator
            // order of contigid2seq).
            std::vector<std::string> refNames;
            refNames.reserve(refMap.size() + shortRefsNameLen.size());

            shortIdx = 0; longIdx = 0; prevExt = 0;
//            for (uint64_t i = 0; i < refMap.size(); i++) {
            while (longIdx < refLengths.size()) {
                if (refIdExtensions[longIdx] == prevExt) {
                    refNames.push_back(refMap[longIdx]);
                    longIdx++;
                } else {
                    refNames.push_back(shortRefsNameLen[shortIdx].first);
                    prevExt++;
                    shortIdx++;
                }
            }
            for (uint64_t i = shortIdx; i < shortRefsNameLen.size(); i++) {
                refNames.push_back(shortRefsNameLen[i].first);
            }

            ar(refNames);

            ar(refIdExtensions);

            class VecHasher {
            public:
                size_t operator()(const std::vector<uint32_t> &vec) const {
                    return XXH64(const_cast<std::vector<uint32_t> &>(vec).data(),
                                 vec.size() * sizeof(decltype(vec.front())), 0);
                }
            };

            spp::sparse_hash_map<std::vector<uint32_t>, uint32_t, VecHasher> eqMap;
            std::vector<uint32_t> eqIDs;
            //std::vector<std::vector<pufferfish::util::Position>> cpos;

            // Build and store equivalence classes
            if (buildEqClses_) {
                std::string eqfile = odir + "/eqtable.bin";
                std::ofstream et(eqfile);
                cereal::BinaryOutputArchive eqAr(et);
                uint64_t offset = 0;
                for (uint64_t idx = 0; idx < cpos_offsets->size(); idx++) {
                    std::vector<uint32_t> tlist;
                    while (offset < (*cpos_offsets)[idx]) {
                        tlist.push_back(contig2pos[offset].transcript_id());
                        offset++;
                    }
                    std::sort(tlist.begin(), tlist.end());
                    tlist.erase(std::unique(tlist.begin(), tlist.end()), tlist.end());
                    auto eqID = static_cast<uint32_t >(eqMap.size());
                    if (eqMap.contains(tlist)) {
                        eqID = eqMap[tlist];
                    } else {
                        eqMap[tlist] = eqID;
                    }
                    eqIDs.push_back(eqID);
                }
                logger_->info("there were {:n}  equivalence classes", eqMap.size());
                eqAr(eqIDs);
                eqIDs.clear();
                eqIDs.shrink_to_fit();
                std::vector<std::vector<uint32_t>> eqLabels;
                eqLabels.reserve(eqMap.size());
                for (auto &kv : eqMap) {
                    eqLabels.push_back(kv.first);
                }
                std::sort(eqLabels.begin(), eqLabels.end(),
                          [&](const std::vector<uint32_t> &l1,
                              const std::vector<uint32_t> &l2) -> bool {
                              return eqMap[l1] < eqMap[l2];
                          });
                eqAr(eqLabels);
            }
            ar(contig2pos);
            {
              std::string fname = odir + "/" + pufferfish::util::CONTIG_OFFSETS;
              std::ofstream bfile(fname, std::ios::binary);
              cpos_offsets->serialize(bfile);
              bfile.close();
            }
        }
    }

    void BinaryGFAReader::deserializeContigTable() {
        // TODO read the file in the same order as you've written it down.
    }

} // namespace pufferfish
