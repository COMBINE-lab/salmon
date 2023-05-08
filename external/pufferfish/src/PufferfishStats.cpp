#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <bitset>
#include <cereal/archives/binary.hpp>
#include <set>

#include "ProgOpts.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"

#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "PufferfishBinaryGFAReader.hpp"

struct PrefixTree {
    std::string txps;
    uint32_t exactCnt;
    uint32_t prefixCnt;
    uint32_t prefixTotalCnt;
    uint32_t numTxps;

    PrefixTree(const std::string &txps, uint32_t exactCnt, uint32_t prefixCnt, uint32_t prefixTotalCnt)
            : txps(txps), exactCnt(exactCnt), prefixCnt(prefixCnt), prefixTotalCnt(prefixTotalCnt) {
        numTxps = std::count(txps.begin(), txps.end(), '-') + 1;
    }
};

struct NodeInfo {
    uint32_t cid;
    std::set<uint32_t> refid;

    NodeInfo(uint32_t cid) : cid(cid) {}
};

void doFindMotifs(std::string &indexDir) {
    PufferfishIndex pfi(indexDir);
    compact::vector<uint64_t, 1> visited(pfi.numContigs()-1);
    visited.clear_mem();
    std::cerr << "\nTotal contigs: " << visited.size() << "\n";
    uint64_t prevV{0};
    std::vector<char> chars{'A', 'C', 'G', 'T'};
    for (uint64_t v = 0; v < visited.size(); v++) {
        if (not visited[v]) {
            if (v - prevV > 1000) {
                std::cerr << "\r" << v;
                prevV = v;
            }
            visited[v] = 1;
            std::vector<CanonicalKmer> cks;
            cks.push_back(pfi.getStartKmer(v));
            auto end = pfi.getEndKmer(v);
            if (cks[0].getCanonicalWord() != end.getCanonicalWord()) {
                cks.push_back(end);
            }
            auto chits = pfi.getRefPos(cks[0]);
            NodeInfo self(chits.contigIdx_); // should be equal to v
            for (auto &r: chits.refRange) {
                self.refid.insert(r.transcript_id());
            }
            uint32_t intersect{0}, x{0}, y{0}, z{0};
            bool found = false;
            for (auto &ck : cks) {
                std::vector<NodeInfo> oneEndNeighbors;
                uint32_t cntr = 0;
                while (not found and cntr < 2) {
                    for (auto c: chars) {
                        auto cck(ck);
                        if (cntr) {
                            cck.shiftFw(c);
                        } else {
                            cck.shiftBw(c);
                        }
                        chits = pfi.getRefPos(cck);
                        if (not chits.empty() and not visited[chits.contigIdx_]) {
                            oneEndNeighbors.emplace_back(chits.contigIdx_);
                            for (auto &r: chits.refRange) {
                                oneEndNeighbors.back().refid.insert(r.transcript_id());
                            }
                        }

                        if (oneEndNeighbors.size() >= 2) {
                            for (auto r : self.refid) {
                                if (oneEndNeighbors[0].refid.find(r) != oneEndNeighbors[0].refid.end() and
                                    oneEndNeighbors[1].refid.find(r) != oneEndNeighbors[1].refid.end()) {
                                    intersect++;
                                }
                            }
                            x = self.refid.size() - intersect;
                            y = oneEndNeighbors[0].refid.size() - intersect;
                            z = oneEndNeighbors[1].refid.size() - intersect;
                            found = true;
                            visited[oneEndNeighbors[0].cid] = 1;
                            visited[oneEndNeighbors[1].cid] = 1;
                            break;
                        }
                    }
                    cntr++;
                }
            }
            if (found) {
                std::cout << intersect << "\t" << x << "\t" << y << "\t" << z << "\n";
            } else {
                std::cout << self.refid.size() << "\t-1\t-1\t-1\n";
            }
        }
    }
    std::cerr << "\n";
}

void doCtabStats(std::string &indexDir) {
    compact::vector<uint64_t> contigOffsets_{16};
    std::vector<std::string> refNames_;
    std::vector<uint32_t> refExt_;
    std::vector<pufferfish::util::Position> contigTable_;
    std::ifstream contigTableStream(indexDir + "/" + pufferfish::util::CTABLE);
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(refNames_);
    contigTableArchive(refExt_);
    refNames_.clear();
    refNames_.shrink_to_fit();
    refExt_.clear();
    refExt_.shrink_to_fit();

    contigTableArchive(contigTable_);
    contigTableStream.close();
    std::string pfile = indexDir + "/" + pufferfish::util::CONTIG_OFFSETS;
    auto bits_per_element = compact::get_bits_per_element(pfile);
    contigOffsets_.set_m_bits(bits_per_element);
    contigOffsets_.deserialize(pfile, false);
    std::cerr << "contigTable size: " << contigTable_.size()
              << " contigOffsets size: " << contigOffsets_.size() << ", bpe: " << bits_per_element << " , total elements in the contig table: " << contigOffsets_[contigOffsets_.size()-1] << "\n";
    std::vector<std::string> txps;
    for (uint64_t i = 0; i < contigOffsets_.size() - 1; i++) {
        uint32_t idx = contigOffsets_[i];
        auto txp = contigTable_[idx].transcript_id();
        std::string txpstr = std::to_string(txp);
        idx++;
        while (idx < contigOffsets_[i + 1]) {
            if (contigTable_[idx].transcript_id() != txp) {
                txp = contigTable_[idx].transcript_id();
                txpstr += ("-" + std::to_string(txp));
            }
            idx++;
        }
        txps.push_back(txpstr);
    }
    contigTable_.clear();
    contigTable_.shrink_to_fit();
    contigOffsets_.clear();
    std::sort(txps.begin(), txps.end());
    auto prevt = txps[0];
    std::vector<PrefixTree> tree;
    tree.emplace_back(prevt, 0, 0, 0);
    for (auto &t : txps) {
        if (t == prevt) {
            tree.back().exactCnt++;
        } else {
            prevt = t;
            tree.emplace_back(t, 1, 0, 0);
        }
    }
    std::cerr << "original cnt= " << txps.size() << " distinct cnt=" << tree.size() << "\n";
    for (uint64_t i = 0; i < tree.size(); i++) {
        std::string prefix = tree[i].txps;
        uint64_t j = i;
        while (j < tree.size()) {
            bool isPrefix = prefix.size() <= tree[j].txps.size() &&
                            std::mismatch(prefix.begin(), prefix.end(),
                                          tree[j].txps.begin(), tree[j].txps.end()).first == prefix.end();
            if (isPrefix) {
                tree[i].prefixCnt++;
                tree[i].prefixTotalCnt += tree[j].exactCnt;
            } else break;
            j++;
        }
    }
    std::cout << "numTxps exactCnt prefixCnt multiPrefixCnt\n";
    for (auto &t:tree) {
        std::cout << t.numTxps << " "
                  << t.exactCnt << " " << t.prefixCnt << " " << t.prefixTotalCnt << "\n";
    }
};

int pufferfishStats(pufferfish::StatsOptions &statsOpts) {
    auto indexDir = statsOpts.indexDir;
    switch (statsOpts.statType) {
        case pufferfish::StatType::ctab:
            doCtabStats(indexDir);
            break;
        case pufferfish::StatType::motif:
            doFindMotifs(indexDir);
            break;
    }
    return 0;
}


