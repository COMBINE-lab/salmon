//
// Created by Fatemeh Almodaresi on 5/17/18.
//

#ifndef PUFFOUTREADER_H
#define PUFFOUTREADER_H

#include "spdlog/spdlog.h"
#include <fstream>
#include <iostream>

#include "PuffoutFilePointer.h"
#include "AlignmentGroup.hpp"
#include "RapMapUtils.hpp"

using QuasiAlignment = rapmap::utils::QuasiAlignment;
using AlnGroupVector = std::vector<AlignmentGroup<QuasiAlignment>>;
using MateStatus = rapmap::utils::MateStatus;


class PuffoutParser {
public:
    PuffoutParser() {
        logger = spdlog::get("jointLog");
    }

    /*
     * will prepare at max as many alignments as miniBatch requires
     */
    bool nextChunkOfAlignments(PuffoutFilePointer *puffoutFilePointer,
                               AlnGroupVector &alignments,
                               std::mutex &iomutex) {
        using rapmap::utils::MateStatus;
        readCnt = 0;
        while (readCnt < alignments.size()) {
            // parse read chunks
            if (!chunk.hasNext()) {// try to read a new chunk from file
                if (!puffoutFilePointer->readChunk(chunk, iomutex)) {// if nothing left to process return false
                    if (readCnt == 0) // and no alignments from previous chunk return false
                        return false;
                    else
                        return true; // if the alignment vector is not empty return true
                }
            }

            while (chunk.hasNext() && readCnt < alignments.size()) {
                uint32_t mcnt{0};
                rLenType llen{0}, rlen{0};
                // load all the alignments for the read
                // if we've reached max available space for alignments
                /* if (alignments.size() == readCnt) {
                     AlignmentGroup<QuasiAlignment> alng;
                     alignments.push_back(alng);
                 }*/
                auto &jointHitGroup = alignments[readCnt];
                jointHitGroup.clearAlignments();
                // read name
                uint8_t readNameLength{0};
                std::string readName;
                chunk.fill(readNameLength);
                //std::cerr << "rl " << (uint64_t)readNameLength << "\n";
                chunk.fill(readName, readNameLength);
                // mapping count
                chunk.fill(mcnt);
                chunk.fill(llen); // left read len
                if (puffoutFilePointer->isMappingPaired()) {
                    chunk.fill(rlen); // right read len
                }
                if (mcnt == 0) std::cerr << readName << "\n";
                //std::cerr << "ra" << readCnt << " " << readName << " " << mcnt << "\n";

                //logger->info("read name {}", readName);
                //logger->info("count {}", mcnt);
                double maxScore{0};
                for (uint64_t mappingCntr = 0; mappingCntr < mcnt; mappingCntr++) {
                    uint32_t puff_id{0}, fraglen{0};
                    rLenType lcnt{0}, rcnt{0};
                    refLenType reflPos{0}, refrPos{0};
                    double reflScore{0}, refrScore{0};
                    bool lori, rori;
                    chunk.fill(puff_id);
                    chunk.fill(lcnt);
                    if (puffoutFilePointer->isMappingPaired()) {
                        chunk.fill(rcnt);
                    }
                    if (lcnt > 0) {
                        chunk.fill(reflScore);
                        chunk.fill(reflPos);
                        lori = reflPos & PuffoutParser::HighBitMask;
                        reflPos = reflPos & PuffoutParser::LowBitsMask;
                    }
                    if (rcnt > 0) {
                        chunk.fill(refrScore);
                        chunk.fill(refrPos);
                        rori = refrPos & PuffoutParser::HighBitMask;
                        refrPos = refrPos & PuffoutParser::LowBitsMask;
                    }
                    if (lcnt > 0 && rcnt > 0) {
                        fraglen = reflPos < refrPos ? refrPos + rlen - reflPos : reflPos + llen - refrPos;
                        QuasiAlignment qaln(puff_id, reflPos, lori, llen, fraglen, true);

                        qaln.mateLen = rlen;
                        qaln.matePos = refrPos;
                        qaln.mateIsFwd = rori;
                        qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;

                        jointHitGroup.addAlignment(qaln);

                    } else if (lcnt > 0) {
                        QuasiAlignment qaln(puff_id, reflPos, lori, llen, 0, false);
                        if (puffoutFilePointer->isMappingPaired())
                            qaln.mateStatus = MateStatus::PAIRED_END_LEFT;
                        else
                            qaln.mateStatus = MateStatus::SINGLE_END;
                        jointHitGroup.addAlignment(qaln);
                    } else if (rcnt > 0) {
                        QuasiAlignment qaln(puff_id, refrPos, rori, rlen, 0, false);
                        qaln.mateStatus = MateStatus::PAIRED_END_RIGHT;
                        jointHitGroup.addAlignment(qaln);
                    }
                    auto &h = jointHitGroup.alignments().back();
                    h.score(reflScore + refrScore);
                    if (h.score() > maxScore) maxScore = h.score();
                    switch (h.mateStatus) {
                        case MateStatus::PAIRED_END_LEFT: {
                            h.format = salmon::utils::hitType(h.pos, h.fwd);
                        }
                            break;
                        case MateStatus::PAIRED_END_RIGHT: {
                            // we pass in !h.fwd here because the right read
                            // will have the opposite orientation from its mate.
                            // NOTE : We will try recording what the mapped fragment
                            // actually is, not to infer what it's mate should be.
                            h.format = salmon::utils::hitType(h.pos, h.fwd);
                        }
                            break;
                        case MateStatus::PAIRED_END_PAIRED: {
                            uint32_t end1Pos = (h.fwd) ? h.pos : h.pos + h.readLen;
                            uint32_t end2Pos =
                                    (h.mateIsFwd) ? h.matePos : h.matePos + h.mateLen;
                            bool canDovetail = false;
                            h.format =
                                    salmon::utils::hitType(end1Pos, h.fwd, h.readLen, end2Pos,
                                                           h.mateIsFwd, h.mateLen, canDovetail);
                        }
                            break;
                        case MateStatus::SINGLE_END: {
                            // do nothing
                        }
                            break;
                        case MateStatus::NOTHING: {
                            std::cerr << "WARNING: shouldn't get a mateStatus of NOTHING\n";
                        }
                        break;
                    }
                    /*switch(jointHitGroup.alignments().back().mateStatus) {
                        case MateStatus::PAIRED_END_PAIRED:
                            std::cerr << "paired\n"; break;
                        case MateStatus::PAIRED_END_LEFT:
                            std::cerr << "left\n"; break;
                        case MateStatus::PAIRED_END_RIGHT:
                            std::cerr << "right\n"; break;
                        case MateStatus::SINGLE_END:
                            std::cerr << "single\n"; break;
                    }*/
                }
                for (auto &h : jointHitGroup.alignments()) {
                    h.score(std::exp(-(maxScore - h.score())));
                }
                // sort alignments based on their TranscriptIDs
                std::sort(jointHitGroup.alignments().begin(), jointHitGroup.alignments().end(),
                          [](const QuasiAlignment & x, const QuasiAlignment & y) -> bool {
                              return x.transcriptID() < y.transcriptID();
                          });

                readCnt++;
            }
        }
        return true;
    }

    uint64_t getReadCnt() { return readCnt; }

private:

    Chunk chunk;
    uint64_t readCnt{0};
    std::shared_ptr <spdlog::logger> logger;
    static constexpr const refLenType HighBitMask = 1u << (sizeof(refLenType) * 8 - 1);
    static constexpr const refLenType LowBitsMask = HighBitMask - 1;
    std::mutex iomutex;

};

#endif //PUFFOUTREADER_H
