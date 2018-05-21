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
    PuffoutParser() {}

    /*
     * will prepare at max as many alignments as miniBatch requires
     */
    bool nextChunkOfAlignments(PuffoutFilePointer &puffoutFilePointer, AlnGroupVector &alignments) {

        // parse read chunks
        readCnt = 0;
        while (chunk.hasNext() || readCnt < alignments.size()) {
            if (!chunk.hasNext()) // if nothing has left from the last chunk fetched from file
                if (!puffoutFilePointer.readChunk(chunk) && readCnt == 0) // try to read a new chunk from file
                    return false; // if nothing left to process return false
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
            uint8_t readNameLength;
            std::string readName;
            chunk.fill(readNameLength);
            chunk.fill(readName, readNameLength);

            // mapping count
            chunk.fill(mcnt);
            chunk.fill(llen); // left read len
            if (puffoutFilePointer.isMappingPaired()) {
                chunk.fill(rlen); // right read len
            }

            for (uint64_t mappingCntr = 0; mappingCntr < mcnt; mappingCntr++) {
                uint32_t puff_id{0}, fraglen{0};
                rLenType lcnt{0}, rcnt{0};
                refLenType reflPos{0}, refrPos{0};
                bool lori, rori;
                chunk.fill(puff_id);
                chunk.fill(lcnt);
                if (puffoutFilePointer.isMappingPaired()) {
                    chunk.fill(rcnt);
                }
                if (lcnt > 0) {
                    chunk.fill(reflPos);
                    lori = reflPos & PuffoutParser::HighBitMask;
                    reflPos = reflPos & PuffoutParser::LowBitsMask;
                }
                if (rcnt > 0) {
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
                    if (puffoutFilePointer.isMappingPaired())
                        qaln.mateStatus = MateStatus::PAIRED_END_LEFT;
                    else
                        qaln.mateStatus = MateStatus::SINGLE_END;
                    jointHitGroup.addAlignment(qaln);
                } else if (rcnt > 0) {
                    QuasiAlignment qaln(puff_id, refrPos, rori, rlen, 0, false);
                    qaln.mateStatus = MateStatus::PAIRED_END_RIGHT;
                    jointHitGroup.addAlignment(qaln);
                }
            }

            readCnt++;
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

};

#endif PUFFOUTREADER_H
