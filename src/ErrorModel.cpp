#include <iostream>
#include <cstdio>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY

#include "ErrorModel.hpp"
#include "Transcript.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "UnpairedRead.hpp"
#include "ReadPair.hpp"

ErrorModel::ErrorModel(double alpha, uint32_t maxExpectedReadLen) :
    maxExpectedLen_(maxExpectedReadLen),
    mismatchLeft_(maxExpectedReadLen, AtomicMatrix<double>(16, 4, alpha)),
    mismatchRight_(maxExpectedReadLen, AtomicMatrix<double>(16, 4, alpha)),
    isEnabled_(true),
    maxLen_(0),
    burnedIn_(false) {}

bool ErrorModel::burnedIn() { return burnedIn_; }
void ErrorModel::burnedIn(bool burnedIn) { burnedIn_ = burnedIn; }

double ErrorModel::logLikelihood(bam_seq_t* read, Transcript& ref,
                                 std::vector<AtomicMatrix<double>>& mismatchProfile){
    using namespace salmon::stringtools;
    bool useQual{false};
    size_t readIdx{0};
    auto transcriptIdx = bam_pos(read);
    size_t transcriptLen = ref.RefLength;
    // if the read starts before the beginning of the transcript,
    // only consider the part overlapping the transcript
    if (transcriptIdx < 0) {
        readIdx = -transcriptIdx;
        transcriptIdx = 0;
    }

    // unsigned version of transcriptIdx
    size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

    if (uTranscriptIdx >= transcriptLen) {
        std::lock_guard<std::mutex> l(outputMutex_);
        std::cerr << "transcript index = " << uTranscriptIdx << ", transcript length = " << transcriptLen << "\n";
        return salmon::math::LOG_0;
    }

    //salmon::stringtools::strand readStrand = (BAM_FREVERSE & read->core.flag) ? salmon::stringtools::strand::reverse :
    //                                            salmon::stringtools::strand::forward;
    salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;
    double logLike = salmon::math::LOG_1;

    uint8_t* qseq = reinterpret_cast<uint8_t*>(bam_seq(read));//bam1_seq(read);
    auto qualStr = reinterpret_cast<uint8_t*>(bam_qual(read));//bam1_qual(read);
    size_t readLen = static_cast<size_t>(bam_seq_len(read));
    size_t basesChecked = 0;
    size_t numInc = std::min(transcriptLen - uTranscriptIdx, readLen - readIdx);
    while (basesChecked < numInc) {
        size_t curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
        size_t prevReadBase = (readIdx > 0) ? samToTwoBit[bam_seqi(qseq, readIdx-1)] : 0;
        size_t refBase = samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)];
        size_t index = prevReadBase + refBase;
        int qval = qualStr[readIdx];
        double qual = (useQual) ? salmon::stringtools::phredToLogProb[qval] : salmon::math::LOG_1;
        logLike += mismatchProfile[readIdx](index, curReadBase) + qual;
        ++readIdx;
        ++uTranscriptIdx;
        ++basesChecked;
    }
    return logLike;
}

double ErrorModel::logLikelihood(const ReadPair& hit, Transcript& ref){
    double logLike = salmon::math::LOG_1;
    if (BOOST_UNLIKELY(!isEnabled_)) { return logLike; }

    if (!hit.isPaired()) {
        if (hit.isLeftOrphan()) {
            return logLikelihood(hit.read1, ref, mismatchLeft_);
        } else {
            return logLikelihood(hit.read1, ref, mismatchRight_);
        }
    }

    bam_seq_t* leftRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read1 : hit.read2;
    bam_seq_t* rightRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read2 : hit.read1;

    size_t leftLen = static_cast<size_t>(bam_seq_len(leftRead));
    size_t rightLen = static_cast<size_t>(bam_seq_len(rightRead));

    // NOTE: Raise a warning in this case?
    if (BOOST_UNLIKELY((leftLen > maxExpectedLen_) or
                       (rightLen > maxExpectedLen_))) {
        return logLike;
    }

    if (leftRead) {
        logLike += logLikelihood(leftRead, ref, mismatchLeft_);
    }

    if (rightRead) {
        logLike += logLikelihood(rightRead, ref, mismatchRight_);
    }
    if (logLike == salmon::math::LOG_0) {
            std::lock_guard<std::mutex> lock(outputMutex_);
            std::cerr << "orphan status: " << hit.orphanStatus << "\n";
            std::cerr << "error likelihood: " << logLike << "\n";
    }

    return logLike;
}

double ErrorModel::logLikelihood(const UnpairedRead& hit, Transcript& ref){
    double logLike = salmon::math::LOG_1;
    if (BOOST_UNLIKELY(!isEnabled_)) { return logLike; }

    bam_seq_t* read = hit.read;
    size_t readLen = static_cast<size_t>(bam_seq_len(read));
    // NOTE: Raise a warning in this case?
    if (BOOST_UNLIKELY(readLen > maxExpectedLen_)) {
        return logLike;
    }
    logLike += logLikelihood(read, ref, mismatchLeft_);

    if (logLike == salmon::math::LOG_0) {
            std::lock_guard<std::mutex> lock(outputMutex_);
            std::cerr << "error log likelihood: " << logLike << "\n";
    }

    return logLike;
}

void ErrorModel::update(const UnpairedRead& hit, Transcript& ref, double p, double mass){
    if (mass == salmon::math::LOG_0) { return; }
    if (BOOST_UNLIKELY(!isEnabled_)) { return; }
    bam_seq_t* leftRead = hit.read;
    update(leftRead, ref, p, mass, mismatchLeft_);
}

void ErrorModel::update(bam_seq_t* read, Transcript& ref, double p, double mass,
                        std::vector<AtomicMatrix<double>>& mismatchProfile) {
    using namespace salmon::stringtools;
    bool useQual{false};
    size_t readIdx{0};
    auto transcriptIdx = bam_pos(read);
    size_t transcriptLen = ref.RefLength;
    // if the read starts before the beginning of the transcript,
    // only consider the part overlapping the transcript
    if (transcriptIdx < 0) {
        readIdx = -transcriptIdx;
        transcriptIdx = 0;
    }
    // unsigned version of transcriptIdx
    size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

    // Only attempt to update the model if the read overlaps the transcript.
    if (uTranscriptIdx < transcriptLen) {
        salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;

        uint32_t numMismatch{0};
        uint8_t* qseq = reinterpret_cast<uint8_t*>(bam_seq(read));
        uint8_t* qualStr = reinterpret_cast<uint8_t*>(bam_qual(read));
        size_t readLen = static_cast<size_t>(bam_seq_len(read));
        size_t basesChecked = 0;
        size_t numInc = std::min(transcriptLen - uTranscriptIdx, readLen - readIdx);
        while (basesChecked < numInc) {
            size_t curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
            size_t prevReadBase = (readIdx > 0) ? samToTwoBit[bam_seqi(qseq, readIdx-1)] : 0;
            size_t refBase = samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)];
            size_t index = prevReadBase + refBase;
            if (curReadBase != refBase) { ++numMismatch; }
            int qval = qualStr[readIdx];
            double qual = (useQual) ? salmon::stringtools::phredToLogProb[qval] : salmon::math::LOG_1;
            mismatchProfile[readIdx].increment(index, curReadBase, mass+p+qual);
            ++readIdx;
            ++uTranscriptIdx;
            ++basesChecked;
        }

        /*
        // DEBUG: print number of mismatches
        if (numMismatch > 3) {
            std::lock_guard<std::mutex> lock(outputMutex_);
            std::cerr << "read of length " << readLen << " had " << numMismatch
                      << " mismatches; dir is " << ((readStrand == strand::forward) ? "forward" : "reverse") << "\n";
           int32_t transcriptStartIdx = bam_pos(read);
           size_t rstart = (transcriptStartIdx < 0) ? -transcriptStartIdx : 0;
           if (transcriptStartIdx < 0) { transcriptStartIdx = 0; }
           if (op == BAM_CSOFT_CLIP) {
               rstart += opLen;
           }
           size_t numInc = std::min(transcriptLen - transcriptStartIdx, readLen - rstart);
           int32_t rinc{1};
            std::cerr << "read name is " << bam_name(read) << "\n";
            std::cerr << "readStart = " << rstart << "\n";
            std::cerr << "txp start = " << transcriptStartIdx << " (txp len = " << ref.RefLength << "), (txp name = " << ref.RefName << ")\n";
            std::cerr << "read is: ";
            int32_t rp = rstart;
            for (int32_t ri = 0; ri < numInc; ri++, rp++) {
                std::cerr << samCodeToChar[bam_seqi(qseq, rp)];
            } std::cerr << "\n";
            std::cerr << "ref is:  ";

            int32_t tp = transcriptStartIdx;
            for (int32_t ri = 0; ri < numInc; ri++, tp+=rinc) {
                std::cerr << samCodeToChar[ref.baseAt(tp, readStrand)];
            } std::cerr << "\n";
        }
        */

        maxLen_ = std::max(maxLen_, readLen);
        if (BOOST_UNLIKELY(maxLen_ > maxExpectedLen_)) {
            std::lock_guard<std::mutex> lock(outputMutex_);
            std::cerr << "Encountered read longer than maximum expected length of "
                << maxExpectedLen_ << ", not applying error model\n";
            isEnabled_ = false;
        }
    }
}

void ErrorModel::update(const ReadPair& hit, Transcript& ref, double p, double mass){
    if (mass == salmon::math::LOG_0) { return; }
    if (BOOST_UNLIKELY(!isEnabled_)) { return; }

    bam_seq_t* leftRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read1 : hit.read2;
    bam_seq_t* rightRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read2 : hit.read1;

    if (leftRead) {
        update(leftRead, ref, p, mass, mismatchLeft_);
    }

    if (rightRead) {
        update(rightRead, ref, p, mass, mismatchRight_);
    }
}
