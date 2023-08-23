#include <cstdio>
#include <iostream>
#include <map>
#include <tuple>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/geometric.hpp>

#include "dbg.hpp"
#include "ONTAlignmentModel.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"



/**
 * GetStateIndex to get the index of the state for a kmer based on the number of indels and mismatches it contains, 
 * as well as whether it is a homopolymer (all mismatches and one nucleotide)
 * Some states will be impossible (ex: k-1 indels and k-1 mismatches would be impossible in a kmer), 
 * so their probabilities will remain at zero, and we can ignore them, but we leave them in the state ordering for convenience/clarity.
*/
static const uint32_t numMismatchStates = ONTAlignmentModel::kmerLength + 1;
  constexpr static uint32_t getStateIndexNoHomopolymer(size_t numIndels, size_t numMisMatch) {
      return numMismatchStates * numIndels + numMisMatch;
  }

uint32_t numStatesNoHomo = getStateIndexNoHomopolymer(ONTAlignmentModel::kmerLength, 0)+1; 

uint32_t getStateIndex(size_t numIndels, size_t numMisMatch, bool homopolymer ) {
      if (homopolymer){
        return getStateIndexNoHomopolymer(numIndels,numMisMatch) + numStatesNoHomo;
      }
      return getStateIndexNoHomopolymer(numIndels,numMisMatch);
  }
//The last two states in order are: the start state and the end state 
uint32_t startStateIdx = getStateIndex(ONTAlignmentModel::kmerLength, 0, true) + 1;
uint32_t endStateIdx = getStateIndex(ONTAlignmentModel::kmerLength, 0, true) + 2;
uint32_t numStatesTotal = getStateIndex(ONTAlignmentModel::kmerLength, 0, true) + 3;



ONTAlignmentModel::ONTAlignmentModel(double alpha, uint32_t readBins)
    : isEnabled_(true)
    , readBins_(readBins)
    , printed(false)
    , transitionProbs_(readBins)
{ 
  for (size_t i = 0; i < readBins; ++i) {
    transitionProbs_[i] = std::move(AtomicMatrix<double>(
        numStatesTotal, numStatesTotal, alpha));
  }

}

//The Column struct represents one column of the alignment
struct Column
{uint8_t ref_base;
uint8_t read_base;
};

//is_homopolymer determines if a set of columns contains only one nucleotide in either column + indels
bool is_homopolymer(Column chunk[]){
  size_t counts[4] = {0};
  size_t numGaps = 0;
  for (size_t i = 0; i < ONTAlignmentModel::kmerLength; i++) {
    size_t symbol = chunk[i].ref_base;
    if (symbol <= 3){
     counts[symbol]+=1;
    }
    else {
      numGaps += 1;
    }
  }
  for (size_t i = 0; i < 4; i++){
    if (counts[i] >= ONTAlignmentModel::kmerLength - 1){
      return true; 
    }
  }
  return false;
}


//processChunk both returns the state in the markov model corresponding to a kmer of the alignment, 
//and updates the chunk so the section of the kmer that overlaps with the next kmer is filled in.
uint32_t processChunk(Column chunk[]){
  size_t numIndels = 0;
  size_t numMismatch = 0;

  for  (size_t i = 0; i < ONTAlignmentModel::kmerLength; i++){
    if (chunk[i].ref_base >3 ||chunk[i].read_base>3 ){
      numIndels ++; 
    }
    else if(chunk[i].ref_base != chunk[i].read_base){
      numMismatch++;
    }
  }
  uint32_t stateIdx = getStateIndex(numIndels, numMismatch, is_homopolymer(chunk));

  for (size_t i = 0; i < ONTAlignmentModel::kmerLength-ONTAlignmentModel::stepSize; i++){
   chunk[i].ref_base = chunk[i+ONTAlignmentModel::stepSize].ref_base;
   chunk[i].read_base = chunk[i+ONTAlignmentModel::stepSize].read_base;
  }
  return stateIdx;
}
//Testing functions (numIndels and numMismatch) useful to check the correct stateIndex is being calculated
//The bases 0-3 are nucleotides, anything else is considered to be an indel consistent with the definition of AlignmentModelChar in AlignmentCommon.hpp
uint32_t numIndels(Column chunk[]){
  size_t numIndels = 0;
  size_t numMismatch = 0;

  for  (size_t i = 0; i < ONTAlignmentModel::kmerLength; i++){
    if (chunk[i].ref_base >3 ||chunk[i].read_base>3 ){
      numIndels ++; 
    }
    else if(chunk[i].ref_base != chunk[i].read_base){
      numMismatch++;
    }
  }
  return numIndels;
}
uint32_t numMismatch(Column chunk[]){
  size_t numIndels = 0;
  size_t numMismatch = 0;

  for  (size_t i = 0; i < ONTAlignmentModel::kmerLength; i++){
    if (chunk[i].ref_base >3 ||chunk[i].read_base>3 ){
      numIndels ++; 
    }
    else if(chunk[i].ref_base != chunk[i].read_base){
      numMismatch++;
    }
  }
  return numMismatch;
}





ONTAlignmentModel::AlnModelProb ONTAlignmentModel::logLikelihood( bam_seq_t* read, bam_seq_t* primary, Transcript& ref,
std::vector<AtomicMatrix<double>>& transitionProbs){
  //Implement mm, returns a pair {likelihood, background likelihood}
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
    std::cerr << "transcript index = " << uTranscriptIdx
              << ", transcript length = " << transcriptLen << "\n";
    return {salmon::math::LOG_0, salmon::math::LOG_1};
  }


  uint32_t* cigar = bam_cigar(read);
  uint32_t cigarLen = bam_cigar_len(read);
  
  const bool usePrimary = bam_seq_len(read) == 0 && primary != nullptr;
  uint8_t* qseq = reinterpret_cast<uint8_t*>(!usePrimary ? bam_seq(read) : bam_seq(primary));
  uint8_t* qualStr = reinterpret_cast<uint8_t*>(!usePrimary ? bam_qual(read) : bam_qual(primary));
  int32_t readLen = !usePrimary ? bam_seq_len(read) : bam_seq_len(primary);

  if (cigarLen == 0 or !cigar) {
    return {salmon::math::LOG_EPSILON, salmon::math::LOG_1};
  }

  salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;
  // foreground likelihood
  double logLike = salmon::math::LOG_1;
  // background likelihood
  double bgLogLike = salmon::math::LOG_1;

  bool advanceInRead{false};
  bool advanceInReference{false};
  uint32_t readPosBin{0};
  uint32_t cigarIdx{0};

  uint32_t prevStateIdx{startStateIdx};
  uint32_t curStateIdx{0};
  double invLen = static_cast<double>(readBins_) / bam_seq_len(read);
  size_t chunk_idx = 0;
  size_t overlap = kmerLength - stepSize;
  Column chunk[kmerLength];
  for (uint32_t cigarIdx = 0; cigarIdx < cigarLen; ++cigarIdx) {

    uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
    enum cigar_op op =
        static_cast<enum cigar_op>(cigar[cigarIdx] & BAM_CIGAR_MASK);
    size_t curReadBase = (BAM_CONSUME_SEQ(op)) ? samToTwoBit[bam_seqi(qseq, readIdx)] : 0;
    size_t curRefBase = (BAM_CONSUME_REF(op)) ? samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)] : 0;
    advanceInRead = false;
    advanceInReference = false;

    for (size_t i = 0; i < opLen; ++i) {
      if (advanceInRead) {
        // Shouldn't happen!
        if (readIdx >= static_cast<decltype(readIdx)>(readLen)) {
          if (logger_) {
            logger_->warn("(in logLikelihood()) CIGAR string for read [{}] "
                          "seems inconsistent. It refers to non-existant "
                          "positions in the read! The use_primary flag is {}.",
                          bam_name(read), usePrimary);
            std::stringstream cigarStream;
            for (size_t j = 0; j < cigarLen; ++j) {
              uint32_t opLen = cigar[j] >> BAM_CIGAR_SHIFT;
              enum cigar_op op =
                  static_cast<enum cigar_op>(cigar[j] & BAM_CIGAR_MASK);
              cigarStream << opLen << opToChr(op);
            }
          }
          return {logLike, bgLogLike};
        }
        curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
        readPosBin = static_cast<uint32_t>((readIdx * invLen));
        advanceInRead = false;
      }
      if (advanceInReference) {
        // Shouldn't happen!
        if (uTranscriptIdx >= transcriptLen) {
          if (logger_) {
            logger_->warn(
                "(in logLikelihood()) CIGAR string for read [{}] "
                "seems inconsistent. It refers to non-existant "
                "positions in the reference! Transcript name "
                "is {}, length is {}, id is {}. Read things refid is {}",
                bam_name(read), ref.RefName, transcriptLen, ref.id,
                bam_ref(read));
          }
          return {logLike, bgLogLike};
        }
        curRefBase = samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)];
        advanceInReference = false;
      }

      setBasesFromCIGAROp_(
          op, curRefBase, curReadBase); //, readStream, matchStream, refStream);
      
      chunk[chunk_idx].ref_base = curRefBase;
      chunk[chunk_idx].read_base = curReadBase;
      
      chunk_idx ++;
      //When we get to the end of our kmer, get the state index in the markov model  
      if (chunk_idx >= kmerLength)
        {
          //Update the index in the next kmer
          chunk_idx = kmerLength-stepSize;
          //Shift overlapping section into first section of kmer, and get state index
          curStateIdx = processChunk(chunk);
        }
      double tp = transitionProbs[readPosBin](prevStateIdx, curStateIdx);
      logLike += tp;
      bgLogLike += transitionProbs[readPosBin](0, 0);
      prevStateIdx = curStateIdx;
      if (BAM_CONSUME_SEQ(op)) {
        ++readIdx;
        advanceInRead = true;
      }
      if (BAM_CONSUME_REF(op)) {
        ++uTranscriptIdx;
        advanceInReference = true;
      }
    }
  }
  //Add transition probability for ending the transcript after this state.
  double tp = transitionProbs[readPosBin](curStateIdx, endStateIdx);
  logLike += tp;
  bgLogLike += transitionProbs[readPosBin](0, 0);

  return {logLike, bgLogLike};
}


double ONTAlignmentModel::logLikelihood(const UnpairedRead& hit, const UnpairedRead& primary, Transcript& ref) {
  using salmon::math::LOG_0;
  using salmon::math::LOG_1;
  constexpr auto dmin = std::numeric_limits<double>::min();
  constexpr double llMin = 1e-10; // Ignore alignment with likelihood below that number

  // If chimeric alignment, doesn't align to this transcript. Return 0
  // probability.
  if(bam_aux_find(hit.read, "SA"))
    return LOG_0;

  double errorllh = LOG_1; // Error and clip Log Likelihood (front and back)

  //Use transition probs to determine error log likelihood
  double logLike = salmon::math::LOG_1;
  double bg = salmon::math::LOG_1;
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return logLike;
  }
  bam_seq_t* read = hit.read;
  size_t readLen = static_cast<size_t>(bam_seq_len(read));
  auto alnLogProb = logLikelihood(read, primary.read, ref, transitionProbs_);
  logLike += alnLogProb.fg;
  bg += alnLogProb.bg;

  if (logLike == salmon::math::LOG_0) {
    std::lock_guard<std::mutex> lock(outputMutex_);
    std::cerr << "error log likelihood: " << logLike << "\n";
  }

  errorllh = logLike - bg;

  return errorllh;

}

void ONTAlignmentModel::update(
    bam_seq_t* read, bam_seq_t* primary,
    Transcript& ref, double p, double mass,
    std::vector<AtomicMatrix<double>>& transitionProbs) {
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
    std::cerr << "transcript index = " << uTranscriptIdx
              << ", transcript length = " << transcriptLen << "\n";
    return;
  }

  // std::stringstream readStream, matchStream, refStream;

  uint32_t* cigar = bam_cigar(read);
  uint32_t cigarLen = bam_cigar_len(read);
  const bool usePrimary = bam_seq_len(read) == 0 && primary != nullptr;
  uint8_t* qseq = reinterpret_cast<uint8_t*>(!usePrimary ? bam_seq(read) : bam_seq(primary));
  uint8_t* qualStr = reinterpret_cast<uint8_t*>(!usePrimary ? bam_qual(read) : bam_qual(primary));
  int32_t readLen = !usePrimary ? bam_seq_len(read) : bam_seq_len(primary);
  
  if (cigarLen == 0 or !cigar) {
    return ;
  }

  salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;
  bool advanceInRead{false};
  bool advanceInReference{false};
  uint32_t readPosBin{0};
  uint32_t cigarIdx{0};

  uint32_t prevStateIdx{startStateIdx};
  uint32_t curStateIdx{0};
  double invLen = static_cast<double>(readBins_) / bam_seq_len(read);
  size_t chunk_idx = 0;
  size_t chunk_size = kmerLength;
  size_t overlap = kmerLength - stepSize;
  Column chunk[chunk_size];

  for (uint32_t cigarIdx = 0; cigarIdx < cigarLen; ++cigarIdx) {
    uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
    enum cigar_op op =
        static_cast<enum cigar_op>(cigar[cigarIdx] & BAM_CIGAR_MASK);
    size_t curReadBase = (BAM_CONSUME_SEQ(op)) ? samToTwoBit[bam_seqi(qseq, readIdx)] : 0;
    size_t curRefBase = (BAM_CONSUME_REF(op)) ? samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)] : 0;
    advanceInRead = false;
    advanceInReference = false;

    for (size_t i = 0; i < opLen; ++i) {
      if (advanceInRead) {
        // Shouldn't happen!
        if (readIdx >= static_cast<decltype(readIdx)>(readLen)) {
          if (logger_) {
            logger_->warn("(in update()) CIGAR string for read [{}] "
                          "seems inconsistent. It refers to non-existant "
                          "positions in the read! Use primary flag is {}",
                          bam_name(read), usePrimary);
            std::stringstream cigarStream;
            for (size_t j = 0; j < cigarLen; ++j) {
              uint32_t opLen = cigar[j] >> BAM_CIGAR_SHIFT;
              enum cigar_op op =
                  static_cast<enum cigar_op>(cigar[j] & BAM_CIGAR_MASK);
              cigarStream << opLen << opToChr(op);
            }
          }
          return;
        }
        curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
        readPosBin = static_cast<uint32_t>((readIdx * invLen));
        advanceInRead = false;
      }
      if (advanceInReference) {
        // Shouldn't happen!
        if (uTranscriptIdx >= transcriptLen) {
          if (logger_) {
            logger_->warn(
                "(in update()) CIGAR string for read [{}] "
                "seems inconsistent. It refers to non-existant "
                "positions in the reference! Transcript name "
                "is {}, length is {}, id is {}. Read things refid is {}",
                bam_name(read), ref.RefName, transcriptLen, ref.id,
                bam_ref(read));
          }
          return;
        }
        curRefBase = samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)];
        advanceInReference = false;
      }

      setBasesFromCIGAROp_(
          op, curRefBase, curReadBase); //, readStream, matchStream, refStream);
      

      chunk[chunk_idx].ref_base = curRefBase;
      chunk[chunk_idx].read_base = curReadBase;
     
      chunk_idx ++;
      if (chunk_idx >= kmerLength)
        {             
          curStateIdx = processChunk(chunk);
          chunk_idx = kmerLength-stepSize;
        }
      
      transitionProbs[readPosBin].increment(prevStateIdx, curStateIdx,mass + p);
      prevStateIdx = curStateIdx;
      if (BAM_CONSUME_SEQ(op)) {
        ++readIdx;
        advanceInRead = true;
      }
      if (BAM_CONSUME_REF(op)) {
        ++uTranscriptIdx;
        advanceInReference = true;
      }
    }
    //Add the transition probability for ending
    curStateIdx = endStateIdx;
    transitionProbs[readPosBin].increment(prevStateIdx, curStateIdx,mass + p);
  }
  return ;
  }

// Update the Markov model. The reads are binned based on their
// length (the length after soft clipping).
void ONTAlignmentModel::update(const UnpairedRead& hit, const UnpairedRead& primary,
                               Transcript& ref, double p, double mass) {
    //logger_->warn("in update()");

  if (mass == salmon::math::LOG_0) {
    return;
  }
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return;
  }

  if(bam_aux_find(hit.read, "SA")) // Chimeric alignment. Ignore
    return;


  const double newMass = mass;
  { update(hit.read, primary.read, ref, p, mass, transitionProbs_);
  }


}

void ONTAlignmentModel::printModel(std::ostream& os) {
  // dbg d(os);

  // d << "Model\n";
  // for(size_t i = 0; i < std::max(errorModel_.size(), frontClipModel_.size()); ++i) {
  //   const auto errorP = errorModel_[i].mass != 0.0 ? (errorModel_[i].sum / errorModel_[i].mass) : 0.0;
  //   const auto fclipP = frontClipModel_[i].mass != 0.0 ? (frontClipModel_[i].sum / frontClipModel_[i].mass) : 0.0;
  //   const auto bclipP = backClipModel_[i].mass != 0.0 ? (backClipModel_[i].sum / backClipModel_[i].mass) : 0.0;
  //   if(errorP == 0.0 && fclipP == 0.0 && bclipP == 0.0) continue;

  //   const auto n = i * binLen;
  //   d << (i * binLen) << " - " << ((i+1) * binLen) 
  //     << ' ' << errorP << ' ' << errorModel_[i].mass.load()
  //     << ' ' << fclipP << ' ' << frontClipModel_[i].mass.load()
  //     << ' ' << bclipP << ' ' << backClipModel_[i].mass.load() << '\n';
  // }
  // d << "--------------\n";
}
