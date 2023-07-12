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

//Calculate number of states based on kmerLength (should probably be moved to the header file)
static const uint32_t kmerLength = 50;
static const uint32_t stepSize = 50;
static const uint32_t numMismatchStates = kmerLength + 1;
constexpr static uint32_t getStateIndexNoHomopolymer(size_t numIndels, size_t numMisMatch) {
    return numMismatchStates * numIndels + numMisMatch;
}

uint32_t numStatesNoHomo = getStateIndexNoHomopolymer(kmerLength, 0)+1; 

uint32_t getStateIndex(size_t numIndels, size_t numMisMatch, bool homopolymer ) {
     if (homopolymer){
      return getStateIndexNoHomopolymer(numIndels,numMisMatch) + numStatesNoHomo;
     }
     return getStateIndexNoHomopolymer(numIndels,numMisMatch);
    }
uint32_t startStateIdx = getStateIndex(kmerLength, 0, true) + 1;
uint32_t endStateIdx = getStateIndex(kmerLength, 0, true) + 2;
uint32_t numStatesTotal = getStateIndex(kmerLength, 0, true) + 3;


ONTAlignmentModel::ONTAlignmentModel(double alpha, uint32_t readBins)
    : isEnabled_(true)
    , readBins_(readBins)
    , printed(false)
    , transitionProbs_(readBins)
    , frontClipModel_(maxReadLen / binLen + 1)
    , backClipModel_(maxReadLen / binLen + 1)
    , transcriptFrontModel_()
    , transcriptBackModel_()
{ 
  for (size_t i = 0; i < readBins; ++i) {
    transitionProbs_[i] = std::move(AtomicMatrix<double>(
        numStatesTotal, numStatesTotal, alpha));
  }

}



struct Column
{uint8_t ref_base;
uint8_t read_base;
};

bool is_homopolymer(Column chunk[]){
  size_t counts[4] = {0};
  size_t numGaps = 0;
  for (size_t i = 0; i < kmerLength; i++) {
    size_t symbol = chunk[i].ref_base;
    if (symbol <= 3){
     counts[symbol]+=1;
    }
    else {
      numGaps += 1;
    }
  }
  for (size_t i = 0; i < 4; i++){
    if (counts[i] >= kmerLength - 1){
      return true;
    }
  }
  return false;
}
//Should return the state index based on the kmer and move the end bases to the middle
uint32_t processChunk(Column chunk[]){
  size_t numIndels = 0;
  size_t numMismatch = 0;

  for  (size_t i = 0; i < kmerLength; i++){
    if (chunk[i].ref_base >3 ||chunk[i].read_base>3 ){
      numIndels ++; 
    }
    else if(chunk[i].ref_base != chunk[i].read_base){
      numMismatch++;
    }
  }
  uint32_t stateIdx = getStateIndex(numIndels, numMismatch, is_homopolymer(chunk));

  for (size_t i = 0; i < kmerLength-stepSize; i++){
   chunk[i].ref_base = chunk[i+stepSize].ref_base;
   chunk[i].read_base = chunk[i+stepSize].read_base;
  }
  return stateIdx;
}
//Testing functions (numIndels and numMismatch)
uint32_t numIndels(Column chunk[]){
  size_t numIndels = 0;
  size_t numMismatch = 0;

  for  (size_t i = 0; i < kmerLength; i++){
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

  for  (size_t i = 0; i < kmerLength; i++){
    if (chunk[i].ref_base >3 ||chunk[i].read_base>3 ){
      numIndels ++; 
    }
    else if(chunk[i].ref_base != chunk[i].read_base){
      numMismatch++;
    }
  }
  return numMismatch;
}

std::string chunkString(Column chunk[]){
  using namespace salmon::stringtools;
  std::string ref_string = "";
  std::string read_string = "";
  
  for (size_t i = 0;i < kmerLength; i++ ){
    if(chunk[i].ref_base > 3){
      ref_string = ref_string +  "-";
    }
    else{
      ref_string = ref_string + twoBitToChar[chunk[i].ref_base];
    }
    if(chunk[i].read_base > 3){
      read_string = read_string +  "-";
    }
    else{
      read_string =  read_string + twoBitToChar[chunk[i].read_base];
    }
  }
 return ("Ref: " + ref_string + " Read: " + read_string );
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

  // std::stringstream readStream, matchStream, refStream;

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
      if (chunk_idx >= kmerLength)
        {
          chunk_idx = kmerLength-stepSize;
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
    //Add the transition probability for ending (was here)
  }
  //double tp = transitionProbs[readPosBin](curStateIdx, endStateIdx);
  //logLike += tp;
  //bgLogLike += transitionProbs[readPosBin](0, 0);
  //logger_->warn("( logLike: {} bgLogLike: {}",logLike, bgLogLike);

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

  ErrorCount counts;
  if(!computeErrorCount(hit.read, primary.read, ref, counts, "logLikelihood")) {
    if(logger_)
      logger_->warn("in logLikelihood() error parsing CIGAR string");
    return LOG_1;
  }

  const uint32_t cigarRLen = alnLen(hit, primary); // Read length minus hard clips
  if(counts.sclips() >= cigarRLen)
    return LOG_0; // Empty alignment!

  const uint32_t alignLen  = cigarRLen - counts.clips(); // Length of aligned part (no soft clip)
  const double   errorRate = (double)counts.ims() / alignLen;
  const int32_t  errorBin  = std::min(alignLen / binLen, (uint32_t)errorModel_.size() - 1);
  const int32_t  frontClipBin   = std::min(cigarRLen / binLen, (uint32_t)frontClipModel_.size() - 1);
  const int32_t  backClipBin    = std::min(cigarRLen / binLen, (uint32_t)backClipModel_.size() - 1);
  const auto&    frontClipAvg   = frontClipModel_[frontClipBin];
  const auto&    backClipAvg    = backClipModel_[backClipBin];

  double errorllh = LOG_1, frontClipllh = LOG_1, backClipllh = LOG_1; // Error and clip Log Likelihood (front and back)

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

  // Likelihood to have so many bases soft clipped based on the
  // average error rate. Don't penalize for having fewer clipped bases
  // than average, only if more.
  // front clips:
  if(frontClipAvg.sum > dmin && frontClipAvg.mass > dmin) {
    using        boost::math::geometric;
    const double  mean          = frontClipAvg.sum / frontClipAvg.mass;
    geometric     clipDist(1.0 / (mean + 1.0));
    const int32_t rmean         = std::round(mean);
    const auto    clips         = counts.fclips();
    const double clipLikelihood =
      clips <= rmean
      ? 1.0
      : (1.0 - cdf(clipDist, clips)) / (1.0 - cdf(clipDist, rmean));
    frontClipllh = clipLikelihood < llMin ? LOG_0 : std::log(clipLikelihood);
  } else {
    if(logger_)
      logger_->warn("read {} (length {}) has no trained clipping model",
                    bam_name(hit.read), cigarRLen);
  }

  // back clips:
  if(backClipAvg.sum > dmin && backClipAvg.mass > dmin) {
    using        boost::math::geometric;
    const double  mean          = backClipAvg.sum / backClipAvg.mass;
    geometric     clipDist(1.0 / (mean + 1.0));
    const int32_t rmean         = std::round(mean);
    const auto    clips         = counts.bclips();
    const double clipLikelihood =
      clips <= rmean
      ? 1.0
      : (1.0 - cdf(clipDist, clips)) / (1.0 - cdf(clipDist, rmean));
    backClipllh = clipLikelihood < llMin ? LOG_0 : std::log(clipLikelihood);
  } else {
    if(logger_)
      logger_->warn("read {} (length {}) has no trained clipping model",
                    bam_name(hit.read), cigarRLen);
  }

  //Transcript clipping model -- Considers soft clips to be "not aligned"/not cover that part of the transcript
   double transcriptFrontllh = LOG_1, transcriptBackllh = LOG_1;
auto transcriptIdx = bam_pos(hit.read);
  auto readAlnStart = transcriptIdx;
  if (transcriptIdx<0){
    transcriptIdx = 0;
   }
   size_t transcriptLen = ref.RefLength;


  int32_t numTranscriptFrontExcluded = static_cast<int32_t>(transcriptIdx );
  //auto alignmentEndIdx = transcriptIdx + counts.fclips() + alignLen - counts.hclips(); 
  if(transcriptIdx < 0){
    transcriptIdx = 0;
  }
    auto alignmentEndIdx = readAlnStart + alignLen; 

  int32_t numTranscriptBackExcluded  = static_cast<int32_t>(transcriptLen - alignmentEndIdx);
  //Update transcript front clip model
  if (numTranscriptFrontExcluded + alignLen + numTranscriptBackExcluded !=  transcriptLen){
    logger_->warn("in loglikelihood, number of bases of read seems inconsistent, transcript length: {}, transcript bases excluded at front: {} , transcript bases aligned: {}, transcript bases excluded at back: {}", 
  transcriptLen, numTranscriptFrontExcluded, alignLen,  numTranscriptBackExcluded);
  }
  
  if(transcriptFrontModel_.sum > dmin && transcriptFrontModel_.mass > dmin) {
    using        boost::math::geometric;
    const double  mean          = transcriptFrontModel_.sum / transcriptFrontModel_.mass;
    geometric     clipDist(1.0 / (mean + 1.0));
    const int32_t rmean         = std::round(mean);
    const auto    clips         = numTranscriptFrontExcluded;
    const double clipLikelihood =
      clips <= rmean
      ? 1.0
      : (1.0 - cdf(clipDist, clips)) / (1.0 - cdf(clipDist, rmean));
     transcriptFrontllh = clipLikelihood < llMin ? LOG_0 : std::log(clipLikelihood);
  } else {
    if(logger_)
      logger_->warn("read {} (length {}) has no trained transcript front clipping model",
                    bam_name(hit.read), cigarRLen);
  }

  if(transcriptBackModel_.sum > dmin && transcriptBackModel_.mass > dmin) {
    using        boost::math::geometric;
    const double  mean          = transcriptBackModel_.sum / transcriptBackModel_.mass;
    geometric     clipDist(1.0 / (mean + 1.0));
    const int32_t rmean         = std::round(mean);
    const auto    clips         = numTranscriptBackExcluded;
    const double clipLikelihood =
      clips <= rmean
      ? 1.0
      : (1.0 - cdf(clipDist, clips)) / (1.0 - cdf(clipDist, rmean));
     transcriptBackllh = clipLikelihood < llMin ? LOG_0 : std::log(clipLikelihood);
  } else {
    if(logger_)
      logger_->warn("read {} (length {}) has no trained transcript back clipping model",
                    bam_name(hit.read), cigarRLen);
  }

//logger_->warn("( errorllh: {} frontClipllh: {} backClipllh: {}  ",errorllh, frontClipllh,backClipllh );
  return errorllh + frontClipllh + backClipllh + transcriptFrontllh + transcriptBackllh;
  //return errorllh + transcriptFrontllh + transcriptBackllh;
  //return errorllh;
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
    /*if (chunk_idx == 0){
      logger_->warn("new alignment start");
    }*/ 
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
            //logger_->warn("(in update()) CIGAR = {}", cigarStream.str());
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
      
     // logger_->warn("chunk index: {}, curRefBase: {} curReadBase: {}",chunk_idx, curRefBase, curReadBase);

      chunk[chunk_idx].ref_base = curRefBase;
      chunk[chunk_idx].read_base = curReadBase;
     
      //chunk[chunk_idx].num_mismatch = chunk_idx == 0 ? 0 : chunk[chunk_idx-1].num_mismatch;
      //chunk[chunk_idx].num_indel = chunk_idx == 0 ? 0 : chunk[chunk_idx-1].num_indel;
       //Match/mismatch case
      /*if (BAM_CONSUME_SEQ(op) && BAM_CONSUME_REF(op)) {
        if(curRefBase != curReadBase){
          chunk[chunk_idx].num_mismatch += 1;
        }
        }
      //Indel case
      else {
        chunk[chunk_idx].num_indel += 1;
      }*/
      chunk_idx ++;
      if (chunk_idx >= kmerLength)
        {
             
        //logger_-> warn("Update Chunk: {},  num_indel: {}, num mismatch: {}, homopolymer: {} ",chunkString(chunk),  numIndels(chunk), numMismatch(chunk), is_homopolymer(chunk));
          curStateIdx = processChunk(chunk);
          //logger_->warn("Log Likelihood Current State Index: ",curStateIdx, "out of: ",  numStatesTotal);
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

// Update the probability model. The reads are binned based on their
// length (the length after soft clipping). For each bin, we assume a
// Binomial distribution B(p,n) (where n is the length of the read).
//
// The probability p in each bin is estimated as the empirical mean of
// the number of errors in the reads (all reads counted the same at
// this point: indels and mutations).
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

  ErrorCount counts;
  if(!computeErrorCount(hit.read, primary.read, ref, counts, "update")) {
    if(logger_)
      logger_->warn("in update() error parsing CIGAR string");
    return;
  }

  // Update error model
  // Not taking p and mass into account. What's up with those?
  const int32_t readLen   = alnLen(hit, primary);
  const int32_t alignLen  = readLen - counts.sclips();
  const double  errorRate = (double)counts.ims() / alignLen;
  const double  clipRateFront  = (double)counts.fclips() / (readLen + counts.hclips());
  const double  clipRateBack   = (double)counts.bclips() / (readLen + counts.hclips());
  if(errorRate > 1.0 || clipRateFront > 1.0 || clipRateBack > 1.0) { // Should not happen
    if (logger_) {
      logger_->warn("(in update()) CIGAR string for read [{}] "
                    "seems inconsistent. It implied an error rate "
                    "greater than 1: {} {} {}",
                    bam_name(hit.read), errorRate, clipRateFront, clipRateBack);
    }
    return;
  }

  const double newMass = mass;
  { update(hit.read, primary.read, ref, p, mass, transitionProbs_);
  }

  // Update front clip model
  { int32_t binIndex = std::min(readLen / binLen, (uint32_t)frontClipModel_.size() - 1);
    auto& bin = frontClipModel_[binIndex];
    salmon::utils::incLoop(bin.mass, newMass);
    salmon::utils::incLoop(bin.sum, (binIndex + 1) * binLen * newMass * clipRateFront);
  }

  // Update back clip model
  { int32_t binIndex = std::min(readLen / binLen, (uint32_t)backClipModel_.size() - 1);
    auto& bin = backClipModel_[binIndex];
    salmon::utils::incLoop(bin.mass, newMass);
    salmon::utils::incLoop(bin.sum, (binIndex + 1) * binLen * newMass * clipRateBack);
  }

  auto transcriptIdx = bam_pos(hit.read);
  auto readAlnStart = transcriptIdx;
  if (transcriptIdx<0){
    transcriptIdx = 0;
   }
   size_t transcriptLen = ref.RefLength;


  int32_t numTranscriptFrontExcluded = static_cast<int32_t>(transcriptIdx );
  //auto alignmentEndIdx = transcriptIdx + counts.fclips() + alignLen - counts.hclips(); 
  if(transcriptIdx < 0){
    transcriptIdx = 0;
  }
    auto alignmentEndIdx = readAlnStart + alignLen; 

  int32_t numTranscriptBackExcluded  = static_cast<int32_t>(transcriptLen - alignmentEndIdx);
  //const double backExcludeRate = numTranscriptBackExcluded / transcriptLen;
  //const double frontExcludeRate = numTranscriptFrontExcluded / transcriptLen;

  //Update transcript front clip model
  if (numTranscriptFrontExcluded + alignLen + numTranscriptBackExcluded !=  transcriptLen){
    logger_->warn("in update, number of bases of read seems inconsistent, transcript length: {}, transcript bases excluded at front: {} , transcript bases aligned: {}, transcript bases excluded at back: {}", 
  transcriptLen, numTranscriptFrontExcluded, alignLen,  numTranscriptBackExcluded);
  }
  
  {
    salmon::utils::incLoop(transcriptFrontModel_.mass, newMass);
    salmon::utils::incLoop(transcriptFrontModel_.sum, newMass * numTranscriptFrontExcluded);

  }
  //Update transcript back clip model

  {
    salmon::utils::incLoop(transcriptBackModel_.mass, newMass);
    salmon::utils::incLoop(transcriptBackModel_.sum, newMass * numTranscriptBackExcluded);
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
