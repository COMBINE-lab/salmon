#include <bitset>
#include <fstream>
#include <iostream>

#include "CLI/Timer.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"

#include "PufferFS.hpp"
#include "PufferfishSparseIndex.hpp"

PufferfishSparseIndex::PufferfishSparseIndex() {}

PufferfishSparseIndex::PufferfishSparseIndex(const std::string& indexDir, pufferfish::util::IndexLoadingOpts opts) {
  if (!puffer::fs::DirExists(indexDir.c_str())) {
    std::cerr << "The index directory " << indexDir << " does not exist!\n";
    std::exit(1);
  }

  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("k", k_));
    infoArchive(cereal::make_nvp("num_kmers", numKmers_));
    infoArchive(cereal::make_nvp("num_sampled_kmers", numSampledKmers_));
    infoArchive(cereal::make_nvp("extension_size", extensionSize_));
    infoArchive(cereal::make_nvp("have_edge_vec", haveEdges_));
    infoArchive(cereal::make_nvp("have_ref_seq", haveRefSeq_));
    infoArchive(cereal::make_nvp("num_decoys", numDecoys_));
    infoArchive(cereal::make_nvp("first_decoy_index", firstDecoyIndex_));

    std::cerr << "k = " << k_ << '\n';
    std::cerr << "num kmers = " << numKmers_ << '\n';
    std::cerr << "num sampled kmers = " << numSampledKmers_ << '\n';
    std::cerr << "extension size = " << extensionSize_ << '\n';
    twok_ = 2 * k_;
    infoStream.close();
  }
  haveEdges_ = opts.try_loading_edges and haveEdges_;
  haveRefSeq_ = opts.try_loading_ref_seqs and haveRefSeq_;
  haveEqClasses_ = opts.try_loading_eqclasses and haveEqClasses_;
  
  // std::cerr << "loading contig table ... ";
  {
    CLI::AutoTimer timer{"Loading contig table", CLI::Timer::Big};
    std::ifstream contigTableStream(indexDir + "/" + pufferfish::util::CTABLE);
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(refNames_);
    contigTableArchive(refExt_);
    contigTableArchive(contigTable_);
    contigTableStream.close();
  }
  {
    CLI::AutoTimer timer{"Loading contig offsets", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::CONTIG_OFFSETS;
    auto bits_per_element = compact::get_bits_per_element(pfile);
    contigOffsets_.set_m_bits(bits_per_element);
    contigOffsets_.deserialize(pfile, false);
  }
  numContigs_ = contigOffsets_.size()-1;

  {
    std::string rlPath = indexDir + "/" + pufferfish::util::REFLENGTH;
    if (puffer::fs::FileExists(rlPath.c_str())) {
      CLI::AutoTimer timer{"Loading reference lengths", CLI::Timer::Big};
      std::ifstream refLengthStream(rlPath);
      cereal::BinaryInputArchive refLengthArchive(refLengthStream);
      refLengthArchive(refLengths_);
    } else {
      refLengths_ = std::vector<uint32_t>(refNames_.size(), 1000);
    }
  }

  {
    std::string rlPath = indexDir + "/" + pufferfish::util::COMPLETEREFLENGTH;
    if (puffer::fs::FileExists(rlPath.c_str())) {
      std::ifstream completeRefLengthStream(rlPath);
      cereal::BinaryInputArchive completeRefLengthArchive(completeRefLengthStream);
      completeRefLengthArchive(completeRefLengths_);
    } else {
      throw std::runtime_error("could not load complete reference lengths!");
    }
  }
  
  if (haveEqClasses_) {
    CLI::AutoTimer timer{"Loading eq table", CLI::Timer::Big};
    std::ifstream eqTableStream(indexDir + "/" + pufferfish::util::EQTABLE);
    cereal::BinaryInputArchive eqTableArchive(eqTableStream);
    eqTableArchive(eqClassIDs_);
    eqTableArchive(eqLabels_);
    eqTableStream.close();
  }
  // std::cerr << "done\n";

  {
    CLI::AutoTimer timer{"Loading mphf table", CLI::Timer::Big};
    std::string hfile = indexDir + "/" + pufferfish::util::MPH;
    std::ifstream hstream(hfile);
    hash_.reset(new boophf_t);
    hash_->load(hstream);
    hstream.close();
  }

  {
    CLI::AutoTimer timer{"Loading contig boundaries", CLI::Timer::Big};
    std::string bfile = indexDir + "/" + pufferfish::util::RANK;
    contigBoundary_.deserialize(bfile, false);
    rankSelDict = rank9sel(&contigBoundary_, (uint64_t)contigBoundary_.size());
  }

  {
    CLI::AutoTimer timer{"Loading sequence", CLI::Timer::Big};
    std::string sfile = indexDir + "/" + pufferfish::util::SEQ;
    seq_.deserialize(sfile, true);
    lastSeqPos_ = seq_.size() - k_;
  }

  if (haveRefSeq_) {
    CLI::AutoTimer timer{"Loading reference sequence", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::REFSEQ;
    refseq_.deserialize(pfile, true);
  }

  {
    std::string rlPath = indexDir + "/" + pufferfish::util::REFACCUMLENGTH;
    if (puffer::fs::FileExists(rlPath.c_str())) {
      CLI::AutoTimer timer{"Loading reference accumulative lengths", CLI::Timer::Big};
      std::ifstream refLengthStream(rlPath);
      cereal::BinaryInputArchive refLengthArchive(refLengthStream);
      refLengthArchive(refAccumLengths_);
    } else {
      refAccumLengths_ = std::vector<uint64_t>(refNames_.size(), 1000);
    }
  }

  if (haveEdges_) {
    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::EDGE;
    edge_.deserialize(pfile, true);
  }


  {
    CLI::AutoTimer timer{"Loading presence vector", CLI::Timer::Big};
    std::string bfile = indexDir + "/" + pufferfish::util::PRESENCE;
    presenceVec_.deserialize(bfile, false);
    presenceRank_ = rank9b(presenceVec_.get(), presenceVec_.size());
    std::cerr << "NUM 1s in presenceVec_ = " << presenceRank_.rank(presenceVec_.size()-1) << "\n\n";
    //presenceRank_ = decltype(presenceVec_)::rank_1_type(&presenceVec_);
    //presenceSelect_ = decltype(presenceVec_)::select_1_type(&presenceVec_);
  }
  {
    CLI::AutoTimer timer{"Loading canonical vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::CANONICAL;
    canonicalNess_.deserialize(pfile, false);
  }
  {
    CLI::AutoTimer timer{"Loading sampled positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::SAMPLEPOS;
    auto bits_per_element = compact::get_bits_per_element(pfile);
    sampledPos_.set_m_bits(bits_per_element);
    sampledPos_.deserialize(pfile, false);
  }

  {
    CLI::AutoTimer timer{"Loading extension vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::EXTENSION;
    auto bits_per_element = compact::get_bits_per_element(pfile);
    auxInfo_.set_m_bits(bits_per_element);
    auxInfo_.deserialize(pfile, false);
    std::string pfileSize = indexDir + "/" + pufferfish::util::EXTENSIONSIZE;
    bits_per_element = compact::get_bits_per_element(pfileSize);
    extSize_.set_m_bits(bits_per_element);
    extSize_.deserialize(pfileSize, false);
  }

  {
    CLI::AutoTimer timer{"Loading direction vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::DIRECTION;
    directionVec_.deserialize(pfile,false);
  }

  firstDecoyEncodedIndex_ = (numDecoys_ > 0) ? getRefId(firstDecoyIndex_) : std::numeric_limits<decltype(firstDecoyEncodedIndex_)>::max();
}

auto PufferfishSparseIndex::getRefPosHelper_(CanonicalKmer& mer, uint64_t pos,
                                             pufferfish::util::QueryCache& qc, bool didWalk)
    -> pufferfish::util::ProjectedHits {
  using IterT = std::vector<pufferfish::util::Position>::iterator;
  if (pos <= lastSeqPos_) {
    uint64_t fk = seq_.get_int(2*pos, 2*k_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      uint64_t rank = 0;
      uint64_t sp = 0;
      uint64_t contigEnd = 0;

      // make sure that the rank vector, from the 0th through k-1st position
      // of this k-mer is all 0s
      auto rankInterval =
          (didWalk) ? contigBoundary_.get_int(pos, (k_ - 1)) : 0;

      if (rankInterval > 0) {
        return {std::numeric_limits<uint32_t>::max(),
                std::numeric_limits<uint64_t>::max(),
                std::numeric_limits<uint32_t>::max(),
                true,
                0,
                k_,
                core::range<IterT>{}};
      }

      // check if the current position is on the same contig
      // as the cached contig or not.
      const bool same_contig = ((qc.contigStart <= pos) and (pos <= qc.contigEnd));
      if (same_contig) {
        rank = qc.prevRank;
        sp = qc.contigStart;
        contigEnd = qc.contigEnd;
      } else {
        // the index of this contig
        rank = rankSelDict.rank(pos);//contigRank_(pos);
        sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank - 1)) + 1;
        contigEnd = rankSelDict.select(rank);
        qc.prevRank = rank;
        qc.contigStart = sp;
        qc.contigEnd = contigEnd;
      }

      // the reference information in the contig table
      auto contigIterRange = contigRange(rank);

      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      // auto clen =
      // cPosInfo_[rank].length();//static_cast<uint64_t>(contigSelect_(rank +
      // 1) + 1 - sp);

      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
              contigIterRange};
      //core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      return {std::numeric_limits<uint32_t>::max(),
              std::numeric_limits<uint64_t>::max(),
              std::numeric_limits<uint32_t>::max(),
              true,
              0,
              k_,
              core::range<IterT>{}};
    }
  }

  return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
}

auto PufferfishSparseIndex::getRefPosHelper_(CanonicalKmer& mer, uint64_t pos,
                                             bool didWalk)
    -> pufferfish::util::ProjectedHits {

  using IterT = std::vector<pufferfish::util::Position>::iterator;
  if (pos <= lastSeqPos_) {
    uint64_t fk = seq_.get_int(2*pos, 2*k_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // make sure that the rank vector, from the 0th through k-1st position
      // of this k-mer is all 0s
      auto rankInterval =
          (didWalk) ? contigBoundary_.get_int(pos, (k_ - 1)) : 0;
      // auto rankEnd = contigRank_(pos + k_ - 1);
      if (rankInterval > 0) {
        return {std::numeric_limits<uint32_t>::max(),
                std::numeric_limits<uint64_t>::max(),
                std::numeric_limits<uint32_t>::max(),
                true,
                0,
                k_,
                core::range<IterT>{}};
      }

      // the index of this contig
      auto rank = rankSelDict.rank(pos);

      // the reference information in the contig table
      auto contigIterRange = contigRange(rank);

      // start position of this contig
      uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank - 1)) + 1;
      uint64_t contigEnd = rankSelDict.select(rank);


      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);

      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
              contigIterRange};
          //core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      return {std::numeric_limits<uint32_t>::max(),
              std::numeric_limits<uint64_t>::max(),
              std::numeric_limits<uint32_t>::max(),
              true,
              0,
              k_,
              core::range<IterT>{}};
    }
  }
  return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
}

auto PufferfishSparseIndex::getRefPos(CanonicalKmer mern, pufferfish::util::QueryCache& qc)
    -> pufferfish::util::ProjectedHits {
  using IterT = std::vector<pufferfish::util::Position>::iterator;
  pufferfish::util::ProjectedHits emptyHit{std::numeric_limits<uint32_t>::max(),
                               std::numeric_limits<uint64_t>::max(),
                               std::numeric_limits<uint32_t>::max(),
                               true,
                               0,
                               k_,
                               core::range<IterT>{}};

  bool didWalk{false};

  auto km = mern.getCanonicalWord();
  CanonicalKmer mer = mern;
  if (!mer.isFwCanonical()) {
    mer.swap();
  }

  // lookup this k-mer
  size_t idx = hash_->lookup(km);

  // if the index is invalid, it's clearly not present
  if (idx >= numKmers_) {
    return emptyHit;
  }

  uint64_t pos{0};
  //auto currRank = (idx == 0) ? 0 : presenceRank_.rank(idx);
  auto currRank = presenceRank_.rank(idx);

  if (presenceVec_[idx] == 1) {
    pos = sampledPos_[currRank];
  } else {
    didWalk = true;
    int signedShift{0};
    int inLoop = 0;

    /*
    do{
            if(inLoop >= 1){
        return emptyHit;
            }
    */

    auto extensionPos = idx - currRank;
    uint64_t extensionWord = auxInfo_[extensionPos];

    if (!canonicalNess_[extensionPos] and mer.isFwCanonical()) {
      mer.swap();
    }

    bool shiftFw = (directionVec_[extensionPos] == 1);
    // + 1 because we encode 1 as 00, 2 as 01, etc.
    int32_t llimit =
        extensionSize_ - static_cast<int32_t>(extSize_[extensionPos] + 1);

    if (shiftFw) {
      for (int32_t i = extensionSize_; i > llimit; --i) {
        uint32_t ssize = 2 * (i - 1);
        int currCode =
            static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
        mer.shiftFw(currCode);
        --signedShift;
      }
    } else {
      for (int32_t i = extensionSize_; i > llimit; --i) {
        uint32_t ssize = 2 * (i - 1);
        int currCode =
            static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
        mer.shiftBw(currCode);
        ++signedShift;
      }
    }

    km = mer.getCanonicalWord();
    idx = hash_->lookup(km);

    if (idx >= numKmers_) {
      return emptyHit;
    }
    // if we didn't find a present kmer after extension, this is a no-go
    if (presenceVec_[idx] != 1) {
      return emptyHit;
    }

    //currRank = (idx == 0) ? 0 : presenceRank_.rank(idx);
    currRank = presenceRank_.rank(idx);

    inLoop++;
    auto sampledPos = sampledPos_[currRank];
    pos = sampledPos + signedShift;
  }
  // end of sampling based pos detection
  return getRefPosHelper_(mern, pos, qc, didWalk);
}

auto PufferfishSparseIndex::getRefPos(CanonicalKmer mern)
    -> pufferfish::util::ProjectedHits {
  using IterT = std::vector<pufferfish::util::Position>::iterator;
  pufferfish::util::ProjectedHits emptyHit{std::numeric_limits<uint32_t>::max(),
                               std::numeric_limits<uint64_t>::max(),
                               std::numeric_limits<uint32_t>::max(),
                               true,
                               0,
                               k_,
                               core::range<IterT>{}};

  bool didWalk{false};

  auto km = mern.getCanonicalWord();
  CanonicalKmer mer = mern;
  if (!mer.isFwCanonical()) {
    mer.swap();
  }

  // lookup this k-mer
  size_t idx = hash_->lookup(km);

  // if the index is invalid, it's clearly not present
  if (idx >= numKmers_) {
    return emptyHit;
  }

  uint64_t pos{0};
  auto currRank = (idx == 0) ? 0 : presenceRank_.rank(idx);

  if (presenceVec_[idx] == 1) {
    pos = sampledPos_[currRank];
  } else {
    didWalk = true;
    int signedShift{0};
    int inLoop = 0;

    /*
    do{
            if(inLoop >= 1){
        return emptyHit;
            }
    */

    auto extensionPos = idx - currRank;
    uint64_t extensionWord = auxInfo_[extensionPos];

    if (!canonicalNess_[extensionPos] and mer.isFwCanonical()) {
      mer.swap();
    }

    bool shiftFw = (directionVec_[extensionPos] == 1);
    // + 1 because we encode 1 as 00, 2 as 01, etc.
    int32_t llimit =
        extensionSize_ - static_cast<int32_t>(extSize_[extensionPos] + 1);

    if (shiftFw) {
      for (int32_t i = extensionSize_; i > llimit; --i) {
        uint32_t ssize = 2 * (i - 1);
        int currCode =
            static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
        mer.shiftFw(currCode);
        --signedShift;
      }
    } else {
      for (int32_t i = extensionSize_; i > llimit; --i) {
        uint32_t ssize = 2 * (i - 1);
        int currCode =
            static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
        mer.shiftBw(currCode);
        ++signedShift;
      }
    }

    km = mer.getCanonicalWord();
    idx = hash_->lookup(km);

    if (idx >= numKmers_) {
      return emptyHit;
    }

    currRank = (idx == 0) ? 0 : presenceRank_.rank(idx);
    inLoop++;

    /*
      }while(presenceVec_[idx] != 1) ;
    */

    // if we didn't find a present kmer after extension, this is a no-go
    if (presenceVec_[idx] != 1) {
      return emptyHit;
    }
    auto sampledPos = sampledPos_[currRank];
    pos = sampledPos + signedShift;
  }
  // end of sampling based pos detection
  return getRefPosHelper_(mern, pos, didWalk);
}
