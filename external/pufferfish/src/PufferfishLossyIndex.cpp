#include <fstream>
#include <iostream>

#include "CLI/Timer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "PufferfishLossyIndex.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"

PufferfishLossyIndex::PufferfishLossyIndex() {}

PufferfishLossyIndex::PufferfishLossyIndex(const std::string& indexDir, pufferfish::util::IndexLoadingOpts opts) {
  if (!puffer::fs::DirExists(indexDir.c_str())) {
    std::cerr << "The index directory " << indexDir << " does not exist!\n";
    std::exit(1);
  }

  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("k", k_));
    infoArchive(cereal::make_nvp("num_kmers", numKmers_));
    infoArchive(cereal::make_nvp("have_edge_vec", haveEdges_));
    infoArchive(cereal::make_nvp("have_ref_seq", haveRefSeq_));
    infoArchive(cereal::make_nvp("num_decoys", numDecoys_));
    infoArchive(cereal::make_nvp("first_decoy_index", firstDecoyIndex_));

    std::cerr << "k = " << k_ << '\n';
    std::cerr << "num kmers = " << numKmers_ << '\n';
    infoStream.close();
    twok_ = 2 * k_;
  } 
  haveEdges_ = opts.try_loading_edges and haveEdges_;
  haveRefSeq_ = opts.try_loading_ref_seqs and haveRefSeq_;
  haveEqClasses_ = opts.try_loading_eqclasses and haveEqClasses_;

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

  {
    CLI::AutoTimer timer{"Loading mphf table", CLI::Timer::Big};
    std::string hfile = indexDir + "/" + pufferfish::util::MPH;
    std::ifstream hstream(hfile);
    hash_.reset(new boophf_t);
    hash_->load(hstream);
    hstream.close();
    hash_raw_ = hash_.get();
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

  {
    CLI::AutoTimer timer{"Loading presence vector", CLI::Timer::Big};
    std::string bfile = indexDir + "/" + pufferfish::util::PRESENCE;
    presenceVec_.deserialize(bfile, false);
    presenceRank_ = rank9b(presenceVec_.get(), presenceVec_.size());//decltype(presenceVec_)::rank_1_type(&presenceVec_);
    //presenceSelect_ = decltype(presenceVec_)::select_1_type(&presenceVec_);
  }

  {
    CLI::AutoTimer timer{"Loading sampled positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::SAMPLEPOS;
    auto bits_per_element = compact::get_bits_per_element(pfile);
    sampledPos_.set_m_bits(bits_per_element);
    sampledPos_.deserialize(pfile, false);
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
  firstDecoyEncodedIndex_ = (numDecoys_ > 0) ? getRefId(firstDecoyIndex_) : std::numeric_limits<decltype(firstDecoyEncodedIndex_)>::max();
}

/**
 * Returns a ProjectedHits object containing all of the reference loci matching this
 * provided Canonical kmer (including the oritentation of the match).  The provided
 * QueryCache argument will be used to avoid redundant rank / select operations if feasible.
 */
auto PufferfishLossyIndex::getRefPos(CanonicalKmer& mer, pufferfish::util::QueryCache& qc)
    -> pufferfish::util::ProjectedHits {
  using IterT = std::vector<pufferfish::util::Position>::iterator;
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);

  if (res < numKmers_ and presenceVec_[res] == 1) {
    uint64_t pos{0};
    auto currRank = presenceRank_.rank(res);
    pos = sampledPos_[currRank];
    // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
    /* 
    uint64_t hashbits = pos & 0xF;
    pos = pos >> 4;
    if ((km & 0xF) != hashbits) { 
        return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
    }
    */
    uint64_t twopos = pos << 1;
    uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = rankSelDict.rank(pos);
      // the reference information in the contig table
      auto contigIterRange = contigRange(rank);
      // start position of this contig
      uint64_t sp = 0;
      uint64_t contigEnd = 0;
      if (rank == qc.prevRank) {
        sp = qc.contigStart;
        contigEnd = qc.contigEnd;
      } else {
        //sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
        //contigEnd = contigSelect_(rank);
        sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank - 1)) + 1;
        contigEnd = rankSelDict.select(rank);
        qc.prevRank = rank;
        qc.contigStart = sp;
        qc.contigEnd = contigEnd;
      }

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

auto PufferfishLossyIndex::getRefPos(CanonicalKmer& mer) -> pufferfish::util::ProjectedHits {
  using IterT = std::vector<pufferfish::util::Position>::iterator;
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);

  if (res < numKmers_ and presenceVec_[res] == 1) {
    uint64_t pos{0};
    auto currRank = presenceRank_.rank(res);
    pos = sampledPos_[currRank];

    uint64_t twopos = pos << 1;
    uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = rankSelDict.rank(pos);
      // the reference information in the contig table
      auto contigIterRange = contigRange(rank);
      // start and end position of this contig
      uint64_t sp =
        (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank - 1)) + 1;

      uint64_t contigEnd = rankSelDict.select(rank);
      //uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
      //uint64_t contigEnd =  contigSelect_(rank + 1);

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
