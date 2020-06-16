#ifndef __SB_MODEL_HPP__
#define __SB_MODEL_HPP__

#include <boost/iostreams/filtering_stream.hpp>

#include "UtilityFunctions.hpp"
//#include "jellyfish/mer_dna.hpp"
//#include "rapmap/Kmer.hpp"
#include "pufferfish/Kmer.hpp"
#include <Eigen/Dense>
#include <cmath>

//using Mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 4>;
using SBMer = combinelib::kmers::Kmer<32,4>;

class SBModel {
public:
  SBModel();

  SBModel(const SBModel&) = default;
  SBModel(SBModel&&) = default;
  SBModel& operator=(const SBModel&) = default;
  SBModel& operator=(SBModel&&) = default;

  bool writeBinary(boost::iostreams::filtering_ostream& out) const;

  inline int32_t contextBefore(bool rc) {
    return rc ? _contextRight : _contextLeft;
  }
  inline int32_t contextAfter(bool rc) {
    return rc ? _contextLeft : _contextRight;
  }

  bool addSequence(const char* seqIn, bool revCmp, double weight = 1.0);
  bool addSequence(const SBMer& sbmer, double weight);

  Eigen::MatrixXd& counts();
  Eigen::MatrixXd& marginals();

  double evaluateLog(const char* seqIn);
  double evaluateLog(const SBMer& sbmer);

  bool normalize();

  bool checkTransitionProbabilities();

  void combineCounts(const SBModel& other);

  void dumpConditionalProbabilities(std::ostream& os);

  int32_t getContextLength();

  template <typename CountVecT>
  bool train(CountVecT& kmerCounts, const uint32_t K);

  inline double evaluate(uint32_t kmer, uint32_t K) {
    std::vector<int32_t> _order{0, 0, 2, 2, 2, 2};
    double p{1.0};
    int32_t SK = static_cast<int32_t>(K);
    for (int32_t pos = 0; pos < SK - _order.back(); ++pos) {
      uint32_t offset = static_cast<uint32_t>(2 * (SK - (pos + 1) - _order[pos]));
      auto idx = _getIndex(kmer, offset, _order[pos]);
      p *= _probs(idx, pos);
    }
    return p;
  }

private:
  inline uint32_t _getIndex(uint32_t kmer, uint32_t offset, uint32_t _order) {
    kmer >>= offset;
    switch (_order) {
    case 0:
      return kmer & 0x3;
    case 1:
      return kmer & 0xF;
    case 2:
      return kmer & 0x3F;
    default:
      return 0;
    }
    return 0;
  }
  bool _trained;

  int32_t _contextLength;
  int32_t _contextLeft;
  int32_t _contextRight;

  Eigen::MatrixXd _probs;
  Eigen::MatrixXd _marginals;

  //Mer _mer;
  SBMer _sbmer;
  std::vector<int32_t> _order;
  std::vector<int32_t> _shifts;
  std::vector<int32_t> _widths;
  constexpr static const double _prior_prob = 1e-10;
};

#endif //__SB_MODEL_HPP__
