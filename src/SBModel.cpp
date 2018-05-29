#include "SBModel.hpp"
#include <sstream>
#include <utility>

SBModel::SBModel() : _trained(false) {
  // Roberts et al. model
  // _order = {0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0};
  //       -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8  9  10 11 12

  // Roberts et al. model (eXpress)
  // _order = {0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  //      -10 -9 -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8  9 10
  //_order = {0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  // 3, 3};

  //_order = {0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1};
  //         -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8  9

  // Simple model
  // _order = {0, 0, 1, 2, 2, 2, 2, 1, 0};
  _order = {0, 1, 2, 2, 2, 2, 2, 2, 2};
  //        -3 -2 -1  0  1  2  3  4  5

  // Short model
  //_order = {0, 1, 2, 2, 2, 2};
  //       -2 -1  0 1  2  3

  // The total length of the contexts we'll consider
  _contextLength = _order.size();
  // The number of bases before the read start position.
  _contextLeft = 3;
  // The number of bases after the read start position.
  _contextRight = 5;

  if ((_contextLeft + _contextRight + 1) != _contextLength) {
    std::cerr
        << "The left (" << _contextLeft << ") and right (" << _contextRight
        << ") "
           "context length (+1) didn't match the size of the model provided "
           "("
        << _contextLength << "); something is wrong!\n";
    std::exit(1);
  }

  _marginals = Eigen::MatrixXd(4, _contextLength);
  _marginals.setZero();

  _shifts.clear();
  _widths.clear();
  _shifts.reserve(_order.size());
  _widths.reserve(_order.size());
  for (int32_t i = 0; i < _contextLength; ++i) {
    _shifts.push_back((2 * _contextLength) - 2 * (i + 1));
    _widths.push_back(2 * (_order[i] + 1));
  }

  // Find the maximum order present in our model
  int32_t maxOrder{0};
  for (auto e : _order) {
    maxOrder = std::max(maxOrder, e);
  }

  // Set k equal to the size of the contexts we'll parse.
  _mer.k(_contextLength);

  // To hold all probabilities the matrix must be 4^{max_order + 1} by
  // context-length
  _probs = Eigen::MatrixXd(constExprPow(4, maxOrder + 1), _contextLength);
  // We have no intial observations
  _probs.setZero();
}

bool SBModel::writeBinary(boost::iostreams::filtering_ostream& out) const {
  auto* mutThis = const_cast<SBModel*>(this);
  out.write(reinterpret_cast<char*>(&(mutThis->_contextLength)),
            sizeof(int32_t));
  out.write(reinterpret_cast<char*>(&(mutThis->_contextLeft)), sizeof(int32_t));
  out.write(reinterpret_cast<char*>(&(mutThis->_contextRight)),
            sizeof(int32_t));
  // write the orders
  out.write(reinterpret_cast<char*>(mutThis->_order.data()),
            _contextLength * sizeof(int32_t));
  // write the shifts
  out.write(reinterpret_cast<char*>(mutThis->_shifts.data()),
            _contextLength * sizeof(int32_t));
  // write the widths
  out.write(reinterpret_cast<char*>(mutThis->_widths.data()),
            _contextLength * sizeof(int32_t));

  // Following adopted from:
  // http://stackoverflow.com/questions/25389480/how-to-write-read-an-eigen-matrix-from-binary-file
  // write all probabilities
  typename Eigen::MatrixXd::Index prows = _probs.rows(), pcols = _probs.cols();
  out.write(reinterpret_cast<char*>(&prows),
            sizeof(typename Eigen::MatrixXd::Index));
  out.write(reinterpret_cast<char*>(&pcols),
            sizeof(typename Eigen::MatrixXd::Index));
  out.write(reinterpret_cast<char*>(mutThis->_probs.data()),
            prows * pcols * sizeof(typename Eigen::MatrixXd::Scalar));

  // write marginal probabilities
  typename Eigen::MatrixXd::Index mrows = _marginals.rows(),
                                  mcols = _marginals.cols();
  out.write(reinterpret_cast<char*>(&mrows),
            sizeof(typename Eigen::MatrixXd::Index));
  out.write(reinterpret_cast<char*>(&mcols),
            sizeof(typename Eigen::MatrixXd::Index));
  out.write(reinterpret_cast<char*>(mutThis->_marginals.data()),
            mrows * mcols * sizeof(typename Eigen::MatrixXd::Scalar));

  return true;
}

double SBModel::evaluateLog(const char* seqIn) {
  double p = 0;
  Mer mer;
  mer.from_chars(seqIn);

  for (int32_t i = 0; i < _contextLength; ++i) {
    uint64_t idx = mer.get_bits(_shifts[i], _widths[i]);
    p += _probs(idx, i);
  }
  return p;
}

double SBModel::evaluateLog(const Mer& mer) {
  double p = 0;

  for (int32_t i = 0; i < _contextLength; ++i) {
    uint64_t idx = mer.get_bits(_shifts[i], _widths[i]);
    p += _probs(idx, i);
  }
  return p;
}

/** inlined member functions

inline int32_t SBModel::contextBefore(bool rc);
inline int32_t SBModel::contextAfter(bool rc);
inline bool SBModel::addSequence(const char* seqIn, bool revCmp, double weight
= 1.0); inline double SBModel::evaluate(const char* seqIn); inline double
SBModel::evaluate(uint32_t kmer, uint32_t K); inline uint32_t
SBModel::_getIndex(uint32_t kmer, uint32_t offset, uint32_t _order);

**/

Eigen::MatrixXd& SBModel::counts() { return _probs; }

Eigen::MatrixXd& SBModel::marginals() { return _marginals; }

void SBModel::dumpConditionalProbabilities(std::ostream& os) {
  typedef jellyfish::mer_dna_ns::mer_base_dynamic<uint64_t> mer64;
  // For each position
  for (int32_t i = 0; i < _contextLength; ++i) {
    mer64 k(_order[i] + 1);
    size_t nbit = 2 * (_order[i] + 1);
    int32_t N = constExprPow(4, _order[i] + 1);
    // header
    for (int32_t j = 0; j < N; ++j) {
      k.set_bits(0, nbit, j);
      std::string s = k.to_str();
      if (s.length() > 1) {
        os << s.substr(0, s.length() - 1) << " -> ";
        os << s.back();
      } else {
        os << "\'\' -> " << s.front();
      }
      if (j < N - 1) {
        os << '\t';
      }
    }
    os << '\n';
    // probs
    for (int32_t j = 0; j < N; ++j) {
      auto p = _probs(j, i);
      os << std::exp(p);
      if (j < N - 1) {
        os << '\t';
      }
    }
    os << '\n';
  }
}

bool SBModel::addSequence(const char* seqIn, bool revCmp, double weight) {
  _mer.from_chars(seqIn);
  if (revCmp) {
    _mer.reverse_complement();
  }
  return addSequence(_mer, weight);
}

bool SBModel::addSequence(const Mer& mer, double weight) {
  for (int32_t i = 0; i < _contextLength; ++i) {
    uint64_t idx = mer.get_bits(_shifts[i], _widths[i]);
    _probs(idx, i) += weight;
  }
  return true;
}

/**
 * Once the _prob matrix has been filled out with observations, calling
 * this function will normalize all counts so that the entries of _prob
 * represent proper transition probabilities.
 * NOTE: The _prob matrix can only be normalized once.
 *
 * returns : true if the matrix was normalized and false otherwise.
 **/
bool SBModel::normalize() {
  if (_trained) {
    return false;
  }

  // now normalize the rest of the sub-contexts in groups
  // each consecutive group of 4 rows shares the same prefix
  for (int32_t pos = 0; pos < _contextLength; ++pos) {
    int32_t numStates = constExprPow(4, _order[pos]);
    size_t rowsPerNode = 4;
    size_t nodeStart = 0;
    for (int32_t i = 0; i < numStates; ++i) {
      // Group the transition probabilities corresponding to
      // the current context.  Normalize them so they are
      // conditional probabilities.
      auto tot = _probs.col(pos).segment(nodeStart, rowsPerNode).sum();
      _probs.col(pos).segment(nodeStart, rowsPerNode) /= tot;

      _marginals(0, pos) += _probs(nodeStart, pos);
      _marginals(1, pos) += _probs(nodeStart + 1, pos);
      _marginals(2, pos) += _probs(nodeStart + 2, pos);
      _marginals(3, pos) += _probs(nodeStart + 3, pos);
      nodeStart += rowsPerNode;
    }
    _marginals.col(pos) /= numStates;
    // std::cerr << "pos = " << pos << ", marginals = " << _marginals.col(pos)
    // << '\n';
  }

  double logSmall = std::log(1e-5);
  auto takeLog = [logSmall](double x) -> double {
    return (x > 0.0) ? std::log(x) : logSmall;
  };
  _probs = _probs.unaryExpr(takeLog);
  _trained = true;
  return true;
}

bool SBModel::checkTransitionProbabilities() {
  if (!_trained) {
    return true;
  }

  // now normalize the rest of the sub-contexts in groups
  // each consecutive group of 4 rows shares the same 2-mer prefix
  for (int32_t pos = 0; pos < _contextLength; ++pos) {
    int32_t numStates = constExprPow(4, _order[pos]);
    size_t rowsPerNode = 4;
    size_t nodeStart = 0;
    for (int32_t i = 0; i < numStates; ++i) {
      auto tot = 0.0;
      for (size_t j = nodeStart; j < nodeStart + rowsPerNode; ++j) {
        tot += std::exp(_probs(j, pos));
      }
      if (tot < 0.98 or tot > 1.02) {
        std::cerr << "Transition probabilites for position " << i << ", rows["
                  << nodeStart << ", " << nodeStart + rowsPerNode
                  << "] = " << tot << '\n';
        return false;
      }
      nodeStart += rowsPerNode;
    }
  }

  return true;
}

void SBModel::combineCounts(const SBModel& other) { _probs += other._probs; }

int32_t SBModel::getContextLength() { return _contextLength; }

template <typename CountVecT>
bool SBModel::train(CountVecT& kmerCounts, const uint32_t K) {
  // The _order of the model *ending at* positions 2, 3, 4, and 5 (0-based)
  std::vector<int32_t> _order{0, 0, 2, 2, 2, 2};
  const auto numKmers = constExprPow(4, K);

  if (!_trained) {
    // For each starting position
    for (int32_t pos = 0; pos < static_cast<int32_t>(K - _order.back()); ++pos) {
      uint32_t offset = 2 * (K - (pos + 1) - _order[pos]);

      // See how frequently sub-contexts starting at this position appear
      for (uint32_t kmer = 0; kmer < numKmers; ++kmer) {
        auto idx = _getIndex(kmer, offset, _order[pos]);
        _probs(idx, pos) += static_cast<double>(kmerCounts[kmer]);
      }
    }

    // Normalize the first column (all 3-mers)
    int32_t startIdx = 0;
    _probs.col(0) /= _probs.col(0).sum();
    _probs.col(1) /= _probs.col(1).sum();
    _probs.col(2) /= _probs.col(2).sum();
    // now normalize the rest of the sub-contexts in groups
    // each consecutive group of 4 rows shares the same 2-mer prefix
    for (int32_t pos = 3; pos < static_cast<int32_t>(K - _order.back()); ++pos) {
      int32_t numStates = constExprPow(4, _order[pos]);
      size_t rowsPerNode = 4;
      size_t nodeStart = 0;
      for (int32_t i = 0; i < numStates; ++i) {
        auto tot = _probs.col(pos).segment(nodeStart, rowsPerNode).sum();
        _probs.col(pos).segment(nodeStart, rowsPerNode) /= tot;
        nodeStart += rowsPerNode;
      }
    }
    _trained = true;
  }
  return true;
}

template bool
SBModel::train<std::array<std::atomic<uint32_t>, constExprPow(4, 6)>>(
    std::array<std::atomic<uint32_t>, constExprPow(4, 6)>& counts,
    const uint32_t K);

template bool SBModel::train<std::vector<double>>(std::vector<double>& counts,
                                                  const uint32_t K);
