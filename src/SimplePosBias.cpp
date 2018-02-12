#include "SimplePosBias.hpp"
#include "SalmonMath.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>

SimplePosBias::SimplePosBias(int32_t numBins, bool logSpace)
    : numBins_(numBins),
      masses_(numBins, (logSpace ? salmon::math::LOG_1 : 1.0)),
      isLogged_(logSpace) {}

// Add a mass of @mass to bin @bin
void SimplePosBias::addMass(int32_t bin, double mass) {
  masses_[bin] = salmon::math::logAdd(masses_[bin], mass);
}

// Compute the bin for @pos on a transcript of length @length,
// and add @mass to the appropriate bin
void SimplePosBias::addMass(int32_t pos, int32_t length, double mass) {
  double step = static_cast<double>(length) / numBins_;
  int bin = std::floor(pos / step);
  int msize = static_cast<int>(masses_.size());
  if (bin >= msize) {
    std::cerr << "bin = " << bin << '\n';
  }
  addMass(bin, mass);
}

// Project, the weights contained in "bins"
// into the vector @out (using spline interpolation)
void SimplePosBias::projectWeights(std::vector<double>& out) {
  auto len = out.size();
  for (size_t p = 0; p < len; ++p) {
    // The fractional sampling factor position p would have
    double fracP = static_cast<double>(p) / len;
    out[p] = std::max(0.001, s_(fracP));
  }
}

// Combine the distribution @other
// with this distribution
void SimplePosBias::combine(const SimplePosBias& other) {
  assert(other.masses_.size() == masses_.size());
  for (size_t i = 0; i < masses_.size(); ++i) {
    masses_[i] = salmon::math::logAdd(masses_[i], other.masses_[i]);
  }
}

// We're finished updating this distribution, so
// compute the cdf etc.
void SimplePosBias::finalize() {
  // convert from log space
  double sum{0.0};
  for (size_t i = 0; i < masses_.size(); ++i) {
    masses_[i] = std::exp(masses_[i]);
    sum += masses_[i];
  }
  // Account for mass at endpoints
  std::vector<double> splineMass(masses_.size() + 2);
  // Duplicate the first and last points as the end knots
  double startKnot = masses_.front() / sum;
  double stopKnot = masses_.back() / sum;
  double splineSum = sum + startKnot + stopKnot;
  splineMass[0] = startKnot;
  for (size_t i = 0; i < masses_.size(); ++i) {
    splineMass[i + 1] = (masses_[i] / splineSum);
    masses_[i] /= sum;
  }
  splineMass.back() = stopKnot;

  std::vector<double> splineBins(splineMass.size());
  splineBins[0] = 0.0;
  for (size_t i = 0; i < masses_.size(); ++i) {
    splineBins[i + 1] = positionBins_[i] - 0.01;
  }
  splineBins.back() = 1.0;
  //s_.set_points(splineBins, splineMass);
  s_ = tk::spline(splineBins, splineMass);
  isLogged_ = false;
  isFinalized_ = true;
}

// Seralize this model.
bool SimplePosBias::writeBinary(
    boost::iostreams::filtering_ostream& out) const {
  auto* mutThis = const_cast<SimplePosBias*>(this);
  // We shouldn't write out a non-finalized model
  if (!mutThis->isFinalized_) {
    auto l = spdlog::get("jointLog");
    l->error("Attempting to write out a non-finalized positional bias model. "
             "This should not happen.  Please report this bug on GitHub.");
    return false;
  }

  uint32_t modelLen = static_cast<uint32_t>(masses_.size());
  out.write(reinterpret_cast<char*>(&modelLen), sizeof(modelLen));
  out.write(reinterpret_cast<char*>(
                const_cast<decltype(masses_)::value_type*>(masses_.data())),
            sizeof(masses_.front()) * modelLen);
  return true;
}
