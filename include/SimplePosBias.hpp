#ifndef SIMPLE_POS_BIAS_HPP
#define SIMPLE_POS_BIAS_HPP

#include "spdlog/spdlog.h"
#include "spline.h"
#include <array>
#include <boost/iostreams/filtering_stream.hpp>
#include <vector>

class SimplePosBias {
public:
  SimplePosBias(int32_t numBins = 20, bool logSpace = true);

  // Add a mass of @mass to bin @bin
  void addMass(int32_t bin, double mass);

  // Compute the bin for @pos on a transcript of length @length,
  // and add @mass to the appropriate bin
  void addMass(int32_t pos, int32_t length, double mass);

  // Project, via linear interpolation, the weights contained in "bins"
  // into the vector @out.
  void projectWeights(std::vector<double>& out);

  // Combine the distribution @other
  // with this distribution
  void combine(const SimplePosBias& other);

  // We're finished updating this distribution, so
  // compute the cdf etc.
  void finalize();

  // Seralize this model.
  bool writeBinary(boost::iostreams::filtering_ostream& out) const;

private:
  int32_t numBins_;
  std::vector<double> masses_;
  bool isLogged_{true};
  bool isFinalized_{false};
  ::tk::spline s_;
  // position bins taken from Cufflinks:
  // https://github.com/cole-trapnell-lab/cufflinks/blob/master/src/biascorrection.cpp
  const std::vector<double> positionBins_{{.02, .04, .06, .08, .10, .15, .2,
                                           .3,  .4,  .5,  .6,  .7,  .8,  .85,
                                           .9,  .92, .94, .96, .98, 1.0}};
};

#endif // SIMPLE_POS_BIAS_HPP
