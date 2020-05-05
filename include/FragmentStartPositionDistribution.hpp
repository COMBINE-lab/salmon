/**
 * Model the bias in read start positions across a set of transcripts.
 * This class is similar to and inspired by the FragmentLengthDistribution
 * class, which was itself modified from lengthdistribution.h ---
 * originally written by Adam Roberts as part of the eXpress software.
 * updated by Rob Patro; 2014, 2015
 */

#ifndef FRAGMENT_START_POSITION_DISTRIBUTION
#define FRAGMENT_START_POSITION_DISTRIBUTION

// #include "tbb/atomic.h"
#include <atomic>
#include <mutex>
#include <string>
#include <vector>

/**
 * The FragmentStartPositionDistribution class keeps track of the observed
 * fragment start position distribution. It is initialized with uniform prior
 * with with parameters specified by the arguments to the constructor. All mass
 * values and probabilities are stored and returned in log space (except in
 * to_string).
 */
class FragmentStartPositionDistribution {
  /**
   * A private vector that stores the observed (logged) mass for each length.
   */
  std::vector<std::atomic<double>> pmf_;
  std::vector<std::atomic<double>> cmf_;
  /**
   * A private double that stores the total observed (logged) mass.
   */
  std::atomic<double> totMass_;
  /**
   * A private double that stores the (logged) sum of the product of observed
   * lengths and masses for quick mean calculations.
   */
  //tbb::atomic<double> sum_;
  /**
   * The number of bins we consider within each transcript.
   */
  size_t numBins_;

  // Mutex for this distribution
  std::mutex fspdMut_;
  std::atomic<bool> isUpdated_;
  std::atomic<bool> allowUpdates_;
  std::atomic<uint32_t> performingUpdate_;

public:
  /**
   * FragmentStartPositionDistribution constructor:
   * @param numBins The number of bins to consider for
   *                each transcript.
   *
   */
  FragmentStartPositionDistribution(uint32_t numBins = 20);

  /**
   * A member function that updates the distribution based on a new length
   * observation.
   * @param len an integer for the observed length.
   * @param mass a double for the mass (logged) to add.
   */
  void addVal(int32_t hitPos, uint32_t txpLen, double mass);
  /**
   * A member function that returns the probability that a hit
   * starts at the specified position within the given transcript length.
   * @param hitPos The position where the fragment begins
   * @param txpLen The length of the transcript
   */
  double operator()(int32_t hitPos, uint32_t txpLen, double effLen);

  /**
   * A member function that computes the probability that a hit
   * starts at the specified position within the given transcript length.
   * The overall log probability is given by logNum - logDenom. The function
   * returns true if the probability is non-zero and false otherwise.
   * @param hitPos The position where the fragment begins
   * @param txpLen The length of the transcript
   * @param logEffLen the log of the effective length of the transcript
   * @param logNum the log of the numerator
   * @param logDenom the log of the denominator
   * @return true if the probaility is non-zero, false otherwise.
   */
  bool logNumDenomMass(int32_t hitPos, uint32_t txpLen, double logEffLen,
                       double& logNum, double& logDenom);

  // Evaluate the CDF between two points
  double evalCDF(int32_t hitPos, uint32_t txpLen);
  // Update the distribution (compute the CDF) and
  // set isUpdated_;
  void update();

  /**
   * An accessor for the (logged) probability of a given length.
   * @param len an integer for the length to return the probability of.
   * @return (logged) probability of observing the given length.
   */
  // double pmf(size_t len) const;
  /**
   * A member function that returns a (logged) cumulative mass for a given
   * length.
   * @param len an integer for the length to return the cmf value of.
   * @return (Logged) cmf value of length.
   */
  // double cmf(size_t len) const;
  /**
   * A member function that returns a vector containing the (logged) cumulative
   * mass function *for the bins*.
   * @return (Logged) cmf of bins.
   */
  // std::vector<double> cmf() const;
  /**
   * An accessor for the (logged) observation mass (including pseudo-counts).
   * @return Total observation mass.
   */
  double totMass() const;
  /**
   * A member function that returns a string containing the current
   * distribution.
   * @return Space-separated string of probabilities ordered from length 0 to
   *         max_val (non-logged).
   */
  std::string toString() const;
  // std::string to_string() const;
  /**
   * A member function that appends the LengthDistribution parameters to the end
   * of the given file.
   * @param outfile the file to append to.
   * @param length_type a string specifying the type of length the distribution
   *        is of (ie. "Fragment" or "Target") to be included in the header.
   */
  // void append_output(std::ofstream& outfile, std::string length_type) const;
};

#endif // FRAGMENT_START_POSITION_DISTRIBUTION
