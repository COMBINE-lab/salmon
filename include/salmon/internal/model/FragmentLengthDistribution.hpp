/**
 *  lengthdistribution.h
 *  express
 *
 *  Created by Adam Roberts on 1/30/13.
 *  Copyright 2013 Adam Roberts. All rights reserved.
 */
// Modifications by Rob Patro; 2014

#ifndef FRAGMENT_LENGTH_DISTRIBUTION
#define FRAGMENT_LENGTH_DISTRIBUTION

#include <atomic>
#include <mutex>
#include <string>
#include <vector>

#include "SpinLock.hpp" // RapMap's with try_lock

/**
 * The LengthDistribution class keeps track of the observed length distribution.
 * It is initialized with a Gaussian prior with parameters specified by the
 * arguments to the constructor. An argument-specified binomial kernel is then
 * added for each observation. All mass values and probabilities are stored and
 * returned in log space (except in to_string).
 */
class FragmentLengthDistribution {
  /**
   * A private vector that stores the (logged) kernel values.
   **/
  std::vector<double> kernel_;
  /**
   * A private vector that stores the observed (logged) mass for each length.
   */
  std::vector<std::atomic<double>> hist_;

  /**
   * A private vector that stores the observed (logged) mass for each length.
   */
  std::vector<double> cachedCMF_;
  std::vector<double> cachedPMF_;
  volatile bool haveCachedCMF_;
  // std::mutex fldMut_;
  SpinLock sl_;

  /**
   * A private double that stores the total observed (logged) mass.
   */
  std::atomic<double> totMass_;
  /**
   * A private double that stores the (logged) sum of the product of observed
   * lengths and masses for quick mean calculations.
   */
  std::atomic<double> sum_;
  /**
   * A private int that stores the minimum observed length.
   */
  std::atomic<size_t> min_;
  /**
   * A size for internal binning of the lengths in the distribution.
   */
  size_t binSize_;

public:
  /**
   * LengthDistribution Constructor.
   * @param alpha double that sets the average pseudo-counts (logged).
   * @param max_val an integer that sets the maximum allowable length.
   * @param prior_mu a double for the mean of the prior gaussian distribution.
            If 0, a uniform distribution is used instead.
   * @param prior_sigma a double for the standard deviation of the prior
   *        gaussian distribution.
   * @param kernel_n a size_t specifying the number of trials in the kernel
   *        binomial distribution. Must be odd.
   * @param kernel_p a double specifying the success probability for the kernel
   *        binomial distribution.
   * @param bin_size a size_t specifying the size of bins to use internally to
   *        reduce the number of parameters in the distribution.
   */
  FragmentLengthDistribution(double alpha, size_t max_val, double prior_mu,
                             double prior_sigma, size_t kernel_n,
                             double kernel_p, size_t bin_size = 1);
  /**
   * An accessor for the maximum allowed length.
   * @return Max allowed length.
   */
  size_t maxVal() const;
  /**
   * An accessor for the minimum observed length (1 initially).
   * @return Minimum observed length.
   */
  size_t minVal() const;
  /**
   * An accessor for the mean length in the distribution.
   * @return Mean observed length.
   */
  double mean() const;
  /**
   * A member function that updates the distribution based on a new length
   * observation.
   * @param len an integer for the observed length.
   * @param mass a double for the mass (logged) to add.
   */
  void addVal(size_t len, double mass);
  /**
   * An accessor for the (logged) probability of a given length.
   * @param len an integer for the length to return the probability of.
   * @return (logged) probability of observing the given length.
   */
  double pmf(size_t len) const;
  /**
   * A member function that returns a (logged) cumulative mass for a given
   * length.
   * @param len an integer for the length to return the cmf value of.
   * @return (Logged) cmf value of length.
   */
  double cmf(size_t len) const;

  /**
   * A member function that caches the cumulative mass function into
   * a vector.  This should be called (once), when the fld will no
   * longer be updated.  It will make future calls to cmf(len)
   * much faster.
   */
  void cacheCMF();

  /**
   * A member function that returns a vector containing the (logged) cumulative
   * mass function *for the bins*.
   * @return (Logged) cmf of bins.
   */
  std::vector<double> cmf() const;
  // do the same thing as above, but use pmf to compute this CMF
  std::vector<double> cmf(const std::vector<double>& pmf) const;

  /**
   * A member function that fills in a a vector containing the (logged)
   * probability mass function *for the bins*, and the min and max values
   * @return (Logged) pmf of bins.
   */
  void dumpPMF(std::vector<double>& pmfOut, size_t& minV, size_t& maxV) const;

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

#endif
