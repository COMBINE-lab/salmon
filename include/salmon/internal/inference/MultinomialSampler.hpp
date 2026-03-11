#ifndef _MULTINOMIAL_SAMPLER_HPP_
#define _MULTINOMIAL_SAMPLER_HPP_

#include <algorithm>
#include <random>
#include <vector>

class MultinomialSampler {
public:
  MultinomialSampler(std::random_device& rd) : gen_(rd()), u01_(0.0, 1.0) {}

  void operator()(std::vector<uint64_t>::iterator sampleBegin, uint32_t n,
                  uint32_t k, std::vector<double>::iterator probsBegin,
                  bool clearCounts = true) {
    uint32_t i, j;
    double u, sum;
    std::vector<double> z(k + 1, 0.0);

    if (clearCounts) {
      for (i = 0; i < k; i++) {
        *(sampleBegin + i) = 0;
      }
    }

    z[0] = 0;
    for (i = 1; i <= k; i++) {
      sum = 0;
      for (j = 0; j < i; j++)
        sum += *(probsBegin + j);
      z[i] = sum;
    }

    // If k is small (<= 100), linear search is usually faster
    if (k <= 100) {
      for (j = 0; j < n; j++) {
        u = u01_(gen_);

        for (i = 0; i < k; i++) {
          if ((z[i] < u) && (u <= z[i + 1])) {
            (*(sampleBegin + i))++;
            break;
          }
        }
      }
    } else { // k is large enough to warrant binary search
      for (j = 0; j < n; j++) {
        u = u01_(gen_);

        // Find the offset of the element to increment
        auto it = std::lower_bound(z.begin(), z.end() - 1, u);
        size_t offset = static_cast<size_t>(std::distance(z.begin(), it));

        if (*it > u and offset > 0) {
          offset -= 1;
        }

        (*(sampleBegin + offset))++;
      }
    }
  }

private:
  std::mt19937 gen_;
  std::uniform_real_distribution<> u01_;
};

#endif //_MULTINOMIAL_SAMPLER_HPP_
