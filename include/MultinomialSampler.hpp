#ifndef _MULTINOMIAL_SAMPLER_HPP_
#define _MULTINOMIAL_SAMPLER_HPP_

#include <random>
#include <vector>

class MultinomialSampler {
 public:
     MultinomialSampler(std::random_device& rd) :
         gen_(rd()), u01_(0.0, 1.0) {}

     void operator()(
             std::vector<int>::iterator sampleBegin,
             uint32_t n,
             uint32_t k,
             std::vector<double>::iterator probsBegin,
	     bool clearCounts = true) {
         int i, j;
         double u, sum;
         std::vector<double> z(k+1, 0.0);

	 if (clearCounts) {
		 for (uint32_t i = 0; i < k; i++) {
		     *(sampleBegin + i) = 0;
		 }
	 }

         z[0] = 0;
         for (i = 1; i <= k; i++) {
             sum = 0;
             for (j = 0; j < i; j++) sum+= *(probsBegin + j);
             z[i] = sum;
         }

         for (j = 0; j < n; j++) {
             u = u01_(gen_);

             for (i = 0; i < k; i++) {
                 if ((z[i] < u) && (u <= z[i+1])) {
                     (*(sampleBegin + i))++;
                 }
             }
         }
     }


private:
   std::mt19937 gen_;
   std::uniform_real_distribution<> u01_;
};

#endif //_MULTINOMIAL_SAMPLER_HPP_

