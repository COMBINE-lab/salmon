#ifndef SALMON_MATH_HPP
#define SALMON_MATH_HPP

#ifndef BOOST_LIKELY
#define BOOST_LIKELY(x) (x)
#endif

#ifndef BOOST_UNLIKELY
#define BOOST_UNLIKELY(x) (x)
#endif

#include <cmath>
#include <cassert>

namespace salmon {

    namespace math {
        constexpr double LOG_0 = HUGE_VAL;
        constexpr double LOG_1 = 0;
        constexpr double LOG_ONEHALF = -0.69314718055994530941;
        constexpr double LOG_ORPHAN_PROB = -2.30258509299404568401;
        constexpr double EPSILON = 0.375e-10;
        const double LOG_EPSILON = log(EPSILON);

        // Taken from https://github.com/adarob/eXpress/blob/master/src/main.h
        inline bool approxEqual(double a, double b, double eps=EPSILON) {
            return std::abs(a-b) <= eps;
        }

        // Taken from https://github.com/adarob/eXpress/blob/master/src/main.h
        inline double logAdd(double x, double y) {
            if (std::abs(x) == LOG_0) { return y; }
            if (std::abs(y) == LOG_0) { return x; }
            if (y > x) { std::swap(x,y); }
            double sum = x + std::log(1 + std::exp(y-x));
            return sum;
        }

        // Taken from https://github.com/adarob/eXpress/blob/master/src/main.h
        inline double logSub(double x, double y) {
            if (std::abs(y) == LOG_0) { return x; }
            if (x <= y) {
                assert(std::fabs(x-y) < 1e-5);
                return LOG_0;
            }
            double diff = x + std::log(1-std::exp(y-x));
            return diff;
        }


    }

}


#endif //SALMON_MATH_HPP
