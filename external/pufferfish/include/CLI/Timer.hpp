#pragma once

// Distributed under the 3-Clause BSD License.  See accompanying
// file LICENSE or https://github.com/CLIUtils/CLI11 for details.

#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <utility>

namespace CLI {

class Timer {
  protected:
    /// This is a typedef to make clocks easier to use
    using clock = std::chrono::steady_clock;

    /// This typedef is for points in time
    using time_point = std::chrono::time_point<clock>;

    /// This is the type of a printing function, you can make your own
    using time_print_t = std::function<std::string(std::string, std::string)>;

    /// This is the title of the timer
    std::string title_;

    /// This is the function that is used to format most of the timing message
    time_print_t time_print_;

    /// This is the starting point (when the timer was created)
    time_point start_;

    /// This is the number of times cycles (print divides by this number)
    size_t cycles{1};

  public:
    /// Standard print function, this one is set by default
    static std::string Simple(std::string title, std::string time) { return title + ": " + time; }

    /// This is a fancy print function with --- headers
    static std::string Big(std::string title, std::string time) {
        return std::string("-----------------------------------------\n") + "| " + title + " | Time = " + time + "\n" +
               "-----------------------------------------";
    }

  public:
    /// Standard constructor, can set title and print function
    Timer(std::string title = "Timer", time_print_t time_print = Simple)
        : title_(std::move(title)), time_print_(std::move(time_print)), start_(clock::now()) {}

    /// Time a function by running it multiple times. Target time is the len to target.
    std::string time_it(std::function<void()> f, double target_time = 1) {
        time_point start = start_;
        double total_time;

        start_ = clock::now();
        size_t n = 0;
        do {
            f();
            std::chrono::duration<double> elapsed = clock::now() - start_;
            total_time = elapsed.count();
        } while(n++ < 100 && total_time < target_time);

        std::string out = make_time_str(total_time / n) + " for " + std::to_string(n) + " tries";
        start_ = start;
        return out;
    }

    /// This formats the numerical value for the time string
    std::string make_time_str() const {
        time_point stop = clock::now();
        std::chrono::duration<double> elapsed = stop - start_;
        double time = elapsed.count() / cycles;
        return make_time_str(time);
    }

    // LCOV_EXCL_START
    std::string make_time_str(double time) const {
        auto print_it = [](double x, std::string unit) {
            char buffer[50];
            std::snprintf(buffer, 50, "%.5g", x);
            return buffer + std::string(" ") + unit;
        };

        if(time < .000001)
            return print_it(time * 1000000000, "ns");
        else if(time < .001)
            return print_it(time * 1000000, "us");
        else if(time < 1)
            return print_it(time * 1000, "ms");
        else
            return print_it(time, "s");
    }
    // LCOV_EXCL_END

    /// This is the main function, it creates a string
    std::string to_string() const { return time_print_(title_, make_time_str()); }

    /// Division sets the number of cycles to divide by (no graphical change)
    Timer &operator/(size_t val) {
        cycles = val;
        return *this;
    }
};

/// This class prints out the time upon destruction
class AutoTimer : public Timer {
  public:
    /// Reimplementing the constructor is required in GCC 4.7
    AutoTimer(std::string title = "Timer", 
    time_print_t time_print = Simple, std::ostream& in=std::cerr) : Timer(title, time_print), ostr(in) {}
    // GCC 4.7 does not support using inheriting constructors.

    /// This desctructor prints the string
    ~AutoTimer() { ostr << to_string() << std::endl; }
  private:
    std::ostream& ostr;
};

} // namespace CLI

/// This prints out the time if shifted into a std::cout like stream.
inline std::ostream &operator<<(std::ostream &in, const CLI::Timer &timer) { return in << timer.to_string(); }
