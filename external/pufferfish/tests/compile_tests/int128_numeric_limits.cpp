#include <cstdint>
#include <limits>

// adopted from https://github.com/gmarcais/Jellyfish/blob/ab0f0f2ac8996ad397db89c9baf8b6c80a4c652a/configure.ac
template<bool> struct StaticAssert; 
template<> struct StaticAssert<true> { static void assert() { } };

int main(int argc, char* argv[]) {
    StaticAssert<std::numeric_limits<__int128>::is_specialized>::assert();
    return 0;
}
