#include "UtilityFunctions.hpp"

// from http://stackoverflow.com/questions/2380962/generate-all-combinations-of-arbitrary-alphabet-up-to-arbitrary-length
std::vector<std::string> getAllWords(int length) {

    int N_LETTERS = 4;
    char alphabet[] = {'A', 'C', 'G', 'T'};
  std::vector<int> index(length, 0);
  std::vector<std::string> words;

  while(true)
  {
    std::string word(length, ' ');
    for (int i = 0; i < length; ++i)
      word[i] = alphabet[index[i]];
    words.push_back(word);

    for (int i = length-1; ; --i)
    {
      if (i < 0) return words;
      index[i]++;
      if (index[i] == N_LETTERS)
        index[i] = 0;
      else
        break;
    }
  }
}

#include <atomic>

SCENARIO("Kmers encode and decode correctly") {
    using salmon::utils::Direction;
    GIVEN("All 6-mers") {
        std::vector<std::string> kmers = getAllWords(6);
        //KmerDist<6, std::atomic<uint32_t>> kh;
        for (auto& k : kmers) {
            auto i = indexForKmer(k.c_str(), 6, Direction::FORWARD);
            auto kp = kmerForIndex(i, 6);
            WHEN("kmer is [" + k + "]") {
                THEN("decodes as [" + kp + "]") {
                    REQUIRE(k == kp);
                }
            }
        }
    }
}


SCENARIO("The next k-mer index function works correctly") {
    using salmon::utils::Direction;
    const uint32_t K = 6;
    std::string s = "ATTCTCCACATAGTTGTCATCGAACCAGTACCCCGTAAGCGCCAACATAT";

    GIVEN("The string " + s) {
        auto idx = indexForKmer(s.c_str(), 6, Direction::FORWARD);
        std::string k = s.substr(0, 6);
        WHEN("kmer is [" + k + "]") {
            auto kp = kmerForIndex(idx, 6);
            THEN("decodes as [" + kp + "]") {
                REQUIRE(k == kp);
            }
        }
        for (size_t i = 0; i < s.size() - K; ++i) {
            idx = nextKmerIndex(idx, s[i+K], 6, Direction::FORWARD);
            k = s.substr(i+1, 6);
            WHEN("kmer is [" + k + "]") {
                auto kp = kmerForIndex(idx, 6);
                THEN("decodes as [" + kp + "]") {
                    REQUIRE(k == kp);
                }
            }
        }
    }

    auto rc = [](std::string s) -> std::string {
        std::string rc;
        for (int32_t i = s.size() - 1; i >= 0; --i) {
            switch(s[i]) {
                case 'A': rc += 'T'; break;
                case 'C': rc += 'G'; break;
                case 'G': rc += 'C'; break;
                case 'T': rc += 'A'; break;
            }
        }
        return rc;
    };

    auto rcs = rc(s);

    GIVEN("The string " + s + " in the reverse complement direction") {
        auto idx = indexForKmer(s.c_str() + s.size() - K - 1, 6,
                                Direction::REVERSE_COMPLEMENT);
        std::string k = rc(s.substr(s.size() - K - 1, 6));
        WHEN("kmer is [" + k + "]") {
            auto kp = kmerForIndex(idx, 6);
            THEN("decodes as [" + kp + "]") {
                REQUIRE(k == kp);
            }
        }
        for (int32_t i = s.size() - K - 2; i >= 0; --i) {
            idx = nextKmerIndex(idx, s[i], 6, Direction::REVERSE);
            k = rc(s.substr(i, 6));
            WHEN("kmer is [" + k + "]") {
                auto kp = kmerForIndex(idx, 6);
                THEN("decodes as [" + kp + "]") {
                    REQUIRE(k == kp);
                }
            }
        }
    }

}


