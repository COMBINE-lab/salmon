#include <random>

std::string generateRandomSequence(size_t length, std::uniform_int_distribution<>& dis, std::mt19937& gen) {
  char nucs[] =  {'A', 'C', 'G', 'T'};
  std::string s(length, 'N');
  for (size_t i = 0; i < length; ++i) {
    s[i] = nucs[dis(gen)];
  }
  return s;
}

SCENARIO("GC sampling works properly") {

    GIVEN("A collection of random transcript sequences") {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::mt19937 gen2(rd());
      std::uniform_int_distribution<> dis(0, 3); 
      std::uniform_int_distribution<> dis2(500, 1500); 
      
      char** txpSeqs = new char*[1000];
      char** names = new char*[1000];
      std::vector<size_t> lengths;
      for (size_t tn = 0; tn < 1000; ++tn) {
	auto l = dis2(gen2);
	lengths.emplace_back(l);
	auto s = generateRandomSequence(l, dis, gen);
	txpSeqs[tn] = new char[s.length() + 1];
	names[tn] = new char[3];
	//names[tn] = "HI\0";
	std::strcpy(txpSeqs[tn], s.c_str());
      }
      std::vector<Transcript> txpsSampled;
      std::vector<Transcript> txpsUnSampled;
      for (size_t tn = 0; tn < 1000; ++tn) {
	auto len = lengths[tn];
	txpsSampled.emplace_back(tn, names[tn], len); 
	txpsUnSampled.emplace_back(tn, names[tn], len); 
	txpsSampled[tn].setSequenceBorrowed(txpSeqs[tn], true, 5);
	txpsUnSampled[tn].setSequenceBorrowed(txpSeqs[tn], true, 1);
      }
      
      for (size_t tn = 0; tn < 1000; ++tn) {
	WHEN("Computing GC content") {
	  auto l = txpsSampled[tn].RefLength;
	  for (size_t i = 0; i < l; ++i) {
	    THEN("Sampled is the same as unsampled") {
	      REQUIRE(txpsSampled[tn].gcAt(i) == txpsUnSampled[tn].gcAt(i));
	    }
	  }
	}
      }

      for (size_t tn = 0; tn < 1000; ++tn) {
	delete [] txpSeqs[tn];
	delete [] names[tn];
      }
      delete txpSeqs;
      delete names;
    } // end GIVEN
}
