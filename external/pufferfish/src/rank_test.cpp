#include "rank9b.hpp"
#include "core/range.hpp"
#include "rank9sel.hpp"
#include "compact_vector/compact_vector.hpp"


int main(int argc, char* argv[]) {

  size_t nkeys{30};

    // Note: the compact_vector constructor does not
    // init mem to 0, so we do that with the clear_mem() function.
    compact::vector<uint64_t, 1> cvec(nkeys);
    cvec.clear_mem();

    cvec[0] = 1;
    cvec[5] = 1;
    cvec[13] = 1;
    cvec[24] = 1;
    cvec[25] = 1;

    rank9b crank(cvec.get(), cvec.size());

    for (size_t i = 0; i < nkeys; ++i) {
      std::cerr << "rank cvec[" << i << "] = " << crank.rank(i) << "\n";
    }
    return 0;
}
