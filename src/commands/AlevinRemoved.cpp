#include <iostream>
#include <memory>

#include "salmon/internal/index/SalmonIndex.hpp"

int salmonBarcoding(int /*argc*/, const char* /*argv*/[],
                    std::unique_ptr<SalmonIndex>& /*index*/) {
  std::cerr
      << "The `salmon alevin` command has been removed from this bulk-only "
         "modernized Salmon.\n"
      << "Use alevin-fry for single-cell analysis:\n"
      << "  https://alevin-fry.readthedocs.io/en/latest/\n"
      << "If you need the legacy `alevin` implementation, use the last "
         "pre-removal Salmon release.\n";
  return 1;
}
