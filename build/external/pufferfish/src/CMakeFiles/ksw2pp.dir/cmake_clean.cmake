file(REMOVE_RECURSE
  "libksw2pp.a"
  "libksw2pp.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/ksw2pp.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
