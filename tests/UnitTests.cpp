#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <unordered_map>
#include <iostream>
#include "catch.hpp"
#include "LibraryFormat.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"

bool verbose=false; // Apparently, we *need* this (OSX)
#include "GCSampleTests.cpp"
#include "LibraryTypeTests.cpp"
//#include "KmerHistTests.cpp"

