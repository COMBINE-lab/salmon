#include "salmon/internal/util/LibraryTypeUtils.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/range/join.hpp>

#include <spdlog/spdlog.h>

namespace salmon::utils {

LibraryFormat parseLibraryFormatStringNew(std::string& fmt) {
  using std::map;
  using std::string;
  using std::stringstream;

  map<string, LibraryFormat> formatMap = {
      {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                           ReadStrandedness::U)},
      {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                            ReadStrandedness::SA)},
      {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                            ReadStrandedness::AS)},
      {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                           ReadStrandedness::U)},
      {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                            ReadStrandedness::SA)},
      {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY,
                            ReadStrandedness::AS)},
      {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                           ReadStrandedness::U)},
      {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                            ReadStrandedness::S)},
      {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME,
                            ReadStrandedness::A)},
      {"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                          ReadStrandedness::U)},
      {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                           ReadStrandedness::S)},
      {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                           ReadStrandedness::A)}};

  // inspired by
  // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
  // first convert the string to upper-case
  for (auto& c : fmt) {
    c = std::toupper(c);
  }

  auto libFmtIt = formatMap.find(fmt);

  if (libFmtIt == formatMap.end()) {
    stringstream errstr;
    errstr << "unknown library format string : " << fmt;
    throw std::invalid_argument(errstr.str());
  }

  return libFmtIt->second;
}

std::vector<ReadLibrary>
extractReadLibraries(boost::program_options::parsed_options& orderedOptions) {
  // The current (default) format for paired end data
  LibraryFormat peFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                         ReadStrandedness::U);
  // The current (default) format for single end data
  LibraryFormat seFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                         ReadStrandedness::U);

  auto isAutoLibType = [](std::string& fmt) -> bool {
    return (fmt.length() == 1 and (fmt.front() == 'a' or fmt.front() == 'A'));
  };

  auto log = spdlog::get("jointLog");

  bool sawLibType{false};
  bool sawPairedLibrary{false};
  bool sawUnpairedLibrary{false};
  bool autoLibType{false};

  std::vector<ReadLibrary> peLibs{peFormat};
  std::vector<ReadLibrary> seLibs{seFormat};

  for (auto& opt : orderedOptions.options) {
    // Update the library type
    if (opt.string_key == "libType") {
      if (!isAutoLibType(opt.value[0])) {
        auto libFmt = parseLibraryFormatStringNew(opt.value[0]);
        if (libFmt.type == ReadType::PAIRED_END) {
          peFormat = libFmt;
          peLibs.emplace_back(libFmt);
        } else {
          seFormat = libFmt;
          seLibs.emplace_back(libFmt);
        }
      } else {
        autoLibType = true;
      }
      sawLibType = true;
    }

    if (opt.string_key == "mates1") {
      if (!sawLibType) {
        log->warn("Encountered a read file (--mates1/-1) before a "
                  "library type specification.  The (--libType/-l) "
                  "option must precede the input files.");
        peLibs.clear();
        return peLibs;
      }
      peLibs.back().addMates1(opt.value);
      if (autoLibType) {
        peLibs.back().enableAutodetect();
      }
      sawPairedLibrary = true;
    }
    if (opt.string_key == "mates2") {
      if (!sawLibType) {
        log->warn("Encountered a read file (--mates2/-2) before a "
                  "library type specification.  The (--libType/-l) "
                  "option must precede the input files.");
        peLibs.clear();
        return peLibs;
      }
      peLibs.back().addMates2(opt.value);
      if (autoLibType) {
        peLibs.back().enableAutodetect();
      }
      sawPairedLibrary = true;
    }
    if (opt.string_key == "unmatedReads") {
      if (!sawLibType) {
        log->warn("Encountered a read file (--unmatedReads/-r) before a "
                  "library type specification.  The (--libType/-l) "
                  "option must precede the input files.");
        seLibs.clear();
        return seLibs;
      }
      seLibs.back().addUnmated(opt.value);
      if (autoLibType) {
        seLibs.back().enableAutodetect();
      }
      sawUnpairedLibrary = true;
    }
  }

  std::vector<ReadLibrary> libs;

  // Allow this temporarily for now, but the mixed-library parsing strategy
  // should eventually be cleaned up instead of relying on this fallback.
  /*
  if (sawPairedLibrary and sawUnpairedLibrary) {
    log->warn("It seems you have specified both paired-end and unpaired read "
              "libraries.  Salmon does not accepted mixed library types, and "
              "different library types should typically not be quantified together "
              "anyway.  Please quantifiy distinct library types separately.");
    return libs;
  }
  */
  (void)sawPairedLibrary;
  (void)sawUnpairedLibrary;

  libs.reserve(peLibs.size() + seLibs.size());
  for (auto& lib : boost::range::join(seLibs, peLibs)) {
    if (lib.format().type == ReadType::SINGLE_END) {
      if (lib.unmated().size() == 0) {
        // Didn't use default single end library type
        continue;
      }
    } else if (lib.format().type == ReadType::PAIRED_END) {
      if (lib.mates1().size() == 0 or lib.mates2().size() == 0) {
        // Didn't use default paired-end library type
        continue;
      }
    }
    libs.push_back(lib);
  }

  size_t numLibs = libs.size();
  if (numLibs == 1) {
    log->info("There is 1 library.");
  } else if (numLibs > 1) {
    log->info("There are {} libraries.", numLibs);
  }
  return libs;
}

LibraryFormat parseLibraryFormatString(std::string& fmt) {
  using std::map;
  using std::string;
  using std::stringstream;

  // inspired by
  // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c

  // first convert the string to upper-case
  for (auto& c : fmt) {
    c = std::toupper(c);
  }
  // split on the delimiter ':', and put the key, value (k=v) pairs into a map
  stringstream ss(fmt);
  string item;
  map<string, string> kvmap;
  while (std::getline(ss, item, ':')) {
    auto splitPos = item.find('=', 0);
    string key{item.substr(0, splitPos)};
    string value{item.substr(splitPos + 1)};
    kvmap[key] = value;
  }

  map<string, ReadOrientation> orientationType = {
      {">>", ReadOrientation::SAME},
      {"<>", ReadOrientation::AWAY},
      {"><", ReadOrientation::TOWARD},
      {"*", ReadOrientation::NONE}};
  map<string, ReadStrandedness> strandType = {{"SA", ReadStrandedness::SA},
                                              {"AS", ReadStrandedness::AS},
                                              {"A", ReadStrandedness::A},
                                              {"S", ReadStrandedness::S},
                                              {"U", ReadStrandedness::U}};
  auto it = kvmap.find("T");
  string typeStr = "";
  if (it != kvmap.end()) {
    typeStr = it->second;
  } else {
    it = kvmap.find("TYPE");
    if (it != kvmap.end()) {
      typeStr = it->second;
    }
  }

  if (typeStr != "SE" and typeStr != "PE") {
    string e = typeStr + " is not a valid read type; must be one of {SE, PE}";
    throw std::invalid_argument(e);
  }

  ReadType type =
      (typeStr == "SE") ? ReadType::SINGLE_END : ReadType::PAIRED_END;
  ReadOrientation orientation = (type == ReadType::SINGLE_END)
                                    ? ReadOrientation::NONE
                                    : ReadOrientation::TOWARD;
  ReadStrandedness strandedness{ReadStrandedness::U};
  // Construct the LibraryFormat class from the key, value map
  for (auto& kv : kvmap) {
    auto& k = kv.first;
    auto& v = kv.second;
    if (k == "O" or k == "ORIENTATION") {
      auto orientIt = orientationType.find(v);
      if (orientIt != orientationType.end()) {
        orientation = orientationType[orientIt->first];
      } else {
        string e =
            v + " is not a valid orientation type; must be one of {>>, <>, ><}";
        throw std::invalid_argument(e);
      }
    }
    if (k == "S" or k == "STRAND") {
      auto strandIt = strandType.find(v);
      if (strandIt != strandType.end()) {
        strandedness = strandType[strandIt->first];
      } else {
        string e =
            v + " is not a valid strand type; must be one of {SA, AS, S, A, U}";
        throw std::invalid_argument(e);
      }
    }
  }
  return LibraryFormat(type, orientation, strandedness);
}

} // namespace salmon::utils
