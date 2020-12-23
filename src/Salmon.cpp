/**
>HEADER
    Copyright (c) 2013 -- 2017 Rob Patro rob.patro@cs.stonybrook.edu

    This file is part of salmon.

    salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#include <boost/thread/thread.hpp>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

#include "GenomicFeature.hpp"
#include "SalmonConfig.hpp"
#include "VersionChecker.hpp"

int help(const std::vector<std::string>& /*opts*/) { 
  fmt::MemoryWriter helpMsg;
  helpMsg.write("salmon v{}\n\n", salmon::version);
  helpMsg.write(
      "Usage:  salmon -h|--help or \n"
      "        salmon -v|--version or \n"
      "        salmon -c|--cite or \n"
      "        salmon [--no-version-check] <COMMAND> [-h | options]\n\n");
  helpMsg.write("Commands:\n");
  helpMsg.write("     index      : create a salmon index\n");
  helpMsg.write("     quant      : quantify a sample\n");
  helpMsg.write("     alevin     : single cell analysis\n");
  helpMsg.write("     swim       : perform super-secret operation\n");
  helpMsg.write("     quantmerge : merge multiple quantifications into a single file\n");

  std::cout << helpMsg.str();
  return 0;
}

int dualModeMessage() {
  auto helpmsg = R"(
    ===============

    salmon quant has two modes --- one quantifies expression using raw reads
    and the other makes use of already-aligned reads (in BAM/SAM format).
    Which algorithm is used depends on the arguments passed to salmon quant.
    If you provide salmon with alignments '-a [ --alignments ]' then the
    alignment-based algorithm will be used, otherwise the algorithm for
    quantifying from raw reads will be used.

    to view the help for salmon's selective-alignment-based mode, use the command

    salmon quant --help-reads

    To view the help for salmon's alignment-based mode, use the command

    salmon quant --help-alignment

    )";
  std::cout << "    salmon v" << salmon::version << helpmsg << "\n";
  return 0;
}

/**
 * Bonus!
 */
int salmonSwim(int /*argc*/, const char* /*argv*/[]) {

  std::cout << R"(
    _____       __
   / ___/____ _/ /___ ___  ____  ____
   \__ \/ __ `/ / __ `__ \/ __ \/ __ \
  ___/ / /_/ / / / / / / / /_/ / / / /
 /____/\__,_/_/_/ /_/ /_/\____/_/ /_/


)";

  return 0;
}

/**
 * Citation
 */
void printCite() {

  std::cout << R"(
If you use salmon in your research, please cite the publication in any
papers, pre-prints or reports.  The proper citation information for salmon
appears below.

Reference:
==========

Rob Patro, Geet Duggal, Michael I. Love, Rafael A. Irizarry, Carl Kingsford.
Salmon provides fast and bias-aware quantification of transcript expression.
Nature Methods. 2017;14(4):417-419. doi: 10.1038/nmeth.4197

bibtex:
=======

@article{Patro2017Salmon,
  doi = {10.1038/nmeth.4197},
  url = {https://doi.org/10.1038%2Fnmeth.4197},
  year  = {2017},
  month = {mar},
  publisher = {{Springer Nature}},
  volume = {14},
  number = {4},
  pages = {417--419},
  author = {Rob Patro and Geet Duggal and Michael I Love and Rafael A Irizarry and Carl Kingsford},
  title = {Salmon provides fast and bias-aware quantification of transcript expression},
  journal = {{Nature Methods}}
}
)";
}

int salmonIndex(int argc, const char* argv[]);
int salmonQuantify(int argc, const char* argv[]);
int salmonAlignmentQuantify(int argc, const char* argv[]);
// TODO : PF_INTEGRATION
int salmonBarcoding(int argc, const char* argv[]);
int salmonQuantMerge(int argc, const char* argv[]);

bool verbose = false;

int main(int argc, char* argv[]) {
  using std::string;
  namespace po = boost::program_options;
  std::setlocale(LC_ALL, "en_US.UTF-8");

  // With no arguments, print help
  if (argc == 1) {
    std::vector<std::string> o;
    help(o); // argc, argv);
    std::exit(1);
  }

  try {

    // subcommand parsing code inspired by :
    // https://gist.github.com/randomphrase/10801888
    po::options_description sfopts("Allowed Options");
    sfopts.add_options()("version,v", "print version string")(
        "no-version-check",
        "don't check with the server to see if this is the latest version")(
        "cite,c", "show citation information")(
        "help,h", "produce help message")("command", po::value<string>(),
                                          "command to run {index, quant, sf}")(
        "subargs", po::value<std::vector<std::string>>(),
        "Arguments for command");

    po::options_description all("Allowed Options");
    all.add(sfopts);

    po::positional_options_description pd;
    pd.add("command", 1).add("subargs", -1);

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv)
                                    .options(all)
                                    .positional(pd)
                                    .allow_unregistered()
                                    .run();
    po::store(parsed, vm);

    if (vm.count("version")) {
      std::cout << "salmon " << salmon::version << "\n";
      std::exit(0);
    }

    if (vm.count("help") and !vm.count("command")) {
      std::vector<std::string> o;
      help(o);
      std::exit(0);
    }

    if (vm.count("cite") and !vm.count("command")) {
      printCite();
      std::exit(0);
    }

    const char* no_version_env_ptr = std::getenv("SALMON_NO_VERSION_CHECK");
    std::string no_version_env = (no_version_env_ptr == nullptr) ? "" : std::string(no_version_env_ptr);
    std::transform(no_version_env.begin(), no_version_env.end(), no_version_env.begin(), 
                   [](unsigned char c){ return std::toupper(c); } // correct
                  );
    bool skip_version_check = vm.count("no-version-check") or (no_version_env == "1") 
                              or (no_version_env == "TRUE") or (no_version_env == "T");

    if (!skip_version_check) {
      std::string versionMessage = getVersionMessage();
      std::cerr << versionMessage;
    }

    // po::notify(vm);

    std::string cmd = vm["command"].as<std::string>();
    std::vector<std::string> opts =
        po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());
    // if there was a help and a command, then add the help back since it was
    // parsed
    if (vm.count("help")) {
      opts.insert(opts.begin(), "--help");
    }

    std::unordered_map<string, std::function<int(int, const char* [])>> cmds(
        {{"index", salmonIndex},
         {"quant", salmonQuantify},
         {"quantmerge", salmonQuantMerge},
         // TODO : PF_INTEGRATION
         {"alevin", salmonBarcoding},
         {"swim", salmonSwim}});

    /*
    //string cmd = vm["command"].as<string>();
    int subCommandArgc = argc - topLevelArgc + 1;
    char** argv2 = new char*[subCommandArgc];
    argv2[0] = argv[0];
    std::copy_n( &argv[topLevelArgc], argc-topLevelArgc, &argv2[1] );
    */

    int32_t subCommandArgc = opts.size() + 1;
    std::unique_ptr<const char* []> argv2(new const char*[subCommandArgc]);
    argv2[0] = argv[0];
    for (int32_t i = 0; i < subCommandArgc - 1; ++i) {
      argv2[i + 1] = opts[i].c_str();
    }

    auto cmdMain = cmds.find(cmd);
    if (cmdMain == cmds.end()) {
      // help(subCommandArgc, argv2);
      return help(opts);
    } else {
      // If the command is quant; determine whether
      // we're quantifying with raw sequences or alignments
      if (cmdMain->first == "quant") {

        if (subCommandArgc < 2) {
          return dualModeMessage();
        }
        // detect mode-specific help request
        if (strncmp(argv2[1], "--help-alignment", 16) == 0) {
          std::vector<char> helpStr{'-', '-', 'h', 'e', 'l', 'p', '\0'};
          const char* helpArgv[] = {argv[0], &helpStr[0]};
          return salmonAlignmentQuantify(2, helpArgv);
        } else if (strncmp(argv2[1], "--help-reads", 12) == 0) {
          std::vector<char> helpStr{'-', '-', 'h', 'e', 'l', 'p', '\0'};
          const char* helpArgv[] = {argv[0], &helpStr[0]};
          return salmonQuantify(2, helpArgv);
        }

        // detect general help request
        if (strncmp(argv2[1], "--help", 6) == 0 or
            strncmp(argv2[1], "-h", 2) == 0) {
          return dualModeMessage();
        }

        // otherwise, detect and dispatch the correct mode
        bool useSalmonAlign{false};
        for (int32_t i = 0; i < subCommandArgc; ++i) {
          if (strncmp(argv2[i], "-a", 2) == 0 or
              strncmp(argv2[i], "-e", 2) == 0 or
              strncmp(argv2[i], "--alignments", 12) == 0 or
              strncmp(argv2[i], "--eqclasses", 11) == 0) {
            useSalmonAlign = true;
            break;
          }
        }
        if (useSalmonAlign) {
          return salmonAlignmentQuantify(subCommandArgc, argv2.get());
        } else {
          return salmonQuantify(subCommandArgc, argv2.get());
        }
      } else {
        return cmdMain->second(subCommandArgc, argv2.get());
      }
    }

  } catch (po::error& e) {
    std::cerr << "Program Option Error (main) : [" << e.what()
              << "].\n Exiting.\n";
    std::exit(1);
  } catch (...) {
    std::cerr << argv[0] << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0]
              << " --help\nExiting.\n";
  }

  return 0;
}