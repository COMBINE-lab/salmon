/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#include <boost/thread/thread.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <string>
#include <memory>
#include <functional>
#include <unordered_map>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "BiasIndex.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "SalmonConfig.hpp"
#include "VersionChecker.hpp"

int help(int argc, char* argv[]) {
    auto helpmsg = R"(
    ===============

    Please invoke salmon with one of the following commands {index, quant, swim}.
    For more information on the options for these particular methods, use the -h
    flag along with the method name.  For example:

    salmon index -h

    will give you detailed help information about the index command.
    )";

    std::cerr << "    Salmon v" << salmon::version << helpmsg << "\n";
    return 1;
}


int dualModeMessage() {
    auto helpmsg = R"(
    ===============

    salmon quant has two modes --- one quantifies expression using raw reads
    and the other makes use of already-aligned reads (in BAM/SAM format).
    which algorithm is used depends on the arguments passed to salmon quant.
    If you provide salmon with alignments '-a [ --alignments ]' then the
    alignment-based algorithm will be used, otherwise the algorithm for
    quantifying from raw reads will be used.

    To view the help for salmon's alignment-based mode, use the command

    salmon quant --help-alignment

    to view the help for salmon's read-based mode, use the command

    salmon quant --help-reads
    )";
    std::cerr << "    Salmon v" << salmon::version << helpmsg << "\n";
    return 1;
}


/**
 * Bonus!
 */
int salmonSwim(int argc, char* argv[]) {

  std::cerr << R"(
    _____       __
   / ___/____ _/ /___ ___  ____  ____
   \__ \/ __ `/ / __ `__ \/ __ \/ __ \
  ___/ / /_/ / / / / / / / /_/ / / / /
 /____/\__,_/_/_/ /_/ /_/\____/_/ /_/


)";

  return 0;

}

int salmonIndex(int argc, char* argv[]);
int salmonQuantify(int argc, char* argv[]);
int salmonAlignmentQuantify(int argc, char* argv[]);

bool verbose = false;

int main( int argc, char* argv[] ) {
  using std::string;
  namespace po = boost::program_options;

  try {

    po::options_description hidden("hidden");
    hidden.add_options()
    ("command", po::value<string>(), "command to run {index, quant, sf}");

    po::options_description sfopts("Allowed Options");
    sfopts.add_options()
    ("version,v", "print version string")
    ("no-version-check", "don't check with the server to see if this is the latest version")
    ("help,h", "produce help message")
    ;

    po::options_description all("Allowed Options");
    all.add(sfopts).add(hidden);

    po::positional_options_description pd;
    pd.add("command", 1);

    size_t topLevelArgc = argc;
    for (size_t i : boost::irange(size_t{1}, static_cast<size_t>(argc))) {
      if (argv[i][0] != '-') {
        topLevelArgc = i+1;
        break;
      }
    }

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(topLevelArgc, argv).options(all).positional(pd).allow_unregistered().run();
    po::store(parsed, vm);

    if (vm.count("version")) {
      std::cerr << "version : " << salmon::version << "\n";
      std::exit(0);
    }

    if (vm.count("help") and !vm.count("command")) {
        std::cout << sfopts << std::endl;
        help(argc, argv);
        std::exit(0);
    }

    if (!vm.count("no-version-check")){
      std::string versionMessage = getVersionMessage();
      std::cerr << versionMessage;
    }

    po::notify(vm);

    std::unordered_map<string, std::function<int(int, char*[])>> cmds({
      {"index", salmonIndex},
      {"quant", salmonQuantify},
      {"swim", salmonSwim}
    });

    string cmd = vm["command"].as<string>();

    int subCommandArgc = argc - topLevelArgc + 1;
    char** argv2 = new char*[subCommandArgc];
    argv2[0] = argv[0];
    std::copy_n( &argv[topLevelArgc], argc-topLevelArgc, &argv2[1] );

    auto cmdMain = cmds.find(cmd);
    if (cmdMain == cmds.end()) {
      help(subCommandArgc, argv2);
    } else {
      // If the command is quant; determine whether
      // we're quantifying with raw sequences or alignemnts
      if (cmdMain->first == "quant") {

        // detect mode-specific help request
        if (strncmp(argv2[1], "--help-alignment", 16) == 0) {
            std::vector<char> helpStr{'-','-','h','e','l','p','\0'};
            char* helpArgv[] = {argv[0], &helpStr[0]};
            salmonAlignmentQuantify(2, helpArgv);
        } else if (strncmp(argv2[1], "--help-reads", 12) == 0) {
            std::vector<char> helpStr{'-','-','h','e','l','p','\0'};
            char* helpArgv[] = {argv[0], &helpStr[0]};
            salmonQuantify(2, helpArgv);
        }

        // detect general help request
        if (strncmp(argv2[1], "--help", 6) == 0 or
            strncmp(argv2[1], "-h", 2) == 0) {
            dualModeMessage();
            std::exit(0);
        }

        // otherwise, detect and dispatch the correct mode
        bool useSalmonAlign{false};
        for (size_t i = 0; i < subCommandArgc; ++i) {
            if (strncmp(argv2[i], "-a", 2) == 0 or
                strncmp(argv2[i], "--alignments", 12) == 0) {
                useSalmonAlign = true;
                break;
            }
        }
        if (useSalmonAlign) {
            salmonAlignmentQuantify(subCommandArgc, argv2);
        } else {
            salmonQuantify(subCommandArgc, argv2);
        }
      } else {
        cmdMain->second(subCommandArgc, argv2);
      }
    }
    delete[] argv2;

  } catch (po::error &e) {
    std::cerr << "Program Option Error (main) : [" << e.what() << "].\n Exiting.\n";
    std::exit(1);
  } catch (...) {
    std::cerr << argv[0] << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0] << " --help\nExiting.\n";
  }

  return 0;
}
