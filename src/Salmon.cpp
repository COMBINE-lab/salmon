/**
>HEADER
    Copyright (c) 2013 -- 2016 Rob Patro rob.patro@cs.stonybrook.edu

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
#include <boost/any.hpp>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

#include "GenomicFeature.hpp"
#include "SalmonConfig.hpp"
#include "VersionChecker.hpp"

int help(std::vector<std::string> opts) { //}int argc, char* argv[]) {
    fmt::MemoryWriter helpMsg;
    helpMsg.write("Salmon v{}\n\n", salmon::version);
    helpMsg.write("Usage:  salmon -h|--help or \n"
                  "        salmon -v|--version or \n"
		  "        salmon -c|--cite or \n"
                  "        salmon [--no-version-check] <COMMAND> [-h | options]\n\n");
    helpMsg.write("Commands:\n");
    helpMsg.write("     cite  Show salmon citation information\n");
    helpMsg.write("     index Create a salmon index\n");
    helpMsg.write("     quant Quantify a sample\n");
    //helpMsg.write("     quantmerge Merge multiple quantifications into a single file\n");
    helpMsg.write("     swim  Perform super-secret operation\n");

    /*
    auto orighelpmsg = R"(
    ===============

    Please invoke salmon with one of the following commands {index, quant, swim}.
    For more information on the options for these particular methods, use the -h
    flag along with the method name.  For example:

    salmon index -h

    will give you detailed help information about the index command.
    )";
    */

    std::cerr << helpMsg.str();
    return 0;
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
    return 0;
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

/**
 * Bonus!
 */
void printCite() {

  std::cerr << R"(
Reference:
==========

Salmon provides accurate, fast, and bias-aware transcript expression estimates using dual-phase inference
Rob Patro, Geet Duggal, Michael I Love, Rafael A Irizarry, Carl Kingsford
bioRxiv 021592; doi: http://dx.doi.org/10.1101/021592

bibtex:
=======

@article {Patro:2016,
  author = {Patro, Rob and Duggal, Geet and Love, Michael I and Irizarry, Rafael A and Kingsford, Carl},
  title = {Salmon provides accurate, fast, and bias-aware transcript expression estimates using dual-phase inference},
  year = {2016},
  doi = {10.1101/021592},
  publisher = {Cold Spring Harbor Labs Journals},
  URL = {http://biorxiv.org/content/early/2016/08/30/021592},
  journal = {bioRxiv}
}
)";

}


int salmonIndex(int argc, char* argv[]);
int salmonQuantify(int argc, char* argv[]);
int salmonAlignmentQuantify(int argc, char* argv[]);

bool verbose = false;

int main( int argc, char* argv[] ) {
  using std::string;
  namespace po = boost::program_options;

  // With no arguments, print help
  if (argc == 1) {
      std::vector<std::string> o;
      help(o);//argc, argv);
      std::exit(1);
  }

  try {
      
    // subcommand parsing code inspired by : https://gist.github.com/randomphrase/10801888
    po::options_description sfopts("Allowed Options");
    sfopts.add_options()
        ("version,v", "print version string")
        ("no-version-check", "don't check with the server to see if this is the latest version")
        ("cite,c", "show citation information")
        ("help,h", "produce help message")
        ("command", po::value<string>(), "command to run {index, quant, sf}")
        ("subargs", po::value<std::vector<std::string>>(), "Arguments for command");
    ;

    po::options_description all("Allowed Options");
    all.add(sfopts);

    po::positional_options_description pd;
    pd.add("command", 1).add("subargs", -1);

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv).options(all).positional(pd).allow_unregistered().run();
    po::store(parsed, vm);

    if (vm.count("version")) {
      std::cerr << "version : " << salmon::version << "\n";
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

    if (!vm.count("no-version-check")){
      std::string versionMessage = getVersionMessage();
      std::cerr << versionMessage;
    }
    
    //po::notify(vm);

    std::string cmd = vm["command"].as<std::string>();
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());
    // if there was a help and a command, then add the help back since it was parsed
    if (vm.count("help")) { opts.insert(opts.begin(), "--help"); }

    std::unordered_map<string, std::function<int(int, char*[])>> cmds({
      {"index", salmonIndex},
      {"quant", salmonQuantify},
      //{"quantmerge", salmonQuantMerge},
      {"swim", salmonSwim}
    });

    /*
    //string cmd = vm["command"].as<string>();
    int subCommandArgc = argc - topLevelArgc + 1;
    char** argv2 = new char*[subCommandArgc];
    argv2[0] = argv[0];
    std::copy_n( &argv[topLevelArgc], argc-topLevelArgc, &argv2[1] );
    */

    int subCommandArgc = opts.size() + 1;
    std::unique_ptr<char*[]> argv2(new char*[subCommandArgc]);
    argv2[0] = argv[0];
    for (size_t i = 0; i < subCommandArgc - 1; ++i) {
        argv2[i+1] = &*opts[i].begin();
    }
    
    auto cmdMain = cmds.find(cmd);
    if (cmdMain == cmds.end()) {
        //help(subCommandArgc, argv2);
        return help(opts);
    } else {
      // If the command is quant; determine whether
      // we're quantifying with raw sequences or alignemnts
      if (cmdMain->first == "quant") {

        // detect mode-specific help request
        if (strncmp(argv2[1], "--help-alignment", 16) == 0) {
            std::vector<char> helpStr{'-','-','h','e','l','p','\0'};
            char* helpArgv[] = {argv[0], &helpStr[0]};
            return salmonAlignmentQuantify(2, helpArgv);
        } else if (strncmp(argv2[1], "--help-reads", 12) == 0) {
            std::vector<char> helpStr{'-','-','h','e','l','p','\0'};
            char* helpArgv[] = {argv[0], &helpStr[0]};
            return salmonQuantify(2, helpArgv);
        }

        // detect general help request
        if (strncmp(argv2[1], "--help", 6) == 0 or
            strncmp(argv2[1], "-h", 2) == 0) {
            return dualModeMessage();
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
            return salmonAlignmentQuantify(subCommandArgc, argv2.get());
        } else {
            return salmonQuantify(subCommandArgc, argv2.get());
        }
      } else {
        return cmdMain->second(subCommandArgc, argv2.get());
      }
    }

  } catch (po::error &e) {
    std::cerr << "Program Option Error (main) : [" << e.what() << "].\n Exiting.\n";
    std::exit(1);
  } catch (...) {
    std::cerr << argv[0] << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0] << " --help\nExiting.\n";
  }

  return 0;
}
