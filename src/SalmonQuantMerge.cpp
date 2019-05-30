/**
>HEADER
    Copyright (c) 2013, 2014, 2015, 2016 Rob Patro rob.patro@cs.stonybrook.edu

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

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <memory>
// C++ string formatting library #include "spdlog/fmt/fmt.h"
// logger includes
#include "spdlog/spdlog.h"

enum class TargetColumn { LEN, ELEN, TPM, NREADS };

class QuantMergeOptions {
public:
  std::vector<std::string> samples;
  std::vector<std::string> names;
  std::string outputName;
  bool genesQuant;
  std::string missingValue;
  std::string outputCol;
  std::shared_ptr<spdlog::logger> log;
  TargetColumn tcol;

  void print() {
    std::string slist;
    for (size_t n = 0; n < samples.size(); ++n) {
      if (n > 0) {
        slist += ", ";
      }
      slist += samples[n];
    }
    log->info("samples: [ {} ]", slist);
    if (names.size() > 0) {
      std::string nlist;
      for (size_t n = 0; n < names.size(); ++n) {
        if (n > 0) {
          nlist += ", ";
        }
        nlist += names[n];
      }
      log->info("sample names : [ {} ]", nlist);
    }
    log->info("output column : {}", outputCol);
    log->info("output file : {}", outputName);
  }
};

bool validateOptions(boost::program_options::variables_map& vm,
                     QuantMergeOptions& qmOpts,
                     std::shared_ptr<spdlog::logger> log) {
  // check that if we have a list of names, it's the same length
  // as the list of samples
  if (vm.count("names")) {
    if (qmOpts.names.size() != qmOpts.samples.size()) {
      log->critical("If you provide names for your samples, they number of "
                    "names must be the "
                    "same as the number of samples.  You provided {} samples "
                    "but {} names",
                    qmOpts.samples.size(), qmOpts.names.size());
      std::exit(1);
    }
  } else {
    for (auto& sp : qmOpts.samples) {
      qmOpts.names.push_back(boost::filesystem::path(sp).filename().string());
    }
  }

  std::transform(qmOpts.outputCol.begin(), qmOpts.outputCol.end(),
                 qmOpts.outputCol.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  if (qmOpts.outputCol == "LEN" or qmOpts.outputCol == "LENGTH") {
    qmOpts.tcol = TargetColumn::LEN;
  } else if (qmOpts.outputCol == "ELEN" or qmOpts.outputCol == "ELENGTH" or
             qmOpts.outputCol == "EFFLEN" or qmOpts.outputCol == "EFFLENGTH" or
             qmOpts.outputCol == "EFFECTIVELEN" or
             qmOpts.outputCol == "EFFECTIVELENGTH") {
    qmOpts.tcol = TargetColumn::ELEN;
  } else if (qmOpts.outputCol == "TPM") {
    qmOpts.tcol = TargetColumn::TPM;
  } else if (qmOpts.outputCol == "NUMREADS" or
             qmOpts.outputCol == "NUM_READS" or qmOpts.outputCol == "NREADS") {
    qmOpts.tcol = TargetColumn::NREADS;
  } else {
    qmOpts.log->critical(
        "I do not understand the requested output column name, {}. "
        "The output column should be one of {{len, elen, tmp, numreads}}.",
        qmOpts.outputCol);
    std::exit(1);
  }

  return true;
}

template <typename LenT>
struct ExpressionRecord {
  uint32_t snum;
  LenT len;
  double effectiveLen;
  double tpm;
  double numReads;
};

template <typename LenT>
bool doMerge(QuantMergeOptions& qmOpts) {
  std::unordered_map<std::string, std::vector<ExpressionRecord<LenT>>> recs;
  for (uint32_t n = 0; n < qmOpts.samples.size(); ++n) {
    auto& sampDir = qmOpts.samples[n];
    auto quantFile = boost::filesystem::path(sampDir) /
      (qmOpts.genesQuant?"quant.genes.sf":"quant.sf");
    if (boost::filesystem::exists(quantFile) and
        boost::filesystem::is_regular_file(quantFile)) {
      qmOpts.log->info("Parsing {}", quantFile.string());
      std::ifstream ifile(quantFile.string());
      std::string header;
      std::getline(ifile, header);
      std::string targetName;
      LenT len;
      double effectiveLen, tpm, numReads;
      while (ifile >> targetName >> len >> effectiveLen >> tpm >> numReads) {
        recs[targetName].push_back({n, len, effectiveLen, tpm, numReads});
      }
      ifile.close();
    } else {
      qmOpts.log->critical("The sample directory {} either doesn't exist, "
                           "or doesn't contain a quant.sf file",
                           sampDir);
      return false;
    }
  }

  auto outputPath =
      boost::filesystem::absolute(boost::filesystem::path(qmOpts.outputName))
          .parent_path();
  if (!boost::filesystem::exists(outputPath)) {
    if (!boost::filesystem::create_directories(outputPath)) {
      qmOpts.log->critical("Couldn't create output path {}",
                           outputPath.string());
      std::exit(1);
    }
  }

  size_t missingValues{0};
  // Now, the path exists
  std::ofstream outFile(qmOpts.outputName);
  if (!outFile.is_open()) {
    qmOpts.log->critical("Couldn't create output file {}", qmOpts.outputName);
    outFile.close();
    std::exit(1);
  }

  outFile << "Name";
  for (size_t n = 0; n < qmOpts.samples.size(); ++n) {
    outFile << '\t' << qmOpts.names[n];
  }
  outFile << '\n';

  for (auto& kv : recs) {
    auto& trecs = kv.second;
    outFile << kv.first;
    uint32_t nextTrec{0};
    for (uint32_t n = 0; n < qmOpts.samples.size(); ++n) {
      if (nextTrec < trecs.size() and trecs[nextTrec].snum == n) {
        switch (qmOpts.tcol) {
        case TargetColumn::LEN:
          outFile << '\t' << trecs[nextTrec].len;
          break;
        case TargetColumn::ELEN:
          outFile << '\t' << trecs[nextTrec].effectiveLen;
          break;
        case TargetColumn::TPM:
          outFile << '\t' << trecs[nextTrec].tpm;
          break;
        case TargetColumn::NREADS:
          outFile << '\t' << trecs[nextTrec].numReads;
          break;
        }
        ++nextTrec;
      } else {
        ++missingValues;
        outFile << '\t' << qmOpts.missingValue;
      }
    }
    outFile << '\n';
  }
  outFile.close();

  if (missingValues > 0) {
    qmOpts.log->warn(
        "There were {} missing entries (recorded as \"{}\") in the output",
        missingValues, qmOpts.missingValue);
  }

  return true;
}

int salmonQuantMerge(int argc, const char* argv[]) {
  using std::cerr;
  using std::vector;
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  QuantMergeOptions qmOpts;
  po::options_description generic("\n"
                                  "basic options");
  generic.add_options()("version,v", "print version string")
    ("help,h", "produce help message")
    ("quants",
      po::value<vector<string>>(&qmOpts.samples)->multitoken()->required(),
     "List of quantification directories.")
    ("names", po::value<vector<string>>(&qmOpts.names)->multitoken(),
     "Optional list of names to give to the samples.")
    ("column,c",
     po::value<string>(&qmOpts.outputCol)->required()->default_value("TPM"),
     "The name of the column that will be merged together into the output "
     "files. "
     "The options are {len, elen, tpm, numreads}")
    ("genes", po::bool_switch(&(qmOpts.genesQuant))->default_value(false),
     "Use gene quantification instead of transcript.")
    ("missing",
     po::value<string>(&qmOpts.missingValue)->required()->default_value("NA"),
     "The value of missing values.")
    ("output,o", po::value<std::string>(&qmOpts.outputName)->required(),
     "Output quantification file.");

  po::options_description all("salmon quantmerge options");
  all.add(generic);

  po::options_description visible("salmon quantmerge options");
  visible.add(generic);

  po::variables_map vm;
  try {
    auto orderedOptions =
        po::command_line_parser(argc, argv).options(all).run();

    po::store(orderedOptions, vm);

    if (vm.count("help")) {
      auto hstring = R"(
quantmerge
==========
Merge multiple quantification results into
a single file.
)";
      std::cerr << hstring << std::endl;
      std::cerr << visible << std::endl;
      std::exit(0);
    }

    po::notify(vm);
    size_t max_q_size = 131072;
    spdlog::set_async_mode(max_q_size);

    auto consoleSink =
        std::make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>();
    auto consoleLog = spdlog::create("mergeLog", {consoleSink});
    qmOpts.log = consoleLog;

    validateOptions(vm, qmOpts, consoleLog);
    qmOpts.print();


    using gene_len_t = double;
    using transcript_len_t = uint32_t;

    if (qmOpts.genesQuant) {
      doMerge<gene_len_t>(qmOpts);
    } else {
      doMerge<transcript_len_t>(qmOpts);
    }
  } catch (po::error& e) {
    std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (const spdlog::spdlog_ex& ex) {
    std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n";
    std::cerr << argv[0] << " quant was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0]
              << " quant --help\nExiting.\n";
    std::exit(1);
  }

  return 0;
}
