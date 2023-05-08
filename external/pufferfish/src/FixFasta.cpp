#include "FastxParser.hpp"
#include "jellyfish/mer_dna.hpp"
#include "clipp.h"
#include "sparsepp/spp.h"
#include "spdlog/spdlog.h"
#include "xxhash.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <unordered_set>

#include "cereal/cereal.hpp"
#include "cereal/archives/json.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
#include "string_view.hpp"
#include "digestpp/digestpp.hpp"
#include "ghc/filesystem.hpp"

using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

bool fixFasta(single_parser* parser,
              // std::string& outputDir,
              spp::sparse_hash_set<std::string>& decoyNames,
              bool keepDuplicates, uint32_t k,
              std::string& sepStr, bool expect_transcriptome, 
              bool noclip_polya, std::mutex& iomutex,
              std::shared_ptr<spdlog::logger> log, std::string outFile,
              std::vector<uint32_t>& refIdExtensions,
              std::vector<std::pair<std::string, uint16_t>>& shortRefs) {
  (void)iomutex;

  ghc::filesystem::path outFilePath{outFile};
  ghc::filesystem::path outDir = outFilePath.parent_path();

  // std::shared_ptr<spdlog::logger> log) {
  // Create a random uniform distribution
  std::default_random_engine eng(271828);
  std::uniform_int_distribution<> dis(0, 3);

  // Hashers for getting txome signature
  digestpp::sha256 seqHasher256;
  digestpp::sha256 nameHasher256;
  digestpp::sha512 seqHasher512;
  digestpp::sha512 nameHasher512;
  digestpp::sha256 decoySeqHasher256;
  digestpp::sha256 decoyNameHasher256;

  // Keep track of if we've seen a decoy sequence yet.
  // The index enforces that all decoy sequences are consecutive, and that
  // they come after all valid (non-decoy) sequences.  If we see a non-decoy
  // sequence after having observed a decoy, then we complain and exit.
  bool sawDecoy{false};
  uint64_t numberOfDecoys{0};
  uint64_t numberOfDuplicateDecoys{0};
  uint64_t firstDecoyIndex{std::numeric_limits<uint64_t>::max()};

  bool firstRecord{true};
  bool hasGencodeSep = (sepStr.find('|') != std::string::npos);
  uint32_t n{0};
  std::vector<std::string> transcriptNames;
  std::unordered_set<std::string> transcriptNameSet;
  std::map<std::string, bool> shortFlag ;

  std::vector<int64_t> transcriptStarts;
  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  uint32_t polyAClipLength{10};
  uint32_t numPolyAsClipped{0};
  uint32_t numNucleotidesReplaced{0};
  std::string polyA(polyAClipLength, 'A');

  //using TranscriptList = std::vector<uint32_t>;
  //using KmerBinT = uint64_t;

  bool haveDecoys = !decoyNames.empty();
  bool clipPolyA = !(noclip_polya);

  struct DupInfo {
    uint64_t txId;
    uint64_t txOffset;
    uint32_t txLen;
  };

  auto update_name_hash = [&nameHasher256, &nameHasher512, &decoyNameHasher256](bool is_decoy, const std::string& processed_name) -> void {
    if (is_decoy) {
      decoyNameHasher256.absorb(processed_name);
    } else {
      nameHasher256.absorb(processed_name);
      nameHasher512.absorb(processed_name);
    }
  };

  auto update_seq_hash = [&seqHasher256, &seqHasher512, &decoySeqHasher256](bool is_decoy, const std::string& seq) -> void {
    if (is_decoy) {
      decoySeqHasher256.absorb(seq.begin(), seq.end());
    } else {
      seqHasher256.absorb(seq.begin(), seq.end());
      seqHasher512.absorb(seq.begin(), seq.end());
    }
  };

  // This used to be set based on the length of Titin (108861).  But, as is almost guaranteed 
  // with any arbitrary cutoff in code, the data will eventually prove you chose an inadequate
  // value for some case.  The new winner for longest transcript is ENST00000674361.1
  // weighing in at 347,561 nucleotides.  It appears in human Gencode 35 (and probably will 
  // persist).  So, we set this value to 400000 to give some headway and avoid producing 
  // warnings for standard human transcriptomes.  This addresses 
  // https://github.com/COMBINE-lab/salmon/issues/591; thanks to @guidohooiveld for reporting 
  // it.
  constexpr size_t tooLong = 400000;

  size_t currIndex{0};
  size_t numDups{0};
  int64_t numShortBeforeFirstDecoy{0};
  std::map<XXH64_hash_t, std::vector<DupInfo>> potentialDuplicates;
  spp::sparse_hash_map<uint64_t, std::vector<std::string>> duplicateNames;
  std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

  // rsdic::RSDicBuilder rsdb;
  std::vector<uint64_t>
      onePos; // Positions in the bit array where we should write a '1'
  // remember the initial lengths (e.g., before clipping etc., of all
  // transcripts)
  std::vector<uint32_t> completeLengths;
  // the stream of transcript sequence
  fmt::MemoryWriter txpSeqStream;
  {
    // ScopedTimer timer;
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();
    bool tooShort{false};
    while (parser->refill(rg)) {
      for (auto& read : rg) { // for each sequence
        tooShort = false;
        std::string& readStr = read.seq;
        readStr.erase(
            std::remove_if(readStr.begin(), readStr.end(),
                           [](const char a) -> bool { return !(isprint(a)); }),
            readStr.end());

        uint32_t readLen = readStr.size();
        uint32_t completeLen = readLen;

        // get the hash to check for collisions before we change anything.
        auto txStringHash =
            XXH64(reinterpret_cast<void*>(const_cast<char*>(readStr.data())),
                  readLen, 0);
        auto& readName = read.name;

        // check if we think this is a gencode transcriptome, and the user has not passed the gencode flag
        if (firstRecord and !hasGencodeSep) {
          constexpr const size_t numGencodeSep{8};
          if ( std::count(readName.begin(), readName.end(), '|') == numGencodeSep ) {
            log->warn("It appears that this may be a GENCODE transcriptome (from analyzing the separators in the FASTA header).  However, "
                      "you have not set \'|\' as a header separator.  If this is a GENCODE transcriptome, consider passing --gencode to the "
                      "pufferfish index command.\n\n");
          }
          firstRecord = false;
        }

        bool isDecoy = (haveDecoys) ? decoyNames.contains(readName) : false;
        // If this is *not* a decoy sequence, make sure that
        // we haven't seen any decoys yet.  Otherwise we are violating
        // the condition that decoys must come last.
        if (!isDecoy and sawDecoy) {
          log->critical("Observed a non-decoy sequence [{}] after having already observed a decoy. "
                        "However, it is required that any decoy target records appear, consecutively, "
                        "at the end of the input fasta file.  Please re-format your input file so that "
                        "all decoy records appear contiguously at the end of the file, after all valid "
                        "(non-decoy) records", readName);
          log->flush();
          spdlog::drop_all();
          std::exit(1);
        }

        // If this was a decoy, add it to the decoy hash
        update_seq_hash(isDecoy, readStr);

        // First, replace non ATCG nucleotides
        for (size_t b = 0; b < readLen; ++b) {
          readStr[b] = ::toupper(readStr[b]);
          int c = jellyfish::mer_dna::code(readStr[b]);
          // Replace non-ACGT bases with pseudo-random bases
          if (jellyfish::mer_dna::not_dna(c)) {
            char rbase = bases[dis(eng)];
            c = jellyfish::mer_dna::code(rbase);
            readStr[b] = rbase;
            ++numNucleotidesReplaced;
          }
        }

        // Now, do Kallisto-esque clipping of polyA tails
        if (clipPolyA) {
          if (readStr.size() > polyAClipLength and
              readStr.substr(readStr.length() - polyAClipLength) == polyA) {
            auto newEndPos = readStr.find_last_not_of("Aa");
            // If it was all As
            if (newEndPos == std::string::npos) {
              log->warn("Entry with header [{}] appeared to be all A's; it "
                        "will be removed from the index!",
                        read.name);
              readStr.resize(0);
            } else {
              readStr.resize(newEndPos + 1);
            }
            ++numPolyAsClipped;
          }
        }

        readLen = readStr.size();
        // If the transcript was completely removed during clipping, don't
        // include it in the index.
        if (readStr.size() > 0) {

          // If we're suspicious the user has fed in a *genome* rather
          // than a transcriptome, say so here.
          if (readStr.size() >= tooLong and !isDecoy and expect_transcriptome) {
            log->warn("Entry with header [{}] was longer than {} nucleotides.  "
                      "This is probably a chromosome instead of a transcript.",
                      read.name, tooLong);
          } else if (readStr.size() <= k) { // <= instead of < because of twopaco!
            log->warn("Entry with header [{}], had length less than equal to "
                      "the k-mer length of {} (perhaps after poly-A clipping)",
                      read.name, k);
            tooShort = true;
          }

          uint32_t txpIndex = n++;

          // The name of the current transcript
          auto& recHeader = read.name;
          auto processedName =
              recHeader.substr(0, recHeader.find_first_of(sepStr));

          update_name_hash(isDecoy, processedName);

          // Add this transcript, indexed by it's sequence's hash value
          // to the potential duplicate list.
          bool didCollide{false};
          auto dupIt = potentialDuplicates.find(txStringHash);
          if (dupIt != potentialDuplicates.end()) {
            auto& dupList = dupIt->second;
            for (auto& dupInfo : dupList) {
              // they must be of the same length
              if (readLen == dupInfo.txLen) {
                bool collision =
                    (readStr.compare(0, readLen,
                                     txpSeqStream.data() + dupInfo.txOffset,
                                     readLen) == 0);
                if (collision) {
                  ++numDups;
                  didCollide = true;
                  duplicateNames[dupInfo.txId].push_back(processedName);
                  continue;
                } // if collision
              }   // if readLen == dupInfo.txLen
            }     // for dupInfo : dupList
          }       // if we had a potential duplicate

          // if this was a duplicate and a decoy sequence
          // then we don't care about the status of the `--keepDuplicates` 
          // flag.  It never really makes sense to keep a duplicate 
          // decoy.
          if (didCollide and isDecoy) {
            // roll back the txp index & skip the rest of this loop
            n--;
            ++numberOfDuplicateDecoys;
            continue;
          }

          // if this was a duplicate and not a decoy, then take proper
          // action based on the status of `--keepDuplicates`
          if (!keepDuplicates and didCollide) {
            // roll back the txp index & skip the rest of this loop
            n--;
            continue;
          }

          // Check for duplicate name
          if (transcriptNameSet.find(processedName) != transcriptNameSet.end()) {
            log->error("In FixFasta, two references with the same name but different sequences: {}. "
                       "We require that all input records have a unique name "
                       "up to the first whitespace (or user-provided separator) character.", processedName);
            std::exit(1);
          }
          // If there was no collision, then add the transcript
          transcriptNameSet.insert(processedName);
          transcriptNames.emplace_back(processedName);

          if(!tooShort) {
              shortFlag[processedName] = false;
          } else {
              numShortBeforeFirstDecoy += sawDecoy ? 0 : 1;
              shortFlag[processedName] = true;
          }
          // nameHasher.process(processedName.begin(), processedName.end());

          // The position at which this transcript starts
          transcriptStarts.push_back(currIndex);
          // The un-molested length of this transcript
          completeLengths.push_back(completeLen);

          if (isDecoy) {
            // if we haven't seen another decoy yet, this is the first decoy
            // index
            if (!sawDecoy) {
              firstDecoyIndex = txpIndex;
            }
            // once we see the first decoy, saw decoy is set to true
            // for the rest of the processing.
            sawDecoy = true;
            ++numberOfDecoys;
            //decoyIndices.push_back(txpIndex);
          }

          // If we made it here, we were not an actual duplicate, so add this
          // transcript
          // for future duplicate checking.
          if (!keepDuplicates or (keepDuplicates and !didCollide)) {
            potentialDuplicates[txStringHash].push_back(
                {txpIndex, currIndex, readLen});
          }

          txpSeqStream << readStr;
          currIndex += readLen;
          onePos.push_back(currIndex);
        } else {
          log->warn("Discarding entry with header [{}], since it had length 0 "
                    "(perhaps after poly-A clipping)",
                    read.name);
        }
      }
      if (n % 10000 == 0) {
        std::cerr << "\r\rcounted k-mers for " << n << " transcripts";
      }
    }
  }
  std::cerr << "\n";
  if (numDups > 0) {
    if (!keepDuplicates) {
      log->warn("Removed {} transcripts that were sequence duplicates of "
                "indexed transcripts.",
                numDups);
      log->warn("If you wish to retain duplicate transcripts, please use the "
                "`--keepDuplicates` flag");
    } else {
      log->warn("There were {} transcripts that would need to be removed to "
                "avoid duplicates.",
                numDups);
    }
  }

  if ((numberOfDecoys + numberOfDuplicateDecoys) != decoyNames.size()) {
    log->critical("The decoy file contained the names of {} decoy sequences, but "
    "{} were matched by sequences in the reference file provided. To prevent unintentional "
    "errors downstream, please ensure that the decoy file exactly matches with the "
    "fasta file that is being indexed.", decoyNames.size(), numberOfDecoys);
    return false;
  }

  if (numberOfDuplicateDecoys > 0) {
    log->warn("There were {} duplicate decoy sequences.", numberOfDuplicateDecoys);
  }

  {
    ghc::filesystem::path dcPath = outDir / ghc::filesystem::path{"duplicate_clusters.tsv"};
    std::ofstream dupClusterStream(dcPath.string());
    dupClusterStream << "RetainedRef" << '\t' << "DuplicateRef" << '\n';
    for (auto kvIt = duplicateNames.begin(); kvIt != duplicateNames.end(); ++kvIt) {
      auto& retainedName = transcriptNames[kvIt->first];
      for (auto& droppedName : kvIt->second) {
        dupClusterStream << retainedName << '\t' << droppedName << '\n';
      }
    }
    dupClusterStream.close();
  }

  /*
  std::ofstream dupClusterStream(outputDir + "duplicate_clusters.tsv");
  {
    dupClusterStream << "RetainedTxp" << '\t' << "DuplicateTxp" << '\n';
    for (auto kvIt = duplicateNames.begin(); kvIt != duplicateNames.end();
  ++kvIt) {
      auto& retainedName = transcriptNames[kvIt->first];
      for (auto& droppedName : kvIt->second) {
        dupClusterStream << retainedName << '\t' << droppedName << '\n';
      }
    }
  }
  dupClusterStream.close();
  */

  log->info("Replaced {:n} non-ATCG nucleotides", numNucleotidesReplaced);
  log->info("Clipped poly-A tails from {:n} transcripts", numPolyAsClipped);

  // Put the concatenated text in a string
  std::string concatText = txpSeqStream.str();
  stx::string_view concatTextView(concatText);
  // And clear the stream
  txpSeqStream.clear();

  std::ofstream ffa(outFile);
  size_t prev1{0};
  size_t numWritten{0};
  uint32_t prevExt{0};
  refIdExtensions.reserve(transcriptNames.size());
  for (size_t i = 0; i < transcriptNames.size(); ++i) {
    size_t next1 = onePos[i];
    size_t len = next1 - prev1;
    if(!shortFlag[transcriptNames[i]]){
        ffa << ">" << transcriptNames[i] << "\n";
        ffa << concatTextView.substr(prev1, len) << "\n";
        refIdExtensions.push_back(prevExt);
        ++numWritten;
    } else {
      shortRefs.emplace_back(transcriptNames[i], len);
      prevExt++;
    }
    prev1 = next1;
  }
  ffa.close();
  std::cerr << "wrote " << numWritten << " cleaned references\n";


  {
    ghc::filesystem::path sigPath = outDir / ghc::filesystem::path{"complete_ref_lens.bin"};
    std::ofstream os(sigPath.string(), std::ios::binary);
    cereal::BinaryOutputArchive lenArchive(os);
    lenArchive(completeLengths);
  }

  // Set the hash info
  std::string seqHash256 = seqHasher256.hexdigest();
  std::string nameHash256 = nameHasher256.hexdigest();
  std::string seqHash512 = seqHasher512.hexdigest();
  std::string nameHash512 = nameHasher512.hexdigest();
  std::string decoySeqHash256 = decoySeqHasher256.hexdigest();
  std::string decoyNameHash256 = decoyNameHasher256.hexdigest();

  {
    ghc::filesystem::path sigPath = outDir / ghc::filesystem::path{"ref_sigs.json"};
    std::ofstream os(sigPath.string());
    cereal::JSONOutputArchive ar(os);
    auto adjustedFirstDecoyIndex = firstDecoyIndex - numShortBeforeFirstDecoy;
    ar( cereal::make_nvp("keep_duplicates", keepDuplicates));
    ar( cereal::make_nvp("num_decoys", numberOfDecoys));
    ar( cereal::make_nvp("first_decoy_index", adjustedFirstDecoyIndex));
    ar( cereal::make_nvp("SeqHash", seqHash256) );
    ar( cereal::make_nvp("NameHash", nameHash256) );
    ar( cereal::make_nvp("SeqHash512", seqHash512) );
    ar( cereal::make_nvp("NameHash512", nameHash512) );
    ar( cereal::make_nvp("DecoySeqHash", decoySeqHash256) );
    ar( cereal::make_nvp("DecoyNameHash", decoyNameHash256) );
  }

  /*  header.setSeqHash256(seqHash256);
  header.setNameHash256(nameHash256);
  header.setSeqHash512(seqHash512);
  header.setNameHash512(nameHash512);
  header.setDecoySeqHash256(decoySeqHash256);
  header.setDecoyNameHash256(decoyNameHash256);
  */
 return true;
}

/**
 * parses the decoy file and constructs a sparse_hash_set that contains the names of the decoys.
 *
 * Note: If there are any errors, this will return an empty hash set.  In the future, using 
 * std::optional or some other mechanism to indicate failure would be better.
 **/
spp::sparse_hash_set<std::string> populateDecoyHashPuff(const std::string& fname, std::shared_ptr<spdlog::logger> log) {
  spp::sparse_hash_set<std::string> dset;
  std::ifstream dfile(fname);

  bool had_duplicate = false;
  std::string dname;
  while (dfile >> dname) {
    auto it = dset.insert(dname);
    if (!it.second) {
      log->critical("The decoy name {} was encountered more than once --- please ensure all decoy names and sequences are unique.", dname);
      had_duplicate = true;
      break;
    }
  }

  if (dset.empty()) {
    log->critical("The decoy file was empty.  If you have no decoys, then you should not pass the `--decoys` option while indexing.");
  }
  if (had_duplicate) {
    dset.clear();
  }
  
  dfile.close();
  return dset;
}

bool extractFasta(std::string& inputTsv, std::string& outFile, uint32_t& numFeats) {
  std::ifstream ifile(inputTsv);
  std::ofstream ofile(outFile);

  digestpp::sha256 seqHasher256;
  digestpp::sha256 nameHasher256;
  digestpp::sha512 seqHasher512;
  digestpp::sha512 nameHasher512;
  digestpp::sha256 decoySeqHasher256;
  digestpp::sha256 decoyNameHasher256;

  size_t seqLen{0};
  std::unordered_set<std::string> names;
  std::unordered_set<std::string> seqs;

  if(ifile.is_open()) {
    while( true ) {
      std::string featStr, seqStr;
      ifile >> featStr >> seqStr;
      if( ifile.eof() ) { break; }

      ofile << ">" << featStr << "\n" << seqStr << "\n";
      ++numFeats;

      nameHasher256.absorb(featStr);
      nameHasher512.absorb(featStr);

      seqHasher256.absorb(seqStr.begin(), seqStr.end());
      seqHasher512.absorb(seqStr.begin(), seqStr.end());

      if (seqLen == 0) { seqLen = seqStr.size(); }
      if (seqLen > 0 && seqLen != seqStr.size()) {
        std::cerr << "CRITICAL ERROR: The length of all the feature sequences should be the same."
                  << std::endl;
        std::exit(1);
      }

      if (seqs.find(seqStr) == seqs.end()) {
        seqs.insert(seqStr);
      } else {
        std::cerr << "CRITICAL ERROR: The feature sequences should be unique"
                  << std::endl;
        std::exit(1);
      }

      if (names.find(featStr) == names.end()) {
        names.insert(featStr);
      } else {
        std::cerr << "CRITICAL ERROR: The feature ids should be unique"
                  << std::endl;
        std::exit(1);
      }
    }
    ifile.close();
    ofile.close();
  }

  { // this block mostly copy paste
    // Set the hash info
    std::string seqHash256 = seqHasher256.hexdigest();
    std::string nameHash256 = nameHasher256.hexdigest();
    std::string seqHash512 = seqHasher512.hexdigest();
    std::string nameHash512 = nameHasher512.hexdigest();
    std::string decoySeqHash256 = decoySeqHasher256.hexdigest();
    std::string decoyNameHash256 = decoyNameHasher256.hexdigest();

    ghc::filesystem::path outFilePath{outFile};
    ghc::filesystem::path outDir = outFilePath.parent_path();
    { // writing the json
      ghc::filesystem::path sigPath = outDir / ghc::filesystem::path{"ref_sigs.json"};
      std::ofstream os(sigPath.string());
      cereal::JSONOutputArchive ar(os);

      uint64_t numShortBeforeFirstDecoy {0}; // forcing to all the features to be of equal length
      uint64_t firstDecoyIndex{std::numeric_limits<uint64_t>::max()};
      auto adjustedFirstDecoyIndex = firstDecoyIndex - numShortBeforeFirstDecoy;

      ar( cereal::make_nvp("keep_duplicates", false));
      ar( cereal::make_nvp("num_decoys", 0));
      ar( cereal::make_nvp("first_decoy_index", adjustedFirstDecoyIndex));
      ar( cereal::make_nvp("SeqHash", seqHash256) );
      ar( cereal::make_nvp("NameHash", nameHash256) );
      ar( cereal::make_nvp("SeqHash512", seqHash512) );
      ar( cereal::make_nvp("NameHash512", nameHash512) );
      ar( cereal::make_nvp("DecoySeqHash", decoySeqHash256) );
      ar( cereal::make_nvp("DecoyNameHash", decoyNameHash256) );
    }

    {
      std::vector<uint32_t> completeLengths(numFeats, seqLen);
      ghc::filesystem::path sigPath = outDir / ghc::filesystem::path{"complete_ref_lens.bin"};
      std::ofstream os(sigPath.string(), std::ios::binary);
      cereal::BinaryOutputArchive lenArchive(os);
      lenArchive(completeLengths);
    }

    {
      ghc::filesystem::path dcPath = outDir / ghc::filesystem::path{"duplicate_clusters.tsv"};
      std::ofstream dupClusterStream(dcPath.string());
      dupClusterStream << "RetainedRef" << '\t' << "DuplicateRef" << '\n';
      dupClusterStream.close();
    }
  }

  return true;
}


int fixFastaMain(std::vector<std::string>& args,
        std::vector<uint32_t>& refIdExtension,
        std::vector<std::pair<std::string, uint16_t>>& shortRefs,
        std::shared_ptr<spdlog::logger> log,
        bool hasFeatures) {
  using namespace clipp;

  uint32_t k{31};
  std::vector<std::string> refFiles;
  std::string outFile;
  std::string decoyFile;
  bool keepDuplicates{false};
  bool printHelp{false};
  bool expect_transcriptome{false};
  bool noclip_polya{false};
  std::string sepStr{" \t"};

  auto cli = (
              option("--help", "-h").set(printHelp, true) % "show usage",
              required("--input", "-i") & values("input", refFiles) % "input FASTA file",
              required("--output", "-o") & value("output", outFile) % "output FASTA file",
              option("--headerSep", "-s") & value("sep_strs", sepStr) %
              "Instead of a space or tab, break the header at the first "
              "occurrence of this string, and name the transcript as the token before "
              "the first separator (default = space & tab)",
              option("--expectTranscriptome").set(expect_transcriptome) % 
              "expect (non-decoy) sequences to be transcripts rather than genomic contigs",
              option("--noClip", "-n").set(noclip_polya) % "Don't clip poly-A tails from the ends of target sequences",
              option("--decoys", "-d") & value("decoys", decoyFile) %
              "Treat these sequences as decoys that may be sequence-similar to some known indexed reference",
              option("--keepDuplicates").set(keepDuplicates) % "Retain duplicate references in the input",
              option("--klen", "-k") & value("k-mer length", k) % "length of the k-mer used to build the cDBG (default = 31)"
              );

  //  if (parse(argc, argv, cli)) {
  if (parse(args, cli)) {
    if (printHelp) {
      std::cout << make_man_page(cli, "fixFasta");
      return 0;
    }

    //auto console = spdlog::stderr_color_mt("ff::console");

    spp::sparse_hash_set<std::string> decoyNames;
    if (!decoyFile.empty()) {
      bool decoyFileExists = ghc::filesystem::exists(decoyFile);
      if (!decoyFileExists) {
        log->error("The decoy file {} does not exist.", decoyFile);
        std::exit(1);
      }
      decoyNames = populateDecoyHashPuff(decoyFile, log);
      if (decoyNames.empty()) {
        return 1;
      }
    }

    bool fix_ok {false};
    if (hasFeatures) {
      uint32_t numFeats{0};
      fix_ok = extractFasta(refFiles[0], outFile, numFeats);
      refIdExtension.resize(numFeats, 0);
    } else {
      size_t numThreads{1};
      std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
      size_t numProd = 1;

      transcriptParserPtr.reset(new single_parser(refFiles, numThreads, numProd));
      transcriptParserPtr->start();
      std::mutex iomutex;
      fix_ok = fixFasta(transcriptParserPtr.get(), decoyNames, keepDuplicates, k, sepStr, expect_transcriptome, 
                        noclip_polya, iomutex, log, outFile, refIdExtension, shortRefs);
      transcriptParserPtr->stop();
    }

    return fix_ok ? 0 : 1;
  } else {
    std::cout << usage_lines(cli, "fixFasta") << '\n';
    return 1;
  }
}
