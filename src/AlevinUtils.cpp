#include "AlevinUtils.hpp"

namespace alevin {
  namespace utils {

    std::unordered_map<uint32_t,uint32_t> lenOffset({
        {8,  0},
          {9,  pow(4,8)},
            {10, 5*pow(4,8)},
              {11, 21*pow(4,8)},
                {12, 85*pow(4,8)}});


    std::unordered_map<std::string, std::string> iupacMap ({
        //https://github.com/brwnj/umitools/blob/master/umitools/umitools.py#L26
        //cell-barcodes should follow IUPAC mapping as following:
        {"A","A"},
          {"T","T"},
            {"C","C"},
              {"G","G"},
                {"R","GA"},
                  {"Y","TC"},
                    {"M","AC"},
                      {"K","GT"},
                        {"S","GC"},
                          {"W","AT"},
                            {"H","ACT"},
                              { "B","GTC"},
                                {"V","GCA"},
                                  {"D","GAT"},
                                    {"N","GATC"}
      });

    std::unordered_map<uint32_t, std::string> ntMap({ {0,"A"}, {1,"C"},
                                                      {2,"T"}, {3,"G"}});

    std::unordered_map<std::string, std::string> revNtMap({ {"C","G"}, {"A","T"},
                                                            {"G","C"}, {"T","A"} });

    std::vector<std::string> nts{"A", "T", "C", "G"};


    template <>
    bool extractUMI<apt::DropSeq>(std::string& read,
                                  apt::DropSeq& pt,
                                  std::string& umi){
      std::cout<<"Incorrect call for umi extract";
      exit(0);
    }
    template <>
    bool extractUMI<apt::Chromium>(std::string& read,
                                   apt::Chromium& pt,
                                   std::string& umi){
      umi = read.substr(pt.barcodeLength, pt.umiLength);
      return true;
    }
    template <>
    bool extractUMI<apt::Gemcode>(std::string& read,
                                   apt::Gemcode& pt,
                                   std::string& umi){
      umi = read.substr(pt.barcodeLength, pt.umiLength);
      return true;
    }
    template <>
    bool extractUMI<apt::Custom>(std::string& read,
                                 apt::Custom& pt,
                                 std::string& umi){
      umi = read.substr(pt.barcodeLength, pt.umiLength);
      return true;
    }
    template <>
    bool extractUMI<apt::InDrop>(std::string& read,
                                 apt::InDrop& pt,
                                 std::string& umi){
      std::cout<<"Incorrect call for umi extract";
      exit(0);
    }


    template <>
    bool extractBarcode<apt::DropSeq>(std::string& read,
                                      apt::DropSeq& pt,
                                      std::string& bc){
      bc = read.substr(0, pt.barcodeLength);
      return true;
    }
    template <>
    bool extractBarcode<apt::Chromium>(std::string& read,
                                       apt::Chromium& pt,
                                       std::string& bc){
      bc = read.substr(0, pt.barcodeLength);
      return true;
    }
    template <>
    bool extractBarcode<apt::Gemcode>(std::string& read,
                                       apt::Gemcode& pt,
                                       std::string& bc){
      bc = read.substr(0, pt.barcodeLength);
      return true;
    }
    template <>
    bool extractBarcode<apt::Custom>(std::string& read,
                                     apt::Custom& pt,
                                     std::string& bc){
      bc = read.substr(0, pt.barcodeLength);
      return true;
    }
    template <>
    bool extractBarcode<apt::InDrop>(std::string& read,
                                     apt::InDrop& pt,
                                     std::string& bc){
      std::string::size_type index = read.find(pt.w1);
      if (index == std::string::npos){
        return false;
      }
      bc = read.substr(0, index);
      if(bc.size()<8 or bc.size()>12){
        return false;
      }
      uint32_t offset = bc.size()+pt.w1.size();
      bc += read.substr(offset, offset+8);
      return true;
    }

    void getIndelNeighbors(
                           const std::string& barcodeSeq,
                           std::unordered_set<std::string>& neighbors){
      size_t barcodeLength = barcodeSeq.size();
      std::string newBarcode;

      for (size_t i=0; i<barcodeLength; i++){
        for(auto j: nts){
          if (i != barcodeLength-1){
            //insert
            newBarcode = barcodeSeq;
            newBarcode.insert(i, j);
            newBarcode.erase(barcodeLength, 1);
            neighbors.insert(newBarcode);

            //deletion
            newBarcode = barcodeSeq;
            newBarcode.erase(i, 1);
            newBarcode.insert(barcodeLength-1, j);
            neighbors.insert(newBarcode);
          }//endif
        }//end-j-for
      }//end-i-for
    }
    //void getIndelNeighbors(
    //                       std::string& barcodeSeq,
    //                       std::unordered_set<std::string>& neighbors){
    //  size_t barcodeLength = barcodeSeq.size();
    //  std::string newBarcode;

    //  for (size_t i=0; i<=barcodeLength; i++){
    //    for(auto j: nts){
    //      if(barcodeLength<12){
    //        //insert
    //        newBarcode = barcodeSeq;
    //        newBarcode.insert(i, j);
    //        neighbors.insert(newBarcode);
    //      }
    //    }//end-j-for

    //    if(barcodeLength>8 and i<barcodeLength){
    //      //deletion
    //      newBarcode = barcodeSeq;
    //      newBarcode.erase(i, 1);
    //      neighbors.insert(newBarcode);
    //    }
    //  }//end-i-for
    //}


    /*
      Finds all 1 edit-distance neighbor of a given barcode.
     */
    void findNeighbors(size_t seqSize,
                       const std::string& barcodeSeq,
                       std::unordered_set<std::string>& neighbors){
      size_t barcodeLength { barcodeSeq.size() };
      std::string newBarcode, nt;

      if(barcodeLength > seqSize){
        std::cout<<"Sequence-Size " << barcodeLength << "greater than specified "
                 << seqSize <<".\nPlease report the issue on Github.\n" ;
        exit(1);
      }

      for (size_t i=0; i<barcodeLength; i++){
        nt = barcodeSeq[i];
        for(auto j: nts){
          if (nt == j){
            continue;
          }
          //snp
          newBarcode = barcodeSeq;
          newBarcode.replace(i, 1, j);
          neighbors.insert(newBarcode);
        }
      }//end-i-for
      getIndelNeighbors(barcodeSeq,
                        neighbors);
    }


    template <typename ProtocolT>
    bool processAlevinOpts(AlevinOpts<ProtocolT>& aopt,
                           SalmonOpts& sopt,
                           boost::program_options::variables_map& vm){
      namespace bfs = boost::filesystem;
      namespace po = boost::program_options;

      //Create outputDirectory
      aopt.outputDirectory = vm["output"].as<std::string>() + "/alevin";
      if (!bfs::exists(aopt.outputDirectory)) {
        bool dirSuccess = boost::filesystem::create_directories(aopt.outputDirectory);
        if (!dirSuccess) {
          fmt::print(stderr,"\nCould not create output directory {}\nExiting Now.",
                     aopt.outputDirectory.string());
          return false;
        }
      }

      if (vm.count("mrna")){
        aopt.mRnaFile = vm["mrna"].as<std::string>();
        if (!bfs::exists(aopt.mRnaFile)) {
          fmt::print(stderr,"\n mRna File {} does not exists\n Exiting Now",
                     aopt.mRnaFile.string());
          return false;
        }
      }

      if (vm.count("rrna")){
        aopt.rRnaFile = vm["rrna"].as<std::string>();
        if (!bfs::exists(aopt.rRnaFile)) {
          fmt::print(stderr,"\nrRna File {} does not exists\n Exiting Now",
                     aopt.rRnaFile.string());
          return false;
        }
      }

      if (vm.count("whitelist")){
        aopt.whitelistFile = vm["whitelist"].as<std::string>();
        if (!bfs::exists(aopt.whitelistFile)) {
          fmt::print(stderr,"\nWhitelist File {} does not exists\n Exiting Now",
                     aopt.whitelistFile.string());
          return false;
        }
      }

      if (vm.count("tgMap")){
        aopt.geneMapFile = vm["tgMap"].as<std::string>();
        if (!bfs::exists(aopt.geneMapFile)) {
          fmt::print(stderr,"\nTranscript to Gene Map File {} does not exists\n Exiting Now",
                     aopt.geneMapFile.string());
          return false;
        }
      }
      else{
        fmt::print(stderr,"\nTranscript to Gene Map File not provided\n Exiting Now");
        return false;
      }

      //create logger
      spdlog::set_async_mode(131072);
      auto logPath = aopt.outputDirectory / "alevin.log";
      auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
      auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>();
      std::vector<spdlog::sink_ptr> sinks{consoleSink, fileSink};
      aopt.jointLog = spdlog::create("alevinLog", std::begin(sinks), std::end(sinks));

      uint32_t barEnd{0}, barcodeLength{0}, umiLength{0};

      if(vm.count("end") and vm.count("barcodeLength") and vm.count("umiLength")){
        barEnd = vm["end"].as<uint32_t>();
        barcodeLength = vm["barcodeLength"].as<uint32_t>();
        umiLength = vm["umiLength"].as<uint32_t>();
      }

      aopt.quiet = vm["quiet"].as<bool>();
      aopt.doEM = vm["em"].as<bool>();
      aopt.useCorrelation = vm["useCorrelation"].as<bool>();
      aopt.noDedup = vm["noDedup"].as<bool>();
      aopt.naive = vm["naive"].as<bool>();
      aopt.noQuant = vm["noQuant"].as<bool>();
      aopt.noSoftMap = vm["noSoftMap"].as<bool>();
      aopt.dumpfq = vm["dumpfq"].as<bool>();
      aopt.dumpBFH = vm["dumpBfh"].as<bool>();
      aopt.txpLevel = vm["txpLevel"].as<bool>();
      aopt.eqClassLevel = vm["eqClassLevel"].as<bool>();
      aopt.nobarcode = vm["noBarcode"].as<bool>();
      aopt.dumpfeatures = vm["dumpFeatures"].as<bool>();
      aopt.dumpCsvCounts = vm["dumpCsvCounts"].as<bool>();
      aopt.dumpBarcodeEq = vm["dumpBarcodeEq"].as<bool>();
      aopt.dumpBarcodeMap = vm["dumpBarcodeMap"].as<bool>();
      aopt.dumpUmiToolsMap = vm["dumpUmitoolsMap"].as<bool>();
      aopt.lowRegionMinNumBarcodes = vm["lowRegionMinNumBarcodes"].as<uint32_t>();
      aopt.maxNumBarcodes = vm["maxNumBarcodes"].as<uint32_t>();
      if(vm.count("iupac")){
        aopt.iupac = vm["iupac"].as<std::string>();
      }

      if (not vm.count("threads")) {
        auto tot_cores = std::thread::hardware_concurrency();
        aopt.numThreads = std::max(1, static_cast<int>(tot_cores/4.0));
        aopt.jointLog->warn("threads flag not specified, Using {}(25%) of total cores",
                       aopt.numThreads);
      }
      else{
        aopt.numThreads = vm["threads"].as<uint32_t>();
      }  // things which needs to be updated for salmonOpts

      //validate customized options for custom protocol
      if (barEnd or barcodeLength or umiLength
          or !aopt.protocol.barcodeLength or !aopt.protocol.umiLength){
        if(!barEnd or !barcodeLength or !umiLength){
          aopt.jointLog->error("\nERROR: if you are using any one"
                               " of (end, umilength, barcodelength) flag\n"
                               "you have to provide all of them explicitly.\n"
                               "You can also use pre-defined single-cell protocols."
                               "Exiting Now.");
          return false;
        }

        aopt.protocol.barcodeLength = barcodeLength;
        aopt.protocol.umiLength = umiLength;
        if (barEnd == 3) {
          aopt.protocol.end = BarcodeEnd::THREE;
        }
        else if (barEnd == 5) {
          aopt.protocol.end = BarcodeEnd::FIVE;
        }
        else{
          aopt.jointLog->error("\nERROR: Wrong value for Barcode-end of read -> {}"
                               "\nExiting now: please provide `5` for barcodes "
                               "starting at 5' end or `3`` for barcode starting "
                               "at 3' end.\n", barEnd);
          return false;
        }
        if (barcodeLength > 60 or umiLength > 12){
          aopt.jointLog->error("\nERROR: umiLength/barcodeLength too long\n"
                               "Exiting Now;");
          return false;
        }
      }

      //validate specified iupac
      if (aopt.iupac.size()>0){
        std::string correctIUPACodes = "ACTGRYMKSWHBVDN";

        if(aopt.iupac.size() != aopt.protocol.barcodeLength){
          aopt.jointLog->error("\nERROR: iupac length: {} and barcode"
                               " length: {} do not match. Please check command "
                               "line flags.\n Exiting Now.", aopt.iupac.size(),
                               aopt.protocol.barcodeLength);
          return false;
        }

        bool allNflag{true};
        for (auto x : aopt.iupac){
          size_t found = correctIUPACodes.find(x);
          if (x != 'N'){
            allNflag = false;
          }
          if (found==std::string::npos){
            aopt.jointLog->error("\nERROR: Wrong IUPAC charachter {} in {}\n"
                                 "\nExiting now: Please check "
                                 "https://www.bioinformatics.org/sms/iupac.html"
                                 "for more details about iupac.",
                                 x, aopt.iupac);
            return false;
          }
        }
        if (allNflag){
          aopt.iupac.clear();
        }
      }

      // code from SalmonAlevin
      sopt.numThreads = aopt.numThreads;
      sopt.quiet = aopt.quiet;
      sopt.quantMode = SalmonQuantMode::MAP;
      bool optionsOK =
        salmon::utils::processQuantOptions(sopt, vm, vm["numBiasSamples"].as<int32_t>());
      if (!optionsOK) {
        if (aopt.jointLog) {
          aopt.jointLog->error("Could not properly process salmon-level options!");
          aopt.jointLog->flush();
          spdlog::drop_all();
        }
        return false;
      }
      return true;
    }

    //template <typename ProtocolT>
    bool sequenceCheck(const std::string& sequence,
                       //AlevinOpts<ProtocolT>& aopt,
                       //std::mutex& iomutex,
                       Sequence seqType){
      //size_t lenSequnce;
      //auto log = aopt.jointLog;
      //lenSequnce = sequence.length();
      if (sequence.length() == 0){
        return false;
      }

      for(auto nt: sequence) {
        if('N' == nt){
          return false;
        }
        switch(nt){
        case 'A':
        case 'C':
        case 'G':
        case 'T':break;
        default:
          {
            //log->error("Unidentified charachter {} --- non(A,C,G,T,N)"
            //           " in barcode sequence", nt);
            return false;
          }
        }//end-switch
      }//end-for
      return true;
    }

    bool checkSetCoverage(std::vector<std::vector<uint32_t>>& tgroup,
                          std::vector<uint32_t> txps){
      // make sparse hash set for constant membership check
      spp::sparse_hash_set<uint32_t> txpSet (txps.begin(), txps.end());
      for (auto& tg: tgroup){
        bool covered = false;
        for(auto txp: tg){
          if (txpSet.contains(txp)){
            covered = true;
            break;
          }
        }
        if (not covered){
          return false;
        }
      }
      return true;
    }

    // from here https://www.geeksforgeeks.org/print-subsets-given-size-set/
    void combinationUtil(std::vector<uint32_t>& arr, int n, int r,
                         int index, std::vector<uint32_t> data,
                         int i, std::vector<std::vector<uint32_t>>& comb) {
      // Current cobination is ready, print it
      if (index == r) {
        comb.emplace_back(data);
        return;
      }

      // When no more elements are there to put in data[]
      if (i >= n)
        return;

      // current is included, put next at next location
      data[index] = arr[i];
      combinationUtil(arr, n, r, index + 1, data, i + 1, comb);

      // current is excluded, replace it with next
      // (Note that i+1 is passed, but index is not
      // changed)
      combinationUtil(arr, n, r, index, data, i + 1, comb);
    }

    bool hasOneGene(const std::vector<uint32_t>& txps, uint32_t& geneId,
                    spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                    const size_t numGenes){
      spp::sparse_hash_set<uint32_t> geneids;
      for (auto& tid: txps){
        uint32_t gid;
        if(txpToGeneMap.contains(tid)){
          gid = txpToGeneMap.at(tid);
        }
        else{
          std::cerr << "Out of Range error for txp to gene Map: " << '\n' << std::flush;
          std::cerr << tid << "\t not found" << std::flush;
          exit(1);
        }
        geneids.insert(gid);
        if (geneids.size() > 1){
          return false;
        }
      }
      if (geneids.size() == 1){
        uint32_t gid = *geneids.begin();
        if (gid > numGenes){
          std::cerr<< "Gene id out of range"
                   << "Please check txp2gene has the write entries"
                   << std::flush;
          exit(1);
        }
        geneId = gid;
        return true;
      }
      return false;
    }

    template
    bool processAlevinOpts(AlevinOpts<apt::DropSeq>& aopt,
                           SalmonOpts& sopt,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::InDrop>& aopt,
                           SalmonOpts& sopt,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::Chromium>& aopt,
                           SalmonOpts& sopt,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::Gemcode>& aopt,
                           SalmonOpts& sopt,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::Custom>& aopt,
                           SalmonOpts& sopt,
                           boost::program_options::variables_map& vm);

    /*
    template bool sequenceCheck(std::string sequence,
                                AlevinOpts<apt::DropSeq>& aopt,
                                std::mutex& iomutex,
                                Sequence seqType);
    template bool sequenceCheck(std::string sequence,
                                AlevinOpts<apt::InDrop>& aopt,
                                std::mutex& iomutex,
                                Sequence seqType);
    template bool sequenceCheck(std::string sequence,
                                AlevinOpts<apt::Chromium>& aopt,
                                std::mutex& iomutex,
                                Sequence seqType);
    template bool sequenceCheck(std::string sequence,
                                AlevinOpts<apt::Gemcode>& aopt,
                                std::mutex& iomutex,
                                Sequence seqType);
    template bool sequenceCheck(std::string sequence,
                                AlevinOpts<apt::Custom>& aopt,
                                std::mutex& iomutex,
                                Sequence seqType);
    */
  }
}
