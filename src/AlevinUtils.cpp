#include "AlevinUtils.hpp"
#include "peglib.h"

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
    std::string* getReadSequence(apt::CITESeq& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      subseq.clear();
      subseq = seq2.substr(protocol.featureStart,
                          protocol.featureLength);
      return &subseq;
    }
    template <>
    std::string* getReadSequence(apt::DropSeq& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string* getReadSequence(apt::Chromium& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string*  getReadSequence(apt::ChromiumV3& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string*  getReadSequence(apt::CELSeq& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string*  getReadSequence(apt::CELSeq2& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string*  getReadSequence(apt::QuartzSeq2& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string*  getReadSequence(apt::Custom& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string*  getReadSequence(apt::CustomGeometry& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      bool ok = protocol.read_geo.extract_read(seq, seq2, subseq);
      return &subseq;
      //subseq = seq2;
    }
    template <>
    std::string*  getReadSequence(apt::Gemcode& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    template <>
    std::string*  getReadSequence(apt::InDrop& protocol,
                         std::string& seq,
                         std::string& seq2,
                         std::string& subseq){
      (void)seq;
      return &seq2;
    }
    // end of read extraction

    template <>
    bool extractUMI<apt::DropSeq>(std::string& read,
                                  std::string& read2,
                                  apt::DropSeq& pt,
                                  std::string& umi){
      (void)read2;
      return (read.length() >= pt.barcodeLength + pt.umiLength) ?
        (umi.assign(read, pt.barcodeLength, pt.umiLength), true) : false;
    }
    template <>
    bool extractUMI<apt::CITESeq>(std::string& read,
                                  std::string& read2,
                                  apt::CITESeq& pt,
                                  std::string& umi){
      (void)read2;
      return (read.length() >= pt.barcodeLength + pt.umiLength) ?
        (umi.assign(read, pt.barcodeLength, pt.umiLength), true) : false;
    }
    template <>
    bool extractUMI<apt::Chromium>(std::string& read,
                                   std::string& read2,
                                   apt::Chromium& pt,
                                   std::string& umi){
      (void)read2;
      return (read.length() >= pt.barcodeLength + pt.umiLength) ?
        (umi.assign(read, pt.barcodeLength, pt.umiLength), true) : false;
    }
    template <>
    bool extractUMI<apt::ChromiumV3>(std::string& read,
                                     std::string& read2,
                                     apt::ChromiumV3& pt,
                                     std::string& umi){
      (void)read2;
      return (read.length() >= pt.barcodeLength + pt.umiLength) ?
        (umi.assign(read, pt.barcodeLength, pt.umiLength), true) : false;
    }
    template <>
    bool extractUMI<apt::Gemcode>(std::string& read,
                                  std::string& read2,
                                   apt::Gemcode& pt,
                                   std::string& umi){
      (void)read2;
      return (read.length() >= pt.barcodeLength + pt.umiLength) ?
        (umi.assign(read, pt.barcodeLength, pt.umiLength), true) : false;
    }
    template <>
    bool extractUMI<apt::Custom>(std::string& read,
                                 std::string& read2,
                                 apt::Custom& pt,
                                 std::string& umi){
      (void)read2;
      if ( pt.end == BarcodeEnd::FIVE ) {
        umi.assign(read, pt.barcodeLength, pt.umiLength);
      } else if (pt.end == BarcodeEnd::THREE ) {
        umi.assign(read, 0, pt.umiLength);
      } else {
        return false;
      }
      return true;
    }
    template <>
    bool extractUMI<apt::CustomGeometry>(std::string& read1,
                                 std::string& read2,
                                 apt::CustomGeometry& pt,
                                 std::string& umi){
      
      return pt.umi_geo.extract_tag(read1, read2, umi);
    }
    template <>
    bool extractUMI<apt::QuartzSeq2>(std::string& read,
                                     std::string& read2,
                                     apt::QuartzSeq2& pt,
                                     std::string& umi){
      (void)read2;
      return (read.length() >= pt.barcodeLength + pt.umiLength) ?
        (umi.assign(read, pt.barcodeLength, pt.umiLength), true) : false;
      return true;
    }
    template <>
    bool extractUMI<apt::CELSeq2>(std::string& read,
                                  std::string& read2,
                                  apt::CELSeq2& pt,
                                  std::string& umi){
      (void)read2;
      return (read.length() >= pt.umiLength) ?
        (umi.assign(read, 0, pt.umiLength), true) : false;
    }
    template <>
    bool extractUMI<apt::CELSeq>(std::string& read,
                                 std::string& read2,
                                 apt::CELSeq& pt,
                                 std::string& umi){
      (void)read2;
      return (read.length() >= pt.umiLength) ?
        (umi.assign(read, 0, pt.umiLength), true) : false;
      return true;
    }
    template <>
    bool extractUMI<apt::InDrop>(std::string& read,
                                 std::string& read2,
                                 apt::InDrop& pt,
                                 std::string& umi){
      (void)read;
      (void)read2;
      std::cout<<"Incorrect call for umi extract";
      exit(1);
    }

    template <>
    bool extractBarcode<apt::DropSeq>(std::string& read,
                                      std::string& read2,
                                      apt::DropSeq& pt,
                                      std::string& bc){
      (void)read2;
      return (read.length() >= pt.barcodeLength) ?
        (bc.assign(read, 0, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::CITESeq>(std::string& read,
                                      std::string& read2,
                                                               apt::CITESeq& pt,
                                                               std::string& bc){
      (void)read2;
      return (read.length() >= pt.barcodeLength) ?
        (bc.assign(read, 0, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::ChromiumV3>(std::string& read,
                                         std::string& read2,
                                                                  apt::ChromiumV3& pt,
                                                                  std::string& bc){
      (void)read2;
      return (read.length() >= pt.barcodeLength) ?
        (bc.assign(read,0, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::Chromium>(std::string& read,
                                       std::string& read2,
                                       apt::Chromium& pt,
                                       std::string& bc){
      (void)read2;
      return (read.length() >= pt.barcodeLength) ?
        (bc.assign(read, 0, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::Gemcode>(std::string& read,
                                      std::string& read2,
                                      apt::Gemcode& pt,
                                      std::string& bc){
      (void)read2;
      return (read.length() >= pt.barcodeLength) ?
        (bc.assign(read, 0, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::Custom>(std::string& read,
                                     std::string& read2,
                                     apt::Custom& pt,
                                     std::string& bc){
      (void)read2;
      if (pt.end == BarcodeEnd::FIVE) {
        return (read.length() >= pt.barcodeLength) ?
          (bc.assign(read, 0, pt.barcodeLength), true) : false;
      } else if (pt.end == BarcodeEnd::THREE) {
        return (read.length() >= (pt.umiLength + pt.barcodeLength)) ?
          (bc.assign(read, pt.umiLength, pt.barcodeLength), true) : false;
      } else {
        return false; 
      }
    }
    template <>
    bool extractBarcode<apt::CustomGeometry>(std::string& read1,
                                            std::string& read2,
                                     apt::CustomGeometry& pt,
                                     std::string& bc){
      return pt.bc_geo.extract_tag(read1, read2, bc);
    }
    template <>
    bool extractBarcode<apt::QuartzSeq2>(std::string& read,
                                         std::string& read2,
                                         apt::QuartzSeq2& pt,
                                         std::string& bc){
      (void)read2;
      return (read.length() >= pt.barcodeLength) ?
        (bc.assign(read, 0, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::CELSeq2>(std::string& read,
                                      std::string& read2,
                                                               apt::CELSeq2& pt,
                                                               std::string& bc){
      (void)read2;
      return (read.length() >= (pt.umiLength + pt.barcodeLength)) ?
        (bc.assign(read, pt.umiLength, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::CELSeq>(std::string& read,
                                     std::string& read2,
                                                              apt::CELSeq& pt,
                                                              std::string& bc){
      (void)read2;
      return (read.length() >= (pt.umiLength + pt.barcodeLength)) ?
        (bc.assign(read, pt.umiLength, pt.barcodeLength), true) : false;
    }
    template <>
    bool extractBarcode<apt::InDrop>(std::string& read, 
                                     std::string& read2, 
                                     apt::InDrop& pt, std::string& bc){
      (void)read2;
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
        exit(64);
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

    void getTxpToGeneMap(spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                         spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                         const std::string& t2gFileName,
                         const std::string& refNamesFile,
                         const std::string& refLengthFile,
                         const std::string& headerFile,
                         std::shared_ptr<spdlog::logger>& jointLog,
                         bool noTgMap){
      size_t  kSize, numDupTxps{0};
      uint64_t numberOfDecoys, firstDecoyIndex;
      spp::sparse_hash_map<std::string, std::vector<uint32_t>> txpIdxMap;
      // reading in the binary file
      std::vector<std::string> txpNames;
      {
        std::ifstream ctstream(refNamesFile);
        cereal::BinaryInputArchive contigTableArchive(ctstream);
        contigTableArchive(txpNames);
        ctstream.close();
      }

      std::vector<uint32_t> txpLengths;
      {
        std::ifstream ctstream(refLengthFile);
        cereal::BinaryInputArchive contigTableArchive(ctstream);
        contigTableArchive(txpLengths);
        ctstream.close();
      }

      {
        std::ifstream hstream(headerFile);
        cereal::JSONInputArchive headerArchive(hstream);
        headerArchive( cereal::make_nvp("num_decoys", numberOfDecoys) );
        headerArchive( cereal::make_nvp("first_decoy_index", firstDecoyIndex) );
        headerArchive( cereal::make_nvp("k", kSize) );
        hstream.close();
      }

      size_t numShort {0};
      {
        // kk this is tricky.
        // firstDecoyIndex is the index of the first decoy *after*
        // removing short transcripts but the reported mappings are
        // in full vector i.e. including small transcripts
        size_t i {0};
        if (numberOfDecoys > 0) {
          while ( i-numShort < firstDecoyIndex ) {
            if (txpLengths[i] <= kSize) { numShort += 1; }
            txpIdxMap[txpNames[i]].emplace_back(i);
            i += 1;
          }
        } else {
          for (size_t i=0; i < txpNames.size(); i++) {
            if (txpLengths[i] <= kSize) { numShort += 1; }
            txpIdxMap[txpNames[i]].emplace_back(i);
          } // end-for
        } // end else
      } // end block for populating txpIdxMap

      for (auto it: txpIdxMap) {
        size_t bucketLen = it.second.size();
        if (bucketLen > 1) { numDupTxps += (bucketLen - 1); }
      }

      jointLog->info("Found {} transcripts(+{} decoys, +{} short and +{}"
                     " duplicate names in the index)",
                     txpNames.size() - numberOfDecoys - numShort - numDupTxps ,
                     numberOfDecoys, numShort, numDupTxps);

      if (noTgMap){ // mapping the protein name to itself
        for (size_t i=0; i<txpNames.size(); i++) {
          geneIdxMap[txpNames[i]] = i;
          txpToGeneMap[i] = i;
        }
      } else {
        std::ifstream t2gFile(t2gFileName);
        uint32_t gid, geneCount{0};
        std::vector<uint32_t> tids;
        std::string tStr, gStr;
        if(t2gFile.is_open()) {
          while( not t2gFile.eof() ) {
            t2gFile >> tStr >> gStr;

            if(not txpIdxMap.contains(tStr)  ){
              continue;
            }

            tids = txpIdxMap[tStr];
            if (geneIdxMap.contains(gStr)){
              gid = geneIdxMap[gStr];
            }
            else{
              gid = geneCount;
              geneIdxMap[gStr] = gid;
              geneCount++;
            }

            for (auto tid: tids) {
              if (txpToGeneMap.find(tid) != txpToGeneMap.end() &&
                  txpToGeneMap[tid] != gid ) {
                jointLog->warn("Dual txp to gene map for txp {}", txpNames[tid]);
              }
              txpToGeneMap[tid] = gid;
            }
          }
          t2gFile.close();
        }
      } // done parsing txp to gene map file

      jointLog->info( "Filled with {} txp to gene entries ", txpToGeneMap.size());
      for ( auto it: txpIdxMap ) {
        for (auto tid: it.second) {
          if (tid < firstDecoyIndex + numShort &&
              txpToGeneMap.find( tid ) == txpToGeneMap.end() ) {
            jointLog->error( "ERROR: Can't find gene mapping for : {} w/ index {}",
                             it.first, tid );
            jointLog->error( "ERROR: "
                             "Txp to Gene Map not found for {}"
                             " transcripts. Exiting",
                             txpIdxMap.size() - txpToGeneMap.size() - numberOfDecoys
                             );
            jointLog->flush();
            exit(1);
          }
        }
      }

      jointLog->info("Found all transcripts to gene mappings");
    }


    struct Desc {
      uint8_t read_num;
      std::vector<std::pair<size_t,size_t>> locs;
      friend std::ostream& operator<<(std::ostream& os, const Desc& de);
    };
    
    std::ostream& operator<<(std::ostream& os, const Desc& de) {
      os << "Geometry description :: [\n";
      os << "\tRead Num : " << static_cast<int32_t>(de.read_num) << "\n\t[";
      for (size_t i = 0; i < de.locs.size(); ++i) {
        os << "  (" << de.locs[i].first << ", " << de.locs[i].second << ")";
      }
      os << "\t]\n";
      os << "]\n";
      return os;
    }
 

    bool parse_geometry_desc(std::string& geo_string, bool can_be_unbounded, peg::parser& parser, std::vector<Desc>& d, std::shared_ptr<spdlog::logger> log, alevin::protocols::TagGeometry& tg) {
      // the first element in the string must be 
      // the read number (currently, only 1 or 2 is )
      // acceptable.

      // currently, only the read geometry can be unbounded 
      // (i.e. can stretch to std::string::npos).  If 
      // can_be_unbounded is false, then reject any parsed
      // geometry tag containing `end`.

      d.clear();
      if ( parser.parse(geo_string.c_str()) ) {
        /*
        if (d.size() > 1) {
          log->error("Though supported in the syntax, the current implementation \n"
                     "of custom tag geometry does not support having a tag \n"
                     "split over more than one read.");    
          return false;
        }
        */

        //auto desc = d[0];
        //tg.read_num = desc.read_num;
        for (auto& desc : d) {
          auto rn = desc.read_num;
          auto& tg_length = (rn == 1) ? tg.length1 : tg.length2;
          auto& tg_largest_index = (rn == 1) ? tg.largest_index1 : tg.largest_index2;
          auto& tg_substr_locs = (rn ==1) ? tg.substr_locs1 : tg.substr_locs2;

          // the length of the tag
          size_t tlen{0};
        
          // the largest index the tag reaches
          size_t largest_index{0};

          for (auto& start_stop : desc.locs) {
            size_t start = start_stop.first;
            if (start <= 0) { log->error("tag string subscript must be strictly positive (tag geometry is 1-indexed)."); return false; }
            size_t stop = start_stop.second;
            // make sure we don't contain `end` if we can't
            if (!can_be_unbounded and stop == std::string::npos) {
              log->error("only the geometry of the read can be unbounded (can contain 'end'), "
                         "but it was specified for a cell barcode / umi.");
              return false;
            }

            // does this piece reach to the end of the read
            // i.e. is it a range of the form `X-end`?
            bool reach_to_end = (stop == std::string::npos);

            if (stop < start) { 
              log->error("geometry contained a range {}-{}; cannot have stop < start.", start, stop);
              return false; 
            }

            // internally, we start indexing from 0
            start -= 1;
            if (!reach_to_end) { stop -= 1; }

            largest_index = (stop > largest_index) ? stop : largest_index;

            size_t len = (reach_to_end) ? std::string::npos : (stop - start) + 1;

            // if the current piece is unbounded, or any piece so far
            // has been unbounded, then the tag length is unbounded
            tlen = (reach_to_end or (tlen == std::string::npos)) ? len : tlen + len;

            if (tg_substr_locs.size() >= alevin::protocols::num_tag_pieces) {
              log->error("Currently, alevin does not support the tag (barcode / umi) being "
                         "split into more than {} pieces.  If the current bound is a "
                         "problem for your protocol, please reach out on GitHub.",
                         alevin::protocols::num_tag_pieces);
              return false;
            }
            tg_substr_locs.push_back(std::make_pair(start, len));
          }
          tg_length = tlen;
          tg_largest_index = largest_index;
        } // for (auto& desc : d) 
        return true;
      } else {
        log->error("parser failure!");
        return false;
      }
      return false;
    }

    template <typename ProtocolT>
    bool processAlevinOpts(AlevinOpts<ProtocolT>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm){
      namespace bfs = boost::filesystem;
      namespace po = boost::program_options;

      // mark in salmon options that we are running
      // in alevin mode
      sopt.alevinMode = true;
      sopt.hardFilter = true;
      if (sopt.initUniform) { aopt.initUniform = true; }

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

      if (vm.count("vbemPrior")){
        aopt.vbemPriorFile = vm["vbemPrior"].as<std::string>();
        aopt.useVBEM = true;

        aopt.vbemNorm = vm["vbemNorm"].as<double>();
        if (aopt.vbemNorm == 0.0) {
          fmt::print(stderr,"\nVBEM Normalization Factor not provided");
          return false;
        }

        if (!bfs::exists(aopt.vbemPriorFile)) {
          fmt::print(stderr,"\nVBEM Prior File {} does not exists\n Exiting Now",
                     aopt.vbemPriorFile.string());
          return false;
        }
      }

      if (vm.count("hash")){
        aopt.bfhFile = vm["hash"].as<std::string>();
        if (!bfs::exists(aopt.bfhFile)) {
          fmt::print(stderr,"\nBfh File {} does not exists\n Exiting Now",
                     aopt.bfhFile.string());
          return false;
        }
      }

      if (not noTgMap) {
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
      }

      //create logger
      spdlog::set_async_mode(131072);
      auto logPath = aopt.outputDirectory / "alevin.log";
      auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
      auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>();
      consoleSink->set_color(spdlog::level::warn, consoleSink->magenta);
      std::vector<spdlog::sink_ptr> sinks{consoleSink, fileSink};
      aopt.jointLog = spdlog::create("alevinLog", std::begin(sinks), std::end(sinks));

      aopt.just_align = vm["rad"].as<bool>();
      aopt.sketch_mode = vm["sketch"].as<bool>();
      if (aopt.sketch_mode and !aopt.just_align) {
        aopt.jointLog->info("currently, --sketch implies --rad. Running in "
                            "alignment-only mode (will write a RAD output).");
        aopt.just_align = true;
      }

      aopt.quiet = vm["quiet"].as<bool>();
      aopt.noEM = vm["noem"].as<bool>();
      aopt.noDedup = vm["noDedup"].as<bool>();
      aopt.noWhitelist = vm["noWhitelist"].as<bool>();
      aopt.naiveEqclass = vm["naiveEqclass"].as<bool>();
      aopt.noQuant = vm["noQuant"].as<bool>();
      aopt.dumpfq = vm["dumpfq"].as<bool>();
      aopt.dumpArborescences = vm["dumpArborescences"].as<bool>();
      aopt.dumpfeatures = vm["dumpFeatures"].as<bool>();
      aopt.dumpMtx = vm["dumpMtx"].as<bool>();
      aopt.dumpBarcodeEq = vm["dumpBarcodeEq"].as<bool>();
      aopt.dumpBFH = vm["dumpBfh"].as<bool>();
      aopt.dumpUmiGraph = vm["dumpUmiGraph"].as<bool>();
      aopt.dumpCellEq = vm["dumpCellEq"].as<bool>();
      aopt.trimRight = vm["trimRight"].as<uint32_t>();
      aopt.numBootstraps = vm["numCellBootstraps"].as<uint32_t>();
      aopt.numGibbsSamples = vm["numCellGibbsSamples"].as<uint32_t>();
      aopt.lowRegionMinNumBarcodes = vm["lowRegionMinNumBarcodes"].as<uint32_t>();
      aopt.maxNumBarcodes = vm["maxNumBarcodes"].as<uint32_t>();
      aopt.freqThreshold = vm["freqThreshold"].as<uint32_t>();
      aopt.umiEditDistance = vm["umiEditDistance"].as<uint32_t>();
      aopt.forceCells = vm["forceCells"].as<uint32_t>();
      aopt.expectCells = vm["expectCells"].as<uint32_t>();

      if (aopt.just_align) {
        aopt.jointLog->info("The --rad flag was passed to alevin. The "
        "reads will be selectively aligned and the output written to a RAD file."
        "Arguments passed that correspond to other processing steps will be ignored");
        if (aopt.sketch_mode) {
          aopt.jointLog->info("The --sketch flag was passed; the alignment will be run "
          "in sketch mode.");
          sopt.mismatchSeedSkip = 7;
        }
      } else {
        if (aopt.sketch_mode) {
          aopt.jointLog->info("The --sketch flag is not meaningful without the "
          "--rad flag.  This flag will be ignored.");
        }
      }

      if (aopt.umiEditDistance > 4 ) {
        aopt.jointLog->error("Too high edit distance collapsing {}, expected <= 4",
                             aopt.umiEditDistance);
        return false;
      }

      if(vm.count("iupac")){
        aopt.iupac = vm["iupac"].as<std::string>();
      }

      if (sopt.numBootstraps>0) {
        aopt.jointLog->error("Do you mean numCellBootstraps ?");
        return false;
      }

      if (aopt.numGibbsSamples > 0 and aopt.numBootstraps > 0) {
        aopt.jointLog->error("Either of --numCellGibbsSamples or --numCellBootstraps "
                             "can be used");
        return false;
      }

      if ( aopt.numBootstraps > 0 and aopt.noEM ) {
        aopt.jointLog->error("cannot perform bootstrapping with noEM option.");
        return false;
      }

      aopt.keepCBFraction = vm["keepCBFraction"].as<double>();
      if ( aopt.keepCBFraction > 0.0 ) {
        if ( vm.count("whitelist") ) {
          aopt.jointLog->error("keepCBFraction and whitelist cannot be used together");
          aopt.jointLog->flush();
          exit(1);
        }

        aopt.jointLog->warn("Force Cells to {} fraction of All possible CB."
                            "This is not recommended way to run the pipeline,"
                            "and it might slow the pipeline",
                            aopt.keepCBFraction);
      } else if (noTgMap) { aopt.keepCBFraction = 1.0;  }

      if (not vm.count("threads")) {
        auto tot_cores = std::thread::hardware_concurrency();
        aopt.numThreads = std::max(1, static_cast<int>(tot_cores/4.0));
        aopt.jointLog->warn("threads flag not specified, Using {}(25%) of total cores",
                       aopt.numThreads);
      }
      else{
        aopt.numThreads = vm["threads"].as<uint32_t>();
      }  // things which needs to be updated for salmonOpts

      // handling of customized barcode and umi geometry 
      bool have_custom_umi_geo = vm.count("umi-geometry");
      bool have_custom_bc_geo = vm.count("bc-geometry");
      bool have_custom_read_geo = vm.count("read-geometry");
      // need both
      bool have_any_custom_geo = have_custom_read_geo or have_custom_umi_geo or have_custom_bc_geo;
      bool have_all_custom_geo = have_custom_read_geo and have_custom_umi_geo and have_custom_bc_geo;
      if ( have_any_custom_geo and !have_all_custom_geo ) {
        aopt.jointLog->error("If you are using either of the umi-geometry or \n"
                             "the barcode-geometry options, then you have to provide both.\n"
                             "Alternatively, you can use pre-defined single-cell protocol flags.\n"
                             "Exiting Now.");
        return false; 
      }

      // validate customized options for custom protocol
      bool haveCustomEnd = vm.count("end");
      bool haveCustomBC= vm.count("barcodeLength");
      bool haveCustomUMI = vm.count("umiLength");

      bool allCustom = (haveCustomEnd and haveCustomBC and haveCustomUMI);
      bool noCustom = !(haveCustomEnd or haveCustomBC or haveCustomUMI);

      if (!noCustom and have_any_custom_geo) {
        aopt.jointLog->warn("Note: the use of --end, --barcodeLength and --umiLength \n"
                            "to describe the barcode and umi geometry incompatible \n"
                            "with the new options `--barcode-geometry` and `--umi-geometry`.\n"
                            "The former are deprecated, please adopt the new options instead\n.");  
        return false;
      }

      // These are all or nothing.  Either the user must provide all 3
      // or none of these options.
      if (!(noCustom or allCustom)) {
        aopt.jointLog->error("If you are using any one"
                             " of (end, umilength, barcodelength) flag\n"
                             "you have to provide all of them explicitly.\n"
                             "You can also use pre-defined single-cell protocols."
                             "Exiting Now.");
        return false;
      }

      if (allCustom) {
        uint32_t barEnd = vm["end"].as<uint32_t>();
        uint32_t barcodeLength = vm["barcodeLength"].as<uint32_t>();
        uint32_t umiLength = vm["umiLength"].as<uint32_t>();

        aopt.jointLog->warn("Note: the use of --end, --barcodeLength and --umiLength "
                            "to describe the barcode and umi geometry is deprecated. "
                            "Please start using the `--barcode-geometry` and `--umi-geometry` "
                            "options instead.");  
        // validate that BC and UMI lengths are OK
        uint32_t maxBC{20};
        uint32_t maxUMI{20};
        if (barcodeLength < 1 or barcodeLength > maxBC) {
          aopt.jointLog->error("Barcode length ({}) was not in the required length range [1, {}].\n"
                               "Exiting now.", barcodeLength, maxBC);
          return false;
        }
        if (umiLength < 1 or umiLength > maxUMI) {
          aopt.jointLog->error("UMI length ({}) was not in the required length range [1, {}].\n"
                               "Exiting now.", umiLength, maxUMI);
          return false;
        }

        // validate that protocol end is OK and set it
        if (barEnd == 3) {
          aopt.protocol.end = BarcodeEnd::THREE;
        }
        else if (barEnd == 5) {
          aopt.protocol.end = BarcodeEnd::FIVE;
        } else{
          aopt.jointLog->error("Wrong value for Barcode-end of read -> {}.\n"
                               "Please provide `5` for barcodes "
                               "starting at 5' end or `3` for barcode starting "
                               "at 3' end.\nExiting now.", barEnd);
          return false;
        }

        // If all validation passed, then set the appropriate variables.
        aopt.protocol.barcodeLength = barcodeLength;
        aopt.protocol.umiLength = umiLength;
        // Since we passed a custom UMI length and need to update the value here.
        if (umiLength != alevin::types::AlevinUMIKmer::k()) {
          alevin::types::AlevinUMIKmer::k(umiLength);
          aopt.jointLog->info("A custom protocol (END, BC length, UMI length) = ({}, {}, {}) "
                              "is being used.  Updating UMI k-mer length accordingly.",
                              barEnd, barcodeLength, umiLength);
        }
      } else if (have_all_custom_geo) {
        std::string  bc_geo_string = vm["bc-geometry"].as<std::string>();
        std::string  umi_geo_string = vm["umi-geometry"].as<std::string>();
        std::string  read_geo_string = vm["read-geometry"].as<std::string>();

        // This describes the parsing expression grammar that 
        // we use for describing barcode geometry.
        peg::parser parser(R"(
         DescriptionList <- Description (','Description){0,1}
         Description <- ReadNumber'['NumberRangeList']'
         ReadNumber  <- [1,2]
         Number      <- [0-9]+
         End         <- 'end'
         NumberRange <- (Number '-' Number) / (Number '-' End)
         NumberRangeList <- NumberRange (','NumberRange)*
        )");
        
        if (!(bool)parser) {
          aopt.jointLog->error("Failed to instantiate the tag geometry parser.\n"
                               "This should not happen. Please report this issue on GitHub.");
          return false;
        }

        // The variable to hold the temporary parsing information
        std::vector<Desc> d;
        
        parser["ReadNumber"] = [&](const peg::SemanticValues& sv) { 
          d.push_back({0,{}}); d.back().read_num = std::stoi(sv.token()); 
        };
        
        parser["Number"] = [&](const peg::SemanticValues& sv) { 
          return static_cast<size_t>(std::stoull(sv.token())); 
        };

        parser["NumberRange"] = [&](const peg::SemanticValues& sv) {
          switch (sv.choice()) {
          case 0:
           {
            d.back().locs.push_back(std::make_pair(
                peg::any_cast<size_t>(sv[0]), peg::any_cast<size_t>(sv[1])));
           }
            break;
          default:
            {
             d.back().locs.push_back(std::make_pair(
               peg::any_cast<size_t>(sv[0]), std::string::npos));
            }
            break;
          }
        };

        // NOTE: this is just for backwards-compatibility.
        // The barcode end is redundant with the new geometry
        // specification.
        aopt.protocol.end = BarcodeEnd::FIVE;

        // parse the cellular barcode geometry
        alevin::protocols::TagGeometry bc_geo;
        if ( parse_geometry_desc(bc_geo_string, false, parser, d, aopt.jointLog, bc_geo) ) {
          // if we are *not* in RAD mode, then the barcode cannot span beyond 
          // read 1
          if ( (!aopt.just_align) and bc_geo.uses_r2() ) {
            aopt.jointLog->error("Currently, barcodes spanning read 1 and read 2 "
                                 "are supported in alignment/sketch mode only "
                                 "(i.e. with the --rad flag or the --sketch flag). "
                                 "Custom barcode geometry without this mode must reside entirely on read1. "
                                 "If you require custom barcode geometry spanning both reads, consider "
                                 "using the alevin-fry pipeline.");
            return false;
          }
          aopt.protocol.set_bc_geo(bc_geo);
        } else {
          aopt.jointLog->error("Failed to parse `--bc-geometry` argument {}, please make sure it is correct.", bc_geo_string);
          return false;
        }
        //std::cerr << "BC GEO\n---------\n";
        //std::cerr << bc_geo << "\n\n";

        // parse the UMI geometry
        alevin::protocols::TagGeometry umi_geo;
        if ( parse_geometry_desc(umi_geo_string, false, parser, d, aopt.jointLog, umi_geo) ) { 
          // if we are *not* in RAD mode, then the UMI cannot span beyond 
          // read 1
          if ( (!aopt.just_align) and bc_geo.uses_r2() ) {
            aopt.jointLog->error("Currently, umis spanning read 1 and read 2 "
                                 "are supported in alignment/sketch mode only "
                                 "(i.e. with the --rad flag or the --sketch flag). "
                                 "Custom umi geometry without this mode must reside entirely on read1. "
                                 "If you require custom umi geometry spanning both reads, consider "
                                 "using the alevin-fry pipeline.");
            return false;
          }
          aopt.protocol.set_umi_geo(umi_geo);
        } else {
          aopt.jointLog->error("Failed to parse `--umi-geometry` argument {}, please make sure it is correct.", umi_geo_string);
        }
        //std::cerr << "UMI GEO\n---------\n";
        //std::cerr << umi_geo << "\n\n";

        // parse the read geometry
        alevin::protocols::TagGeometry read_geo;
        if ( parse_geometry_desc(read_geo_string, true, parser, d, aopt.jointLog, read_geo) ) {
          aopt.protocol.set_read_geo(read_geo);
        } else {
          aopt.jointLog->error("Failed to parse `--read-geometry` argument {}, please make sure it is correct.", read_geo_string);
        }

        // validate that BC and UMI lengths are OK
        uint32_t maxBC{31};
        uint32_t maxUMI{31};
        // the barcode length must be in [1,31]
        if ((bc_geo.length() < 1) or (bc_geo.length() > maxBC)) {
          aopt.jointLog->error("Barcode length ({}) was not in the required length range [1, {}].\n"
                               "Exiting now.", bc_geo.length(), maxBC);
          return false;
        }
        // if it's OK, set the barcode kmer length
        alevin::types::AlevinCellBarcodeKmer::k( static_cast<uint16_t>(bc_geo.length()) );

        // the UMI length must be in [1,31]
        if ((umi_geo.length() < 1) or (umi_geo.length() > maxUMI)) {
          aopt.jointLog->error("UMI length ({}) was not in the required length range [1, {}].\n"
                               "Exiting now.", umi_geo.length(), maxUMI);
          return false;
        }
        // if it's OK, set the umi kmer length
        alevin::types::AlevinUMIKmer::k( static_cast<uint16_t>(umi_geo.length()) );

      } // new custom barcode geometry

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
            aopt.jointLog->error("\nERROR: Wrong IUPAC character {} in {}\n"
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

      // enable validate mappings
      sopt.validateMappings = true;
      sopt.hitFilterPolicyStr = "BOTH";
      bool optionsOK =
        salmon::utils::processQuantOptions(sopt, vm, vm["numBiasSamples"].as<int32_t>());
      if (!vm.count("minScoreFraction")) {
        if (noTgMap) {
          int32_t l = vm["featureLength"].as<size_t>();
          sopt.minScoreFraction = salmon::utils::compute_1_edit_threshold(l, sopt);
          aopt.jointLog->info("set CITE-seq minScoreFraction parameter to : {}",
                              sopt.minScoreFraction);
        } else {
          sopt.minScoreFraction = alevin::defaults::minScoreFraction;
        }

        sopt.consensusSlack = alevin::defaults::consensusSlack;
        sopt.jointLog->info(
            "Using default value of {} for minScoreFraction in Alevin\n"
            "Using default value of {} for consensusSlack in Alevin",
            sopt.minScoreFraction, sopt.consensusSlack);
      }

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
      return (sequence.length() > 0) and (sequence.find_first_not_of("ACGTacgt") == std::string::npos);
    }

    bool recoverBarcode(std::string& sequence){
      size_t pos = sequence.find_first_of("Nn");
      if (pos == std::string::npos) { return false; }

      // Randomly assigning 'A' to first base with 'N'
      sequence[pos] = 'A';
      return sequenceCheck(sequence);
    }

    void readWhitelist(bfs::path& filePath,
                       TrueBcsT& trueBarcodes) {
      std::ifstream whiteFile(filePath.string());
      std::string whtBc;
      if(whiteFile.is_open()) {
        while(getline(whiteFile, whtBc)) {
          trueBarcodes.insert(whtBc);
        }
        whiteFile.close();
      }
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
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::CITESeq>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::InDrop>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::ChromiumV3>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::Chromium>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::Gemcode>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::Custom>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::CustomGeometry>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::CELSeq>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::CELSeq2>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
    template
    bool processAlevinOpts(AlevinOpts<apt::QuartzSeq2>& aopt,
                           SalmonOpts& sopt, bool noTgMap,
                           boost::program_options::variables_map& vm);
  }
}
