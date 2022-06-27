#ifndef __SINGLE_CELL_PROTOCOLS_HPP__
#define __SINGLE_CELL_PROTOCOLS_HPP__

#include <string>
#include <ostream>

#include "AlevinOpts.hpp"
#include "AlevinTypes.hpp"
#include "pufferfish/itlib/static_vector.hpp"

namespace alevin{
  namespace protocols {

    static constexpr size_t num_tag_pieces{16};
    struct TagGeometry {
      // uint32_t read_num{0};
      // tuples are read_num, start_pos, length
      itlib::static_vector<std::pair<uint32_t, size_t>, num_tag_pieces> substr_locs1{};
      itlib::static_vector<std::pair<uint32_t, size_t>, num_tag_pieces> substr_locs2{};
      // the total length of the tag on read 1 
      size_t length1{0};
      // the total length of the tag on read 2
      size_t length2{0};
      // the largest index on read 1
      size_t largest_index1{0};
      // the largest index on read 2
      size_t largest_index2{0};

      inline bool unbounded1() const { return length1 == std::string::npos; }
      inline bool unbounded2() const { return length2 == std::string::npos; }

      inline bool uses_r1() const { return !substr_locs1.empty(); }
      inline bool uses_r2() const { return !substr_locs2.empty(); }

      size_t length() const { return length1 + length2; }

      // Given the geometry of the tag, fill in the tag from 
      // read 1 (`from1`) and read 2 (`from2`), placing the constructed
      // tag in `to`.
      //
      // *assumption*: `to` is large enough to hold the tag
      // *returns*: true if the tag was written completely, and false otherwise
      inline bool extract_tag(std::string& from1, std::string& from2, std::string&to) {
        // if anything is too short, just ignore the whole thing
        if (uses_r1() and (from1.length() < largest_index1)) { return false; }
        if (uses_r2() and (from2.length() < largest_index2)) { return false; }

        // will point to the next place to 
        // begin filling the output string
        auto fill_it = to.begin();
        
        // grab anything from read 1
        auto f1b = from1.begin();
        for (auto& st_len : substr_locs1) {
          auto f1 = f1b + st_len.first;
          fill_it = std::copy(f1, f1 + st_len.second, fill_it);
        }
        
        // grab anything from read 2
        auto f2b = from2.begin();
        for (auto& st_len : substr_locs2) {
          auto f2 = f2b + st_len.first;
          fill_it = std::copy(f2, f2 + st_len.second, fill_it);
        }
        return true;
      }

      inline bool extract_read(std::string& from1, std::string& from2, std::string&to) {
        // if anything is too short, just ignore the whole thing
        if (uses_r1() and !unbounded1() and (from1.length() < largest_index1)) { return false; }
        if (uses_r2() and !unbounded2() and (from2.length() < largest_index2)) { return false; }
        
        // since the read extraction doesn't have a 
        // fixed size, we'll append rather than 
        // overwrite.
        to.clear();
        // grab anything from read 1
        for (auto& st_len : substr_locs1) {
          to.append(from1, st_len.first, st_len.second);
        }
        // grab anything from read 2
        for (auto& st_len : substr_locs2) {
          to.append(from2, st_len.first, st_len.second);
        }
        return true;
      }

    };

    std::ostream& operator<<(std::ostream& os, const TagGeometry& tg);

    struct Rule{
      Rule(){}
      Rule(uint32_t barcodeLength_,
           uint32_t umiLength_,
           BarcodeEnd end_,
           uint32_t maxValue_,
           ReadsToUse readsToUse_):
        barcodeLength(barcodeLength_),
        umiLength(umiLength_),
        end(end_),
        maxValue(maxValue_),
        readsToUse(readsToUse_) {
        alevin::types::AlevinUMIKmer::k(umiLength);
      }
      // NOTE: these functions are duplicated 
      // with those in `CustomGeometry` below, and 
      // due to semantics have slightly different 
      // implementations. See if the design can be 
      // unified later.
      void set_umi_geo(TagGeometry& g) { umiLength = g.length(); };
      void set_bc_geo(TagGeometry& g) { barcodeLength = g.length(); };
      void set_read_geo(TagGeometry& g) { (void)g; };
      uint32_t barcode_length() const { return barcodeLength; }
      uint32_t umi_length() const { return umiLength; }
      ReadsToUse get_reads_to_use() const { return readsToUse; }

      uint32_t barcodeLength, umiLength, maxValue;
      BarcodeEnd end;
      ReadsToUse readsToUse;
    };


    struct DropSeq : Rule{
      //Drop-Seq starts from 5 end with 12 length
      //barcode and 8 length umi & iupac can be
      //changed
      DropSeq(): Rule(12, 8, BarcodeEnd::FIVE, 16777216, ReadsToUse::USE_SECOND){}
    };

    struct InDropV2 : Rule{
        //InDropV2 starts from 5end with variable
        //length barcodes where barcode1 varies from 8 to 11 bp
        // followed by w1 sequence, 8 bp barcode2 and 6bp UMI
      InDropV2(): Rule(20, 6, BarcodeEnd::FIVE, 22347776, ReadsToUse::USE_FIRST){}

      std::string w1;
      std::size_t w1Length, maxHammingDist = 2, bc2Len = 8;
      void setW1(std::string& w1_){
        w1 = w1_;
        w1Length = w1.length();
      }
      std::size_t w1Pos = 0, bc2EndPos;
    };

    struct CITESeq : Rule{
      CITESeq(): Rule(16, 10, BarcodeEnd::FIVE, 4294967295, ReadsToUse::USE_SECOND){
        featureLength = 15;
        featureStart = 10;
      }

      size_t featureLength, featureStart;
      void setFeatureLength(size_t length) { featureLength = length; }
      void setFeatureStart(size_t startIdx) { featureStart = startIdx; }
    };

    struct Chromium5V2 : Rule{
      // fix barcodeLength and umiLength
      Chromium5V2(): Rule(16, 12, BarcodeEnd::FIVE, 4294967295, ReadsToUse::USE_BOTH){}
    };

    struct ChromiumV3 : Rule{
      ChromiumV3(): Rule(16, 12, BarcodeEnd::FIVE, 4294967295, ReadsToUse::USE_SECOND){}
    };

    struct Chromium : Rule{
      Chromium(): Rule(16, 10, BarcodeEnd::FIVE, 4294967295, ReadsToUse::USE_SECOND){}
    };

    struct Gemcode : Rule{
      Gemcode(): Rule(14, 10, BarcodeEnd::FIVE, 268435456, ReadsToUse::USE_SECOND){}
    };

    struct QuartzSeq2 : Rule{
      QuartzSeq2(): Rule(15, 8, BarcodeEnd::FIVE, 1073741824, ReadsToUse::USE_SECOND){}
    };

    struct CELSeq : Rule{
      // WEHI SCORE's CEL-Seq2 starts from 5' end with a 8 bp barcode
      // and a 6 bp UMI.
      CELSeq(): Rule(8, 6, BarcodeEnd::FIVE, 65536, ReadsToUse::USE_SECOND){}
    };
    struct CELSeq2 : Rule{
      // WEHI SCORE's CEL-Seq2 starts from 5' end with a 8 bp barcode
      // and a 6 bp UMI.
      CELSeq2(): Rule(6, 6, BarcodeEnd::FIVE, 4096, ReadsToUse::USE_SECOND){}
    };

    struct SplitSeqV2 : Rule{
        SplitSeqV2(): Rule(24, 10, BarcodeEnd::FIVE, 4294967295, ReadsToUse::USE_SECOND){}
        std::size_t const bcLen = 8, bc1Pos = 10, bc2Pos = 48, bc3Pos = 78;
    };

    struct SplitSeqV1 : Rule{
        SplitSeqV1(): Rule(24, 10, BarcodeEnd::FIVE, 4294967295, ReadsToUse::USE_SECOND){}
        std::size_t const bcLen = 8, bc1Pos = 10, bc2Pos = 48, bc3Pos = 86;
    };

    //dummy class
    struct Custom : Rule{
      Custom() : Rule(0,0,BarcodeEnd::FIVE,0, ReadsToUse::USE_SECOND){}
    };
    struct SciSeq3 : Rule{
      SciSeq3() : Rule(21, 8, BarcodeEnd::FIVE, 1073741824, ReadsToUse::USE_SECOND){}
      std::string anchorSeq = "CAGAGC";
      std::size_t anchorSeqLen = anchorSeq.length();
      std::size_t anchorPos = 0;
      u_int16_t const maxHairpinIndexLen = 10;
      u_int16_t const rtIdxLen = 10; // rev transcription index length
    };

    // for the new type of specification
    struct CustomGeometry {
      // vector of offset, length pairs
      TagGeometry umi_geo;
      TagGeometry bc_geo;
      TagGeometry read_geo;

      void set_umi_geo(TagGeometry& g) { umi_geo = g; umiLength = umi_geo.length(); }
      void set_bc_geo(TagGeometry& g) { bc_geo = g; barcodeLength = bc_geo.length(); }
      void set_read_geo(TagGeometry& g) { read_geo = g; }
      uint32_t barcode_length() const { return barcodeLength; }
      uint32_t umi_length() const { return umiLength; }
      ReadsToUse get_reads_to_use() const { return readsToUse; }

      // These values are set only when `set_umi_geo` and 
      // `set_bc_geo` are called.  See if this design can 
      // be better integrated with `Rule` later.
      uint32_t barcodeLength, umiLength, maxValue;
      BarcodeEnd end;
      ReadsToUse readsToUse = (alevin::defaults::isFivePrimeLibrary) ? ReadsToUse::USE_BOTH : ReadsToUse::USE_SECOND;
    };

  }
}

#endif
