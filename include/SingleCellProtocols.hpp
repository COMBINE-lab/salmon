#ifndef __SINGLE_CELL_PROTOCOLS_HPP__
#define __SINGLE_CELL_PROTOCOLS_HPP__

#include <string>
#include <ostream>

#include "AlevinOpts.hpp"
#include "AlevinTypes.hpp"
#include "pufferfish/chobo/static_vector.hpp"

namespace alevin{
  namespace protocols {

    static constexpr size_t num_tag_pieces{16};
    struct TagGeometry {
      // uint32_t read_num{0};
      // tuples are read_num, start_pos, length
      chobo::static_vector<std::pair<uint32_t, size_t>, num_tag_pieces> substr_locs1{};
      chobo::static_vector<std::pair<uint32_t, size_t>, num_tag_pieces> substr_locs2{};
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
           uint32_t maxValue_):
        barcodeLength(barcodeLength_),
        umiLength(umiLength_),
        end(end_),
        maxValue(maxValue_){
        alevin::types::AlevinUMIKmer::k(umiLength);
      }
      // NOTE: these do nothing but satisfy 
      // template requirements right now
      void set_umi_geo(TagGeometry& g) { (void)g; };
      void set_bc_geo(TagGeometry& g) { (void)g; };
      void set_read_geo(TagGeometry& g) { (void)g; };
      uint32_t barcode_length() const { return barcodeLength; }
      uint32_t umi_length() const { return umiLength; }

      uint32_t barcodeLength, umiLength, maxValue;
      BarcodeEnd end;
    };


    struct DropSeq : Rule{
      //Drop-Seq starts from 5 end with 12 length
      //barcode and 8 length umi & iupac can be
      //changed
      DropSeq(): Rule(12, 8, BarcodeEnd::FIVE, 16777216){}
    };

    struct InDrop : Rule{
        //InDrop starts from 5end with variable
        //length barcodes so provide the full
        // length of the barcod eincluding w1.
        // UMI length is 6
      InDrop(): Rule(42, 6, BarcodeEnd::FIVE, 22347776){}

      std::string w1;
      void setW1(std::string& w1_){
        w1 = w1_;
      }
    };

    struct CITESeq : Rule{
      CITESeq(): Rule(16, 10, BarcodeEnd::FIVE, 4294967295){
        featureLength = 15;
        featureStart = 10;
      }

      size_t featureLength, featureStart;
      void setFeatureLength(size_t length) { featureLength = length; }
      void setFeatureStart(size_t startIdx) { featureStart = startIdx; }
    };

    struct ChromiumV3 : Rule{
      ChromiumV3(): Rule(16, 12, BarcodeEnd::FIVE, 4294967295){}
    };

    struct Chromium : Rule{
      Chromium(): Rule(16, 10, BarcodeEnd::FIVE, 4294967295){}
    };

    struct Gemcode : Rule{
      Gemcode(): Rule(14, 10, BarcodeEnd::FIVE, 268435456){}
    };

    struct QuartzSeq2 : Rule{
      QuartzSeq2(): Rule(15, 8, BarcodeEnd::FIVE, 1073741824){}
    };

    struct CELSeq : Rule{
      // WEHI SCORE's CEL-Seq2 starts from 5' end with a 8 bp barcode
      // and a 6 bp UMI.
      CELSeq(): Rule(8, 6, BarcodeEnd::FIVE, 65536){}
    };
    struct CELSeq2 : Rule{
      // WEHI SCORE's CEL-Seq2 starts from 5' end with a 8 bp barcode
      // and a 6 bp UMI.
      CELSeq2(): Rule(6, 6, BarcodeEnd::FIVE, 4096){}
    };

    //dummy class
    struct Custom : Rule{
      Custom() : Rule(0,0,BarcodeEnd::FIVE,0){}
    };
    

    // for the new type of specification
    struct CustomGeometry {
      // vector of offset, length pairs
      TagGeometry umi_geo;
      TagGeometry bc_geo;
      TagGeometry read_geo;

      void set_umi_geo(TagGeometry& g) { umi_geo = g; umiLength = umi_geo.length1 + umi_geo.length2; }
      void set_bc_geo(TagGeometry& g) { bc_geo = g; barcodeLength = bc_geo.length1 + bc_geo.length2; }
      void set_read_geo(TagGeometry& g) { read_geo = g; }
      uint32_t barcode_length() const { return barcodeLength; }
      uint32_t umi_length() const { return umiLength; }

      // These values do nothing in this class except
      // maintain template compat ... fix this design later.
      uint32_t barcodeLength, umiLength, maxValue;
      BarcodeEnd end;
    };

  }
}

#endif
