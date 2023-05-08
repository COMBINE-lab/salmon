#ifndef _UTIL__H
#define _UTIL__H

#include "core/range.hpp"
#include "string_view.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>
#include <limits>

#include "CanonicalKmer.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
//#include "jellyfish/mer_dna.hpp"
#include "Kmer.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"
#include "itlib/small_vector.hpp"
#include "parallel_hashmap/phmap.h"
#include "compact_vector/compact_vector.hpp"

#ifdef PUFFERFISH_SALMON_SUPPORT
#include "LibraryFormat.hpp"
#endif

#ifndef __DEFINE_LIKELY_MACRO__
#define __DEFINE_LIKELY_MACRO__
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif
#endif

#define PUFF_DEBUG_VERBOSE 0
#ifdef PUFF_DEBUG_VERBOSE
#  define VERB(x) x
#else
#  define VERB(x)
#endif // PUFF_DEBUG_VERBOSE

namespace pufferfish {

    namespace util {
        using namespace std;

        constexpr const char MPH[] = "mphf.bin";
        constexpr const char CTABLE[] = "ctable.bin";
        constexpr const char CONTIG_OFFSETS[] = "ctg_offsets.bin";
        constexpr const char EQTABLE[] = "eqtable.bin";
        constexpr const char REFLENGTH[] = "reflengths.bin";
        constexpr const char COMPLETEREFLENGTH[] = "complete_ref_lens.bin";
        constexpr const char REFACCUMLENGTH[] = "refAccumLengths.bin";
        constexpr const char REFNAME[] = "reflengths.bin";
        constexpr const char RANK[] = "rank.bin";
        constexpr const char SEQ[] = "seq.bin";
        constexpr const char POS[] = "pos.bin";
        constexpr const char REFSEQ[] = "refseq.bin";
        constexpr const char EDGE[] = "edge.bin";
        constexpr const char PRESENCE[] = "presence.bin";
        constexpr const char CANONICAL[] = "canonical.bin";
        constexpr const char SAMPLEPOS[] = "sample_pos.bin";
        constexpr const char EXTENSION[] = "extension.bin";
        constexpr const char EXTENSIONSIZE[] = "extensionSize.bin";
        constexpr const char DIRECTION[] = "direction.bin";
		    constexpr const char INFO[] = "info.json";

        static constexpr int8_t rc_table[128] = {
                78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 15
                78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 31
                78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 787
                78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 63
                78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78, // 79
                78, 78, 78, 78, 65, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, // 95
                78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78, // 101
                78, 78, 78, 78, 65, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78  // 127
        };

      struct IndexLoadingOpts {
        // These are the optional components for the index
        // this says whether or not we will even try to load
        // these things.  Obviously, if they were not built 
        // during index construction, we cannot load them.
        bool try_loading_eqclasses{false};
        bool try_loading_edges{false};
        bool try_loading_ref_seqs{true};
      };

        enum ReadEnd : uint8_t {
            LEFT, RIGHT
        };

      enum class HitFilterPolicy : uint8_t {
            FILTER_AFTER_CHAINING = 0, FILTER_BEFORE_CHAINING,
            FILTER_BEFORE_AND_AFTER_CHAINING, DO_NOT_FILTER
      };

      // encapsulates policy choices about what types of mappings
      // should be allowed (e.g. orphans, dovetails, etc.)
      struct MappingConstraintPolicy {
        bool noOrphans;
        bool noDiscordant;
        bool noDovetail;
        // after merging chains for paired-end reads 
        // only chains having this threshold score 
        // *with respect to the best chain on the same target*
        // will be passed to the next stage of mapping.
        double post_merge_chain_sub_thresh{0.9};
        double inv_post_merge_chain_sub_thresh{1.0 / post_merge_chain_sub_thresh};
        double orphan_chain_sub_thresh{1.0};

        double postMergeChainSubThresh() const { return post_merge_chain_sub_thresh; }
        void setPostMergeChainSubThresh(double t) {
          post_merge_chain_sub_thresh = t;
          inv_post_merge_chain_sub_thresh = std::nexttoward(1.0/post_merge_chain_sub_thresh, std::numeric_limits<long double>::infinity());
        }

        double orphanChainSubThresh() const { return orphan_chain_sub_thresh; }
        void setOrphanChainSubThresh(const double t) {
          orphan_chain_sub_thresh = t;
        }

      };

      struct pair_hash
      {
        template <class T1, class T2>
        std::size_t operator() (const std::pair<T1, T2> &pair) const
        {
          return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
        }
      };


      enum class BestHitReferenceType : uint8_t { NON_FILTERED, FILTERED, BOTH, UNKNOWN };

      template <typename K, typename V, typename H>
      class CachedVectorMap {
      private:
        phmap::flat_hash_map<K, uint32_t, H> index_map_;
        std::vector<V> cache_;
        uint32_t next_avail_{0};
      public:
        CachedVectorMap(){}

        V& operator[](const K& k) {
          auto it = index_map_.find(k);
          if (it == index_map_.end()) {
            auto idx = next_avail_;
            ++next_avail_;
            index_map_[k] = idx;
            if (idx >= cache_.size()) {
              cache_.emplace_back(V());
              return cache_.back();
            } else {
              cache_[idx].clear();
              return cache_[idx];
            }
          } else {
            return cache_[it->second];
          }
        }

        V& cache_index(uint32_t ci) {
          return cache_[ci];
        }

        decltype(index_map_.begin()) key_begin() { return index_map_.begin(); }
        decltype(index_map_.end()) key_end() { return index_map_.end(); }

        size_t size() const { return index_map_.size(); }

        void clear() {
          next_avail_ = 0;
          index_map_.clear();
        }

        class iterator {

            typedef iterator self_type;
            typedef std::pair<K,V*> value_type;
            typedef value_type& reference;
            typedef value_type* pointer;
//            typedef std::input_iterator_tag iterator_category;
//            typedef int64_t difference_type;

        public:
            explicit iterator(CachedVectorMap &vmIn): vm(vmIn) {
              key = vm.index_map_.begin();
              if (key != vm.index_map_.end()) {
                setKV();
              }
            }

            reference operator*() {
                return kv;
            }

            pointer operator->() { return &operator*(); }

            iterator& operator++() {
                if (++key != vm.index_map_.end())
                    setKV();
                return *this;
            }

            iterator operator++(int) {
                auto tmp = *this;
                ++*this;
                return tmp;
            }

            bool operator==(const self_type& itr) {
                if (key == itr.key and key == vm.index_map_.end()) return true;
                if (key != itr.key) return false;
                if (kv.first != itr.kv.first) return false;
                if (kv.second->size() != itr.kv.second->size()) return false;
                for (uint64_t i = 0; i < kv.second->size(); i++) {
                    if ((*kv.second)[i] != (*itr.kv.second)[i]) return false;
                }
                return true;
            }

            bool operator!=(const self_type& itr) {
                return !((*this) == itr);
            }

            void set2End() {key = vm.index_map_.end();}

          private:
            CachedVectorMap &vm;
            value_type kv;
            decltype(index_map_.begin()) key;

            void setKV() {
                kv.first = key->first;
                kv.second = key->second >= vm.cache_.size()?nullptr:&vm.cache_[key->second];
            }

          };

        iterator begin() {iterator it_(*this); return it_;}
        iterator end() {iterator it_(*this); it_.set2End(); return it_;}
      };



// Adapted from
// https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
        inline void reverseRead(std::string &seq,
                                std::string &qual,
                                std::string &readWork,
                                std::string &qualWork) {

            // NOTE: Should we explicitly check here that qual and seq
            // have the same length? This should have been a parsing 
            // error upstream, but consider a defensive check here too.
            readWork.resize(seq.length(), 'A');
            qualWork.resize(qual.length(), 'I');
            int32_t end = seq.length() - 1, start = 0;
            //readWork[end] = '\0';
            //qualWork[end] = '\0';
            while (LIKELY(start < end)) {
                readWork[start] = (char) rc_table[(int8_t) seq[end]];
                readWork[end] = (char) rc_table[(int8_t) seq[start]];
                qualWork[start] = qual[end];
                qualWork[end] = qual[start];
                ++start;
                --end;
            }
            // If odd # of bases, we still have to complement the middle
            if (start == end) {
                readWork[start] = (char) rc_table[(int8_t) seq[start]];
                // but don't need to mess with quality
                // qualWork[start] = qual[start];
            }
            //std::swap(seq, readWork);
            //std::swap(qual, qualWork);
        }

        // Adapted from
        // https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
        // Don't modify the qual
        inline void reverseRead(std::string &seq,
                                std::string &readWork) {

            readWork.resize(seq.length(), 'A');
            int32_t end = seq.length() - 1, start = 0;
            //readWork[end] = '\0';
            //qualWork[end] = '\0';
            while (LIKELY(start < end)) {
                readWork[start] = (char) rc_table[(int8_t) seq[end]];
                readWork[end] = (char) rc_table[(int8_t) seq[start]];
                ++start;
                --end;
            }
            // If odd # of bases, we still have to complement the middle
            if (start == end) {
                readWork[start] = (char) rc_table[(int8_t) seq[start]];
                // but don't need to mess with quality
                // qualWork[start] = qual[start];
            }
            //std::swap(seq, readWork);
            //std::swap(qual, qualWork);
        }

        inline std::string reverseComplement(std::string &seq) {
            std::string work;
            reverseRead(seq, work);
            return work;
        }

        struct cmpByPair {
            bool operator()(std::pair<uint64_t, bool> a,
                            std::pair<uint64_t, bool> b) const {
                return (a.first != b.first) ? (a.first < b.first) : (a.second < b.second);
            }
        };

// We need a wraper that will provide a "type" field
        template<typename>
        struct Void {
            typedef void type;
        };

// The "default" type T does not have a key
        template<typename T, typename Sfinae = void>
        struct has_key : std::false_type {
        };

// Use decltype to access the key_type field of the
// provided type T if it exists.  If T has a key_type trait
// then this will succeed, and we can instantiate this struct.
// If T doesn't have a key_type trait, then SFINAE means that that
// struct won't be instantiated and we'll fall back to the default
// above.
        template<typename T>
        struct has_key<T, typename Void<decltype(typename T::key_type{})>::type>
                : std::true_type {
        };

/*
Print a map-like container
 */
        template<typename T, bool b>
        struct ContainerPrinter {
            static string str(T &container) {
                ostringstream s;
                /*
                 Loop over and print out the map using the new c++11
                 style for loop.  The new-style loop will work for any
                 type that has begin() and end() methods that return an iterator,
                 OR if there are _free_ functions begin() and end() which
                 consume the type and return an iterator.
                */
                s << "[";
                for (auto &e : container) {
                    s << e.first << " : " << e.second << ", ";
                }
                s.seekp(-2, ios_base::end);
                s << "]";
                return s.str();
            }
        };

/*
Print a list-like container
 */
        template<typename T>
        struct ContainerPrinter<T, false> {
            static string str(T &container) {
                ostringstream s;
                s << "[";
                for (auto &e : container) {
                    s << e << ", ";
                }
                s.seekp(-2, ios_base::end);
                s << "]";
                return s.str();
            }
        };

/*
Compile-time selection between list-like and map-like printing.
 */
        template<typename T>
        string str(T &container) {
            return ContainerPrinter<T, has_key<T>::value>::str(container);
        }

      struct CIGARGenerator {
        // TODO: @fataltes --- think about just replacing this
        // with CIGAR Op class or some such.
        std::vector<uint32_t> cigar_counts;
        std::string cigar_types;
        int32_t begin_softclip_len{0}, end_softclip_len{0};
        bool beginOverhang{false}, endOverhang{false};
        bool allowOverhangSoftClip{false};

        void clear() {
          cigar_counts.clear(); cigar_types.clear();
          begin_softclip_len = end_softclip_len = 0;
          beginOverhang = endOverhang = false;
        }

        void add_item(uint32_t count, char type) {
          cigar_counts.push_back(count);
          cigar_types.push_back(type);
        }

        void get_approx_cigar(int32_t readLen, std::string& cigar) {
          if (begin_softclip_len > 0) {
            cigar += std::to_string(begin_softclip_len);
            cigar += beginOverhang and allowOverhangSoftClip ? "I" : "S";
          }
          cigar += std::to_string(readLen - (begin_softclip_len + end_softclip_len));
          cigar += "M";
          if (end_softclip_len > 0) {
            cigar += std::to_string(end_softclip_len);
            cigar += endOverhang and allowOverhangSoftClip ? "I" : "S";
          }
        }

        std::string get_cigar(uint32_t readLen, bool &cigar_fixed) {
          cigar_fixed = false;
          std::string cigar = "";
          if (cigar_counts.size() != cigar_types.size() or cigar_counts.size() == 0) {
            return "!";
          }
          if (cigar_counts.size() == 0) {
            return cigar;
          }

          uint32_t cigar_length = 0;
          uint32_t count = cigar_counts[0];

          if (cigar_counts.size() == 1) {
            if (count != readLen) {
              count = readLen;
              cigar_fixed = true;
            }
            cigar += std::to_string(count);
            cigar += cigar_types[0];
            return cigar;
          }

          char type = cigar_types[0];
          if (type == 'I' or type == 'M')
            cigar_length += count;
          for (size_t i = 1; i < cigar_counts.size(); i++) {
            if (cigar_types[i] == 'I' or cigar_types[i] == 'M')
              cigar_length += cigar_counts[i];
            if (type == cigar_types[i]) {
              count += cigar_counts[i];
            } else {
              cigar += std::to_string(count);
              cigar += type;
              count = cigar_counts[i];
              type = cigar_types[i];
            }
            if (i == cigar_counts.size() - 1) {
              cigar += std::to_string(count);
              cigar += type;
              if (cigar_length < readLen) {
                cigar_fixed = true;
                count = readLen - cigar_length;
                cigar += std::to_string(count);
                cigar += 'I';
              } else if (cigar_length > readLen) {
                cigar_fixed = true;
                count = cigar_length - readLen;
                cigar += std::to_string(count);
                cigar += 'I';
              }
            }
          }
          return cigar;
        }
      };

        //Mapped object contains all the information
        //about mapping the struct is a bit changed from
        //quasi mapping
        enum class MateStatus : uint8_t {
            SINGLE_END = 0,
            PAIRED_END_LEFT = 1,
            PAIRED_END_RIGHT = 2,
            PAIRED_END_PAIRED = 3
        };

        // Reporting merge (join) status in paired end reads
        enum class MergeResult : uint8_t {
            HAD_NONE,
            HAD_EMPTY_INTERSECTION,
            HAD_CONCORDANT,
            HAD_DISCORDANT,
            HAD_ONLY_LEFT,
            HAD_ONLY_RIGHT,
        };

        //required for edge extension

        enum class Direction : bool {
            FORWARD = 0, BACKWORD = 1
        };

        enum class Type : bool {
            APPEND = 0, PREPEND = 1
        };


        struct extension {
            char c;
            Direction dir;
        };


        struct ContigCecheBlock {
            uint32_t cpos;
            uint32_t blockLen;
            std::string cseq;


            ContigCecheBlock(uint32_t cposIn, uint32_t lenIn, std::string cseqIn) :
                    cpos(cposIn), blockLen(lenIn), cseq(cseqIn) {}

        };

        struct UniMemInfo {
            uint32_t cid;
            bool cIsFw;
            uint32_t rpos;
            uint32_t memlen;
            uint32_t cpos;
            uint64_t cGlobalPos;
            uint32_t clen;
            ReadEnd readEnd;


            UniMemInfo(uint32_t cidIn, bool cIsFwIn, uint32_t rposIn, uint32_t memlenIn, uint32_t cposIn,
                       uint64_t cGlobalPosIn, uint32_t clenIn, ReadEnd rendIn = ReadEnd::LEFT) :
                    cid(cidIn), cIsFw(cIsFwIn), rpos(rposIn), memlen(memlenIn), cpos(cposIn), cGlobalPos(cGlobalPosIn),
                    clen(clenIn), readEnd(rendIn) {}

            UniMemInfo() : cid(0), cIsFw(false), rpos(3), memlen(1), cpos(0), cGlobalPos(0), clen(0), readEnd(LEFT) {}
        };

        struct MemInfo {
            std::vector<UniMemInfo>::iterator memInfo;
            size_t tpos;
            bool isFw;
            uint32_t extendedlen;
            uint32_t rpos;

            //MemInfo(std::vector<UniMemInfo>::iterator uniMemInfoIn, size_t tposIn, bool isFwIn = true) :
            //        memInfo(uniMemInfoIn), tpos(tposIn), isFw(isFwIn) {}
            MemInfo(std::vector<UniMemInfo>::iterator uniMemInfoIn, size_t tposIn, bool isFwIn = true) :
                    memInfo(uniMemInfoIn), tpos(tposIn), isFw(isFwIn) {
                extendedlen = uniMemInfoIn->memlen;
                rpos = uniMemInfoIn->rpos;
            }

            MemInfo(std::vector<UniMemInfo>::iterator uniMemInfoIn, size_t tposIn, uint32_t extendedlenIn,
                    uint32_t rposIn, bool isFwIn = true) :
                    memInfo(uniMemInfoIn), tpos(tposIn), isFw(isFwIn), extendedlen(extendedlenIn), rpos(rposIn) {}
            bool operator==(const MemInfo& mi) {
                return memInfo == mi.memInfo and tpos == mi.tpos and isFw == mi.isFw and extendedlen == mi.extendedlen and rpos == mi.rpos;
            }
            bool operator!=(const MemInfo& mi) {
                return !((*this) == mi);
            }
        };

        struct MemCluster {
            // second element is the transcript position
            std::vector<MemInfo> mems;
            bool isFw;
            bool isVisited = false;
            double coverage{0};
//            std::vector<std::pair<std::string, std::string>> alignableStrings; //NOTE we don't need it [cigar on the fly]
            int score;
            std::string cigar;
            bool perfectChain = false;
            uint32_t readLen{0};
            uint32_t openGapLen{0};
            uint16_t softClipStart{0};
            uint64_t queryChainHash{0};
            
            //bool isValid = true;
            MemCluster(bool isFwIn, uint32_t readLenIn) : isFw(isFwIn), readLen(readLenIn) {}
            /*MemCluster(bool isFwIn, MemInfo memIn): isFw(isFwIn) {
              mems.push_back(memIn);
              }*/
            /*MemCluster() {
              mems.push_back(std::make_pair(memInfoPtr, tposIn));
              }*/
            MemCluster(const MemCluster &other) = default;

            MemCluster &operator=(const MemCluster &other) = default;

            bool operator==(const MemCluster& mc) {
                if (!(isFw == mc.isFw and score == mc.score and coverage == mc.coverage
                and cigar == mc.cigar and perfectChain == mc.perfectChain
                and readLen == mc.readLen and openGapLen == mc.openGapLen 
                and softClipStart == mc.softClipStart )) return false;
                if (mems.size() != mc.mems.size()) return false;
                for (uint64_t i = 0; i < mems.size(); i++) {
                    if (mems[i] != mc.mems[i]) return false;
                }
                return true;
            }

            bool operator!=(const MemCluster& mc) {
                return !((*this) == mc);
            }

            MemCluster() {}

            // Add the new mem to the list and update the coverage, designed for clustered Mems
            void addMem(std::vector<UniMemInfo>::iterator uniMemInfo, size_t tpos, uint32_t extendedlen, uint32_t rpos,
                        bool isFw) {
                if (mems.empty())
                    coverage = extendedlen;
                else if (tpos > mems.back().tpos + mems.back().extendedlen) {
                    coverage += extendedlen;
                } else { // they overlap
                    coverage += (uint32_t) std::max(
                            (int) (tpos + extendedlen) - (int) (mems.back().tpos + mems.back().extendedlen), 0);
                }
                mems.emplace_back(uniMemInfo, tpos, extendedlen, rpos, isFw);
            }

          inline int64_t getReadLastHitPos() const { return mems.empty() ? 0 : static_cast<int64_t>(mems.back().rpos); }

          inline int64_t getTrLastHitPos() const {
            return mems.empty() ? 0 : static_cast<int64_t>(mems.back().tpos);
          }

          inline int64_t getTrLastMemLen() const {
            //return mems.empty()?0:mems.back().memInfo->memlen;
            return mems.empty() ? 0 : static_cast<int64_t>(mems.back().extendedlen);
          }

          inline int64_t getTrFirstHitPos() const {
            return mems.empty() ? 0 : static_cast<int64_t>(mems[0].tpos) - openGapLen;
            //return mems.empty() ? 0 : (mems[0].isFw ? mems[0].tpos-mems[0].rpos : mems[0].tpos - (readLen-mems[0].rpos-mems[0].extendedlen));
          }

          // returns the approximate read start position (approximate because
          // alignment / selective alignment hasn't been computed yet, and we 
          // don't know where the read necessarily starts).  If this value is 
          // less than 0, then the read overhangs the start of the transcript.
          inline int64_t approxReadStartPos() const {
            if (mems.empty()) { return 0; }
            auto& m = mems.front();
            return isFw ? (static_cast<int64_t>(m.tpos) - static_cast<int64_t>(m.rpos)) :
              (static_cast<int64_t>(m.tpos) - ((static_cast<int64_t>(readLen) - static_cast<int64_t>(m.rpos + m.extendedlen))));
          }

            inline int64_t firstRefPos() const { return getTrFirstHitPos(); }

            inline int64_t lastRefPos() const { return getTrLastHitPos(); }

            inline size_t lastMemLen() const { return getTrLastMemLen(); }

            void calcCoverage() {
                if (mems.size() == 0) {
                    std::cerr << "Shouldn't happen! Cluster is empty.\n";
                    return;
                }
                // we keep prev to take care of overlaps while calculating the coverage
                auto lstart = mems.begin();
                size_t offset = 0;
                auto prev = lstart;
                //coverage = mems[0].memInfo->memlen;
                coverage = mems[0].extendedlen;
                for (auto &&mem : mems) {
                    ++offset;
                    //coverage += std::max((int)(mem.tpos+mem.memInfo->memlen) - (int)(prev->tpos+prev->memInfo->memlen), 0);
                    coverage += std::max((int) (mem.tpos + mem.extendedlen) - (int) (prev->tpos + prev->extendedlen),
                                         0);
                    prev = lstart + offset;
                }
            }
        };


        struct JointMems {
            uint32_t tid;
            std::vector<pufferfish::util::MemCluster>::iterator leftClust;
            std::vector<pufferfish::util::MemCluster>::iterator rightClust;
            size_t fragmentLen;
            size_t rmemMaxLen{0}, lmemMaxLen{0};
            int32_t alignmentScore{0};
            int32_t mateAlignmentScore{0};
            MateStatus mateStatus;
            bool recovered{false};

            bool isLeftAvailable() {
                return mateStatus == MateStatus::PAIRED_END_PAIRED ||
                       mateStatus == MateStatus::PAIRED_END_LEFT;
            }

            bool isRightAvailable() {
                return mateStatus == MateStatus::PAIRED_END_PAIRED ||
                       mateStatus == MateStatus::PAIRED_END_RIGHT;
            }

            bool isOrphan() { return !isLeftAvailable() || !isRightAvailable(); }

          // NOTE: needed for vector compatibility, should not be used.
          JointMems() {
            std::cerr << "JointMems default constructor called; should not happen!";
          }
            JointMems(uint32_t tidIn,
                      std::vector<pufferfish::util::MemCluster>::iterator leftClustIn,
                      std::vector<pufferfish::util::MemCluster>::iterator rightClustIn,
                      size_t fragmentLenIn,
                      MateStatus mateStatusIn = MateStatus::PAIRED_END_PAIRED) :
                    tid(tidIn), leftClust(leftClustIn), rightClust(rightClustIn),
                    fragmentLen(fragmentLenIn), mateStatus(mateStatusIn) {
                /*if (leftClust and rightClust) {
                    if (!leftClust->mems.size() and !rightClust->mems.size()) {
                        std::cerr << "ERROR: Both read end mem lists cannot be empty.\n";
                        std::exit(1);
                    }
                }*/
                //if (leftClust->mems.size() == 0) mateStatus = MateStatus::PAIRED_END_RIGHT;
                //if (rightClust->mems.size() == 0) mateStatus = MateStatus::PAIRED_END_LEFT;
            }

            double coverage() {
                return (isLeftAvailable() ? leftClust->coverage : 0) + (isRightAvailable() ? rightClust->coverage : 0);
            }

            double leftCoverage() { return isLeftAvailable() ? leftClust->coverage : 0; }

            double rightCoverage() { return isRightAvailable() ? rightClust->coverage : 0; }

            // FIXME : what if the mapping is not orphan? who takes care of this function not being called from outside?
            std::vector<pufferfish::util::MemCluster>::iterator orphanClust() {
              return isLeftAvailable() ? leftClust : rightClust;
            }

        };

      /*
    struct QuasiAlignment {
  	QuasiAlignment() :
		tid(std::numeric_limits<uint32_t>::max()),
		pos(std::numeric_limits<int32_t>::max()),
		fwd(true),
		fragLen(std::numeric_limits<uint32_t>::max()),
		readLen(std::numeric_limits<uint32_t>::max()),
		isPaired(false)
#ifdef PUFFERFISH_SALMON_SUPPORT
        ,format(LibraryFormat::formatFromID(0))
#endif // PUFFERFISH_SALMON_SUPPORT
        {}

        QuasiAlignment(uint32_t tidIn, int32_t posIn,
                bool fwdIn, uint32_t readLenIn,
                uint32_t fragLenIn = 0,
                bool isPairedIn = false) :
            tid(tidIn), pos(posIn), fwd(fwdIn),
            fragLen(fragLenIn), readLen(readLenIn), 
            isPaired(isPairedIn)
#ifdef PUFFERFISH_SALMON_SUPPORT
        ,format(LibraryFormat::formatFromID(0))
#endif // PUFFERFISH_SALMON_SUPPORT
        {}
        QuasiAlignment(QuasiAlignment&& other) = default;
        QuasiAlignment& operator=(QuasiAlignment&) = default;
        QuasiAlignment& operator=(QuasiAlignment&& o) = default;
        QuasiAlignment(const QuasiAlignment& o) = default;
        QuasiAlignment(QuasiAlignment& o) = default;

      inline void setChainScore(double chainScoreIn) {
        chainScore_ = chainScoreIn;
      }

      inline double chainScore() const {
        return chainScore_;
      }

      inline uint32_t transcriptID() const { return tid; }
      inline double score() const { return score_; }
      inline void score(double scoreIn) { score_ = scoreIn; }
      inline int32_t alnScore() const { return alnScore_; }
      inline void alnScore(int32_t alnScoreIn) { alnScore_ = alnScoreIn; }
      inline uint32_t fragLength() const { return fragLen; }
      inline int32_t hitPos() { return std::min(pos, matePos); }

// Some convenience functions to allow salmon interop
#ifdef RAPMAP_SALMON_SUPPORT
      inline uint32_t fragLengthPedantic(uint32_t txpLen) const {
        if (mateStatus != rapmap::utils::MateStatus::PAIRED_END_PAIRED
            or fwd == mateIsFwd) {
          return 0;
        }
        int32_t p1 = fwd ? pos : matePos;
        int32_t sTxpLen = static_cast<int32_t>(txpLen);
        p1 = (p1 < 0) ? 0 : p1;
        p1 = (p1 > sTxpLen) ? sTxpLen : p1;
        int32_t p2 = fwd ? matePos + mateLen : pos + readLen;
        p2 = (p2 < 0) ? 0 : p2;
        p2 = (p2 > sTxpLen) ? sTxpLen : p2;

        return (p1 > p2) ? p1 - p2 : p2 - p1;
      }

      double logProb{HUGE_VAL};
      double logBias{HUGE_VAL};
      inline LibraryFormat libFormat() { return format; }
      LibraryFormat format;
#endif // RAPMAP_SALMON_SUPPORT
       bool hasMultiPos{false};
       chobo::small_vector<int32_t> allPositions;
       chobo::small_vector<int32_t> oppositeStrandPositions;

        // Only 1 since the mate must have the same tid
        // we won't call *chimeric* alignments here.
        uint32_t tid;
        // Left-most position of the hit
        int32_t pos;
        // left-most position of the mate
        int32_t matePos;
        // Is the read from the forward strand
        bool fwd;
        // Is the mate from the forward strand
        bool mateIsFwd;
        // The fragment length (template length)
        // This is 0 for single-end or orphaned reads.
        uint32_t fragLen;
        // The read's length
        uint32_t readLen;
        // The mate's length
        uint32_t mateLen;
        // Is this a paired *alignment* or not
        bool isPaired;
        MateStatus mateStatus;
        // numeric score associated with this mapping
        double score_{1.0};
        // actual ``alignment'' score associated with this mapping.
        int32_t alnScore_{0};
        // If one or both of the reads is a complete match (no mismatch, indels), say what kind.
        FragmentChainStatus chainStatus;
        double chainScore_{std::numeric_limits<double>::lowest()};
      //int32_t queryOffset{-1};
      //MateStatus completeMatchType{MateStatus::NOTHING};
    };
    */

      enum class PuffAlignmentMode : uint8_t { SCORE_ONLY, APPROXIMATE_CIGAR,  EXACT_CIGAR};

      struct AlignmentConfig {
        int32_t refExtendLength{20};
        bool fullAlignment{false};
        int16_t matchScore;
        int16_t gapExtendPenalty;
        int16_t gapOpenPenalty;
        int16_t mismatchPenalty;
        double minScoreFraction{0.0};
        bool mimicBT2{false};
        bool mimicBT2Strict{false};
        bool allowOverhangSoftclip{false};
        bool allowSoftclip{false};
        bool useAlignmentCache{true};
        bool noDovetail{false};
        uint32_t maxFragmentLength{1000};
        PuffAlignmentMode alignmentMode{PuffAlignmentMode::SCORE_ONLY};
        bool bestStrata{false};
        bool decoyPresent{false};
      };

        struct QuasiAlignment {
            QuasiAlignment() :
                    tid(std::numeric_limits<uint32_t>::max()),
                    pos(std::numeric_limits<int32_t>::max()),
                    matePos(std::numeric_limits<int32_t>::max()),
                    fwd(true),
                    mateIsFwd(true),
                    fragLen(std::numeric_limits<uint32_t>::max()),
                    readLen(std::numeric_limits<uint32_t>::max()),
                    mateLen(std::numeric_limits<uint32_t>::max()),
                    isPaired(false),
#ifdef PUFFERFISH_SALMON_SUPPORT
                    format(LibraryFormat::formatFromID(0)),
#endif // PUFFERFISH_SALMON_SUPPORT
                    score(std::numeric_limits<int32_t>::min()),
                    mateScore(std::numeric_limits<int32_t>::min())
          {}

            QuasiAlignment(uint32_t tidIn, int32_t posIn,
                           bool fwdIn, uint32_t readLenIn, std::string cigarIn, //NOTE can we make it uint32?
                           uint32_t fragLenIn = 0,
                           bool isPairedIn = false) :
                    tid(tidIn), pos(posIn), fwd(fwdIn),
                    mateIsFwd(true),
                    fragLen(fragLenIn), readLen(readLenIn),
                    mateLen(std::numeric_limits<uint32_t>::max()),
                    isPaired(isPairedIn), cigar(cigarIn),
#ifdef PUFFERFISH_SALMON_SUPPORT
                    format(LibraryFormat::formatFromID(0)),
#endif // PUFFERFISH_SALMON_SUPPORT
                    score(std::numeric_limits<int32_t>::min()),
              mateScore(std::numeric_limits<int32_t>::min())
          {}

            QuasiAlignment(QuasiAlignment &&other) = default;

            QuasiAlignment &operator=(QuasiAlignment &) = default;

            QuasiAlignment &operator=(QuasiAlignment &&o) = default;

            QuasiAlignment(const QuasiAlignment &o) = default;

            QuasiAlignment(QuasiAlignment &o) = default;

          // Some convenience functions to allow salmon interop
          inline uint32_t transcriptID() const { return tid; }
          //inline int32_t alnScore() const { return alnScore_; }
          //inline void alnScore(int32_t alnScoreIn) { alnScore_ = alnScoreIn; }
          inline uint32_t fragLength() const { return fragLen; }
          inline int32_t hitPos() { return std::min(pos, matePos); }

          // Some convenience functions to allow salmon interop
#ifdef PUFFERFISH_SALMON_SUPPORT
          inline uint32_t fragLengthPedantic(uint32_t txpLen) const {
            if (mateStatus != pufferfish::util::MateStatus::PAIRED_END_PAIRED
                or fwd == mateIsFwd) {
              return 0;
            }
            int32_t p1 = fwd ? pos : matePos;
            int32_t sTxpLen = static_cast<int32_t>(txpLen);
            p1 = (p1 < 0) ? 0 : p1;
            p1 = (p1 > sTxpLen) ? sTxpLen : p1;
            int32_t p2 = fwd ? matePos + mateLen : pos + readLen;
            p2 = (p2 < 0) ? 0 : p2;
            p2 = (p2 > sTxpLen) ? sTxpLen : p2;

            return (p1 > p2) ? p1 - p2 : p2 - p1;
          }

          double estAlnProb_{0.0};
          double logProb{HUGE_VAL};
          inline LibraryFormat libFormat() { return format; }
          LibraryFormat format;

          inline double estAlnProb() const { return estAlnProb_; }
          inline void estAlnProb(double scoreIn) { estAlnProb_ = scoreIn; }
#endif // PUFFERFISH_SALMON_SUPPORT

/*
#ifdef RAPMAP_SALMON_SUPPORT
        inline uint32_t transcriptID() const { return tid; }
        inline double score() { return 1.0; }
        inline uint32_t fragLength() const { return fragLen; }

        inline uint32_t fragLengthPedantic(uint32_t txpLen) const {
            if (mateStatus != rapmap::utils::MateStatus::PAIRED_END_PAIRED
                or fwd == mateIsFwd) {
                return 0;
            }
            int32_t p1 = fwd ? pos : matePos;
            p1 = (p1 < 0) ? 0 : p1;
            p1 = (p1 > txpLen) ? txpLen : p1;
            int32_t p2 = fwd ? matePos + mateLen : pos + readLen;
            p2 = (p2 < 0) ? 0 : p2;
            p2 = (p2 > txpLen) ? txpLen : p2;

            return (p1 > p2) ? p1 - p2 : p2 - p1;
        }

        inline int32_t hitPos() { return std::min(pos, matePos); }
        double logProb{HUGE_VAL};
        double logBias{HUGE_VAL};
        inline LibraryFormat libFormat() { return format; }
        LibraryFormat format;
        double estAlnProb;
        double logProb;


#endif // RAPMAP_SALMON_SUPPORT
*/
            // Only 1 since the mate must have the same tid
            // we won't call *chimeric* alignments here.
            uint32_t tid;
            // Left-most position of the hit
            int32_t pos;
            // left-most position of the mate
            int32_t matePos;
            // Is the read from the forward strand
            bool fwd;
            // Is the mate from the forward strand
            bool mateIsFwd;
            // The fragment length (template length)
            // This is 0 for single-end or orphaned reads.
            uint32_t fragLen;
            // The read's length
            uint32_t readLen;
            // The mate's length
            uint32_t mateLen;
            // Is this a paired *alignment* or not
            bool isPaired;


            std::string cigar;
            std::string mateCigar;

            int32_t score;
            int32_t mateScore;

            MateStatus mateStatus;
            bool active = true;
            uint32_t numHits = 0;
        };

// from https://github.com/cppformat/cppformat/issues/105
        class FixedBuffer : public fmt::Buffer<char> {
        public:
            FixedBuffer(char *array, std::size_t size)
                    : fmt::Buffer<char>(array, size) {}

        protected:
            void grow(std::size_t size) {
                (void) size;
                throw std::runtime_error("buffer overflow");
            }
        };

        class FixedWriter : public fmt::Writer {
        private:
            FixedBuffer buffer_;
        public:
            FixedWriter(char *array, std::size_t size)
                    : fmt::Writer(buffer_), buffer_(array, size) {}
        };


// For the time being, assume < 4B contigs
// and that each contig is < 4B bases
        struct Position {
            // std::string transcript_id;
            uint32_t transcript_id_;
            uint32_t pos_;

            // bool orien;
            Position() {
                transcript_id_ = std::numeric_limits<decltype(transcript_id_)>::max();
                pos_ = std::numeric_limits<decltype(pos_)>::max();
            }

            Position(uint32_t tid, uint32_t tpos, bool torien) {
                transcript_id_ = tid;
                pos_ = tpos;
                setOrientation(torien);
                // orien = torien;
            }

            //The most significant bit carry
            //the orientation information

            void setOrientation(bool orientation) {
                if (orientation) {
                    pos_ |= 1 << 31;
                } else {
                    pos_ &= 0x7FFFFFFF;
                }
            }

            inline uint32_t transcript_id() { return transcript_id_; }

            inline uint32_t pos() { return (pos_ & 0x7FFFFFFF); }

            inline bool orientation() { return (pos_ & 0x80000000); }

            template<class Archive>
            void serialize(Archive &ar) {
                ar(transcript_id_, pos_);
            }

            void update(uint32_t tid, uint32_t tpos, bool torien) {
                transcript_id_ = tid;
                pos_ = tpos;
                setOrientation(torien);
            }

        private:
            // uint32_t orientMask_
        };

//struct HitPos
        struct HitQueryPos {
            HitQueryPos(uint32_t queryPosIn, uint32_t posIn, bool queryFwdIn) :
                    queryPos(queryPosIn), pos(posIn), queryFwd(queryFwdIn) {}

            uint32_t queryPos, pos;
            bool queryFwd;

        };

        struct QueryCache {
            uint64_t prevRank{std::numeric_limits<uint64_t>::max()};
            uint64_t contigStart{std::numeric_limits<uint64_t>::max()};
            uint64_t contigEnd{std::numeric_limits<uint64_t>::max()};
        };

        struct ContigPosInfo {
            size_t offset_;
            uint32_t length_;

            inline size_t offset() { return offset_; }

            inline uint32_t length() { return length_; }

            template<class Archive>
            void serialize(Archive &ar) { ar(offset_, length_); }
        };

        struct PackedContigInfo {
            size_t fileOrder;
            size_t offset;
            uint32_t length;

            PackedContigInfo(size_t fileOrder, size_t offset, uint32_t length) : fileOrder(fileOrder),
                                                                                 offset(offset),
                                                                                 length(length) {}
        };

        struct PackedContigInfoVec {
          uint64_t useq_len_{0};
          std::unique_ptr<std::vector<uint64_t>> data_{nullptr};

          PackedContigInfoVec() {}

          PackedContigInfoVec(uint64_t useq_len,
                              size_t ncontig) {
            useq_len_ = useq_len;
            data_.reset(new std::vector<uint64_t>);
            data_->reserve(ncontig);
          }

          void add(uint64_t o) { data_->push_back(o); }

          PackedContigInfo operator[](uint64_t i) const {
            uint32_t len = (i < data_->size() - 1)
                               ? ((*data_)[i + 1] - (*data_)[i])
                               : (useq_len_ - (*data_)[i]);
            return PackedContigInfo{i, (*data_)[i], len};
          }

          void clear() {
            useq_len_ = 0;
            data_->clear();
            data_->shrink_to_fit();
            data_.reset(nullptr);
          }

          size_t size() { return data_->size(); }

          struct PackedContigInfoVecIterator {
            const PackedContigInfoVec* pci_{nullptr};
            uint64_t it{std::numeric_limits<uint64_t>::max()};
	    std::pair<uint64_t, PackedContigInfo> p_{std::numeric_limits<uint64_t>::max(), {size_t(0), size_t(0), uint32_t(0)}};

            PackedContigInfoVecIterator& operator++() {
              ++it;
              return *this;
            }

            bool operator==(const PackedContigInfoVecIterator& o) const {
              return it == o.it;
            }

            bool operator!=(const PackedContigInfoVecIterator& other) const {
              return !(*this == other);
            }

            std::pair<uint64_t, PackedContigInfo>& operator*() {
              p_.second = (*pci_)[it];
              p_.first = it;
              return p_;
            }
          };

          const PackedContigInfoVecIterator find(uint64_t idx) const {
            PackedContigInfoVecIterator i;
            i.pci_ = this;
            i.it = (idx >= data_->size()) ? data_->size() : idx;
            return i;
          }

          const PackedContigInfoVecIterator begin() const {
            PackedContigInfoVecIterator i;
            i.pci_ = this;
            i.it = 0;
            return i;
          }

          const PackedContigInfoVecIterator end() const {
            PackedContigInfoVecIterator i;
            i.pci_ = this;
            i.it = data_->size();
            return i;
          }
        };

        struct RefPos {
            uint32_t pos;
            bool isFW;
        };

        struct HitCounters {
            std::atomic<uint64_t> numMapped{0};
            std::atomic<uint64_t> numMappedAtLeastAKmer{0};
            std::atomic<uint64_t> numOfOrphans{0};
            std::atomic<uint64_t> peHits{0};
            std::atomic<uint64_t> seHits{0};
            std::atomic<uint64_t> trueHits{0};
            std::atomic<uint64_t> totHits{0};
            std::atomic<uint64_t> numReads{0};
            std::atomic<uint64_t> tooManyHits{0};
            std::atomic<uint64_t> lastPrint{0};
            std::atomic<uint64_t> totAlignment{0};
            std::atomic<uint64_t> correctAlignment{0};
            std::atomic<uint64_t> maxMultimapping{0};
            std::atomic<uint64_t> numDovetails{0};
            std::atomic<uint64_t> tooShortReads{0};

            std::atomic<uint64_t> skippedAlignments_notAlignable{0};
            std::atomic<uint64_t> skippedAlignments_byCache{0};
            std::atomic<uint64_t> skippedAlignments_byCov{0};
            std::atomic<uint64_t> totalAlignmentAttempts{0};
            std::atomic<uint64_t> cigar_fixed_count{0};
        };

        struct ContigBlock {
            ContigBlock() : 
              contigIdx_(std::numeric_limits<uint64_t>::max()),
              globalPos_(std::numeric_limits<uint64_t>::max()),
              contigLen_(std::numeric_limits<uint32_t>::max()),
              isDummy_(false) {}

            ContigBlock(uint64_t idIn, uint64_t cposIn, uint32_t len, std::string seqIn, bool isDummyIn = false) :
                    contigIdx_(idIn), globalPos_(cposIn), contigLen_(len), seq(seqIn), isDummy_(isDummyIn) {}

            uint64_t contigIdx_;
            uint64_t globalPos_;
            uint32_t contigLen_;

            std::string seq;

            bool isDummy_;

            bool isDummy() { return isDummy_; }

            std::string substrSeq(size_t s, size_t len) {
                //std::cerr << contigLen_ << "\t" << s << "\t" << len << "\n" ;
                if (s + len <= contigLen_) {
                    return seq.substr(s, len);
                } else
                    return "";
            }


        };


// Structure to hold a list of "projected" (i.e. reference) hits
// for a k-mer


        struct ProjectedHits {
            uint32_t contigIdx_;
            // The relative position of the k-mer inducing this hit on the
            // contig
            uint64_t globalPos_;

            uint32_t contigPos_;
            // How the k-mer inducing this hit maps to the contig
            // true for fw, false for rc
            bool contigOrientation_;
            uint32_t contigLen_;
            uint32_t k_;
            core::range<std::vector<pufferfish::util::Position>::iterator> refRange;

            inline bool empty() { return refRange.empty(); }

            inline uint32_t contigID() const { return contigIdx_; }

            //inline uint64_t getGlobalPos() const { return globalPos_; }
            inline RefPos decodeHit(pufferfish::util::Position &p) {
                // true if the contig is fowrard on the reference
                bool contigFW = p.orientation();
                // we are forward with respect to the reference if :
                // (1) contigFW and contigOrientation_
                // (2) !contigFW and !contigOrientation_
                // we are reverse complement with respect to the reference if :
                // (3) configFW and !contigOrientation_
                // (4) !configFW and contigOrientation_

                // if we're in the forward orientation, then our position is
                // just the contig offset plus or relative position
                uint32_t rpos;//{0};
                bool rfw;//{false};
                if (contigFW and contigOrientation_) {
                    // kmer   :          AGC
                    // contig :      ACTTAGC
                    // ref    :  GCA[ACTTAGC]CA
                    rpos = p.pos() + contigPos_;
                    rfw = true;
                } else if (contigFW and !contigOrientation_) {
                    // kmer   :          GCT
                    // contig :      ACTTAGC
                    // ref    :  GCA[ACTTAGC]CA
                    rpos = p.pos() + contigPos_;
                    rfw = false;
                } else if (!contigFW and contigOrientation_) {
                    // kmer   :          AGT
                    // contig :      GCTAAGT
                    // ref    :  GCA[ACTTAGC]CA
                    rpos = p.pos() + contigLen_ - (contigPos_ + k_);
                    rfw = false;
                } else {// if (!contigFW and !contigOrientation_) {
                    // kmer   :          ACT
                    // contig :      GCTAAGT
                    // ref    :  GCA[ACTTAGC]CA
                    rpos = p.pos() + contigLen_ - (contigPos_ + k_);
                    rfw = true;
                }

                return {rpos, rfw};
            }
        };

        struct AlignmentResult {
            AlignmentResult(bool isFwIn, int32_t scoreIn, std::string cigarIn, uint32_t openGapLenIn, 
                            uint16_t softclip_start_in ) :
                    isFw(isFwIn), score(scoreIn), cigar(cigarIn), openGapLen(openGapLenIn), softclip_start(softclip_start_in) {}

            AlignmentResult() : isFw(true), score(0), cigar(""), openGapLen(0), softclip_start(0) {}
            bool isFw;
            int32_t score;
            std::string cigar;
            uint32_t openGapLen;
            uint16_t softclip_start{0};
        };

      void joinReadsAndFilterSingle( pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& leftMemClusters,
                                     //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> &leftMemClusters,
                                     std::vector<pufferfish::util::JointMems> &jointMemsList,
                                     uint32_t perfectCoverage,
                                     double coverageRatio);

      pufferfish::util::MergeResult joinReadsAndFilter(
                                                       pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& leftMemClusters,
                                                       pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& rightMemClusters,
                                                       std::vector<pufferfish::util::JointMems> &jointMemsList,
                                                       uint32_t maxFragmentLength,
                                                       uint32_t perfectCoverage,
                                                       double coverageRatio,
                                                       uint64_t firstDecoyIndex,
                                                       const pufferfish::util::MappingConstraintPolicy& mpol,
                                                       pufferfish::util::HitCounters& hctr);

      char * getRefSeqOwned(compact::vector<uint64_t, 2> &refseq, uint64_t refAccPos, uint32_t len);

        char complement(char &c);

        std::string revcomp(std::string s);

        bool isRevcomp(std::string s);

        std::vector<std::pair<uint64_t, bool>> explode(const stx::string_view str,
                                                       const char &ch);

        bool is_number(const std::string &s);

        std::vector<std::string> tokenize(const std::string &s, char delim);

// Avoiding un-necessary stream creation + replacing strings with string view
// is a bit > than a 2x win!
// implementation from : https://marcoarena.wordpress.com/tag/string_view/
        std::vector<stx::string_view> split(stx::string_view str, char delims);

        std::vector<extension> getExts(uint8_t e);
    }
}

#endif // _UTIL__H
