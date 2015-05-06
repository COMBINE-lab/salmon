#ifndef ITERATIVE_OPTIMIZER_HPP
#define ITERATIVE_OPTIMIZER_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <map>
#include <vector>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <sstream>
#include <exception>
#include <random>
#include <queue>
#include "btree_map.h"

/** Boost Includes */
#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>

//#include <Eigen/Core>
#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/task_scheduler_init.h"

//#include "nnls.h"
#include "BiasIndex.hpp"
#include "poisson_solver.hpp"
#include "matrix_tools.hpp"
#include "ezETAProgressBar.hpp"
//#include "CompatibilityGraph.hpp"


template <typename ReadHashT, typename TranscriptHashT>
class IterativeOptimizer {

private:
    /**
    * Typedefs
    */
    // A KmerMapT is a map from a kmer (encoded as an integer) to the set
    // of transcripts where it occurs
    typedef uint32_t TranscriptIDT;
    typedef uint64_t KmerIDT;
    typedef double KmerQuantityT;
    typedef double PromiscutityT;
    //typedef std::unordered_map< KmerIDT, std::vector<TranscriptIDT> > KmerMapT;
    typedef tbb::concurrent_unordered_map< uint64_t, tbb::concurrent_vector<uint32_t> > KmerMapT;

    struct TranscriptGeneVectors;
    typedef std::vector<TranscriptGeneVectors> KmerIDMap;


    typedef std::tuple<TranscriptIDT, std::vector<KmerIDT>> TranscriptKmerSet;
    typedef std::string *StringPtrT;
    typedef uint64_t TranscriptScoreT;
    typedef jellyfish::invertible_hash::array<uint64_t, atomic::gcc, allocators::mmap> HashArrayT;
    typedef size_t ReadLengthT;

    // Necessary forward declaration
    struct TranscriptDataT;
    typedef std::tuple<TranscriptScoreT, TranscriptIDT> HeapPair;
    typedef typename boost::heap::fibonacci_heap<HeapPair>::handle_type HandleT;

    struct TranscriptGeneVectors {
        tbb::concurrent_vector<uint32_t> transcripts;
        tbb::concurrent_vector<uint32_t> genes;
    };

    struct TranscriptDataT {
        TranscriptIDT id;
        StringPtrT header;
        std::map<KmerIDT, KmerQuantityT> binMers;
        KmerQuantityT mean;
        size_t length;
    };

    struct TranscriptInfo {
        btree::btree_map<KmerIDT, KmerQuantityT> binMers;
        KmerQuantityT mean;
        ReadLengthT length;
        ReadLengthT effectiveLength;
    };

    // This struct represents a "job" (transcript) that needs to be processed
    struct TranscriptJob {
        StringPtrT header;
        StringPtrT seq;
        TranscriptIDT id;
    };

    struct TranscriptResult {
        TranscriptDataT *data;
        TranscriptKmerSet *ks;
    };

    struct BinmerUpdates {
        std::vector<KmerIDT> zeroedBinmers;
        std::vector<KmerIDT> updatedBinmers;
    };

    size_t _merLen;
    ReadHashT &_readHash;
    TranscriptHashT &_transcriptHash;
    BiasIndex& biasIndex_;

    // The number of occurences above whcih a kmer is considered promiscuous
    size_t _promiscuousKmerCutoff {15};

    // Map each kmer to the set of transcripts it occurs in
    //KmerMapT _transcriptsForKmer;
    KmerIDMap _transcriptsForKmer;

    // The actual data for each transcript
    std::vector<TranscriptInfo *> _transcripts;

    TranscriptGeneMap& _transcriptGeneMap;

    tbb::concurrent_unordered_set<uint64_t> _genePromiscuousKmers;

    std::vector<PromiscutityT> kmerPromiscuities_;
    std::vector<PromiscutityT> kmerBiases_;


    /******** REFACTOR ME *************/
    void extractKmerTranscriptCompatibilityGraph() {
        std::vector<int> kmerComponents(_transcriptsForKmer.size(), -1);
        std::vector<int> transcriptComponents(_transcripts.size(), -1);

        int component = 0;
        size_t tstart = 0;
        size_t numAssigned = 0;
        std::queue<size_t> toProcess;

        while (numAssigned < transcriptComponents.size()) {
            // Find the first unassigned transcript
            while ( transcriptComponents[tstart] != -1 and tstart < transcriptComponents.size()) { 
                ++tstart; 
            }

            if ( tstart == transcriptComponents.size() ) { std::cerr << "no unassigned transcripts\n"; }

            toProcess.push(tstart);
            while( !toProcess.empty() ) {
                auto tid = toProcess.front();
                toProcess.pop();
                if ( transcriptComponents[tid] == -1 ) {
                    transcriptComponents[tid] = component;    
                    ++numAssigned;

                    // Find all kmers in this transcript
                    for( auto binmer : _transcripts[tid]->binMers ) {

                        auto existingComponent = kmerComponents[binmer.first];

                        if ( existingComponent == -1 ) {
                            kmerComponents[binmer.first] = component;
                            for( auto ts : _transcriptsForKmer[binmer.first].transcripts ) {
                                if ( transcriptComponents[ts] == -1 ) {
                                    toProcess.push(ts);
                                } else {
                                    assert(transcriptComponents[ts] == component);
                                }
                            }
                        }
                    }

                } else {
                    assert( transcriptComponents[tid] == component );
                }                    
            }
            ++component;
        }

        btree::btree_map<int, size_t> components;
        for( auto tc : transcriptComponents ) {
            assert(tc != -1);
            components[tc] += 1;
            //components.insert(tc);
        }
        for( auto kc : kmerComponents ) {
            assert(kc != -1);
        }
        for ( auto cit : components ) {
            std::cerr << "component " << cit.first << " contains " << cit.second << " transcripts\n";
        }
        std::cerr << "there were " << components.size() << " components\n";
    }


    /******** END REFACTOR ME *************/


    inline double _idf( uint64_t k ) {
        double df = _transcriptsForKmer[k].transcripts.size();
        return (df > 0.0) ? std::log(_transcripts.size() / df) : 0.0;
    }

    // Should the given kmer be considered?
    inline bool _considered( uint64_t mer ) {
        // The kmer is only considered if it exists in the transcript set
        // (i.e. it's possible to cover) and it's less prmiscuous than the
        // cutoff.
        return ( _transcriptHash.atIndex(mer) > 0 and
                 _transcriptHash.atIndex(mer) < _promiscuousKmerCutoff );
    }

    KmerQuantityT _weight( KmerIDT k ) {
        //return 1.0 / (_transcriptHash[k]);
        return 1.0 / (_transcriptHash.atIndex(k) );
    }

    KmerQuantityT _computeMedian( TranscriptInfo* ti ) {

        KmerQuantityT median = 0.0;
        auto& binMers = ti->binMers;
        auto len = binMers.size();
        if ( len > 0 ) {
            if ( len % 2 == 0 ) {
                auto it = binMers.begin();
                for ( size_t i = 0; i < (len / 2); ++i, ++it ) {}
                median = it->second; ++it;
                median += it->second;
                median /= 2.0;
            } else {
                auto it = binMers.begin();
                for ( size_t i = 0; i < (len / 2) + 1; ++i, ++it ) {}
                median = it->second;
            }
        }
        return median;
    }

    /**
     * Computes the sum of kmer counts within the transcript given by ti, but clamping
     * all non-zero counts to the given quantile.  For example, if quantile was 0.25, and
     * x and y represented the 1st and 3rd quantile of kmer counts, then every nonzero count c 
     * would be transformed as c = max(x, min(y,c));
     * 
     * @param  ti       [description]
     * @param  quantile [description]
     * @return          [description]
     */
    KmerQuantityT _computeSumQuantile( TranscriptInfo* ti, double quantile ) {
        using namespace boost::accumulators;
        typedef accumulator_set<double, stats<tag::p_square_quantile> > accumulator_t;
        KmerQuantityT sum = 0.0;
        
        accumulator_t accLow(quantile_probability = quantile);
        accumulator_t accHigh(quantile_probability = 1.0-quantile);
        for ( auto binmer : ti->binMers ) {
            if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
                accLow(binmer.second);
                accHigh(binmer.second);
            }        
        }

        auto cutLow = p_square_quantile(accLow);
        auto cutHigh = p_square_quantile(accHigh);

        for ( auto binmer : ti->binMers ) {
            if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
                sum += std::min( cutHigh, std::max( cutLow, binmer.second ) );
            }
        }
        return sum;
    }

    KmerQuantityT _computeSum( TranscriptInfo* ti ) {
        KmerQuantityT sum = 0.0;
        for ( auto binmer : ti->binMers ) {
            if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
                sum += kmerBiases_[binmer.first] * binmer.second;
            }
        }
        return sum;
    }

    bool _discard(TranscriptInfo* ti) {
        if ( ti->mean == 0.0 ) { 
            return false; 
        } else {
            ti->mean = 0.0;
            ti->binMers.clear();
            return true;
        }
    }

    KmerQuantityT _computeMean( TranscriptInfo* ti ) {
        //return (ti->effectiveLength > 0.0) ? _computeSumQuantile(ti, 0.25) / ti->effectiveLength : 0.0;
        return (ti->effectiveLength > 0.0) ? (_computeSum(ti) / ti->effectiveLength) : 0.0;
    }

    KmerQuantityT _computeWeightedMean( TranscriptInfo* ti ) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::count, tag::weighted_mean>, double> acc;

        for ( auto binmer : ti->binMers ) {
          if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
            acc(binmer.second, weight=kmerBiases_[binmer.first] * kmerPromiscuities_[binmer.first]);
          }
        }

        auto nnz = count(acc);
        
        if ( nnz < ti->effectiveLength ) {
            acc(0.0, weight=ti->effectiveLength-nnz);
        }
       
        auto sum = sum_of_weights(acc);
        return sum > 0.0 ? weighted_mean(acc) : 0.0;
    }

    double _effectiveLength( TranscriptInfo* ts ) {
        double length = 0.0;
        for ( auto binmer : ts->binMers ) {
            length += _weight(binmer.first);
        }
        return length;
    }

    void normalizeTranscriptMeans_(){
        auto sumMean = 0.0;
        for ( auto ti : _transcripts ) { sumMean += ti->mean; }
    
        // compute the new mean for each transcript
        tbb::parallel_for( size_t(0), size_t(_transcripts.size()),
            [this, sumMean]( size_t tid ) -> void { this->_transcripts[tid]->mean /= sumMean; });

    }

    double averageCount(TranscriptInfo* ts){
        if ( ts->binMers.size() == 0 ) { return 0.0; }
        double sum = 0.0;
        for ( auto binmer : ts->binMers ) {
            sum += kmerBiases_[binmer.first] * binmer.second;
        }
        return sum / ts->binMers.size();

    }

    /**
     * This function should be run only after <b>after</b> an EM loop.
     * It estimates kmer specific biases based on how much a kmer's count
     * deviates from it's corresponding transcript's mean 
     */
    void computeKmerFidelities_() {

        std::vector<double> transcriptFidelities(_transcripts.size(), 0.0);
        tbb::parallel_for( size_t(0), _transcripts.size(),
            [this, &transcriptFidelities]( TranscriptIDT tid ) {
                double sumDiff = 0.0;
                auto ts = this->_transcripts[tid];
                for ( auto& b : ts->binMers ) {
                    auto diff = (this->kmerBiases_[b.first] * b.second) - ts->mean;
                    //auto diff = (this->kmerBiases_[b.first] * b.second) - this->averageCount(ts);
                    sumDiff += diff*diff;
                }
                transcriptFidelities[tid] = std::sqrt(sumDiff / ts->binMers.size());
            });
        
        auto totalFidelity = std::accumulate(transcriptFidelities.begin(), transcriptFidelities.end(), 0.0);
        auto maxFidelity = *std::max_element(transcriptFidelities.begin(), transcriptFidelities.end());
        auto averageFidelity = totalFidelity / transcriptFidelities.size();

        std::cerr << "max fidelity = " << maxFidelity << "\n";

        tbb::parallel_for( size_t(0), _transcriptsForKmer.size(),
            [this, &transcriptFidelities, averageFidelity, maxFidelity]( KmerIDT kid ) -> void {
                double sumT = 0.0; double sumK = 0.0;
                double confidence = 0.0;
                for( auto tid : this->_transcriptsForKmer[kid].transcripts ) {
                    sumT += this->_transcripts[tid]->mean;
                    sumK += this->_transcripts[tid]->binMers[kid];   
                    confidence += transcriptFidelities[tid];// / averageFidelity;
                }

                confidence /= this->_transcriptsForKmer[kid].transcripts.size();
                double alpha = 0.5; //std::min( 1.0, 10.0*averageFidelity / confidence);
                double bias = sumK > 0.0 ? (sumT / sumK) : 0.0;
                double prevBias = this->kmerBiases_[kid];
                this->kmerBiases_[kid] = alpha * bias + (1.0 - alpha) * prevBias;
            }
        );

    }

    size_t _initialize( const std::vector<std::string> &transcriptFiles ) { 

        char** fnames = new char*[transcriptFiles.size()];
        size_t z{0};
        size_t numFnames{0};
        for ( auto& s : transcriptFiles ){
            // Ugly, yes?  But this is not as ugly as the alternatives.
            // The char*'s contained in fname are owned by the transcriptFiles
            // vector and need not be manually freed.
            fnames[numFnames] = const_cast<char*>(s.c_str());
            ++numFnames;
        }

        // Create a jellyfish parser
        jellyfish::parse_read parser( fnames, fnames+numFnames, 1000);

        // So we can concisely identify each transcript
        TranscriptIDT transcriptIndex {0};

        size_t numActors = 12;
        std::vector<std::thread> threads;

        _transcripts.resize( _transcriptGeneMap.numTranscripts(), nullptr );

        tbb::parallel_for( size_t{0}, _transcriptGeneMap.numTranscripts(),
            [this](size_t tid) -> void {
             auto procRead = new TranscriptInfo {
                btree::btree_map<KmerIDT, KmerQuantityT>(),
                0.0,
                0,
                0
                };
              this->_transcripts[tid] = procRead;
            }
        );

        _transcriptsForKmer.resize( _transcriptHash.size() );
        kmerBiases_.resize(_transcriptHash.size(), 1.0);

        bool done {false};
        std::atomic<size_t> numRes {0};
        std::atomic<size_t> nworking{numActors-1};

        // Start the thread that will print the progress bar
        threads.push_back( std::thread( [&numRes, &nworking, this] () {
            auto numJobs = this->_transcriptGeneMap.numTranscripts();
            ez::ezETAProgressBar show_progress(numJobs);
            //boost::progress_display show_progress( numJobs );
            size_t lastCount = numRes;
            show_progress.start();
            while ( lastCount < numJobs and nworking > 0 ) {
                auto diff = numRes - lastCount;
                if ( diff > 0 ) {
                    show_progress += static_cast<unsigned int>(diff);
                    lastCount += diff;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(1010));
            }
        }) );

        std::cerr << "Processing transcripts\n";

        // Start the desired number of threads to parse the transcripts
        // and build our data structure.
        for (size_t i = 0; i < numActors - 1; ++i) {

            threads.push_back( std::thread( 
                [&numRes, &parser, &nworking, this]() -> void {

                    // Each thread gets it's own stream
                    jellyfish::parse_read::read_t* read;
                    jellyfish::parse_read::thread stream = parser.new_thread();

                    while ( (read = stream.next_read()) ) {
                        ++numRes;
                        
                        std::string header(read->header, read->hlen);
                        size_t d = std::min(header.length()-1, header.find(' '));
                        std::string tname( header.substr(0,d) );

                        auto biases = this->biasIndex_.getBiases(tname);

                        // Strip the newlines from the read and put it into a string (newSeq)
                        std::string seq(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );
                        auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
                        auto readLen = std::distance( seq.begin(), newEnd );
                        std::string newSeq(seq.begin(), seq.begin() + readLen);

                        // The set of kmers in this transcript
                        //auto ts = new TranscriptKmerSet {transcriptIndex, {}};
                        // We know how many kmers there will be, so allocate the space
                        /*
                        auto procRead = new TranscriptInfo {
                            btree::btree_map<KmerIDT, KmerQuantityT>(),
                            0.0,
                            readLen,
                            0
                        };
                        */
                        
                        // Lookup the ID of this transcript in our transcript -> gene map
                        auto transcriptIndex = this->_transcriptGeneMap.findTranscriptID(tname); 
                        bool valid = ((transcriptIndex != this->_transcriptGeneMap.INVALID) and 
                                      (readLen > this->_merLen));
                        
                        if ( not valid ) {
                            //this->_transcripts[transcriptIndex] = procRead;
                            continue; 
                        }
                        
                        auto procRead = this->_transcripts[transcriptIndex];
                        procRead->length = readLen;
                        auto geneIndex = this->_transcriptGeneMap.gene(transcriptIndex);

                        size_t numKmers {readLen - this->_merLen + 1};
                        // The binary representation of the kmers in this transcript
                        //btree::btree_map<KmerIDT, KmerQuantityT> binMers;

                        // Iterate over the kmers
                        ReadLengthT effectiveLength(0);
                        float weightedLen = 0.0;
                        size_t coverage = _merLen;
                        auto INVALID = this->_transcriptHash.INVALID;
                        for ( auto offset : boost::irange( size_t(0), numKmers ) ) { 
                            // the kmer and it's uint64_t representation
                            try {
                                auto mer = newSeq.substr( offset, this->_merLen );                            
                            auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), _merLen );
                            auto rcBinMer = jellyfish::parse_dna::reverse_complement( binMer, _merLen );
                            binMer = std::min(binMer, rcBinMer);

                            auto binMerId = this->_transcriptHash.id(binMer);
                            // Only count and track kmers which should be considered
                            //if ( this->_considered(binMer) ) {
                            if ( binMerId != INVALID and this->_considered(binMerId) ) {
                                weightedLen += (1.0 / this->_transcriptHash.atIndex(binMerId));
                                effectiveLength += coverage;
                                coverage = 1;
                                //binMers[binMer] += 1.0;
                                //procRead->binMers[binMerId] += 1.0;
                                procRead->binMers[binMerId] += 1.0 / biases[offset];
                                
                                this->_transcriptsForKmer[binMerId].transcripts.push_back(transcriptIndex);
                                this->_transcriptsForKmer[binMerId].genes.push_back(geneIndex);
                            } else {
                                coverage += (coverage < _merLen) ? 1 : 0;
                            }

                            } catch ( std::exception& e ) { 
                                std::cerr << "trying to take kmer starting at " << offset << 
                                ", length of read is " << readLen << " merLen is " << this->_merLen << "\n";
                                std::cerr << "read is " << newSeq << "\n";
                                std::cerr << e.what() << "\n";
                                std::exit(1);
                            }

                        }
                        procRead->effectiveLength = numKmers;

                        /*
                        auto procRead = new TranscriptInfo {
                            binMers, // the set of kmers
                            0.0, // initial desired mean
                            readLen, // real transcript length
                            effectiveLength // effective transcript length
                        };*/

                        //this->_transcripts[transcriptIndex] = procRead;
                    }
                    nworking--;
            }) );

        }

        // Wait for all of the threads to finish
        for ( auto& thread : threads ){ thread.join(); }

        numRes = 0;
        tbb::task_scheduler_init tbb_init;

        std::cerr << "\n\nDetermining gene-promiscuous kmers ... ";
        
        //tbb::parallel_for_each( _transcriptsForKmer.begin(), _transcriptsForKmer.end(),
        tbb::parallel_for( size_t(0), size_t(_transcriptsForKmer.size()),
            [&numRes, this]( size_t idx ) { //KmerMapT::value_type& kt ) {

                /*
                //auto& transcripts = kt.second;
                */
                
                // Uniqify the transcripts
                auto& transcripts = this->_transcriptsForKmer[idx].transcripts;
                std::sort(transcripts.begin(), transcripts.end());
                auto it = std::unique(transcripts.begin(), transcripts.end()); 
                transcripts.resize( it - transcripts.begin() );
                
                bool filterGenePromiscuous = false;

                if (filterGenePromiscuous) {
                    auto numGenes = this->_transcriptsForKmer[idx].genes.size();
                    if (numGenes > 1) {
                        auto geneId = this->_transcriptsForKmer[idx].genes[0];
                        for ( auto i : boost::irange(size_t(1), this->_transcriptsForKmer[idx].genes.size()) ) {
                          auto id = this->_transcriptsForKmer[idx].genes[i];
                          if( id != geneId) { 
                            this->_genePromiscuousKmers.insert(idx);
                            break; 
                          }
                        }
                    }
                }
    

                // auto geneId = this->_transcriptGeneMap.gene(transcripts[0]);
                // for ( auto i : boost::irange(size_t(1), transcripts.size()) ) {
                //     auto id = this->_transcriptGeneMap.gene(transcripts[i]);
                //     if( id != geneId) { this->_genePromiscuousKmers.insert(kt.first); break; }
                // }

                /*
                std::unordered_set<uint32_t> geneSet;
                for ( auto t : transcripts ){
                    geneSet.insert(this->_transcriptGeneMap.gene(t));
                }
                if ( geneSet.size() > 1 ){ this->_genePromiscuousKmers.insert(kt.first); }
                */
                ++numRes;
            }
        );

        auto gpfrac = _genePromiscuousKmers.size() / static_cast<double>(_transcriptsForKmer.size());
        std::cerr << "There were " << _genePromiscuousKmers.size() << " (" << 100 * gpfrac << "%) gene promiscuous kmers\n";

        std::cerr << "Computing promiscuity rates\n";
        kmerPromiscuities_.resize(_transcriptsForKmer.size());
        tbb::parallel_for( size_t{0}, kmerPromiscuities_.size(),
            [this]( KmerIDT kid ) -> void { this->kmerPromiscuities_[kid] = this->_weight(kid); }
        );

        /**
         * gene-promiscuous kmers can never be added to a transcript's counts, so
         * it's unfair to consider them in the transcripts effective length. 
         */
        std::for_each( _genePromiscuousKmers.begin(), _genePromiscuousKmers.end(),
            [this]( KmerIDT kmerId ) { 
                for ( auto tid : _transcriptsForKmer[kmerId].transcripts ) {
                    _transcripts[tid]->effectiveLength -= 1.0;
                }
            });

        std::cerr << "done\n";

        std::cerr << "computing number of mapped (usable) reads\n";
        size_t mappedReads = 0;
        for ( auto kidx : boost::irange(size_t(0), _transcriptsForKmer.size()) ) {
            if ( _genePromiscuousKmers.find(kidx) == _genePromiscuousKmers.end() ) {
                mappedReads += _readHash.atIndex(kidx);
            }
        }

        //std::cerr << "computing kmer compatibility graph\n";
        //extractKmerTranscriptCompatibilityGraph();
        //std::cerr << "done\n";
        return mappedReads;
    }

    void _dumpCoverage( const std::string &cfname ) {
        typedef std::string* StringPtr;

        size_t numTrans = _transcripts.size();
        size_t numProc = 0;
        std::ofstream ofile(cfname);

        ofile << "# num_transcripts\n";
        ofile << "# transcript_name_{1} num_kmers_{1} count_1 count_2 ... count_{num_kmers}\n";
        ofile << "# ... \n";
        ofile << "# transcript_name_{num_transcripts} num_kmers_{num_transcripts} count_1 count_2 ... count_{num_kmers_{num_transcripts}}\n";

        ofile << _transcripts.size() << "\n";

        std::cerr << "Dumping coverage statistics to " << cfname << "\n";

        boost::lockfree::queue<StringPtr> covQueue(_transcripts.size());
        
        tbb::parallel_for( size_t{0}, _transcripts.size(),
            [this, &covQueue] (size_t index) -> void {

                auto td = this->_transcripts[index];
                
                std::stringstream ostream;
                ostream << this->_transcriptGeneMap.transcriptName(index) << " " << td->binMers.size();
                for ( auto bm : td->binMers ) {
                    ostream << " " << bm.second;
                }
                ostream << "\n";
                std::string* ostr = new std::string(ostream.str());
                while(!covQueue.push(ostr));
            }
        );


                        
        ez::ezETAProgressBar pb(_transcripts.size());
        pb.start();

        std::string* sptr = nullptr;
        while ( numProc < numTrans ) {
            while( covQueue.pop(sptr) ) {
                ofile << (*sptr);
                ++pb;
                ++numProc;
                delete sptr;
            }
        }

        ofile.close();

    }



    std::vector<double> _estimateIsoformNNLS( const size_t geneID, const size_t mappedReads ) {

        typedef std::vector<std::vector<double>> IsoOptMatrix;
        typedef std::vector<double> IsoOptVector;
    
        auto transcripts = _transcriptGeneMap.transcriptsForGene( geneID );
        size_t numTranscripts = transcripts.size();

        std::unordered_map<uint64_t, size_t> kmerIDs;
        std::unordered_set<uint64_t> geneKmers;

        // reindex the kmers for this gene
        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kc : _transcripts[ transcripts[i] ]->binMers ) {
                auto kmer = kc.first;
                auto kmerIt = geneKmers.find(kmer);
                
                if ( kmerIt == geneKmers.end() ) {
                    kmerIDs[kmer] = geneKmers.size();
                    geneKmers.insert(kmer);
                }
            }            
        }


        size_t numKmers = geneKmers.size();
        if ( numKmers == 0 ) {
            return IsoOptVector();
        }


        // Create a numTranscripts x numKmers matrix to hold the "rates"
        IsoOptVector x(numTranscripts, 0.0);
        IsoOptMatrix A(numTranscripts, IsoOptVector(numKmers, 0.0) );
        IsoOptVector counts(numKmers, 0.0);
        double numReads = 0.0;

        std::vector<double> numMappedReads(transcripts.size(), 0.0);

        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kc : _transcripts[ transcripts[i] ]->binMers ) {
                auto kmer = kc.first; auto count = kc.second;
                size_t j = kmerIDs[kmer];
                A[i][j] = count;
                numMappedReads[i] += count;
                if ( counts[j] == 0 ) {
                    auto promIt = _genePromiscuousKmers.find(kmer);
                    double normFact = (promIt == _genePromiscuousKmers.end() ) ? 1.0 : 0.0;
                    counts[j] = _readHash.atIndex(kmer) * normFact;
                    numReads += counts[j];
                }
            }
        }


        std::unique_ptr<IsoOptMatrix> collapsedA;
        std::unique_ptr<IsoOptVector> collapsedCounts;
        matrix_tools::collapseIntoCategories(A, counts, collapsedA, collapsedCounts);

        std::vector<double> cc(collapsedCounts->size(), 0.0);
        for( size_t i = 0; i < cc.size(); ++i ) { cc[i] = (*collapsedCounts)[i]; }

        auto nr = collapsedA->size();
        auto nc = (*collapsedA)[0].size();

        typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
        typedef Eigen::Triplet<double> T;

        auto M = std::unique_ptr<SpMat>(new SpMat(nr,nc));
        std::vector<T> trips;

        for( size_t r = 0; r < nr; ++r ) {
            for( size_t c = 0; c < nc; ++c ) {
                if ( (*collapsedA)[r][c] > 0 ) {
                    trips.push_back( T(r,c,(*collapsedA)[r][c]) );
                }
            }
        }
        M->setFromTriplets(trips.begin(), trips.end());
                
        /*
        auto dupCols = matrix_tools::markDependentColumns(*M);
        if (dupCols.size() > 0) {
            std::vector<T> triplets;

            auto mCol = M->cols() - dupCols.size();
            size_t dupPtr = 0;
            for (size_t k = 0; k < M->outerSize(); ++k) {
                // If this is a colum we skip, then continue this loop and
                // move to the next column to skip
                if ( dupPtr < dupCols.size() && k == dupCols[dupPtr] ) {
                    ++dupPtr;
                    continue;
                }
                // The column index is the column index of A minus the number of
                // columns we've skipped so far
                size_t colInd = k - dupPtr;
                assert(colInd < mCol);
                // Add all rows of this column to the new matrix
                for (SpMat::InnerIterator it(*M, k); it; ++it) {
                    triplets.push_back( {it.row(), colInd, it.value()} );
                }
            }
            assert(M->cols() > dupCols.size());
            auto B = std::unique_ptr<SpMat>(new SpMat(M->rows(), M->cols() - dupCols.size()));
            B->setFromTriplets(triplets.begin(), triplets.end());
            std::swap(M, B);
        }
       */ 
       
        auto r = matrix_tools::nnlsSolve(*M, cc);
        for( auto i : boost::irange(size_t{0}, r.size()) ){
            r[i] *= numMappedReads[i];
        }

        return r;
    }



    std::vector<double> _estimateIsoformAbundance( const size_t geneID, const size_t mappedReads
      //std::vector<std::vector<uint64_t>> &transcriptKmers,
      //std::vector<double>& transcriptLengths,
      //size_t mappedReads,
      //tbb::concurrent_unordered_map<uint64_t, uint32_t>& genePromiscuousKmers
      ) {

        typedef std::vector<std::vector<double>> IsoOptMatrix;
        typedef std::vector<double> IsoOptVector;
	
        auto transcripts = _transcriptGeneMap.transcriptsForGene( geneID );
        size_t numTranscripts = transcripts.size();

        std::unordered_map<uint64_t, size_t> kmerIDs;
        std::unordered_set<uint64_t> geneKmers;
        // reindex the kmers for this gene
        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kc : _transcripts[ transcripts[i] ]->binMers ) {
                auto kmer = kc.first;
                auto kmerIt = geneKmers.find(kmer);
                if ( kmerIt == geneKmers.end() ) {
                    kmerIDs[kmer] = geneKmers.size();
                    geneKmers.insert(kmer);
                }
            }            
        }


	size_t numKmers = geneKmers.size();
        if ( numKmers == 0 ) {
            return IsoOptVector();
        }


	// Create a numTranscripts x numKmers matrix to hold the "rates"
        IsoOptVector x(numTranscripts, 0.0);
        IsoOptMatrix A(numTranscripts, IsoOptVector(numKmers, 0.0) );
        IsoOptVector counts(numKmers, 0.0);
        double numReads = 0.0;

        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kc : _transcripts[ transcripts[i] ]->binMers ) {
                auto kmer = kc.first; auto count = kc.second;
                size_t j = kmerIDs[kmer];
                A[i][j] = mappedReads;//1.0;
                if ( counts[j] == 0 ) {
                    auto promIt = _genePromiscuousKmers.find(kmer);
                    double normFact = (promIt == _genePromiscuousKmers.end() ) ? 1.0 : 0.0;
                    counts[j] = _readHash.atIndex(kmer) * normFact;
                    numReads += counts[j];
                }
            }
        }


        std::unique_ptr<IsoOptMatrix> collapsedA;
        std::unique_ptr<IsoOptVector> collapsedCounts;
        matrix_tools::collapseIntoCategories(A, counts, collapsedA, collapsedCounts);
        
        for (size_t i = 0; i < collapsedA->size(); ++i) {
            auto tlen = _transcripts[ transcripts[i] ]->length;
            for (size_t j = 0; j < (*collapsedA)[0].size(); ++j) {
                (*collapsedA)[i][j] = (*collapsedA)[i][j] / ( tlen * mappedReads / 1000000000.0 );
                //(*collapsedA)[i][j] = (*collapsedA)[i][j] * tlen / 1000.0 * static_cast<double>(mappedReads) / 1000000.0;  
            }
        }

        solve_likelihood( *collapsedA, *collapsedCounts, x, false );
        return x;
    }

public:
    /**
     * Construct the solver with the read and transcript hashes
     */
    IterativeOptimizer( ReadHashT &readHash, TranscriptHashT  &transcriptHash, TranscriptGeneMap& transcriptGeneMap,
                        BiasIndex& biasIndex ) : 
        _readHash(readHash), _transcriptHash(transcriptHash), _merLen(transcriptHash.kmerLength()), 
        _transcriptGeneMap(transcriptGeneMap), biasIndex_(biasIndex) {}


    KmerQuantityT optimize( const std::vector<std::string> &transcriptFiles, const std::string &outputFile, size_t numIt, double minMean ) {

        _initialize(transcriptFiles);

        // Holds results that we will use to update the index map
        size_t numActors = 1;

        KmerQuantityT globalError {0.0};
        bool done {false};
        std::atomic<size_t> numJobs {0};
        std::atomic<size_t> completedJobs {0};
        std::vector<KmerIDT> kmerList( _transcriptsForKmer.size(), 0 );
        size_t idx = 0;

        tbb::task_scheduler_init tbb_init;


        std::cerr << "Computing initial coverage estimates ... ";

        std::default_random_engine generator;
        std::normal_distribution<double> distribution(5.0,2.0);

        tbb::parallel_for( size_t(0), size_t(_transcriptGeneMap.numTranscripts()),
        [this, &distribution, &generator]( size_t tid ) -> void {
            auto transcriptData = this->_transcripts[tid];

            for ( auto & kv : transcriptData->binMers ) {
                auto kmer = kv.first;
                if ( this->_genePromiscuousKmers.find(kmer) == this->_genePromiscuousKmers.end() ){
                    // count is the number of times kmer appears in transcript (tid)
                    auto count = kv.second;
                    kv.second = count * this->_readHash.atIndex(kmer) * this->_weight(kmer);
                }
            }
            transcriptData->mean = this->_computeMean( transcriptData );
            //transcriptData->mean = distribution(generator);
            //this->_computeWeightedMean( transcriptData );
        }
        );

        //normalizeTranscriptMeans_();

        std::cerr << "done\n";

        size_t outerIterations = 1;
        for ( size_t oiter = 0; oiter < outerIterations; ++oiter ) {

        for ( size_t iter = 0; iter < numIt; ++iter ) {

            auto reqNumJobs = _transcriptsForKmer.size();

            std::cerr << "iteraton: " << iter << "\n";

            globalError = 0.0;
            numJobs = 0;
            completedJobs = 0;

            auto pbthread = std::thread( 
                [&completedJobs, reqNumJobs]() -> bool {
                    auto prevNumJobs = 0;
                    ez::ezETAProgressBar show_progress(reqNumJobs);
                    //boost::progress_display show_progress( reqNumJobs );
                    show_progress.start();
                    while ( prevNumJobs < reqNumJobs ) {
                        if ( prevNumJobs < completedJobs ) {
                            show_progress += completedJobs - prevNumJobs;
                        }
                        prevNumJobs = completedJobs.load();

                        boost::this_thread::sleep_for(boost::chrono::seconds(1));
                    }
                    show_progress.done();
                    return true;
                }
            );

            tbb::parallel_for( size_t(0), size_t(_transcriptsForKmer.size()),
                // for each kmer    
                [&kmerList, &completedJobs, this]( size_t kid ) {

                        auto kmer = kid;
                        if ( this->_genePromiscuousKmers.find(kmer) == this->_genePromiscuousKmers.end() ){
                        
                            auto &transcripts = this->_transcriptsForKmer[kmer].transcripts;
                            if ( transcripts.size() > 1 ) {

                                // for each transcript containing this kmer
                                double totalMass = 0.0;
                                for ( auto tid : transcripts ) {
                                    totalMass += this->_transcripts[tid]->mean;
                                }


                                if ( totalMass > 0.0 ) {
                                    double norm = 1.0 / totalMass;
                                    for ( auto tid : transcripts ) {
                                        if ( this->_transcripts[tid]->mean > 0.0 ) {
                                            this->_transcripts[tid]->binMers[kmer] =
                                            this->_transcripts[tid]->mean * norm * kmerBiases_[kmer] * this->_readHash.atIndex(kmer);
                                        }
                                    }
                                }

                            }
                        }
                        ++completedJobs;
                    }
            );

            pbthread.join();
            //tp.wait();

            // reset the job counter
            completedJobs = 0;

            double delta = 0.0;
            double norm = 1.0 / _transcripts.size();

            std::vector<KmerQuantityT> prevMeans( _transcripts.size(), 0.0 );
            tbb::parallel_for( size_t(0), size_t(_transcripts.size()),
                [this, &prevMeans]( size_t tid ) -> void { prevMeans[tid] = this->_transcripts[tid]->mean; });

            std::cerr << "\ncomputing new means ... ";
            size_t discard = 0;
            // compute the new mean for each transcript
            tbb::parallel_for( size_t(0), size_t(_transcripts.size()),
                [this, iter, numIt, norm, minMean, &discard]( size_t tid ) -> void {
                        // this thread claimed myJobID;
                        auto ts = this->_transcripts[tid];
                        auto tsNorm = 1.0;//(ts->effectiveLength > 0.0) ? 1.0 / std::sqrt(ts->effectiveLength) : 1.0;
                        //ts->mean = tsNorm * this->_computeWeightedMean( ts );
                        //ts->mean = tsNorm * this->averageCount( ts );
                        ts->mean = tsNorm * this->_computeMean( ts );
                }
            );

            //normalizeTranscriptMeans_();
            for( auto tid : boost::irange(size_t{0}, prevMeans.size()) ){
                delta += std::abs( _transcripts[tid]->mean - prevMeans[tid] );
            }

            std::cerr << "done\n";
            std::cerr << "total variation in mean = " << delta << "\n";
            std::cerr << "discarded " << discard << " transcripts in this round whose mean was below " << minMean << "\n";
        }

        std::cerr << "end of outer iteration " << oiter << " recomputing biases\n";
        computeKmerFidelities_();
    }

        std::cerr << "Writing output\n";
        ez::ezETAProgressBar pb(_transcripts.size());
        pb.start();

        std::ofstream ofile( outputFile );
        size_t index = 0;
        ofile << "Transcript" << '\t' << "Length" << '\t' << "Effective Length" << '\t' << "Weighted Mapped Reads" << '\n';
        for ( auto ts : _transcripts ) {
            ofile << _transcriptGeneMap.transcriptName(index) << '\t' << ts->length << '\t' <<
                    ts->effectiveLength << '\t' << _computeSum(ts) << "\n";
            ++index;
            ++pb;
        }
        ofile.close();

        auto writeCoverageInfo = true;
        if ( writeCoverageInfo ) {
            std::string cfname("transcriptCoverage.txt");
            _dumpCoverage( cfname );
        }

    }

    KmerQuantityT optimizePoisson( const std::vector<std::string> &transcriptFiles, const std::string &outputFile ) {

        auto mappedReads = _initialize(transcriptFiles);

        // Holds results that we will use to update the index map
        size_t numActors = 1;

        KmerQuantityT globalError {0.0};
        bool done {false};
        std::atomic<size_t> numJobs {0};
        std::atomic<size_t> completedJobs {0};
        std::mutex nnlsmutex;
        std::vector<KmerIDT> kmerList( _transcriptsForKmer.size(), 0 );
        size_t idx = 0;


        std::vector<double> abundances( _transcriptGeneMap.numTranscripts(), 0.0 );

        _transcriptGeneMap.needReverse();
        std::atomic<size_t> numGenesProcessed{0};

        auto abundanceEstimateProgress = std::thread( [&numGenesProcessed, this] () -> void {
            auto numJobs = this->_transcriptGeneMap.numGenes();
            ez::ezETAProgressBar pbar(numJobs);
            pbar.start();
            size_t lastCount = numGenesProcessed;
            while ( lastCount < numJobs ) {
                auto diff = numGenesProcessed - lastCount;
                if ( diff > 0 ) {
                    pbar += diff;
                    lastCount += diff;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(1010));
            }
        });
        std::cerr << "Processing abundances\n";

        // process each gene and solve the estimation problem
        tbb::parallel_for( size_t {0}, _transcriptGeneMap.numGenes(), 
            [this, mappedReads, &abundances, &numGenesProcessed]( const size_t & geneID ) {
                //auto abundance = this->_estimateIsoformAbundance(geneID, mappedReads);
                auto abundance = this->_estimateIsoformNNLS(geneID, mappedReads);
                auto transcripts = this->_transcriptGeneMap.transcriptsForGene(geneID);

                for ( size_t i = 0; i < abundance.size(); ++i ) {
                    abundances[ transcripts[i] ] = abundance[i];
                }
                ++numGenesProcessed;
            }
        );

        /*
        for( size_t geneID = 0; geneID < _transcriptGeneMap.numGenes(); ++geneID ) {
                //auto abundance = this->_estimateIsoformAbundance(geneID, mappedReads);
                std::cerr << "computing abundance for " << geneID << "  . . . ";
                auto abundance = this->_estimateIsoformNNLS(geneID, mappedReads);
                std::cerr << "done, size = " << abundance.size() << "\n";
                auto transcripts = this->_transcriptGeneMap.transcriptsForGene(geneID);

                for ( size_t i = 0; i < abundance.size(); ++i ) {
                    abundances[ transcripts[i] ] = abundance[i];
                }
                ++numGenesProcessed;
        }
        */
       
        abundanceEstimateProgress.join();

        std::cerr << "Writing output\n";
        ez::ezETAProgressBar pb(_transcripts.size());
        pb.start();

        std::ofstream ofile( outputFile );
        size_t index = 0;
        ofile << "Transcript" << '\t' << "Length" << '\t' << "Effective Length" << '\t' << "Abundance" << '\n';
        for ( auto ts : _transcripts ) {
            ofile << _transcriptGeneMap.transcriptName(index) << '\t' << ts->length << '\t' <<
                    ts->effectiveLength << '\t' << abundances[index] << "\n";
            ++index;
            ++pb;
        }
        ofile.close();

        auto writeCoverageInfo = false;
        if ( writeCoverageInfo ) {
            std::string cfname("transcriptCoverage.txt");
            _dumpCoverage( cfname );
        }

    }


};

#endif // ITERATIVE_OPTIMIZER_HPP
