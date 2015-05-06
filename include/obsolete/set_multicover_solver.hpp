#ifndef SET_MULTICOVER_SOLVER_HPP
#define SET_MULTICOVER_SOLVER_HPP

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <unordered_set>

/** Boost Includes */
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

template <typename ReadHashT, typename TranscriptHashT>
class SetMulticoverSolver {

    private:
        /**
        * Typedefs
        */
        // A KmerMapT is a map from a kmer (encoded as an integer) to the set
        // of transcripts where it occurs
        typedef uint32_t TranscriptIDT;
        typedef uint64_t KmerIDT;
        typedef std::unordered_map< KmerIDT, std::unordered_set<TranscriptIDT> > KmerMapT;
        typedef std::tuple<TranscriptIDT, std::vector<KmerIDT>> TranscriptKmerSet;
        typedef std::string* StringPtrT;
        typedef uint64_t TranscriptScoreT;
        typedef jellyfish::invertible_hash::array<uint64_t, atomic::gcc, allocators::mmap> HashArrayT;
        // Necessary forward declaration
        struct TranscriptDataT;
        typedef std::tuple<TranscriptScoreT,TranscriptIDT> HeapPair;
        typedef typename boost::heap::fibonacci_heap<HeapPair>::handle_type HandleT;

        size_t _merLen;
        ReadHashT& _readHash;
        TranscriptHashT& _transcriptHash;
        std::unordered_map<KmerIDT, size_t> _merHash;

        struct TranscriptDataT {
            uint32_t id;
            uint32_t numTimesUsed;
            HandleT heapHandle;
            StringPtrT header;
            std::unordered_set<uint64_t> binMers;
            bool active;
            double cost;
            double benefit;
        };
    
        // This struct represents a "job" (transcript) that needs to be processed 
        struct TranscriptJob{
        	StringPtrT header;
        	StringPtrT seq;
        	uint32_t id;
        };

        struct TranscriptResult {
        	TranscriptDataT* data;
        	TranscriptKmerSet* ks;
        };

        struct BinmerUpdates {
        	std::vector<KmerIDT> zeroedBinmers;
        	std::vector<KmerIDT> updatedBinmers;
        };

        // The number of occurences above whcih a kmer is considered promiscuous
        size_t _promiscuousKmerCutoff{10};

        // Map each kmer to the set of transcripts it occurs in
        KmerMapT _transcriptsForKmer;

        // The actual data for each transcript
        std::vector<TranscriptDataT*> _transcripts;

        // Is a given kmer considered part of the universe?
        inline bool _considered( const char* mer ) {
        	// The kmer is only considered if it exists in the transcript set 
        	// (i.e. it's possible to cover) and it's less prmiscuous than the
        	// cutoff.
        	return ( _readHash[mer] > 0 and 
        		     _transcriptHash[mer] > 0 and 
        		     _transcriptHash[mer] < _promiscuousKmerCutoff );
        }

        // Process the transcript
        BinmerUpdates _processTranscript( uint32_t transcriptID, uint64_t &rem, size_t& numTimesUsed ) {
        	auto transcript = _transcripts[transcriptID];
        	std::unordered_map<uint64_t, std::unordered_set<uint64_t>> binMerCounts;
        	size_t minCount = std::numeric_limits<size_t>::max();
        	// Both the zeroed kmers and all affected kmers start out empty
        	BinmerUpdates updates{ {}, {} };

            // loop over all kmers (remember, only considered kmers are stored here)
            for ( auto binMer : transcript->binMers ) {
            	// See if this kmer is still active
            	auto merIt = _merHash.find(binMer);
            	// If it is, then
            	if ( merIt != _merHash.end() ) {
            		auto currentCount = merIt->second;
            		// Is it the minimum count kmer?
            		minCount = (currentCount < minCount) ? currentCount : minCount;
            		// Insert it into the map with it's current count
            		binMerCounts[currentCount].insert(binMer);
            	}
            }

            // The number of kmers we're getting rid of
            rem -= binMerCounts[minCount].size();

            // For every kmer in this transcript
            for ( auto countMer : binMerCounts ) {
            	auto origCount = countMer.first;
            	// Get rid of minCount copies of it
            	for ( auto binMer : countMer.second ) {
            		updates.updatedBinmers.push_back(binMer);
            		_merHash[binMer] -= minCount;
            		// If this is a kmer with minimum count, then get rid of it
            		if ( origCount == minCount ) {
            			_merHash.erase(binMer);
            			updates.zeroedBinmers.push_back(binMer);
            		}
            	}
            }
            
            numTimesUsed = minCount;
            return updates;
        }

        void _createUniverse( const std::string& transcriptFile, uint64_t& numRequiredCover ) {

        	// open up the transcript file for reading
            const char *fnames[] = { transcriptFile.c_str() };
            jellyfish::parse_read parser( fnames, fnames + 1, 100);
            jellyfish::parse_read::thread stream = parser.new_thread();

            // Read pointer
            jellyfish::read_parser::read_t *read;

            // So we can concisely identify each transcript
            uint32_t transcriptIndex {0};

            // Holds the tasks to be processed
            boost::lockfree::fifo< TranscriptJob * > workQueue;
            // Holds the processed transcript data
            boost::lockfree::fifo< TranscriptDataT * > processedReadQueue;
            // Holds results that we will use to update the index map
            boost::lockfree::fifo< TranscriptKmerSet * > q;

            size_t numActors = 26;
            boost::threadpool::pool tp(numActors);

            auto readResults = [&q, &tp, this] () -> void {
                size_t numRes = 0;
                // While tasks remain to be processed
                while (tp.active() + tp.pending() > 1) {
                    TranscriptKmerSet *ks;
                    // Process existing results
                    while ( q.dequeue(ks) ) {

                        uint32_t transcriptID = std::get<0>(*ks);
                        // For each kmer in this transcrpt, add this
                        // transcript to all kmers it contains in the map
                        for ( uint64_t kmer : std::get<1>(*ks) ) {
                            this->_transcriptsForKmer[kmer].insert(transcriptID);
                        }

                        // no memory leaks!
                        delete ks;
                        ++numRes;
                        std::cerr << "#res = " << numRes << "\n";
                    }
                    boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
                }
            };

            tp.schedule(readResults);

            bool done {false};
            std::atomic<size_t> numJobs {0};

            for (size_t i = 0; i < numActors - 1; ++i) {

                auto task = [&workQueue, &q, &done, &processedReadQueue, &numJobs, this]() -> void {
                    TranscriptJob *tj;
                    while ( !done || numJobs > 0 ) {
                        while (workQueue.dequeue(tj) ) {
                            auto seq = tj->seq;
                            auto newEnd  = std::remove( seq->begin(), seq->end(), '\n' );
                            auto readLen = std::distance( seq->begin(), newEnd );
                            StringPtrT newSeq = new std::string( seq->begin(), seq->begin() + readLen); 
                            auto transcriptIndex = tj->id;

                            // The set of kmers in this transcript
                            auto ts = new TranscriptKmerSet {transcriptIndex, {}};
                            
                            size_t offset {0};
                            // We know how many kmers there will be, so allocate the space
                            size_t numKmers {readLen - this->_merLen+1};
                            std::get<1>(*ts).reserve(numKmers);
                            // The binary representation of the kmers in this transcript
                            std::unordered_set<uint64_t> binMers;
                            //binMers.reserve(numKmers);

                            // Iterate over the kmers
                            while ( offset < numKmers )  {
                                auto mer = newSeq->substr( offset, this->_merLen );
                                // only count and track kmers which should be considered
                                if ( this->_considered(mer.c_str()) ) {
                                    auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), _merLen );
                                    //binMers.push_back(binMer);
                                    binMers.insert(binMer);
                                    std::get<1>(*ts).push_back(binMer);
                                }

                                ++offset;
                            }

                            // Create the processed read
                            auto procRead = new TranscriptDataT {
                            	transcriptIndex,
                            	0, // Number of times used in the cover
                            	HandleT(0), // The priority heap handle, we don't know it yet
                            	tj->header, // header to idenfity this transcript in the ouptut
                            	binMers, // the set of kmers
                            	false, // not initially active
                            	binMers.size(), // the cost is the number of considered kmers
                            	0 // we don't know the initial benefit
                            };

                            // Push it onto the processed read queue
                            processedReadQueue.enqueue(procRead);
                            // Trnascript set enqueue
                            q.enqueue(ts);

                            // Free what we no longer need
                            delete tj->seq;
                            delete tj;
                            --numJobs;
                        }
                        // if we're not done, but there is no work, then sleep a bit
                        boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
                    }

                };

                tp.schedule(task);
            }


            // Read the input data
            while ( (read = stream.next_read()) ) {

                StringPtrT header = new std::string( read->header, read->hlen );
                StringPtrT oseq = new std::string(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );
                // Create a job to process this read and throw it on the work queue
                auto job = new TranscriptJob { header, oseq, transcriptIndex };
                workQueue.enqueue(job);
                ++transcriptIndex;
                ++numJobs;
            } // end while
            done = true;

            // Wait for all reads to be processed
            tp.wait();

            // Now we know how many transrcipts there will be
            _transcripts.resize(transcriptIndex, nullptr);

            // Pop each transcript off the queue and put it in the right
            // index in the transcript array
            TranscriptDataT* td;
            while ( processedReadQueue.dequeue(td) ) {
            	auto id = td->id;
            	_transcripts[id] = td;
            }

            // iterate over the jellyfish hash and count the number of occurences of
            // each kmer that appears in the set of reads
            auto it = _transcriptHash.iterator_all();
            while ( it.next() ) {
            	auto mer = it.get_dna_str();
            	auto binMer = jellyfish::parse_dna::mer_string_to_binary(mer, _merLen);
            	auto count = _readHash[binMer];//it.get_val();
            	// we only care about considered kmers
            	if ( _considered(mer) ) {
            		_merHash[binMer] = count;//.map( binMer, count );
            		numRequiredCover += 1;
            	}
            }
        }

    public:
    	/**
    	* Construct the solver with the read and transcript hashes
    	*/
    	SetMulticoverSolver( ReadHashT& readHash, TranscriptHashT&  transcriptHash ) : 
    	_readHash(readHash), _transcriptHash(transcriptHash), _merLen(readHash.get_mer_len()),
    	_merHash( std::unordered_map<KmerIDT,size_t>() ) {}
    		/*
    	_merHash( HashArrayT(transcriptHash.get_size(), // Size of hash. Will be rounded up to a power of 2
                     2*_merLen,      // Key length in bits (2 * mer length)
                     30,      // Value length in bits
                     126,     // Max # of reprobe
                     jellyfish::quadratic_reprobes // Reprobe policy
                     ) ) {}
	    */

        bool operator()(const std::string& transcriptFile, const std::string& outputFile ) {

            typedef size_t TranscriptCountT;

            // Hold the number of kmers that remain to be processed
            uint64_t numRemainingKmers = 0;

            // Create the universe and initialize all necesssary data
            // structures for solving the cover problem
            _createUniverse( transcriptFile, numRemainingKmers );

            boost::heap::fibonacci_heap<HeapPair> pq;

            // Function to compute the cost of a given transcript.  Right now it's
            // 1 (all transcripts have uniform cost), but this could be changed in
            // the future.
            auto cost = []( TranscriptDataT* transcript ) -> double { return 1.0; };

            // The reward of a given transcript.  For each required kmer that a transcript
            // covers, it's score is higher.
            auto reward = [this]( TranscriptDataT* transcript ) -> TranscriptScoreT {
                TranscriptScoreT benefit{0};
                // Iterate over the considered kmers in the transcript
                for ( auto binMer : transcript->binMers ) {
                	auto currentIt = _merHash.find(binMer);
                	benefit += ( currentIt != _merHash.end() ) ? 1 : 0;//currentIt->second : 0; 
                }
                return benefit;
            };

            //std::atomic<size_t> transcriptsToProcess = _transcripts.size();
            // Score each transcript and add it to the queue
            for ( auto tdata : _transcripts ) { 
            	auto benefit = reward(tdata);
                //auto c = cost(tdata);
                tdata->benefit = benefit;
                if ( tdata->benefit > 0 && tdata->cost > 0 ) {
                	tdata->heapHandle = pq.push( std::make_tuple( tdata->benefit / tdata->cost , tdata->id) );
                	tdata->active = true;
                }
            }

            //tp.wait();
            std::cerr << "numRemainingKmers = " << numRemainingKmers << "\n";
            // While not everything has been covered the required number of times
            boost::progress_display show_progress( numRemainingKmers );
            while ( numRemainingKmers > 0 and !pq.empty() ) {
                // The score of this best transcript
                double score;
                // The best transcript's ID
                TranscriptIDT t;

                // Get the best transcript to use and its score
                std::tie(score, t) = pq.top();
                // Update the cover --- add a copy of this transcript 
                // and update the number of things remaining to be covered
                auto preSize = numRemainingKmers;
                // Number of times this transcript is used in this round
                size_t numCopies{0};
                // Gets the set of kmers whose remaining count became zeros during
                // this round
                auto binmerUpdates = _processTranscript(t, numRemainingKmers, numCopies);
                /*
                std::cerr << "best transcript is " << t << "\n";
                std::cerr << "before, used " << _transcripts[t]->numTimesUsed << " times; numCopies = " << numCopies;
                */
                _transcripts[t]->numTimesUsed += numCopies;
                /*
                std::cerr << " after, used " << _transcripts[t]->numTimesUsed << " times\n";
                std::cerr << "got rid of " << binmerUpdates.zeroedBinmers.size() << " kmers\n";
                std::cerr << "with total weight " << binmerUpdates.updatedBinmers.size() * numCopies << "\n";
                */
                if( binmerUpdates.zeroedBinmers.size() == 0 ) { std::cerr << "zeroedBinmers == 0\n"; std::abort(); }
                show_progress += (preSize - numRemainingKmers);

                // Update all of the transcripts containing a changed kmer (including
                // the guy at the top of the heap)
                std::unordered_map<TranscriptDataT*, uint64_t> affectedTranscripts;
                for ( auto binMer : binmerUpdates.updatedBinmers ) {
                    // For every transcript containing this kmer, its scorew ill be decreased
                    // by numCopies kmers
                    for ( auto transcriptID : _transcriptsForKmer[binMer] ){
                    	affectedTranscripts[ _transcripts[transcriptID] ] += numCopies;
                    }
                }
                //std::cerr << "\n";
                // update the affected kmers
                for ( auto transcriptDecrease : affectedTranscripts ) {
                	
                	auto transcript = transcriptDecrease.first;
                	//std::cerr << "Updating transcript " << transcript->id << "\n";
                	// Remove these kmers from the transcript
                	
                	size_t numErased = 0;
                	/*
                	auto br = reward(transcript);
                	std::cerr << "it contains " << transcriptDecrease.second.size() << " affected kmers; [";
                	for ( auto binMer : transcriptDecrease.second ) { std::cerr << " " << binMer; }
                	std::cerr << " ]\n";
	                */

                	for ( auto binMer : binmerUpdates.zeroedBinmers ) {
                	    //for ( auto binMer : changedKmers ) { 
                		/*
                		auto it = _merHash.find(binMer);
                		if( transcript->binMers.find(binMer) == transcript->binMers.end() ) {
                			std::cerr << "Transcript should have contained " << binMer << " but it didn't!\n";
                			std::abort();
                		}
                		std::cerr << "erasing " << binMer << " from " << transcript->id << "\n";
                		if ( it != _merHash.end() ) {
                	    */
              			auto erased = transcript->binMers.erase(binMer); 
               			numErased += erased;
                	}

                	if ( transcript->active ) {
                	transcript->benefit -= numErased;
                	std::get<0>( *transcript->heapHandle ) = transcript->benefit / transcript->cost;

                	 /*              	
                	auto decrease = transcriptDecrease.second;
                	auto benefit = std::get<0>( *transcript->heapHandle );
                	std::get<0>( *transcript->heapHandle ) -= decrease;
                	*/

                	/*
                	if (  std::get<0>( *transcript->heapHandle ) != reward(transcript) ) {
                		std::cerr << "original = " << benefit << ", decrease = " << decrease << "\n";
                		std::cerr << "calculated reward explicitlly as " << reward(transcript) <<
                		"; a decrease of " << benefit-reward(transcript) << 
                		" but updated it to" << std::get<0>( *transcript->heapHandle ) << "\n";
                		std::abort();
                	}
                	*/
                	if ( transcript->benefit == 0 ) {
                		if ( transcript->active ) {
                			pq.erase( transcript->heapHandle ); 
                		}
                		transcript->active = false;
                	} else {
                		pq.update( transcript->heapHandle );
                	}
                    }
                	//std::cerr << "zeroed: " << binmerUpdates.zeroedBinmers.size() << ", affected: " << binmerUpdates.updatedBinmers.size() << "\n";

                	    /*
                        auto benefit = reward(transcript);
                        auto c = cost(transcript);
                        // Compute the new score and update it in the heap
                        // If the reward is zero, then drop this guy
                        if ( benefit < 1e-5 ) { 
                        	pq.erase( transcript->heapHandle ); 
                        } else {
                   	         // Otherwise, update him
                        	std::get<0>( *transcript->heapHandle ) = benefit;
                        	pq.update( transcript->heapHandle );
                        }
                        */
                } // done updating affected transcripts

            } // done solving the cover problem
            
            std::cerr << "numRemainingKmers: " << numRemainingKmers << "\n";

            std::ofstream ofile( outputFile );
            for ( auto ts : _transcripts ) {
            	ofile << *(ts->header) << '\t' << ts->numTimesUsed << "\n";
            }
            ofile.close();
        }

};

#endif // MULTI_SETCOVER_SOLVER_HPP