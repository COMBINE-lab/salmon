#ifndef TRANSCRIPTSEGMENTER_HPP
#define TRANSCRIPTSEGMENTER_HPP

#include <cassert>
#include <cstdint>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <sstream>
#include <limits>
#include <exception>
#include <ext/vstring.h>
#include <ext/vstring_fwd.h>

#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/progress.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/range/irange.hpp>

#include <Eigen/Core>
//#include <Eigen/Sparse>

#include "g2log.h"

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#include "boost/lockfree/fifo.hpp"
#include "threadpool.hpp"
#include "graph_utils.hpp"
#include "transcript_gene_map.hpp"
#include "poisson_solver.hpp"

//#include "tsnnls/tsnnls.h"
//#include "shotgun/common.h"
//extern "C"
//{
#include "nnls.h"
//}

#include "matrix_tools.hpp"

template <typename HashT>
class TranscriptSegmenter {

private:
    // typedefs to make things less noisy later on
    //typedef __gnu_cxx::__vstring string;
    typedef std::string string;
    typedef string *StringPtrT;
    typedef std::mutex mutex;
    typedef tbb::concurrent_unordered_map< uint64_t, tbb::concurrent_unordered_set<uint32_t> > KmerMapT;

    // The final result type of our analysis; mapping each transcript to an abundance
    struct ResultT {
        StringPtrT name;
        double abundance;
    };

    // This struct represents a "job" (transcript) that needs to be processed
    struct TranscriptJob {
        StringPtrT header;
        StringPtrT seq;
        uint32_t index;
    };


    size_t errCtr {0};

    // If a connected subgraph contains more than this many transcripts, then
    // it will be partitioned into smaller subgraphs
    const size_t tooManyVertices {
        1500
    };

    // Filename where we'll read the transcripts from
    std::string _transcriptFile;

    // Map from each kmer to the transcripts it occurs in
    KmerMapT _indexMap;
    KmerMapT _kmerGeneMap;

    // Jellyfish hash of hte reads
    HashT _readHash;
    //HashT _transcriptHash;

    TranscriptGeneMap &_transcriptGeneMap;

    /**
    * Private member functions
    **/
    bool _considered( const char *mer ) {
        return true;
    }
    /**
    * Given a transcript graph G, and the vector of transcripts, partition G into
    */
    template <typename GraphT>
    void _partitionIntoSubgraphs( GraphT &G,
                                  std::vector< TranscriptJob * > &transcripts,
                                  std::unordered_map<size_t, GraphT> &subgraphMap) {

        LOG(DEBUG) << "in _partitionIntoSubgraphs\n";

        // Create a metis graph from our boost graph
        graphdef mG;
        int wgtflag {1}; // we have edge weights

        LOG(DEBUG) << "calling boostToMetis\n";

        utils::graph::boostToMetis( G, &mG );

        LOG(DEBUG) << "done boost to metis\n";
        LOG(DEBUG) << "num verts = " << mG.nvtxs << ", num edges = " << mG.nedges << "\n";

        // Hold the partitions
        idxtype *part = new idxtype[mG.nvtxs];

        LOG(DEBUG) << "alloctated new partition vector\n";

        // Initialize the clustering options
        Options opt;
        initOptions(&opt);
        opt.matchType = MATCH_SHEMN;

        LOG(DEBUG) << "calling mlmcl\n";
        LOG(DEBUG) << "num verts = " << mG.nvtxs << ", num edges = " << mG.nedges << "\n";
        LOG(DEBUG) << "num edges 2 (should = num edges) = " << mG.xadj[mG.nvtxs] << "\n";

        // Call the clustering procedure
        mlmcl(&mG.nvtxs, mG.xadj, mG.adjncy, mG.vwgt,
              mG.adjwgt, &wgtflag, part, opt );

        LOG(DEBUG) << "done mlmcl\n";

        size_t numClust {0};
        // Create the subgraphs
        for ( size_t i = 0; i < mG.nvtxs; ++i ) {
            auto cluster = part[i];
            if ( subgraphMap.find(cluster) == subgraphMap.end() ) {
                subgraphMap[cluster] = G.create_subgraph();
            }
        }
        // For each vertex, add it to the appropriate subgraph
        for ( size_t i = 0; i < mG.nvtxs; ++i ) {
            auto cluster = part[i];
            boost::add_vertex(i, subgraphMap[cluster]);
        }

    }

    /**
    * Given the matrix A (m x n) and the right-hand-side vector b (m x 1) solve
    * for the solution vector x ( n x 1 ) satisfying
    * min_{x} || Ax - b ||
    * subject to x[i] >= 0 for all i
    */
    template< typename MatT, typename VecT >
    std::vector<double> _solveNonNegativeLeastSquares( MatT &A, VecT &b ) {
        std::cerr << "in least squares problem\n";

        Eigen::MatrixXd ADense(A);


        int    mda = A.rows();
        int    m = A.rows(), n = A.cols();
        std::vector<double> x( A.cols(), 1.0 );
        double rnorm;
        std::vector<double> w(n, 0.0);
        std::vector<double> zz(m, 0.0);
        std::vector<int> indx(n, 0);
        int mode;

        nnls( ADense.data(), mda, m, n, &b[0], &x[0], &rnorm, &w[0], &zz[0], &indx[0], &mode);
        //auto adata = ADense.data();
        //nnls( &adata, m, n, &b[0], &x[0], &rnorm, &w[0], &zz[0], &indx[0]);
        std::vector<double> myX( A.cols(), 0.0 );

        for ( size_t i = 0; i < myX.size(); ++i ) {
            myX[i] = x[i];
            std::cerr << "x[" << i << "]  = " << x[i] << "\n";
        }

        return myX;

        /*
        // Create a new TAUCS matrix structure
        taucs_ccs_matrix* mat;
        mat = new taucs_ccs_matrix;//(taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));

        mat->n = A.cols();
        mat->m = A.rows();
        mat->flags = TAUCS_DOUBLE;

        std::cerr << "inited TAUCS matrix\n";

        // Compress this matrix so we can steal / share the data
        A.makeCompressed();

        std::cerr << "compressed matrix\n";

        mat->colptr = A.outerIndexPtr();
        mat->rowind = A.innerIndexPtr();
        mat->values.d = A.valuePtr();

        double residualNorm{0.0};

        std::cerr << "called NNLS solver\n";
        auto x = t_snnls_fallback( mat, &b[0], &residualNorm, 0.0, 0);
        char *errString;
        tsnnls_error(&errString);
        std::cerr << "NNLS solver returned\n";
        std::cerr << "ERROR STRING " << errString << "\n";

        if( x == NULL ) {
          std::cerr << "AHHH, x is STILL NULL\n";

          std::stringstream ss;
          ss << "errorA_" << errCtr << ".mtx";
          std::ofstream Afile(ss.str());

          Afile << "SPARSE\n";
          Afile << A.rows() << '\t' << A.cols() << '\n' << A.nonZeros() << '\n';
          for (int k=0; k < A.outerSize(); ++k) {
            for (typename MatT::InnerIterator it(A,k); it; ++it) {
              Afile << it.row()+1 << '\t' << it.col()+1 << '\t' << it.value() << '\n';
            }
          }
          Afile.close();

          ss.clear();
          ss << "errorb_"<< errCtr << ".mat";
          std::ofstream bfile("errorb.mtx");
          bfile << "# error right hand side\n";
          bfile << "# rows: " << b.size() << "\n";
          bfile << "# columns: 1\n";
          for (size_t i = 0; i < b.size(); ++i) {
              bfile << b[i] << '\n';
          }
          bfile.close();
          ++errCtr;

          delete mat;
          return std::vector<double>();
          //std::abort();
        } else {

          std::vector<double> myX( A.cols(), 0.0 );

          for( size_t i = 0; i < myX.size(); ++i ) { myX[i] = x[i]; std::cerr << "x[" << i << "]  = " << x[i] << "\n"; }

          free(x);
          delete mat;

          return myX;
        }
        */

        /*
        // Try shotgun LASSO solver
        shotgun_data prob;
        std::cerr << "trying to reserve " << A.nonZeros() << " nonzeros in prob rows/cols\n";
        prob.A_cols.resize( A.nonZeros(), sparse_array() );
        prob.A_rows.resize( A.nonZeros(), sparse_array() );

        std::cerr << "done\n";

        auto nnz = A.nonZeros();
        size_t ctr{0};
        for (int k=0; k < A.outerSize(); ++k) {
          std::cerr << "column " << k << " of " << A.outerSize() << "\n";
          for (typename MatT::InnerIterator it(A,k); it; ++it) {
            auto val = it.value();
            auto i = it.row();   // row index
            auto j = it.col();   // col index (here it is equal to k)
            //std::cerr << "adding " << i << ", " << j << " : " << val << "\n";
            try {
            prob.A_cols[j].add(i, val);
            prob.A_rows[i].add(j, val);
            ++ctr;
            if ( ctr > nnz ) { std::cerr << "only reserved space for " << nnz <<
              " elements, but you're pushing on the " << ctr << "th element\n";
            }
            } catch ( std::exception& e ) {
              std::cerr << "died trying to add " << i << ", " << j << " : " << val << "\n";
              std::cerr << e.what();
            }
          }
        }

        prob.nx = A.cols();
        prob.ny = A.rows();

        prob.y.reserve(b.size());
        for ( auto e : b ) { prob.y.push_back(e); };
        double lambda = 0.0;
        int K = 10;
        int maxiter = 5000;
        int verbose = 10;
        double threshold = 1e-8;
        assert( prob.ny == prob.y.size() );
        std::cerr << "before solve lasso\n";
        solveLasso(&prob, lambda, K, threshold, maxiter, verbose);
        std::cerr << "after solve lasso\n";

        std::vector<double> myX(prob.x);
        return myX;
        */

    }



    template <typename GraphT, typename TranscriptHashT>
    std::vector<ResultT> _estimateViaNNLS( GraphT &G, TranscriptHashT &transcriptHash, std::vector< TranscriptJob * > &transcriptJobs ) {

        typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
        typedef Eigen::Triplet<double> T;

        struct TranscriptColumn {
            size_t localID;
            StringPtrT name;
            std::vector<uint64_t> kmers;
        };

        struct LocalTranscriptJob {
            size_t localID;
            TranscriptJob *globalJob;
        };

        auto numJobs = boost::num_vertices(G);
        size_t numActors = (26 > numJobs + 1) ? numJobs + 1 : 26; //min(26, boost::num_vertices(G)+1);
        boost::threadpool::pool tp(numActors);

        // Holds the transcripts to be processed
        boost::lockfree::fifo< LocalTranscriptJob * > workQueue;

        // Holds the results of processed transcripts
        boost::lockfree::fifo< TranscriptColumn * > resultQueue;

        // The matrix that will hold the left-hand side of our problem
        SpMat *A = nullptr;
        // The vector that will hold the right-hand side of our problem
        std::vector<double> b;

        auto merLen = _readHash.get_mer_len();

        // Collect the processed transcripts and build a matrix out of them
        auto collectResults = [&resultQueue, &A, &b, &tp, &transcriptHash, this, merLen, numJobs]() -> void {
            // Char array to hold a kmer string
            char merStr[merLen];
            // Triplets we'll use to make our matrix
            std::vector<T> triplets;
            // Map from a kmer (encoded in binary) to its row ID
            auto kmerSet = std::unordered_map<uint64_t, uint32_t>();

            // The number of rows and columns in our matrix
            size_t numRows{0};
            size_t numCols{0};
            size_t resultsProcessed{0};

            // While tasks remain to be processed
            TranscriptColumn *tc;
            while ( resultsProcessed < numJobs ) {
                // Pull a result off the queue
                while ( resultQueue.dequeue(tc) ) {
                    std::cerr << "unqueueing result " << tc->name << "\n";
                    // Get the column id
                    auto columnID = tc->localID;
                    // Update the max column id
                    numCols = ( numCols > (columnID + 1)) ? numCols : columnID + 1 ;

                    // Map from a row ID to it's value (likely just a "1")
                    std::unordered_map<uint32_t, float> rowIDs;
                    // The number of non-zeros in this column
                    rowIDs.reserve(tc->kmers.size());

                    // For each kmer, map it to a local id
                    for ( auto & kmer : tc->kmers ) {
                        auto weight = ( 1.0 / (0.1 + transcriptHash[merStr]));
                        // If this is a new kmer, then give it the next available id
                        if ( kmerSet.find(kmer) == kmerSet.end() ) {
                            // Update the number of rows
                            kmerSet[kmer] = numRows;
                            ++numRows;

                            // Put the count of this kmer in the right-hand side
                            jellyfish::parse_dna::mer_binary_to_string( kmer, merLen, &merStr[0] );
                            b.push_back( this->_readHash[merStr] * weight );
                        };
                        // Assign the local id to this column
                        rowIDs[ kmerSet[kmer] ] += weight;
                    }

                    // Push this column's values onto our triplets array
                    for (auto & kv : rowIDs ) {
                        triplets.push_back( {kv.first, columnID, kv.second} );
                    }
                    ++resultsProcessed;
                }

                // if we're not done, but there is no work, then sleep a bit
                boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
            }

            // Build a proper matrix from our triplets
            /*size_t maxRow{0};
            int minRow = 2000000;
            for ( auto& t : triplets ) { std::cerr << "pushing " << t.row() << ", " << t.col() << ", " << t.value() << "\n";
            maxRow = max(maxRow,t.row());
            minRow = min(minRow,t.row());
            }
            */
            std::cerr << "numCols = " << numCols << "\n";
            std::cerr << "kmerSet.size = " << kmerSet.size() << ", numRows = " << numRows << "\n";
            /*std::cerr << "maxRow = " << maxRow << "\n";
            std::cerr << "minRow = " << minRow << "\n";*/
            if ( numRows > 0 && numCols > 0 ) {
                A = new SpMat(numRows, numCols);
                A->setFromTriplets(triplets.begin(), triplets.end());
            }
            std::cerr << "done building matrix\n";
        };

        tp.schedule(collectResults);

        // Read Hash
        bool done {false};

        std::atomic<size_t> jobsRemaining {numJobs};
        // Launch the worker jobs
        for ( size_t i = 0; i < numActors - 1; ++i ) {
            //
            auto columnBuilder = [&workQueue, &resultQueue, &done, &jobsRemaining, &transcriptHash, merLen, this]() -> void {
                LocalTranscriptJob *ltj;
                while ( jobsRemaining > 0 ) {
                    while (workQueue.dequeue(ltj) ) {

                        TranscriptJob *tj = ltj->globalJob;
                        StringPtrT seq = tj->seq;
                        auto readLen = seq->length();
                        auto transcriptIndex = tj->index;
                        auto localIndex = ltj->localID;

                        size_t offset {0};
                        auto transcriptColumn = new TranscriptColumn {
                            localIndex,
                            tj->header,
                            std::vector<uint64_t>()
                        };

                        if ( transcriptColumn->name->length() != tj->header->length() ) {
                            std::cerr << "transcriptColumn name is " << *transcriptColumn->name << ", but "  <<
                                      "tj->header is " << *tj->header << "\n";
                        }

                        std::cerr << "trying to allocate kmers vector of size " << readLen - merLen + 1 << "\n";
                        transcriptColumn->kmers.reserve( readLen - merLen + 1 );

                        // Iterate over the kmers
                        while ( offset < readLen - merLen + 1 )  {
                            auto mer = seq->substr( offset, merLen );
                            if ( this->_considered(mer.c_str()) ) {// transcriptHash[mer.c_str()] < 10 ) {
                                auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
                                transcriptColumn->kmers.push_back(binMer);
                            }
                            ++offset;
                        }

                        // Enqueue the column
                        resultQueue.enqueue(transcriptColumn);
                        std::cerr << "processed job " << *transcriptColumn->name << " enqueueing result\n";
                        --jobsRemaining;
                        std::cerr << "there are " << jobsRemaining << " jobs remaining\n";
                        // Delete the original transcript job
                        // delete tj;
                        // Delete the local job as well
                        // delete ltj;

                    }
                    // if we're not done, but there is no work, then sleep a bit
                    boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
                }

            };

            tp.schedule(columnBuilder);
        } // end launch of workers

        std::vector<StringPtrT> columnNames( boost::num_vertices(G), nullptr );
        // for each vertex in the graph
        auto se = boost::vertices(G);
        for ( auto &v = se.first; v != se.second; ++v ) {
            // Get the global id of the vertex
            auto globalIdx = G.local_to_global( *v );
            // Create a local job from the global job
            auto localJob = new LocalTranscriptJob { *v, transcriptJobs[globalIdx] };
            // Add this transcripts name to the index
            columnNames[*v] = transcriptJobs[globalIdx]->header;//std::string(header.begin(), he);
            std::cerr << "pusing job " << *(transcriptJobs[globalIdx]->header) << ", for vertex " << *v << "\n";
            // build column for this vertex
            workQueue.enqueue(localJob);
        } // end launch of jobs
        done =  true;
        tp.wait();

        std::cerr << "done waiting\n";
        std::vector<ResultT> results;

        if ( A != nullptr ) {
            auto dupCols = matrix_tools::markDependentColumns(*A);

            if (dupCols.size() > 0) {
                std::cerr << "Matrix had " << dupCols.size() << " duplicate columns\n";
                std::vector<T> triplets;
                std::vector<StringPtrT> newColumnNames;

                size_t dupPtr = 0;
                for (size_t k = 0; k < A->outerSize(); ++k) {
                    // If this is a colum we skip, then continue this loop and
                    // move to the next column to skip
                    if ( dupPtr < dupCols.size() && k == dupCols[dupPtr] ) {
                        ++dupPtr;
                        continue;
                    }
                    // The column index is the column index of A minus the number of
                    // columns we've skipped so far
                    size_t colInd = k - dupPtr;
                    // Add all rows of this column to the new matrix
                    for (SpMat::InnerIterator it(*A, k); it; ++it) {
                        triplets.push_back( {it.row(), colInd, it.value()} );
                    }
                    newColumnNames.push_back(columnNames[colInd]);
                }
                SpMat *B = new SpMat(A->rows(), A->cols() - dupCols.size());
                B->setFromTriplets(triplets.begin(), triplets.end());
                std::swap(A, B);
                std::swap(columnNames, newColumnNames);
                delete B;
            }




            // At this point we've built the matrix, so solve the problem
            auto x = _solveNonNegativeLeastSquares(*A, b);
            results.reserve(x.size());
            for ( size_t i = 0; i < x.size(); ++i ) {
                results.push_back( {columnNames[i], x[i]} );
            }
            delete A;
        }

        return results;
    }

    template <typename GraphT, typename TranscriptHashT>
    std::vector<ResultT> _computeAbundance( GraphT &G, TranscriptHashT &transcriptHash, std::vector< TranscriptJob * > transcriptJobs ) {
        auto numVerts = boost::num_vertices(G);
        // If the graph is non-trivial
        //if ( numVerts > 1 ) {

        // If the graph is "too big", then we cluster it
        if ( numVerts > tooManyVertices ) {
            std::unordered_map<size_t, GraphT> subgraphMap;
            _partitionIntoSubgraphs(G, transcriptJobs, subgraphMap);
            // Now G has some set of "children" subgraphs
            auto se = G.children();

            // Iterate over each child subgraph and build it's linear system
            std::vector<ResultT> res;
            for ( auto & kv : subgraphMap ) {
                // recursively decompose
                auto sgres = _computeAbundance(kv.second, transcriptHash, transcriptJobs);//_estimateViaNNLS(kv.second, transcriptJobs);
                res.insert(res.end(), sgres.begin(), sgres.end());
            }
            return res;

        } else {
            std::cerr << "Component has only " << boost::num_vertices(G) << " vertices; solving NNLS directly\n";
            // Otherwise, solve the nnls problem directly
            auto res = _estimateViaNNLS(G, transcriptHash, transcriptJobs);
            return res;
        }
        /*
        } else {
        // The graph has only one vertex; its kmer set doesn't overlap with
        // any others and so its abundance can be uniquely estimated
        auto singletonVertex = boost::vertices(G).first;
        auto globalIdx = G.local_to_global( *singletonVertex );
        //auto res = _estimateAbundanceTrivial( globalIdx );
        return std::vector<;
        }
        */
    }

    template <typename hash_t>
    std::vector<double> _estimateIsoformAbundanceByCount( const size_t geneID, 
      std::vector<std::vector<uint64_t>> &transcriptKmers,
      std::vector<double>& transcriptLengths,
      size_t mappedReads,
      hash_t& _transcriptHash
      ) {

        auto transcripts = _transcriptGeneMap.transcriptsForGene( geneID );
        size_t numTranscripts = transcripts.size();
        std::vector<double> exps(numTranscripts, 0.0);

        for ( auto ti : boost::irange(0ul, numTranscripts) ){
            auto globalTranscriptID = transcripts[ti];

            for ( auto kmer : transcriptKmers[ globalTranscriptID ] ) {
                exps[ti] += _readHash[kmer] / (1.0 + _transcriptHash[kmer]);
            }
        }

        return exps;
    }

    std::vector<double> _estimateIsoformAbundance( const size_t geneID, 
      std::vector<std::vector<uint64_t>> &transcriptKmers,
      std::vector<double>& transcriptLengths,
      size_t mappedReads,
      tbb::concurrent_unordered_map<uint64_t, uint32_t>& genePromiscuousKmers
      ) {

        typedef std::vector<std::vector<double>> IsoOptMatrix;
        typedef std::vector<double> IsoOptVector;

        auto transcripts = _transcriptGeneMap.transcriptsForGene( geneID );
        size_t numTranscripts = transcripts.size();

        std::unordered_map<uint64_t, size_t> kmerIDs;
        std::unordered_set<uint64_t> geneKmers;
        // reindex the kmers for this gene
        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kmer : transcriptKmers[ transcripts[i] ] ) {                
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
        IsoOptMatrix A( numTranscripts, IsoOptVector(numKmers, 0.0) );
        IsoOptVector counts(numKmers, 0.0);
        double numReads = 0.0;

        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kmer : transcriptKmers[ transcripts[i] ] ) {
                size_t j = kmerIDs[kmer];
                A[i][j] = 1.0;
                if ( counts[j] == 0 ) {

                    auto promIt = genePromiscuousKmers.find(kmer);
                    double normFact = 1.0;
                    if ( promIt != genePromiscuousKmers.end() ) {
                        normFact /= 4.0*promIt->second;
                    }

                    counts[j] = _readHash[kmer] * normFact;
                    numReads += counts[j];
                }
            }
        }

        
        std::unique_ptr<IsoOptMatrix> collapsedA;
        std::unique_ptr<IsoOptVector> collapsedCounts;
        matrix_tools::collapseIntoCategories(A, counts, collapsedA, collapsedCounts);
        
        for (size_t i = 0; i < collapsedA->size(); ++i) {
            auto tlen = transcriptLengths[ transcripts[i] ];
          for (size_t j = 0; j < (*collapsedA)[0].size(); ++j) {
            (*collapsedA)[i][j] = (*collapsedA)[i][j] / ( tlen * mappedReads / 1000000000.0 );
            //(*collapsedA)[i][j] *= numReads;
          }
        }

        //transcriptLengths[transcripts[i]] / 1000.0 * mappedReads / 1000000.0;

        solve_likelihood( *collapsedA, *collapsedCounts, x, false );

        //solve_likelihood( A, counts, x, false );
        //for ( auto& e : x ) { e *= numReads; }
        return x;
    }


public:
    TranscriptSegmenter( const std::string &transcriptFile, HashT &readHash, TranscriptGeneMap &t2g ) :
        _transcriptFile (transcriptFile), _indexMap( KmerMapT() ), _readHash(readHash), _transcriptGeneMap(t2g) {}

    template <typename TranscriptHashT>
    bool findOverlappingSegments( int merLen, TranscriptHashT &transcriptHash ) {
        // We'll need the reverse map (gene -> [transcripts])
        _transcriptGeneMap.needReverse();

        // Holds the set of kmers present in each transcript
        typedef std::tuple<uint32_t, std::vector<uint64_t>> TranscriptKmerSet;

        // Parser for the reads
        const char *fnames[] = { _transcriptFile.c_str() };
        jellyfish::parse_read parser( fnames, fnames + 1, 100);
        jellyfish::parse_read::thread stream = parser.new_thread();

        // Read pointer
        jellyfish::read_parser::read_t *read;

        // So we can concisely identify each transcript
        uint32_t transcriptIndex {0};

        // Holds the tasks to be processed
        boost::lockfree::fifo< TranscriptJob * > workQueue;

        // Holds the processed transcripts (removed newlines and index of each)
        boost::lockfree::fifo< TranscriptJob * > processedReadQueue;

        // Holds results that we will use to update the index map
        //boost::lockfree::fifo< TranscriptKmerSet * > q;

        std::atomic<size_t> numJobsProcessed {0};
        std::atomic<size_t> numRes {0};
        size_t numActors = 26 + 1;
        boost::threadpool::pool tp(numActors);

        std::vector<double> transcriptLengths( _transcriptGeneMap.numTranscripts(), 0.0 );

        auto printProgressBar = [&numRes, this] () {
            auto numJobs = this->_transcriptGeneMap.numTranscripts();
            boost::progress_display show_progress( numJobs );
            size_t lastCount = numRes;
            while ( lastCount < numJobs ) {
                auto diff = numRes - lastCount;
                if ( diff > 0 ) {
                    show_progress += diff;
                    lastCount += diff;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
            }
        };
        std::cerr << "Processing transcripts\n";
        tp.schedule(printProgressBar);


        /**
        * Read results from the result queue, and build the index map
        */
        /*
        auto readResults = [&q, &tp, &numRes, this] () -> void {
            // While tasks remain to be processed
            while (tp.active() + tp.pending() > 1) {
                TranscriptKmerSet *ks;
                // Process existing results
                while ( q.dequeue(ks) ) {

                    uint32_t transcriptID = std::get<0>(*ks);
                    // For each kmer in this transcrpt, add this
                    // transcript to all kmers it contains in the map
                    for ( uint64_t & kmer : std::get<1>(*ks) ) {
                        //this->_indexMap[kmer].insert(transcriptID);
                        this->_kmerGeneMap[kmer].insert( this->_transcriptGeneMap.gene(transcriptID) );
                    }

                    // no memory leaks!
                    //delete ks;
                    ++numRes;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
            }
        };

        tp.schedule(readResults);
        */

        // Holds, for each transcript, the kmers present in it
        std::vector< std::vector<uint64_t> > transcriptKmers(_transcriptGeneMap.numTranscripts(), std::vector<uint64_t>());

        /**
        * Spawn off the appropriate number of worker threads.  Each will process reads from the
        * transcript file.
        */
        for (size_t i = 0; i < numActors; ++i) {

            auto task = [&workQueue, &transcriptHash, &processedReadQueue, &numJobsProcessed,
            &transcriptKmers, &numRes, &transcriptLengths, merLen, this]() -> void {
                TranscriptJob *tj;
                size_t numTranscripts = this->_transcriptGeneMap.numTranscripts();

                while ( numJobsProcessed < numTranscripts ) {
                    while (workQueue.dequeue(tj) ) {
                        // Get the sequence
                        auto seq = tj->seq;
                        // Strip the newlines from it
                        auto newEnd  = std::remove( seq->begin(), seq->end(), '\n' );
                        auto readLen = std::distance( seq->begin(), newEnd );
                        // Copy the newline-stripped read into newSeq
                        StringPtrT newSeq = new string( seq->begin(), seq->begin() + readLen);
                        // We no longer need the original sequence
                        delete seq;

                        // Lookup the ID of this transcript in our transcript -> gene map
                        auto transcriptIndex = this->_transcriptGeneMap.findTranscriptID( *(tj->header) );
                        assert( transcriptIndex != this->_transcriptGeneMap.INVALID );

                        /*
                        // Create the processed read
                        auto procRead = new TranscriptJob{
                          tj->header,
                          newSeq,
                          transcriptIndex
                        };
                        // Push it onto the processed read queue
                        processedReadQueue.enqueue(procRead);
                        */

                        // The set of kmers in this transcript
                        //auto ts = new TranscriptKmerSet{transcriptIndex, {}};

                        // We know how many kmers there will be, so allocate the space
                        //std::get<1>(*ts).reserve(readLen-merLen+1);
                        transcriptKmers[transcriptIndex].reserve(readLen - merLen + 1);
                        size_t offset {0};

                        // Iterate over the kmers and add them to the kmer list for this transcript
                        while ( offset < readLen - merLen + 1 )  {
                            auto mer = newSeq->substr( offset, merLen );
                            auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
                            if ( this->_considered(mer.c_str()) ) {                                
                                this->_kmerGeneMap[binMer].insert( this->_transcriptGeneMap.gene(transcriptIndex) );
                                transcriptKmers[transcriptIndex].push_back(binMer);//std::get<1>(*ts).push_back(binMer);
                            }
                            ++offset;
                        }
                        transcriptLengths[transcriptIndex] = readLen - merLen + 1;

                        ++numJobsProcessed;
                        ++numRes;
                        delete tj;
                        delete newSeq;
                        //q.enqueue(ts);
                        //delete tj;
                    }
                    // if we're not done, but there is no work, then sleep a bit
                    boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
                }

            };

            tp.schedule(task);
        }


        // Read the input data from the trnascript file
        while ( (read = stream.next_read()) ) {
            StringPtrT header = new string( read->header, read->hlen );
            StringPtrT oseq = new string(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );
            auto job = new TranscriptJob { header, oseq, transcriptIndex };
            workQueue.enqueue(job);
        }// end while

        tp.wait();
        std::cerr << "done waiting\n";

        /*
        // Now that we have all of the processed reads, along with their indices
        // in the processed read queue; put them in a vector
        std::vector< TranscriptJob* > transcripts( _transcriptGeneMap.numTranscripts() , nullptr);
        // Dump the queue into the vector, where each element is in its appropriate
        // position / index
        TranscriptJob* tj;
        while ( processedReadQueue.dequeue(tj) ) { transcripts[tj->index] = tj; }
        std::cerr << "done processing read queue\n";
        */

        tbb::task_scheduler_init tbb_init;

        size_t numKmers = 0;
        std::atomic<size_t> mappedReads{0};
        
        tbb::concurrent_unordered_map<uint64_t, uint32_t> genePromiscuousKmers;
        tbb::parallel_for( _kmerGeneMap.range(),
        [&genePromiscuousKmers, &mappedReads, this]( KmerMapT::range_type & range ) {
            for ( auto kg = range.begin(); kg != range.end(); ++kg ) {
                auto kmer = kg->first;
                auto &geneSet = kg->second;
                if ( geneSet.size() > 1 ) {
                    genePromiscuousKmers[kmer] = geneSet.size();//.insert(kmer);
                }
                mappedReads += this->_readHash[kmer];
            }
        }
        );
        std::cerr << "I counted " << mappedReads << " mapped reads\n";
        std::cerr << " out of " << _kmerGeneMap.size() << " kmers, " << genePromiscuousKmers.size() << " were gene non-unique\n";
        std::cerr << 100.0 * (_kmerGeneMap.size() - genePromiscuousKmers.size()) / static_cast<double>(_kmerGeneMap.size()) <<
                  "% of the kmers are gene unique\n";
        
        std::vector<double> abundances( _transcriptGeneMap.numTranscripts(), 0.0 );

        std::atomic<size_t> numGenesProcessed{0};

        auto abundanceEstimateProgress = [&numGenesProcessed, this] () {
            auto numJobs = this->_transcriptGeneMap.numGenes();
            boost::progress_display show_progress( numJobs );
            size_t lastCount = numGenesProcessed;
            while ( lastCount < numJobs ) {
                auto diff = numGenesProcessed - lastCount;
                if ( diff > 0 ) {
                    show_progress += diff;
                    lastCount += diff;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
            }
        };
        std::cerr << "Processing abundances\n";
        tp.schedule(abundanceEstimateProgress);

        // process each gene
        tbb::parallel_for( size_t {0}, _transcriptGeneMap.numGenes(), size_t {1},
        [&transcriptKmers, &transcriptLengths, &abundances, &numGenesProcessed, 
            &mappedReads, &transcriptHash, &genePromiscuousKmers, this]( const size_t & geneID ) {
            auto abundance = this->_estimateIsoformAbundance(geneID, transcriptKmers, transcriptLengths, mappedReads,
                genePromiscuousKmers);
            //auto abundance = this->_estimateIsoformAbundanceByCount(geneID, transcriptKmers, transcriptLengths, mappedReads, transcriptHash);

            auto transcripts = this->_transcriptGeneMap.transcriptsForGene(geneID);
            for ( size_t i = 0; i < abundance.size(); ++i ) {
                //_transcriptGeneMap.transcriptName( transcripts[i] );
                abundances[ transcripts[i] ] = abundance[i];
            }
            ++numGenesProcessed;
        }
        );
        //for ( size_t geneID = 0; geneID < _transcriptGeneMap.numGenes(); ++geneID ) {
        //auto abundance = _estimateIsoformAbundance(geneID, transcriptKmers);
        //}

        std::ofstream afile("abundances.txt");
        for ( size_t ti = 0; ti < _transcriptGeneMap.numTranscripts(); ++ti ) {
            /*double effectiveLength = 0.0;
            for ( auto kmer : transcriptKmers[ti] ) {
                effectiveLength += (1.0 / (1.0 + transcriptHash[kmer]));
            }*/
            afile << _transcriptGeneMap.transcriptName(ti) << '\t' << transcriptLengths[ti] << '\t' << abundances[ti] << '\n';
        }
        afile.close();

        std::ofstream gfile("abundances.gene.txt");
        for ( size_t gi = 0; gi < _transcriptGeneMap.numGenes(); ++gi ) {
          auto transcripts = _transcriptGeneMap.transcriptsForGene(gi);
          double exp = 0.0;
          for ( auto t : transcripts ) { exp += abundances[t]; }

          gfile << _transcriptGeneMap.nameFromGeneID(gi) << '\t' << exp << '\n';
        }
        gfile.close();
        std::exit(1);

        // Our graph will have weighted edges
        struct EdgePropT {
            float weight;
        };
        typedef boost::property< boost::edge_index_t, size_t, EdgePropT > edge_property;

        // We'll be taking subgraphs, so make sure we implement this concept
        typedef typename boost::subgraph <
        boost::adjacency_list < boost::hash_setS,
              boost::vecS,
              boost::undirectedS,
              boost::no_property,
              edge_property
              > > GraphT;
        typedef typename boost::graph_traits<GraphT>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<GraphT>::vertices_size_type VertexIndex;

        // Create our transcript overlap graph
        auto numTranscripts = _transcriptGeneMap.numTranscripts();
        auto numGenes = _transcriptGeneMap.numGenes();
        GraphT G( numGenes );
        std::vector<VertexIndex> rank( numGenes );
        std::vector<Vertex> parent( numGenes );

        typedef VertexIndex *Rank;
        typedef Vertex *Parent;

        boost::disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);
        boost::initialize_incremental_components(G, ds);
        boost::incremental_components(G, ds);

        std::cerr << "Building transcript overlap graph\n";
        boost::progress_display showGraphProgress(_kmerGeneMap.size());
        for ( auto & kv : _kmerGeneMap ) {
            ++showGraphProgress;

            using std::tie;
            typename GraphT::edge_descriptor e;
            bool exists {false};

            auto kmer = kv.first;
            for ( auto u = kv.second.begin(); u != kv.second.end(); ++u ) {//& u : kv.second ) {
                auto v = u; ++v;
                for ( ; v != kv.second.end(); ++v ) {//& u : kv.second ) {
                    //for ( auto& v : kv.second ) {
                    auto vertU = *u;
                    auto vertV = *v;
                    if ( vertU > vertV ) {
                        std::swap(vertU, vertV);
                    }

                    tie(e, exists) = boost::edge( vertU, vertV, G );
                    if ( !exists ) {
                        auto e = boost::add_edge(vertU, vertV, G);
                        G[e.first].weight = 1.0f;
                        ds.union_set(vertU, vertV);
                    } else {
                        G[ e ].weight += 1.0;
                    }

                }
            }

        }

        typedef boost::component_index<VertexIndex> Components;

        // NOTE: Because we're using vecS for the graph type, we're
        // effectively using identity_property_map for a vertex index map.
        // If we were to use listS instead, the index map would need to be
        // explicitly passed to the component_index constructor.
        Components components(parent.begin(), parent.end());

        /**
        * Create a graph (appropriately indexed) for each
        * connected component of our original graph
        */
        // Vector to hold the new indices of each vertex
        // for vertex i in the original graph, it's new index is relabel[i]
        std::vector<size_t> relabel( boost::num_vertices(G), 0 );
        std::vector<size_t> relabelInv( boost::num_vertices(G), 0 );
        std::vector<size_t> offsets;
        size_t vnum {0};

        std::cerr << "building the connected component subgraphs \n";
        std::vector< GraphT > componentGraphs;
        size_t compID {0};
        //std::ofstream ofile("nnls_abundance.txt");
        // Iterate through the component indices
        BOOST_FOREACH(VertexIndex current_index, components) {

            // The offset index of this subgraph
            offsets.push_back(vnum);
            // The graph for this component
            auto compG = G.create_subgraph();//componentGraphs.back();

            // For each vertex in this component, add it to the graph
            // and calculate its label

            auto compIter = components[current_index];
            for ( auto &vert = compIter.first; vert != compIter.second; ++vert ) {
                boost::add_vertex(*vert, compG);
                relabel[*vert] = vnum;
                relabelInv[vnum] = *vert;
                ++vnum;
            }

            std::cerr << "subgraph " << compID << "; order = " << boost::num_vertices(compG) << "; size = " <<
                      boost::num_edges(compG) << "\n";

            // process the graph
            /*
            auto res = _computeAbundance( compG, transcriptHash, transcripts );
            for ( auto& r : res ) {
              ofile << *(r.name) << '\t' << r.abundance << '\n';
            }
            */
            //componentGraphs.push_back(compG);
            ++compID;
        }



        //ofile.close();


        std::ofstream ofile("sharedkmers.gra");

        // Write number of verts, edges and format (edge weights)
        ofile << boost::num_vertices(G) << " " << boost::num_edges(G)  << " 1\n";

        // We have to iterate over the vertices in numerical order
        // for this format
        for ( size_t i = 0; i < boost::num_vertices(G); ++i ) {
            // Iterate over all edges adjacent to vertex i
            auto se = boost::out_edges(i, G);
            for ( auto &e = se.first; e != se.second; ++e ) {
                // Output the pair (target, weight)
                ofile << boost::target(*e, G) << " " << G[*e].weight << " ";
            }
            ofile << "\n";
        }
        ofile.close();

        return true;
    } // end findOverlappingSegments

};


#endif //TRANSCRIPTSEGMENTER_HPP
