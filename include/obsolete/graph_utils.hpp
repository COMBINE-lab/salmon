#ifndef GRAPH_UTILS_HPP
#define GRAPH_UTILS_HPP

#include <type_traits>
#include <boost/graph/adjacency_list.hpp>

#include "mlrmcl/metis.h"

namespace utils {
	namespace graph {
	
	template <typename GraphT>//, typename IndexT, typename WeightT>
	bool boostToMetis( GraphT& G, graphdef* metisG) {

		std::cerr << "calling InitGraph\n";
	    // Initalize the graph structure
		InitGraph(metisG);
		std::cerr << "done\n";

  	    // Check if this graph has weighted edges or not
		bool hasEdgeWeights = !(std::is_same< typename GraphT::edge_property_type, typename boost::no_property >::value);
	    // Number of vertices in this graph
		auto numVerts = metisG->nvtxs = boost::num_vertices(G);

	    // Allocate the arrays for the vertex adj. offsets and adj. structure
		metisG->xadj = new idxtype[numVerts+1]; metisG->xadj[0] = 0;
		metisG->adjncy = new idxtype[ 2 * boost::num_edges(G) + numVerts ];

	    // If the original graph has edge weights, then allocate the memory for the
	    // corresponding array.
		if ( hasEdgeWeights ) {
			metisG->adjwgt = new idxtype[ 2 * boost::num_edges(G) + numVerts];
		}

		size_t edgeNum{0};
	    
	    // For each vertex in order
		for ( size_t i = 0; i < numVerts; ++i ) {
			bool selfLoop{false};

			// Keep track of the maximum edge weight (this is what
			// we'll set the self-loop edge weight to).
			idxtype maxWeight = std::numeric_limits<idxtype>::min();

		    // Iterate over all edges adjacent to vertex i
			auto se = boost::out_edges(i,G);
			for ( auto& e = se.first; e != se.second; ++e ) {
				auto tgt = boost::target(*e,G);
				if ( tgt == i ) { selfLoop = true; }

                // Output the pair (target, weight)
				metisG->adjncy[edgeNum] = tgt;
				if ( hasEdgeWeights ) { 
					auto weight = G[*e].weight; 
					metisG->adjwgt[edgeNum] = weight;
					maxWeight = (weight > maxWeight) ? weight : maxWeight;
				}
				++edgeNum;
			} 

		    // Ensure self-loops
			if ( !selfLoop ) { 
				metisG->adjncy[edgeNum] = i; 
				if ( hasEdgeWeights ) { metisG->adjwgt[edgeNum] = maxWeight; }
				++edgeNum; 
			}

			if ( i < numVerts - 1 ) { metisG->xadj[i+1] = edgeNum; }
		}
	    // # of edges

		// Apparently, this is necessary
		metisG->xadj[metisG->nvtxs] = edgeNum;

		metisG->nedges = edgeNum;
		return true;
	} // end of boostToMetis

	} // end of graph
} // end of utils

#endif // GRAPH_UTILS_HPP