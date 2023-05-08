#ifndef PUFFER_GRAPH_H
#define PUFFER_GRAPH_H

#include "Util.hpp"
#include "sparsepp/spp.h"
#include <algorithm>
#include <iostream>
//#include <map>

namespace pufg {

enum class EdgeType : uint8_t {
  PLUS_PLUS = 0,
  PLUS_MINUS = 1,
  MINUS_PLUS = 2,
  MINUS_MINUS = 3
};

inline bool fromSign_(EdgeType et) {
  return (et == EdgeType::PLUS_PLUS or et == EdgeType::PLUS_MINUS);
}
inline bool toSign_(EdgeType et) {
  return (et == EdgeType::PLUS_PLUS or et == EdgeType::MINUS_PLUS);
}

inline EdgeType typeFromBools_(bool sign, bool toSign) {
  if (sign and toSign) {
    return EdgeType::PLUS_PLUS;
  } else if (sign and !toSign) {
    return EdgeType::PLUS_MINUS;
  } else if (!sign and toSign) {
    return EdgeType::MINUS_PLUS;
  } else {
    return EdgeType::MINUS_MINUS;
  }
}

struct edgetuple {
  edgetuple(bool fSign, uint64_t cId, bool tSign) : contigId(cId) {
    t = typeFromBools_(fSign, tSign);
  }

  edgetuple() {}

  bool baseSign() {
    return (t == EdgeType::PLUS_PLUS or t == EdgeType::PLUS_MINUS);
  }
  bool neighborSign() {
    return (t == EdgeType::PLUS_PLUS or t == EdgeType::MINUS_PLUS);
  }

  EdgeType t;
  uint64_t contigId;
};

inline bool operator==(const edgetuple& e1, const edgetuple& e2) {
  return (e1.t == e2.t and e1.contigId == e2.contigId);
}

inline bool operator!=(const edgetuple& e1, const edgetuple& e2) {
  return !(e1 == e2);
}

class Node {
  template <typename E>
  constexpr typename std::underlying_type<E>::type to_index(E e) {
    return static_cast<typename std::underlying_type<E>::type>(e);
  }
  int8_t getCountPlus_(std::vector<edgetuple>& elist) {
    int8_t r{0};
    for (auto& et : elist) {
      r += (et.baseSign() ? 1 : 0);
    }
    return r;
  }
  int8_t getCountMinus_(std::vector<edgetuple>& elist) {
    int8_t r{0};
    for (auto& et : elist) {
      r += (!et.baseSign() ? 1 : 0);
    }
    return r;
  }

public:
  Node() {}
  Node(uint64_t idIn) { id = idIn; }

  int8_t getIndegP() { return getCountPlus_(in_edges); }

  int8_t getOutdegP() { return getCountPlus_(out_edges); }

  int8_t getIndegM() { return getCountMinus_(in_edges); }

  int8_t getOutdegM() { return getCountMinus_(out_edges); }

  size_t getRealOutdeg() {
    // outgoing + and incoming -
    spp::sparse_hash_set<std::string> distinctNeighbors;
    for (auto& e : in_edges) {
      if (!e.baseSign()) {
        std::string name = std::to_string(e.contigId) + (e.neighborSign() ? "-" : "+");
        distinctNeighbors.insert(name);
      }
    }
    for (auto& e : out_edges) {
      if (e.baseSign()) {
        std::string name = std::to_string(e.contigId) + (e.neighborSign() ? "+" : "-");
        distinctNeighbors.insert(name);
      }
    }
    return distinctNeighbors.size();
  }
  size_t getRealIndeg() {
    // outgoing - and incoming +
    spp::sparse_hash_set<std::string> distinctNeighbors;
    for (auto& e : in_edges) {
      if (e.baseSign()) {
        std::string name = std::to_string(e.contigId) + (e.neighborSign() ? "+" : "-");
        distinctNeighbors.insert(name);
      }
    }
    for (auto& e : out_edges) {
      if (!e.baseSign()) {
        std::string name = std::to_string(e.contigId) + (e.neighborSign() ? "-" : "+");
        distinctNeighbors.insert(name);
      }
    }
    return distinctNeighbors.size();
  }

  edgetuple& getOnlyRealIn() {
    if (getIndegP() > 0) {
      for (auto& e : in_edges)
        if (e.baseSign()) {
          return e;
        }
    } else {
      for (auto& e : out_edges)
        if (!e.baseSign()) {
          return e;
        }
    }
    // should not get here
    return in_edges.front();
  }

  edgetuple& getOnlyRealOut() {
    if (getOutdegP() > 0) {
      for (auto& e : out_edges)
        if (e.baseSign()) {
          return e;
        }
    } else { // The real outgoing edge should be an incoming edge to negative if
             // it's not an outgoing edge from positive
      for (auto& e : in_edges)
        if (!e.baseSign()) {
          return e;
        }
    }
    // should not get here
    return out_edges.front();
  }

  uint64_t& getId() { return id; }

  void insertNodeTo(uint64_t nodeId, bool sign, bool toSign) {
    edgetuple ekey = {sign, nodeId, toSign};
    if (std::find(out_edges.begin(), out_edges.end(), ekey) ==
        out_edges.end()) {
      out_edges.emplace_back(ekey);
    }
  }

  void removeEdgeTo(uint64_t nodeId) {
    out_edges.erase(std::remove_if(out_edges.begin(), out_edges.end(),
                                   [&nodeId](edgetuple& etup) -> bool {
                                     return etup.contigId == nodeId;
                                   }),
                    out_edges.end());
  }

  void insertNodeFrom(uint64_t nodeId, bool sign, bool fromSign) {
    edgetuple ekey = {sign, nodeId, fromSign};
    if (std::find(in_edges.begin(), in_edges.end(), ekey) == in_edges.end()) {
      in_edges.emplace_back(ekey);
    }
  }

  void removeEdgeFrom(uint64_t nodeId) {
    in_edges.erase(std::remove_if(in_edges.begin(), in_edges.end(),
                                  [&nodeId](edgetuple& etup) -> bool {
                                    return etup.contigId == nodeId;
                                  }),
                   in_edges.end());
  }

  bool checkExistence(bool bSign, uint64_t toId, bool toSign) {
    edgetuple ekey = {bSign, toId, toSign};
    return (std::find(out_edges.begin(), out_edges.end(), ekey) !=
            out_edges.end());
  }

  std::vector<edgetuple>& getPredecessors() { return in_edges; }
  std::vector<edgetuple>& getSuccessors() { return out_edges; }

private:
  // TODO: Make the actual node IDs into numbers instead of strings
  // uint64_t id_;
  uint64_t id;
  std::vector<edgetuple> out_edges;
  std::vector<edgetuple> in_edges;
};

class Graph {
private:
  spp::sparse_hash_map<uint64_t, Node> Vertices;
  // std::vector<Node> Vertices;
  // std::vector<std::string> NodeNames;
  // std::vector<std::pair<Node,Node> > Edges ;

public:
  // case where I do
  spp::sparse_hash_map<uint64_t, Node>& getVertices() { return Vertices; }

  bool getNode(uint64_t nodeId) {
    // return (Vertices.find(nodeId) == Vertices.end());
    if (Vertices.find(nodeId) == Vertices.end())
      return true;
    else
      return false;
  }

  bool addEdge(uint64_t fromId, bool fromSign, uint64_t toId,
               bool toSign) {
//		  if (fromId == 70545 or fromId == 70547 or toId == 70545 or toId == 70547)
//				  std::cerr << "ADD EDGE FROM " << fromId << " " << fromSign << " TO " << toId << " " << toSign << "\n";
    // case 1 where the from node does not exist
    // None of the nodes exists
    if (Vertices.find(fromId) == Vertices.end()) {
      Node newNode(fromId);
      Vertices[fromId] = newNode;
    }
    if (Vertices.find(toId) == Vertices.end()) {
      Node newNode(toId);
      Vertices[toId] = newNode;
    }
    auto& fromNode = Vertices[fromId];
    auto& toNode = Vertices[toId];

    if (!fromNode.checkExistence(fromSign, toNode.getId(), toSign)) {
      fromNode.insertNodeTo(toNode.getId(), fromSign, toSign);
      toNode.insertNodeFrom(fromNode.getId(), toSign, fromSign);
      return true;
    }

    return false;

    // Edges.emplace_back(fromNode,toNode) ;
  }

  bool merge(uint64_t tobeMerged, bool tobeMergedOri, uint64_t remaining, bool remainingOri, bool isMergeRight) {
		if (!Vertices.contains(tobeMerged) or !Vertices.contains(remaining)) return false;
		Node& tobeMergedNode = Vertices[tobeMerged];

		for (auto& in : tobeMergedNode.getPredecessors()) {
			if (isMergeRight != in.baseSign()){ 
					if (in.contigId == remaining and ( (in.neighborSign() == remainingOri and in.baseSign() == tobeMergedOri) or (in.neighborSign() != remainingOri and in.baseSign() != tobeMergedOri) ) ) continue;
					std::cerr << "BUG BUG!! IN GRAPH MERGING Predecessors. isMergeRight=" <<isMergeRight << " " << tobeMerged << "," << tobeMergedOri << " " << remaining << "," << remainingOri << " buggi: " << in.contigId << "," <<in.neighborSign() << "\n";
			}
	        addEdge(in.contigId, in.neighborSign(), remaining, isMergeRight?remainingOri == tobeMergedOri:remainingOri!=tobeMergedOri);
		}
		for (auto& out : tobeMergedNode.getSuccessors()) {			
			if (isMergeRight == out.baseSign()) { 
					if (out.contigId == remaining and ( (out.neighborSign() == remainingOri and out.baseSign() == tobeMergedOri) or (out.neighborSign() != remainingOri and out.baseSign() != tobeMergedOri) ) ) continue;
					std::cerr << "BUG BUG!! IN GRAPH MERGING Successors. isMergeRight=" << isMergeRight << " " << tobeMerged << "," << tobeMergedOri << " " << remaining << "," << remainingOri << "\n";
			}
	        addEdge(remaining, isMergeRight?remainingOri != tobeMergedOri:remainingOri==tobeMergedOri, out.contigId, out.neighborSign());					
		}
		removeNode(tobeMerged);
		return true;

  }

  bool removeNode(uint64_t id) {
	Node& n = Vertices[id];
    for (auto& in : n.getPredecessors()) {
      Node& from = Vertices[in.contigId];
      from.removeEdgeTo(id);
    }
    for (auto& out : n.getSuccessors()) {
      Node& to = Vertices[out.contigId];
      to.removeEdgeFrom(id);
    }
	Vertices.erase(n.getId());
    return false;
  }

  /*
      bool removeNode_orig(std::string id) {
          if (Vertices.find(id) == Vertices.end()) return false;
          Node& n = Vertices[id];
          for (auto & in : n.getIn()) {
              for (auto & out : n.getOut()) {
                  addEdge(in.contigId, in.neighborSign(), out.contigId,
  out.neighborSign());
              }
          }
          for (auto& in : n.getIn()) {
              Node& from = Vertices[in.contigId];
              from.removeEdgeTo(id);
          }
          for (auto& out : n.getOut()) {
              Node& to = Vertices[out.contigId];
              to.removeEdgeFrom(id);
          }
    return false;
  }
  */

  /*
  void buildGraph(std::string gfa_file){
      std::ifstream file(gfa_file);
      while(std::get_line(file, ln)){
          char firstC = ln[0];
          if(firstC != 'L') continue ;

      }
  }*/
};
}

#endif // end puffer
