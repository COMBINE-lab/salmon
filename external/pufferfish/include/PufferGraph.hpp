#ifndef PUFFER_GRAPH_H
#define PUFFER_GRAPH_H

#include "Util.hpp"
#include <map>

namespace puffergraph {
class Node {
public:
  Node() {}
  Node(std::string idIn) {
    id = idIn;
    indegp = 0;
    outdegp = 0;
    indegm = 0;
    outdegm = 0;
  }

  uint8_t getIndegP() { return indegp; }
  uint8_t getOutdegP() { return outdegp; }
  uint8_t getIndegM() { return indegm; }
  uint8_t getOutdegM() { return outdegm; }
  uint8_t getOutdeg() { return (outdegp + outdegm); }
  uint8_t getIndeg() { return (indegp + indegm); }
  uint8_t getRealOutdeg() { return (outdegp + indegm); }
  uint8_t getRealIndeg() { return (indegp + outdegm); }

  // std::vector<std::string> getPositiveToNodes()
  std::string getId() { return id; }

  std::vector<std::string> getToNeighbor() {
    std::vector<std::string> result;
    for (auto i : out) {
      if (i.baseSign) {
        result.push_back(i.contigId);
      }
    }
    for (auto i : in) {
      if (!i.baseSign) {
        result.push_back(i.contigId);
      }
    }

    return result;
  }

  std::vector<std::string> getFromNeighbor() {
    std::vector<std::string> result;
    for (auto i : out) {
      if (!i.baseSign) {
        result.push_back(i.contigId);
      }
    }
    for (auto i : in) {
      if (i.baseSign) {
        result.push_back(i.contigId);
      }
    }

    return result;
  }

  void insertNodeTo(std::string nodeId, bool sign, bool toSign) {
    if (toSign)
      outdegp++;
    else
      outdegm++;

    out.emplace_back(sign, nodeId, toSign);
  }
  void insertNodeFrom(std::string nodeId, bool sign, bool fromSign) {
    if (fromSign)
      indegp++;
    else
      indegm--;

    in.emplace_back(sign, nodeId, fromSign);
  }

private:
  std::string id;
  uint8_t indegp;
  uint8_t outdegp;
  uint8_t indegm;
  uint8_t outdegm;
  struct edgetuple {
    edgetuple(bool fSign, std::string cId, bool tSign)
        : baseSign(fSign), contigId(cId), neighborSign(tSign) {}

    bool baseSign;
    std::string contigId;
    bool neighborSign;
  };
  std::vector<edgetuple> in;
  std::vector<edgetuple> out;
};

class Graph {
private:
  std::map<std::string, Node> Vertices;
  // std::vector<std::pair<Node,Node> > Edges ;

public:
  // case where I do
  std::map<std::string, Node> getVertices() { return Vertices; }

  bool getNode(std::string nodeId) {
    if (Vertices.find(nodeId) == Vertices.end())
      return true;
    else
      return false;
  }

  void addEdge(std::string fromId, bool fromSign, std::string toId,
               bool toSign) {
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

    fromNode.insertNodeTo(toNode.getId(), fromSign, toSign);
    toNode.insertNodeFrom(fromNode.getId(), toSign, fromSign);

    // Edges.emplace_back(fromNode,toNode) ;
  }

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
