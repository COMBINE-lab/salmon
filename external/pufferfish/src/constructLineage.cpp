//
// Created by Fatemeh Almodaresi on 2020-06-19.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <sstream>

constexpr uint64_t noRank = 1;

struct taxInfo {
    uint64_t parent{noRank}, lowestSpecies{noRank}, highestSpecies{noRank}, genus{noRank}, family{noRank}, order{
            noRank}, kingdom{noRank};
    std::vector<uint64_t> children;
    std::string rank{"NA"};

    taxInfo(uint64_t p, std::string &r) : parent(p), rank(r) {}

    taxInfo() = default;

    void copyLineage(taxInfo &p) {
        lowestSpecies = p.lowestSpecies;
        highestSpecies = p.highestSpecies;
        genus = p.genus;
        family = p.family;
        order = p.order;
        kingdom = p.kingdom;
    }

};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Missing taxa file as input\n";
        std::exit(3);
    }


// parse file
    std::cerr << argv[1] << "\n";
    std::ifstream taxaFile(argv[1]);

    std::vector<taxInfo> tlist;
    std::string col1, col2, col3;
    uint64_t taxa{0}, parent{0};
    std::string rank;
    if (taxaFile.good()) {
        std::string s;
        if (!getline( taxaFile, s )) {
            std::cerr << "File is empty\n";
            std::exit(3);
        }

        std::istringstream ss( s );
        if (!getline( ss, col1, ',' ))  {
            std::cerr << "Expect 3 columns of taxa, parent, rank separated by comma\n";
            std::exit(3);
        }
        if (!getline( ss, col2, ',' ))  {
            std::cerr << "Expect 3 columns of taxa, parent, rank separated by comma\n";
            std::exit(3);
        }
        if (!getline( ss, col3, ',' ))  {
            std::cerr << "Expect 3 columns of taxa, parent, rank separated by comma\n";
            std::exit(3);
        }
    }
    uint64_t lineCntr{0};
    while (taxaFile.good()) {
        std::string s;
        if (!getline( taxaFile, s )) break;
        std::istringstream ss( s );
        if (!getline( ss, rank, ',' ))  {
            std::cerr << "Expect 3 columns of taxa, parent, rank separated by comma\n";
            std::exit(3);
        }
        taxa = std::stoull(rank);

        if (!getline( ss, rank, ',' ))  {
            std::cerr << "Expect 3 columns of taxa, parent, rank separated by comma\n";
            std::exit(3);
        }
        parent = std::stoull(rank);

        if (!getline( ss, rank, ',' ))  {
            std::cerr << "Expect 3 columns of taxa, parent, rank separated by comma\n";
            std::exit(3);
        }

        lineCntr++;
        if (tlist.size() <= taxa) {
            tlist.resize(taxa + 1);
        }
        tlist[taxa] = taxInfo(parent, rank);
    }
    taxaFile.close();
    std::cerr << "Lines: " << lineCntr << "\n";


    // set children
    for (uint64_t i = 0; i < tlist.size(); i++) {
        if (tlist[i].parent != i) {
            tlist[tlist[i].parent].children.push_back(i);
        }
    }

// calculate main seven level lineage via BFS
    std::queue<uint64_t> q;
    q.push(1);
    while (not q.empty()) {
        auto curIdx = q.front();
        q.pop();
        auto par = tlist[curIdx];
        for (auto c : par.children) {
            auto &child = tlist[c];
            child.copyLineage(par);
            if (child.rank == "species") {
                child.lowestSpecies = c;
                if (child.highestSpecies ==
                    noRank) { // The node is the first species observed in the path of root to node
                    child.highestSpecies = c;
                }
            } else if (child.rank == "genus") {
                child.genus = c;
            } else if (child.rank == "family") {
                child.family = c;
            } else if (child.rank == "order") {
                child.order = c;
            } else if (child.rank == "kingdom") {
                child.kingdom = c;
            }
            q.push(c);
        }
    }

    // store lineage
    std::ofstream linFile(std::string(argv[1]) + ".lineage", std::ios::out);
    linFile << col1 << ","
            << col2 << ","
            << col3 << ","
            << "lowestSpecies" << ","
            << "highestSpecies" << ","
            << "genus" << ","
            << "family" << ","
            << "order" << ","
            << "kingdom"
            << "\n";
    for (uint64_t i = 0; i < tlist.size(); i++) {
        linFile << i << ","
                << tlist[i].parent << ","
                << tlist[i].rank << ","
                << tlist[i].lowestSpecies << ","
                << tlist[i].highestSpecies << ","
                << tlist[i].genus << ","
                << tlist[i].family << ","
                << tlist[i].order << ","
                << tlist[i].kingdom
                << "\n";
    }
    linFile.close();

    return EXIT_SUCCESS;
}