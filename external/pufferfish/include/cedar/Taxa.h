#ifndef __TAXA_H__
#define __TAXA_H__

#include <string>
#include <vector>
#include <set> // std::set
#include "iostream"

#define NO_PARENT -1
typedef uint16_t rLenType;
typedef uint32_t refLenType;

enum class ReadEnd : uint8_t {
    LEFT, RIGHT
};
// inclusive for begin 
// exclusive for end
struct Interval {
    uint32_t begin;
    uint32_t end;
    Interval() {}
    Interval(uint32_t inBegin, uint32_t inEnd) : begin(inBegin), end(inEnd) {}

    bool operator==(const Interval& rhs) const {
        return begin == rhs.begin and end == rhs.end;
    }
};

enum class Rank : uint8_t {
    STRAIN, FORMA, VARIETAS, // strain sub group
    SUBSPECIES, SPECIES, SPECIESSUBGROUP, SPECIESGROUP, // species sub group 
    SUBGENUS, GENUS, SUPERGENUS, // genus sub group
    SUBTRIBE, TRIBE, SUBFAMILY, FAMILY, SUPERFAMILY, // family sub group
    PARVORDER, INFRAORDER, SUBORDER, ORDER, SUPERORDER, // order sub group
    COHORT, INFRACLASS, SUBCLASS, CLASS, SUPERCLASS, // class sub group
    SUBPHYLUM, PHYLUM, SUPERPHYLUM, // phylum sub group
    SUBKINGDOM, KINGDOM, SUPERKINGDOM , DOMAINE, LIFE, // general sub group
    UNCLASSIFIED
};


class TaxaNode {
    public:
        TaxaNode() { 
            id = static_cast<uint64_t>(-1);
            rank = Rank::STRAIN; 
            score = 0; 
            notIncorporatedChildrenCounter = 0;
						paired = false;
        }
        TaxaNode(uint64_t inId) : id(inId) {paired = false;}
        TaxaNode(uint64_t inId, uint64_t inPid, Rank inRank) : 
            id(inId), parentId(inPid), rank(inRank),
            score(0), notIncorporatedChildrenCounter(0) 
             {
                // we go into the loop just for 1
                if (id == parentId) {
                    rank = Rank::LIFE;
                    parentId = NO_PARENT;
                    //std::cerr << id << "\n";
                }
								paired = false;
            }
        bool isRoot() { return rank == Rank::LIFE; }//TODO not easy with the new design return children.size(); }
        bool isRipe() { 
            return !notIncorporatedChildrenCounter;
        } // ripe if zero
        void addInterval(uint64_t begin, uint64_t len, ReadEnd readEnd);
        void updateIntervals(TaxaNode* child, ReadEnd readEnd);
        void updateScore();
        /**
         * Sorts intervals
         * Merge intervals if possible
         * Calculates score
         **/
        void cleanIntervals(ReadEnd readEnd);
        bool addChild(TaxaNode* child);
        void setOneMoreChildAsProcessed() {
            notIncorporatedChildrenCounter--;
        }
        uint64_t getParentId() {return parentId;}
        uint64_t getId() {return id;}
        void setId(uint64_t inId) {id = inId;}
        Rank getRank() {return rank;}
        uint64_t getScore() {return score;}
        std::set<uint64_t>& getActiveChildren() {return activeChildren;}
        std::vector<Interval>& getIntervals(ReadEnd readEnd) {return readEnd == ReadEnd::LEFT?lintervals:rintervals;}
        refLenType getPos(ReadEnd re) {if(re == ReadEnd::LEFT) return lPos; return rPos;}
        bool isFw(ReadEnd re) {if (re == ReadEnd::LEFT) return leftFw; return rightFw;}
        void setPos(refLenType pos, ReadEnd re) {re == ReadEnd::LEFT?lPos = pos:rPos = pos;}
        void setFw(bool isFw, ReadEnd re) {re == ReadEnd::LEFT?leftFw = isFw:rightFw = isFw;}

    void setScore(uint32_t scoreIn) { score = scoreIn; }
        bool isConcordant() {
            return lintervals.size() && rintervals.size();
        }
        bool compareIntervals(TaxaNode& other);
        void reset();
        static Rank str2rank(std::string rankstr) {
            if (rankstr == "no rank") return Rank::STRAIN;
            else if (rankstr == "varietas") return Rank::VARIETAS;
            else if (rankstr == "subspecies") return Rank::SUBSPECIES;
            else if (rankstr == "species") return Rank::SPECIES;
            else if (rankstr == "species subgroup") return Rank::SPECIESSUBGROUP;
            else if (rankstr == "species group") return Rank::SPECIESGROUP;
            else if (rankstr == "subgenus") return Rank::SUBGENUS;
            else if (rankstr == "genus") return Rank::GENUS;
            else if (rankstr == "supergenus") return Rank::SUPERGENUS;
            else if (rankstr == "subfamily") return Rank::SUBFAMILY;
            else if (rankstr == "family") return Rank::FAMILY;
            else if (rankstr == "superfamily") return Rank::SUPERFAMILY;
            else if (rankstr == "subtribe") return Rank::SUBTRIBE;
            else if (rankstr == "tribe") return Rank::TRIBE;
            else if (rankstr == "forma") return Rank::FORMA;
            else if (rankstr == "cohort") return Rank::COHORT;
            else if (rankstr == "parvorder") return Rank::PARVORDER;
            else if (rankstr == "suborder") return Rank::SUBORDER;
            else if (rankstr == "order") return Rank::ORDER;
            else if (rankstr == "infraorder") return Rank::INFRAORDER;
            else if (rankstr == "superorder") return Rank::SUPERORDER;
            else if (rankstr == "subclass") return Rank::SUBCLASS;
            else if (rankstr == "class") return Rank::CLASS;
            else if (rankstr == "infraclass") return Rank::INFRACLASS;
            else if (rankstr == "superclass") return Rank::SUPERCLASS;
            else if (rankstr == "subphylum") return Rank::SUBPHYLUM;
            else if (rankstr == "phylum") return Rank::PHYLUM;
            else if (rankstr == "superphylum") return Rank::SUPERPHYLUM;
            else if (rankstr == "subkingdom") return Rank::SUBKINGDOM;
            else if (rankstr == "kingdom") return Rank::KINGDOM;
            else if (rankstr == "superkingdom") return Rank::SUPERKINGDOM;
            else if (rankstr == "domain") return Rank::DOMAINE;
            else if (rankstr == "life") return Rank::LIFE; 
            std::cerr << "ERROR: Not a valid rank: --" << rankstr << "--\n";
            std::exit(1);
        }
        static std::string rank2str(Rank r) {
            switch (r) {
                case Rank::STRAIN:
                    return "no rank";
                case Rank::VARIETAS:
                    return "varietas";
                case Rank::SUBSPECIES:
                    return "subspecies";
                case Rank::SPECIES:
                    return "species";
                case Rank::SPECIESSUBGROUP:
                    return "species subgroup";
                case Rank::SPECIESGROUP:
                    return "species group";
                case Rank::SUBGENUS:
                    return "subgenus";
                case Rank::GENUS:
                    return "genus";
                case Rank::SUPERGENUS:
                    return "supergenus";
                case Rank::SUBFAMILY:
                    return "subfamily";
                case Rank::FAMILY:
                    return "family";
                case Rank::SUPERFAMILY:
                    return "superfamily";
                case Rank::SUBTRIBE:
                    return "subtribe";
                case Rank::TRIBE:
                    return "tribe";
                case Rank::FORMA:
                    return "forma";
                case Rank::COHORT:
                    return "cohort";
                case Rank::PARVORDER:
                    return "parvorder";
                case Rank::SUBORDER:
                    return "suborder";
                case Rank::ORDER:
                    return "order";
                case Rank::INFRAORDER:
                    return "infraorder";
                case Rank::SUPERORDER:
                    return "superorder";
                case Rank::SUBCLASS:
                    return "subclass";
                case Rank::CLASS:
                    return "class";
                case Rank::INFRACLASS:
                    return "infraclass";
                case Rank::SUPERCLASS:
                    return "superclass";
                case Rank::SUBPHYLUM:
                    return "subphylum";
                case Rank::PHYLUM:
                    return "phylum";
                case Rank::SUPERPHYLUM:
                    return "superphylum";
                case Rank::SUBKINGDOM:
                    return "subkingdom";
                case Rank::KINGDOM:
                    return "kingdom";
                case Rank::SUPERKINGDOM:
                    return "superkingdom";
                case Rank::DOMAINE:
                    return "domain";
                case Rank::LIFE:
                    return "life";
                case Rank::UNCLASSIFIED:
                    return "unclassified";
            }
            return "really!! none?";
        }
				void setPaired() { paired = true; }
				bool isPaired() { return paired; }
    private:
        // fix for all the mappings (coming from taxonomy tree structure)
        uint64_t id;
        uint64_t parentId;
        Rank rank;
        // change per read mapping
        uint32_t score{0};
        uint64_t notIncorporatedChildrenCounter;
        std::set<uint64_t> activeChildren; //it's a set because, we might add the same child multiple times
        bool leftFw{true};
        bool rightFw{false};
        refLenType lPos{0};
        refLenType rPos{0};
        std::vector<Interval> lintervals;
        std::vector<Interval> rintervals;
        bool isCleaned = false;
				bool paired;
};

struct TaxaInfo {
    uint64_t cnt;
    uint64_t subTreeCnt;
    Rank rank;
    TaxaInfo() {};
    TaxaInfo(uint64_t cntIn, Rank rankIn) : cnt(cntIn), subTreeCnt(cntIn), rank(rankIn) {}
    void increment() {
        cnt++;
        subTreeCnt++;
    }
};

#endif
