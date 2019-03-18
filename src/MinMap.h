//
// Created by Ales Varabyou on 3/16/19.
// ====================================
// This class contains a map for storing
// kmers and their positions
//
// two types of maps can be constructed: transcriptomic, which can have kmers that span introns; and genomic
// right now the class will store vectors of coordinates for each kmer

#ifndef HAIRPIN_MINMAP_H
#define HAIRPIN_MINMAP_H

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>

// this class contains a single contiguous stretch of coordinates
class EPair{
public:
    EPair() = default;
    EPair(uint32_t start, uint32_t end){
        this->start=start;
        this->end=end;
    }
    ~EPair() = default;

    uint32_t getStart(){return this->start;}
    uint32_t getEnd(){return this->end;}

private:
    uint32_t start{};
    uint32_t end{};
};

// this class contains the coordinates spanned by a sequence such as a kmer
class EVec{
public:
    EVec() = default;
    EVec(uint8_t chrID,uint8_t strand){
        this->chrID=chrID;
        this->strand=strand;
        this->size=0;
    }
    ~EVec() = default;

    void _push_back(uint32_t start, uint32_t end){
        coords.push_back(new EPair(start,end));
        size++;
    };
    EPair* _get(int pos){return coords[pos];}
    uint8_t _getChr(){return this->chrID;}
    uint8_t _getStrand(){return this->strand;}
    int _getSize(){return this->size;}
    uint32_t _getStart(int pos){return this->coords[pos]->getStart();}
    uint32_t _getEnd(int pos){return this->coords[pos]->getEnd();}

private:
    std::vector<EPair*> coords{};
    uint8_t chrID{};
    uint8_t strand{};
    int size{};
};

class MinMap {
public:
    MinMap();
    ~MinMap();

    typedef std::map<std::string,std::vector< EVec> > KmerMap;
    KmerMap::iterator insert(std::string key, int chrID,int strand); // insert key and initiate a new EVec with given a given chrID and strand information

private:
    std::map<std::string,std::vector<EVec> > minmap{};


};

#endif //HAIRPIN_MINMAP_H
