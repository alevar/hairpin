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
    void incStart(){this->start++;}

    friend std::ostream& operator<<(std::ostream& os, const EPair& ep){
        os<<ep.start<<"-"<<ep.end;
        return os;
    }

    bool operator< (const EPair& ep) const{
        return (this->start < ep.start || this->end < ep.end);
    }
    bool operator> (const EPair& ep) const{
        return (this->start > ep.start);
    }
    bool operator== (const EPair& ep) const{
        return (this->start == ep.start && this->end == ep.end);
    }


private:
    uint32_t start{};
    uint32_t end{};
};

class EVec{
public:
    EVec(){this->size=0;};
    EVec(uint8_t chrID,uint8_t strand){
        this->chrID=chrID;
        this->strand=strand;
        this->size=0;
    }
    ~EVec() = default;

    void _push_back(uint32_t start, uint32_t end){
        coords.emplace_back(start,end);
        size++;
    }
    void _pop_back(){
        this->coords.pop_back();
        this->size--;
    }

    EPair _get(int pos){return coords[pos];}
    uint8_t _getChr(){return this->chrID;}
    uint8_t _getStrand(){return this->strand;}
    int _getSize(){return this->size;}
    uint32_t _getStart(int pos){return this->coords[pos].getStart();}
    uint32_t _getEnd(int pos){return this->coords[pos].getEnd();}
    void _incStart(){this->coords[0].incStart();}
    void _eraseFirst(){
        this->coords.erase(this->coords.begin());
        this->size--;
    }
    void _erase(){
        this->coords.erase(this->coords.begin(),this->coords.end());
        this->size=0;
    }

    typedef std::vector<EPair>::iterator iterator;
    typedef std::vector<EPair>::const_iterator const_iterator;
    typedef std::vector<EPair>::reference reference;
    iterator begin() {return coords.begin();}
    const_iterator begin() const { return coords.begin();}
    iterator end() {return coords.end();}
    const_iterator end() const { return coords.end();}

    friend std::ostream& operator<<(std::ostream& os, const EVec& ev){
        os<<(int)ev.chrID<<":"<<(int)ev.strand<<"@";
        for(auto it_ep: ev){
            os<<it_ep<<",";
        }
        return os;
    }

    bool operator< (const EVec& ev) const{
        return (this->chrID < ev.chrID || this->strand < ev.strand || this->coords <  ev.coords);
    }
    bool operator> (const EVec& ev) const{
        return (this->chrID > ev.chrID || this->strand > ev.strand || this->coords >  ev.coords);
    }
    bool operator==(const EVec& ev) const{
        return (this->chrID == ev.chrID && this->strand == ev.strand && this->coords ==  ev.coords);
    }

    void flipStrand(){
        if(this->strand=='-'){
            this->strand='+';
        }
        else{
            this->strand='-';
        }
    }

private:
    std::vector<EPair> coords{};
    uint8_t chrID{};
    uint8_t strand{};
    int size{};
};

class MinMap {
public:
    MinMap()=default;
    ~MinMap()=default;

    typedef std::map<std::string,std::vector<EVec> > KmerMap;

    KmerMap::iterator _find(std::string key){
        return minmap.find(key);
    }

    // this method inserts a new kmer if exists, or updates with a new kmer if does not exist
    KmerMap::iterator _insert(std::string key, uint8_t chrID,uint8_t strand){
        std::pair<KmerMap::iterator,bool> mm_it = minmap.insert(std::pair<std::string,std::vector<EVec>>(key,{}));
        if(mm_it.second){ // the kmer previously did not exist
            numKmer++;
        }
        mm_it.first->second.emplace_back(chrID,strand);
        return mm_it.first;
    } // insert key and initiate a new EVec with given a given chrID and strand information
    KmerMap::iterator _insert(std::string key, uint8_t chrID,uint8_t strand, uint32_t start, uint32_t end){ // initiate with the coordinate pair
        std::pair<KmerMap::iterator,bool> mm_it = minmap.insert(std::pair<std::string,std::vector<EVec>>(key,{}));
        if(mm_it.second){ // the kmer previously did not exist
            numKmer++;
        }
        mm_it.first->second.emplace_back(chrID,strand);
        mm_it.first->second.back()._push_back(start,end);
        return mm_it.first;
    }
    KmerMap::iterator _insert(std::string key,EVec ev){
        std::pair<KmerMap::iterator,bool> mm_it = minmap.insert(std::pair<std::string,std::vector<EVec>>(key,{}));
        if(mm_it.second){ // the kmer previously did not exist
            numKmer++;
        }
        mm_it.first->second.push_back(ev);
        return mm_it.first;
    }
    KmerMap::iterator _insert(std::string key,std::vector<EVec> vev){
        std::pair<KmerMap::iterator,bool> mm_it = minmap.insert(std::make_pair(key,vev));
        if(mm_it.second){ // the kmer previously did not exist
            numKmer++;
        }
        return mm_it.first;
    }

    // this method adds coordinates to the kmer at a given map given an iterator
    void _add(KmerMap::iterator mm_it,uint32_t start,uint32_t end){
        mm_it->second.back()._push_back(start,end);
    }
    int _getNumMultimappers(std::string key){ // this method returns the number of multimappers the given kmer has
        return minmap[key].size();
    }
    int _getNumKmers(){
        return numKmer;
    }

    typedef KmerMap::iterator iterator;
    typedef KmerMap::const_iterator const_iterator;
    typedef KmerMap::reference reference;
    iterator begin() {return minmap.begin();}
    const_iterator begin() const { return minmap.begin();}
    iterator end() {return minmap.end();}
    const_iterator end() const { return minmap.end();}

    friend std::ostream& operator<<(std::ostream& os, const MinMap& mm){
        for(auto it_mm: mm){
            os<<it_mm.first<<"\t";
            for(auto it_ev: it_mm.second){
                os<<it_ev<<";";
            }
            os<<std::endl;
        }
        return os;
    }

private:
    KmerMap minmap{};
    int numKmer{}; // number of kmers currently stored
};

#endif //HAIRPIN_MINMAP_H
