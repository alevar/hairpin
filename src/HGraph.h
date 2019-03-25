//
// Created by Ales Varabyou on 3/16/19.
// ====================================
// this class contains the implementation
// of a graph populated with kmers from sequenced reads
//

#ifndef HAIRPIN_GRAPH_H
#define HAIRPIN_GRAPH_H

#include <iostream>
#include <fstream>
#include <map>

#include "HDB.h"
#include "MinMap.h"

class Vertex;

class Edge{
public:
    Edge()=default;
    Edge(std::set<Vertex>::iterator vit_1,std::set<Vertex>::iterator vit_2, bool known){
        this->origin=vit_1;
        this->destination=vit_2;
        this->known=known;
        this->weight=1;
    }

    std::set<Vertex>::iterator getOrigin() {return origin;}
    std::set<Vertex>::iterator getDestination() {return destination;}
    void incWeight(){this->weight++;}

//    bool operator< (const Edge& et) const{
//        return (*origin < *et.origin || *destination < *et.destination);
//    }
//    bool operator> (const Edge& vt) const{
//        return (*this->origin > *vt.origin || *this->destination > *vt.destination);
//    }
//    bool operator==(const Edge& vt) const{
//        return (*this->origin == *vt.origin && *this->destination == *vt.destination);
//    }

private:
    std::set<Vertex>::iterator origin;
    std::set<Vertex>::iterator destination;
    bool known{}; // set to true if the edge is based on the transcriptome
    int weight=0; // number of reads supporting the edge
};

class Vertex{
public:
    Vertex(){this->weight=1;}

    ~Vertex()=default;

    int getWeight() {return this->weight;}
    void incWeight(){this->weight++;}
    void addOutEdge(std::set<Edge>::iterator eit){this->outEdges.emplace_back(eit);}
    void addInEdge(std::set<Edge>::iterator eit){this->inEdges.emplace_back(eit);}
    int getOutDegree(){return static_cast<int>(this->outEdges.size());}
    int getInDegree(){return static_cast<int>(this->inEdges.size());}

private:
    int weight=0;
    std::vector<std::set<Edge>::iterator> outEdges{};
    std::vector<std::set<Edge>::iterator> inEdges{};
};

class VCoords{
public:
    VCoords()=default;
    VCoords(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length){
        this->chrID=chrID;
        this->strand=strand;
        this->pos=pos;
        this->length=length;
    }
    ~VCoords()=default;

    bool operator< (const VCoords& vc) const{
        return (this->chrID < vc.chrID || this->strand < vc.strand || this->pos <  vc.pos);
    }
    bool operator> (const VCoords& vc) const{
        return (this->chrID > vc.chrID || this->strand > vc.strand || this->pos >  vc.pos);
    }
    bool operator==(const VCoords& vc) const{
        return (this->chrID == vc.chrID && this->strand == vc.strand && this->pos ==  vc.pos && this->length == vc.length);
    }

    int operator- (const VCoords& vc) const{
        if(this->chrID == vc.chrID && this->strand == vc.strand){
            return vc.pos - this->pos;
        }
        return -1;
    }

    struct coord_hash { // in case it is implemented as an unordered_map
        uint64_t operator()(const std::vector<std::pair<int,int> > &coords ) const
        {
            uint64_t resHash=1;
            for (auto const& coord: coords) {
                resHash=resHash^std::hash<uint64_t>()(coord.first)^std::hash<uint64_t>()(coord.second);
            }
            return resHash;
        }
    };

private:
    uint8_t chrID;
    uint8_t strand;
    uint32_t pos;
    uint8_t length;
};

class VMap{
public:
    VMap()=default;
    ~VMap()=default;

    std::map<VCoords,Vertex>::iterator _insert(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length){
        this->vm_it=this->vertices.insert(std::make_pair(VCoords(chrID,strand,pos,length),Vertex()));
        if(!this->vm_it.second){ // already exists
            this->vm_it.first->second.incWeight();
        }
        return this->vm_it.first;
    }
    void _addOutEdge(Vertex& vt){

    }
    void _addInEdge(Vertex& vt){

    }

    int getSize(){return static_cast<int>(this->vertices.size());}

private:
    std::map<VCoords,Vertex> vertices;
    std::pair<std::map<VCoords,Vertex>::iterator,bool> vm_it;
};

class HGraph {
public:
    HGraph();
    explicit HGraph(HDB* hdb);
    ~HGraph();

    void add_read(std::string& read);
    void print_stats();

    void to_sam(std::string out_sam_fname);
    void sort_graph();
    void parse_graph();

private:

    std::map<VCoords,Vertex>::iterator add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length);

    HDB* hdb;

    struct Stats{
        int numReads=0;
        int numReadsIgnored=0;
        int kmerlen;
    } stats;

    MinMap::iterator trans_it;
    HDB::GenMap::iterator genom_it;

    // below are graph declarations
    typedef std::tuple<uint8_t,uint8_t,uint32_t,uint8_t> vt_coords;
    VMap vertices;
    std::pair<std::map<VCoords,Vertex>::iterator,bool> vertex_it;

    std::vector<Edge> edges;
    std::pair<std::vector<Edge>::iterator,bool> edge_it;

};

#endif //HAIRPIN_GRAPH_H
