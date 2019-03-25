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
#include <set>

#include "HDB.h"
#include "MinMap.h"

class Vertex;

class Edge{
public:
    Edge()=default;

    Vertex* getOrigin() {return origin;}
    Vertex* getDestination() {return destination;}

private:
    Vertex* origin;
    Vertex* destination;
    bool known{}; // set to true if the edge is based on the transcriptome
    bool weight=0; // number of reads supporting the edge
};

class Vertex{
public:
    Vertex() = default;
    explicit Vertex(uint8_t chrID,uint8_t strand, uint32_t pos, uint8_t length){
        this->chrID=chrID;
        this->strand=strand;
        this->pos=pos;
        this->length=length;
        this->weight=1;
    }
    ~Vertex()=default;

    void addEdge(Vertex *v){

    }

    uint8_t _getChr(){return this->chrID;}
    uint8_t _getStrand(){return this->strand;}
    uint32_t _getPos(){return this->pos;}
    uint8_t _getLength(){return this->length;}
    int getWeight() {return this->weight;}

    bool operator< (const Vertex& vt) const{
        return (this->chrID < vt.chrID || this->strand < vt.strand || this->pos <  vt.pos);
    }
    bool operator> (const Vertex& vt) const{
        return (this->chrID > vt.chrID || this->strand > vt.strand || this->pos >  vt.pos);
    }
    bool operator==(const Vertex& vt) const{
        return (this->chrID == vt.chrID && this->strand == vt.strand && this->pos ==  vt.pos && this->length == vt.length);
    }

private:
    uint8_t chrID{};
    uint8_t strand{};
    uint32_t pos{};
    uint8_t length{};
    int weight=0;
    std::vector<Edge> edges{};
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

    std::set<Vertex>::iterator add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length);

    HDB* hdb;

    struct Stats{
        int numReads=0;
        int numReadsIgnored=0;
        int kmerlen;
        int numVertices=0;
    } stats;

    MinMap::iterator trans_it;
    HDB::GenMap::iterator genom_it;

    // below are graph declarations
    std::set<Vertex> vertices;
    std::pair<std::set<Vertex>::iterator,bool> vertex_it;
    std::set<Edge> edges;

};


#endif //HAIRPIN_GRAPH_H
