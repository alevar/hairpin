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
    explicit Vertex(EVec * coords){
        this->coordinates=coords;
        this->weight=1;
    }
    ~Vertex()=default;

    void addEdge(Vertex *v){

    }

    // TODO: comparison operator so that verticies can be stored in a set
    //

    int getWeight() {return this->weight;}

private:
    EVec* coordinates{};
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

    HDB* hdb;

    struct Stats{
        int numReads=0;
        int numReadsIgnored=0;
        int kmerlen;
    } stats;

    MinMap::iterator trans_it;
    HDB::GenMap::iterator genom_it;

    // below are graph declarations
    std::set<Vertex> verices;
    std::set<Edge> edges;

};


#endif //HAIRPIN_GRAPH_H
