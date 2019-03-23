//
// Created by Ales Varabyou on 3/16/19.
// ====================================
// this class contains the implementation
// of a graph populated with kmers from sequenced reads
//

#ifndef HAIRPIN_GRAPH_H
#define HAIRPIN_GRAPH_H

#include <iostream>
#include <boost/graph/adjacency_list.hpp>

#include "HDB.h"
#include "MinMap.h"

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
    struct Vertex{
        uint32_t start; // start position of the current vertex
        uint8_t length; // number of bases covered by the current node
        uint8_t chrID;
        uint8_t strand;
    };
    struct Edge{
        uint16_t weight; // might have to be uint32_t
        bool known; //whether the edge was inesrted based on the transcriptomic match
    };

    //Define the graph using those classes
    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, Vertex, Edge > Graph;
    //Some typedefs for simplicity
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

    Graph graph;

    // Example code for adding vertices and edges to the graph
//    Vertex u = boost::add_vertex(g);
//    Vertex v = boost::add_vertex(g);
//    // Create an edge conecting those two vertices
//    Edge e; bool b;
//    boost::tie(e,b) = boost::add_edge(u,v,g);
//    // Set the properties of a vertex and the edge
//    g[u].foo = 42;
//    g[e].blah = "Hello world";

    HDB* hdb;

    struct Stats{
        int numReads=0;
        int numReadsIgnored=0;
        int kmerlen;
    } stats;

    MinMap::iterator trans_it;
    HDB::GenMap::iterator genom_it;

};


#endif //HAIRPIN_GRAPH_H
