//
// Created by Ales Varabyou on 3/16/19.
//

#include "HGraph.h"

HGraph::HGraph() = default;

HGraph::HGraph(HDB* hdb){
    this->hdb=hdb;
    this->stats.kmerlen=this->hdb->getKmerLen();
}

HGraph::~HGraph() = default;

void HGraph::add_read(std::string &read) {
    this->stats.numReads++;
    if(read.length()<this->stats.kmerlen){ // do not attempt to do anything with reads that are shorter than the kmer length chosen for the database
        this->stats.numReadsIgnored++;
        return;
    }
    std::string kmer;
    int start; // start of the first kmer
    int length; // length of the current match
    //TODO: how to deal with multimappers? need multiple starts to be evaluated
    for(int i=0;i<read.length()-this->stats.kmerlen;i++){
        kmer=read.substr(i,this->stats.kmerlen);
        // for each kmer position output create a new node and keep track of all the nodes created by the current read
        // increment the length for each extension on each vertex
        // create edges respectively for each if no match is found
    }
}

// print general statistics about the quantification process
void HGraph::print_stats() {
    std::cout<<"Total number of reads"<<stats.numReads<<std::endl;
    std::cout<<"Number of reads ignored:" <<stats.numReadsIgnored<<std::endl;
}

// this function takes the parsed graph and output the SAM-formatted file for further analysis
void HGraph::to_sam(std::string out_sam_fname) {

}

// topological sort of the graph
void HGraph::sort_graph() {

}

// parse graph and evaluate gaps and assign splice junctions and mismatches and gaps
void HGraph::parse_graph() {

}
