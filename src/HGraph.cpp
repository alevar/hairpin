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
    std::vector<Vertex> cur_matches; // each element in the current vector contains a single one of the multimappers

    //TODO: how to deal with multimappers? need multiple starts to be evaluated
    for(int i=0;i<read.length()-this->stats.kmerlen;i++){
        kmer=read.substr(i,this->stats.kmerlen);
        // lookup in the transcriptome database
        this->trans_it=this->hdb->find_trans(kmer);
        if(this->trans_it!=this->hdb->trans_end()){ // match to transcriptome found
//
        }
        else{ // search in the genom map
            this->genom_it=this->hdb->find_genom(kmer);
            if(this->genom_it!=this->hdb->genom_end()){ // match to genome found
                if(!cur_matches.empty()){

                }
                else{ // nothing was there before, so need to add new elements based on the results of the database search
                    for (auto &item : this->genom_it->second) { // for each multimapper
                        this->add_vertex(std::get<0>(item),std::get<1>(item),std::get<2>(item),(uint8_t)this->hdb->getKmerLen());
                    }
                }
            }
        }
    }
}

// print general statistics about the quantification process
void HGraph::print_stats() {
    std::cout<<"Total number of reads: "<<this->stats.numReads<<std::endl;
    std::cout<<"Number of reads ignored: " <<this->stats.numReadsIgnored<<std::endl;
    std::cout<<"Number of vertices: "<<this->stats.numVertices<<std::endl;
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

std::set<Vertex>::iterator HGraph::add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length) {
    this->vertex_it = this->vertices.insert(Vertex(chrID,strand,pos,length));
    if(this->vertex_it.second){ // new vertex inserted
        this->stats.numVertices++;
    }
    return this->vertex_it.first;
}

