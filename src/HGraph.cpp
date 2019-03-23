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
    std::vector<std::vector<std::pair<uint32_t,uint8_t>>> cur_matches; // each element in the current vector contains a single one of the multimappers

    //TODO: how to deal with multimappers? need multiple starts to be evaluated
    for(int i=0;i<read.length()-this->stats.kmerlen;i++){
        kmer=read.substr(i,this->stats.kmerlen);
        // lookup in the transcriptome database
        this->trans_it=this->hdb->find_trans(kmer);
        if(this->trans_it!=this->hdb->trans_end()){ // match to transcriptome found
//            std::cout<<kmer<<"\t"<<this->trans_it->first<<"\t"<<this->trans_it->second[0]._getStart(0)<<std::endl;

        }
        else{ // search in the genom map
            this->genom_it=this->hdb->find_genom(kmer);
            if(this->genom_it!=this->hdb->genom_end()){ // match to genome found

            }
            std::cout<<"searched in the genome: "<<kmer<<std::endl;
        }
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
