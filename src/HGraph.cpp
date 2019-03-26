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
    std::vector<std::map<VCoords,Vertex>::iterator> cur_matches; // each element in the current vector contains a single one of the multimappers
    std::vector<std::map<VCoords,Vertex>::iterator> prev_matches; // vertices inserted from the previous kmer evaluation

    std::map<VCoords,Vertex>::iterator cur_vertex_it;

    for(int i=0;i<read.length()-this->stats.kmerlen;i++){
        kmer=read.substr(i,this->stats.kmerlen);
        // lookup in the transcriptome database
        this->trans_it=this->hdb->find_trans(kmer);
        if(this->trans_it!=this->hdb->trans_end()){ // match to transcriptome found

        }
        else{ // search in the genom map
            this->genom_it=this->hdb->find_genom(kmer);
            if(this->genom_it!=this->hdb->genom_end()){ // match to genome found
                if(!cur_matches.empty()){ // need to evaluate the predecessors, if a match did not occur near one of the positions
                    // subtract the new vertex from the other and see if the difference is 1
                    // if not 1 - this is a candidate for an edge


                    for (auto &item : this->genom_it->second) { // for each multimapper
                        cur_vertex_it=this->add_vertex(std::get<0>(item),std::get<1>(item),std::get<2>(item),(uint8_t)this->stats.kmerlen);
                        for (int prev_vi=0;prev_vi<cur_matches.size();prev_vi++){ // see if a match is found in the current list of vertices
                            if(cur_vertex_it==cur_matches[prev_vi]){ //points to the same node
                                continue;
                            }
                            int dist=cur_vertex_it->first-cur_matches[prev_vi]->first;
                            if(dist == 1){ //only one base difference - can replace the entry with the new
                                cur_matches[prev_vi]=cur_vertex_it;
                                break; // found match can evaluate the next multimapper
                            }
                            else if(dist > 1){ // means no close match was found and this is a candidate for an edge
//                              this->add_edge(cur_matches[prev_vi],cur_vertex_it);
                                // add edges
                                int edge_inc = cur_matches[prev_vi]->second.addOutEdge(cur_vertex_it);
                                this->stats.numEdges = this->stats.numEdges + edge_inc; // increment the number of edges
                                cur_vertex_it->second.addInEdge(cur_matches[prev_vi]);
                                // also add respective edges to the vertices using addOutEdge and addInEdge
                                cur_matches[prev_vi]=cur_vertex_it;
                            }
                            else{ // anything else than one should be ignored and inserted back into the cur_matches;
//                                cur_matches.emplace_back(cur_vertex_it);
                            }
                        }
                    }
                }
                else{ // nothing was there before, so need to add new elements based on the results of the database search
                    for (auto &item : this->genom_it->second) { // for each multimapper
                        cur_vertex_it=this->add_vertex(std::get<0>(item),std::get<1>(item),std::get<2>(item),(uint8_t)this->stats.kmerlen);
                        cur_matches.emplace_back(cur_vertex_it);
                    }
                }
            }
        }
    }
}

// print general statistics about the quantification process
void HGraph::print_stats() {
    std::cout<<"Number of reads: "<<this->stats.numReads<<std::endl;
    std::cout<<"Number of reads ignored: " <<this->stats.numReadsIgnored<<std::endl;
    std::cout<<"Number of vertices: "<<this->vertices.getSize()<<std::endl;
    std::cout<<"Number of edges: "<<this->stats.numEdges<<std::endl;
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

std::map<VCoords,Vertex>::iterator HGraph::add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length) {
    return this->vertices._insert(chrID,strand,pos,length);
}

