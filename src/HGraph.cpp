//
// Created by Ales Varabyou on 3/16/19.
//

#include "HGraph.h"

HGraph::HGraph() = default;

HGraph::HGraph(HDB* hdb){
    this->hdb=hdb;
    this->stats.kmerlen=this->hdb->getKmerLen();
}

HGraph::HGraph(HDB* hdb,int max_intron,int min_intron, std::string out_fname){
    this->hdb=hdb;
    this->stats.kmerlen=this->hdb->getKmerLen();
    this->maxIntron=max_intron;
    this->minIntron=min_intron;
    this->out_fname=out_fname;
}

HGraph::~HGraph() = default;

// TODO: need to account for reverse strand reads - for now need to make sure that simulations are with respect to the same strand

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
                    for (auto &item : this->genom_it->second) { // for each multimapper
                        cur_vertex_it=this->add_vertex(std::get<0>(item),std::get<1>(item),std::get<2>(item),(uint8_t)this->stats.kmerlen);
                        // there should be a priority when deciding which one to choose
                        // 1 base difference has a priority

                        // subtract the new vertex from the other and see if the difference is 1
                        // if not 1 - this is a candidate for an edge
                        bool first=false; // most parsimonious case found - 1 base difference - should skip other cases
                        std::vector<int> second; // second cases where no match was found but an edge might exist
                        std::vector<int> third; // third case - might not even be relevant, but we should still include these coordinates just in case

                        for (int prev_vi=0;prev_vi<cur_matches.size();prev_vi++){ // see if a match is found in the current list of vertices
                            if(cur_vertex_it==cur_matches[prev_vi]){ //points to the same node
                                continue;
                            }
                            int dist=cur_vertex_it->first-cur_matches[prev_vi]->first;
                            if(dist == 1){ //only one base difference - can replace the entry with the new
                                cur_matches[prev_vi]=cur_vertex_it;
                                first=true;
                                break; // found match can evaluate the next multimapper
                            }
                            else if(dist > this->minIntron && dist < this->maxIntron){ // means no close match was found and this is a candidate for an edge
                                second.emplace_back(prev_vi);
                            }
                            else{ // anything else than one should be ignored and inserted back into the cur_matches;
                                // add as a candidate for a new node
                                third.emplace_back(prev_vi);
                            }
                        }
                        if(!first){ // did not find the best case and need to add cases two and three
                            for(int si=0;si<second.size();si++){ // first output edges
                                // add as a candidate for a potential splice site
                                // add edges
                                int edge_inc = cur_matches[si]->second.addOutEdge(cur_vertex_it);
                                this->stats.numEdges = this->stats.numEdges + edge_inc; // increment the number of edges
                                cur_vertex_it->second.addInEdge(cur_matches[si]);
                                // replace old entry with the new one
                                cur_matches[si]=cur_vertex_it;
                            }
                            if(!third.empty()){
                                cur_matches.emplace_back(cur_vertex_it);
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

std::map<VCoords,Vertex>::iterator HGraph::add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length) {
    return this->vertices._insert(chrID,strand,pos,length);
}

// topological sort of the graph
void HGraph::sort_graph() {

}

// parse graph and evaluate gaps and assign splice junctions and mismatches and gaps
void HGraph::parse_graph() {
    // this is a toy example for prsing the data
    std::string edges_fname(this->out_fname);
    edges_fname.append(".edges");
    std::ofstream edges_fp(edges_fname.c_str());

    int counter=0; // sets id
    for(auto vit: this->vertices){
        std::cout<<"hello\t"<<vit.second.getOutDegree()<<std::endl;
        if(vit.second.getOutDegree()>0) {
            std::cout<<"hello2"<<std::endl;
            for(auto eit: vit.second.getEdges()) {
                counter++;
                edges_fp << counter << "\t" << vit.first.getChr() << "\t" << vit.first.getStrand() << "\t"
                         << vit.first.getStart() + vit.first.getLength() << "\t" << eit.first->first.getStart()+eit.first->first.getLength() << "\t"
                         << eit.second.getWeight() << "\t" << std::endl;
            }
        }
    }

    edges_fp.close();
}

// this function takes the parsed graph and output the SAM-formatted file for further analysis
void HGraph::to_sam(std::string out_sam_fname) {

}
