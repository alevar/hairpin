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
                        int secondFound=-1; // the index of the first element to which an edge is added
                        std::vector<int> seconds; // second cases where no match was found but an edge might exist
                        bool third=false; // third case - should the new entry be created?

                        for (int prev_vi=0;prev_vi<cur_matches.size();prev_vi++){ // see if a match is found in the current list of vertices
                            if(cur_vertex_it==cur_matches[prev_vi]){ //points to the same node
                                continue;
                            }
                            int dist=cur_vertex_it->first-cur_matches[prev_vi]->first;
                            if (dist==0){
                                continue; //same as the previous comparison
                            }
                            else if(dist == 1){ //only one base difference - can replace the entry with the new
                                cur_matches[prev_vi]=cur_vertex_it;
                                first=true;
                                break; // found match can evaluate the next multimapper
                            }
                            else if((dist-this->stats.kmerlen) > this->minIntron && (dist-this->stats.kmerlen) < this->maxIntron){ // means no close match was found and this is a candidate for an edge
                                if(secondFound==-1){
                                    secondFound=prev_vi;
                                }
                                else{
                                    seconds.emplace_back(prev_vi);
                                }
                            }
                            else{ // anything else than one should be ignored and inserted back into the cur_matches;
                                // add as a candidate for a new node
                                third=true;
                            }
                        }
                        if(!first){ // did not find the best case and need to add cases two and three
                            if(secondFound>=0){
                                // process original second
                                int edge_inc = cur_matches[secondFound]->second.addOutEdge(cur_vertex_it);
                                this->stats.numEdges = this->stats.numEdges + edge_inc; // increment the number of edges
                                cur_vertex_it->second.addInEdge(cur_matches[secondFound]);
                                cur_matches[secondFound]=cur_vertex_it;

                                for (int second_idx : seconds) { // first output edges
                                    // add as a candidate for a potential splice site
                                    // add edges
                                    edge_inc = cur_matches[second_idx]->second.addOutEdge(cur_vertex_it);
                                    this->stats.numEdges = this->stats.numEdges + edge_inc; // increment the number of edges
                                    cur_vertex_it->second.addInEdge(cur_matches[second_idx]);
                                }
                                for (int second_idx : seconds) { // remove redundant entries
                                    cur_matches.erase(cur_matches.begin()+second_idx);
                                }
                            }
                            if(!third && secondFound == -1){
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

void HGraph::write_intron_gff() {
    std::string edges_fname(this->out_fname);
    edges_fname.append(".intron.gff");
    std::ofstream edges_fp(edges_fname.c_str());

    int counter=0; // sets id
    for(auto vit: this->vertices){
        if(vit.second.getOutDegree()>0) {
            for(auto eit: vit.second.getOutEdges()) {
                counter++;
                edges_fp << this->hdb->getContigFromID(vit.first.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                         << vit.first.getStart() + vit.first.getLength() << "\t" << eit.first->first.getStart() << "\t"
                         << "." << "\t" << vit.first.getStrand() << "\t" << "." << "\t" <<  eit.second.getWeight() << std::endl;
            }
        }
    }

    edges_fp.close();
}

// parse graph and evaluate gaps and assign splice junctions and mismatches and gaps
void HGraph::parse_graph() {
    // this is a toy example for prsing the data
    std::string edges_fname(this->out_fname);
    edges_fname.append(".edges");
    std::ofstream edges_fp(edges_fname.c_str());

    int counter=0; // sets id
    for(auto vit: this->vertices){
        if(vit.second.getOutDegree()>0) {
            for(auto eit: vit.second.getOutEdges()) {
                counter++;
                edges_fp << this->hdb->getContigFromID(vit.first.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                         << vit.first.getStart() + vit.first.getLength() << "\t" << eit.first->first.getStart() << "\t"
                         << "." << "\t" << vit.first.getStrand() << "\t" << "." << "\t" <<  eit.second.getWeight() << std::endl;
            }
        }
    }

    edges_fp.close();
}

// this function takes the parsed graph and output the SAM-formatted file for further analysis
void HGraph::to_sam(std::string out_sam_fname) {

}
