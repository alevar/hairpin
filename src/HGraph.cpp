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

void HGraph::add_edge(std::map<VCoords,Vertex>::iterator prev,std::map<VCoords,Vertex>::iterator next){

    int edge_inc = prev->second.addOutEdge(next);
    next->second.addInEdge(prev);
    this->stats.numEdges = this->stats.numEdges + edge_inc; // increment the number of edges

    // now add to the main cluster
    this->emap_it=this->emap.insert(std::make_pair(std::make_pair(prev,next),Aggregate_edge_props(prev,next)));
}

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
        kmer=read.substr(i, static_cast<unsigned long>(this->stats.kmerlen));
        // lookup in the transcriptome database
        this->trans_it=this->hdb->find_trans(kmer);
        if(this->trans_it!=this->hdb->trans_end()){ // match to transcriptome found
            // TODO: do the transcriptomic search now
            if(!cur_matches.empty()) { // TODO: should the same cur_matches be used for the transcriptomic and genomic? At the moment it seems to be the case since there could be cases when a transcriptomic match follows a genomic and a potential splice junction can exist
//                for (auto &ev : this->trans_it->second) { // for each multimapper in the transcriptomic hit
//                    for (auto &ep : ev){ // since each transcriptomic hit can have several coordinates spanning a splice junction we need to add a Vertex for each chunk within the coordinate vector
//                        cur_vertex_it=this->add_vertex(ev._getChr(),ev._getStrand(),ep.getStart(),ep.getLength());
//
//                        bool first=false; // most parsimonious case found - 1 base difference - should skip other cases
//                        int secondFound=-1; // the index of the first element to which an edge is added
//                        std::vector<int> seconds; // second cases where no match was found but an edge might exist
//                        bool third=false; // third case - should the new entry be created?
//
//                        // the logic here will be a little more involved
//
//                        // we are operating on chunks of a kmer whenever there is a splice junction
//
//                        // 1. if this is the first chunk in the coordinate vector - do the same evaluation as before
//                        // 2. if the chunk is not the first - then need to figure out which multimapper it is associated with
//                        //      this need to be handled correctly since the logic is likely to be different from that of the genomic lookup
//
//                        // here a different type of lookup is possible
//                        // if a kmer spans a junction - and the first vertex has a distance of 1
//                        // then the junction spanned should also be compatible with the previously observed kmer
//
//                        int counter=0;
//                        for (int prev_vi=0;prev_vi<cur_matches.size();prev_vi++){
//                            if(cur_vertex_it==cur_matches[prev_vi]){continue;} // first case
//
//                            int dist=cur_vertex_it->first-cur_matches[prev_vi]->first;
//                            if (dist==0){continue;} //same as first case
//
//                        }
//                    }
//                }
            }
            else{ // nothing was there before, so need to add new elements based on the results of the database search
                for (auto &ev : this->trans_it->second) { // for each multimapper
                    std::map<VCoords,Vertex>::iterator prev_vertex,first_vertex;
                    int counter=0;
                    for (auto &ep : ev){ // for each pair in the coordinate vector
                        cur_vertex_it=this->add_vertex(ev._getChr(),ev._getStrand(),ep.getStart(),ep.getLength());
                        if(counter>0){
                            if(prev_vertex->first.getEnd()==210675&& cur_vertex_it->first.getStart()==211590){
                                std::cout<<"found add_edge"<<std::endl;
                            }
                            this->add_edge(prev_vertex,cur_vertex_it);
                        } // if the current piece indicates a splice junction - need to insert an edge - add edge
                        else{first_vertex=cur_vertex_it;} // remeber the first vertex
                        counter++;
                        prev_vertex=cur_vertex_it;
                    }
                    cur_matches.emplace_back(first_vertex); // only remember the position of the first vertex for the next comparrisons to the cur_matches
                }
            }
        }
        else{ // search in the genom map
//            std::cerr<<"searching in genome"<<std::endl;
            this->genom_it=this->hdb->find_genom(kmer);
            if(this->genom_it!=this->hdb->genom_end()){ // match to genome found
                if(!cur_matches.empty()){ // need to evaluate the predecessors, if a match did not occur near one of the positions
                    for (auto &item : this->genom_it->second) { // for each multimapper
                        cur_vertex_it=this->add_vertex(std::get<0>(item),std::get<1>(item),std::get<2>(item),(uint8_t)this->stats.kmerlen);
                        // there should be a priority when deciding which one to choose
                        // 1 base difference has a priority

                        bool first=false; // most parsimonious case found - 1 base difference - should skip other cases
                        int secondFound=-1; // the index of the first element to which an edge is added
                        std::vector<int> seconds; // second cases where no match was found but an edge might exist
                        bool third=false; // third case - should the new entry be created?

                        for (int prev_vi=0;prev_vi<cur_matches.size();prev_vi++){ // see if a match is found in the current list of vertices
                            if(cur_vertex_it==cur_matches[prev_vi]){continue;} //points to the same node

                            int dist=cur_vertex_it->first-cur_matches[prev_vi]->first;
                            if (dist==0){continue;} //same as the previous comparison

                            else if(dist == 1){ //only one base difference - can replace the entry with the new
                                cur_matches[prev_vi]=cur_vertex_it;
                                first=true;
                                break; // found match can evaluate the next multimapper
                            }
                            else if((dist-this->stats.kmerlen) > this->minIntron && (dist-this->stats.kmerlen) < this->maxIntron){ // means no close match was found and this is a candidate for an edge
                                if(secondFound==-1){secondFound=prev_vi;}
                                else{seconds.emplace_back(prev_vi);}
                            }
                            else{third=true;} // anything else than one should be ignored and inserted back into the cur_matches;
                        }
                        if(!first){ // did not find the best case and need to add cases two and three
                            if(secondFound>=0){
                                // process original second
                                this->add_edge(cur_matches[secondFound],cur_vertex_it);
                                cur_matches[secondFound]=cur_vertex_it;

                                for (int second_idx : seconds) {this->add_edge(cur_matches[second_idx],cur_vertex_it);} // first output edges - add as a candidate for a potential splice site
                                for (int second_idx : seconds) {cur_matches.erase(cur_matches.begin()+second_idx);}// remove redundant entries
                            }
                            if(!third && secondFound == -1){cur_matches.emplace_back(cur_vertex_it);}
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
    std::cout<<"Number of edges: "<<this->emap.size()<<std::endl;
}

std::map<VCoords,Vertex>::iterator HGraph::add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length) {
    return this->vertices._insert(chrID,strand,pos,length);
}

// TODO: need to implement add_vertex for transcriptome - in that case it should somehow unite vertices based on the splice junction
//      since no longer vertices have the same length - they now need to be able to form the same vertex even if of different lengths and with different start sites
//      potentially we can add them either based on the start or based on the end if splice junction, in which case the vertex can be extended in either direction
//      however, the problem with that approach is that we can not modify vertex coordinates since they are of const type as the map key

// topological sort of the graph
void HGraph::sort_graph() {

}

void HGraph::write_intron_gff() {
    std::string edges_fname(this->out_fname);
    edges_fname.append(".intron.gff");
    std::ofstream edges_fp(edges_fname.c_str());

    int counter=0;
    for(auto eit : this->emap){
        edges_fp << this->hdb->getContigFromID(eit.second.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                 << eit.second.getStart() << "\t" << eit.second.getEnd() << "\t"
                 << "." << "\t" << eit.second.getStrand() << "\t" << "." << "\t" <<  eit.second.getWeight() << std::endl;

        counter++;
    }

    edges_fp.close();
}

// parse graph and evaluate gaps and assign splice junctions and mismatches and gaps
void HGraph::parse_graph() {
    // this is a toy example for prsing the data
    std::string edges_fname(this->out_fname);
    edges_fname.append(".intron.gff");
    std::ofstream edges_fp(edges_fname.c_str());

    int counter=0;
    for(auto eit : this->emap){ // cycles through all the edges in ascending order
        edges_fp << this->hdb->getContigFromID(eit.second.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                 << eit.second.getStart() << "\t" << eit.second.getEnd() << "\t"
                 << "." << "\t" << eit.second.getStrand() << "\t" << "." << "\t" <<  eit.second.getWeight() << std::endl;

        counter++;
    }

    edges_fp.close();
}

// this function takes the parsed graph and output the SAM-formatted file for further analysis
void HGraph::to_sam(std::string out_sam_fname) {

}
