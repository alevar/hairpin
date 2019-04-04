//
// Created by Ales Varabyou on 3/16/19.
//

#include "HGraph.h"

HGraph::HGraph() = default;

HGraph::HGraph(HDB* hdb){
    this->hdb=hdb;
    this->stats.kmerlen=this->hdb->getKmerLen();
}

HGraph::HGraph(HDB* hdb,int max_intron,int min_intron,int min_mismatch,std::string out_fname){
    this->hdb=hdb;
    this->stats.kmerlen=this->hdb->getKmerLen();
    this->maxIntron=max_intron;
    this->minIntron=min_intron;
    this->out_fname=out_fname;
    this->minMismatch=min_mismatch;
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

int HGraph::coord_distance(const CoordVec& cv1,const CoordVec& cv2){ // compute the distance between two coordinate tuples
    if(std::get<0>(cv1) == std::get<0>(cv2) && std::get<1>(cv1) == std::get<1>(cv2)){
        return std::get<2>(cv1)-std::get<2>(cv2);
    }
    return -1;
}

void HGraph::add_read(std::string &read) {
    this->stats.numReads++;
    if(read.length()<this->stats.kmerlen){ // do not attempt to do anything with reads that are shorter than the kmer length chosen for the database
        this->stats.numReadsIgnored++;
        return;
    }
    std::string kmer;

    std::map<VCoords,Vertex>::iterator cur_vertex_it;

    std::map<std::map<VCoords,Vertex>::iterator,uint8_t,map_vit_cmp> extends_final; // holds final extension which should no longer be modified
    std::map<CoordVec,PosPair> extends; // groups of vertices that form stretches - will form edges between them afterwards; second element is the current length of the extendion

    for(int i=0;i<(read.length()-this->stats.kmerlen)+1;i++){
    kmer=read.substr(i, static_cast<unsigned long>(this->stats.kmerlen));

        this->genom_it=this->hdb->find_genom(kmer);
        if(this->genom_it!=this->hdb->genom_end()){ // match to genome found
            this->stats.numKmerMatchedGenom++;
            if(!extends.empty()){ // need to evaluate the predecessors, if a match did not occur near one of the positions

                std::vector<std::pair<CoordVec,PosPair> > case_one; // coordinates for which a distance of 1 was observed
                std::vector<std::pair<CoordVec,PosPair> >::iterator case_one_it;
                std::vector<std::pair<std::pair<CoordVec,PosPair>,std::pair<CoordVec,PosPair> > > case_two; // coordinates for which no distance of 1 was observed - thus needs to create a new entry - contains the closest coordinate vector found for here
                std::vector<std::pair<CoordVec,PosPair> > case_three; // coordinates for which no distance of 1 was observed and no closest match was observed either

                std::set<int> seen_idx; // which indices of the previous multimappers have been observed - is used to detect previous entries for which no extension was found
                for (auto &item : this->genom_it->second) { // for each multimapper

                    std::pair<CoordVec,PosPair> tmp_case_one; // temporary assignments
                    std::vector<std::pair<CoordVec,PosPair> > tmp_case_two;
                    bool tmp_case_one_flag=false,tmp_case_two_flag=false,tmp_case_three_flag=false; // flag set when created

                    for (auto cv : extends){
                        int dist=coord_distance(item,cv.first)-(cv.second.first-this->stats.kmerlen); // accounts for the length
                        if (dist==0){continue;}
                        else if (dist==1){ // do an extend
                            tmp_case_one=cv; // remember that this position has already been extended and needs not be evaluated in future multimappers
                            tmp_case_one_flag = true;
                            break; // great - found it - now can proceed to the next one
                        }
                        else{
                            if((dist-this->stats.kmerlen) >= this->minIntron && (dist-this->stats.kmerlen) <= this->maxIntron){
                                tmp_case_two_flag=true;
                                tmp_case_two.emplace_back(cv);
                            }
                            else{tmp_case_three_flag = true;} // case three evaluation - case three - no closest previous position observed
                        }
                    }
                    // preliminary case evaluation is performed here
                    if (tmp_case_one_flag){ // if for a given kmer match a 1-based distance was found - great - can clear cases two and three
                        case_one.emplace_back(tmp_case_one);
                        continue;
                    }
                    if (tmp_case_two_flag) { // no distance 1 was found - need to remember current multimapper value
                        for(auto cv: tmp_case_two) {
                            case_two.emplace_back(std::make_pair(std::make_pair(item, std::make_pair(this->stats.kmerlen,i)), cv));
                        }
                    }
                    if (tmp_case_three_flag){ case_three.emplace_back(std::make_pair(item,std::make_pair(this->stats.kmerlen,i)));}
                }
                // also need to know if any given previous entry has never been observed in the current round of evaluation
                // if such is the case - such entry needs to be pushed into the final_extends immediately and removed from the current search space

                // second case evaluation is performed here now that all current multimappers have been observed
                if(case_one.size()==this->genom_it->second.size()){ // all case ones
                    for (auto cv : case_one){
                        extends[cv.first].first++;
//                        extends[cv.first].second++; // update the position of the last kmer in the extend
                    }
                }
                else{ // add the case three entries
                    // first need to evaluate case two:
                    //    push any case_ones to the extends final if need be
                    //    and then increment them if they appear
                    for (auto cv : case_two){
                        case_one_it=std::find(case_one.begin(),case_one.end(),cv.second); // find the corresponding case_one - deal with it and remove from case ones
                        if (case_one_it != case_one.end()){ // found the case one entry
                            cur_vertex_it = this->add_vertex(std::get<0>(cv.second.first), std::get<1>(cv.second.first), std::get<2>(cv.second.first),(uint8_t)cv.second.second.first);
                            extends_final.insert(std::make_pair(cur_vertex_it,case_one_it->second.second)); // need to move it to extends final first
                            extends[case_one_it->first].first++; // increment the current length in case_one
//                            extends[case_one_it->first].second++; // update the position of the last kmer in the extend
                            extends.insert(cv.first); // add the new value to the extends
                            case_one.erase(case_one_it); // remove the entry from the case_one
                        }
                        else{
                            cur_vertex_it = this->add_vertex(std::get<0>(cv.second.first), std::get<1>(cv.second.first), std::get<2>(cv.second.first),(uint8_t)cv.second.second.first);
                            extends_final.insert(std::make_pair(cur_vertex_it,cv.second.second.second));
                            extends.insert(cv.first);
                        }
                    }
                    for (auto cv : case_one){
                        extends[cv.first].first++;
//                        extends[cv.first].second++; // update the position of the last kmer in the extend
                    } // now process the remaining entries in case one
//                    for (auto cv : case_three){extends.insert(cv);} // should this remove any from extends?
                }
            }
            else{ // nothing was there before, so need to add new elements based on the results of the database search
                for (auto &item : this->genom_it->second) { // for each multimapper
                    extends.insert(std::make_pair(item,std::make_pair(this->stats.kmerlen,i))); // initialize starting diagonals
                }
            }
        }
        else{this->stats.numKmerUnmatched++;}
    }
    // the read has been evaluated
    // now we can proceed to generate vertices and edges based on the final extensions
    //merge extends with the extends_final
    for (auto ex : extends){
        cur_vertex_it = this->add_vertex(std::get<0>(ex.first), std::get<1>(ex.first), std::get<2>(ex.first),(uint8_t)ex.second.first);
        extends_final.insert(std::make_pair(cur_vertex_it,ex.second.second));
    }
    std::map<std::map<VCoords,Vertex>::iterator,uint8_t>::iterator ci1,ci2;
    int counter=0;
    for (ci1=extends_final.begin();ci1!=std::prev(extends_final.end());ci1++){
        // create the first vertex
        for(ci2=std::next(ci1);ci2!=extends_final.end();ci2++){

            int dist=ci2->first->first - ci1->first->first - ci1->first->first.getLength();

            int kmerDist = (ci2->second - (ci1->second + (ci1->first->first.getLength()-this->stats.kmerlen)) ); // number of kmers separating the two vertices in the given read

            if (dist > 0 && dist <= this->maxIntron && // is not backward and is not greater than the maximum intron length; no minum threshold to permit edges to form for errors
                    kmerDist <= this->stats.kmerlen && kmerDist >= this->minMismatch){ // follows the expected number of missed kmers
                this->add_edge(ci1->first,ci2->first); // form an edge
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
    std::cout<<"Number of kmers matched to transcriptome: "<<this->stats.numKmerMatchedTrans<<std::endl;
    std::cout<<"Number of kmers matched to genome: "<<this->stats.numKmerMatchedGenom<<std::endl;
    std::cout<<"Number of kmers not matched to genome or transcriptome: "<<this->stats.numKmerUnmatched<<std::endl;
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

    int counter=0;
    std::string sub_seq;
    for(auto eit : this->emap){
        // first get fasta sequence
        this->getGenomeSubstr(eit.first,10,sub_seq);
        edges_fp << this->hdb->getContigFromID(eit.second.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                 << eit.second.getStart() << "\t" << eit.second.getEnd() << "\t"
                 << "." << "\t" << eit.second.getStrand() << "\t" << "." << "\t" <<  "weight="<<eit.second.getWeight()
                 <<";start="<< sub_seq.substr(0,20) <<";end=" << sub_seq.substr(sub_seq.length()-20,20) << std::endl;

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
    std::string sub_seq;
    for(auto eit : this->emap){
        // first get fasta sequence
        this->getGenomeSubstr(eit.first,this->stats.kmerlen-2,sub_seq);
        edges_fp << this->hdb->getContigFromID(eit.second.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                 << eit.second.getStart() << "\t" << eit.second.getEnd() << "\t"
                 << "." << "\t" << eit.second.getStrand() << "\t" << "." << "\t" <<  "weight="<<eit.second.getWeight()
                 <<";start="<< sub_seq.substr(0,this->stats.kmerlen) <<";end=" << sub_seq.substr(sub_seq.length()-(this->stats.kmerlen-2),this->stats.kmerlen) << std::endl;

        counter++;
    }

    edges_fp.close();
}

// this function takes the parsed graph and output the SAM-formatted file for further analysis
void HGraph::to_sam(std::string out_sam_fname) {

}

int HGraph::getGenomeSubstr(Edge et,int overhang, std::string& sub_seq){ // get substring from the genome based on the edges

    FastaReader fastaReader(this->hdb->getGenomeFname());
    FastaRecord cur_contig;

    while (fastaReader.good()) {
        fastaReader.next(cur_contig);

        uint8_t contigID = this->hdb->getIDFromContig(cur_contig.id_); //this is the id of the new contig

        if (contigID == et.first->first.getChr()){ // correct contig found
            int start = et.first->first.getEnd()-overhang;
            int length = (et.second->first.getStart()+overhang)-start;
            sub_seq=cur_contig.seq_.substr(start, length);

            if(et.first->first.getStrand()=='-'){ // need to reverse complement
                char *rev_sub_seq = new char[sub_seq.length()+1];
                strcpy(rev_sub_seq, sub_seq.c_str());
                reverseComplement(rev_sub_seq, length);
                sub_seq = rev_sub_seq;
                delete[] rev_sub_seq;
            }

            return 1;
        }
    }
    return 0;
}

int HGraph::getGenomeSubstr(uint8_t chrID,int8_t strand, uint32_t start, uint8_t length, std::string& sub_seq){
    return 0;
}