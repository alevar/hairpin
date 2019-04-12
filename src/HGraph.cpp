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

//            std::cout<<kmer<<std::endl;

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

//                            std::cout<<"\t1: "<<kmer<<std::endl;

                        }
                        else{
                            cur_vertex_it = this->add_vertex(std::get<0>(cv.second.first), std::get<1>(cv.second.first), std::get<2>(cv.second.first),(uint8_t)cv.second.second.first);
                            extends_final.insert(std::make_pair(cur_vertex_it,cv.second.second.second));

//                            std::cout<<"\t2: "<<kmer<<std::endl;

//                            if (cur_vertex_it->first.getEnd() == 70567){
//                                std::cout<<read<<std::endl;
//                            }

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
        else{
            this->stats.numKmerUnmatched++;
//            std::cout<<"  "<<kmer<<std::endl;
        }
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
                    kmerDist <= this->stats.kmerlen && kmerDist >= this->stats.kmerlen/2){ // follows the expected number of missed kmers; the division by two here is due to the fact that an overhang can not be in the second half, since otherwise there would not be enough space on the receiving end
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
    for(const auto& eit : this->emap){
        // first get fasta sequence
        this->getGenomeSubstr(eit.first,this->stats.kmerlen-2,sub_seq);
        edges_fp << this->hdb->getContigFromID(eit.second.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                 << eit.second.getStart() << "\t" << eit.second.getEnd() << "\t"
                 << "." << "\t" << eit.second.getStrand() << "\t" << "." << "\t" <<  "weight="<<eit.second.getWeight()
                 <<";start="<< sub_seq.substr(0,this->stats.kmerlen) <<";end=" << sub_seq.substr(sub_seq.length()-(this->stats.kmerlen),this->stats.kmerlen) << std::endl;

        counter++;
    }

    edges_fp.close();
}

//TODO: Ruleset for evaluation of potential splice junctions
//        1. INCORRECT DUE TO ALTERNATIVe SPLICING: Each vertex may only be allowed to have a single outgoing or ingoing edge
//        2. INCORRECT DUE TO ALTERNATIVe SPLICING: There may not be parallel edges - that is no two edges may overlap by coordinates
//        3. No two edges may exist with less than k bases between the end of the first edge and the start of the next edge
//        4. DONE: donor/acceptor sites are weighted (canonical/semi-canonical/non-canonical)
//        5. Vertices at each end of the junction shall have identical weights corresponding to the evaluated edge
//        6. IN PROGRESS: Number of distinct start positions of reads that cover a splice junction can also be taken into account when computing the likelihood score
//        7. IN PROGRESS: We can also create a threshold for the number of unique start sites preceding a splice junction
//        8. DONE: If two contradictory junctions exist with otherwise identical scores - the shorter one shall be preferred
//        9. DONE: the number of kmers that span a vertex-edge-vertex should be no less than readlength-(kmerlen*2)
//
//
//TODO: Output of the edge parser:
//  The edge iterator is a pointer to a <key,value> pair from the edge map, where the value is of the type Aggregate_edge_props and contains a function validate();
//  Whenever the parser decides that the splice junction is valid it should use this function to set the known flag to true
//
//TODO: IN PROGRESS: The parser has to make sure that the total number of bases contained in connected vertices is above a threshold
//      this is important to ensure spurious/incomplete multimappers are not present in the final graph
//
//
//TODO: All the criteria in the ruleset must be normalized to a 0-1 scale
//   for instance the weight of canonical vs semi vs noncanonical splice sites is already enforced to be in the correct range
//   length of the junction can also be waited, where if many potential precise splice junctions exist - the shortest one (above the minimum
//      intron length threshold is given the precedence by asigning the highest value and the rest are weighted with respect to this shortest length
//   afterwards, when all the individual scores are computed, we can put them through a series of sigmoid function with cutoffs
//      and then be able to chose the optimal alignment

uint8_t HGraph::getEdgeChr(const std::pair<Edge,Aggregate_edge_props>& eit){
    return eit.second.getChr();
}

uint8_t HGraph::getEdgeStrand(const std::pair<HGraph::Edge, Aggregate_edge_props> &eit) {
    return eit.second.getStrand();
}

uint8_t HGraph::getEdgeStart(const std::pair<HGraph::Edge, Aggregate_edge_props> &eit) {
    return eit.second.getStart();
}

uint8_t HGraph::getEdgeEnd(const std::pair<HGraph::Edge, Aggregate_edge_props> &eit) {
    return eit.second.getEnd();
}

// this function identifies all possible precise splice junctions given one potential splice junction (edge)
// the evaluation is based entirely on the donor/acceptor pairs
void HGraph::evaluate_sj(const std::pair<Edge,Aggregate_edge_props>& eit,const std::pair<std::string,double>& donor,const std::pair<std::string,double>& acceptor, SJS& sj_map){
    std::string sub_seq, start_seq, end_seq;
    this->getGenomeSubstr(eit.first,this->stats.kmerlen-2,sub_seq); // first get fasta sequence
    start_seq = sub_seq.substr(0,this->stats.kmerlen); // get sequence for the potential start of the splice junction
    std::transform(start_seq.begin(), start_seq.end(), start_seq.begin(), ::toupper);
    end_seq = sub_seq.substr(sub_seq.length()-(this->stats.kmerlen),this->stats.kmerlen); // get sequence for the potential end of the splice junction
    std::transform(end_seq.begin(), end_seq.end(), end_seq.begin(), ::toupper);

//    std::cout<<start_seq<<std::endl;
//    std::cout<<end_seq<<std::endl;

    std::map<int,std::string> starts,ends; // holds positions and other information from donor and acceptor sites

    size_t pos = start_seq.find(donor.first, 0); // find canonical donor on the donor site
    while(pos != std::string::npos){
        starts.insert(std::make_pair(pos,start_seq.substr(pos,start_seq.length()-pos-2)));
        pos = start_seq.find(donor.first,pos+1);
    }

    pos = end_seq.find(acceptor.first, 0); // find canonical acceptor on the acceptor site
    while(pos != std::string::npos){
        ends.insert(std::make_pair(pos,end_seq.substr(pos+2,end_seq.length()-pos)));
        pos = end_seq.find(acceptor.first,pos+1);
    }

    for (const auto& start_it : starts){ // now figure out which one is the true junction
        for (const auto& end_it : ends){

            int numBases_after_donor=(start_seq.length()-2)-start_it.first;
            int numBases_before_acceptor=end_it.first;

            if (start_it.first - numBases_before_acceptor >= 0 &&
                end_it.first + numBases_before_acceptor <= this->stats.kmerlen){ // possible now compare the actual strings

                std::string donor_overhang=start_seq.substr(start_it.first,numBases_after_donor);
                std::string acceptor_overhang=end_seq.substr((end_it.first+2)-numBases_before_acceptor,numBases_before_acceptor);

                std::string donor_overhang_on_acceptor=end_seq.substr(end_it.first+2,numBases_after_donor);
                std::string acceptor_overhang_on_donor=start_seq.substr(start_it.first-numBases_before_acceptor,numBases_before_acceptor);

                if (donor_overhang_on_acceptor == donor_overhang && acceptor_overhang_on_donor == acceptor_overhang){
                    sj_map.insert(std::make_pair(std::make_tuple(getEdgeChr(eit),
                                                                 getEdgeStrand(eit),
                                                                 eit.second.getStart() - numBases_after_donor,
                                                                 eit.second.getEnd()+1 + numBases_before_acceptor),
                                                 std::make_tuple(donor.second,acceptor.second)));
                }
            }
        }
    }
}

// this function runs splice junction verification for all donor/acceptor pairs
void HGraph::evaluate_donor_acceptor(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm){
    for (const auto& dit : this->donors) { // iterate over all the donors
        for (const auto& ait : this->acceptors) { // iterate over all acceptors
            evaluate_sj(eit,dit,ait,sm);
        }
    }
}

// given all precise splice junctions this function enforces any specified constraints
// currently supports only the minimum and the maximum intron lengths
// but can be extended to deal with more if needed
void HGraph::enforce_constraints(SJS& sm){
    int length = 0;
    for (auto sm_it = sm.cbegin(); sm_it != sm.cend();){
        length = std::get<3>(sm_it->first) - std::get<2>(sm_it->first);
        if (length < this->minIntron || length > this->maxIntron){
            sm_it = sm.erase(sm_it);
        }
        else{
            ++sm_it;
        }
    }
}

// This function computes the number of bases in some exon. taking into account implicit edges (overlaps in vertex coordinates)
// computes untill no vertex overlaps are found
int HGraph::compute_maximal_exon_length(SJS& sm){
    return 1;
}

// This function computes the number of bases in some exon, taking into account implicit edges (overlaps n vertex coordinates)
// only comptes till the first validated edge
int HGraph::compute_minimal_exon_length(SJS& sm){
    return 1;
}

// This function computes the number of bases that belong to a given junction
// this takes into account implicit edges - overlaps between vertices
int HGraph::compute_maximal_clique_length(SJS& sm){
    return 1;
}

// This function computes the total number of bases in the vertices directly connected to the given edge
int HGraph::compute_minimal_clique_length(SJS& sm){
    return 1;
}

// This is a recursive implementation which return the total number of bases in the vertices connected by edges to the current edge
// no implicit edge evaluation is done here
// previous vertices only
int HGraph::_get_read_length_before(std::map<VCoords,Vertex>::iterator vit){
    int cur_length = vit->first.getLength();

    if (vit->second.getInDegree()>0){ // need to follow the exon structure for reads that span multiple introns
        // TODO: this is a little problematic at the moment, since some of the "raw" edges might have already been removed
        //      and as such need not be evaluated. Ideally, whenever an edge is refined - it has to be changed in the main graph
        //      however, for the time being, we shall ignore this case and see how it works without taking that into account

        // for each in edge - compute the length using recursion and return max
        int max_length=0;
        for (auto ep_it : vit->second.getInEdges()){ // for each incoming edge
            int tmp_length = this->_get_read_length_before(ep_it.first);
            if (tmp_length > max_length){max_length = tmp_length;}
        }
        cur_length=cur_length+max_length;
    }
    return cur_length;
}

// This is a recursive implementation which return the total number of bases in the vertices connected by edges to the current edge
// no implicit edge evaluation is done here
// next vertices only
int HGraph::_get_read_length_after(std::map<VCoords,Vertex>::iterator vit){
    int cur_length = vit->first.getLength(); // get length of the current vertex

    if (vit->second.getOutDegree()>0){ // need to follow the exon structure for reads that span multiple introns
        // TODO: this is a little problematic at the moment, since some of the "raw" edges might have already been removed
        //      and as such need not be evaluated. Ideally, whenever an edge is refined - it has to be changed in the main graph
        //      however, for the time being, we shall ignore this case and see how it works without taking that into account

        // for each in edge - compute the length using recursion and return max
        int max_length=0;
        for (auto ep_it : vit->second.getOutEdges()){ // for each incoming edge
            int tmp_length = this->_get_read_length_after(ep_it.first);
            if (tmp_length > max_length){max_length = tmp_length;}
        }
        cur_length=cur_length+max_length;
    }
    return cur_length;
}

// This function makes sure that the total number of bases in a given read
void HGraph::enforce_read_length(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm){
    // first get the number of bases in the "raw" edge

    int total_length=0;

    // first get the maximum start over vertices preceeding the edge
    int max_length = 0;
    for (auto vit : eit.second.getPrevs()){
        // get current vertex length
        int cur_length = this->_get_read_length_before(vit);
        if (cur_length>max_length){max_length=cur_length;}
    }

    total_length = total_length + max_length;

    // second get the maximum end over the vertices after the edge
    max_length=0;
    for (auto vit : eit.second.getNexts()){
        // get current vertex length
        int cur_length = this->_get_read_length_after(vit);
        if (cur_length>max_length){max_length=cur_length;}
    }

    total_length = total_length + max_length;
    std::cout<<"total length: "<<total_length<<std::endl;

    // TODO: evaluate modified edges
    // see if they need to be removed or not
    if(((total_length - this->stats.kmerlen) + 1) < this->minKmers){
        // erase contents of the current map
        sm.clear();
    }
}

// This function computes the number of distinct vertex start sites for a given edge
void HGraph::_get_starts(std::map<VCoords,Vertex>::iterator vit,std::vector<int>& starts){
    if (vit->second.getInDegree()>0){ // need to follow the exon structure for reads that span multiple introns
        // TODO: this is a little problematic at the moment, since some of the "raw" edges might have already been removed
        //      and as such need not be evaluated. Ideally, whenever an edge is refined - it has to be changed in the main graph
        //      however, for the time being, we shall ignore this case and see how it works without taking that into account

        // for each in edge - find starts and add them to the vector of starts
        for (auto ep_it : vit->second.getInEdges()){ // for each incoming edge
            this->_get_starts(ep_it.first,starts);
        }
    }
    else{
        starts.emplace_back(vit->first.getStart());
    }
}

// This function computes the number of distinct vertex end sites for a given edge
void HGraph::_get_ends(std::map<VCoords,Vertex>::iterator vit,std::vector<int>& ends){
    if (vit->second.getOutDegree()>0){ // need to follow the exon structure for reads that span multiple introns
        // TODO: this is a little problematic at the moment, since some of the "raw" edges might have already been removed
        //      and as such need not be evaluated. Ideally, whenever an edge is refined - it has to be changed in the main graph
        //      however, for the time being, we shall ignore this case and see how it works without taking that into account

        // for each in edge - find starts and add them to the vector of starts
        for (auto ep_it : vit->second.getOutEdges()){ // for each incoming edge
            this->_get_ends(ep_it.first,ends);
        }
    }
    else{
        ends.emplace_back(vit->first.getEnd());
    }
}

// evaluates the number of unique ends and starts that are connected to a given edge
void HGraph::enforce_unique_start_end(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm){

    std::vector<int>starts{},ends{};

    // first get the maximum start over vertices preceeding the edge
    for (auto vit : eit.second.getPrevs()){
        this->_get_starts(vit,starts);
    }

    // second get the maximum end over the vertices after the edge
    for (auto vit : eit.second.getNexts()){
        this->_get_ends(vit,ends);
    }

    // now compute the entropy and compare to the threshold
    if (starts.size() == 1 && ends.size() == 1){ // TODO: also needs to take into account the weight and compute perform a more informative comparison
        sm.clear();
    }

}

// This function edits the graph based on the "raw" edge and the "parsed" edges
void HGraph::edit_graph(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm){

}

// This function tests for an overlap between two vertices
// this method only works in one direction where vit2 is greater than vit1
bool HGraph::overlap(std::map<VCoords,Vertex>::iterator vit1,std::map<VCoords,Vertex>::iterator vit2){
    return vit1->first.getChr() == vit2->first.getChr() &&
           vit1->first.getStrand() == vit2->first.getStrand() &&
           vit1->first.getStart() <= vit2->first.getStart() &&
           vit1->first.getStart() <= vit2->first.getEnd();
}

void add_bases_to_set(int start, int end, std::set<int>& cb){
    for (int i=start;i<end;i++){
        cb.insert(i);
    }
}

void HGraph::remove_vertex(std::map<VCoords,Vertex>::iterator vit){
    std::map<Edge,Aggregate_edge_props,edge_cmp>::iterator cur_emap_it;
    for (auto eit : vit->second.getInEdges()){ // for all incoming edges
        // first find iterator to the edge
        cur_emap_it = this->emap.find(std::make_pair(eit.first,vit));
        if (cur_emap_it != this->emap.end()) {
            // now we need to remove this pair of vertices from the edge
            cur_emap_it->second.remove_vertex_pair(eit.first, vit);

            if (cur_emap_it->second.isInEmpty() && cur_emap_it->second.isOutEmpty()){
                this->emap.erase(cur_emap_it); // remove the edge all together
            }
        }
        eit.first->second.remove_out_edge(vit); // now need to remove the association from the previous vertex
    }
    for (auto eit : vit->second.getOutEdges()){ // for all outgoing edges
        // first find iterator to the edge
        cur_emap_it = this->emap.find(std::make_pair(vit,eit.first));
        // now we need to remove this pair of vertices from the edge
        if (cur_emap_it != this->emap.end()) {
            cur_emap_it->second.remove_vertex_pair(vit,eit.first);

            if (cur_emap_it->second.isInEmpty() && cur_emap_it->second.isOutEmpty()){
                this->emap.erase(cur_emap_it); // remove the edge all together
            }
        }
        eit.first->second.remove_in_edge(vit);
    }
    this->vertices.remove_vt(vit); // now need to remove the actual vertex
}

// evaluates connected bases and removes related vertices and edges from the graph
void HGraph::remove_vertices(std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& vts){
    for (auto vit : vts){
        this->remove_vertex(vit);
    }
}

// This function runs a breadth first search starting with a given iterator
// stores all vertices that belong to current connected stretch into vector
// returns the iterator to the start of the next clique
std::map<VCoords,Vertex>::iterator HGraph::taylor_bfs(std::map<VCoords,Vertex>::iterator cur_vit, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited){

    std::vector<std::map<VCoords,Vertex>::iterator> queue;

    std::map<VCoords,Vertex>::iterator next_vit;

    // Mark the current node as visited and enqueue it
    visited.insert(cur_vit);
    queue.emplace_back(cur_vit);

    auto last_vit = cur_vit;

    while(!queue.empty()){
        // deque a vertex from queue and print it
        cur_vit = queue.front();
        queue.erase(queue.begin());

        // Get all adjacent vertices of the dequeued vertex.
        // If an adjacent has not been visited, then mark it visited and enqueue it

        // check if explicit edges exist
        for (auto eit : cur_vit->second.getOutEdges()){
            if(visited.find(eit.first) == visited.end()){ // element did not previously exist
                visited.insert(eit.first);
                queue.emplace_back(eit.first);
                if (last_vit->first < eit.first->first){ // is this vertex a better candidate for the last vertex in the sequence?
                    last_vit = eit.first;
                }
            }
        }

        // check if implicit edges exist
        next_vit = cur_vit;
        next_vit++;
        if (this->overlap(cur_vit,next_vit)){
            if (visited.find(next_vit) == visited.end()){
                visited.insert(next_vit);
                queue.emplace_back(next_vit);
                if(last_vit->first < next_vit->first){ // is this vertex a better candidate for the last vertex in the sequence?
                    last_vit = next_vit;
                }
            }
        }
    }
    return last_vit;
}

// this function computes the total number of bases in a given clique
int HGraph::clique_length(std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& vts){
    std::set<int> connected_bases; // all the bases in a single stretch of connected vertices (implicit and explicit edges included)
    for (auto vit : vts){
        add_bases_to_set(vit->first.getStart(),vit->first.getEnd(),connected_bases); // add current bases
    }
    return connected_bases.size();
}

// This function runs BFS variants and checks BFS-specific constraints
void HGraph::enforce_bfs_constraints(){
    std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp> cur_vertices;
    int cur_clique_length = 0;

    auto cur_vit = this->vertices.begin();
    std::map<VCoords,Vertex>::iterator next_vit;

    while (true){
        next_vit = this->taylor_bfs(cur_vit,cur_vertices);
        next_vit++;

        cur_clique_length = clique_length(cur_vertices);
//        std::cout<<"evaluated: "<<"\t"<<cur_clique_length<<std::endl;
        if (cur_clique_length < 151){
            this->remove_vertices(cur_vertices);
        }
        // clear the vertex sequence
        cur_vertices.clear();

        cur_vit = next_vit;

        if (cur_vit == this->vertices.end()){
            break;
        }
    }
}

// This is a helper function for the Depth First Search on explicit edges
void HGraph::_helper_topological_sort_explicit(std::map<VCoords,Vertex>::iterator cur_vit,
                                      std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited,
                                      std::vector<std::map<VCoords,Vertex>::iterator>& stack){

    visited.insert(cur_vit);
    for (auto eit : cur_vit->second.getOutEdges()){
        if(visited.find(eit.first) == visited.end()){ // element did not previously exist
            this->_helper_topological_sort_explicit(eit.first,visited,stack);
        }
    }
    stack.emplace_back(cur_vit);

}

// This function runs a depth first search algorithm to find the longest path within the sequence of connected vertices
// taking into account explicit edges only
// this permits filtering out partial matches that are candidate for removal as noise-multimappers
// (incomplete multimapper smaller than the read_length or the expected length given mismatches)
void HGraph::topological_sort_explicit(std::map<VCoords,Vertex>::iterator cur_vit,
                              std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited,
                              std::vector<std::map<VCoords,Vertex>::iterator>& stack){

    _helper_topological_sort_explicit(cur_vit,visited,stack);
}

void HGraph::shortest_path(){
    std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp> cur_vertices;
    std::vector<std::map<VCoords,Vertex>::iterator> stack;
    std::map<std::map<VCoords,Vertex>::iterator,int,Vertex::vmap_cmp> dist;

    auto cur_vit = this->vertices.begin();
    std::map<VCoords,Vertex>::iterator next_vit;

    while (true){
        if (cur_vit->second.getInDegree() == 0){ // check that the has no incoming edges
            this->topological_sort_explicit(cur_vit,cur_vertices,stack); // topologically sort the vertices returned by BFS

            // search for minimum distance
            for(int i=0;i<stack.size();i++){ // Initialize distances to all vertices as infinite and distance; to source as 0
                if(i+1==stack.size()){
                    dist.insert(std::make_pair(stack[i],stack[i]->first.getLength()));
                }
                else{
                    dist.insert(std::make_pair(stack[i],INF));
                }
            }

            while (!stack.empty()){ // Process vertices in topological order
                // Get the next vertex from topological order
                auto u = stack.back();
                stack.pop_back();

                // Update distances of all adjacent vertices
                std::map<VCoords,Vertex>::iterator i;
                if (dist[u] != INF){
                    int eit_idx = 0; // position_within_outEdges;
                    for (auto ei : u->second.getOutEdges()) {
                        i = ei.first;
//
                        if (dist[i] > dist[u] + (int)i->first.getLength()) {
                            dist[i] = dist[u] + (int)i->first.getLength();
                        }
                    }
                }
            }
            // final distances are in dist

            // search for maximum distance

            cur_vertices.clear(); // clear the contents of the visited vertices
            stack.clear(); // clear the contents of the stack
        }

        cur_vit++;

        if (cur_vit == this->vertices.end()){
            break;
        }
    }
}

// This is a helper function for the Depth First Search on explicit edges
// TODO: needs to store the farthest distance somehow - perhaps a different algorithm
//      also would be helpful to know the minimum distance in addition to the maximum distance
void HGraph::_helper_taylor_dfs_explicit(std::map<VCoords,Vertex>::iterator cur_vit, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited){
    visited.insert(cur_vit); // mark the current vertex as visited

    // iterate over all explicit edges
    for (auto eit : cur_vit->second.getOutEdges()){
        // if the vertex at the other end of the edge is not in visited
        if(visited.find(eit.first) == visited.end()){ // element did not previously exist
            visited.insert(eit.first);
            this->_helper_taylor_dfs_explicit(eit.first,visited);
        }
    }
}

// This function runs a depth first search algorithm to find the longest path within the sequence of connected vertices
// taking into account explicit edges only
// this permits filtering out partial matches that are candidate for removal as noise-multimappers
// (incomplete multimapper smaller than the read_length or the expected length given mismatches)
void HGraph::taylor_dfs_explicit(std::map<VCoords,Vertex>::iterator cur_vit, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited){
    _helper_taylor_dfs_explicit(cur_vit,visited);
}

// This function runs DFS variants and checks DFS-specific constraints
void HGraph::enforce_dfs_constraints(){
    std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp> cur_vertices;
    int cur_clique_length = 0;

    auto cur_vit = this->vertices.begin();
    std::map<VCoords,Vertex>::iterator next_vit;

    while (true){
        this->taylor_dfs_explicit(cur_vit,cur_vertices);

        cur_vit++;

        if (cur_vit == this->vertices.end()){
            break;
        }
    }
}

// This function parses through the vertices and removes anything that does not pass vertex-specific constraints
// such as the total length of a clique (all connected vertices)
void HGraph::parse_vertices(){

    enforce_bfs_constraints();

//    enforce_dfs_constraints(); // TODO: see if needed after topological sort
}

// This function selects the best possible combination of donor-acceptor sites
// for a given set of potential precise junctions
void HGraph::filter_best_donor_acceptor(SJS& sm){

    std::vector<SJS::iterator> bad_its;
    SJS::iterator prev_it;
    double best_score = 0.0; //(int) ((~((unsigned int) 0)) >> 1);

    double cur_score;

    for (auto sm_it = sm.begin(); sm_it != sm.end();sm_it++){
        cur_score = std::get<0>(sm_it->second) * std::get<1>(sm_it->second);
        if (cur_score < best_score){
            bad_its.emplace_back(sm_it);
        }
        else if (cur_score > best_score){
            if (sm_it != sm.begin()){
                bad_its.emplace_back(prev_it);
            }
            prev_it = sm_it;
            best_score = cur_score;
        }
        else{}
    }
    for (auto bad_it : bad_its){
        sm.erase(bad_it);
    }
}

// This function returns the most plausible match given an expanded set of edges
void HGraph::remove_overlapping_edges(SJS& sm){
    if(sm.size() > 1 ){
        filter_best_donor_acceptor(sm);
    }

    if(sm.size() > 1){ // only if an edge was actually expanded

        auto prev_it = sm.cbegin();
        int min_intron_len = 0; //(int) ((~((unsigned int) 0)) >> 1);

        int length = 0;
        for (auto sm_it = sm.cbegin(); sm_it != sm.cend();){
            length = std::get<3>(sm_it->first) - std::get<2>(sm_it->first);
            if (length < min_intron_len){
                sm_it = sm.erase(sm_it);
            }
            else{
                sm_it = sm.erase(prev_it);
                prev_it = sm_it;
                min_intron_len = length;
                ++sm_it;
            }
        }
    }
}

// this function writes out the dot-file for the current graph
void HGraph::make_dot(){
    std::string dot_fname(this->out_fname);
    dot_fname.append(".dot");
    std::ofstream dot_fp(dot_fname.c_str());

    for (const auto& vit : this->vertices){
        if (!vit.second.getOutEdges().empty()){ // write with edges
            for (auto eit : vit.second.getOutEdges()){
                dot_fp<<"\""<<vit.first.getStart()<<"-"<<vit.first.getEnd()<<"\""
                      <<"->"
                        <<"\""<<eit.first->first.getStart()<<"-"<<eit.first->first.getEnd()<<"\""<<std::endl;
            }
        }
        else{
            dot_fp<<"\""<<vit.first.getStart()<<"-"<<vit.first.getEnd()<<"\""<<std::endl;
        }
    }

    dot_fp.close();
}

// parse graph and evaluate gaps and assign splice junctions and mismatches and gaps
void HGraph::parse_graph() {
    // first need to parse the vertices, and remove any stretches that do not pass vertex-specific constraints
    std::cerr<<"\tparsing vertices"<<std::endl;
    this->parse_vertices();

    std::cerr<<"\tmaking dot"<<std::endl;
    make_dot();

    std::cerr<<"\tparsing edges"<<std::endl;

    // next iterate over the edges and output ready introns
    std::string edges_fname(this->out_fname);
    edges_fname.append(".intron.parsed.gff");
    std::ofstream edges_fp(edges_fname.c_str());

    int counter=0;

    for(const auto& eit : this->emap){
        SJS sjs_map;
        evaluate_donor_acceptor(eit,sjs_map);
        if (!sjs_map.empty()){
            counter++;
        }

        enforce_constraints(sjs_map);
//        enforce_read_length(eit,sjs_map); // TODO: currently does not quite work - introduces false negatives -  should revisit - can be done as a second vertex pass - needs to be implemented with dfs but without implicit edge evaluation
//        enforce_unique_start_end(eit,sjs_map); // TODO: currently does not quite work - introduces false negatives - should revisit - instead needs to compute the distribution of starts and ends, and see if it falls within expected thresholds

        remove_overlapping_edges(sjs_map); // TODO: currently introduces false negatives - needs work

        for(auto eit2 : sjs_map) {
            // TODO: need to remove the old edge and include the new edges in the main graph here
            //      so that upon the next edge evaluation everything is correctly accounted for.
            edges_fp << this->hdb->getContigFromID(eit.second.getChr()) << "\t" << "hairpin" << "\t" << "intron" << "\t"
                     << std::get<2>(eit2.first) << "\t" << std::get<3>(eit2.first) << "\t"
                     << "." << "\t" << eit.second.getStrand() << "\t" << "." << "\t" << "weight=" << eit.second.getWeight() << std::endl;
        }
    }
    edges_fp.close();
}

void HGraph::generate_sam_header(std::ofstream& sam_fp,std::string& cl){ // generate the header from the database information
    sam_fp << "@HD\tVN:1.0\tSO:unsorted"<<std::endl;
    this->hdb->generate_sq(sam_fp);
    sam_fp << "@PG\tID:hairpin\tPN:hairpin\tVN:1.0\tCL:\"" << cl << "\"" << std::endl;
}

// this function takes the parsed graph and output the SAM-formatted file for further analysis
void HGraph::to_sam(std::string& cl) {
    std::string sam_fname(this->out_fname);
    sam_fname.append(".sam");
    std::ofstream sam_fp(sam_fname.c_str());

    this->generate_sam_header(sam_fp,cl);

    sam_fp<<"body"<<std::endl;

    sam_fp.close();
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
