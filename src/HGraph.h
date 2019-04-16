//
// Created by Ales Varabyou on 3/16/19.
// ====================================
// this class contains the implementation
// of a graph populated with kmers from sequenced reads
//

#ifndef HAIRPIN_GRAPH_H
#define HAIRPIN_GRAPH_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <map>

#include "FastaTools.h"

#include "HDB.h"
#include "MinMap.h"

#define INF INT_MAX

class Vertex;
class VCoords;
class VMap;
class Edge;
class Edge_props;
class EMap;

class Edge_props{
public:
    Edge_props(){this->weight=1;}
    void setKnown(){} // set the current edge as known/annotated
    bool isKnown(){} // is the current edge known/annotated
    void incWeight(){this->weight++;}
    int getWeight(){return this->weight;}

private:
    bool known{};
    int weight{};
};

class VCoords{
public:
    VCoords()=default;
    VCoords(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length){
        this->chrID=chrID;
        this->strand=strand;
        this->pos=pos;
        this->length=length;
    }
    ~VCoords()=default;

    bool operator< (const VCoords& vc) const{
        uint32_t this_end = this->getEnd(), vc_end = vc.getEnd();
        return std::tie(this->chrID,this->strand,this->pos,this_end) < std::tie(vc.chrID,vc.strand,vc.pos,vc_end);
    }
    bool operator> (const VCoords& vc) const{
        return (this->chrID > vc.chrID || this->strand > vc.strand || this->pos >  vc.pos || this->getEnd() > vc.getEnd());
    }
    bool operator==(const VCoords& vc) const{
        return (this->chrID == vc.chrID && this->strand == vc.strand && this->pos ==  vc.pos && this->length == vc.length);
    }

    int operator- (const VCoords& vc) const{
        if(this->chrID == vc.chrID && this->strand == vc.strand){
            return this->pos-vc.pos;
        }
        return -1;
    }

    int findIntronLength(const VCoords& vc1, const VCoords& vc2) const{ // non functional method
        if(vc1.chrID == vc2.chrID && vc1.strand == vc2.strand){
            return vc2.pos - vc1.pos;
        }
        return -1;
    }

    uint8_t getChr() const {return this->chrID;}
    uint8_t getStrand() const {return this->strand;}
    uint32_t getStart() const {return this->pos;}
    uint32_t getEnd() const {return this->pos+(uint32_t)this->length;}
    uint8_t getLength() const {return this->length;}

    struct coord_hash { // in case it is implemented as an unordered_map
        uint64_t operator()(const std::vector<std::pair<int,int> > &coords ) const{
            uint64_t resHash=1;
            for (auto const& coord: coords) {
                resHash=resHash^std::hash<uint64_t>()(coord.first)^std::hash<uint64_t>()(coord.second);
            }
            return resHash;
        }
    };

private:
    uint8_t chrID;
    uint8_t strand;
    uint32_t pos;
    uint8_t length;
};

class Vertex{
public:
    Vertex(){this->weight=1;}

    ~Vertex()=default;

    int getWeight() {return this->weight;}
    int getOutEdgeWeight(uint32_t end){
        for(auto vt: this->outEdges){
            if(vt.first->first.getEnd()==end){
                return vt.second.getWeight();
            }
        }
        return -1;
    }
    int getInEdgeWeight(uint32_t start){
        for(auto vt: this->inEdges){
            if(vt.first->first.getStart()==start){
                return vt.second.getWeight();
            }
        }
        return -1;
    }
    void incWeight(){this->weight++;}

    int addOutEdge(std::map<VCoords,Vertex>::iterator vit){
        this->e_it = this->outEdges.insert(std::make_pair(vit,Edge_props()));
        if(!this->e_it.second){ // not inserted - need to increment the weight
            this->e_it.first->second.incWeight();
            return 0;
        }
        return 1;
    }
    int addInEdge(std::map<VCoords,Vertex>::iterator vit){
        this->e_it = this->inEdges.insert(std::make_pair(vit,Edge_props()));
        if(!this->e_it.second){ // not inserted - need to increment the weight
            this->e_it.first->second.incWeight();
            return 0;
        }
        return 1;
    }

    int getOutDegree() const {return static_cast<int>(this->outEdges.size());}
    int getInDegree(){return static_cast<int>(this->inEdges.size());}

    struct vmap_cmp {
        bool operator()(const std::map<VCoords,Vertex>::iterator& vit_1, const std::map<VCoords,Vertex>::iterator& vit_2) const {
            return vit_1->first < vit_2->first;
        }
    };

    std::map<std::map<VCoords,Vertex>::iterator,Edge_props, vmap_cmp> getOutEdges() const {return this->outEdges;}
    std::map<std::map<VCoords,Vertex>::iterator,Edge_props, vmap_cmp> getInEdges(){return this->inEdges;}

    void remove_out_edge(std::map<VCoords,Vertex>::iterator vit){
        auto emap_it_found = this->outEdges.find(vit);
//        if (emap_it_found == this->outEdges.end()){
//            std::cerr<<"out edge not found"<<std::endl;
//        }
        this->outEdges.erase(emap_it_found);
    }
    void remove_in_edge(std::map<VCoords,Vertex>::iterator vit){
        auto emap_it_found = this->inEdges.find(vit);
//        if (emap_it_found == this->inEdges.end()){
//            std::cerr<<"in edge not found"<<std::endl;
//        }
        this->inEdges.erase(emap_it_found);
    }

private:
    int weight=0;

    typedef std::map<std::map<VCoords,Vertex>::iterator,Edge_props,vmap_cmp> EMap;
    typedef std::pair<std::map<std::map<VCoords,Vertex>::iterator,Edge_props,vmap_cmp>::iterator,bool> EMap_it;
    std::map<std::map<VCoords,Vertex>::iterator,Edge_props,vmap_cmp> outEdges{};
    std::map<std::map<VCoords,Vertex>::iterator,Edge_props,vmap_cmp> inEdges{};
    EMap_it e_it;
};

class VMap{
public:
    VMap()=default;
    ~VMap()=default;

    std::map<VCoords,Vertex>::iterator _insert(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length){
        this->vm_it=this->vertices.insert(std::make_pair(VCoords(chrID,strand,pos,length),Vertex()));
        if(!this->vm_it.second){ // already exists
            this->vm_it.first->second.incWeight();
        }
        return this->vm_it.first;
    }

    int getSize(){return static_cast<int>(this->vertices.size());}

    typedef std::map<VCoords,Vertex>::iterator iterator;
    typedef std::map<VCoords,Vertex>::const_iterator const_iterator;
    typedef std::map<VCoords,Vertex>::reference reference;
    iterator begin() {return this->vertices.begin();}
    const_iterator cbegin() const { return this->vertices.cbegin();}
    iterator end() {return this->vertices.end();}
    const_iterator cend() const { return this->vertices.cend();}

    bool _exists(VMap::iterator vit){
        return this->vertices.find(vit->first) != this->vertices.end();
    }
    bool _exists(VMap::const_iterator vit){
        return this->vertices.find(vit->first) != this->vertices.end();
    }

    void remove_vt(std::map<VCoords,Vertex>::iterator vit){
        VCoords searching_vit = vit->first;
//        if (this->vertices.find(vit->first) == this->vertices.end()){
//            std::cerr<<"'VCoords not found"<<std::endl;
//        }
        // let's do a check in the out and in edges here as well
//        if (!vit->second.getOutEdges().empty()){
//            std::cerr<<"something's wrong out"<<std::endl;
//        }
//        if (!vit->second.getInEdges().empty()){
//            std::cerr<<"something's wrong in"<<std::endl;
//        }
        this->vertices.erase(vit);
//        if (this->vertices.find(vit->first) != this->vertices.end()){
//            std::cerr<<"the vertex was not removed"<<std::endl;
//        }
    }

private:
    std::map<VCoords,Vertex> vertices;
    std::pair<std::map<VCoords,Vertex>::iterator,bool> vm_it;
};

class Aggregate_edge_props{ // TODO: potentially entirely get rid of storing edges in the vertices and only keep track of them here if this approach works
public:
    Aggregate_edge_props()=default;
    Aggregate_edge_props(std::map<VCoords,Vertex>::iterator prev, std::map<VCoords,Vertex>::iterator next){
        this->addNext(next);
        this->addPrev(prev);
    }

    void setKnown(){this->known=true;} // set the current edge as known/annotated
    bool isKnown(){return this->known;} // is the current edge known/annotated
    void validate(){
        this->valid=true;
    } // for an unknown edge - modify the edge and the vertices to the inferred ones; refine coordinates and set the validation flag
    int getWeight() const {
        int total=0;
        for(auto vt : this->nexts){ // there is a problem here - we can not go to the exact edge in the vertex which is needed to compute the weight of a given edge
            int wt = vt->second.getInEdgeWeight(prevs[0]->first.getStart());
            total=total+wt;
            if(wt<0){
                std::cerr<<"something is wrong when getting edges"<<std::endl;
                exit(-1);
            }
        }
        return total;
    }
    uint32_t getStart() const {return prevs[0]->first.getEnd();}
    uint32_t getEnd() const {return nexts[0]->first.getStart();}
    uint8_t getChr() const {return prevs[0]->first.getChr();}
    uint8_t getStrand() const {return prevs[0]->first.getStrand();}

    void addNext(std::map<VCoords,Vertex>::iterator vt){this->nexts.emplace_back(vt);}
    void addPrev(std::map<VCoords,Vertex>::iterator vt){this->prevs.emplace_back(vt);}
    std::vector<std::map<VCoords,Vertex>::iterator> getNexts() const {return this->nexts;}
    std::vector<std::map<VCoords,Vertex>::iterator> getPrevs() const {return this->prevs;}

    bool _exists_pair(std::map<VCoords,Vertex>::iterator prev,std::map<VCoords,Vertex>::iterator next){
        auto nit = std::find(this->nexts.begin(),this->nexts.end(),next);
        auto pit = std::find(this->prevs.begin(),this->prevs.end(),prev);
        if (nit != this->nexts.end() && pit != this->prevs.end()) { // pair exists
            return true;
        }
        return false;
    }

    void remove_vertex_pair(std::map<VCoords,Vertex>::iterator prev, std::map<VCoords,Vertex>::iterator next){
        auto nit = std::find(this->nexts.begin(),this->nexts.end(),next);
        if (nit != this->nexts.end()) {
            this->nexts.erase(nit);
        }
        auto pit = std::find(this->prevs.begin(),this->prevs.end(),prev);
        if (pit != this->prevs.end()){
            this->prevs.erase(pit);
        }

//        auto nit = this->nexts.erase(std::remove(this->nexts.begin(), this->nexts.end(), next), this->nexts.end());
//        auto pit = this->prevs.erase(std::remove(this->prevs.begin(), this->prevs.end(), prev), this->prevs.end());
    }
    int getOutSize(){return this->nexts.size();} // return the number of incoming vertices
    int getInSize(){return this->prevs.size();} // return the number of outgoing vertices
    bool isOutEmpty(){return this->nexts.empty();} // return the number of incoming vertices
    bool isInEmpty(){return this->prevs.empty();} // return the number of outgoing vertices


private:
    bool known{};
    bool valid{};

    typedef std::vector<std::map<VCoords,Vertex>::iterator> VP; // vertex pointers
    VP nexts{},prevs{};
};

class HGraph {
public:
    HGraph();
    explicit HGraph(HDB* hdb);
    HGraph(HDB* hdb,int max_intron,int min_intron,int min_mismatch,std::string out_fname);
    ~HGraph();

    void add_read(std::string& read);
    void print_stats();

    void add_edge(std::map<VCoords,Vertex>::iterator vit1,std::map<VCoords,Vertex>::iterator vit2);

    void _helper_tay_dfs(std::map<VCoords,Vertex>::iterator cur_vit, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited);
    void tay_dfs(std::map<VCoords,Vertex>::iterator cur_vit, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited);

    void to_sam(std::string& cl);
    void sort_graph();
    void parse_graph();
    void write_intron_gff();

private:

    int minMismatch=0;

    std::map<VCoords,Vertex>::iterator add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length);

    typedef std::map<VCoords,Vertex>::iterator VIT;
    typedef std::pair<VIT,VIT> Edge;
    struct edge_cmp { // the comparator in this case simply compares if the splice site coordinates are identical
        bool operator()(const Edge& prev, const Edge& next) const {
            uint32_t prev_end = prev.first->first.getEnd(), next_end = next.first->first.getEnd();
            uint32_t prev_start = prev.second->first.getStart(), next_start = next.second->first.getStart();
            return std::tie(prev_end,prev_start) < std::tie(next_end,next_start);
        }
    };

    std::map<Edge,Aggregate_edge_props,edge_cmp> emap,final_emap;
    std::pair<std::map<Edge,Aggregate_edge_props,edge_cmp>::iterator,bool> emap_it,final_emap_it;

    HDB* hdb;
    std::string out_fname="";

    int maxIntron=10;
    int minIntron=500000;
    int minKmers=65; // TODO: minimum number of kmers that have to be matched from a single read - needs to be computed during runtime as a function of the read length
    struct Stats{
        int numReads=0;
        int numReadsIgnored=0;
        int numEdges=0;
        int kmerlen;
        int numKmerUnmatched=0;
        int numKmerMatchedTrans=0;
        int numKmerMatchedGenom=0;
    } stats;

    MinMap::iterator trans_it;
    HDB::GenMap::iterator genom_it;

    VMap vertices,final_vertices;
    std::pair<std::map<VCoords,Vertex>::iterator,bool> vertex_it,final_vertex_it;

    int getGenomeSubstr(Edge et,int overhang, std::string& sub_seq);
    int getGenomeSubstr(uint8_t chrID,int8_t strand, uint32_t start, uint8_t length, std::string& sub_seq);

    typedef std::tuple<uint8_t,uint8_t,uint32_t> CoordVec;
    int coord_distance(const CoordVec& cv1,const CoordVec& cv2);

    struct map_vit_cmp { // this comparator can be applied whenever there is a need to build a sorted object (eg. std::map) with std::map<VCoords,Vertex>::iterator as a key
        bool operator()(const std::map<VCoords,Vertex>::iterator& prev, const std::map<VCoords,Vertex>::iterator& next) const {
            uint32_t prev_chr = prev->first.getChr(), prev_strand = prev->first.getStrand(), prev_start = prev->first.getStart(), prev_end = prev->first.getEnd();
            uint32_t next_chr = next->first.getChr(), next_strand = next->first.getStrand(), next_start = next->first.getStart(), next_end = next->first.getEnd();
            return std::tie(prev_chr,prev_strand,prev_start,prev_end) < std::tie(next_chr,next_strand,next_start,next_end);
        }
    };

    typedef std::pair<uint8_t,uint8_t> PosPair; // first -> length of the associated Coords Vector; second - position of the kmer within the current read

    void generate_sam_header(std::ofstream& sam_fp,std::string& cl);

    std::map<std::string,double> donors = {
            {"GT",1.0},
            {"GC",0.5},
            {"AT",0.5},
    };
    std::map<std::string,double> acceptors = {
            {"AG",1.0},
            {"CA",0.5},
            {"AC",0.5}
    };

    typedef std::tuple<uint8_t,uint8_t,uint32_t,uint32_t> SJ; // defines the type of a splice junction
    struct sjs_cmp { // Comparator for the SJS map
        bool operator()(const SJ& prev, const SJ& next) const {
//            return std::tie(std::get<0>(prev),std::get<1>(prev),std::get<2>(prev),std::get<3>(prev)) < std::tie(std::get<0>(next),std::get<1>(next),std::get<2>(next),std::get<3>(next));
            return std::get<0>(prev) < std::get<0>(next) || std::get<1>(prev) < std::get<1>(next) ||
                   std::get<2>(prev) < std::get<2>(next) || std::get<3>(prev) < std::get<3>(next);
        }
    };

    bool overlap(std::map<VCoords,Vertex>::iterator vit1,std::map<VCoords,Vertex>::iterator vit2);
    int clique_length(std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& vts);
    std::map<VCoords,Vertex>::iterator taylor_bfs(std::map<VCoords,Vertex>::iterator, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& vts);

    void _helper_taylor_dfs_explicit(std::map<VCoords,Vertex>::iterator cur_vit, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited,std::vector<int>& lens,int cur_len);
    void taylor_dfs_explicit(std::map<VCoords,Vertex>::iterator cur_vit, std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited,std::vector<int>& lens);

    void _helper_topological_sort_explicit(std::map<VCoords,Vertex>::iterator cur_vit,
            std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited,
            std::vector<std::map<VCoords,Vertex>::iterator>& stack);
    void topological_sort_explicit(std::map<VCoords,Vertex>::iterator cur_vit,
            std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& visited,
            std::vector<std::map<VCoords,Vertex>::iterator>& stack);
    void shortest_path();

    void remove_vertex(std::map<VCoords,Vertex>::iterator vit);
    void remove_vertex_dfs(std::map<VCoords,Vertex>::iterator vit);
    void remove_vertices(std::set<std::map<VCoords,Vertex>::iterator,Vertex::vmap_cmp>& vts);
    void enforce_bfs_constraints();
    void enforce_dfs_constraints();
    void parse_vertices();

    typedef std::map<SJ,std::tuple<double,double>,sjs_cmp> SJS; // defines a type for a map of splice junctions

    void edit_graph(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm);

    void evaluate_sj(const std::pair<Edge,Aggregate_edge_props>& eit,const std::pair<std::string,double>& donor,const std::pair<std::string,double>& acceptor,SJS& sjs, std::string& sub_seq);

    void filter_best_donor_acceptor(SJS& sm);
    void remove_overlapping_edges(SJS& sm);

    void enforce_constraints(SJS& sm);

    int _get_read_length_before(std::map<VCoords,Vertex>::iterator vit);
    int _get_read_length_after(std::map<VCoords,Vertex>::iterator vit);
    void enforce_read_length(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm);

    void _get_starts(std::map<VCoords,Vertex>::iterator vit,std::vector<int>& starts);
    void _get_ends(std::map<VCoords,Vertex>::iterator vit,std::vector<int>& ends);
    void enforce_unique_start_end(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm);

    int compute_minimal_clique_length(SJS& sm);
    int compute_maximal_clique_length(SJS& sm);
    int compute_minimal_exon_length(SJS& sm);
    int compute_maximal_exon_length(SJS& sm);

    void evaluate_donor_acceptor(const std::pair<Edge,Aggregate_edge_props>& eit, SJS& sm);
    uint8_t getEdgeChr(const std::pair<Edge,Aggregate_edge_props>& eit);
    uint8_t getEdgeStrand(const std::pair<Edge,Aggregate_edge_props>& eit);
    uint8_t getEdgeStart(const std::pair<Edge,Aggregate_edge_props>& eit);
    uint8_t getEdgeEnd(const std::pair<Edge,Aggregate_edge_props>& eit);

    void make_dot();

    void testAllEdges();
    void testEdge(std::map<Edge,Aggregate_edge_props,edge_cmp>::iterator eit);

};

#endif //HAIRPIN_GRAPH_H
