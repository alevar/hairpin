//
// Created by Ales Varabyou on 3/16/19.
// ====================================
// this class contains the implementation
// of a graph populated with kmers from sequenced reads
//

#ifndef HAIRPIN_GRAPH_H
#define HAIRPIN_GRAPH_H

#include <iostream>
#include <fstream>
#include <map>

#include "HDB.h"
#include "MinMap.h"

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
        return (this->chrID < vc.chrID || this->strand < vc.strand || (this->pos <  vc.pos) || this->getEnd() < vc.getEnd()); // the coordinate test should return equal if
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
    uint32_t getEnd() const {return this->pos+this->length;}
    uint8_t getLength() const {return this->length;}

    struct coord_hash { // in case it is implemented as an unordered_map
        uint64_t operator()(const std::vector<std::pair<int,int> > &coords ) const
        {
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

    int getOutDegree(){return static_cast<int>(this->outEdges.size());}
    int getInDegree(){return static_cast<int>(this->inEdges.size());}

    struct vmap_cmp {
        bool operator()(const std::map<VCoords,Vertex>::iterator& vit_1, const std::map<VCoords,Vertex>::iterator& vit_2) const {
            return vit_1->first < vit_2->first;
        }
    };

    std::map<std::map<VCoords,Vertex>::iterator,Edge_props, vmap_cmp> getOutEdges(){return this->outEdges;}
    std::map<std::map<VCoords,Vertex>::iterator,Edge_props, vmap_cmp> getInEdges(){return this->inEdges;}

private:
    int weight=0;

    typedef std::map<std::map<VCoords,Vertex>::iterator,Edge_props,vmap_cmp> EMap;
    typedef std::pair<std::map<std::map<VCoords,Vertex>::iterator,Edge_props>::iterator,bool> EMap_it;
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
    const_iterator begin() const { return this->vertices.begin();}
    iterator end() {return this->vertices.end();}
    const_iterator end() const { return this->vertices.end();}

private:
    std::map<VCoords,Vertex> vertices;
    std::pair<std::map<VCoords,Vertex>::iterator,bool> vm_it;
};

class Aggregate_edge_props{ // TODO: potentially entirely get rifd of storing edges in the vertices and only keep track of them here if this approach works
public:
    Aggregate_edge_props()=default;
    Aggregate_edge_props(std::map<VCoords,Vertex>::iterator prev, std::map<VCoords,Vertex>::iterator next){
        this->addNext(next);
        this->addPrev(prev);
    }

    void setKnown(){} // set the current edge as known/annotated
    bool isKnown(){} // is the current edge known/annotated
    int getWeight(){
        // TODO: needs to compute weight across all vertices that belong to the edge
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
    uint32_t getStart(){return prevs[0]->first.getEnd();}
    uint32_t getEnd(){return nexts[0]->first.getStart();}
    uint8_t getChr(){return prevs[0]->first.getChr();}
    uint8_t getStrand(){return prevs[0]->first.getStrand();}

    void addNext(std::map<VCoords,Vertex>::iterator vt){this->nexts.emplace_back(vt);}
    void addPrev(std::map<VCoords,Vertex>::iterator vt){this->prevs.emplace_back(vt);}
    std::vector<std::map<VCoords,Vertex>::iterator> getNexts(){return this->nexts;}
    std::vector<std::map<VCoords,Vertex>::iterator> getPrevs(){return this->prevs;}


private:
    bool known{};

    typedef std::vector<std::map<VCoords,Vertex>::iterator> VP; // vertex pointers
    VP nexts{},prevs{};
};

class HGraph { // TODO: add interface to substring the database to query for splice junctions
public:
    HGraph();
    explicit HGraph(HDB* hdb);
    HGraph(HDB* hdb,int max_intron,int min_intron,std::string out_fname);
    ~HGraph();

    void add_read(std::string& read);
    void print_stats();

    void add_edge(std::map<VCoords,Vertex>::iterator vit1,std::map<VCoords,Vertex>::iterator vit2);

    void to_sam(std::string out_sam_fname);
    void sort_graph();
    void parse_graph();
    void write_intron_gff();

private:

    std::map<VCoords,Vertex>::iterator add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length);

    typedef std::map<VCoords,Vertex>::iterator VIT;
    typedef std::pair<VIT,VIT> Edge;
    struct edge_cmp { // the comparator in this case simply compares if the splice site coordinates are identical
        bool operator()(const Edge& prev, const Edge& next) const {
            return prev.first->first.getEnd() < next.first->first.getEnd()
            || prev.second->first.getStart() < next.second->first.getStart();
        }
    };

    std::map<Edge,Aggregate_edge_props,edge_cmp> emap;
    std::pair<std::map<Edge,Aggregate_edge_props,edge_cmp>::iterator,bool> emap_it;

    HDB* hdb;
    std::string out_fname="";

    int maxIntron=10;
    int minIntron=500000;
    struct Stats{
        int numReads=0;
        int numReadsIgnored=0;
        int numEdges=0;
        int kmerlen;
    } stats;

    MinMap::iterator trans_it;
    HDB::GenMap::iterator genom_it;

    VMap vertices;
    std::pair<std::map<VCoords,Vertex>::iterator,bool> vertex_it;

};

#endif //HAIRPIN_GRAPH_H
