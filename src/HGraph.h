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
        return (this->chrID < vc.chrID || this->strand < vc.strand || this->pos <  vc.pos);
    }
    bool operator> (const VCoords& vc) const{
        return (this->chrID > vc.chrID || this->strand > vc.strand || this->pos >  vc.pos);
    }
    bool operator==(const VCoords& vc) const{
        return (this->chrID == vc.chrID && this->strand == vc.strand && this->pos ==  vc.pos && this->length == vc.length);
    }

    int operator- (const VCoords& vc) const{
        if(this->chrID == vc.chrID && this->strand == vc.strand){
            return vc.pos - this->pos;
        }
        return -1;
    }

    uint8_t getChr() const {return this->chrID;}

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
    void incWeight(){this->weight++;}

    int addOutEdge(std::map<VCoords,Vertex>::iterator vit){
        this->e_it = this->inEdges.insert(std::make_pair(vit,Edge_props()));
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
    void test(){
        for(auto vit: this->outEdges){
            std::cout<<vit.first->first.getChr()<<std::endl;
        }
    }

    int getOutDegree(){return static_cast<int>(this->outEdges.size());}
    int getInDegree(){return static_cast<int>(this->inEdges.size());}

private:
    int weight=0;

    struct vmap_cmp {
        bool operator()(const std::map<VCoords,Vertex>::iterator& vit_1, const std::map<VCoords,Vertex>::iterator& vit_2) const {
            return vit_1->first < vit_2->first;
        }
    };

    typedef std::map<std::map<VCoords,Vertex>::iterator,Edge_props,vmap_cmp> EMap;
    typedef std::pair<std::map<std::map<VCoords,Vertex>::iterator,Edge_props>::iterator,bool> EMap_it;
    EMap outEdges{};
    EMap inEdges{};
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
    void _addOutEdge(Vertex& vt){}
    void _addInEdge(Vertex& vt){}

    int getSize(){return static_cast<int>(this->vertices.size());}

private:
    std::map<VCoords,Vertex> vertices;
    std::pair<std::map<VCoords,Vertex>::iterator,bool> vm_it;
};

class HGraph {
public:
    HGraph();
    explicit HGraph(HDB* hdb);
    HGraph(HDB* hdb,int max_intron,int min_intron);
    ~HGraph();

    void add_read(std::string& read);
    void print_stats();

    void to_sam(std::string out_sam_fname);
    void sort_graph();
    void parse_graph();

private:

    std::map<VCoords,Vertex>::iterator add_vertex(uint8_t chrID,uint8_t strand,uint32_t pos,uint8_t length);
//    TODO: std::map<Edge,Edge_props>::iterator add_edge(std::map<VCoords,Vertex>::iterator vt_1,std::map<VCoords,Vertex>::iterator vt_2);

    HDB* hdb;

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
