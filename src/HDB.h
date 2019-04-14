//
// Created by Ales Varabyou on 3/16/19.
// ====================================
// this class contains the implementation
// of a database required for the hairpin operation
// this class allows building and loading the database
//

#ifndef HAIRPIN_DB_H
#define HAIRPIN_DB_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>

#include "gff.h"
#include "GFaSeqGet.h"
#include "FastaTools.h"
#include "MinMap.h"
#include "gdna.h"

class HDB {
public:
    HDB(std::string gtf_fname, std::string genome_fname);
    HDB();
    ~HDB();

    void init(std::string genome_fname);
    void init(std::string gtf_fname,std::string genome_fname);

    void make_db(std::string out_fname, int kmerlen);

    typedef std::vector<std::tuple<uint8_t,uint8_t,uint32_t>> GenVec;
    typedef std::map<std::string,GenVec> GenMap;

    MinMap::iterator find_trans(std::string& kmer);
    GenMap::iterator find_genom(std::string& kmer);

    typedef GenMap::iterator genom_iterator;
    typedef GenMap::const_iterator const_genom_iterator;
    typedef GenMap::reference genom_reference;
    genom_iterator genom_begin() {return genom_map.begin();}
    const_genom_iterator genom_begin() const { return genom_map.begin();}
    genom_iterator genom_end() {return genom_map.end();}
    const_genom_iterator genom_end() const { return genom_map.end();}

    typedef MinMap::iterator trans_iterator;
    typedef MinMap::const_iterator const_trans_iterator;
    typedef MinMap::reference trans_reference;
    trans_iterator trans_begin() {return trans_map.begin();}
    const_trans_iterator trans_begin() const { return trans_map.begin();}
    trans_iterator trans_end() {return trans_map.end();}
    const_trans_iterator trans_end() const { return trans_map.end();}

    void save_db();
    void load_db(std::string fb_fname_base);

    std::string getGenomeFname(){return this->genome_fname_;}

    int getKmerLen();

    std::string getContigFromID(uint8_t id){return id_to_contig[id].first;}
    uint8_t getIDFromContig(std::string id){return contig_to_id[id];}

    std::string getContigStr(uint8_t id){
        return this->contigs[id];
    }

    void generate_sq(std::ofstream& fp);

private:
    bool transcriptomeBuild=true; // whether a transcriptome index is requested
    GffReader gtfReader_;

    std::string gtf_fname_;
    std::string genome_fname_;

    FILE* gtf_fhandle_;

    typedef std::map<std::string, std::vector< int >* > ContigTransMap;
    ContigTransMap contigTransMap_;

    MinMap kmers;

    void transcript_map();

    int kmerlen{};
    std::string out_fname{};

    void get_exonic_sequence(GffObj& p_trans, FastaRecord& rec, uint8_t contigID, int& lastPosition);
    void process_contig(std::string seq, uint8_t chrID, uint8_t strand, uint8_t rev_strand, bool checkTrans);

    void process_kmers(MinMap& mm); // process individual kmers

    MinMap trans_map;
    GenMap genom_map;
    std::pair<GenMap::iterator,bool> genom_map_it;
    MinMap::KmerMap::iterator trans_map_it;

    std::set<EVec> kmer_coords; // genomic positions encountered in trans_map
    std::set<EVec> kmer_coords_genom; // genomic positions encountered in genom_map
    std::pair<std::set<EVec>::iterator,bool> kmer_coords_exist,kmer_coords_genom_exist;

    std::map<std::string, uint8_t > contig_to_id;
    std::map<uint8_t, std::pair<std::string,int> > id_to_contig;
    std::pair<std::map<std::string,uint8_t>::iterator,bool> contig_exists;
    uint8_t maxID=0;

    void save_trans_db();
    void save_genom_db();
    void save_db_info();
    void save_contig_info();

    void load_trans_db(std::ifstream& stream);
    void load_genom_db(std::ifstream& stream);
    void load_db_info(std::ifstream& stream);
    void load_contig_info(std::ifstream& fp);

    void load_contigs();

    std::map<uint8_t,std::string> contigs;
};

#endif //HAIRPIN_DB_H
