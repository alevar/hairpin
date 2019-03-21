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
    ~HDB();

    void make_db(std::string out_fname, int kmerlen);

    void save_trans_db();
    void save_genom_db();
    void save_db_info();
    void save_contig_info();

    void load_trans_db();
    void load_genom_db();
    void load_db_info();
    void load_contig_info();

private:
    GffReader gtfReader_;

    std::string gtf_fname_;
    std::string genome_fname_;

    FILE* gtf_fhandle_;

    typedef std::map<std::string, std::vector< int >* > ContigTransMap;
    ContigTransMap contigTransMap_;

    MinMap kmers;

    void transcript_map();

    HDB(); // Don't want anyone calling the constructor w/o options

    int kmerlen{};
    std::string out_fname{};

    void get_exonic_sequence(GffObj& p_trans, FastaRecord& rec, uint8_t contigID, int& lastPosition);
    void process_contig(std::string seq, uint8_t chrID, uint8_t strand, uint8_t rev_strand, bool checkTrans); // TODO: reverse complement each kmer

    void process_kmers(MinMap& mm); // process individual kmers

    MinMap trans_map;
    std::map<std::string,std::tuple<uint8_t,uint8_t,uint32_t> > genom_map;
    MinMap::KmerMap::iterator trans_map_it,genom_map_it;

    std::set<EVec> kmer_coords; // genomic positions encountered in trans_map
    std::set<EVec> kmer_coords_genom; // genomic positions encountered in genom_map
    std::pair<std::set<EVec>::iterator,bool> kmer_coords_exist,kmer_coords_genom_exist;

    std::map<std::string, uint8_t > contig_to_id;
    std::map<uint8_t, std::pair<std::string,int> > id_to_contig;
    std::pair<std::map<std::string,uint8_t>::iterator,bool> contig_exists;
    uint8_t maxID=0;
};

#endif //HAIRPIN_DB_H
