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

    void load_trans_db();
    void load_genom_db();

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

    void get_exonic_sequence(GffObj& p_trans, FastaRecord& rec);
    void get_contig_sequence(FastaRecord& rec); // TODO: reverse complement each kmer

    MinMap trans_map,genom_map;
    MinMap::KmerMap::iterator trans_map_it,genom_map_it;

    std::set<EVec> kmer_coords; // genomic positions encountered
    std::pair<std::set<EVec>::iterator,bool> kmer_coords_exist;
};

#endif //HAIRPIN_DB_H
