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

class HDB {
public:
    HDB(std::string gtf_fname, std::string genome_fname);
    ~HDB();

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
};

#endif //HAIRPIN_DB_H
