//
// Modified by Ales Varabyou on 3/16/19.
//  Created by Harold Pimentel on 10/26/11.
//
#include "HDB.h"

HDB::HDB() = default;

HDB::HDB(std::string gtf_fname, std::string genome_fname){
    gtf_fname_ = gtf_fname;
    gtf_fhandle_ = fopen(gtf_fname_.c_str(), "r");
    if (gtf_fhandle_ == nullptr){
        std::cerr << "FATAL: Couldn't open annotation: " << gtf_fname_<< std::endl;
        exit(1);
    }
    std::cout << "Reading the annotation file: " << gtf_fname_ << std::endl;
    gtfReader_.init(gtf_fhandle_, true); //load recognizable transcript features only
    gtfReader_.readAll();

    genome_fname_ = genome_fname;

    // Make a map from the GffObj
    transcript_map();
}

HDB::~HDB(){
    ContigTransMap::iterator it;

    for (it = contigTransMap_.begin(); it != contigTransMap_.end(); ++it) {
        delete it->second;
    }
}

void HDB::transcript_map(){
    GffObj *p_gffObj;
    const char *p_contig_name;
    std::vector<int> *p_contig_vec;

    for (int i = 0; i < gtfReader_.gflst.Count(); ++i){
        p_gffObj = gtfReader_.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0) {
            continue;
        }

        p_contig_name = p_gffObj->getRefName();
        std::string contig_name(p_contig_name);

        // Check if the current contig exists in the map
        // If it doesn't, add it
        if (contigTransMap_.find(contig_name) == contigTransMap_.end()){
            p_contig_vec = new std::vector<int>;
            contigTransMap_[contig_name] = p_contig_vec;
        }
        else{
            p_contig_vec = contigTransMap_[contig_name];
        }

        p_contig_vec->push_back(i);
    }
}
