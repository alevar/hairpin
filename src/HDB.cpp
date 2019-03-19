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

void HDB::make_trans_db(std::string out_fname,int kmerlen) {

    this->out_fname=out_fname;
    this->kmerlen=kmerlen;

    std::vector<int> *p_contig_vec;

    FastaReader fastaReader(genome_fname_);
    FastaRecord cur_contig;

    // file for the transcriptomic map
    std::string trans_map_fname(this->out_fname);
    trans_map_fname.append(".trans");
    std::ofstream trans_map(trans_map_fname.c_str());

    while (fastaReader.good()) {

        fastaReader.next(cur_contig);

        // If this contig isn't in the map, then there are no transcripts
        // associated with it. Skip it.
        if (contigTransMap_.find(cur_contig.id_) == contigTransMap_.end()){
            continue;
        }

        p_contig_vec = contigTransMap_[cur_contig.id_];

        for (int trans_idx : *p_contig_vec) {
            GffObj *p_trans = gtfReader_.gflst.Get(trans_idx);
            HDB::get_exonic_sequence(*p_trans, cur_contig);
        }
    }
    trans_map.close();
}

void HDB::make_genom_db() { //excludes parsing kmers in the regions not spanned by the

}

void HDB::get_exonic_sequence(GffObj &p_trans, FastaRecord &rec) {
    
    GList<GffExon>& exon_list = p_trans.exons;

    std::string exon_seq;
    size_t length;
    int cur_len=0; // current length left from the previous exon
    int cur_pos=0;
    std::string sub_seq;

    EVec *cur_coords = new EVec((uint8_t)p_trans.gseq_id,(uint8_t)p_trans.strand);

    for (int i = 0; i < exon_list.Count(); ++i) {
        GffExon& cur_exon = *(exon_list.Get(i));
        length = cur_exon.end - cur_exon.start + 1;

        if (length>1){ // sanity check for 0 and 1 baes exons
            for(int j=0;j<length;j++){ // iterate over all kmers in the given exon
                if ((length-j)+cur_len<this->kmerlen){ // not enough coordinates - need to look at the next exon
                    cur_coords->_push_back((uint32_t)(cur_exon.start+j),(uint32_t)cur_exon.end);
                    cur_len+=(cur_exon.end-(cur_exon.start+j)); // save the length that has been seen thus far
                    break;
                }
                else{ // otherwise we have all the bases we need and can evaluate their uniqueness
                    if (cur_len!=0){ // some information is left from previous exons
                        for (int g=this->kmerlen-cur_len;g<this->kmerlen;g++){ // build new sequences using past coordinates
                            sub_seq="";
                            cur_coords->_push_back((uint32_t)(cur_exon.start+j),(uint32_t)(cur_exon.start+j+g));
                            this->kmer_coords_exist=this->kmer_coords.insert(cur_coords);
                            if (this->kmer_coords_exist.second){
                                for (int d=1;d<cur_coords->_getSize();d++){
                                    sub_seq+=rec.seq_.substr(cur_coords->_getStart(d)-1,cur_coords->_getEnd(d)-cur_coords->_getStart(d));
                                }
                                if(sub_seq.length()>this->kmerlen){
                                    std::cerr<<"1: "<<sub_seq.length()<<" "<<p_trans.getID()<<std::endl;
                                    for(int y=0;y<cur_coords->_getSize();y++){
                                        std::cerr<<cur_coords->_getStart(y)<<"-"<<cur_coords->_getEnd(y)<<";";
                                    }
                                    std::cerr<<std::endl;
                                }
                                this->trans_map._insert(sub_seq,cur_coords);
                            }
                            cur_len-=1;
                            cur_coords->_pop_back(); // need to test
                            cur_coords->_incStart(); //need to test
                            if (cur_coords->_getStart(0)==cur_coords->_getEnd(0)){
                                cur_coords->_eraseFirst(); // delete the first element // need to test
                            }
                        }
                        std::transform(sub_seq.begin(), sub_seq.end(), sub_seq.begin(), ::toupper);
                        // add new coordinates first
                    }
                    // need to resume from the current index
                    if ((length-j)+cur_len<this->kmerlen){ // not enough coordinates - need to look at the next exon
                        cur_coords->_push_back((uint32_t)(cur_exon.start+j),(uint32_t)(cur_exon.end));
                        cur_len+=(cur_exon.end-(cur_exon.start+j));
                        break;
                    }
                    else{
                        cur_coords->_push_back((uint32_t)(cur_exon.start+j),(uint32_t)(cur_exon.start+j+this->kmerlen));
                        this->kmer_coords_exist=this->kmer_coords.insert(cur_coords);
                        if (this->kmer_coords_exist.second){ // was successfully inserted
                            // get sequence
                            sub_seq=rec.seq_.substr(cur_exon.start+j-1,this->kmerlen);
                            if(sub_seq.length()>this->kmerlen){
                                std::cerr<<"2: "<<sub_seq.length()<<" "<<p_trans.getID()<<" ";
                                for(int y=0;y<cur_coords->_getSize();y++){
                                    std::cerr<<cur_coords->_getStart(y)<<"-"<<cur_coords->_getEnd(y)<<";";
                                }
                                std::cerr<<std::endl;
                            }
                            std::transform(sub_seq.begin(), sub_seq.end(), sub_seq.begin(), ::toupper);
                            // add new coordinates here
                        }
                    }
                    // since we went into this conditional, that means no previously encountered exons are left
                    // and we can reset some things
                    cur_coords->_erase(); // need to test
                    cur_len=0;
                }
            }
        }
        exon_seq += rec.seq_.substr(cur_exon.start - 1, length);
    }
}
