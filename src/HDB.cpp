//
// Modified by Ales Varabyou on 3/16/19.
// Based partially on gtf_to_fasta by Harold Pimentel from 10/26/11.
//
#include "HDB.h"

HDB::HDB() = default;

HDB::HDB(std::string gtf_fname, std::string genome_fname){
    this->gtf_fname_ = gtf_fname;
    this->gtf_fhandle_ = fopen(this->gtf_fname_.c_str(), "r");
    if (gtf_fhandle_ == nullptr){
        std::cerr << "FATAL: Couldn't open annotation: " << this->gtf_fname_<< std::endl;
        exit(1);
    }
    std::cout << "Reading the annotation file: " << this->gtf_fname_ << std::endl;
    this->gtfReader_.init(this->gtf_fhandle_, true); //load recognizable transcript features only
    this->gtfReader_.readAll();

    this->genome_fname_ = genome_fname;

    // Make a map from the GffObj
    this->transcript_map();
}

HDB::~HDB(){
    ContigTransMap::iterator it;

    for (it = contigTransMap_.begin(); it != contigTransMap_.end(); ++it) {
        delete it->second;
    }
}

void HDB::init(std::string genome_fname){
    this->genome_fname_ = genome_fname;
    this->transcriptomeBuild=false;
}

void HDB::init(std::string gtf_fname, std::string genome_fname){
    this->gtf_fname_ = gtf_fname;
    this->gtf_fhandle_ = fopen(this->gtf_fname_.c_str(), "r");
    if (gtf_fhandle_ == nullptr){
        std::cerr << "FATAL: Couldn't open annotation: " << this->gtf_fname_<< std::endl;
        exit(1);
    }
    std::cout << "Reading the annotation file: " << this->gtf_fname_ << std::endl;
    this->gtfReader_.init(this->gtf_fhandle_, true); //load recognizable transcript features only
    this->gtfReader_.readAll();

    this->genome_fname_ = genome_fname;

    // Make a map from the GffObj
    this->transcript_map();
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

void HDB::make_db(std::string out_fname, int kmerlen) {

    this->out_fname=out_fname;
    this->kmerlen=kmerlen;

    std::vector<int> *p_contig_vec;

    FastaReader fastaReader(genome_fname_);
    FastaRecord cur_contig;

    // file for the transcriptomic map
    std::string trans_map_fname(this->out_fname);
    trans_map_fname.append(".trans");
    std::ofstream trans_map(trans_map_fname.c_str());

    int lastPosition=0; // this tells what was the last position encountered. Passed by reference and modified. It is used by the get_exonic_sequence function to process non-transcriptomic regions

    while (fastaReader.good()) {

        fastaReader.next(cur_contig);

        // first make sure that the contig info is added to the contig id map along with the length
        this->contig_exists = this->contig_to_id.insert(std::make_pair(cur_contig.id_,this->maxID+1));
        if (contig_exists.second){ //contig info was successfully inserted
            this->maxID++;
            this->id_to_contig.insert(std::make_pair(this->maxID,std::make_pair(cur_contig.id_,cur_contig.seq_.length())));
        }
        uint8_t contigID = this->contig_exists.first->second; //this is the id of the new contig

        // now proceed to safelly process the information from the contig
        if (contigTransMap_.find(cur_contig.id_) == contigTransMap_.end() || !this->transcriptomeBuild){ // no transcript found on this contig or just the genome index requested - just add to the genomic map
            HDB::process_contig(cur_contig.seq_,contigID,'+','-', false);
        }
        else{
            p_contig_vec = contigTransMap_[cur_contig.id_];

            for (int trans_idx : *p_contig_vec) {
                GffObj *p_trans = gtfReader_.gflst.Get(trans_idx);
                HDB::get_exonic_sequence(*p_trans, cur_contig,contigID, lastPosition);
            }
            HDB::process_contig(cur_contig.seq_,contigID,'+','-', true);
        }
    }
    trans_map.close();
}

// this method creates a kmer map for the genome
// the method ignores any sequence spanned by transcriptome
// additionally, the method ignores any kmeres with Ns
void HDB::process_contig(std::string seq, uint8_t chrID, uint8_t strand, uint8_t rev_strand, bool checkTrans){
    if(checkTrans){
        for(uint32_t i=0;i<(seq.length()-this->kmerlen)+1;i++){ // iterate over the kmers in the current sequence and put them all into the map
            if(i%1000000==0){
                std::cout<<"process_contig:"<<i<<std::endl;
            }
            // now check if coordinates have been seen when parsing the transcriptome
            EVec cur_coords_genom(chrID,'+');
            cur_coords_genom._push_back((uint32_t) i, (uint32_t) (i+this->kmerlen));
            this->kmer_coords_exist = this->kmer_coords.insert(cur_coords_genom);
            if (!this->kmer_coords_exist.second) { // kmer was not inserted
                continue;
            }

            std::string cur_kmer=seq.substr(i,this->kmerlen);
            std::transform(cur_kmer.begin(), cur_kmer.end(), cur_kmer.begin(), ::toupper);
            size_t n_pos=cur_kmer.rfind('N');
            if(n_pos!=std::string::npos){ // found n and need to now skip by the kmerlen
                i+=n_pos;
                continue;
            }
            this->genom_map_it=this->genom_map.insert(std::pair<std::string,GenVec>(cur_kmer,{}));
            this->genom_map_it.first->second.push_back(std::make_tuple(chrID,strand,i));
            char *rev_cur_kmer = new char[this->kmerlen];
            strcpy(rev_cur_kmer, cur_kmer.c_str());
            reverseComplement(rev_cur_kmer,this->kmerlen);
            this->genom_map_it=this->genom_map.insert(std::pair<std::string,GenVec>(rev_cur_kmer,{}));
            this->genom_map_it.first->second.push_back(std::make_tuple(chrID,rev_strand,i));
            delete [] rev_cur_kmer;
        }
    }
    else {
        for(uint32_t i=0;i<(seq.length()-this->kmerlen)+1;i++){ // iterate over the kmers in the current sequence and put them all into the map
            if(i%1000000==0){
                std::cout<<"process_contig:"<<i<<std::endl;
            }
            std::string cur_kmer=seq.substr(i,this->kmerlen);
            std::transform(cur_kmer.begin(), cur_kmer.end(), cur_kmer.begin(), ::toupper);
            size_t n_pos=cur_kmer.find_first_of('N');
            if(n_pos!=std::string::npos){ // found n and need to now skip by the kmerlen
                i+=n_pos;
                continue;
            }
            this->genom_map_it=this->genom_map.insert(std::pair<std::string,GenVec>(cur_kmer,{}));
            this->genom_map_it.first->second.push_back(std::make_tuple(chrID,strand,i));
            char *rev_cur_kmer = new char[this->kmerlen];
            strcpy(rev_cur_kmer, cur_kmer.c_str());
            reverseComplement(rev_cur_kmer,this->kmerlen);
            this->genom_map_it=this->genom_map.insert(std::pair<std::string,GenVec>(rev_cur_kmer,{}));
            this->genom_map_it.first->second.push_back(std::make_tuple(chrID,rev_strand,i));
            delete [] rev_cur_kmer;
        }
    }
}

// should eventually take the code that is currently hosted in the get_exonic_sequence
void HDB::process_kmers(MinMap& mm){

}

void HDB::get_exonic_sequence(GffObj &p_trans, FastaRecord &rec, uint8_t contigID, int& lastPosition) {
    
    GList<GffExon>& exon_list = p_trans.exons;

    std::string exon_seq;
    size_t length;
    int cur_len=0; // current length left from the previous exon
    int cur_pos=0;
    std::string sub_seq;
    std::string genomic_sub_seq;

    EVec cur_coords(contigID,(uint8_t)p_trans.strand);

    bool exhaustedExon=false;

    for (int i = 0; i < exon_list.Count(); ++i) {
        GffExon& cur_exon = *(exon_list.Get(i));
        length = (cur_exon.end+1) - cur_exon.start;

        for (int j = 0; j < length; j++) { // iterate over all kmers in the given exon
            if ((length - j) + cur_len < this->kmerlen) { // not enough coordinates - need to look at the next exon
                cur_coords._push_back((uint32_t) (cur_exon.start + j), (uint32_t) cur_exon.end);
                cur_len += ((cur_exon.end+1) - (cur_exon.start + j)); // save the length that has been seen thus far
                break;
            } else { // otherwise we have all the bases we need and can evaluate their uniqueness
                if (cur_len != 0) { // some information is left from previous exons
                    for (int g = this->kmerlen - cur_len; g < this->kmerlen; g++) { // build new sequences using past coordinates
                        sub_seq = "";
                        cur_coords._push_back((uint32_t) (cur_exon.start + j), (uint32_t) (cur_exon.start + j + g-1));
                        this->kmer_coords_exist = this->kmer_coords.insert(cur_coords);
                        if (this->kmer_coords_exist.second) { // these coordinates have not been previously seen
                            for (int d = 0; d < cur_coords._getSize(); d++) {
                                sub_seq += rec.seq_.substr(cur_coords._getStart(d) - 1,
                                                           (cur_coords._getEnd(d)+1) - cur_coords._getStart(d));
                            }
                            if (sub_seq.length() != this->kmerlen) {
                                std::cerr << "1: " << j << " " << g << "\t" << sub_seq.length() << " "
                                          << p_trans.getID() << std::endl;
                                for (int y = 0; y < cur_coords._getSize(); y++) {
                                    std::cerr << cur_coords._getStart(y) << "-" << cur_coords._getEnd(y) << ";";
                                }
                                std::cerr << std::endl;
                            }
                            std::transform(sub_seq.begin(), sub_seq.end(), sub_seq.begin(), ::toupper);
                            this->trans_map._insert(sub_seq, cur_coords);
                            // do reverse_complement
                            char *rev_cur_kmer = new char[this->kmerlen];
                            strcpy(rev_cur_kmer, sub_seq.c_str());
                            reverseComplement(rev_cur_kmer, this->kmerlen);
//                            cur_coords.flipStrand(); // change strand information for the reverse complement
                            this->trans_map._insert(std::string(rev_cur_kmer), cur_coords);
//                            cur_coords.flipStrand(); // reset the strand information back
                            delete[] rev_cur_kmer;
                        }
                        cur_len -= 1;
                        cur_coords._pop_back();
                        cur_coords._incStart();
                        if (cur_coords._getStart(0)-1 == cur_coords._getEnd(0)) {
                            cur_coords._eraseFirst(); // delete the first element // need to test
                        }
                        if((length - j) + cur_len < this->kmerlen){ // check again, since we are decreasing the cur_len with each iteration, thus it is possible that the length may not be sufficient
                            cur_coords._push_back((uint32_t) (cur_exon.start + j), (uint32_t) cur_exon.end);
                            cur_len += ((cur_exon.end+1) - (cur_exon.start + j)); // save the length that has been seen thus far
                            exhaustedExon=true; //set flag to exit from both for loops
                            break;
                        }
                    }
                    if(exhaustedExon){
                        exhaustedExon=false;
                        break;
                    }
                    // add new coordinates first
                }
                // need to resume from the current index
                if ((length - j) + cur_len < this->kmerlen) { // not enough coordinates - need to look at the next exon
                    cur_coords._push_back((uint32_t) (cur_exon.start + j), (uint32_t) cur_exon.end);
                    cur_len += ((cur_exon.end+1) - (cur_exon.start + j));
                    break;
                } else {
                    cur_coords._push_back((uint32_t) (cur_exon.start + j),(uint32_t) (cur_exon.start + j + this->kmerlen));
                    this->kmer_coords_exist = this->kmer_coords.insert(cur_coords);
                    if (this->kmer_coords_exist.second) { // was successfully inserted
                        // get sequence
                        sub_seq = rec.seq_.substr((cur_exon.start + j) - 1, this->kmerlen);
                        if (sub_seq.length() != this->kmerlen) {
                            std::cerr << "2: " << sub_seq.length() << " " << p_trans.getID() << " ";
                            for (int y = 0; y < cur_coords._getSize(); y++) {
                                std::cerr << cur_coords._getStart(y) << "-" << cur_coords._getEnd(y) << ";";
                            }
                            std::cerr << std::endl;
                        }
                        std::transform(sub_seq.begin(), sub_seq.end(), sub_seq.begin(), ::toupper);
                        this->trans_map._insert(sub_seq, cur_coords);
                        // do reverse_complement
                        char *rev_cur_kmer = new char[this->kmerlen];
                        strcpy(rev_cur_kmer, sub_seq.c_str());
                        reverseComplement(rev_cur_kmer, this->kmerlen);
//                        cur_coords.flipStrand(); // change strand information for the reverse complement
                        this->trans_map._insert(std::string(rev_cur_kmer), cur_coords);
//                        cur_coords.flipStrand(); // reset the strand information back
                        delete[] rev_cur_kmer;
                        // add new coordinates here
                    }
                }
                // since we went into this conditional, that means no previously encountered exons are left
                // and we can reset some things
                cur_coords._erase(); // need to test
                cur_coords = EVec(contigID, (uint8_t) p_trans.strand);
                cur_len = 0;
            }
        }
        exon_seq += rec.seq_.substr(cur_exon.start - 1, length);
    }
}

void HDB::save_db_info(){ // this method saves the general info about the database such as kmerlen, num kmers in transcriptome, num kmers in genome etc.
    std::string dbinfo_fname(this->out_fname);
    dbinfo_fname.append("/db.info");
    std::ofstream dbinfo_fp(dbinfo_fname.c_str());

    dbinfo_fp<<"kmerlen"<<"\t"<<(int)this->kmerlen<<std::endl;
    dbinfo_fp<<"genome"<<"\t"<<this->genome_fname_<<std::endl;
    dbinfo_fp<<"num trans kmer"<<"\t"<<(int)this->trans_map._getNumKmers()<<std::endl;
    dbinfo_fp<<"num genome kmer"<<"\t"<<(int)this->genom_map.size()<<std::endl;

    dbinfo_fp.close();
}

void HDB::save_contig_info(){ // this method saves the map of contigs to IDs
    std::string contiginfo_fname(this->out_fname);
    contiginfo_fname.append("/contig.info");
    std::ofstream contiginfo_fp(contiginfo_fname.c_str());

    for(auto it: this->id_to_contig){
        contiginfo_fp<<(int)it.first<<"\t"<<it.second.first<<"\t"<<it.second.second<<std::endl;
    }

    contiginfo_fp.close();
}

void HDB::save_trans_db(){

    std::string trans_fname(this->out_fname);
    trans_fname.append("/db.kmer.trans");
    std::ofstream trans_fp(trans_fname.c_str());

    trans_fp<<this->trans_map;

    trans_fp.close();
}

void HDB::save_genom_db(){
    std::string genom_fname(this->out_fname);
    genom_fname.append("/db.kmer.genom");
    std::ofstream genom_fp(genom_fname.c_str());

    uint8_t chrID,strand;
    uint32_t pos;

    for(auto it: this->genom_map){
        genom_fp<<it.first<<'\t';
        sort(it.second.begin(),it.second.end());
        for(auto vit: it.second) {
            std::tie(chrID,strand,pos) = vit;
            genom_fp << (int) chrID << ":" << (int) strand << "@" << pos<<";";
        }
        genom_fp<<std::endl;
    }
    genom_fp.close();
}

void HDB::load_genom_db(std::ifstream& stream) {
    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss(""), sub_ss("");
    std::string kmer,pos_string,chrID,strand,pos;
    while(std::getline(stream,line)) { // iterate over all kmers
        GenVec gv;
        ss.str(line);
        ss.clear();

        std::getline(ss, kmer, '\t');
        std::getline(ss, pos_string, '\t');

        ss.str(pos_string);
        ss.clear();

        while(std::getline(ss,pos_string,';')){ // iterate over all positions for a given kmer
            sub_ss.str(pos_string);
            sub_ss.clear();
            std::getline(sub_ss,chrID,':');
            std::getline(sub_ss,strand,'@');
            std::getline(sub_ss,pos,';');
            gv.push_back(std::make_tuple((uint8_t)std::stoi(chrID),(uint8_t)std::stoi(strand),(uint32_t)std::stoi(pos)));
        }
        this->genom_map.insert(std::make_pair(kmer,gv));
    }
}

void HDB::load_trans_db(std::ifstream& stream) {
    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss(""), sub_ss(""), sub_sub_ss("");
    std::string kmer,pos_string,chrID,strand, pos_vec,start,end;
    while(std::getline(stream,line)) { // iterate over all kmers
        std::vector<EVec> vev;
        ss.str(line);
        ss.clear();

        std::getline(ss, kmer, '\t');
        std::getline(ss, pos_string, '\t');

        ss.str(pos_string);
        ss.clear();

        while(std::getline(ss,pos_string,';')){ // iterate over all positions for a given kmer
            sub_ss.str(pos_string);
            sub_ss.clear();
            std::getline(sub_ss,chrID,':');
            std::getline(sub_ss,strand,'@');
            std::getline(sub_ss,pos_vec,';');

            EVec ev((uint8_t)std::stoi(chrID),(uint8_t)std::stoi(strand));
            sub_ss.str(pos_vec);
            sub_ss.clear();
            while(std::getline(sub_ss,pos_vec,',')){
                sub_sub_ss.str(pos_vec);
                sub_sub_ss.clear();
                std::getline(sub_sub_ss,start,'-');
                std::getline(sub_sub_ss,end,',');
                ev._push_back((uint32_t)std::stoi(start),(uint32_t)std::stoi(end));
            }
            vev.push_back(ev);
        }
        this->trans_map._insert(kmer,vev);
    }
}

// This function preloads all the contigs into memory for efficient lookup of splice junctions
void HDB::load_contigs(){
    FastaReader fastaReader(this->getGenomeFname());
    FastaRecord cur_contig;

    while (fastaReader.good()) {
        fastaReader.next(cur_contig);

        uint8_t contigID = this->getIDFromContig(cur_contig.id_); //this is the id of the new contig
        this->contigs.insert(std::make_pair(contigID,cur_contig.seq_));
    }
}

void HDB::load_db_info(std::ifstream& stream){
    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss("");
    std::string kmerlen,prename,genome_fn;
    while(std::getline(stream,line)) {
        ss.str(line);
        ss.clear();

        std::getline(ss, prename, '\t');
        if (std::strcmp(prename.c_str(), "kmerlen") == 0) {
            std::getline(ss, kmerlen, '\t');
            this->kmerlen=std::stoi(kmerlen);
        }
        if (std::strcmp(prename.c_str(), "genome") == 0) { // save the genome filename
            std::getline(ss, genome_fn, '\t');
            this->genome_fname_=genome_fn;
            // now test that the file still exists
            std::ifstream genome_fp;
            genome_fp.open(this->genome_fname_.c_str(),std::ios::in);
            if(!genome_fp.good()){
                std::cerr<<"genome fasta no longer can be found"<<std::endl;
                exit(-1);
            }
            genome_fp.close();
        }
        if(!genome_fn.empty() && !kmerlen.empty()){
            return;
        }
    }
    std::cerr<<"kmerlen unknown from the databse"<<std::endl;
    exit(1);
}

void HDB::load_contig_info(std::ifstream& stream){

    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss("");
    std::string contigName,contigID,contigLength;

    while (std::getline(stream,line)){
        ss.str(line);
        ss.clear();

        std::getline(ss,contigID,'\t');
        std::getline(ss,contigName,'\t');
        std::getline(ss,contigLength,'\t');
        this->contig_to_id.insert(std::make_pair(contigName,std::stoi(contigID)));
        this->id_to_contig.insert(std::make_pair(std::stoi(contigID),std::make_pair(contigName,std::stoi(contigLength))));
    }
}

void HDB::save_db(){
    if(this->transcriptomeBuild){
        std::cout<<"saving transcriptome"<<std::endl;
        HDB::save_trans_db();
    }
    HDB::save_genom_db();
    HDB::save_db_info();
    HDB::save_contig_info();
}

void HDB::load_db(std::string db_fname_base){
    // first check that everything is there
    // then load the componenets
    if(db_fname_base.rfind('/')==db_fname_base.length()-1){
        db_fname_base.pop_back();
    }
    std::string genom_fname(db_fname_base);
    genom_fname.append("/db.kmer.genom");
    std::string trans_fname(db_fname_base);
    trans_fname.append("/db.kmer.trans");
    std::string dbinfo_fname(db_fname_base);
    dbinfo_fname.append("/db.info");
    std::string contiginfo_fname(db_fname_base);
    contiginfo_fname.append("/contig.info");

    std::ifstream trans_fp;
    trans_fp.open(trans_fname.c_str(),std::ios::in);
    if(!trans_fp.good()){
        this->transcriptomeBuild=false;
    }
    if(this->transcriptomeBuild){
        HDB::load_trans_db(trans_fp);
    }
    trans_fp.close();

    std::ifstream genom_fp;
    genom_fp.open(genom_fname.c_str(),std::ios::in);
    if(!genom_fp.good()){
        std::cerr<<"FATAL: Couldn't open genome kmer data: "<<genom_fname<<std::endl;
        exit(1);
    }
    HDB::load_genom_db(genom_fp);
    genom_fp.close();

    std::ifstream dbinfo_fp;
    dbinfo_fp.open(dbinfo_fname.c_str(),std::ios::in);
    if(!dbinfo_fp.good()){
        std::cerr<<"FATAL: Couldn't open database info file: "<<dbinfo_fname<<std::endl;
        exit(1);
    }
    HDB::load_db_info(dbinfo_fp);
    dbinfo_fp.close();

    std::ifstream contiginfo_fp;
    contiginfo_fp.open(contiginfo_fname.c_str(),std::ios::in);
    if(!contiginfo_fp.good()){
        std::cerr<<"FATAL: Couldn't open contig info file: "<<contiginfo_fname<<std::endl;
        exit(1);
    }
    HDB::load_contig_info(contiginfo_fp);
    contiginfo_fp.close();

    this->load_contigs(); // preload current contigs into memory for efficient lookup
}

int HDB::getKmerLen() {
    return this->kmerlen;
}

MinMap::iterator HDB::find_trans(std::string &kmer) {
    return this->trans_map._find(kmer);
}

HDB::GenMap::iterator HDB::find_genom(std::string& kmer){
    return this->genom_map.find(kmer);
}

void HDB::generate_sq(std::ofstream &fp) {
    for (const auto& contig_it : this->id_to_contig){
        fp << "@SQ\tSN:" << contig_it.second.first << "\tLN:" << contig_it.second.second <<std::endl;
    }
}
