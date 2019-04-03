#include <iostream>
#include <string>
#include <cstring>

#include "src/arg_parse.h"
#include "src/HDB.h"
#include "src/HGraph.h"

void print_help(){
    std::cout<<"help page"<<std::endl;
}

void process_reads_single(std::string readsFP, HGraph& hg){
    FastaReader fastaReader(readsFP);
    FastaRecord read;

    while (fastaReader.good()) {
        fastaReader.next(read);
        hg.add_read(read.seq_);
    }
}

int hairpin_quant(int argc,char* argv[]){
    enum Opt_Quant {HDB_FP= 'x',
        OUTPUT= 'o',
        READ1 = '1',
        READ2 = '2',
        UNPAIR= 'u',
        REF   = 'r',
        MAX   = 'm',
        MIN   = 'i',
        MM    = 'l'};

    ArgParse args_quant("hairpin_align");
    args_quant.add_string(Opt_Quant::HDB_FP,"hdb","","");
    args_quant.add_string(Opt_Quant::OUTPUT,"output","","");
    args_quant.add_string(Opt_Quant::READ1,"input1","","");
    args_quant.add_string(Opt_Quant::READ2,"input2","","");
    args_quant.add_string(Opt_Quant::UNPAIR,"unpaired","","");
    args_quant.add_string(Opt_Quant::REF,"reference","","");
    args_quant.add_int(Opt_Quant::MAX,"max_intron",30000,"");
    args_quant.add_int(Opt_Quant::MIN,"min_intron",20,"");
    args_quant.add_int(Opt_Quant::MM,"min_mismatch",10,"");

    args_quant.parse_args(argc,argv);

    HDB hdb;
    std::cerr<<"Loading the database"<<std::endl;
    hdb.load_db(args_quant.get_string(Opt_Quant::HDB_FP));

    HGraph hg(&hdb,args_quant.get_int(Opt_Quant::MAX),args_quant.get_int(Opt_Quant::MIN),args_quant.get_int(Opt_Quant::MM),args_quant.get_string(Opt_Quant::OUTPUT));
    std::cerr<<"processing reads"<<std::endl;
    process_reads_single(args_quant.get_string(Opt_Quant::UNPAIR),hg);
    std::cerr<<"parsing the graph"<<std::endl;
    hg.parse_graph();
//    hg.write_intron_gff();
    hg.print_stats();

    // when parsing a read - need to set the minimum number of kmers that need ot be mapped from that read
    // if fewer than n reads are mapped - remove any additions to the graph
    // to do so efficiently - for each read keep pointer to each the node where a kmer has been inserted and the associated edges
    // and remove them if fewer than n kmers have been matched from a read.

    return 0;
}

int hairpin_build(int argc,char* argv[]){
    enum Opt_Build {HDB_FP   = 'o',
                    REF      = 'r',
                    GFF      = 'a',
                    KMERLEN  = 'k'};

    ArgParse args_build("hairpin_build");
    args_build.add_string(Opt_Build::HDB_FP,"hdb","","");
    args_build.add_string(Opt_Build::REF,"reference","","");
    args_build.add_string(Opt_Build::GFF,"annotation","","");
    args_build.add_int(Opt_Build::KMERLEN,"kmerlen",31,"");

    args_build.parse_args(argc,argv);

    HDB hdb;
    if(args_build.get_string(Opt_Build::GFF).length()==0){
        hdb.init(args_build.get_string(Opt_Build::REF));
    }
    else {
        hdb.init(args_build.get_string(Opt_Build::GFF), args_build.get_string(Opt_Build::REF));
    }
    std::cerr<<"building the database:\t"<<std::endl;
    hdb.make_db(args_build.get_string(Opt_Build::HDB_FP), args_build.get_int(Opt_Build::KMERLEN));
    std::cerr<<"saving the database:\t"<<std::endl;
    hdb.save_db();

    return 0;
}

int main(int argc, char* argv[]) {

    if(strcmp(argv[1],"build") == 0){
        std::cerr<<"building index"<<std::endl;
        int argc_build=argc-1;
        char* argv_build[argc_build];
        memcpy(argv_build, argv+1, argc_build*sizeof(char*));
        hairpin_build(argc_build,argv_build);
    }
    else if(strcmp(argv[1],"quant") ==0 ){
        std::cerr<<"quantifying"<<std::endl;
        int argc_quant=argc-1;
        char* argv_quant[argc_quant];
        memcpy(argv_quant, argv+1, argc_quant*sizeof(char*));
        hairpin_quant(argc_quant,argv_quant);
    }
    else if (strcmp(argv[1],"help") == 0 || strcmp(argv[1],"--help") == 0){
        print_help();
    }
    else{
        std::cout<<"Unrecognized Mode: please consult the help manual"<<std::endl;
        print_help();
    }

    return 0;
}