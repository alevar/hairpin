#include <iostream>
#include <string>
#include <cstring>

#include "src/arg_parse.h"
#include "src/HDB.h"

void print_help(){
    std::cout<<"help page"<<std::endl;
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

    HDB hdb(args_build.get_string(Opt_Build::GFF),args_build.get_string(Opt_Build::REF));
    std::cout<<"building the database"<<std::endl;
    hdb.make_db(args_build.get_string(Opt_Build::HDB_FP), args_build.get_int(Opt_Build::KMERLEN));
    std::cout<<"saving the database"<<std::endl;
    hdb.save_trans_db();
    hdb.save_genom_db();
    hdb.save_db_info();
    hdb.save_contig_info();

    return 0;
}

int hairpin_align(int argc,char* argv[]){
    enum Opt_Align {HDB   = 'x',
                    OUTPUT= 'o',
                    READ1 = '1',
                    READ2 = '2',
                    UNPAIR= 'u'};

    ArgParse args_align("hairpin_align");
    args_align.add_string(Opt_Align::HDB,"hdb","","");
    args_align.add_string(Opt_Align::OUTPUT,"output","","");
    args_align.add_string(Opt_Align::READ1,"input1","","");
    args_align.add_string(Opt_Align::READ2,"input2","","");
    args_align.add_string(Opt_Align::UNPAIR,"unpaired","","");

    args_align.parse_args(argc,argv);

    // when parsing a read - need to set the minimum number of kmers that need ot be mapped from that read
    // if fewer than n reads are mapped - remove any additions to the graph
    // to do so efficiently - for each read keep pointer to each the node where a kmer has been inserted and the associated edges
    // and remove them if fewer than n kmers have been matched from a read.

    return 0;
}

int main(int argc, char* argv[]) {

    if(strcmp(argv[1],"build") == 0){
        std::cout<<"building index"<<std::endl;
        int argc_build=argc-1;
        char* argv_build[argc_build];
        memcpy(argv_build, argv+1, argc_build*sizeof(char*));
        hairpin_build(argc_build,argv_build);
    }
    else if(strcmp(argv[1],"quant") ==0 ){
        std::cout<<"quantifying"<<std::endl;
        int argc_quant=argc-1;
        char* argv_quant[argc_quant];
        memcpy(argv_quant, argv+1, argc_quant*sizeof(char*));
        hairpin_build(argc_quant,argv_quant);
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