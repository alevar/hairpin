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

    HDB hdb(args_build.get_string(Opt_Build::REF),args_build.get_string(Opt_Build::GFF));

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

//    gffMapper(args_align.get_string(Opt::TLST_FP),args_align.get_string(Opt::IN_AL),args_align.get_string(Opt::MULTI),args_align.get_string(Opt::GLST_FP),args_align.get_flag(Opt::MULTI_FLAG));
    return 0;
}

int main(int argc, char* argv[]) {

    if(strcmp(argv[1],"build") == 0){
        std::cout<<"building index"<<std::endl;
    }
    else if(strcmp(argv[1],"quant") ==0 ){
        std::cout<<"quantifying"<<std::endl;
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