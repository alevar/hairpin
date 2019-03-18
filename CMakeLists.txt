cmake_minimum_required(VERSION 3.13)
project(hairpin)

set(CMAKE_CXX_STANDARD 11)

add_executable(hairpin main.cpp src/HGraph.cpp src/HGraph.h src/HDB.cpp src/HDB.h src/MinMap.cpp src/MinMap.h src/arg_parse.cpp src/arg_parse.h src/codons.h src/codons.cpp src/FastaTools.cpp src/FastaTools.h src/GBase.cpp src/GBase.h src/gdna.cpp src/gdna.cpp src/GFaSeqGet.cpp src/GFaSeqGet.h src/gff.cpp src/gff.h src/GHash.hh src/GList.hh src/GVec.hh src/tokenize.cpp src/tokenize.h)