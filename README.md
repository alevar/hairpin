# hairpin
Fast and reference-free novel isofom assembly and quantification.

Hairpin is a tool that leverages the speed of pseudo-alignment techniques via exact-matching and dyrected acyclic graph representation of all potential mappings of reads in the transcritome landscape to accurately resolve ambiguitites via assesment of the local structure within the graph.

# Installation
version 3.13 of CMAKE required\
support for c++14 standard required\
`cmake ./`\
`make`

Due to the short development cycle, there are several memory allocation and access bugs still present in the code. As such, certain combinations of parameters and datasets might result in Segmentation Faults. Most of such mistakes can be avoided by changing paramters such as the kmer size.

# Building index
`mkdir -p drosophila_4_genomeDB`\
`hairpin build -o drosophila_4_genomeDB/ -r ./data/drosophila_4.fasta -k 22`

# Performin Pseudo-alignment
`mkdir -p drosophila_4_out`\
`hairpin quant -x ./drosophila_4_genomeDB/ -o ./drosophila_4_out/output -u ./data/reads_drosophila_4_pos_6_1_151/sample_01_1.fasta -s`

# A note on SegFaults
Please try a different kmer length - the code is under very active development and in 20 days a large project contains many bugs to be detected and fixed with a less hasty implementation.
