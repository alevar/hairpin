# hairpin
Fast and reference-free novel isofom assembly and quantification

# Installation
version 3.13 of CMAKE required
support for c++14 standard required
`cmake ./`
`make`

In our testing, the current development version does not work unless compiled within CLion. Parser branch should contains a slightly outdated version, but more suitable for testing in a separate environment.

# Building index
`mkdir -p drosophila_4_genomeDB`
`hairpin build -o drosophila_4_genomeDB/ -r ./data/drosophila_4.fasta -k 22`

# Performin Pseudo-alignment
`mkdir -p drosophila_4_out`
`hairpin quant -x ./drosophila_4_genomeDB/ -o ./drosophila_4_out/output -u ./data/reads_drosophila_4_pos_6_1_151/sample_01_1.fasta -s`

# A note on SegFaults
Please try a different kmer length - the code is under very active development and in 20 days a large project contains many bugs to be detected and fixed with a less hasty implementation