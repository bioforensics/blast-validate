# BLAST-based validation of metagenomic sequence assignments

----
## Publication

Preprint available at <https://doi.org/10.1101/181636>

> Bazinet, A. L., Ondov, B. D., Sommer, D. D., & Ratnayake, S. (2017). BLAST-based validation of metagenomic sequence assignments. [_bioRxiv_:181636](https://doi.org/10.1101/181636).


----
## Installation
1. git clone --recurse-submodules https://github.com/bioforensics/blast-validate.git

2. Install [KronaTools](https://github.com/marbl/Krona/wiki/KronaTools) (located in blast-validate/scripts/Krona/KronaTools)
```
cd blast-validate/scripts/Krona/KronaTools/
./install.pl --prefix <install path>
./updateTaxonomy.sh
```
see krona wiki link above for more detailed instructions

3. Install [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

4. Optional: download and install the [ART read simulator](https://omictools.com/art-tool) (needed for parameter optimization)


----
## Parameter optimization
### (example using *B. anthracis*)

1.) Simulate reads with ART.

Example:

     art_illumina -ss HS25 -i B_anthracis_Ames.fa -p -l 250 -f 10 -m 868 -s 408 -o ba_ART_sim -1 /usr/local/packages/art_bin_MountRainier/Illumina_profiles/HiSeq2500L250R1.txt -2 /usr/local/packages/art_bin_MountRainier/Illumina_profiles/HiSeq2500L250R2.txt

     -p = paired-end  
     -l = length of reads to simulate (250 bp)  
     -f = fold-coverage (10x)  
     -m = mean fragment length  
     -s = fragment length standard deviation  
     -o = output prefix  
     -1 / -2 = custom quality profiles for R1 and R2

**************************************************

2.) Convert FASTQ files to FASTA format and combine R1 and R2 into a single file.

Example:

     perl scripts/convertFastqToFasta.pl < ba_ART_sim1.fq > ba_ART_sim1.fa; perl scripts/convertFastqToFasta.pl < ba_ART_sim2.fq > ba_ART_sim2.fa; cat ba_ART_sim1.fa ba_ART_sim2.fa > ba_ART_sim.fasta

**************************************************

3.) BLAST all reads against the NCBI nt database.

Example:

     blastn -query ba_ART_sim.fasta -task blastn -db /your/path/to/NCBI/nt -outfmt 7 > ba_ART_sim.blast

**************************************************

4.) Run parameter sweep.

Example:

     perl scripts/blastValidate.pl -e 0,-1,-2,-4,-8,-16,-32,-64,-128 -b 0,1,2,4,8,16,32,64,128 -x -t 1392 ba_ART_sim.fasta ba_ART_sim.blast

**************************************************

5.) (Optional) Perform steps 1-4 for a near neighbor genome (e.g., a *B. cereus* genome)

**************************************************

6.) Evaluate results of parameter sweep.

Example:

     perl scripts/evaluate_parameter_sweep.pl 1392 /path/to/B_anthracis/LCA/files [/path/to/near neighbor (B. cereus)/LCA/files]


----
## Validate reads provisionally assigned to a target taxon
### (example using NYC subway data and *B. anthracis*)

1.) BLAST reads against the NCBI nt database.

Example:

     blastn -query subway_analysis/SRR1748708_Kraken_Ba.fasta -task blastn -db /your/path/to/NCBI/nt -outfmt 7 > subway_analysis/SRR1748708_Kraken_Ba.blast

**************************************************

2.) Run the BLAST-validate script using parameter values previously determined to be optimal.

Example:

     perl scripts/blastValidate.pl -e -64 -b 8 -x -t 1392 subway_analysis/SRR1748708_Kraken_Ba.fasta subway_analysis/SRR1748708_Kraken_Ba.blast
