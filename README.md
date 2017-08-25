**************************************************
**************************************************

Find nearest neighbor to target taxon: see http://localhost/mediawiki/index.php/Mash

**************************************************
**************************************************

Simulate reads with ART, example:

/nbacc/local/packages/art_bin_MountRainier/art_illumina -ss HS25 -i B_anthracis_Ames.fa -p -l 250 -f 10 -m 868 -s 408 -o ba_ART_sim -1 /nbacc/local/packages/art_bin_MountRainier/Illumina_profiles/HiSeq2500L250R1.txt -2 /nbacc/local/packages/art_bin_MountRainier/Illumina_profiles/HiSeq2500L250R2.txt

-p = paired-end
-l = length of reads to simulate (250 bp)
-f = fold-coverage (10x)
-m = mean fragment length
-s = fragment length standard deviation
-o = output prefix
-1 / -2 = custom quality profiles for R1 and R2

**************************************************
**************************************************

Run BLAST jobs on cluster, example:

dist_blast.pl cbot_ART_sim.fasta cbot_ART_sim 100 -task blastn

- first argument is path to reads in FASTA format
- second argument is the output prefix
- third argument is the number of files to break the input into (hence the number of jobs to run)
- last arguments are parameters to be passed to BLAST

**************************************************
**************************************************

Run parameter sweep, example:

varyParameters.pl -c 0_100_10 -i 0_100_10 -b 6_8_1 "qsub -V -cwd -b y blastValidate.pl" -x -t 1392 <fasta> <blast results>

OR

varyParameters.pl -e -2,-4,-8,-16,-24 -b 5,6,7 "qsub -V -cwd -b y blastValidate.pl" -x -t 632 yp_ART_sim.fasta yp_ART_sim

(the quotes ensure the varied parameters don't get passed directly to qsub)

**************************************************
**************************************************

Evaluate parameter sweep, example:

scripts/evaluate_parameter_sweep.pl 1392 /nbacc/project/fy16/BLAST_validation_experiments/ART/B_anthracis /nbacc/project/fy16/BLAST_validation_experiments/ART/B_cereus

**************************************************
**************************************************