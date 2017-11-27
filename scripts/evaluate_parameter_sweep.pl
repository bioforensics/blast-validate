#!/usr/bin/perl
#
# This script evaluates the results of a blastValidate parameter sweep using either the "false negatives" approach only, 
# or additionally the "near neighbor" approach.
#
# As input, takes 1) the NCBI taxonomy ID of the "target taxon" (the taxon whose classification is being evaluated);
#                 2) the absolute path to a directory of target taxon LCA files (produced by the blastValidate script); and
#                 3) optionally, the path to a directory of near neighbor LCA files (also produced by the blastValidate script)
#
# Currently this script uses the following optimality criterion:
# first maximize precision, then maximize sensitivity for parameter value combinations yielding sensitivity > 0.
#
# For context, see the following publication:
# "BLAST-based validation of metagenomic sequence assignments"
#
# usage: perl evaluate_parameter_sweep.pl taxon-id absolute-path-to-taxon-lca-files [absolute-path-to-near-neighbor-lca-files]
#

use lib (`ktGetLibPath`);
use KronaTools;

if(@ARGV < 2 || @ARGV > 3) {
    print "usage: perl evaluate_parameter_sweep.pl taxon-id absolute-path-to-taxon-lca-files [absolute-path-to-near-neighbor-lca-files]\n";
    exit;
}

$taxonid = $ARGV[0];
$taxon_lcafiles_path = $ARGV[1];
$near_neighbor_lcafiles_path = "";
if(@ARGV == 3) {
    $near_neighbor_lcafiles_path = $ARGV[2];
}

%taxon_lcafiles = ();
%near_neighbor_lcafiles = ();

# for determining "optimal" results
$optimalFNstring = "\n";
$optimalNNstring = "\n";
$FNprec = 0.0;
$FNsens = 0.0;
$FNeval = 100;
$FNbitscore = 0;
$NNprec = 0.0;
$NNsens = 0.0;
$NNeval = 100;
$NNbitscore = 0;

# store lca files for taxon
foreach $lcafile (glob ($taxon_lcafiles_path . "/*.lca")) {
    $taxon_lcafiles{$lcafile} = 1;
}

if($near_neighbor_lcafiles_path ne "") {
    # store lca files for near neighbor
    foreach $lcafile (glob ($near_neighbor_lcafiles_path . "/*.lca")) {
	$near_neighbor_lcafiles{$lcafile} = 1;
    }
}

# write results to file
$outputfile = "param_sweep_eval";
open OUTPUT, ">$outputfile";

if($near_neighbor_lcafiles_path ne "") { # we will have an extra column for precision based on false positives
    print OUTPUT "TP\tFN\tFP\tsens (TP/read count)\tprec (TP/(TP+FN))\tprec (TP/(TP+FP))\tminlen\tminid\tmaxeval\tminbitscorediff\n";
} else {
    print OUTPUT "TP\tFN\tsens (TP/read count)\tprec (TP/(TP+FN))\tminlen\tminid\tmaxeval\tminbitscorediff\n";
}

# loop through lca files produced by blastValidate script
# example file name: ba_ART_sim.fasta.c_i_l0_b16_e-80_x.lca
foreach $lcafile (keys %taxon_lcafiles) {
    if($lcafile =~ /.*\/(.*)\.c([0-9]*?)_i([0-9]*?)_l([0-9]*?)_b([0-9]*)/) {
	$prefix = $1;
	$coverage = $2;
	$identity = $3;
	$length = $4;
	$bitscorediff = $5;
	$evalue = "";

	if($lcafile =~ /.*b[0-9]*_e(.*?)_/ || $lcafile =~ /.*b[0-9]*_e(.*?)\./) {
	    $evalue = $1;
	}

	# make sure corresponding near neighbor lca file exists, if it's supposed to
	$nnlcafile = "";
	if($near_neighbor_lcafiles_path ne "") {
	    $nnprefix = "";
	    foreach $nnlcafiletemp (keys %near_neighbor_lcafiles) {
		if($nnlcafiletemp =~ /(.*\/.*)\.c.*/) {
		    $nnprefix = $1;
		} else {
		    print "couldn't get near neighbor lca file prefix from file: $nnlcafiletemp  exiting...\n";
		    exit;
		}
	    }
	    $nnlcafile = $nnprefix . ".c" . $coverage . "_i" . $identity . "_l" . $length . "_b" . $bitscorediff;
	    if($evalue ne "") {
		$nnlcafile .= "_e" . $evalue;
	    }

	    $foundnnlcafile = 0;
	    $nnlcafile1 = $nnlcafile . "_x.lca";
	    $nnlcafile2 = $nnlcafile . ".lca";
	    foreach $nnlcafiletemp (keys %near_neighbor_lcafiles) {
		if($nnlcafiletemp =~ /$nnlcafile1/ || $nnlcafiletemp =~ /$nnlcafile2/) {
		    $nnlcafile = $nnlcafiletemp;
		    $foundnnlcafile = 1;
		    last;
		}
	    }

	    if($foundnnlcafile == 0) {
		print "didn't find near neighbor lca file with prefix: $nnlcafile  skipping...\n";
		next;
	    }
	}

	# parse lca file(s)
	# example: # query	taxon	name	kingdom	best-hit-subj-len
	# NC_003997.3-348480/1	1	root	(root)

	$read_count = 0;
	$taxon_count = 0;
	$taxon_lineage_count = 0;
	$no_hits_count = 0;
	$taxon_count_in_nn = 0;

	open LCA, $lcafile;
	print "lca file: $lcafile\n";
	while(<LCA>) {
	    chomp;
	    $line = $_;
	    if(($line =~ /(.*?)\t(.*?)\t(.*?)\t.*$/ || $line =~ /(.*?)\t(.*?)\t(.*)$/) && $line !~ /query\ttaxon\tname/) {
		$query = $1;
		$taxon = $2;
		$taxname = $3;
		chomp($query, $taxon, $taxname);

		$read_count++;
		if($taxon eq $taxonid || taxContains($taxonid, $taxon)) { # adding the second clause in case taxon was not provided to blastValidate script
		    $taxon_count++;
		} elsif(taxContains($taxon, $taxonid)) {
		    $taxon_lineage_count++;
		} elsif($taxon eq "-") {
		    $no_hits_count++;
		} else {
		    #print "false negative: $line\n";
		}
	    }
	}
	close LCA;

	if($nnlcafile ne "") {
	    open LCA, $nnlcafile;
	    while(<LCA>) {
		chomp;
		$line = $_;
		if(($line =~ /(.*?)\t(.*?)\t(.*?)\t.*$/ || $line =~ /(.*?)\t(.*?)\t(.*)$/) && $line !~ /query\ttaxon\tname/) {
		    $query = $1;
		    $taxon = $2;
		    $taxname = $3;
		    chomp($query, $taxon, $taxname);

		    if($taxon eq $taxonid || taxContains($taxonid, $taxon)) { # adding the second clause in case taxon was not provided to blastValidate script
			$taxon_count_in_nn++;
		    }
		}
	    }
	    close LCA;
	}

	$tp = $taxon_count;
	$fn = $read_count - $taxon_count - $taxon_lineage_count - $no_hits_count;
	chomp($tp);
	chomp($fn);

	if($coverage eq "") {
	    $coverage = 0;
	}

	if($identity eq "") {
	    $identity = 0;
	}

	if($evalue eq "") {
	    $evalue = "inf";
	}

	# defining sensitivity so that it's useful/not misleading
	if($tp > 0) {
	    $sensitivity = $tp / $read_count;
	} else {
	    $sensitivity = 0;
	}
	# (as opposed to defining sensitivity this more standard way)
	#if($tp + $fn > 0) {
	#$sensitivity = $tp / ($tp + $fn);
	#} else {
	#$sensitivity = 0;
	#}

	# precision based on false negatives can always be computed
	if($tp + $fn > 0) {
	    $precisionfn = $tp / ($tp + $fn);
	} else {
	    $precisionfn = 1;
	}
	
	if($nnlcafile ne "") {
	    $fp = $taxon_count_in_nn;
	    chomp($fp);
	    
	    # compute precision based on false positives (the "near neighbor" approach)
	    if($tp + $fp > 0) {
		$precisionfp = $tp / ($tp + $fp);
	    } else {
		$precisionfp = 1;
	    }

	    $outputstring = "$tp\t$fn\t$fp\t$sensitivity\t$precisionfn\t$precisionfp\t$length\t$identity\t$evalue\t$bitscorediff\n";
	    print OUTPUT $outputstring;

	    # debugging
	    #print "precisionfp: $precisionfp\n";
	    #print "sens: $sensitivity\n";
	    #print "NNsens: $NNsens\n";
	    #print "bitscorediff: $bitscorediff\n";
	    #print "NNbitscore: $NNbitscore\n";
	    #print "evalue: $evalue\n";
	    #print "NNeval: $NNeval\n";

	    # possibly update optimal near neighbor results
	    if($sensitivity > 0 && $precisionfp > $NNprec) {
		# print "better precision: $outputstring\n";
		$NNprec = $precisionfp;
		$NNsens = $sensitivity;
		$NNbitscore = $bitscorediff;
		$NNeval = $evalue;
		$optimalNNstring = $outputstring;
	    } elsif($sensitivity > 0 && $precisionfp == $NNprec && $sensitivity > $NNsens) {
		# print "better sensitivity: $outputstring\n";
		$NNsens = $sensitivity;
		$NNbitscore = $bitscorediff;
		$NNeval = $evalue;
		$optimalNNstring = $outputstring;
	    } elsif($sensitivity > 0 && $precisionfp == $NNprec && $sensitivity == $NNsens) {
		if($bitscorediff > $NNbitscore) {
		    # print "better bit score: $outputstring\n";
		    $NNbitscore = $bitscorediff;
		    $optimalNNstring = $outputstring;
		} elsif($bitscorediff == $NNbitscore && $evalue < $NNeval) {
		    # print "better evalue: $outputstring\n";
		    $NNeval = $evalue;
		    $optimalNNstring = $outputstring;
		}
	    }
	} else {
	    $outputstring = "$tp\t$fn\t$sensitivity\t$precisionfn\t$length\t$identity\t$evalue\t$bitscorediff\n";
	    print OUTPUT $outputstring;
	}

	# possibly update optimal false negatives results
	if($sensitivity > 0 && $precisionfn > $FNprec) {
	    $FNprec = $precisionfn;
	    $FNsens = $sensitivity;
	    $FNbitscore = $bitscorediff;
	    $FNeval = $evalue;
	    $optimalFNstring = $outputstring;
	} elsif($sensitivity > 0 && $precisionfn == $FNprec && $sensitivity > $FNsens) {
	    $FNsens = $sensitivity;
	    $FNbitscore = $bitscorediff;
	    $FNeval = $evalue;
	    $optimalFNstring = $outputstring;
	} elsif($sensitivity > 0 && $precisionfn == $FNprec && $sensitivity == $FNsens) {
	    if($bitscorediff > $FNbitscore) {
		$FNbitscore = $bitscorediff;
		$optimalFNstring = $outputstring;
	    } elsif($bitscorediff == $FNbitscore && $evalue < $FNeval) {
		$FNeval = $evalue;
		$optimalFNstring = $outputstring;
	    }
	}
    } else {
	print "could not parse lca file: $lcafile  exiting...\n";
	exit;
    }
}

close OUTPUT;

# print results for optimal parameter value combinations
# example:
# TP	FN	sens (TP/read count)	prec (TP/(TP+FN))	minlen	minid	maxeval	minbitscorediff
# 20272	0	0.104943831857949	1	0	0	-64	8

# check to see if we're using the near neighbor approach
if($near_neighbor_lcafiles_path ne "") {
    # print "optimal results for near neighbor approach: $optimalNNstring\n";
    @chunks = split("\t", $optimalNNstring);
    print "\noptimal results for near neighbor approach:\n";
    print "\ttrue positives: " . $chunks[0] . "\n";
    print "\tfalse negatives: " . $chunks[1] . "\n";
    print "\tfalse positives: " . $chunks[2] . "\n";
    print "\tsensitivity (TP/read count): " . $chunks[3] . "\n";
    print "\tprecision (TP/(TP+FN)): " . $chunks[4] . "\n";
    print "\tprecision (TP/(TP+FP)): " . $chunks[5] . "\n";
    print "\tminimum alignment length: " . $chunks[6] . "\n";
    print "\tminimum alignment percent identity: " . $chunks[7] . "\n";
    print "\tmaximum E-value: " . $chunks[8] . "\n";
    print "\tminimum bit score difference: " . $chunks[9] . "\n";
}

# we always have the false negatives approach
# print "optimal results for false negatives approach: $optimalFNstring\n";
@chunks = split("\t", $optimalFNstring);
if($near_neighbor_lcafiles_path ne "") {
    print "\noptimal results for false negatives approach:\n";
    print "\ttrue positives: " . $chunks[0] . "\n";
    print "\tfalse negatives: " . $chunks[1] . "\n";
    print "\tsensitivity (TP/read count): " . $chunks[3] . "\n";
    print "\tprecision (TP/(TP+FN)): " . $chunks[4] . "\n";
    print "\tminimum alignment length: " . $chunks[6] . "\n";
    print "\tminimum alignment percent identity: " . $chunks[7] . "\n";
    print "\tmaximum E-value: " . $chunks[8] . "\n";
    print "\tminimum bit score difference: " . $chunks[9] . "\n";
} else {
    print "\noptimal results for false negatives approach:\n";
    print "\ttrue positives: " . $chunks[0] . "\n";
    print "\tfalse negatives: " . $chunks[1] . "\n";
    print "\tsensitivity (TP/read count): " . $chunks[2] . "\n";
    print "\tprecision (TP/(TP+FN)): " . $chunks[3] . "\n";
    print "\tminimum alignment length: " . $chunks[4] . "\n";
    print "\tminimum alignment percent identity: " . $chunks[5] . "\n";
    print "\tmaximum E-value: " . $chunks[6] . "\n";
    print "\tminimum bit score difference: " . $chunks[7] . "\n";
}

exit;
