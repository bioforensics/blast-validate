#! /usr/bin/perl
#
# This script filters BLAST hits based on input parameters and then classifies the sequence
# using the LCA algorithm, reporting the target taxon if the classification is a child of the
# target taxon, or the exact classification otherwise.
#
# For context, see the following publication:
# "BLAST-based validation of metagenomic sequence assignments"
#
# Run script without any arguments to see detailed usage.
#

use strict;

use Getopt::Long;

use lib (`ktGetLibPath`);
use KronaTools;

my $minCov = 0;
my $minId = 0;
my $maxEval = -256;
my $minLength = 0;
my $taxId = 0;
my $bitWin = 0;
my $exclude = 0;
my $organism;

GetOptions
(
 "c=s" => \$minCov,
 "l=s" => \$minLength,
 "i=s" => \$minId,
 "t=i" => \$taxId,
 "b=s" => \$bitWin,
 "x" => \$exclude,
 "o=s" => \$organism,
 "e=s" => \$maxEval
);

if ( @ARGV < 1 )
{
    print STDERR "blastValidate.pl [options] <fasta> [<blast>]\n";
    print STDERR "   -o str  organism (\"ba\", \"cb\", \"yp\")\n";
    print STDERR "   -c %    min coverage [0] or <value1>,<value2>,...\n";
    print STDERR "   -l int  min alignment length or <value1>,<value2>,...\n";
    print STDERR "   -i %    min % ID [0]  or <value1>,<value2>,...\n";
    print STDERR "   -e num  log10 of max e-value [inf] <value1>,<value2>,...\n";
    print STDERR "   -t int  target taxon ID [0]\n";
    print STDERR "   -b num  bit score window for ties [0]  or <value1>,<value2>,...\n";
    print STDERR "   -x      exclude hits to the \"other sequences\" clade\n";
    print STDERR "\n";
    print STDERR "Filters BLAST hits based on input parameters and then classifies the sequence using\n";
    print STDERR "the LCA algorithm, reporting the target taxon if the classification is a child of\n";
    print STDERR "the target taxon, or the exact classification otherwise.\n";
    print STDERR "\n";
    print STDERR "Organism defaults (-c, -i, -b, and -e can be overridden):\n";
    print STDERR "   ba (B. anthracis): -c 0 -i 0 -b 8 -e -64 -x -t 1392\n";
    print STDERR "   cb (C. botulinum): -c 0 -i 0 -b 8 -e -64 -x -t 1491\n";
    print STDERR "   yp (Y. pestis):    -c 0 -i 0 -b 8 -e -64 -x -t 632\n";
    print STDERR "\n";
    print STDERR "   Defaults were chosen to maximize sensitivity while retaining perfect\n";
    print STDERR "   precision, using a ground truth from simulated reads.\n";
    print STDERR "\n";
    print STDERR "<blast>  defaults to <fasta>.blastn\n";
    print STDERR "Output is written to <fasta>.<params>.lca, where <params> describes parameters.\n";
    
    exit;
}

if ( $organism eq 'ba')
{
    if ( ! defined $minCov )
    {
	$minCov = 0;
    }
    
    if ( ! defined $minId )
    {
	$minId = 0;
    }
    
    if ( ! defined $bitWin )
    {
	$bitWin = 8;
    }

    if ( ! defined $maxEval )
    {
	$maxEval = -64;
    }
    
    $taxId = 1392;
    $exclude = 1;
}
elsif ( $organism eq 'cb')
{
    if ( ! defined $minCov )
    {
	$minCov = 0;
    }
    
    if ( ! defined $minId )
    {
	$minId = 0;
    }
    
    if ( ! defined $bitWin )
    {
	$bitWin = 8;
    }

    if ( ! defined $maxEval )
    {
	$maxEval = -64;
    }
    
    $taxId = 1491;
    $exclude = 1;
}
elsif ( $organism eq 'yp')
{
    if ( ! defined $minCov )
    {
	$minCov = 0;
    }
    
    if ( ! defined $minId )
    {
	$minId = 0;
    }
    
    if ( ! defined $bitWin )
    {
	$bitWin = 8;
    }

    if ( ! defined $maxEval )
    {
	$maxEval = -64;
    }
    
    $taxId = 632;
    $exclude = 1;
}
elsif ( defined $organism )
{
    print STDERR "ERROR: unrecognized organism (\"$organism\") given to -o\n";
    exit 1;
} else {
    if ( ! defined $minCov )
    {
	$minCov = 0;
    }
    
    if ( ! defined $minId )
    {
	$minId = 0;
    }
    
    if ( ! defined $bitWin )
    {
	$bitWin = 0;
    }

    if ( ! defined $maxEval )
    {
	$maxEval = -256;
    }
}

my ($fasta, $blast) = @ARGV;

if ( ! defined $blast )
{
    $blast = "$fasta.blastn";
}

my $params = "$fasta.c$minCov\_i$minId\_l$minLength\_b$bitWin" . (defined $maxEval ? "_e$maxEval" : "") . ($exclude ? "_x" : "");

my %lengths = getSequenceLengths($fasta);

print STDERR "Loading taxonomy...\n";

setOption('noRank', 1);

loadTaxonomy();

my %idsFiltered;
my %noHits;
my %minLen;

my $filtered = "$params.filt";

my $filter = ! -e "$filtered" || (stat($blast))[9] > (stat("$filtered"))[9];

print STDERR "Read in blast output...\n";
open BLAST, $blast or die $blast;
my @lines = <BLAST>;
close BLAST;

print STDERR " start search with $bitWin(bitwin), $minLength(len), $minCov(cov), $minId(id), $maxEval(eval) \n";
my @lengths = split /,/, $minLength;
foreach my $l ( @lengths ) {
    my @covs = split /,/, $minCov;
    foreach my $c ( @covs ) {
	my @ids = split /,/, $minId;
	foreach my $i ( @ids ) {
	    if (defined $maxEval && $maxEval ne '') { 
		my @evs = split /,/, $maxEval;
		foreach my $e ( @evs ) {
		    print STDERR " run blast_validate with $bitWin(bitwin), $l(len), $c(cov), $i(id), $e(eval) \n";
		    validate_blast(\@lines, $bitWin, $l, $c, $i, $e);
		}
	    } else {
		print STDERR " run blast_validate with $bitWin(bitwin), $l(len), $c(cov), $i(id), $maxEval(eval) \n";
		validate_blast(\@lines, $bitWin, $l, $c, $i, $maxEval);
	    }
	}
    }
}
print STDERR "done.\n";
exit;


## BEGIN SUBROUTINES ##

sub classifyBlastArr
{
    # taxonomically classifies BLAST results based on LCA (or random selection)
    # of 'best' hits.
    #
    # Options used: bitScore, factor, include, percentIdentity, random, score
    
    my # parameters
	(
	 $blastArr, # array with tabular BLAST results
	 
	 # hash refs to be populated with results (keyed by query ID)
	 #
	 $taxIDs,
	 $scores
	) = @_;
    
    # open BLAST, "<$fileName" or ktDie("Could not open $fileName\n");
    
    my $lastQueryID;
    my $topScore;
    my $topEVal;
    my $ties;
    my $taxID;
    my %lcaSet;
    my $totalScore;
    my $zeroEVal;
    
    for my $line ( @{$blastArr} )
    {
	
	chomp $line;
	# print STDERR $line . " threshold " . getOption('threshold') . " \n";
	
	if ( $line =~ /^#/ )
	{
	    if ( $line =~ /Query: ([\S]+)/ )
	    {
		# Initialize taxID and score in case this query has no hits
		$taxIDs->{$1} = -1;
		$scores->{$1} = 0;
	    }
	    next;
	}
	
	my
	    (
	     $queryID,
	     $hitID,
	     $identity,
	     $length,
	     $mismatches,
	     $gaps,
	     $queryStart,
	     $queryEnd,
	     $subjectStart,
	     $subjectEnd,
	     $eVal,
	     $bitScore
	    ) = split /\t/, $line;

	# print "query ID: $queryID  last query ID: $lastQueryID\n";
	
	if ( $queryID ne $lastQueryID )
	{
	    if
		(
		 ! defined $lastQueryID &&
		 ! defined $taxIDs->{$queryID} &&
		 getOption('include')
		)
	    {
		ktWarn("-i specified but array does not contain comment lines. Queries with no hits will not be included for this file."); 
	    }
	    
	    if (  $ties )
	    {
		# add the chosen hit from the last queryID
		
		if ( ! getOption('random') )
		{
		    $taxID = taxLowestCommonAncestor(keys %lcaSet)
		}
		
		$taxIDs->{$lastQueryID} = $taxID;
		$scores->{$lastQueryID} = $totalScore / $ties;
	    }
	    
	    $ties = 0;
	    $totalScore = 0;
	    %lcaSet = ();
	}
	
	if ( ! defined $hitID )
	{
	    last; # EOF
	}
	
	my $acc = getAccFromSeqID($hitID);
	
	if ( ! defined $acc )
	{
	    $lastQueryID = $queryID;
	    next;
	}
	
	if # this is a 'best' hit if...
	    (
	     $queryID ne $lastQueryID || # new query ID (including null at EOF)
	     $bitScore >= $topScore - getOption('threshold') || # within score threshold
	     getOption('factor') && $eVal <= getOption('factor') * $topEVal # within e-val factor
	    )
	{
	    # add score for average
	    #
	    if ( getOption('percentIdentity') )
	    {
		$totalScore += $identity;
	    }
	    elsif ( getOption('bitScore') )
	    {
		$totalScore += $bitScore;
	    }
	    else
	    {
		if ( $eVal > 0 )
		{
		    $totalScore += (log $eVal) / log 10;
		}
		else
		{
		    $totalScore += $KronaTools::minEVal;
		    $zeroEVal = 1;
		}
	    }

	    $ties++;
	    
	    if # use this hit if...
		(
		 ! getOption('random') || # using LCA
		 $queryID ne $lastQueryID || # new query ID
		 int(rand($ties)) == 0 # randomly chosen to replace other hit
		)
	    {
		my $newTaxID = getTaxIDFromAcc($acc);
		
		if ( ! $newTaxID || ! taxIDExists($newTaxID) )
		{
		    $newTaxID = 1;
		}
		
		if ( getOption('random') )
		{
		    $taxID = $newTaxID;
		}
		else
		{
		    $lcaSet{$newTaxID} = 1;
		}
	    }
	}
	
	if ( $queryID ne $lastQueryID )
	{
	    $topScore = $bitScore;
	    $topEVal = $eVal;
	}
	
	$lastQueryID = $queryID;
    }
    
    if ( $zeroEVal )
    {
	ktWarn("Array had e-values of 0. Approximated log[10] of 0 as $KronaTools::minEVal.");
    }
}	

# takes parsed blast and bitWin string
sub validate_blast {
    my # parameters
	(
	 $blast_lines, # array with tabular BLAST results
	 $bitWin,
	 $minLength,
	 $minCov,
	 $minId,
	 $maxEval,
	) = @_;

    my @filterArr;
    my @idOrder;
    print STDERR " filtering hits by minLength, minCov, minId...\n";
    my $queryIdLast;

    for my $line (@{$blast_lines}) {
	chomp $line;

	if ( $line =~ /^\#/ ) {
	    push @filterArr, $line;
	    
	    if ( $line =~ /Query: ([\S]+)/ ) {
		push @idOrder, $1;
	    }
	    elsif ( $line eq "# 0 hits found" ) {
		$noHits{$idOrder[-1]} = 1;
	    }
	    
	    next;
	}

	my
	    (
	     $queryId,
	     $hitId,
	     $identity,
	     $length,
	     $mismatches,
	     $gaps,
	     $queryStart,
	     $queryEnd,
	     $subjectStart,
	     $subjectEnd,
	     $eVal,
	     $bitScore,
	     $slen
	    ) = split /\t/, $line;

	if ( $idsFiltered{$queryId} )
	{
	    next;
	}

	if ( ! defined $minLen{$queryId} )
	{
	    $minLen{$queryId} = $slen;
	}

	if
	    (
	     $queryId eq $queryIdLast ||
	     $filter &&
	     $identity >= $minId &&
	     ($queryEnd - $queryStart + 1) / $lengths{$queryId} * 100 >= $minCov &&
	     (! defined $maxEval || log($eVal)/log(10) <= $maxEval) &&
	     $length >= $minLength
	    )
	{
	    my $taxon = getTaxIDFromAcc(getAccFromSeqID($hitId));
	    if ( ! $exclude || ! taxContains(28384, $taxId) )
	    {
		my $taxName = getTaxName($taxon);
		push @filterArr, "$line\t$taxName";
		
	    }
	}
	else
	{
	    $idsFiltered{$queryId} = 1;
	}
	
	$queryIdLast = $queryId;
    }

    print STDERR " applying LCA algorithm to hits...\n";
    my @bits = split /,/, $bitWin;
    foreach my $bit ( @bits ) {
	my %taxIds;
	my %scores;

	setOption('threshold', $bit);

	# push a null string onto the end of filterArr so Krona code (which previously read from a file) works properly
	push @filterArr, "";
	classifyBlastArr(\@filterArr, \%taxIds, \%scores);

	my $lca = "$fasta.c$minCov\_i$minId\_l$minLength\_b$bit\_e" . (defined $maxEval ? $maxEval : "none") . ".lca";

	#print STDERR "output $bit to $lca...\n";
	
	open LCA, ">$lca";
	
	print LCA "# Minimum query coverage:     $minCov\%\n";
	print LCA "# Minimum alignment length:   $minLength\n";
	print LCA "# Minimum identity:           $minId\%\n";
	print LCA "# Maximum e-value (log10):    " . (defined $maxEval ? $maxEval : "none") . "\n";
	print LCA "# Target taxon:               $taxId\n";
	print LCA "# LCA bit score range:        $bit\n";
	print LCA "# Exclude \"other sequences\":  " . ($exclude ? 'yes' : 'no') . "\n";
	print LCA "#\n";
	print LCA "# query	taxon	name	kingdom	best-hit-subj-len\n";
	
	foreach my $id ( @idOrder )
	{
	    if ( $noHits{$id} )
	    {
		print LCA "$id\t-\tno hits\n";
	    }
	    elsif ( $taxIds{$id} == -1 )
	    {
		#print STDERR "$id\t-\tfiltered\n";
		print LCA "$id\t-\tfiltered\n";
	    }
	    elsif ( $taxId && taxContains($taxId, $taxIds{$id}) )
	    {
		print LCA "$id\t$taxId\t" . getTaxName($taxId) . "\n";
	    }
	    else
	    {
		my $kingdom = $taxIds{$id};
		
		while ( $kingdom > 1 && getTaxRank($kingdom) !~ /kingdom/ )
		{
		    $kingdom = getTaxParent($kingdom);
		}
		
		print LCA "$id\t$taxIds{$id}\t" . getTaxName($taxIds{$id}) . "\t(" . getTaxName($kingdom) . ")\t$minLen{$id}\n";
	    }
	}
	
	close LCA;
    }
}

sub getSequenceLengths
{
    my ($fileName) = @_;
    
    open FILE1, "<$fileName" or die "Could not open $fileName";
    
    my $tag;
    my $length;
    my %lengths;
    
    while (my $line = <FILE1>)
    {
	if ( $line =~ /^#/ )
	{
	    next;
	}
	
	if ( $line =~ /> *([^ ]*).*\n/)
	{
	    if ( defined $tag )
	    {
		$lengths{$tag} = $length;
	    }
	    
	    $tag = $1;
	    $length = 0;
	}
	elsif ( defined $tag )
	{
	    chomp $line;
	    
	    $length += length $line;
	}
    }
    
    $lengths{$tag} = $length; # final length
    
    close FILE1;
    
    return %lengths;
}
