#!/usr/bin/perl
#
# Jonathan Klassen
# v1.1 April 29, 2020
#
# Derived from Clustcompare script blastp_to_BBH_list.pl v2.0
# Script to create a BLAST homology network table that can be opened in Cytoscape from an all-vs-all BLAST
# Assumes the following command to make the BLAST output formal:
# blastp -query [file] -db [file] -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" > AMPs_vs_self.blastp

# v1.1 fixes a bug that produced an incomplete node output file

use strict;
use warnings;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "BLAST_parser_for_cytoscape.pl

DEPENDANCIES: none

USAGE:
-a	minimum percent identify to declare a match						DEFAULT: 0			e.g., perl BLAST_parser_for_cytoscape.pl -i infile -a 70
-b	minimum percent sequence overlap to declare a match					DEFAULT: 0			e.g., perl BLAST_parser_for_cytoscape.pl -i infile -b 70
-c	minimum score to declare a match							DEFAULT: 0			e.g., perl BLAST_parser_for_cytoscape.pl -i infile -c 50
-d	maximum evalue to declare a match							DEFAULT: 0			e.g., perl BLAST_parser_for_cytoscape.pl -i infile -d 1e-5
-h	displays this usage statement (also using --help)
-i	input list of blastp output files (modified outfmt 7) produced by mult_blastp.pl	REQUIRED			e.g., perl BLAST_parser_for_cytoscape.pl -i infile
-o	output tab-delimited table of nodes in the BLAST homology network			DEFAULT: BLAST_nodes.tsv	e.g., perl BLAST_parser_for_cytoscape.pl -i infile -o Cytoscape_nodes.tsv
-p	output tab-delimited table of edges in the BLAST homology network			DEFAULT: BLAST_edges.tsv	e.g., perl BLAST_parser_for_cytoscape.pl -i infile -p Cytoscape_edges.tsv
-q	run quietly, i.e., no STDOUT (Y or N)							DEFAULT: N			e.g., perl BLAST_parser_for_cytoscape.pl -i infile -q Y

OUTPUT FILES:
	a table listing all nodes and the edges that connect them based on BLAST homology, compatible with Cytoscape
";

# input arguements

my %options = ();
GetOptions( 
	"a=i"  => \$options{id},
	"b=i"  => \$options{overlap},
	"c=i"  => \$options{score},
	"d=s"  => \$options{evalue},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"o=s"  => \$options{nodes_outfile},
	"p=s"  => \$options{edges_outfile},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input file is not specified:\n, $usage" if (!$options{infile});
my $infile = $options{infile};

# defaut arguments

unless ($options{id}){      $options{id} = 0};
unless ($options{overlap}){ $options{overlap} = 0};
unless ($options{score}){ $options{score} = 0};
unless ($options{evalue}){ $options{evalue} = 5};
unless ($options{nodes_outfile}){ $options{nodes_outfile} = "BLAST_nodes.tsv"};
unless ($options{edges_outfile}){ $options{edges_outfile} = "BLAST_edges.tsv"};
unless ($options{quiet}){   $options{quiet} = "N"};


# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
BLAST_parser_for_cytoscape.pl	Jonathan Klassen	v1.1	Apr 29, 2020

parameters used:
	percent sequence id threshold = $options{id}
	minimum score = $options{score}
	maximum evalue = $options{evalue}
	percent sequence overlap threshold = $options{overlap}
	input file = $options{infile}
	output nodes file = $options{nodes_outfile}
	output edges file = $options{edges_outfile}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads list of input files
#############################################################################

my %best_hits;
my %found;
open (INFILE, $infile) or die "Cannot open $infile";
print "Parsing $infile\n";
while (<INFILE>){
	next if (/^#/);
	my @line = split /\t/;
	next if (scalar @line != 14); 
	$found{$line[0]} = "y";
#	next if ($line[0] eq $line[1]); # no self hits
	next if ($line[2] < $options{id}); # no low %ID hits
	next if ($line[11] < $options{score}); # no low score hits
	next if ($line[10] > $options{evalue}); # no high evalue hits
	my $max_length = 0; # no low % overlap hits
	if ($line[12] >= $line[13]){
		$max_length = $line[12];
	}
	else {
		$max_length = $line[13];
	}
	next if ($line[3] / $max_length * 100 < $options{overlap});
	$best_hits{$line[0]}{$line[1]}{id} = $line[2]; # keep %ID if all thresholds passed
	$best_hits{$line[0]}{$line[1]}{score} = $line[11]; # keep score if all thresholds passed
	$best_hits{$line[0]}{$line[1]}{evalue} = $line[10]; # keep evalue if all thresholds passed
	$best_hits{$line[0]}{$line[1]}{aligned} = $line[3] / $max_length * 100; # keep %aligned if all thresholds passed
}
close INFILE;


############################################################################
# Creates output table of all BBHs
############################################################################

print "Creating output files\n";

# find all possible nodes
open (OUTFILE, ">$options{nodes_outfile}") or die "Cannot open $options{nodes_outfile}";
print OUTFILE "Sequence_name\n";

foreach my $hit (sort keys %found){
	print OUTFILE $hit, "\n";
}
close OUTFILE;

# find all BBHs
open (OUTFILE, ">$options{edges_outfile}") or die "Cannot open $options{edges_outfile}";
print OUTFILE "Sequence_name_1\tSequence_name_2\tavg_\%id\tavg_score\tavg_evalue\tavg_\%aligned\n";

my %found_bbh;
foreach my $hit1 (sort keys %best_hits){
	foreach my $hit2 (sort keys %{$best_hits{$hit1}}){
		next unless ($best_hits{$hit1}{$hit2} and $best_hits{$hit2}{$hit1}); # require bidirectional hits
		next if ($best_hits{$hit1}{$hit2} eq $best_hits{$hit2}{$hit1}); # no self hits
		next if ($found_bbh{$hit1}{$hit2} or $found_bbh{$hit2}{$hit1});
		my $avg_id = ($best_hits{$hit1}{$hit2}{id} + $best_hits{$hit2}{$hit1}{id}) / 2;
		my $avg_score = ($best_hits{$hit1}{$hit2}{score} + $best_hits{$hit2}{$hit1}{score}) / 2;
		my $avg_evalue = ($best_hits{$hit1}{$hit2}{evalue} + $best_hits{$hit2}{$hit1}{evalue}) / 2;
		my $avg_aligned = ($best_hits{$hit1}{$hit2}{aligned} + $best_hits{$hit2}{$hit1}{aligned}) / 2;
		print OUTFILE "$hit1\t$hit2\t$avg_id\t$avg_score\t$avg_evalue\t$avg_aligned\n";
		$found_bbh{$hit1}{$hit2} = $found_bbh{$hit2}{$hit1} = "y";
	}
}
close OUTFILE;

