# BLAST_parser_for_cytoscape
Script to parse an all-vs-all BLAST into a homology network that can be visualized using Cytoscape

v1.1 April 29, 2020

Derived from Clustcompare script blastp_to_BBH_list.pl v2.0
Script to create a BLAST homology network table that can be opened in Cytoscape from an all-vs-all BLAST
Assumes the following command to make the BLAST output formal:
blastp -query [file] -db [file] -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" > AMPs_vs_self.blastp

v1.1 fixes a bug that produced an incomplete node output file

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
