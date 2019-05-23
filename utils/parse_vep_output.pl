#!/usr/bin/env perl
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra

use warnings;
use strict;

my %SNPseverity = (
	"3_prime_UTR_variant" => 2,
	"5_prime_UTR_variant" => 2,
	"coding_sequence_variant" => 3,
	"downstream_gene_variant" => 1,
	"intergenic_variant" => 1,
	"intron_variant" => 2,
	"mature_miRNA_variant" => 2,
	"missense_variant" => 5,
	"non_coding_transcript_exon_variant" => 2,
	"non_coding_transcript_variant" => 2,
	"splice_acceptor_variant" => 5,
	"splice_donor_variant" => 5,
	"splice_region_variant" => 4,
	"start_lost" => 5,
	"stop_gained" => 5,
	"stop_lost" => 5,
	"stop_retained_variant" => 3,
	"synonymous_variant" => 3,
	"upstream_gene_variant" => 2
);

my $in_file = $ARGV[0];

open(IN_FILE, $in_file) or die "Can't open ".$in_file."\n";

while(<IN_FILE>){
	next if $_ =~ /^#/;
	
	chomp $_;
	my @cols = split("\t", $_);
	
	my $severity = (exists($SNPseverity{$cols[6]})) ? $SNPseverity{$cols[6]} : 0;
	if ($cols[6] =~ /,/){
		foreach my $c (split(",", $cols[6])){
			if (!exists($SNPseverity{$c})){ print "ERROR - no severity code for term - ".$c."\n"; }
			else{ $severity = ($SNPseverity{$c} > $severity) ? $SNPseverity{$c} : $severity; }
		}
	}
	
	my $symbol = "\\N";
	my $biotype = "\\N";
	my $impact = "\\N";
	my $sift_val = "\\N";
	my $sift_text = "\\N";
	my @extras = split(";", pop @cols);
	
	foreach my $ann (@extras){
		my ($k, $v) = split("=",$ann);
		if ($k eq "SYMBOL") { $symbol = $v; }
		if ($k eq "BIOTYPE"){ $biotype = $v; }
		if ($k eq "IMPACT") { $impact = $v; }
		if ($k eq "SIFT")   { ($sift_text, $sift_val) = $v =~ /(\w+)\((.*)\)/; }
	}
	
	print join("\t", @cols)."\t";
	print join("\t", $symbol, $severity, $biotype, $impact, $sift_text, $sift_val)."\t";
	print join("\;", @extras)."\n";
	
}

## Close infile
close(IN_FILE);