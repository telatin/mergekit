#!/usr/bin/env perl
#ABSTRACT - Check quality of merged paired ends

use 5.012;
use warnings;
use Getopt::Long;
use FASTX::Reader;
use FindBin qw($Bin);
use lib "$Bin/lib";
use Local::Align qw(align);
use Data::Dumper;
use Term::ANSIColor;
use File::Basename;
use utf8;

my $opt_visualqual = 1;
binmode STDOUT, ":utf8";
my $usage=<<END;
   USAGE:
     16S.pl file.fa

END

my ($opt_debug, $opt_max_seq);

my $_opt = GetOptions(
	'd|debug'        => \$opt_debug,
	'm|max=i'        => \$opt_max_seq,
);

if (not defined $ARGV[0]) {
	say $usage;
	exit;
}

my %regions = (
	V1   => [68, 99],
	V2   => [136, 242],
	V3   => [338, 533],
	V4   => [576, 682],
	V5   => [821, 879],
	V6   => [970, 1046],
	V7   => [1117, 1294],
	V8   => [1435, 1465],
);
my $reference_16S = load_ref();

my $query = FASTX::Reader->new( {filename => "$ARGV[0]" });

while (my $seq = $query->getRead() )  {
	last if (defined $opt_max_seq and $query->{counter} > $opt_max_seq);
	my ($top, $middle, $bottom, $score) = align($seq->{seq}, $reference_16S);
	my $start = 0;
	my $len = 0;
	my $end = 0;
	my $ref = '';
	if ($top=~/^([-]*)([A-Z].*?[A-Z])([-]*)$/) {
		$start = length($1);
		my $query_match = $2;
		my $ref_match   = substr($bottom, $start, length($query_match));
		$ref_match=~s/[-]//g;
		$len = length($ref_match);
		$ref = substr($reference_16S, $start, $len);
		$end = $start + $len;
	}
	next if (not defined $start or not defined $len);
	my @r = ();
	foreach my $region (sort keys %regions) {
		my $s = $regions{$region}[0];
		my $e = $regions{$region}[1];

		if ( ($s > $start) and ($e < $end) ) {
			 #included
			push(@r, "$region");
		} elsif (( $e > $start ) and ( $s < $start) ) {
			push(@r, "$region");
		} elsif ( ( $s < $end ) and ( $e > $end ) ) {
			push(@r, "$region");
		}
	}
	say $seq->{name}, "\t", "regions=",join(",", @r), "\tregion=$start-$end;score=$score";
	say join("\n", $top, $middle, $bottom) if ($opt_debug);

}




sub load_ref   {
 return 'AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAGCTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAATGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA';
}

