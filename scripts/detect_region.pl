#!/usr/bin/env perl
#ABSTRACT - Detect the hypervariable regions of input sequences, in either FASTA or FASTQ format

use 5.012;
use warnings;
use Getopt::Long;
use FASTX::Reader;
use FindBin qw($Bin);
use lib "$Bin/lib";
use Local::Align qw(align);
use Data::Dumper;
use File::Basename;
use JSON::PP;

my $usage=<<END;
   USAGE:
     16S.pl [options] file.fa

   OPTIONS:
     -m, --max  INT       Read 'm' number of sequences and then stop

     -t, --thr  FLOAT     Percentage of coverage of target region
                          (0-100), if below the coverage percentage will
                          be reported

     -j, --json           Prints output in JSON format
     -p                   Pretty print JSON
     -f                   Add parameters to JSON
END

my ($opt_debug, $opt_max_seq, $opt_json, $opt_pretty, $opt_full);
my $opt_ths = 80;
my $opt_min_score = 100;
my %output;

my $_opt = GetOptions(
	'd|debug'        => \$opt_debug,
	't|threshold=f'  => \$opt_ths,
	'm|max=i'        => \$opt_max_seq,
	's|min-score=f'  => \$opt_min_score,
	'j|json'         => \$opt_json,
	'p|pp'           => \$opt_pretty,
	'f|full'         => \$opt_full,
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

my %regions_counter = ();

my $reference_16S = load_ref();

if ($opt_full) {
	$output{reference_regions} = \%regions;
	$output{params}{max_seqs} = $opt_max_seq;
}

my $query = FASTX::Reader->new( {filename => "$ARGV[0]" });

while (my $seq = $query->getRead() )  {
	last if (defined $opt_max_seq and $query->{counter} > $opt_max_seq);
	
	my ($top, $middle, $bottom, $score) = align($seq->{seq}, $reference_16S);
	my $start = 0;
	my $len = 0;
	my $end = 0;
	my $ref = '';
	my $data;
	$data->{align_score} = $score;
	next if ($score < $opt_min_score);
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
		my $l = $e - $s + 1;
		my $span = 0;
		if ( ($s > $start) and ($e < $end) ) {
			 #included
			$regions_counter{$region}++;
			push(@r, "$region");
			$span = 100.00;

		} elsif (( $e > $start ) and ( $s < $start) ) {
			$span =  100 * ( $e - $start ) / $l ; 
			if ($span > $opt_ths) {
				push(@r, "$region");
				$regions_counter{$region}++;
			} else {
				push(@r, "$region($span)");
			}

		} elsif ( ( $s < $end ) and ( $e > $end ) ) {
			$span =  100 * ( $end - $s ) / $l ; 
			if ($span > $opt_ths) {
				$regions_counter{$region}++;
				push(@r, "$region($span)");
			}
		}
		$data->{regions}->{$region} = sprintf("%.2f", $span) + 0 if ($span);

	}
	if (not $opt_json) {
		say $seq->{name}, "\t", "regions=",join(",", @r), "\tregion=$start-$end;score=$score";
		say join("\n", $top, $middle, $bottom) if ($opt_debug);
	} else {
		$data->{detected_regions} = join(',', @r);		
		$output{input_seqs}{ $seq->{name} } = $data;

	}

}

foreach my $r (keys %regions_counter) {
	$regions_counter{$r} /= ($query->{counter} - 1) if ($query->{counter} > 1);
}
$output{global_seqs}{hit_ratios} = \%regions_counter;
$output{global_seqs}{parsed_seqs} = $query->{counter} - 1;
if ($opt_json) {
  my $json = JSON::PP->new->ascii->allow_nonref;
  $json = JSON::PP->new->ascii->pretty->allow_nonref if ($opt_pretty);
  say $json->encode( sort \%output );
}



sub load_ref   {
 return 'AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAGCTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAATGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA';
}

