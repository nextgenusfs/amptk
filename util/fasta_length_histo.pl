#!/usr/bin/env perl

#script takes a multi-fasta file and prints out stats and length histogram.
#majority of this script from here: http://seqanswers.com/forums/showthread.php?t=15856

use warnings;
use strict;
use Statistics::Descriptive;
use Carp;
use English qw(-no_match_vars);

my $stat = Statistics::Descriptive::Full->new();
my (%distrib);
if (not defined $ARGV[0]) {
    die "Need fasta file input\nusage: fasta_length_histo.pl file.fasta\n";
}
my $fastaFile = $ARGV[0];
my @shortname = split(/\.fa/,$fastaFile);
my $outputFile = $shortname[0];
$outputFile .= "_histogram\.txt";
open (FASTA, "<$fastaFile");
$/ = ">";

my $junkFirstOne = <FASTA>;

while (<FASTA>) {
	chomp;
	my ($def,@seqlines) = split /\n/, $_;
	my $seq = join '', @seqlines;
	$stat->add_data(length($seq));
}
my $start = $stat->min();
my $end = $stat->max();
my @bins = ($start..$end);
%distrib = $stat->frequency_distribution(\@bins);

print "Total reads:\t" . $stat->count() . "\n";
print "Total nt:\t" . $stat->sum() . "\n";
print "Mean length:\t" . $stat->mean() . "\n";
print "Median length:\t" . $stat->median() . "\n";
print "Mode length:\t" . $stat->mode() . "\n";
print "Max length:\t" . $stat->max() . "\n";
print "Min length:\t" . $stat->min() . "\n";
print "Histogram data saved to file:  $outputFile\n";

#check if output file exists
if (-f $outputFile) {
    unlink $outputFile
        or croak "Cannot delete $outputFile: $!";
}
my $OUTFILE;
open $OUTFILE, '>>', $outputFile
    or croak "Cannot open $outputFile: $OS_ERROR";
print { $OUTFILE } "Length\t# Seqs\n";
foreach (sort {$a <=> $b} keys %distrib) {
	print { $OUTFILE } "$_\t$distrib{$_}\n";
}