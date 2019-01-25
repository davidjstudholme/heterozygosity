#!/usr/bin/perl -w

use strict;
use warnings;

my $file = shift or die;

die unless $file =~ m/.major-allele-freq.csv$/;

open (FILE, "<$file") or die $!;

my $heterozygous_count = 0;
my $homozygous_count = 0;
my $ambiguous_count = 0;

while(<FILE>) {
    chomp;
    if (m/(\S+)\t(\d+)\t([\.\d]+)$/) {
	#warn "'$1' '$2' '$3'\n";
	
	my $allele_freq = $3;
	if($allele_freq >= 0.98) {
	    $homozygous_count++;
	} elsif ($allele_freq >= 0.48 and $allele_freq <= 0.52) {
	    $heterozygous_count++;
	} else {
	    $ambiguous_count++
	}
    }
}

print "$file\t$heterozygous_count\t$homozygous_count\t$ambiguous_count\n";
warn "$file\t$heterozygous_count\t$homozygous_count\t$ambiguous_count\n";
