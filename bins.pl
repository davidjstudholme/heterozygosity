#!/usr/bin/perl -w

use strict;
use warnings;
use Statistics::Descriptive;
#use diagnostics;

my $usage = "Usage: $0 <samtools pileup file>";

my $file = shift or die "$usage\n";

### Read the data into memory
my %scaffold2pos2freq2depth;
open(FILE,"<$file") or die "Failed to open file '$file'\n$!\n";
warn "parsing file '$file'\n";
while (<FILE>) {
    chomp;
    my @fields = split /\t/;
    #warn "@fields\n";
    my ($scaffold, $pos, $ref_base, $reads, $alignment ) = @fields; 
    
    my ($proportion1, $proportion2) = (0, 0);
    if ($reads) {
	($proportion1, $proportion2) = get_proportions(\@fields);
    }
    $scaffold2pos2freq2depth{$scaffold}{$pos} = [$proportion1, $reads];
}
warn "Finished parsing file '$file'\n";


my $window_size = 5000;
my $step_size = 1000;
my $x = 0.10; ### wriggle room
my $min_depth = 5;

print "\"scaffold\"\t";
print "\"midpoint\"\t";
print "\"start\"\t";
print "\"end\"\t";
print "\"mean_depth\"\t";
print '"bin50"';
print "\t";
print '"bin66"';
print "\t";
print '"bin75"';
print "\t";
print '"bin100"';
print "\n";

foreach my $scaffold (keys %scaffold2pos2freq2depth) {
    my @positions =  sort {$a<=>$b} keys %{ $scaffold2pos2freq2depth{$scaffold} };
    my $length = $positions[-1];
    warn "length of $scaffold: $length\n";
    
    for (my $start_pos = 1; $start_pos + $window_size < $length; $start_pos += $step_size) {
	
	### Calculate stats for this window
	my $previous_position = 0;
	my $previous_scaffold = '';
	my @depths;
	my @proportions;
	my $bin50 = 0;
	my $bin66 = 0;
	my $bin75 = 0;
	my $bin100 = 0;
	foreach my $pos ( $start_pos .. ($start_pos+$window_size) ) {
	    $scaffold2pos2freq2depth{$scaffold}{$pos} = [0,0] unless defined $scaffold2pos2freq2depth{$scaffold}{$pos};
	    my ($proportion1, $reads) = @{ $scaffold2pos2freq2depth{$scaffold}{$pos} };
	    push @depths, $reads;
	    
	    if ($reads > $min_depth) {
		
		$bin50++ if $proportion1 > 0.50 - $x and $proportion1 < 0.50 + $x;
		$bin66++ if $proportion1 > 0.67 - $x and $proportion1 < 0.67 + $x;
		$bin75++ if $proportion1 > 0.75 - $x and $proportion1 < 0.75 + $x;
		
		### We want to calculate the mean and median frequencies of the major allele for non-homozygous sites
		if ($proportion1 > 1.0 - (1*$x)) {
		    $bin100++;
		} else {
		    push @proportions, $proportion1;
		    #warn "$proportion1:\t$scaffold $pos $ref_base $reads $alignment \n";
		}
	    }
	}    
	my $mid_pos = 0.5*($start_pos + ($start_pos + $window_size - 1));
	my $end_pos = $start_pos + $window_size - 1;
	my $depth_stat = Statistics::Descriptive::Full->new();
	$depth_stat->add_data(@depths); 
	my $depth_mean = $depth_stat->mean();
	my $depth_median = $depth_stat->median();
	my $depth_sum = $depth_stat->sum();
	
	print "$scaffold\t";
	print "$mid_pos\t";
	print "$start_pos\t";
	print "$end_pos\t";
	print "$depth_mean\t";
	print $bin50;
	print "\t";
	print $bin66;
	print "\t";
	print $bin75;
	print "\t";
	print $bin100;
	print "\n";
	
    }
    ### Free up some RAM
    delete $scaffold2pos2freq2depth{$scaffold}
}


exit;


sub get_proportions {
    my $fields = shift or die;
    my ($chromosome, $pos, $ref_base, $reads, $alignment ) = @{$fields}; 
 
    
    #warn "Reference base = $ref_base" unless $ref_base =~ m/^[ACGTN]$/i;
    $ref_base = 'N' unless $ref_base =~ m/^[ACGTN]$/i;


    ### Remove special symbols
    $alignment =~ s/\^.//gi; 
    $alignment =~ s/\$//gi; 
    #$alignment =~ s/\*//gi 
    
    ### Count the dots and commas
    my @agree_hits = ($alignment =~ m/[\,\.]/g );
    #warn "Done counting dots and commas\n";
    
    ### Count the deletions
    my @deletion;
    while  ($alignment =~ m/\-(\d+)/) {
	my $n = $1;
	if ( $alignment =~ m/(\-$n[\w]{$n})/) {
	    my $deletion = $1;
	    #warn "\n$_\n$deletion\n";
	    push @deletion, $deletion;
	    $alignment =~ s/$deletion//;
	}
    }
    #warn "Done counting deletions\n";
    
    ### Count the insertions
    my @insertion;
    while  ($alignment =~ m/\+(\d+)/) {
	my $n = $1;
	if ( $alignment =~ m/(\+$n[\w]{$n})/) {
	    my $insertion = "\\$1";
	    #warn "\n$_\n$insertion\n";
	    push @insertion, $insertion;
	    $alignment =~ s/$insertion//;
	}
    }
    #warn "Done counting insertions\n";
    
    ### Count substitutions
    my @a = ($alignment =~ m/a/gi );
    my @c = ($alignment =~ m/c/gi );	
    my @g = ($alignment =~ m/g/gi );
    my @t = ($alignment =~ m/t/gi );
    my @n = ($alignment =~ m/n/gi );
    my @star = ($alignment =~ m/\*/gi );
    #warn "Done counting substitutions\n";
    
    ### Do a sanity check that everything adds up!
    my $depth = @agree_hits + @a + @c + @g + @t + @n + @star;
    my $length = length($alignment);
    
    die unless $depth == $reads;
    
    ### Identify the change. Which is the most abundant base?
    my %base2counts;
    $base2counts{A} = @a;
    $base2counts{C} = @c;
    $base2counts{G} = @g;
    $base2counts{T} = @t;
    $base2counts{'*'} = @star;
    $base2counts{'ins'} = @insertion;
    $base2counts{'del'} = @deletion;
    $base2counts{'ref'} = @agree_hits;
    my %counts2bases;
    foreach my $base(keys %base2counts) {
	my $count = $base2counts{$base};
	$counts2bases{$count} = $base;
    }
    my @sorted_counts = sort {$b<=>$a} keys %counts2bases;
    my $highest_count = shift @sorted_counts;
    my $highest_base = $counts2bases{$highest_count};
    my $second_highest_count = shift @sorted_counts;
    $second_highest_count = 0 unless defined $second_highest_count;
    my $second_highest_base = $counts2bases{$second_highest_count};
    $second_highest_base = '' unless defined $second_highest_base;
    
    ### What is the proportion of the most abundant base?
    my $proportion1 = $highest_count / $depth;
    my $proportion2 = $second_highest_count / $depth;
    
    return ($proportion1, $proportion2); 
}

