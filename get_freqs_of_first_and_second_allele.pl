#!/usr/bin/perl -w

use strict;
use warnings;

my $usage = "Usage: $0 <samtools pileup file> <min depth>";

my $file = shift or die "$usage\n";
my $min_depth = shift or die "$usage\n";

open(FILE,"<$file") or die "Failed to open file '$file'\n$!\n";
warn "parsing file '$file'\n";
while (<FILE>) {
    chomp;
    my @fields = split /\t/;
    my ($chromosome, $pos, $ref_base, $reads, $alignment ) = @fields; 
    
    
    #warn "@fields\n";
    
    if ($ref_base =~ m/^[ACGTN]$/i and
	$reads >= $min_depth and
	$alignment =~ m/^[ACGT\^\$\,\.]+$/i ) {
	#warn "Alignment: '$alignment'\n";
	
	### Parse the alignment string ...
	
	### Remove special symbols for read ends and starts
	$alignment =~ s/\^.//gi; 
	$alignment =~ s/\$//gi; 
	
	### Count the dots and commas
	my @agree_hits = ($alignment =~ m/[\,\.]/g );
	#warn "Done counting dots and commas\n";
	
	### Count the deletions
	my @deletion;
	while  ($alignment =~ m/\-(\d+)/) {
	    my $n = $1;
	    if ( $alignment =~ m/(\-$n[ACGTNacgtn]{$n})/) {
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
	    if ( $alignment =~ m/(\+$n[ACGTNacgtn]{$n})/) {
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
	
	if ($depth == $reads) {
	    
	    ### Identify the change. Which is the most abundant base?
	    my %base2counts;
	    $base2counts{a} = @a;
	    $base2counts{c} = @c;
	    $base2counts{g} = @g;
	    $base2counts{t} = @t;
	    $base2counts{ lc($ref_base) } = @agree_hits;
	    #$base2counts{'*'} = @star;
	    #$base2counts{'ins'} = @insertion;
	    #$base2counts{'del'} = @deletion;
	    #$base2counts{'ref'} = @agree_hits;
	    my %counts2bases;
	    foreach my $base(keys %base2counts) {
		$base = lc($base);
		if ($base eq 'ref') {
		    $base = $ref_base;
		}
		
		my $count = $base2counts{$base};
		$counts2bases{$count}{$base}++;
		#warn "base '$base' => '$count'\n";
	    }
	    
	    ### Generate an ordered list of alleles, ordered by abundance
	    my @sorted_counts;
	    if (keys %counts2bases) {
		@sorted_counts = sort {$b<=>$a} keys %counts2bases;
	    } else {
		warn "No eligible bases found for position $chromosome: $pos\n";
	    }
	    my @ordered_bases;
	    foreach my $count (@sorted_counts) {
		foreach my $base (keys %{  $counts2bases{$count} }) {
		    push @ordered_bases, $base;
		}
	    }
	    
	    my $most_abundant_base = shift @ordered_bases;
	    my $secondmost_abundant_base = shift @ordered_bases;
	    
  	    ### What is the proportion of the first-most abundant base?
	    my $proportion_firstmost;
	    if (defined $most_abundant_base) {
		$proportion_firstmost = $base2counts{$most_abundant_base} / $depth;
		

		unless ($proportion_firstmost) {
		    warn "Most abundant base is $most_abundant_base: $proportion_firstmost\n";
		    warn "base2counts =\n";
		    foreach my $key (keys %base2counts) {
			warn "\t$key => $base2counts{$key}\n";
		    }
		    die "Alignment = '$alignment'\n";
		}
		
	    } else {
		warn "At $chromosome $pos I cannot determine the most abundant base !\n\t$alignment\n";
		$proportion_firstmost = 'NA';
	    }
	    


	    ### What is the proportion of the second-most abundant base?
	    my $proportion_secondmost;
	    if (defined $secondmost_abundant_base) {
		$proportion_secondmost = $base2counts{$secondmost_abundant_base} / $depth;
	    } else {
		$proportion_secondmost = 0;
		#No second base!
	    }
	    #warn "$proportion_secondmost\n";
	    
	  	    
	    ### Only print this if the first and second most abundant are A, C, G ot T
	    #warn "Two most abundant alleles are: $most_abundant_base and $secondmost_abundant_base\n";
	    if ($most_abundant_base =~ m/^[ACGT]$/i and  $secondmost_abundant_base =~ m/^[ACGT]$/i) {
		print "$chromosome\t$pos\t$proportion_firstmost\t$proportion_secondmost\t$depth\n";
	    }
	    
	} else {
	    warn "$_\nlength($alignment)=$length\n$depth != $reads\n";
	}
	
    }
    
}
