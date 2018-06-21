#!/usr/bin/perl
# A program to a) parse the vcf header b) select variants progressively

use strict;
use Data::Dumper;
#use Acme::Tools qw(minus);

# make sure the user passed the required arguments
if (scalar @ARGV != 2 ) {
   print STDERR "Usage: progressiveSelection.pl <vcf file> <samples in region-segment file>\n";
 	  exit(1);
}
chomp(@ARGV);
my ($vcffile, $sfile) = @ARGV;

open(my $isfh, "<", $sfile) || die "ERROR: failed to read regions file: $!";

#my %reg_sam;
#my @target_samples;
my %regionIDX; # {segment} -> [] = sample ID blocks
my %regions; # {chr} -> {start} = end

while(my $line = <$isfh>){
	#print "Sample information loaded\n";
	chomp $line;
	$line =~ s/\r//g;
	my @data  = split(/\t/, $line);
	#$reg_sam{$data[0]} = $data[1]; 
	push(@{$regionIDX{$data[0]}}, $data[1]); 
	my ($chr, $start, $end) = $data[0] =~ m/(.+):(\d+)-(\d+)/;
	#@target_samples = split(",", $data[1])
	#print Dumper $regionIDX{$data[0]};
	$regions{$chr}->{$start} = $end;
}

close $isfh;

# Now open vcf and go to the region and sample to check if the variant is homozygous in target samples

open(my $ifh, "gunzip -c $vcffile |") || die "ERROR: failed to read vcf file: $!";
open (my $OUTFILE, "> filtered_run2/filtered_run2.vcf");

my @samples;
my %sampleIDX; # {sample} = index of gtype

while (my $line = <$ifh>) {
	chomp $line;
	my @breaks = split (/\t/, $line);
	if ($breaks[0] =~ /^#CHROM/) { 
	      	@samples = @breaks[9..scalar(@breaks)];
		for(my $x = 9; $x < scalar(@breaks); $x++){
			if(length($breaks[$x]) > 10){
				($breaks[$x]) = $breaks[$x] =~ /^(.{10})/;
				print "Shortened $breaks[$x]\n";
			}
			$sampleIDX{$breaks[$x]} = $x;
		}
		print $OUTFILE "$line\n";
		last;
      	}
	print $OUTFILE "$line\n";
}
print "Sample information loaded and indexed\n";

my $mean_dp=4337;
my $std_dp=1555;
my $min_dp= (3*$std_dp) - $mean_dp;
my $max_dp= $mean_dp + (3*$std_dp);
my @depth;

#$mean_dp = average(@depth);
#$std_dp = stdev(@depth);

print "Calculated mean and standard deviation of DP\n";

# Now select the variants

my $filter1_count = 0; # Selected variants
my $filter2_count = 0; # Discarded variants
my $filter3_count = 0; # Filter for singletons within regions
my $filter4_count = 0; # Filter for hetsOnly
my $filter5_count = 0; # AN filter 
my $filter6_count = 0; # MQ0F filter
my $filter7_count = 0; # MQSB filter
my $filter8_count = 0; # DP filter


while(my $line = <$ifh>){
	chomp $line;
	my @segs = split(/\t/, $line);
	#next if($segs[0] ne $chr || ($segs[1] < $start || $segs[1] > $end));		# skip if not in region of interest
	my $hapUCSC = returnRegion(\%regions, $segs[0], $segs[1]);
		# filter 1: variant should be homozygous in target animals and non heterozygous in non-target samples
	if(exists($regionIDX{$hapUCSC})){
		my @groups = @{$regionIDX{$hapUCSC}}; 

			
		my ($ac) = $segs[7] =~ /AC=(\d{1,4})\;/;
		if($ac == 1){
			print "$segs[0]\t$segs[1]\tsingleton\n";
			$filter3_count++;
			next;
		}
		
		my $valid = isHomOnlyOne(\@segs, \%sampleIDX, \@groups);
		my $hetall = isHetOnly(\@segs, \%sampleIDX, \@groups);
			my ($an) = $segs[7] =~ /AN=(\d{1,4})\;/;
			my ($dp) = $segs[7] =~ /DP=(\d{1,6})\;/;
			push @depth, $dp;
			my ($mq0f) = $segs[7] =~ /MQ0F=(\d{1,4}\.?\d*)\;/;
			my ($mqb) = $segs[7] =~ /MQB=(\d+\.?\d*)\;/;
			my ($mqsb) = $segs[7] =~ /MQSB=(\d+\.?\d*)\;/;
			
						
		if($valid && ! $hetall && $an > 337 && $mq0f <= 0.1 && $mqsb >= 0.95 && $dp > $min_dp && $dp < $max_dp){
		#if($valid && ! $hetall && $dp > ($mean_dp - (3*$std_dp))  && $dp < ($mean_dp + (3*$std_dp))){	
			print "$segs[0]\t$segs[1]\tkept\n";
			$filter1_count++;
			print {$OUTFILE} join("\t", @segs) . "\n";
		}else{
			print "$segs[0]\t$segs[1]\tnotkept\n";
			$filter2_count++ if (!$valid);
			$filter4_count++ if ($hetall);
			$filter5_count++ if ($an < 337);							#AN > 337    (2% reduction from total max AN due to CNV)
			$filter6_count++ if ($mq0f >= 0.1);							#MQ0F filter away >= 0.1
			$filter7_count++ if ($mqsb <= 0.95);							#MQSB filter away <= 0.95
			$filter8_count++ if ($dp < $min_dp && $dp > $max_dp);					#DP filtered on +/- 3 stdevs   (remove repeats, condensed duplications, faulty calls, etc)
		}	
	}else{
		$filter2_count++;
	}	
}

print "Kept $filter1_count out of " . ($filter1_count + $filter2_count) . " variant sites\n";
print "Removed $filter3_count singletons\n";
print "Removed $filter4_count sites that were heterozygous in all other regions\n";
print "Removed $filter5_count sites with AN < 337\n";
print "Removed $filter6_count sites with MQ0F >= 0.1\n";
print "Removed $filter7_count sites with MQSB <= 0.95\n";
print "Removed $filter8_count sites with $min_dp > DP > $max_dp\n";



close $ifh;
close $OUTFILE;


exit;
sub average{
        my @data = @_;
        my $total = 0;
        foreach my $v (@data) {
                $total += $v;
        }
        my $average = $total / scalar(@data);
        return ($average);
}

sub stdev{
        my @data = @_;
        if(scalar(@data) == 1){
                return 0;
        }
        my $average = average(@data);
        my $sqtotal = 0;
        foreach my $v (@data) {
                $sqtotal += ($average-$v) ** 2;
        }
        my $std = ($sqtotal / scalar(@data)-1) ** 0.5;
        return $std;
}

sub isHetOnly{
	my ($segs, $sampidx, $groups) = @_;
	
	my %exclude; # container for indicies that we want to ignore
	foreach my $g (@{$groups}){
		foreach my $s (split(/,/, $g)){
			$exclude{$sampidx->{$s}} = 1;
		}
	}

	my $hetCount = 0;
	for(my $x = 9; $x < scalar(@{$segs}); $x++){
		if(! exists($exclude{$x})){
			my @gtarray = split(/:/, $segs->[$x]);
			if($gtarray[0] =~ /[123]/){
				$hetCount++;
			}
		}
	}
	return ($hetCount >= scalar(@{$segs}) - 9 - scalar(keys(%exclude)))? 1 : 0;
}

sub isHomOnlyOne{
	my ($segs, $sampidx, $groups) = @_;
	
	my @gtsum; # container for boolean homozygous group
	foreach my $g (@{$groups}){
		my @gtarray; 
		my $homozygous = 1; 
		foreach my $s (split(/,/, $g)){
			my $idx = $sampidx->{$s};
			my @gtsegs = split(/:/, $segs->[$idx]);
			push(@gtarray, $gtsegs[0]);
		}
		# Implement test for single genotype groups
		if(scalar(@gtarray) == 1){
			my @sgtsegs = split(/[\/|]/, $gtarray[0]);
			# test if this is not homozygous
			$homozygous = ($sgtsegs[0] eq $sgtsegs[1])? 1 : 0;
			if($homozygous && $sgtsegs[0] eq "0"){
				# if homozygous reference, check to see if this is a singleton site
				my ($ac) = $segs->[7] =~ /AC=(\d{1,4})\;/;
				if($ac == 1){
					return 0;
				}
			}
		}
		for(my $x = 1; $x < scalar(@gtarray); $x++){
			my @sgtsegs = split(/[\/|]/, $gtarray[$x]);
			if($sgtsegs[0] ne $sgtsegs[1]){
				$homozygous = 0;
				last;
			}
			if($gtarray[$x] ne $gtarray[$x - 1] || $homozygous == 0){
				$homozygous = 0;
			}
		}
		push(@gtsum, ($homozygous)? $gtarray[-1] : "0");
		# If the variant site is homozygous, store the homozygous variant site in the array
	}
	my %count; # If count is greater than 1 or less than zero, it's not a good region
	for(my $x = 0; $x < scalar(@gtsum); $x++){
		if($gtsum[$x] ne "0"){
			$count{$gtsum[$x]} += 1;
		}
	}
	# Take the max value of count for the variant site and keep it only if it is "1"
	my $final = maxKeys(\%count);
	return ($final == 1)? 1 : 0; 
}

sub maxKeys{
	my ($href) = @_;
	my $max = -1;
	foreach my $k (keys(%{$href})){
		if($href->{$k} > $max){
			$max = $href->{$k};
		}
	}
	return $max;
}

sub returnRegion{
	my ($href, $chr, $pos) = @_;
	if(!exists($href->{$chr})){
		return "-1";
	}
	foreach my $starts (sort{$b cmp $a}keys(%{$href->{$chr}})){
		if($starts < $pos && $href->{$chr}->{$starts} >= $pos){
			return "$chr:$starts-" . $href->{$chr}->{$starts};
		}
	}
	return "-1";
}
