#!/usr/bin/perl -w
# June 16, 2011
use Getopt::Long;
use strict;
my $usage = qq( Usage:
    join2d.pl  -i <index column> -j <index column> -s <sum column> <input file> 
    This script takes a single file containing three columns--two index columns 
    and one target varaible-- and turns it into a m x n table.
    sum column = 0 causes counter behavior.
);

my $verbose = "";
my $rank =0;
my $index =1;
my $index2 =2;
my $sumcol =0;
my $binomial =0;
my $delim ="\t";
my %key1index =();
my %key2index =();
#my %hash = ();
my $hashref ; # = \%hash;
my $count = 0;
if ( (@ARGV > 0) && ($ARGV[0] =~ /-h/) ) { print $usage;  }
if ( ! GetOptions (
                   "verbose!"       => \$verbose,
                   "rank!"       => \$rank,
                   "j:i"            => \$index,
                   "i:i"            => \$index2,
                   "s:i"            => \$sumcol,
                   "d:s"            => \$delim,
                   "binomial!"            => \$binomial,
                  )
   ) { print $usage; }
if($sumcol == 0) {$count = 1; $sumcol=1;}

$index--;
$index2--;
$sumcol--;

if($#ARGV+1 ==0) {die "$usage";}
my @filenames = @ARGV;
my $i= 0;
my $target=1;
my $key1; my $key2;
foreach my $file (@filenames) 
	{
	print STDERR "Processing file $file\n";
	open my $FILE, "<$file";
	while(<$FILE>)
		{
		chomp;
		my @a = split(/$delim/);
#		print "$a[$index]\t$a[$sumcol]\n";
#	if($binomial) {my @w = split(/ /,$a[$index]); $key = "$w[0] $w[1]";}
	 $key1 =$a[$index]; 
	 $key2 =$a[$index2]; 
	if($count==0) {$target = $a[$sumcol];} 
	     else   {$target = 1; ;}
#		if(!exists($hashref->{$key1}->{$key2})) 
{$key1index{$key1} =1; $key2index{$key2} = 1;}

if(exists($hashref->{$key1}->{$key2} )) 
	{
	if($count ==0) {print STDERR "Duplicate ($key1, $key2) already equal to  $hashref->{$key1}->{$key2} \n";}
	         else {$hashref->{$key1}->{$key2} += $target;}
	} 
	else
	{
	$hashref->{$key1}->{$key2} = $target;
	}

		}
	$i++;
	}
# Print first (header) row:
print "0\t";
foreach my $key1 (sort {$a <=> $b} keys %key1index)
	{
	print "$key1\t";
	}
print "\n";

# Print data table:
foreach my $key2 (sort {$a <=> $b} keys %key2index) 
{
print "$key2\t";
foreach my $key1 (sort {$a <=> $b} keys %key1index)
	{
# print "$key1\t";
if(!(exists($hashref->{$key1}->{$key2}))) { $hashref->{$key1}->{$key2} = 0;}
 
print "$hashref->{$key1}->{$key2}\t";	

	}
print "\n";
	}
