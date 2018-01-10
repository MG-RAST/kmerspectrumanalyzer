#!/usr/bin/perl -w

use Getopt::Long;
use strict;
my $usage = qq( Usage:
    multijoin.pl  -i <index column> -s <sum column> <input files> 
    This script creates a table of numbers from column s.
    Rows are indexed by column i and columns are given by input files.
);

my $verbose = "";
my $rank =0;
my $index =0;
my $sumcol =0;
my $sep ="\t";
my $binomial =0;
my $sortme =0;
my %keyindex =();
my @hash = ();
if ( (@ARGV > 0) && ($ARGV[0] =~ /-h/) ) { print $usage;  }
if ( ! GetOptions (
                   "verbose!"       => \$verbose,
                   "rank!"       => \$rank,
                   "sep:s"       => \$sep,
                   "i:i"            => \$index,
                   "s:i"            => \$sumcol,
                   "binomial!"            => \$binomial,
                   "sort!"            => \$sortme,
                  )
   ) { print $usage; }
$index--;
$sumcol--;

if($#ARGV+1 ==0) {die "put filenames to parse on command line";}
my @filenames = @ARGV;
print "#Index\t";
foreach my $file (@filenames)
	{
	if( !(-e $file) ) {die "Can't find file $file \n$usage";}
	print "$file\t";
	}
print "\n";
my $i= 0;
foreach my $file (@filenames) 
	{
	print STDERR "Processing file $file\n";
	open my $FILE, "<$file";
	while(<$FILE>)
		{
		chomp;
                if ( $_ =~ s/$sep/$sep/g == 0 )   # if at first you don't find the separator
			{$sep = " ";}             # change the separator
		my @a = split(/$sep/);

#                if ( $#a + 1 < 1 ) { die "Can't find fields in file $file!"; } 
	my $key;
	if($binomial) {my @w = split(/ /,$a[$index]); $key = "$w[0] $w[1]";}
	else { $key =$a[$index]; }
	
		if(!exists($keyindex{$key})) {$keyindex{$key} =1;  }

		if(exists($hash[$i]{$key} ) ) { $hash[$i]{$key}  += $a[$sumcol] ; } 
		else { $hash[$i]{$key}   = $a[$sumcol] ;    } 
		}
	$i++;
	}
if($sortme)
{
print STDERR "Sorting and outputting\n";
foreach my $key (sort  {$a <=> $b} keys %keyindex)
	{
	print "$key\t";
	for ($i=0; $i< $#filenames+1; $i++)
		{
		if(!defined($hash[$i]{$key})) {$hash[$i]{$key} = "0";}
		print "$hash[$i]{$key}\t";	
		}
	print "\n";
	}
}else
{
print STDERR "Outputting\n";
foreach my $key ( keys %keyindex )
	{
	print "$key\t";
	for ($i=0; $i< $#filenames+1; $i++)
		{
		if(!defined($hash[$i]{$key})) {$hash[$i]{$key} = "0";}
		print "$hash[$i]{$key}\t";	
		}
	print "\n";
	}
}
