#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#########################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$type,$popName);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$popName,
				"o:s"=>\$fOut,
				"t:s"=>\$type,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $type and $popName);
$type="DH" if($type=~/BC\d+/i);     #马立祥
#
# read in loc file 
#
my @marker_info ;
my $number_of_individual;

open (LOC , $fIn ) or die $!;
while (<LOC>) {
	chomp;
	next if (/^;/ || /^\s*$/ || /=/) ;
	s/\ba\b/A/g;
	s/\bb\b/B/g;
	s/-/U/g;
	s/\bh\b/X/g;
	my ($marker , @genotype ) = split /\s+/ , $_ ;
	next if (@genotype < 5);
	my $genotype = join("\t",@genotype) ;
	push @marker_info , "$marker\t$genotype";
	$number_of_individual = scalar @genotype ;
}
close (LOC) ;

my $number_of_loci = scalar (@marker_info) ;


#
# output MSTmap loc 
#
open (OUT, ">", $fOut) or die $!;

my $head = <<"Headend" ;
population_type $type
population_name $popName
distance_function kosambi
cut_off_p_value 2.0
no_map_dist 15.0
no_map_size 0
missing_threshold 1.00
estimation_before_clustering no
detect_bad_data no
objective_function COUNT
number_of_loci $number_of_loci
number_of_individual $number_of_individual
Headend

print OUT $head ,"\n";
print OUT "locus_name\t";
map {print OUT "i$_\t" } (1..$number_of_individual);
print OUT "\n" ;
foreach (@marker_info) {
	print OUT $_ ,"\n";
}
close (OUT) ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage: 将loc文件转换为MSTmap所需输入格式，参考http://alumni.cs.ucr.edu/~yonghui/mstmap.html
  Options:
  -help		USAGE
  -i		input loc file, forced  
  -k		population name  
  -o		output file 
  -t		population type: RIL{2,3,4,5,6..}, DH, 
  -h		help

USAGE
	print $usage;
	exit;
}
