#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

##############.......... .....SMOOTH a statistical method for successful removal of genotyping errors from high-density genetic linkage data
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------


my ($fIn,$fOut,$threshold,$windows_size,$min_mis);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"D:s"=>\$threshold,
				"M:s"=>\$min_mis,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
#
# command line options 
#
$threshold ||= 0.95;
$windows_size ||= 15 ;
$min_mis ||= 0.85 ;

#
# read ordered loc 
#
open ( LOC , $fIn ) or die $!;
my %genotype;
my @marker;my @head ;
while (<LOC>) {
	chomp;
	next if ( /^\s*$/ );
	if (/=/) {
		push @head ,$_;
	}
	my ($marker,@genotype) = split /\s+/ , $_ ;
	next if (@genotype < 6) ;
	push @marker , $marker;
	for (my $i =0 ; $i < @genotype ;$i++) {
		$genotype{$i}{$marker} = $genotype[$i] ;
	}
}
close (LOC) ;
#
# score weights 
#
my %weight = (
	1  => 0.998,
	2  => 0.981,
	3  => 0.934,
	4  => 0.857,
	5  => 0.758,
	6  => 0.647,
	7  => 0.537,
	8  => 0.433,
	9  => 0.342,
	10 => 0.265,
	11 => 0.202,
	12 => 0.151,
	13 => 0.112,
	14 => 0.082,
	15 => 0.059,
);

my %genotype_new;
my $count ;
my @pos = sort {$a <=> $b} keys %genotype ;
foreach my $pos (sort {$a <=> $b}  keys %genotype) {
	my @genotype;my @genotype_new;
	for (my $i=0; $i<@marker ; $i++) {
		push @genotype , $genotype{$pos}{$marker[$i]}
	}

	$count += filter(\@genotype ,\%weight, \@genotype_new ) ;
#	print @genotype_new;<stdin>;
	for (my $i=0; $i<@marker ; $i++) {
		$genotype_new{$pos}{$marker[$i]} = $genotype_new[$i] ;
	}
}

my $sum = @marker * @pos ;
#
# output 
#
my %miss_num;
my %singleton_num = (
	"before"=> $count,
	"after" => 0 ,
	);
$miss_num{before} = 0;
$miss_num{after}  = 0;
open (Nloc,">","$fOut.loc") or die $!;
print Nloc join ("\n" , @head ),"\n\n";
open (Sta,">","$fOut.loc.detail") or die $!;

for (my $i=0; $i<@marker ; $i++) {
	print Nloc $marker[$i] ,"\t" ;
	print Sta $marker[$i],"\t" ;
	foreach my $pos (sort {$a <=> $b} @pos) {

		print Nloc $genotype_new{$pos}{$marker[$i]} ,"\t" ;
		$miss_num{before} ++ if ($genotype{$pos}{$marker[$i]} eq "-") ; 
		$miss_num{after} ++ if ($genotype_new{$pos}{$marker[$i]} eq "-") ; 

		if ($genotype_new{$pos}{$marker[$i]} ne $genotype{$pos}{$marker[$i]}) {
			print Sta $genotype{$pos}{$marker[$i]},"->",$genotype_new{$pos}{$marker[$i]} ,"\t";
		}else{
			print Sta $genotype{$pos}{$marker[$i]},"\t";
		}
	}
	print Nloc "\n";
	print Sta "\n";
}
close (Nloc) ;
close (Sta) ;

open (Diff,">","$fOut\_Diff.xls") or die $!;
print Diff "\tSum\tmiss\tsingl\tmiss\%\tsingl\%\n";
my $misretio = sprintf ("%.4f", $miss_num{before} / $sum * 100 );
my $sngtonretio =sprintf ("%.4f",  $singleton_num{before} / $sum * 100);
print Diff "before\t$sum\t$miss_num{before}\t$singleton_num{before}\t$misretio\%\t$sngtonretio\%\n" ;
$misretio =sprintf ("%.4f", $miss_num{after} / $sum * 100) ;
$sngtonretio = sprintf ("%.4f", $singleton_num{after} / $sum * 100) ;
print Diff "after\t$sum\t$miss_num{after}\t$singleton_num{after}\t$misretio\%\t$sngtonretio\%\n" ;
close (Diff) ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub filter {#\@array1 \%weight \@array2 
	my ($array1,$weight,$array2) = @_ ;
	my %Di;
	my $count=0;
	for (my $i=0; $i<@{$array1} ;$i++) {

		my $left_size  = $windows_size < $i ? $windows_size : $i ;
		my $right_size = $windows_size < @{$array1} - $i -1 ? $windows_size : @{$array1} - $i -1 ;

		my %count;
		for (my $j=$i-$left_size; $j<$i+$right_size ;$j++) {
			next if ($array1->[$j] eq "-" ) ;
			$count{ $array1->[$j] } ++  ;
		}

		my @type = keys %count ;
		if ($array1->[$i] ne "-") {
			my %stand;
			my $value;
			my $numerator;
			my $denominator;
			$stand{$array1->[$i]} = 0 ;
			foreach (@type) {
				 $stand{$_} = 1 if ($_ ne $array1->[$i]) ;
			}

			my %probility ;
			for (my $j=$i-$left_size; $j<$i+$right_size ;$j++) {
				next if ($j == $i) ;
				next if ($array1->[$j] eq "-" ) ;
				$numerator += $weight{abs($i-$j)} * $stand{$array1->[$j]} ;
				$probility{$array1->[$j]} += $weight{abs($i-$j)} * $stand{$array1->[$j]} ;
				$denominator += $weight{abs($i-$j)};
			}

			if (!defined $numerator || !defined $denominator || $denominator == 0) {
				$value = "X" ;
			}else{
				$value = $numerator / $denominator ;
			}

			if ( $value ne "X") {
				$Di{$i} = abs ($stand{$array1->[$i]} - $value );
				if ($Di{$i} > $threshold ) {
					my $max_type = (sort {$probility{$b} <=> $probility{$a}} keys %probility )[0];
					$array2->[$i] = $max_type ;
					$count++;
				}else{
					$array2->[$i] = $array1->[$i];
				}
			}else{
				$Di{$i} = "x" ;
				$array2->[$i] = $array1->[$i];
			}
		}else{

			my %stand;
			my $value;
			my $numerator;
			my $denominator;
			my %probility;
			my %value_type ;
			my %weight_type;

			foreach (@type) {
				$stand{$_} = 1 ;
			}
			my $min_loc = 0 ;my $sum_loc = 0 ;
			for (my $j=$i-$left_size; $j<$i+$right_size ;$j++) {
				next if ($j == $i) ;
				$sum_loc += $weight{abs($i-$j)} ;
				next if ($array1->[$j] eq "-" ) ;
				$min_loc += $weight{abs($i-$j)} ;
				$numerator += $weight{abs($i-$j)} * $stand{$array1->[$j]} ;
				$weight_type{$array1->[$j]} += $weight{abs($i-$j)} * $stand{$array1->[$j]} ;
			}

			foreach (@type) {
				$value_type{$_} = $weight_type{$_} / $numerator ;
			}

			my $max_probality = "X";
			foreach my $type (@type) {
				if ($value_type{$type} > $threshold ){
					$max_probality = $type ;
					last;
				}
			}
			$max_probality = "X" if ( $min_loc/$sum_loc < $min_mis ) ;
			if ($max_probality ne "X") {
				$array2->[$i] = $max_probality ;
			}else {
				$array2->[$i] = $array1->[$i] ;
			}
		}
	}
	return $count;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: wangml <wangml\@biomarker.com.cn> 
		

Usage: 根据具有map顺序的loc文件进行基因型纠错
  Options:
  -i <file>   input loc file  forced
  -o <file>   output file stem, xxx.loc,xxx.loc.detail, forced
  -D <float>  threshold values to change a loci, defaut [0.95]
  -M <float>  minimim miss ratio, [0.85]
  -h          Help

USAGE
	print $usage;
	exit;
}




