#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
###################### ........loc.... 
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);

$fIn = Cwd::abs_path($fIn);

#
# loc 
#
my %marker_info ;
my @marker ;
my %mis_marker ;

open (LOC,$fIn) or die $!;
while (<LOC>) {
	chomp;
	next if (/^\s*$/ || /=/) ;

	my ($marker,@genotype) = split ;
	next if (@genotype < 8) ;
	for (my $i=0 ;$i<@genotype ;$i++) {
		$marker_info{$i}{$marker} = $genotype[$i];
		$mis_marker{$i}{$marker} = 1 if ($genotype[$i] eq "-") ;
	}
	push @marker , $marker ;
}
close (LOC) ;

#
# indi and marker numbers 
#

my $loc_num = scalar keys %{$marker_info{0}};
my $ind_num = scalar keys %marker_info;

#
# find noise genotypes：singleton
#
my %sus_marker ;
for (my $i=0; $i< keys %marker_info ;$i++) {
	for (my $j=0; $j < @marker ;$j++) {
		next if ($marker_info{$i}{$marker[$j]} eq "-" ) ;

		my $before = $marker_info{$i}{$marker[$j]} ;
		my $after  = $marker_info{$i}{$marker[$j]} ;
		my $sign=0;
		for (my $k =$j-1;$k>=0 ;$k--) {
			next if ($marker_info{$i}{$marker[$k]} eq "-" );
			if ($j == @marker - 1 ) {
				if ($sign eq 0) {
					$sign = $marker_info{$i}{$marker[$k]} ;
				}elsif ($sign eq "a" || $sign eq "b" || $sign eq "h") {
					$before = $sign ;
					$after  = $marker_info{$i}{$marker[$k]} ;
					last;
				}
			}else {
				$before = $marker_info{$i}{$marker[$k]} ;
				last;
			}
		}

		$sign =0 ;
		for (my $k =$j+1;$k< @marker ;$k++) {
			next if ($marker_info{$i}{$marker[$k]} eq "-" );
			if ($j == 0 ) {
				if ($sign eq 0) {
					$sign = $marker_info{$i}{$marker[$k]} ;
				}elsif ($sign eq "a" || $sign eq "b" || $sign eq "h") {
					$before = $sign;
					$after  = $marker_info{$i}{$marker[$k]} ;
					last;
				}
			}else {
				$after = $marker_info{$i}{$marker[$k]} ;
				last;
			}
		}
		
		if ($before eq $after && $marker_info{$i}{$marker[$j]} ne $before) {
			$sus_marker{$marker[$j]}{$i} = 1 ;
		}elsif ($before ne $after && $marker_info{$i}{$marker[$j]} ne $before && $marker_info{$i}{$marker[$j]} ne $after) {
			$sus_marker{$marker[$j]}{$i} = 1 ;
		}
	}
}

#
# statistic: missing observations 
#
my %mis;
my $mis;

for (my $i=0; $i< keys %marker_info ;$i++) {
	for (my $j=0; $j < @marker ;$j++) {
		if (defined $mis_marker{$i}{$marker[$j]} ) {
			$mis{indi}{$i}++ ;
			$mis{marker}{$marker[$j]}++ ;
			$mis++;
		}
	}
}

#
# statistic: suspicious 
#

my %sus;
my $sus;
for (my $i=0; $i< keys %marker_info ;$i++) {
	for (my $j=0; $j < @marker ;$j++) {
		if (defined $sus_marker{$marker[$j]}{$i} ) {
			$sus{indi}{$i}++ ;
			$sus{marker}{$marker[$j]}++ ;
			$sus++;
		}
	}
}



#
# output 
#
open ( STAD, ">" ,"$fOut.detail" ) or die $!;
print STAD "marker\ttotal\tmiss\tnoise_loc\tmiss\%\tnoise_loc\%\n";
for (my $i=0; $i<@marker ;$i++) {
	$mis{marker}{$marker[$i]} = 0 if (!defined $mis{marker}{$marker[$i]} ) ;
	$sus{marker}{$marker[$i]} = 0 if (!defined $sus{marker}{$marker[$i]} ) ;
	print STAD $marker[$i],"\t",$ind_num,"\t",$mis{marker}{$marker[$i]},"\t",$sus{marker}{$marker[$i]},"\t",$mis{marker}{$marker[$i]}/$ind_num,"\t",$sus{marker}{$marker[$i]}/$ind_num,"\n";
}

print STAD "\n\nIndi\ttotal\tmiss\tnoise_loc\tmiss\%\tnoise_loc\n";
for (my $i=0; $i< keys %marker_info ;$i++) {
	$mis{indi}{$i} = 0 if (!defined $mis{indi}{$i} ) ;
	$sus{indi}{$i} = 0 if (!defined $sus{indi}{$i} ) ;
	print STAD $i,"\t",$loc_num,"\t",$mis{indi}{$i},"\t",$sus{indi}{$i},"\t",$mis{indi}{$i}/$ind_num,"\t",$sus{indi}{$i}/$ind_num,"\n";
}
close (STAD) ;

##############

open ( STA, ">" ,"$fOut.xls" ) or die $!;
print STA "#total\tmiss\tnoise\tmiss\%\tnoise\%\n" ;

$sus = 0 if (!defined $sus ) ;
$mis = 0 if (!defined $mis ) ;

print STA $loc_num*$ind_num,"\t",$mis,"\t",$sus,"\t",$mis/($loc_num*$ind_num)*100,"\%\t",$sus/($loc_num*$ind_num)*100 ,"\%\n\n";
print STA "nloc = $loc_num\n";
print STA "nind = $ind_num\n";

my @nos_marker = (sort {$sus{marker}{$b} <=> $sus{marker}{$a}} keys %{$sus{marker}} )[0..4];
my @mis_marker = (sort {$mis{marker}{$b} <=> $mis{marker}{$a}} keys %{$mis{marker}} )[0..4];
my @nos_indi = (sort {$sus{indi}{$b} <=> $sus{indi}{$a}} keys %{$sus{indi}} )[0..4];
my @mis_indi = (sort {$mis{indi}{$b} <=> $mis{indi}{$a}} keys %{$mis{indi}} )[0..4];

print STA "\n\#most noise markers:\n" ;
map { print STA "\#$_\t$sus{marker}{$_}\n" } @nos_marker ;
print STA "\n\#most miss  markers:\n" ;
map { print STA "\#$_\t$mis{marker}{$_}\n" } @mis_marker ;

print STA "\n\#most noise indis:\n" ;
map { print STA "\#$_\t$sus{indi}{$_}\n" } @nos_indi ;
print STA "\n\#most miss  indis:\n" ;
map { print STA "\#$_\t$mis{indi}{$_}\n" } @mis_indi ;

close (STA) ;
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

Usage: 对近交群体排好序的loc文件统计杂点及缺失情况并发现可疑marker和可疑个体
  Options:
  -help			USAGE
  -i			input ordered loc file, forced 
  -o			output file 

USAGE
	print $usage;
	exit;
}



