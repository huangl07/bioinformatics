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
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($locfile,$pwdfile,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"loc:s"=>\$locfile,
				"pwd:s"=>\$pwdfile,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($locfile and $pwdfile and $fOut );


$locfile = Cwd::abs_path($locfile) ;
$pwdfile = Cwd::abs_path($pwdfile) ;


#
# read in pwd 
#
open (PWD,$pwdfile) or die $!;
my %linph ;
my %LOD ;
my $pre_LOD = 0 ;

my $ini_marker1 ;
my $ini_marker2 ;
while (<PWD>) {
	chomp;
	next if (/^\s*$/ || /^;/);

	my ($marker1 ,$marker2 ,$recom , $LOD , $linph ) = split ;

	$LOD{$marker1}{$marker2} = $LOD ;
	$LOD{$marker2}{$marker1} = $LOD ;
	
	$linph{$marker1}{$marker2} = $linph ;
	$linph{$marker2}{$marker1} = $linph ;

	if ($LOD > $pre_LOD) {
		$ini_marker1 = $marker1 ;
		$ini_marker2 = $marker2 ;
		$pre_LOD = $LOD ;
	}
}
close (PWD) ;


##\\\\\\\\\\\\\\\\\\

my %abs_linph ;
my %redefind_ratio;

my @linph = ($ini_marker1,$ini_marker2) ;

$abs_linph{$ini_marker1}{linq} = 0 ;
$abs_linph{$ini_marker2}{linq} = 0 if ($linph{$ini_marker1}{$ini_marker2} == 0) ;
$abs_linph{$ini_marker2}{linq} = 1 if ($linph{$ini_marker1}{$ini_marker2} == 1) ;

open (LIN,">$fOut.dlp.detail") or die $!;
print LIN "initial Marker:\n",$ini_marker1,"\t",$abs_linph{$ini_marker1}{linq},"\n",$ini_marker2 ,"\t", $abs_linph{$ini_marker2}{linq},"\n";

while (1) {

	my @temp = keys %abs_linph ;
	my $cur_marker = Find_MAXLOD_marker(\@temp , \%LOD ) ;
	last if ( $cur_marker eq "0" ) ;

	my %temp_linph ;
	foreach my $marker2 ( keys %abs_linph ) {
		if ($linph{$cur_marker}{$marker2} eq "0") {
			push @{$temp_linph{ $abs_linph{$marker2}{linq}}} , $marker2 ;
		}else{
			push @{$temp_linph{ 1 - $abs_linph{$marker2}{linq}}} , $marker2 ;
		}
	}
	@{$temp_linph{"0"}} = () if (!defined $temp_linph{"0"} ) ;
	@{$temp_linph{"1"}} = () if (!defined $temp_linph{"1"} ) ;

	if ( @{$temp_linph{"0"}} >= @{$temp_linph{"1"}} ) {
		$abs_linph{$cur_marker}{linq} = 0 ;
		$abs_linph{$cur_marker}{"sup"} = \@{$temp_linph{"0"}} ;
		$abs_linph{$cur_marker}{"opp"} = \@{$temp_linph{"1"}} ;
	}else{
		$abs_linph{$cur_marker}{linq} = 1 ;
		$abs_linph{$cur_marker}{"sup"} = \@{$temp_linph{"1"}} ;
		$abs_linph{$cur_marker}{"opp"} = \@{$temp_linph{"0"}} ;
	}
	my $temp_ratio = scalar(@{$abs_linph{$cur_marker}{"sup"}}) / (scalar (@{$abs_linph{$cur_marker}{"opp"}}) + scalar(@{$abs_linph{$cur_marker}{"sup"}}) ) ;
	print LIN $cur_marker,"\t",$abs_linph{$cur_marker}{linq} , "\t" , $temp_ratio , "\n" ;
	print LIN "Support ",$cur_marker,":\t",join ( "\t" , @{$abs_linph{$cur_marker}{"sup"}}) ,"\n";
	print LIN "against ",$cur_marker,":\t",join ( "\t" , @{$abs_linph{$cur_marker}{"opp"}}) ,"\n";

	foreach my $temp_marker ( @{$abs_linph{$cur_marker}{"opp"}} ) {
		$redefind_ratio{$cur_marker}{$temp_marker} = 1;
		$redefind_ratio{$temp_marker}{$cur_marker} = 1;
	}
}
close (LIN) ;

#
# loc file 
#
my %marker_info ;
my $head ;

open (LOC,$locfile) or die $!;
open (NLOC,">$fOut.DH.loc") or die $!;
while (<LOC>) {
	chomp;

	if (/=/) {
		print NLOC $_,"\n";
		next ;
	}

	next if (/^\s*$/) ;
	my ($marker , @genotype ) = split ;
	next if (@genotype < 8) ;

	print NLOC $marker,"\t";
	if ( $abs_linph{$ini_marker1}{linq} eq "0" ) {
		print NLOC join( "\t" , @genotype ),"\n";
	}else{
		my @temp_gen;
		for (my $i = 0 ; $i < @genotype ; $i ++ ) {
			if ( $genotype[$i] eq "b" ) {
				$temp_gen[$i] = "a" ;
			}elsif ($genotype[$i] eq "a"){
				$temp_gen[$i] = "b" ;
			}else{
				$temp_gen[$i] = "-" ;
			}
		}
		print NLOC join("\t" , @temp_gen ),"\n" ;
	}

}
close (LOC) ;
close (NLOC) ;

#
# redefined pwd
#
open (PWD,$pwdfile) or die $!;
open (NPWD,">$fOut.DH.pwd") or die $!;

while (<PWD>) {
	chomp;
	next if (/^\s*$/ || /^;/);

	my ($marker1 ,$marker2 ,$recom , $LOD , $linph ) = split ;
	if (defined $redefind_ratio{$marker2}{$marker1}) {
		print NPWD $marker1 ,"\t",$marker2 ,"\t", "0.4999" ,"\t" ,0.0001 ,"\t;pre: $recom   $LOD\n" ;
	}else {
		print NPWD $marker1 , "\t" ,$marker2 ,"\t", $recom ,"\t" ,$LOD ,"\n" ;
	}
}
close (PWD);
close (NPWD);

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub Find_MAXLOD_marker {#\@ \%
	my ($array ,$LOD) = @_ ;
	if ( @{$array} == scalar (keys %{$LOD})) {
		return 0 ;
	}else{
		my %temp ;
		my $max_marker ;
		map {$temp{$_} ++ } @{$array} ;
		my $cur_LOD = Cal_LODS ( $array , $LOD ) ;
		foreach my $marker ( keys %{$LOD} ) {
			my @temp = @{$array};
			next if ( defined $temp{$marker} ) ;
			push @temp , $marker ;
			my $temp_LOD = Cal_LODS ( \@temp , $LOD ) ;
			if ($temp_LOD >= $cur_LOD ) {
				$max_marker = $marker ;
				$cur_LOD = $temp_LOD  ;
			}
		}
		return $max_marker ;
	}
}

sub Cal_LODS {#\@  \%
	my ($array ,$LOD) = @_ ;
	my $sum_LOD ;
	map {$sum_LOD += $LOD->{$array->[$_-1]}->{$array->[$_] } } 1..@{$array} - 1 ;
	return $sum_LOD ;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage: 按连锁相统一DH群体的个体分型，作为MSTmap的输入
  Options:
  -help               USAGE
  -loc                input loc file, forced 
  -pwd                input pwd file, forced 
  -o                  output file stem, xxx.DH.loc, xxx.DH.pwd, xxx.dlp.detail


USAGE
	print $usage;
	exit;
}

