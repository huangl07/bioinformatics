#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$popt,$homo,$seg,$mis);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"popt:s"=>\$popt,
	"homo:s"=>\$homo,
	"seg:s"=>\$seg,
	"mis:s"=>\$mis,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $popt);
close In;
$homo||="aa";
$mis||=0.3;
$seg||=0.05;
open In,$fIn;
open Out,">$fOut.marker.stat";
print Out "#MarkerID\ttype\tNind\tNmiss\tGeno\tNGeno\tSegretion\n";
open Geno,">$fOut.filtered.marker";
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ );
	if (/^#/) {
		print Geno "$_\n";
		next;
	}
	my ($id,$type,@info)=split(/\s+/,$_);
	$stat{total}++;
	if ($type eq "aaxbb" && $popt eq "BC1") {
		for (my $i=0;$i<@info;$i++) {
			if ($info[$i] eq "bb") {
				$info[$i] = "--";
			}
		}
	}
	if ($type eq "abxcc") {
		$type = "lmxll";
		for (my $i=0;$i<@info;$i++) {
			if ($info[$i] eq "ac") {
				$info[$i]="ll"
			}elsif ($info[$i] eq "bc") {
				$info[$i]="lm"
			}elsif ($info[$i] eq "ab"||$info[$i] eq "cc") {
				$info[$i]="--"
			}
		}
	}elsif ($type eq "ccxab") {
		$type = "nnxnp";
		for (my $i=0;$i<@info;$i++) {
			if ($info[$i] eq "ac") {
				$info[$i]="np"
			}elsif ($info[$i] eq "bc") {
				$info[$i]="nn"
			}elsif ($info[$i] eq "ab"|| $info[$i] eq "cc") {
				$info[$i]="--"
			}
		}
	}
	my %gstat;
	my $miss=0;
	for (my $i=0;$i<@info;$i++){
		if ($info[$i] eq "ff" || $info[$i] eq "gg" ||$info[$i] eq "mm" || $info[$i] eq "pp") {
			$info[$i] = "--";
		}elsif (($info[$i] eq "ab"||$info[$i] eq "cd") && $type eq "abxcd") {
			$info[$i] = "--";
		}
		if ($info[$i] eq "--"){
			$miss++;
			next;
		}
		$gstat{$info[$i]}++;
	}
	my $flag;
	my $order;
	my($p)=&SegregationX2($type,\@info,\$flag,\$order);
	print Out "$id\t$type\t",scalar @info,"\t",$miss,"\t",$order,"\t",$flag,"\n";
	if ($miss / scalar @info > $mis) {
		$stat{miss}++;
		next;
	}
	if (($popt eq "CP" && $type eq "aaxbb") || ($popt ne "CP" && $type ne "aaxbb")) {
		$stat{type}++;
		next;
	}
	if ($p < $seg) {
		$stat{seg}++;
		next;
	}
	$stat{filtered}++;
	print Geno "$id\t$type\t",join("\t",@info),"\n";
}
close Out;
close In;
close Geno;
open Out,">$fOut.filtered.stat";
print Out "#type\tnum\tpercent\n";
foreach my $seq (sort keys %stat) {
	print Out $seq,"\t",$stat{$seq},"\t",$stat{$seq}/$stat{total},"\n";
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub SegregationX2{#.....
	my($type,$genotype,$flag,$order)=@_;
	my %Data;
	foreach  (@$genotype) {
		$Data{$_}++;
	}
	if($type eq "hkxhk"){
		my ($hh,$hk,$kk,$missingData);
		$hh=exists $Data{"hh"}?$Data{"hh"}:0;
		$hk=exists $Data{"hk"}?$Data{"hk"}:0;
		$kk=exists $Data{"kk"}?$Data{"kk"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$hh+$hk+$kk+$missingData;
		my $valid=$hh+$hk+$kk;
		my $genotypeOrder="hh:hk:kk";
		my $theoretical_segregation="1:2:1";
		my $segregation="$hh:$hk:$kk";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "efxeg"){
		my ($eg,$ef,$ee,$fg,$missingData);
		$eg=exists $Data{"eg"}?$Data{"eg"}:0;
		$ef=exists $Data{"ef"}?$Data{"ef"}:0;
		$ee=exists $Data{"ee"}?$Data{"ee"}:0;
		$fg=exists $Data{"fg"}?$Data{"fg"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$eg+$ef+$ee+$fg+$missingData;
		my $valid=$eg+$ef+$ee+$fg;
		my $genotypeOrder="ee:ef:eg:fg";
		my $theoretical_segregation="1:1:1:1";
		my $segregation="$ee:$ef:$eg:$fg";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "abxcd"){
		my ($ac,$ad,$bc,$bd,$missingData);
		$ac=exists $Data{"ac"}?$Data{"ac"}:0;
		$ad=exists $Data{"ad"}?$Data{"ad"}:0;
		$bc=exists $Data{"bc"}?$Data{"bc"}:0;
		$bd=exists $Data{"bd"}?$Data{"bd"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ac+$ad+$bc+$bd+$missingData;
		my $valid=$ac+$ad+$bc+$bd;
		my $genotypeOrder="ac:ad:bc:bd";
		my $theoretical_segregation="1:1:1:1";
		my $segregation="$ac:$ad:$bc:$bd";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "lmxll"){
		my ($ll,$lm,$missingData);
		$ll=exists $Data{"ll"}?$Data{"ll"}:0;
		$lm=exists $Data{"lm"}?$Data{"lm"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ll+$lm+$missingData;
		my $valid=$ll+$lm;
		my $genotypeOrder="lm:ll";
		my $theoretical_segregation="1:1";
		my $segregation="$lm:$ll";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "nnxnp"){
		my ($nn,$np,$missingData);
		$nn=exists $Data{"nn"}?$Data{"nn"}:0;
		$np=exists $Data{"np"}?$Data{"np"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$nn+$np+$missingData;
		my $valid=$nn+$np;
		my $genotypeOrder="nn:np";
		my $theoretical_segregation="1:1";
		my $segregation="$nn:$np";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "abxcc"){
		my ($ac,$bc,$missingData);
		$ac=exists $Data{"ac"}?$Data{"ac"}:0;
		$bc=exists $Data{"bc"}?$Data{"bc"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ac+$bc+$missingData;
		my $valid=$ac+$bc;
		my $genotypeOrder="ac:bc";
		my $theoretical_segregation="1:1";
		my $segregation="$ac:$bc";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "ccxab"){
		my ($ac,$bc,$missingData);
		$ac=exists $Data{"ac"}?$Data{"ac"}:0;
		$bc=exists $Data{"bc"}?$Data{"bc"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ac+$bc+$missingData;
		my $valid=$ac+$bc;
		my $genotypeOrder="ac:bc";
		my $theoretical_segregation="1:1";
		my $segregation="$ac:$bc";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "aaxbb"){
		my ($aa,$bb,$ab,$missingData);
		$aa=exists $Data{"aa"}?$Data{"aa"}:0;
		$bb=exists $Data{"bb"}?$Data{"bb"}:0;
		$ab=exists $Data{"ab"}?$Data{"ab"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$aa+$bb+$ab+$missingData;
		my $valid=$aa+$bb+$ab;
		my $genotypeOrder="ab:aa:bb";
		if($popt=~/BC\d+/i){
			$genotypeOrder=$homo=~/aa/i?"aa:ab":"bb:ab";
		}elsif ($popt=~/Ri\d+/i || $popt=~/DH/i) {
			$genotypeOrder="aa:bb";
		}
		my $theoretical_segregation;
		if($popt=~/f2/i || $popt=~/F2/i){
			$theoretical_segregation="2:1:1";
		}elsif ($popt=~/Ri\d+/i || $popt=~/DH/i) {
			$theoretical_segregation="1:1";
		}elsif($popt=~/BC\d+/i){
			$theoretical_segregation="1:1"; 
		}elsif ($popt = "CP") {
			my $segregation="$ab:$aa:$bb";
			$$flag="$segregation";
			$$order="$genotypeOrder\t$segregation";
			return("-");
		}else{
			warn "please input true group!\n";
			exit(0);
		}
		my $segregation="$ab:$aa:$bb";
		if($popt=~/bc\d+/i){
			$segregation=$homo=~/aa/i?"$aa:$ab":"$bb:$ab";
		}elsif ($popt=~/Ri\d+/i || $popt=~/DH/i) {
			$segregation="$aa:$bb";
		}
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}else{#.......
		warn "wrong grouptype! $type\n";
		exit(0);
	}

	sub Segregation {#
		my ($theoretical_segregation,$segregation,$all)=@_;
		my @a=split ":",$theoretical_segregation;
		my @b=split ":",$segregation;
		return "0.01" if (scalar @a != scalar @b || $all == 0) ;
		my @theoretical;
		my $a_sum=0;
		$a_sum+=$_ foreach (@a);
		push @theoretical,$_/$a_sum*$all foreach (@a);
		my $df=scalar @a -1;
		my $X2=0;
		if ($df == 1) {
			for (my $i=0;$i<@a ;$i++) {
				$X2+=X2df2($b[$i],$theoretical[$i]);
			}
		}else{
			for (my $i=0;$i<@a ;$i++) {
				$X2+=X2df1($b[$i],$theoretical[$i]);
			}
		}
		my $p=0;
		$p=Statistics::Distributions::chisqrprob($df,$X2);
		return int($p*10000)/10000;
	}

	sub X2df1 {#
		my ($A,$T)=@_;
		return ($A-$T)**2/$T;
	}

	sub X2df2 {#
		my ($A,$T)=@_;
		return (abs($A-$T)-0.5)**2/$T;
	}
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -input	<file>	input file name
  -output	<file>	input keys of file name
  -popt	<str>	population type CP/BCi/Fi/Rix/
  -homo	<str>	if BC homotype,default=aa
  -seg	<num>	p value of segragation,default 0.05
  -mis	<num>	missing values  default 0.3
  -h         Help

USAGE
        print $usage;
        exit;
}
