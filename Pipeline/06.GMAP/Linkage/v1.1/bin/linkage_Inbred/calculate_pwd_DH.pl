#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$type);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"t:s"=>\$type,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $type);

my @markers ;
my %marker_info;

#
# read loc 
#
open (LOC,$fIn) or die $!;
my $head ;
while (<LOC>) {
	chomp;
	if (/=/) {
		$head .= "$_\n" ;
	}
	next if (/=/ || /^\s*$/) ;
	my ($marker,@genotype ) = split ;
	next if (@genotype < 8);
	push @markers ,$marker ;
	$marker_info{$marker} = \@genotype ;
}
close (LOC) ;

#
# calculate pair wise data 
#
my %linph ;
my $ini_marker1 ;
my $ini_marker2 ;
my %LOD;
my $pre_LOD = 0;
open (TT,">$fOut.pwd") or die $!;
for (my $i = 0; $i < @markers ;$i ++ ) {
	for (my $j = $i+1; $j < @markers ;$j ++) {
		my $array1 = $marker_info{$markers[$i]} ;
		my $array2 = $marker_info{$markers[$j]} ;
		my $comb ;
		my $LOD  ;
		my $linph ;
		if($type=~/[Dd][Hh]/){
			&calcomb( $array1,$array2,\$comb,\$linph);
			$linph{$markers[$i]}{$markers[$j]} = $linph ;
			$linph{$markers[$j]}{$markers[$i]} = $linph ;
		}
		if($type=~/BC\d+/i){
#			my @array1_com=@$array1;
#			my @array2_com=@$array2;
#			for (my $x=0;$x<@{$array1};$x++) {
#				{
#					$array1_com[$x]=~ s/\ba\b/aa/g; 
#					$array1_com[$x]=~ s/\bb\b/bb/g; 
#					$array1_com[$x]=~ s/\bh\b/ab/g; 
#					$array1_com[$x]=~ s/\bc\b/a-/g; 
#					$array1_com[$x]=~ s/\bd\b/-b/g; 
#					$array1_com[$x]=~ s/-/--/g;     
#				}
#				{
#					$array2_com[$x]=~ s/\ba\b/aa/g; 
#					$array2_com[$x]=~ s/\bb\b/bb/g; 
#					$array2_com[$x]=~ s/\bh\b/ab/g; 
#					$array2_com[$x]=~ s/\bc\b/a-/g; 
#					$array2_com[$x]=~ s/\bd\b/-b/g; 
#					$array2_com[$x]=~ s/-/--/g;
#				}
#			}
			$comb=Gamete_BC1(\@$array1,\@$array2,\$linph);
			$linph{$markers[$i]}{$markers[$j]} = $linph ;
			$linph{$markers[$j]}{$markers[$i]} = $linph ;
		}
		my $gen1 ;
		my $gen2 ;

		for (my $i = 0 ;$i<@{$array1} ;$i++) {
			if($type=~/[Dd][Hh]/){
				$gen1 -> [$i] = "aa" if ($array1->[$i] eq "a") ;
				$gen1 -> [$i] = "bb" if ($array1->[$i] eq "b") ;
				$gen1 -> [$i] = "--" if ($array1->[$i] eq "-") ;
				$gen2 -> [$i] = "aa" if ($array2->[$i] eq "a") ;
				$gen2 -> [$i] = "bb" if ($array2->[$i] eq "b") ;
				$gen2 -> [$i] = "--" if ($array2->[$i] eq "-") ;
			}else{															#马立祥
				$gen1 -> [$i] = "XX" if ($array1->[$i] eq "a") ;		#马立祥
				$gen1 -> [$i] = "XX" if ($array1->[$i] eq "b") ;		#马立祥
				$gen1 -> [$i] = "ab" if ($array1->[$i] eq "h") ;		#马立祥
				$gen1 -> [$i] = "--" if ($array1->[$i] eq "-") ;		#马立祥
				$gen2 -> [$i] = "XX" if ($array2->[$i] eq "a") ;		#马立祥
				$gen2 -> [$i] = "XX" if ($array2->[$i] eq "b") ;		#马立祥
				$gen2 -> [$i] = "ab" if ($array2->[$i] eq "h") ;		#马立祥
				$gen2 -> [$i] = "--" if ($array2->[$i] eq "-") ;		#马立祥
			}
		}																#马立祥
		$LOD = &calculatemLOD( "aaxbb" , "aaxbb" , $gen1 , $gen2  ) ;
		$LOD{$markers[$i]}{$markers[$j]} = $LOD ;
		$LOD{$markers[$j]}{$markers[$i]} = $LOD ;

		if ($LOD > $pre_LOD) {
			$ini_marker1 = $markers[$i] ;
			$ini_marker2 = $markers[$j] ;
			$pre_LOD = $LOD ;
		}
		print TT $markers[$i],"\t",$markers[$j],"\t",$comb,"\t", $LOD ,"\t",$linph , "\n" ;
	}
}
close (TT) ;

#
# unify liankge phase and rewrite loc file 
#
my $comand  = "perl $Bin/unify_linkage_phase_DH.pl -loc $fIn -pwd $fOut.pwd -o $fOut " ;
`$comand ` ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub calcomb {#\@ \@ \$theta \$lin
	my ($array1,$array2,$theta,$linph) = @_;
	my $bd = 0 ; my $ac = 0;
	for (my $i=0; $i<@{$array1} ;$i++) {
		next if ( $array1->[$i] eq "-" || $array2->[$i] eq "-" ) ;
		if ($type=~/[Dd][Hh]/) {
			if ( $array2->[$i] ne $array1->[$i] ) {
				$bd++ ;
			}else{
				$ac++ ;
			}
		}else{
			if ( ($array2->[$i] eq "h" && $array1->[$i] ne "h") || ($array2->[$i] ne "h" && $array1->[$i] eq "h")) {
				$bd++ ;
			}else{
				$ac++ ;
			}
		}
	}
	$$theta = $bd / ($ac+$bd) ;

	if ($$theta == 0.5) {
		$$theta = 0.499;
	}
	if ($$theta > 0.5 ) {
		$$theta = 1- $$theta;
		$$linph = 1;
	}else{
		$$linph = 0;
	}
}

sub calculatemLOD {#
	my ($type1,$type2,$indi1,$indi2) = @_;
	my ($aPgeno,$aMgeno) = split /x/,$type1;
	my ($bPgeno,$bMgeno) = split /x/,$type2;

	return -1 if ((substr($aPgeno,0,1) eq substr($aPgeno,1,1) and substr($aMgeno,0,1) ne substr($aMgeno,1,1) and substr($bPgeno,0,1) ne substr($bPgeno,1,1) and substr($bMgeno,0,1) eq substr($bMgeno,1,1)) or (substr($aPgeno,0,1) ne substr($aPgeno,1,1) and substr($aMgeno,0,1) eq substr($aMgeno,1,1) and substr($bPgeno,0,1) eq substr($bPgeno,1,1) and substr($bMgeno,0,1) ne substr($bMgeno,1,1))) ;

	die "different population size\n" if (@{$indi1} != @{$indi2});

	my (%pairs,%R,%C) = ();
	
	my $T = 0;
	for (my $p=0;$p<@{$indi2} ;$p++) {
		if ($indi1->[$p] eq "--" or $indi2->[$p] eq "--") {  ## skip uninformative indi 
			next;
		}else{
			$R{"$indi1->[$p]"}++;
			$C{"$indi2->[$p]"}++;
			$pairs{"$indi1->[$p]x$indi2->[$p]"}++;
			$T++;
		}

	}

	my ($bG,$G);
	foreach my $key1 (keys %pairs) {
		my ($r,$c)=split /x/, $key1;
		my $R=$R{"$r"};
		my $C=$C{"$c"};
		my $E=$R*$C/$T;
		$bG+=$pairs{$key1}*log($pairs{$key1}/$E);
	}
	if (!defined $bG) {
			$G=0;
	}else{
		$G=2*$bG;
	}
	
	my ($G1,$e);
	if (scalar keys %R < 2 or scalar keys %C < 2) {
		print "lack of freedom to calculate mLOD\n";
	}

	if (scalar keys %R == 2 and scalar keys %C == 2) {
		$G1=$G;
	}else{
		my $df = 2 ;			#DHqun
		$e = exp(-$G /(2 * ($df -1)));
		$G1 = ((4-$e)*$e-3) * ($df-1) + $G ;
	}
	my $mlod=$G1/ ( 2 * log (10) );
	return $mlod;
}

sub Gamete_BC1{ #
	my ($r_1,$r_2,$linph)=@_;
	my %Stat;
	my $Sum=0;
	for (my $i=0;$i<@$r_1;$i++) {
		next if ($$r_1[$i] eq "-" || $$r_2[$i] eq "-");
		$Stat{"$$r_1[$i]$$r_2[$i]"}++;
		$Sum++;
	}
	
	$Stat{aa}||=0;
	$Stat{ab}||=0;
	$Stat{ah}||=0;
	$Stat{ba}||=0;
	$Stat{bb}||=0;
	$Stat{bh}||=0;
	$Stat{hh}||=0;
	$Stat{ha}||=0;
	$Stat{hb}||=0;
	$Stat{"a-"}||=0;
	$Stat{"b-"}||=0;
	$Stat{"h-"}||=0;
	$Stat{"-a"}||=0;
	$Stat{"-b"}||=0;
	$Stat{"-h"}||=0;
	my $r=0.00001;
	#my $maxLod=0;
	my $maxR=0;
	my $t=$r;
	#my $Lods5 = ($Stat{aa}+$Stat{hh})*log((1-$t)/2)+($Stat{ah}+$Stat{ha})*log($t/2);#+($Stat{ah}+$Stat{bh}+$Stat{ha}+$Stat{hb})*lob($t*(1-$t)/2)+$Stat{hh}*log((1-2*$t+2*$t*$t)/2);
	my $maxLod = ($Stat{aa}+$Stat{hh}+$Stat{bb})*log((1-$t)/2)+($Stat{ah}+$Stat{ha}+$Stat{bh}+$Stat{hb})*log($t/2)-$Sum*log(0.25);

	while (1) {
		$t=$r;
		my $Lods = ($Stat{aa}+$Stat{hh}+$Stat{bb})*log((1-$t)/2)+($Stat{ah}+$Stat{ha}+$Stat{bh}+$Stat{hb})*log($t/2)-$Sum*log(0.25);#+($Stat{ah}+$Stat{bh}+$Stat{ha}+$Stat{hb})*lob($t*(1-$t)/2)+$Stat{hh}*log((1-2*$t+2*$t*$t)/2);
		if ($Lods == 0){
			$r+=0.00001;
			next;
		}
		#$Lods=($Lods-$Lods5)/log(10);
		if ( $Lods >= $maxLod) {
			$maxR=$r;
			$maxLod=$Lods;
		}
		$r+=0.00001;
		last if ($Lods < $maxLod || $r == 1);
	}
#	if ($maxR == 0 || $maxR > 0.5) {
#		$maxR=0.5;
#		$maxLod=0.001;
#	}

	if ($maxR == 0.5) {
		$maxR = 0.499;
	}
	if ($maxR > 0.5 ) {
		$maxR = 1- $maxR;
		$$linph = 1;
	}else{
		$$linph = 0;
	}

	return $maxR;
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

Usage: 计算DH群体pwd并将原来的loc按统一连锁相输出
  Options:
  -i	<file>	Input loc file, forced
  -o	<file>	Output file stem, forced    
  -t	popt
USAGE
	print $usage;
	exit;
}



