#!/usr/bin/perl -w
BEGIN{ #add absolute path of required moduals to array \@INC
	push @INC,"/share/nas1/macx/perlModuals/statistics/Regression/blib/lib";
}
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Regression;
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fMap,$fPwd,$fKey,$dOut,$mapFunction,$episilon, $rThreshold, $lodThreshold);
GetOptions(
				"help|?" =>\&USAGE,
				"m:s"=>\$fMap,
				"p:s"=>\$fPwd,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				
				"mapFunction:s"=>\$mapFunction,
				"r:s"=>\$rThreshold,
				"lod:s"=>\$lodThreshold,
				
				"episilon:s"=>\$episilon,
				) or &USAGE;
&USAGE unless ($fMap and $fPwd and $fKey);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

$mapFunction||='Kosambi';
$episilon||=0.0001;
$rThreshold||=0.4;
$lodThreshold||=1;

#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my (%Map,%pwd) = ();
my %mapFunction=(
	"Kosambi"=>\&Kosambi,
	"Haldane"=>\&Haldane,
	);
my %inverseMapFunction=(
	"Kosambi"=>\&inverseKosambi,
	"Haldane"=>\&inverseHaldane,
	);

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------

loadMap(\%Map,$fMap);
loadPwd(\%pwd,$fPwd);

#-------------------------------------------------------------------
# Process
#-------------------------------------------------------------------

my ($adjusted_order,$adjusted_dist) = ();

my @original_order = sort {$Map{$a} <=> $Map{$b}} keys %Map;
my $distance = weighted_linear_regression(\%pwd,\@original_order);

print join("\t",@{$distance}),"\n";

my ($min_dist) = sort {$a <=> $b} @{$distance};
print $min_dist,"\n";

if (defined $min_dist && $min_dist < -$episilon) {

	($adjusted_order,$adjusted_dist) = adjustOrder(\@original_order,$distance);


	open (OUT,">$dOut/$fKey.map") or die $!;
	print OUT "group 1\n\n";

	my $accum_dist = 0;

	print OUT $adjusted_order->[0],"\t",sprintf("%.3f",$accum_dist),"\n";

	for (my $i=1;$i<@{$adjusted_order} ;$i++) {

		$accum_dist += $adjusted_dist->[$i-1];
		print OUT $adjusted_order->[$i],"\t",sprintf("%.3f",$accum_dist),"\n";
	}

	close (OUT) ;
}else{

	`cp $fMap $dOut/$fKey.map`;
}

#print scalar @{$adjusted_order},"\n",scalar @{$adjusted_dist},"\n";
#-------------------------------------------------------------------
# Print
#-------------------------------------------------------------------

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub adjustOrder{#
	my ($ref_seq,$ref_dis)=@_;

	## calculate acccumutive map distance 
	my @accum_dist = ();
	$accum_dist[0] = 0;
	map {$accum_dist[$_] = $accum_dist[$_-1] + $ref_dis->[$_-1]} 1..@{$ref_seq}-1;

	## init order 
	my @order = ();
	map {$order[$_] = $_} 0..@{$ref_seq}-1;

	## adjust order according to accumutive distance 
	for (my $i=0;$i<@{$ref_seq} ;$i++) {
		for (my $j=$i+1;$j<@{$ref_seq} ;$j++) {

			if ($accum_dist[$i] > $accum_dist[$j]) {
				($order[$i],$order[$j]) = ($order[$j],$order[$i]);
				($accum_dist[$i],$accum_dist[$j]) = ($accum_dist[$j],$accum_dist[$i]);
			}
		}
	}

	## adjust map order 
	my (@newOrder,@newDist) = ();

	map {$newOrder[$_] = $ref_seq->[$order[$_]]}0..@{$ref_seq}-1;
	map {$newDist[$_] = $accum_dist[$_+1] - $accum_dist[$_]} 0..@{$ref_dis}-1;

	return (\@newOrder,\@newDist);

}

sub weighted_linear_regression {
	my ($pairWise,$markerList)=@_;
	
	my $variablesNum=@{$markerList}-1;
	my $obj=Statistics::Regression->new($variablesNum);
	
	for (my $i=0;$i<@{$markerList}-1;$i++) {
			
			my @coefficientVector=();
			my $const=0;

			for (my $j=0;$j<@{$markerList}-1;$j++) {

					my ($min,$max)=sort {$a <=> $b} ($i,$j);
					my $sum=0;
					
					for(my $h=0;$h<=$min;$h++){
						for (my $k=$max+1;$k<@{$markerList};$k++) {

							next if (!defined $pairWise->{$markerList->[$h]}{$markerList->[$k]});
							next unless ($pairWise->{$markerList->[$h]}{$markerList->[$k]}{'r'} <= $rThreshold && $pairWise->{$markerList->[$h]}{$markerList->[$k]}{'lod'} >= $lodThreshold) ;
#							next if (!defined $pairWise->{$markerList->[$h]}{$markerList->[$k]}) ;
							$sum+=$pairWise->{$markerList->[$h]}{$markerList->[$k]}{"lod"}**2;
						}
					}
					
					push @coefficientVector,$sum;
					
					if ($j<=$i) {
						for (my $k=$i+1;$k<@{$markerList};$k++) {
							next if (!defined $pairWise->{$markerList->[$j]}{$markerList->[$k]}) ;
							next unless ($pairWise->{$markerList->[$j]}{$markerList->[$k]}{'r'} <= $rThreshold && $pairWise->{$markerList->[$j]}{$markerList->[$k]}{'lod'} >= $lodThreshold) ;
							$const+=$pairWise->{$markerList->[$j]}{$markerList->[$k]}{"lod"}**2*(&{$mapFunction{$mapFunction}}($pairWise->{$markerList->[$j]}{$markerList->[$k]}{"r"}));
						}

					}
			}

#			print join("\t",@coefficientVector,$const),"\n";

			$obj->include($const,\@coefficientVector);
	}

#	print "==================================================\n";
	
	my @result=$obj->theta;
	return \@result;
}

sub weighted_linear_regression_optimised {#
	my ($pairWise,$markerList) = @_;

	my $variablesNum = @{$markerList}-1;
	my $obj = Statistics::Regression->new($variablesNum);

	my (%coef,%const) = ();

	for (my $i=0;$i<@{$markerList}-1 ;$i++) {
		for (my $j=$i+1;$j<@{$markerList} ;$j++) {			
			next if (not defined $pairWise->{$markerList->[$i]}{$markerList->[$j]}) ;
			for (my $k = $i;$k < $j ;$k++) {
				for (my $h = $i;$h < $j ;$h++) {
					$coef{$k}{$h} += $pairWise->{$markerList->[$i]}{$markerList->[$j]}{'lod'}**2;
#					$const{$k} += $pairWise->{$markerList->[$i]}{$markerList->[$j]}{'lod'}**$weight*(&{$mapFunction{$mapFunction}}($pairWise->{$markerList->[$i]}{$markerList->[$j]}{'r'}));
				}
				$const{$k} += $pairWise->{$markerList->[$i]}{$markerList->[$j]}{'lod'}**2*(&{$mapFunction{$mapFunction}}($pairWise->{$markerList->[$i]}{$markerList->[$j]}{'r'}));
			} 
		}
	}

	print scalar @{$markerList},"\t",scalar keys %coef,"\n";
	print Dumper %{$pairWise->{$markerList->[202]}};

	foreach my $equation (sort {$a <=> $b} keys %coef) {
		my @coef = map {$coef{$equation}{$_}} sort {$a <=> $b} keys %{$coef{$equation}};
		$obj->include($const{$equation},\@coef);
	}
	my @result=$obj->theta;
	return \@result;
}

sub loadMap {#
	my ($mref,$fIn) = @_;

	my $order = 0;

	open (IN,"$fIn") or die $!;
	$/="\n";
	while (<IN>) {
		chomp;
		next if (/^$/ or /^\#/ or /^group/ or /^;/) ;

		my ($marker,$cm) = split /\s+/,$_;

		$mref->{$marker} = ++$order;
	}
	close (IN) ;
}

sub loadPwd {#
	my ($pref,$fIn) = @_;

	open (IN,"$fIn") or die $!;
	$/="\n";
	while (<IN>) {
		chomp;
		next if (/^$/ or /^name/ or /^;/) ;
		last if (/^locus/) ;

		next unless (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) ;

		$pref->{$1}{$2}{'r'} = $3;
		$pref->{$2}{$1}{'r'} = $3;
		$pref->{$1}{$2}{'lod'} = $4;
		$pref->{$2}{$1}{'lod'} = $4;

	}
	
	close (IN) ;
}
sub Haldane {#
	my $r=shift;
	if ($r>=0.5) {
		$r = 0.4999;
	}
	my $result=-100*(1/2)*log(1-2*$r);
	return $result;
}

sub Kosambi {#
	my $r=shift;
	if ($r>=0.5) {
		$r = 0.4999;
	}
	my $result=100*(1/4)*log((1+2*$r)/(1-2*$r));
	return $result;
}


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription:

Usage:
  Options:
  -m	<file>		Map file,joinmap format, forced
  -p	<file>		Pwd file,joinmap format, forced

  -k	<str>		Key of output file,forced
  -d	<str>		Directory where output file produced,optional,default [./]
  
  -r	<double>	Recombination threshold, optional, default [0.4]
  -lod	<double>	Lod threshold, optional, default [1] 
  
  -episilon	<double>	Tolerace of negative distance, optional, default [0.0001] 

  -h			Help

USAGE
	print $usage;
	exit;
}