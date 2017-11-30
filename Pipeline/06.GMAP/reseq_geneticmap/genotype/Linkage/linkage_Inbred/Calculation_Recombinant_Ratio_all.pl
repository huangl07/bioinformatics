#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Math::Cephes::Polynomial;
my $BEGIN_TIME=time();
my $version="1.0.0";
$Script=~s/\.pl//g;
my @Times=localtime();
my $year=$Times[5]+1990;
my $month=$Times[4]+1;
my $day=$Times[3];
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$Type,$fOut,$OnlyFirst);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"p:i"=>\$OnlyFirst,
				"t:s"=>\$Type,
				) or &USAGE;
&USAGE unless ($fIn and $Type and $fOut);
&USAGE unless ($Type ne 'F2' or $Type ne 'BC1' or $Type =~ /[Rr][Ii]\d+/);


$OnlyFirst = defined $OnlyFirst ? $OnlyFirst: 0;
$Type||="F2";

#------------------------------------------------------------------
#  Global variable
#------------------------------------------------------------------

my (%Marker,%Data) = ();

#------------------------------------------------------------------
#  Get genotype or loc file 
#------------------------------------------------------------------

my $file_type = `grep ^nloc $fIn` ? 'loc' : 'genotype';
print "file type : $file_type\n";

open (IN,"$fIn") or die $!;
<IN> if ($file_type eq 'genotype') ; ## skip first line of genotype file 
while (<IN>) {

	next if (/^$/ || /^;/ || /^nloc/ || /^name/ || /^popt/ || /^nind/) ;
	
	my $marker = ''; ## marker ID 
	my $type = '';   ## cross_type 
	my @indi = ();   ## progeny genotypes 
	
	if ($file_type eq 'loc') {
		($marker,@indi) = split; 
	}else{
		($marker,$type,@indi) = split;
	}

	my $line=join("\t",@indi);
	{
		$line=~ s/\ba\b/aa/g;
		$line=~ s/\bb\b/bb/g;
		$line=~ s/\bh\b/ab/g;
		$line=~ s/\bc\b/a-/g;
		$line=~ s/\bd\b/-b/g;
		$line=~ s/-/--/g;
	}
	@indi=split(/\t/,$line);

	@{$Marker{'ORI'}{$marker}{"arr"}} = @indi;
	$Marker{'ORI'}{$marker}{"type"} = 'abxab'; ## for calculate mLOD 

	$line=join("\t",@indi);
	{
		$line=~ s/\baa\b/a/g;
		$line=~ s/\bbb\b/b/g;
		$line=~ s/\bab\b/h/g;
		$line=~ s/\ba-/c/g;
		$line=~ s/-b\b/d/g;
		$line=~ s/--/-/g;
	}
	@indi=split(/\t/,$line);
	
	@{$Marker{'MOD'}{$marker}{"arr"}} = @indi;
	$Marker{'MOD'}{$marker}{"type"} = 'abxab';

}
close IN;

#print Dumper %Marker;die;

#------------------------------------------------------------------
#  calculate recombination raction and modified LOD score 
#------------------------------------------------------------------

my @Marker=sort keys %{$Marker{'MOD'}};
my $EndMarkerCurr=$OnlyFirst?1:$#Marker;
if ($OnlyFirst ==0) {
	$EndMarkerCurr=$#Marker;
}else{
	$EndMarkerCurr=$OnlyFirst;
}
if ($EndMarkerCurr > scalar @Marker) {
	$EndMarkerCurr=scalar @Marker;
}

if ($Type =~ /[Rr][Ii]/ && $Type !~ /[Rr][Ii]2/){
	my ($Grade) = $Type =~/^[Rr][Ii](\d+)$/;
	$Grade||=10;
	my %Equa;

	&Parameter_RiLs($Grade,\%Equa);
	
#	print Dumper %Equa;die;

	# \\\\\\\\ debug:检验概率推导是否与黄龙结果一致 \\\\\\\\

#	for (my $r=0.00001;$r<=0.5;$r+=0.00001) {	
#		print $r;
#		for(1..6){
#			print "\t",$Equa{$_}->eval($r);
#		}
#		print "\n";
#	}
#	die;


	## 经检验，与黄龙计算结果一致

	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	open Out,">$fOut";
	for (my $i=0;$i<$EndMarkerCurr;$i++) {
		for (my $j=$i+1;$j<@Marker;$j++) {

#			print ">$Marker[$i]\t$Marker[$j]\n";
		
			my $Ref = Gamete_RiLs(\@{$Marker{'MOD'}{$Marker[$i]}{"arr"}},\@{$Marker{'MOD'}{$Marker[$j]}{"arr"}},\%Equa);

#			print "calculate recombination fraction done\n";
			my $mlod = calculate_mlod(\%{$Marker{'ORI'}},[$Marker[$i],$Marker[$j]]);

#			print "calculate modified LOD score done\n";
			if ($Ref >= 0.5) {	
				print Out $Marker[$i],"\t",$Marker[$j],"\t",sprintf("%.16e","0.4999"),"\t",sprintf("%.16e","0.001"),"\t","; original estimate: ",sprintf("%.16e",$Ref),"\t",sprintf("%.16e",$mlod),"\n";

			}else{
				
				print Out $Marker[$i],"\t",$Marker[$j],"\t",sprintf("%.16e",$Ref),"\t",sprintf("%.16e",$mlod),"\n";
			}
		}	
	}
	close Out;
}elsif ($Type eq "F2" || $Type =~/[Rr][Ii]2/) {

    #print "ok\n";die;

	open Out,">$fOut";
	for (my $i=0;$i<$EndMarkerCurr;$i++) {
		for (my $j=$i+1;$j<@Marker;$j++) {

#			print ">$Marker[$i]\t$Marker[$j]\n";
			my $Ref = Gamete_F2(\@{$Marker{'MOD'}{$Marker[$i]}{"arr"}},\@{$Marker{'MOD'}{$Marker[$j]}{"arr"}});
			my $mlod = calculate_mlod(\%{$Marker{'ORI'}},[$Marker[$i],$Marker[$j]]);
			if ($Ref >= 0.5) {
				
				print Out $Marker[$i],"\t",$Marker[$j],"\t",sprintf("%.16e","0.4999"),"\t",sprintf("%.16e","0.001"),"\t","; original estimate: ",sprintf("%.16e",$Ref),"\t",sprintf("%.16e",$mlod),"\n";

			}else{
				
				print Out $Marker[$i],"\t",$Marker[$j],"\t",sprintf("%.16e",$Ref),"\t",sprintf("%.16e",$mlod),"\n";
			}
		}
	}
	close Out;
}elsif ($Type eq "BC1") {

	open Out,">$fOut";
	for (my $i=0;$i<$EndMarkerCurr;$i++) {
		for (my $j=$i+1;$j<@Marker;$j++) {
			my ($Ref,$Lod);
			($Ref,$Lod) = Gamete_BC1(\@{$Marker{'MOD'}{$Marker[$i]}{"arr"}},\@{$Marker{'MOD'}{$Marker[$j]}{"arr"}});
			my $mlod = calculate_mlod(\%{$Marker{'ORI'}},[$Marker[$i],$Marker[$j]]);
			if ($Ref >= 0.5) {
				$Ref = 0.4999;
				$mlod = 0.0001;
			}
			print Out $Marker[$i],"\t",$Marker[$j],"\t",$Ref,"\t",$mlod,"\n";
		}
	}
	close Out;
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
### \\\\\\\\\\\ Root finding for nonlinear sets of equations \\\\\\\\\\\\\\\\\\\\\\\\\

sub bisectionMethod {# $eqution is an object of Math::Cephes::Polynomial
	my ($equation,$leftEnd,$rightEnd,$stepSize) = @_;
	
	my @root = ();
	$stepSize ||= 0.1;
	my $precision = 1e-8;
	my $root;

	my %interval = ();
	my $i = $leftEnd;
	my $time = 1;

	while (scalar keys %interval == 0 && $time++ < 100000) {
		
		$i = $leftEnd;
		while ($i < $rightEnd) {

			push @root,$i if ($equation->eval($i) == 0) ;
		
			$interval{$i."to".($i+$stepSize)} = 1 if ($equation->eval($i)*$equation->eval($i+$stepSize) < 0) ;

			$i+=$stepSize;
	
		}

		$stepSize/=10; 
			
	}

#	print join("\t",keys %interval),"\n";

#	<STDIN>;

	if (scalar keys %interval != 0) {
	
		foreach my $interval (keys %interval) {

			my @root_interval = split /to/,$interval;

			while (1) {

				if ($root_interval[1] - $root_interval[0] <= $precision) {

					push @root,($root_interval[0]+$root_interval[1])/2;
					last;
				}

				my $midPoint = ($root_interval[0]+$root_interval[1])/2;

#				print "midPoint: ",$midPoint,"\n";
#				<STDIN>;

				if ($equation->eval($midPoint) == 0) {
					
					push @root,$midPoint;
					last;
				}else{

					if ($equation->eval($midPoint)*$equation->eval($root_interval[0]) > 0) {

						$root_interval[0] = $midPoint;
						$root = $midPoint;

					}else{
						
						$root_interval[1] = $midPoint;
						$root = $midPoint;
					}
				}

			}
		}
	}else{

		print "root is not at interval [$leftEnd,$rightEnd].\n";

	}

#	print "roots :",join("\t",@root),"\n";
#	<STDIN>;

	return \@root;

}

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
sub Gamete_F2{
	my ($r_1,$r_2)=@_;
	my %Stat;
	my $Sum=0;
	for (my $i=0;$i<@$r_1;$i++) {
		next if ($$r_1[$i] eq "-" or  $$r_2[$i] eq "-");
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
	
	my %prob = (
		'1' => [0.25,-0.5,0.25],
		'2' => [0,0,0.25],
		'3' => [0,0.5,-0.5],
		'4' => [0.5,-1,1],
		);

	my %groupByProb = (
		'1' => ['aa','bb'],
		'2' => ['ab','ba'],
		'3' => ['ah','bh','ha','hb'],
		'4' => ['hh'],
		);
	
	my %factor = () ;

	foreach my $nr (sort {$a <=> $b} keys %groupByProb) {

		my $count = 0;

		foreach my $phenoType (@{$groupByProb{$nr}}) {

			$count += $Stat{$phenoType};
		}

		$factor{$nr}{'count'} = $count;

		$factor{$nr}{'poly'} = Math::Cephes::Polynomial->new([map {$_/0.25} @{$prob{$nr}}]);
	}

	my $numerator = Math::Cephes::Polynomial->new([0]);

	foreach my $factor (keys %factor) {
				
		my @temp=();

		&numberByVector($factor{$factor}{'count'},&derivativePolynomial($factor{$factor}{'poly'}->coef),\@temp);
				
		my $term = Math::Cephes::Polynomial->new(\@temp);

		foreach my $restFactor (keys %factor) {

			next if ($restFactor eq $factor) ;

			$term  = $term->mul($factor{$restFactor}{'poly'});

		}

		$numerator = $numerator->add($term);
 	}

	my ($flag,$roots) = $numerator->rts();
	my @realRoots = ();
	
	for (my $i=2;$i<$numerator->{'n'} ;$i++) {

#		print $i,"\t",$roots->[$i]->r,"\t",$roots->[$i]->i,"\n";

		next if ($roots->[$i]->i != 0) ;

		push @realRoots,$roots->[$i]->r;
		
	} 

	my ($rec) = sort {$a <=> $b} @realRoots;

	return $rec;
}
sub Gamete_RiLs{
	my ($r_1,$r_2,$prob)=@_;
	my %Stat;
	my $Sum=0;
	for (my $i=0;$i<@$r_1;$i++) {
		next if ($$r_1[$i] eq "-" && $$r_2[$i] eq "-");
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

	my %groupByProb = (
		'1' => ['aa','bb'],
		'2' => ['ab','ba'],
		'3' => ['ah','bh','ha','hb'],
		'4' => ['hh'],
		);

	my %factor = () ;

	foreach my $nr (sort {$a <=> $b} keys %groupByProb) {

		my $count = 0;

		foreach my $phenoType (@{$groupByProb{$nr}}) {

			$count += $Stat{$phenoType};
		}

		$factor{$nr}{'count'} = $count;

		$factor{$nr}{'poly'} = $prob->{$nr};
	}

#	print Dumper %factor;die;

	my $numerator = Math::Cephes::Polynomial->new([0]);

	foreach my $factor (keys %factor) {
				
		my @temp=();

		&numberByVector($factor{$factor}{'count'},&derivativePolynomial($factor{$factor}{'poly'}->coef),\@temp);
				
		my $term = Math::Cephes::Polynomial->new(\@temp);

		foreach my $restFactor (keys %factor) {

			next if ($restFactor eq $factor) ;

			$term  = $term->mul($factor{$restFactor}{'poly'});

		}

		$numerator = $numerator->add($term);
 	}

#	print Dumper $numerator;die;

	my $coef = $numerator->coef;

	foreach  (@{$coef}) {

		if ($_ == 0) {
			shift @{$coef};
		}else{
			last;
		}
	}

	$numerator = Math::Cephes::Polynomial->new($coef);

#	print Dumper $numerator;die;

	my $rec;
	my ($flag,$roots) = $numerator->rts();

#	print $flag,"\n";die;  ### 破模块，50多次的多项式求根都解不了？shit!!!
	
	if ($flag == 0) {

		my @realRoots = ();
	    my @trivilRoots = ();
		for (my $i=0;$i<$numerator->{'n'} ;$i++) {

#			print $i,"\t",$roots->[$i]->r,"\t",$roots->[$i]->i,"\n";

#            <STDIN>;

			next if ($roots->[$i]->i != 0) ;
            push @realRoots,$roots->[$i]->r;

		} 
        
        ($rec) = sort {$a <=> $b} @realRoots;

	}else{

		my @root = @{&bisectionMethod($numerator,0,1,0.1)};

		if (@root != 0) {
			($rec) = sort {$a <=> $b} @root;
		}else{
			$rec = 0.4999;
		}
	}

	return $rec;

}
sub Gamete_BC1{ #
	my ($r_1,$r_2)=@_;
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
	my $maxLod=0;
	my $maxR=0;
	my $t=$r;
	#my $Lods5 = ($Stat{aa}+$Stat{hh})*log((1-$t)/2)+($Stat{ah}+$Stat{ha})*log($t/2);#+($Stat{ah}+$Stat{bh}+$Stat{ha}+$Stat{hb})*lob($t*(1-$t)/2)+$Stat{hh}*log((1-2*$t+2*$t*$t)/2);
	my $Lods5 = ($Stat{aa}+$Stat{hh}+$Stat{bb})*log((1-$t)/2)+($Stat{ah}+$Stat{ha}+$Stat{bh}+$Stat{hb})*log($t/2)-$Sum*log(0.25);
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
	return $maxR,$maxLod;
}
sub Parameter_RiLs{
	my ($Grade,$result)=@_;

	my %prob = ();
	$prob{'1'}{'1'} = Math::Cephes::Polynomial->new([0]); ## aa||bb
	$prob{'1'}{'2'} = Math::Cephes::Polynomial->new([0]); ## ab||ba
	$prob{'1'}{'3'} = Math::Cephes::Polynomial->new([0]); ## ah||ha||bh||hb
	$prob{'1'}{'4'} = Math::Cephes::Polynomial->new([1]); ## hh->CC
	$prob{'1'}{'5'} = Math::Cephes::Polynomial->new([0]); ## hh->CR
	
	for (my $g=2;$g<=$Grade ;$g++) {

		$prob{$g}{'1'} = Math::Cephes::Polynomial->new([0]);
		$prob{$g}{'2'} = Math::Cephes::Polynomial->new([0]);
		$prob{$g}{'3'} = Math::Cephes::Polynomial->new([0]);
		$prob{$g}{'4'} = Math::Cephes::Polynomial->new([0]);
		$prob{$g}{'5'} = Math::Cephes::Polynomial->new([0]);

		my $multiplier = Math::Cephes::Polynomial->new([1]);
		$prob{$g}{'1'} = $prob{$g}{'1'}->add($prob{$g-1}{'1'});
		$multiplier = Math::Cephes::Polynomial->new([0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'3'});
		$prob{$g}{'1'} = $prob{$g}{'1'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0.25,-0.5,0.25]);
		$multiplier = $multiplier->mul($prob{$g-1}{'4'});
		$prob{$g}{'1'} = $prob{$g}{'1'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0,0,0.25]);
		$multiplier = $multiplier->mul($prob{$g-1}{'5'});
		$prob{$g}{'1'} = $prob{$g}{'1'}->add($multiplier);


		$prob{$g}{'2'} = $prob{$g}{'2'}->add($prob{$g-1}{'2'});
		$multiplier = Math::Cephes::Polynomial->new([0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'3'});
		$prob{$g}{'2'} = $prob{$g}{'2'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0.25,-0.5,0.25]);
		$multiplier = $multiplier->mul($prob{$g-1}{'5'});
		$prob{$g}{'2'} = $prob{$g}{'2'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0,0,0.25]);
		$multiplier = $multiplier->mul($prob{$g-1}{'4'});
		$prob{$g}{'2'} = $prob{$g}{'2'}->add($multiplier);

		
		$multiplier = Math::Cephes::Polynomial->new([0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'3'});
		$prob{$g}{'3'} = $prob{$g}{'3'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0,0.5,-0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'5'});
		$prob{$g}{'3'} = $prob{$g}{'3'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0,0.5,-0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'4'});
		$prob{$g}{'3'} = $prob{$g}{'3'}->add($multiplier);

		$multiplier = Math::Cephes::Polynomial->new([0,0,0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'5'});
		$prob{$g}{'4'} = $prob{$g}{'4'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0.5,-1,0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'4'});
		$prob{$g}{'4'} = $prob{$g}{'4'}->add($multiplier);
		

		$multiplier = Math::Cephes::Polynomial->new([0.5,-1,0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'5'});
		$prob{$g}{'5'} = $prob{$g}{'5'}->add($multiplier);
		$multiplier = Math::Cephes::Polynomial->new([0,0,0.5]);
		$multiplier = $multiplier->mul($prob{$g-1}{'4'});
		$prob{$g}{'5'} = $prob{$g}{'5'}->add($multiplier);
	
	}

	$result->{'1'} = $prob{$Grade}{'1'};
	$result->{'2'} = $prob{$Grade}{'2'};
	$result->{'3'} = $prob{$Grade}{'3'};

	$result->{'4'} = Math::Cephes::Polynomial->new([0]);
	$result->{'4'} = $result->{'4'}->add($prob{$Grade}{'4'});
	$result->{'4'} = $result->{'4'}->add($prob{$Grade}{'5'});

	$result->{'5'} = Math::Cephes::Polynomial->new([0]);
	$result->{'5'} = $result->{'5'}->add($prob{$Grade}{'1'});
	$result->{'5'} = $result->{'5'}->add($prob{$Grade}{'2'});
	$result->{'5'} = $result->{'5'}->add($prob{$Grade}{'3'});

	$result->{'6'} = Math::Cephes::Polynomial->new([0]);
	$result->{'6'} = $result->{'6'}->add($prob{$Grade}{'4'});
	$result->{'6'} = $result->{'6'}->add($prob{$Grade}{'5'});
	$result->{'6'} = $result->{'6'}->add($prob{$Grade}{'3'});
	$result->{'6'} = $result->{'6'}->add($prob{$Grade}{'3'});

}

### \\\\\\\\\\\\\\\\\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
sub add {#
	my ($v1,$v2,$r)=@_;
	if (@{$v1} != @{$v2}) {
		print "different size for two vector\n";
		exit(-1);
	}

	@{$r} = map {$v1->[$_]+$v2->[$_]}0..@{$v1}-1;
}

sub numberByVector {#
	my ($number,$v,$r)=@_;

	@{$r} = map {$number*$_} @{$v};

}

sub multiple {#
	my ($v1,$v2,$r)=@_;
	my %result=();
	for (my $i=0;$i<@{$v1} ;$i++) {
		for (my $j=0;$j<@{$v2} ;$j++) {
			$result{$i+$j}+=$v1->[$i]*$v2->[$j];
		}
	}
	
	@{$r} = map {$result{$_}} (sort {$a <=> $b} keys %result);

}

sub derivativePolynomial {#
	my ($polynomial) = @_;
	if (@{$polynomial} == 1) {
		return (0);
	}
	my @result = map {$_*$polynomial->[$_]} 1..@{$polynomial}-1;
	return \@result;
}

 #------------------------------------------------------------------
 #calculate mlod
 #------------------------------------------------------------------

 sub calculate_mlod {#
	
	my ($refData,$markerPair)=@_;
	my (@classfication,%Ocell,$T,$iLOD,$mLOD,$G1)=();
	
	$iLOD=0;
	
	my ($aPgeno,$aMgeno)=split /x/,$refData->{$markerPair->[0]}->{"type"};
	my ($bPgeno,$bMgeno)=split /x/,$refData->{$markerPair->[1]}->{"type"};
	
	return -1 if ((substr($aPgeno,0,1) eq substr($aPgeno,1,1) and substr($aMgeno,0,1) ne substr($aMgeno,1,1) and substr($bPgeno,0,1) ne substr($bPgeno,1,1) and substr($bMgeno,0,1) eq substr($bMgeno,1,1)) or (substr($aPgeno,0,1) ne substr($aPgeno,1,1) and substr($aMgeno,0,1) eq substr($aMgeno,1,1) and substr($bPgeno,0,1) eq substr($bPgeno,1,1) and substr($bMgeno,0,1) ne substr($bMgeno,1,1))) ;
	
	for (my $i=0;$i<@{$markerPair} ;$i++) {
		$classfication[$i]=&classification($refData,$markerPair->[$i]);
	}
	
	foreach my $genoLocusA (sort {$a cmp $b} keys %{$classfication[0]->{$markerPair->[0]}}) {
		
		foreach my $genoLocusB (sort {$a cmp $b} keys %{$classfication[1]->{$markerPair->[1]}}) {
			
			for (my $i=0;$i<@{$refData->{$markerPair->[0]}->{"arr"}} ;$i++) {
				
				$Ocell{join("-",($genoLocusA,$genoLocusB))}++ if ($genoLocusA eq $refData->{$markerPair->[0]}->{"arr"}->[$i] and $genoLocusB eq $refData->{$markerPair->[1]}->{"arr"}->[$i]) ;
			
			}
			
			$T+=$Ocell{join("-",($genoLocusA,$genoLocusB))} if (exists $Ocell{join("-",($genoLocusA,$genoLocusB))}) ;
		}
	}
	
	foreach my $genoLocusA (sort {$a cmp $b} keys %{$classfication[0]->{$markerPair->[0]}}) {
		
		foreach my $genoLocusB (sort {$a cmp $b} keys %{$classfication[1]->{$markerPair->[1]}}) {
			
			next if (not exists $Ocell{join("-",($genoLocusA,$genoLocusB))}) ;
			
			my $R=0;
			my $C=0;
			
			foreach (keys %{$classfication[1]->{$markerPair->[1]}}) {
				
				next if (not exists $Ocell{join("-",($genoLocusA,$_))}) ;
				$R+=$Ocell{join("-",($genoLocusA,$_))};
			
			}
			foreach (keys %{$classfication[0]->{$markerPair->[0]}}) {
				
				next if (not exists $Ocell{join("-",($_,$genoLocusB))}) ;
				$C+=$Ocell{join("-",($_,$genoLocusB))};

			} 
			
			my $E=$R*$C/$T;
			$iLOD+=2*$Ocell{join("-",($genoLocusA,$genoLocusB))}*log($Ocell{join("-",($genoLocusA,$genoLocusB))}/$E);
		}

	}
	
	my $df=((scalar keys %{$classfication[0]->{$markerPair->[0]}})-1)*((scalar keys %{$classfication[1]->{$markerPair->[1]}})-1);
	
	if ($df==1) {
		$G1=$iLOD;
	}else{
		my $e=exp(-$iLOD/(2*($df-1)));
		$G1=((4-$e)*$e-3)*($df-1)+$iLOD;
	}
	
	$mLOD=$G1/(2*log(10));
	return $mLOD;
}

sub classification {#
	
	my ($refData,$marker)=@_;
	
	my %classfication=();
	
	my ($pGeno,$mGeno)=split /x/,$refData->{$marker}->{"type"};
	
	for (my $i=0;$i<length($pGeno) ;$i++) {
		
		for (my $j=0;$j<length($mGeno) ;$j++) {
			
			my ($penoA,$penoB)=sort {$a cmp $b} (substr($pGeno,$i,1),substr($mGeno,$j,1));
			$classfication{$marker}{$penoA.$penoB}=1;
		
		}
	}
	
	return \%classfication; 
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: ltcLod.pl (计算二倍体CP,F2,RILs的重组率和mLOD)
Version: $version
Contact: Li Tiancheng for (CP) (litc\@biomarker.com.cn)
		 Huang Long for (F2,BC1,RILs) (huangl\@biomarker.com.cn)
Description:
	CP contains:nnxnp,abxcd,efxeg,hkxhk,lmxll
	F2,RILs:aaxbb，offspring have five type a，b，c，d，h
	BC1：aaxbb，offspring have two type a，h
	
	Options:
		-i	<file>	input file
		-o	<file>	output file
		-t	<str>	type of population  RILs,CP,F2
				RILs need give the generation like Ri10 default Ri10
		-p <int>	only calcualte the first n markers and the other markers, default off(-t=CP).

		-h	Help

USAGE
	print $usage;
	exit;
}
