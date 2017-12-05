#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use lib "$Bin/exlib/lib/perl5/";
use Math::Cephes::Matrix;
use Math::Cephes::Polynomial;
use Graph;
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fKey,$dOut,$filter,$FI_threshold,$mapFunction,$imput,$OnlyFirst);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				"p:i"=>\$OnlyFirst,
				
				"FILTER"=>\$filter,
				"t:s"=>\$FI_threshold,

				"MAPFUNCTION:s"=>\$mapFunction,
				
				"IMPUT"=>\$imput,
				) or &USAGE;
&USAGE unless ($fIn and $fKey);

$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

$OnlyFirst = defined $OnlyFirst ? $OnlyFirst: 0;
$FI_threshold||=1;
$mapFunction||='Kosambi';
my %mapFunction=(
	"Kosambi"=>\&Kosambi,
	"Haldane"=>\&Haldane,
	);
my %inverseMapFunction=(
	"Kosambi"=>\&inverseKosambi,
	"Haldane"=>\&inverseHaldane,
	);
# ------------------------------------------------------------------
# Load CP Data
# ------------------------------------------------------------------
my %Data;
my $rGraph = Graph->new();

if (`grep '{..}' -P $fIn`) {
	print "file type: loc with linkage phase\n";
	&LoadCP_fromLOC($fIn,\%Data);
}else{
	print "file type: loc with no info of linkage phase or genotype file\n";
	&LoadCP_fromBMK($fIn,\%Data);
}

# ------------------------------------------------------------------
# Calculate and Output
# ------------------------------------------------------------------
open (OUT,">","$dOut/$fKey.pwd.detail") or die $!;
open (PWD,">","$dOut/$fKey.pwd") or die $!;
open (LOG,">","$dOut/$fKey.log") or die $! if ($filter or ($filter and $imput)) ;
open (FILTER,">","$dOut/$fKey.filtered.pwd") or die $! if ($filter) ;
print LOG ">FILTER\n" if ($filter) ;

my @Marker=sort {$Data{$a}{"order"} <=> $Data{$b}{"order"}} keys %Data;
my $EndMarkerCurr=$OnlyFirst ? $OnlyFirst : $#Marker;


print PWD ";",GetTime(),"\n";

my @EndCurrLM = grep {$Data{$Marker[$_]}{'type'} eq 'lmxll'} 0..$EndMarkerCurr ;
my @EndCurrNP = grep {$Data{$Marker[$_]}{'type'} eq 'nnxnp'} 0..$EndMarkerCurr ;
my @AllLM = grep {$Data{$Marker[$_]}{'type'} eq 'lmxll'} 0..$#Marker ;
my @AllNP = grep {$Data{$Marker[$_]}{'type'} eq 'nnxnp'} 0..$#Marker ;

my $allPairs = @Marker*(@Marker-1)/2;
my $notInformativePairs = @EndCurrLM*@AllNP + @EndCurrNP*@AllLM - @EndCurrLM*@EndCurrNP;
my $informativePairs =$allPairs - $notInformativePairs;

print PWD ";$informativePairs pairs out of $allPairs, $notInformativePairs not informative\n";
print PWD ";name = ",basename($fKey),"\n";

for (my $i=0;$i<$EndMarkerCurr ;$i++) {
	for (my $j=$i+1;$j<@Marker ;$j++) {

#		print ">$Marker[$i]\t$Marker[$j]\n";

		my $type1=$Data{$Marker[$i]}{"type"};
		my $type2=$Data{$Marker[$j]}{"type"};

		my $phase1 = $Data{$Marker[$i]}{"phase"};
		my $phase2 = $Data{$Marker[$j]}{"phase"};

		my ($RecCC,$LodCC,$FiCC,$RecCR,$LodCR,$FiCR,$RecRC,$LodRC,$FiRC,$RecRR,$LodRR,$FiRR,$result) = (defined $phase1 and defined $phase2)? &calculatePairWiseData($type1,$type2,\@{$Data{$Marker[$i]}{"arr"}},\@{$Data{$Marker[$j]}{"arr"}},$phase1,$phase2):&calculatePairWiseData($type1,$type2,\@{$Data{$Marker[$i]}{"arr"}},\@{$Data{$Marker[$j]}{"arr"}});
#		my $mLOD=calculatemLOD(\%Data,[$Marker[$i],$Marker[$j]]);
		my $mLOD = calculatemLOD($type1,$type2,\@{$Data{$Marker[$i]}{"arr"}},\@{$Data{$Marker[$j]}{"arr"}});
		if ($mLOD == -2) {
			print "$Marker[$i] or $Marker[$j] is exteame segregation distortion,calculate mLOD failed\n";
			$mLOD = 0.001;
#			next;
		}

#		print Dumper %{$result};
		my $valid = 0;
		if (defined $result->{'CC'}{'valid'}) {
			$valid = $result->{'CC'}{'valid'};
		}elsif(defined $result->{'CR'}{'valid'}){
			$valid = $result->{'CR'}{'valid'};
		}elsif(defined $result->{'RC'}{'valid'}){
			$valid = $result->{'RC'}{'valid'};
		}else{
			$valid = $result->{'RR'}{'valid'};
		}

		print OUT $Marker[$i],"\t",$Marker[$j],"\t",$type1,"\t",$type2,"\t",$valid,"\t",join("\t",map {sprintf("%.16e",$_)} ($RecCC,$LodCC,$FiCC,$RecCR,$LodCR,$FiCR,$RecRC,$LodRC,$FiRC,$RecRR,$LodRR,$FiRR,$mLOD)),"\n";
		next if ($mLOD == -1) ;
			
		my %r=();
		foreach my $lp (keys %{$result}) {
			$r{$result->{$lp}{'r'}}{$lp}=1;
		}
		my ($minR) = sort {$a <=> $b} keys %r;
		my @lp = keys %{$r{$minR}};

		if ($minR <= 0.49) {
			print PWD join("\t",($Marker[$i],$Marker[$j],(map {sprintf("%.16e",$_)}($minR,$mLOD,$result->{$lp[0]}{'lod'},$result->{$lp[0]}{'fi'})),$result->{$lp[0]}{'valid'},join(",",@lp))),"\n";
		}else{
			print PWD join("\t",($Marker[$i],$Marker[$j],(map {sprintf("%.16e",$_)}(0.49999,0.00001,0,0)),$valid,join(",",@lp)))," ;  original estimate: ", join("\t",map {sprintf("%.16e",$_)}($minR,$mLOD,$result->{$lp[0]}{'lod'},$result->{$lp[0]}{'fi'})),"\n";
		}

		if ($filter) {
			if ($result->{$lp[0]}{'fi'} < $FI_threshold) {
				if ($minR < 0.49) {
					print LOG join("\t",($Marker[$i],$Marker[$j],(map {sprintf("%.16e",$_)}($minR,$mLOD,$result->{$lp[0]}{'lod'},$result->{$lp[0]}{'fi'})),$result->{$lp[0]}{'valid'},join(",",@lp))),"\n";
				}else{
					print LOG join("\t",($Marker[$i],$Marker[$j],(map {sprintf("%.16e",$_)}(0.49999,0.00001,0,0)),$result->{$lp[0]}{'valid'},join(",",@lp)))," ;  original estimate: ", join("\t",map {sprintf("%.16e",$_)}($minR,$mLOD,$result->{$lp[0]}{'lod'},$result->{$lp[0]}{'fi'})),"\n";
				} 
				
			}else{
				if ($minR < 0.49) {
					print FILTER join("\t",($Marker[$i],$Marker[$j],(map {sprintf("%.16e",$_)}($minR,$mLOD,$result->{$lp[0]}{'lod'},$result->{$lp[0]}{'fi'})),$result->{$lp[0]}{'valid'},join(",",@lp))),"\n";
				}else{
					print FILTER join("\t",($Marker[$i],$Marker[$j],(map {sprintf("%.16e",$_)}(0.49999,0.00001,0,0)),$result->{$lp[0]}{'valid'},join(",",@lp)))," ;  original estimate: ", join("\t",map {sprintf("%.16e",$_)}($minR,$mLOD,$result->{$lp[0]}{'lod'},$result->{$lp[0]}{'fi'})),"\n";
				} 
	#				print join("\t",($Marker[$i],$Marker[$j])),"\n";
	
				my $weight = &{$mapFunction{$mapFunction}}($minR);	
				$rGraph->add_weighted_edge($Marker[$i],$Marker[$j],$weight) if ($imput) ;
			}
		}
	}
}

my $graphType = $rGraph->is_directed();

#print $graphType;die;


if ($imput) {
	
	my $APSP = $rGraph->APSP_Floyd_Warshall();

#	print Dumper @{$APSP};die;

	print "      ";
	foreach my $v ( $APSP->vertices ) { printf "%-9s ", "$v" } print "\n";
	foreach my $u ( $APSP->vertices ) {
		print "$u: ";
		foreach my $v ( $APSP->vertices ) {
			my $w = $APSP->get_attribute("weight", $u, $v);
			if (defined $w) {
				my $p = $APSP->get_attribute("path",   $u, $v);
				printf "(%-5s)=%d ", "@$p", $w
			} else {
				printf "%-9s ", "-"
			}
		}
		print "\n"
	}

}

close (OUT) ;
close (PWD);
close (LOG) if ($filter or ($filter and $imput)) ;
close (FILTER) if ($filter) ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub determine_lp_for_loci_pair {#
	my ($phase1,$phase2) = @_;

	my ($p1,$m1) = $phase1=~/(.)(.)/;
	my ($p2,$m2) = $phase2=~/(.)(.)/;
	
	my @plp = ();
	my @mlp = ();

	if ($p1 eq '-' or $p2 eq '-') {
		@plp = ('C','R');
	}elsif($p1 ne $p2){
		@plp = ('R');
	}else{
		@plp = ('C');
	}

	if ($m1 eq '-' or $m2 eq '-') {
		@mlp = ('C','R');
	}elsif($m1 ne $m2){
		@mlp = ('R');
	}else{
		@mlp = ('C');
	}

	my @possibleLP = ();
	foreach my $i (@plp) {
		foreach my $j (@mlp) {
			push @possibleLP,$i.$j;	
		}
	}

#	print $phase1,"\t",$phase2,"\t",join("\t",@possibleLP),"\n";

	my %result = ();
	foreach(@possibleLP){
		if ($_ eq 'CC') {
			$result{$_} = "01X01";
		}elsif($_ eq 'CR'){
			$result{$_} = "01X10";
		}elsif($_ eq 'RC'){
			$result{$_} = "10X01";
		}else{
			$result{$_} = "10X10";
		}			
	} 

	return %result;
}

### \\\\\\\\calculate recombination frequency,linkage LOD score,Fisher information in each linkagr phase configugation \\\\\\\\\\

sub calculatePairWiseData {# for CP°°population 
	
	my ($type1,$type2,$progenyGeno1,$progenyGeno2,$phase1,$phase2)=@_;

	if (@{$progenyGeno1} != @{$progenyGeno2}) {
		print "different population size for two loci,please check!\n";
		exit(-1);
	}

	my $nind = @{$progenyGeno1};
	my %linkagePhase = ();
	if (defined $phase1 and defined $phase2) {
		%linkagePhase = determine_lp_for_loci_pair($phase1,$phase2);
	}else{
		%linkagePhase=(
		"CC" => "01X01",
		"CR" => "01X10",
		"RC" => "10X01",
		"RR" => "10X10",
		);
	}

	my @type1=$type1=~/(\w+)x(\w+)/;
	my @type2=$type2=~/(\w+)x(\w+)/;

	my @p1=$type1[0]=~/(.)(.)/;
	my @m1=$type1[1]=~/(.)(.)/;

	my @p2=$type2[0]=~/(.)(.)/;
	my @m2=$type2[1]=~/(.)(.)/;


	my %result=();
	foreach my $lp (sort {$a cmp $b} keys %linkagePhase) {

		my ($plp1,$plp2,$mlp1,$mlp2)=$linkagePhase{$lp}=~/(.)(.)X(.)(.)/;

		##  gametes produced by parents
		my %Gametes=();
		foreach my $pallel1 (@p1) {
			foreach my $pallel2 (($p2[$plp1],$p2[$plp2])) { 
				push @{$Gametes{'P'}},join('',($pallel1,$pallel2));
			}
		}

		foreach my $mallel1 (@m1) {
			foreach my $mallel2 (($m2[$mlp1],$m2[$mlp2])) {
				push @{$Gametes{'M'}},join('',($mallel1,$mallel2));
			}
		}

		## assign gametes into classification {recombination,non_recombination}
		my %gameteType=();
		$gameteType{'P'}{$Gametes{'P'}->[0]}{'nr'}++;
		$gameteType{'P'}{$Gametes{'P'}->[3]}{'nr'}++;
		$gameteType{'P'}{$Gametes{'P'}->[1]}{'r'}++;
		$gameteType{'P'}{$Gametes{'P'}->[2]}{'r'}++;


		$gameteType{'M'}{$Gametes{'M'}->[0]}{'nr'}++;
		$gameteType{'M'}{$Gametes{'M'}->[3]}{'nr'}++;
		$gameteType{'M'}{$Gametes{'M'}->[1]}{'r'}++;
		$gameteType{'M'}{$Gametes{'M'}->[2]}{'r'}++;

#		print "> linkage phase : $lp \n";
#
#		print "============= gamete type ============\n";
#		print Dumper %gameteType;

		## initialize conditional probability of gamete's calssification 
		my %unit=(
			'r' => [0,0.5],
			'nr' => [0.5,-0.5],
			);
	
		## calculate probability of each gamete
		my %probGametes=();
		foreach my $parent (('P','M')) {

			foreach my $Gametes (keys %{$gameteType{$parent}}) {
			
				$gameteType{$parent}{$Gametes}{'nr'}||=0;
				$gameteType{$parent}{$Gametes}{'r'}||=0;

				my (@arr1,@arr2,@sum)=();

				&numberByVector($gameteType{$parent}{$Gametes}{'nr'},$unit{'nr'},\@arr1);
				&numberByVector($gameteType{$parent}{$Gametes}{'r'},$unit{'r'},\@arr2);
				&add(\@arr1,\@arr2,\@sum);

				$probGametes{$parent}{$Gametes}=\@sum;
				
			}
		}
#		print "============= prob gamete ============\n";
#		print Dumper %probGametes;

		## calculate probability for all possible zygote
		my %zygot_com=();
		foreach my $pGametes (keys %{$gameteType{'P'}}) {
			foreach my $mGametes (keys %{$gameteType{'M'}}) {
				my @multiple=();
				&multiple($probGametes{'P'}{$pGametes},$probGametes{'M'}{$mGametes},\@multiple);
				$zygot_com{$pGametes."X".$mGametes}=\@multiple;

			}
		}

#		print "============ zygote gamete =============\n";
#		print Dumper %zygot_com;

		#calculate probability for all possible phenotype 

		my %locus_1_classification = map {$_,1} (join("",sort {$a cmp $b} ($p1[0],$m1[0])),join("",sort {$a cmp $b} ($p1[0],$m1[1])),join("",sort {$a cmp $b} ($p1[1],$m1[0])),join("",sort {$a cmp $b} ($p1[1],$m1[1])));
		my %locus_2_classification = map {$_,1} (join("",sort {$a cmp $b} ($p2[0],$m2[0])),join("",sort {$a cmp $b} ($p2[0],$m2[1])),join("",sort {$a cmp $b} ($p2[1],$m2[0])),join("",sort {$a cmp $b} ($p2[1],$m2[1])));


		my %combination = ();

		foreach my $pheno_1 (sort {$a cmp $b} keys %locus_1_classification) {
			foreach my $pheno_2 (sort {$a cmp $b} keys %locus_2_classification) {

				my @allel1 = $pheno_1 =~/(.)(.)/;
				my @allel2 = $pheno_2 =~/(.)(.)/;
				
				for (my $i=0;$i<@allel1 ;$i++) {
					for (my $j=0;$j<@allel2 ;$j++) {
						
						my $zygot_com = join("X",(join("",($allel1[$i],$allel2[$j])),join("",($allel1[1-$i],$allel2[1-$j]))));
						$combination{$pheno_1.$pheno_2}{$zygot_com}=$zygot_com{$zygot_com} if (exists $zygot_com{$zygot_com}) ;
					
					}
				}

			}
		}

#		print "======== combination =========\n";
#		print Dumper %combination;

		my %probType=();
		foreach my $gtype (keys %combination) {
			my @sum = map {0}0..2;
			foreach my $gametesCom (keys %{$combination{$gtype}}) {
				&add($combination{$gtype}{$gametesCom},\@sum,\@sum);
			}

			if ($sum[0]!=0 and $sum[1] !=0 and $sum[0] == 0-$sum[1] and $sum[2] == 0) {

				$probType{$gtype}{'nr'}=1;
			}elsif ($sum[0] ==0 and $sum[1] !=0 and $sum[2] == 0) {

				$probType{$gtype}{'r'}=1;
			}elsif ($sum[0] ==0 and $sum[1] ==0 and $sum[2] != 0) {

				$probType{$gtype}{'r'}=2;
			}elsif ($sum[0] !=0 and $sum[1] !=0 and $sum[2] != 0 and $sum[0] == $sum[2] and $sum[0] == 0-$sum[1]/2) {

				$probType{$gtype}{'nr'}=2;
			}elsif ($sum[0] ==0 and $sum[1] !=0 and $sum[2] != 0 and $sum[1] == 0-$sum[2]) {

				$probType{$gtype}{'r'}=1;
				$probType{$gtype}{'nr'}=1;
			}elsif($sum[0] !=0 and $sum[1] ==0 and $sum[2] == 0){
				next;
			}else{
				$probType{$gtype}{'unsolve'}=\@sum;

			}

		}

#		print "======== prob type =======\n";
#		print Dumper %probType;

		## confirm the form of likelihood function
		my %count=();
		my $valid_observation = 0;
		for (my $i=0;$i<@{$progenyGeno1} ;$i++) {

			next unless ($progenyGeno1->[$i] ne '--' and $progenyGeno2->[$i] ne '--') ;
			$valid_observation++;

			my $gtype = join('',($progenyGeno1->[$i],$progenyGeno2->[$i]));

			if (defined $probType{$gtype}) {

				foreach (keys %{$probType{$gtype}}) {

					if ($_ eq 'r' or $_ eq 'nr') {

						$count{$_}+=$probType{$gtype}{$_};

					}else{
						$count{$_}{join("|",@{$probType{$gtype}{$_}})}+=1;
					}
					
					
				}
			}

		}
#		print "======= count ==========\n";
#		print Dumper %count;

		$result{$lp}{'valid'} = sprintf("%.4f",$valid_observation/$nind);
		## maximize the likelihood function and get the estimation of r 
		if (not exists $count{'unsolve'}) { # the estimation can be calculated by an explicit formula:n1/(n1+n2)
			
			$count{'r'}||=0;
			$count{'nr'}||=0;
			my $valid=$count{'r'}+$count{'nr'};
			
			if ($valid!=0) {
				$result{$lp}{'r'} = $count{'r'}/$valid;
				if ($result{$lp}{'r'} > 0.5) {
					
					$result{$lp}{'r'} = 0.4999;
					
					$result{$lp}{'fi'} = 0;
				}
				if ($count{'r'} == 0) {

					$result{$lp}{'fi'} = 'INF';
					
					$result{$lp}{'lod'} = $count{'nr'}*log10(1-$result{$lp}{'r'}) + $valid*log10(2);

				}elsif($count{'nr'} == 0){

					$result{$lp}{'lod'} = $count{'r'}*log10($result{$lp}{'r'}) + $valid*log10(2);
					$result{$lp}{'fi'} = 'INF';

				}else{

					$result{$lp}{'lod'} = $count{'r'}*log10($result{$lp}{'r'}) + $count{'nr'}*log10(1-$result{$lp}{'r'}) + $valid*log10(2);

					$result{$lp}{'fi'} = $valid/($result{$lp}{'r'}*(1-$result{$lp}{'r'})*$nind);
				}
			}else {

				$result{$lp}{'r'} = 0;
				$result{$lp}{'lod'} = 0;
				$result{$lp}{'fi'} = 0;
			}
		}else{

			my %factor = ();
			my $valid = 0;

			foreach my $flag (keys %count) {

				if ($flag eq 'r' or $flag eq 'nr') {

					$factor{$flag}{'poly'} = Math::Cephes::Polynomial->new([0,1]) if ($flag eq 'r') ;

					$factor{$flag}{'poly'} = Math::Cephes::Polynomial->new([1,-1]) if ($flag eq 'nr') ;

					$factor{$flag}{'count'} = $count{$flag};

					$valid += $count{$flag};

				}else{

					foreach my $joinCoef (keys %{$count{$flag}}) {

						$factor{$joinCoef}{'poly'} = Math::Cephes::Polynomial->new([map {$_/0.25} (split /\|/,$joinCoef)]);

						$factor{$joinCoef}{'count'} = $count{$flag}{$joinCoef};

						$valid += $count{$flag}{$joinCoef};
					}
				}
			}

			my $numerator = Math::Cephes::Polynomial->new([0]);

			foreach my $factor (keys %factor) {
				
				my @temp;

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
			for (my $i=0;$i<$numerator->{'n'} ;$i++) {
				next if ($roots->[$i]->i != 0) ;

				push @realRoots,$roots->[$i]->r;

				$result{$lp}{'r'} = $roots->[$i]->r;
			}

			if (@realRoots!=0) {

				($result{$lp}{'r'}) = sort {$a <=> $b} @realRoots;
			
			}else{
				# debugged by macx 2012-6-28
				if (not exists $factor{'r'}) {
					$result{$lp}{'r'}||=0;
				}else{
					$result{$lp}{'r'}||=0.49999999; 
				}
			}
		 
			if ($result{$lp}{'r'} == 0) {

				$result{$lp}{'fi'} = 'INF';
			}else{

				$result{$lp}{'fi'} = ($valid/$result{$lp}{'r'} + $valid/(1-$result{$lp}{'r'})-4*$valid + $valid*(4*$result{$lp}{'r'}-2)**2/(1-2*$result{$lp}{'r'}+2*$result{$lp}{'r'}**2))/$nind;
				
			}
			
		
			$result{$lp}{'lod'} = 0;

			foreach my $factor (keys %factor) {

#				next if ($factor{$factor}{'poly'}->eval($result{$lp}{'r'} == 0)) ;

				if ($factor{$factor}{'poly'}->eval($result{$lp}{'r'}) == 0) {  ## debug by macx 2012-3-27
					
					$result{$lp}{'lod'} += $factor{$factor}{'count'}*log10(2);

					next;

				}

				$result{$lp}{'lod'} += $factor{$factor}{'count'}*log10($factor{$factor}{'poly'}->eval($result{$lp}{'r'})/$factor{$factor}{'poly'}->eval(0.5));

			}
			
		}

#		print "$result{$lp}{'r'}\t$result{$lp}{'lod'}\n";

	}

#	print join("\t",($result{'CC'}{'r'},$result{'CC'}{'lod'},$result{'CC'}{'fi'},$result{'CR'}{'r'},$result{'CR'}{'lod'},$result{'CR'}{'fi'},$result{'RC'}{'r'},$result{'RC'}{'lod'},$result{'RC'}{'fi'},$result{'RR'}{'r'},$result{'RR'}{'lod'},$result{'RR'}{'fi'})),"\n";

	$result{'CC'}{'r'}=0.49999 unless (defined $result{'CC'}{'r'}) ; 
	$result{'CC'}{'lod'}=0.00001 unless (defined $result{'CC'}{'lod'}) ;  
	$result{'CC'}{'fi'}=0.00001 unless (defined $result{'CC'}{'fi'}) ;  
	$result{'CR'}{'r'}=0.49999 unless (defined $result{'CR'}{'r'}) ;  
	$result{'CR'}{'lod'}=0.00001 unless (defined $result{'CR'}{'lod'}) ; 
	$result{'CR'}{'fi'}=0.00001 unless (defined $result{'CR'}{'fi'}) ;
	$result{'RC'}{'r'}=0.49999 unless (defined $result{'RC'}{'r'}) ; 
	$result{'RC'}{'lod'}=0.00001 unless (defined $result{'RC'}{'lod'}) ;
	$result{'RC'}{'fi'}=0.00001 unless (defined $result{'RC'}{'fi'}) ;
	$result{'RR'}{'r'}=0.49999 unless (defined $result{'RR'}{'r'}) ; 
	$result{'RR'}{'lod'}=0.00001 unless (defined $result{'RR'}{'lod'}) ; 
	$result{'RR'}{'fi'}=0.00001 unless (defined $result{'RR'}{'fi'}) ;

	return ($result{'CC'}{'r'},$result{'CC'}{'lod'},$result{'CC'}{'fi'},$result{'CR'}{'r'},$result{'CR'}{'lod'},$result{'CR'}{'fi'},$result{'RC'}{'r'},$result{'RC'}{'lod'},$result{'RC'}{'fi'},$result{'RR'}{'r'},$result{'RR'}{'lod'},$result{'RR'}{'fi'},\%result);

}

### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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

### \\\\\\\\\\\\\  º∆À„mLOD \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

################ Method introduced from joinmap official website : http://www.kyazma.nl/index.php/mc.JoinMap/sc.FAQ/#Q9
#	For each pair of loci a contingency table is produced of the genotypes. The dimensions of the table depend on the population type and the segregation types of the loci (unknown genotypes are ignored). From this table the G statistic (-2 times the logarithm of the likelihood ratio for the Poisson distribution, see e.g. Fienberg, 1979, The analysis of cross-classified categorical data, MIT Press) is calculated to test for independence; the expected number (E) in each cell is calculated from the R-total (R), the column-total (C) and the grand-total (T): 
#
#	E = R*C/T .
#
#
#	The G statistic then is a summation (SUM) over all cells (O is the observed number, ln() is the natural logarithm): 
#
#	G = 2 * SUM [ O*ln(O/E) ] .
#
#
#	The G statistic (Gd) has an approximate chi-square distribution with the number of Rs in the table minus 1 multiplied by the number of columns minus 1 as the degrees of freedom (d). When the loci have different numbers of genotypes in their segregation, this would present a problem, because this number affects the degrees of freedom in the G test and in the testing of linkage one would need to take account of the degrees of freedom. In order to remove this problem and to ensure the comparability of data coming from different population or segregation types the G statistic with d degrees of freedom, Gd, is transformed approximately to a G statistic, G1, that would have been obtained if there was just a single degree of freedom (as if in a backcross). This approximate transformation is an empirically determined formula (exp() is the exponential function): 
#
#	e = exp( -Gd/(2*(d-1)) ) ,
#
#	G1 = ((4-e)*e - 3)*(d-1) + Gd . 
#
#
#	Because in genetics one is used to LOD scores, which are likelihood ratio statistics using the 10-base logarithm instead of the natural logarithm multiplied by -2, the modified LOD score (mLOD) is simply derived from G1: 
#
#	mLOD = G1 / (2*ln(10)) .
#
#	It can be shown that for the case of two loci each segregating in two genotypes (e.g. population types BC1 or DH1) when there is no segregation distortion this mLOD equals the 'normal' LOD score. The modified LOD score is not sensitive to segregation distortion, in contrast to the 'normal' LOD score.
#
#####################################################

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
		return -2;
	}

	if (scalar keys %R == 2 and scalar keys %C == 2) {
		$G1=$G;
	}else{
		my $df = (scalar(keys %R) -1) * (scalar(keys %C) -1);
		$e = exp(-$G /(2 * ($df -1)));
		$G1=((4-$e)*$e-3)*($df-1)+$G;
	}
	my $mlod=$G1/(2*log(10));

	return $mlod;	
}

## \\\\\\\\\\\\\\\\\\\\\\\\\\\ Floyd-wallshall algorithm \\\\\\\\\\\\\\\\\

# APSP_Floyd_Warshall
#
#       $APSP = $G->APSP_Floyd_Warshall
#
#       Returns the All-pairs Shortest Paths graph of the graph $G
#       computed using the Floyd-Warshall algorithm and the attribute
#       'weight' on the edges.
#       The returned graph has an edge for each shortest path.
#       An edge has attributes "weight" and "path"; for the length of
#       the shortest path and for the path (an anonymous list) itself.


sub APSP_Floyd_Warshall {
    my $G = shift;
    my @V = $G->vertices;
    my @E = $G->edges;
    my (%V2I, @I2V);
    my (@P, @W);
    # Compute the vertex <-> index mappings.
    @V2I{ @V     } = 0..$#V;
    @I2V[ 0..$#V ] = @V;
    # Initialize the predecessor matrix @P and the weight matrix @W.
    # (The graph is converted into adjacency-matrix representation.)
    # (The matrix is a list of lists.)
    foreach my $i ( 0..$#V ) { $W[ $i ][ $i ] = 0 }
    while ( my ($u, $v) = splice(@E, 0, 2) ) {
        my ( $ui, $vi ) = ( $V2I{ $u }, $V2I{ $v } );
        $P[ $ui ][ $vi ] = $ui unless $ui == $vi;
        $W[ $ui ][ $vi ] = $G->get_attribute( 'weight', $u, $v );
    }
    # Do the O(N**3) loop.
    for ( my $k = 0; $k < @V; $k++ ) {
        my (@nP, @nW); # new @P, new @W
        for ( my $i = 0; $i < @V; $i++ ) {
            for ( my $j = 0; $j < @V; $j++ ) {
                my $w_ij    = $W[ $i ][ $j ];
                my $w_ik_kj = $W[ $i ][ $k ] + $W[ $k ][ $j ]
                    if defined $W[ $i ][ $k ] and
                       defined $W[ $k ][ $j ];
                # Choose the minimum of w_ij and w_ik_kj.
                if ( defined $w_ij ) {
                    if ( defined $w_ik_kj ) {
                        if ( $w_ij <= $w_ik_kj ) {
                          $nP[ $i ][ $j ] = $P[ $i ][ $j ];                          
						  $nW[ $i ][ $j ] = $w_ij;
                        } else {
                          $nP[ $i ][ $j ] = $P[ $k ][ $j ];
                          $nW[ $i ][ $j ] = $w_ik_kj;
                        }
                    } else {
                        $nP[ $i ][ $j ] = $P[ $i ][ $j ];
                        $nW[ $i ][ $j ] = $w_ij;
                    }
                } elsif ( defined $w_ik_kj ) {
                    $nP[ $i ][ $j ] = $P[ $k ][ $j ];
                    $nW[ $i ][ $j ] = $w_ik_kj;
                }
            }
        }

        @P = @nP; @W = @nW; # Update the predecessors and weights.
    }
    # Now construct the APSP graph.
    my $APSP = (ref $G)->new;
    $APSP->directed( $G->directed ); # Copy the directedness.
    # Convert the adjacency-matrix representation
    # into a Graph (adjacency-list representation).
    for ( my $i = 0; $i < @V; $i++ ) {
        my $iv = $I2V[ $i ];
        for ( my $j = 0; $j < @V; $j++ ) {
            if ( $i == $j ) {
                $APSP->add_weighted_edge( $iv, 0, $iv );
                $APSP->set_attribute("path", $iv, $iv, [ $iv ]);
                next;
            }
            next unless defined $W[ $i ][ $j ];
            my $jv = $I2V[ $j ];
            $APSP->add_weighted_edge( $iv, $W[ $i ][ $j ], $jv );
            my @path = ( $jv );
            if ( $P[ $i ][ $j ] != $i ) {
                my $k = $P[ $i ][ $j ];  # Walk back the path.               
				while ( $k != $i ) {
                    push @path, $I2V[ $k ];
                    $k = $P[ $i ][ $k ]; # Keep walking.
                }
            }
            $APSP->set_attribute( "path", $iv, $jv,[ $iv, reverse @path ] );
        }
    }
    return $APSP;
}

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

sub LoadCP_fromLOC {#
	my ($file,$refData)=@_;
	open (IN,"<",$file) or die $!;
	my $curr=1;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^;/ or /^nloc/ or /^nind/ or /^name/ or /^popt/) ;
		#1	marker001	<abxcd>	{00}	"(ac,ad,bc,bd)"	ad	bc	bd	bd	ac	bd	bc	ad	ad	ad	bc	ac	bd	bc	ad	bc	ad	bd	bd	ac	ad	bc	bc	ac	bd	ac	bd	bc	bd	bd	bd	ac	bd	bc	bd	bd	bc	bd	bc	ac	ad	ac	ac	bc	ad	ac	ac	bd	ad	bc	ac	ad	ad	bc	bd	bc	bd	ad	ac	ac	bc	ac	bc	bc	ac	bd	ad	bc	bc	ac	ac	ad	bd	bd	ad	ac	ac	ac	ad	bd	bc	bc	ad	bd	ac	bd	ad	bc	bd	bd	ad	bc	bd	bc	bd	bd	ad	bc	bc	ad
		my ($ID,$type,$phase,@indi)=split /\s+/,$_;
		$type=~s/[<>]//g;
		$phase=~s/\{|\}//g;
		$$refData{$ID}{"type"}=$type;
		$$refData{$ID}{"phase"}=$phase;
		@{$$refData{$ID}{"arr"}}=@indi;
		$$refData{$ID}{"order"}=$curr++;
	}
	close (IN) ;
}

sub LoadCP_fromBMK {#
	my ($file,$refData)=@_;
	open (IN,"<",$file) or die $!;
	my $curr=1;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^;/ or /[Tt]ype/ or /^nloc/ or /^nind/ or /^name/ or /^popt/) ;
		#Marker10000     nnxnp   np -- -- nn np np nn nn -- np nn nn np nn nn nn nn nn np np nn nn nn -- nn nn np np np nn nn nn nn np -- nn nn np nn -- -- nn np nn nn nn
		my ($ID,$type,@indi)=split /\s+/,$_;
		$type=~s/[<>]//g;
		$$refData{$ID}{"type"}=$type;
		$$refData{$ID}{'phase'}=undef;
		@{$$refData{$ID}{"arr"}}=@indi;
		$$refData{$ID}{"order"}=$curr++;
	}
	close (IN) ;
}

sub log10 {#
	my $l=shift;
	return log($l)/log(10);
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

sub inverseHaldane {#
	my $r=shift;
	my $result=(1/2)*(1-exp((-1)*$r/50));
	return $result;
}

sub inverseKosambi {#
	my $r=shift;
	my $result=(1/2)*(exp($r/25)-1)/(exp($r/25)+1);
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
Contact: Ma Chouxian <macx\@biomarker.com.cn> 
Description:

	This program is used to calcualte the recombination frequency and LOD score in linkage analysis for CP population.All output is equal to joinMap\'result.
	
	The program also give the result of mLOD calculated according to joinMap\'s method.

	Meanwhile, it shows the estimation quality of r using a statistic : Fisher information.
	

Usage:
  Options:
  -i <file>  Input loc format file or bmk genetype format file, forced
  -k	<str>	Key of output file,forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -p <int>   only calcualte the first n markers and the other markers, default off.
  -h         Help

USAGE
	print $usage;
	exit;
}
