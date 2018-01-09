#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fPWD,$fGenotype,$fKey,$dOut,$rThreshold,$lodThreshold,$log);
GetOptions(
				"help|?" =>\&USAGE,
				"p:s"=>\$fPWD,
				"g:s"=>\$fGenotype,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				"r:s"=>\$rThreshold,
				"LOD:s"=>\$lodThreshold,
				) or &USAGE;
&USAGE unless ($fPWD and $fGenotype and $fKey);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$fKey="$dOut/$fKey";
$rThreshold||=0.4;
$lodThreshold||=1;
#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my (%pwd,%G,$nind) = ();

my $nloc = 0;   ## number of loci with linkage phase
#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------

my @pwdMarker = loadPWD(\%pwd,$fPWD);

my $phase=loadGtype(\%G,$fGenotype);
if ($phase == 1) {
	open In,$fGenotype;
	open Out,">$fKey.loc";
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		print Out $_,"\n";
	}
	close Out;
	close In;
}else{
	&checkDataSet(\@pwdMarker,[keys %G]);

	#-------------------------------------------------------------------
	# determine linkage phase
	#-------------------------------------------------------------------

	open ($log,">$fKey.dlp.log") or die $!;

	my @determined = anchorFirstPair(\%pwd,\%G);
	$G{$determined[0]}{'lp'} = singleDetermine($G{$determined[0]}{'type'});
	my @temp = pairDetermine($G{$determined[0]}{'lp'},$G{$determined[1]}{'type'},$pwd{$determined[0]}{$determined[1]}{'lp'});
	#if (@temp == 1) {
		$G{$determined[1]}{'lp'} = $temp[0];
		$nloc = 2;
	#}else{

	#	die "unsuccessful in linkage phase determination\n";
	#}

	##\\\\\\\\\\\\\\\\\\\\\\  DEBUG \\\\\\\\\\\\\\\\\\\\\\\\\\\

	print $log "-----------------------------------\n";

	print $log "start with : @determined\t$pwd{$determined[0]}{$determined[1]}{'lp'} (rec = $pwd{$determined[0]}{$determined[1]}{'r'}  lod = $pwd{$determined[0]}{$determined[1]}{'lod'})\n";

	print $log "\n--- LOCI --- TYPE --- LP -----\n";

	foreach my $dLoci (@determined) {
		print $log $dLoci,"\t",$G{$dLoci}{'type'},"\t{",$G{$dLoci}{'lp'},"}\n";
	}

	print $log "\n\n";

	#die;

	## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	my @undetermined = su([keys %G],\@determined);

	my %undetermined = map {$_,1} @undetermined;
	my %buffer = map {$_,1} @undetermined;
	my $loop_control_flag = scalar keys %buffer;
	my $failure_flag = scalar keys %buffer;



	while (1) {

		if (scalar keys %buffer == 0) {
			last;
		}elsif(scalar keys %buffer < $loop_control_flag){

			%undetermined = %buffer;
			if (scalar keys %buffer == $failure_flag) {

				warn "determine linkage phase failed for loci:\n";
				warn join("\t",keys %buffer),"\n";

				open (ERROR,">$fKey.dlp.udl") or die $!;
				
				print ERROR $_,"\n" for(keys %buffer);
				
				close (ERROR) ;

				last;
			}
			$failure_flag = scalar keys %buffer;
		}

	#	print scalar keys %buffer,"\n";
	#	<STDIN>;

		while (1) {

			last if (scalar keys %undetermined == 0) ;

			## choose the next loci 

			my $nextLoci = chooseNextLoci(\@determined,\%undetermined,\%pwd);

			print $log "next loci : $nextLoci \t $G{$nextLoci}{'type'} \n";

		#	die;
			
			## determine all the possible linkage phase for the next loci

			my %possibleLP = PossibleLP(\@determined,$nextLoci,\%G,\%pwd);

			print $log $_,"\t",join("\t",keys %{$possibleLP{$_}}),"\n" for(keys %possibleLP);

			my @maxiLikelihoodLP = std_dlp(\%possibleLP);

			if (@maxiLikelihoodLP == 1) {

				$G{$nextLoci}{'lp'} = $maxiLikelihoodLP[0];
				$nloc++;

				push @determined,$nextLoci;
				delete $undetermined{$nextLoci};
				delete $buffer{$nextLoci};
				
				print $log "\n--- LOCI --- TYPE --- LP -----\n";

				foreach my $dLoci (@determined) {
					print $log $dLoci,"\t",$G{$dLoci}{'type'},"\t{",$G{$dLoci}{'lp'},"}\n";
				}

				print $log "\n\n";
			
			}else{

				print $log "Sorry,It can't determine linkage Phase of $nextLoci currently \n";
				
				delete $undetermined{$nextLoci};

			}
			
		}
	}

	close ($log);

	#-------------------------------------------------------------------
	# Print
	#-------------------------------------------------------------------
	open (OUT,">$fKey.loc")|| die $!;

	print OUT ";",GetTime(),"\n";
	print OUT "name = ",basename($fKey),"\n";
	print OUT "popt = CP\n";
	print OUT "nloc = ",$nloc,"\n";
	print OUT "nind = $nind\n\n";

	foreach my $marker (keys %G) {

		if (exists $G{$marker}{'lp'}) {

			print OUT join("\t",($marker,$G{$marker}{'type'},"{$G{$marker}{'lp'}}",@{$G{$marker}{'indi'}})),"\n";
		
		}
	}

	close (OUT) ;
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub std_dlp {#
	my ($pos_lp) = @_;
	
	my %count = ();

	foreach my $lp (keys %{$pos_lp}) {

		my @lp = split /,/,$lp;

		foreach (@lp) {

			$count{$_} += scalar keys %{$pos_lp->{$lp}};
		}
	}

	## \\\\\\\\ DEBUG \\\\\\\\\\\\\\\\\

	print $log $_,"  <",$count{$_},">\t" for(sort {$count{$b} <=> $count{$a}} keys %count);
	print $log "\n";

#	<STDIN>;

	##\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	my ($maxFreq,$submaxFreq) = sort {$count{$b} <=> $count{$a}} keys %count;

#	$count{$submaxFreq}||=0;
	
	if (defined $submaxFreq) {

		my $chiStat = log($count{$maxFreq}/$count{$submaxFreq});
		my $chiProb = Statistics::Distributions::chisqrprob (1,$chiStat);

		print $log "ChiValue = $chiStat, p = $chiProb\n";


	}else{

		print $log "ChiValue = infinity, p = 0\n";
	}
	

	return grep {$count{$_} == $count{$maxFreq}} keys %count;

} 

sub PossibleLP {#
	
	my ($determined,$LociToDetermine,$ref_G,$ref_pwd) = @_;
	my %result = ();

	foreach my $dLoci (@{$determined}) {

#		print "lp with $dLoci ($ref_G->{$dLoci}{'type'}) :$ref_pwd->{$dLoci}{$LociToDetermine}{'lp'}\n";
		
#		my $curLP = join(",",pairDetermine($ref_G->{$dLoci}{'lp'},$ref_G->{$LociToDetermine}{'type'},$ref_pwd->{$dLoci}{$LociToDetermine}{'lp'})) if (exists $ref_pwd->{$dLoci}{$LociToDetermine}) ;

#		print $curLP,"\n";

#		$result{$curLP}{$dLoci} = 1;

		if (exists $ref_pwd->{$dLoci}{$LociToDetermine}) {

			my $curLP = join(",",pairDetermine($ref_G->{$dLoci}{'lp'},$ref_G->{$LociToDetermine}{'type'},$ref_pwd->{$dLoci}{$LociToDetermine}{'lp'})) if ($ref_pwd->{$dLoci}{$LociToDetermine}{'r'} <= $rThreshold and $ref_pwd->{$dLoci}{$LociToDetermine}{'lod'} >= $lodThreshold) ;
			$result{$curLP}{$dLoci."($ref_G->{$dLoci}{'type'})"."{$ref_pwd->{$dLoci}{$LociToDetermine}{'lp'}}"."{$ref_G->{$dLoci}{'lp'}}"."{$ref_pwd->{$dLoci}{$LociToDetermine}{'lod'}}"} = 1 if (defined $curLP) ;
			
		}
	}

	return %result;
}

sub singleDetermine {#
	my ($cross_type) = @_;

	my ($p1,$p2,$m1,$m2) = $cross_type =~/(.)(.)x(.)(.)/;

	my ($plp,$mlp) = ();

	$plp = '-' if ($p1 eq $p2) ;
	$mlp = '-' if ($m1 eq $m2) ;
	
	$plp = '0' if (not defined $plp) ;
	$mlp = '0' if (not defined $mlp) ;
	
	return $plp.$mlp;
}

sub pairDetermine {#
	my ($dlp,$cross_type,$pairlp) = @_;

#	print join("\t",($dlp,$cross_type,$pairlp)),"\n";
#	<STDIN>;

	my ($p1,$p2,$m1,$m2) = $cross_type =~/(.)(.)x(.)(.)/;
	my ($plp,$mlp) = ();

	my @lp = split /,/,$pairlp;

#	print "@lp\n";

	my ($dplp,$dmlp) = $dlp =~/(.)(.)/;

	my @dplp = ($dplp eq '-')?(0,1):($dplp);
	my @dmlp = ($dmlp eq '-')?(0,1):($dmlp);

#	print "@dplp\n";
#	print "@dmlp\n";

	my (%plp,%mlp) = ();

	$plp{'-'}=1 if ($p1 eq $p2) ;
	$mlp{'-'}=1 if ($m1 eq $m2) ;

	foreach my $lp (@lp) {

		my ($pPhase,$mPhase) = $lp =~/(.)(.)/;
		if (not exists $plp{'-'}) {

			if ($pPhase eq 'R') {

				$plp{1-$_} = 1 for(@dplp);
			}else{

				$plp{$_} = 1 for(@dplp);
			}	
		}


		if (not exists $mlp{'-'}) {

			if ($mPhase eq 'R') {

				$mlp{1-$_} = 1 for(@dmlp);
			}else{

				$mlp{$_} = 1 for(@dmlp);
			}	
		}
	}
		
#	print "P possible phase : ",join(",",@plp),"\n";
#	print "M possible phase : ",join(",",@mlp),"\n";

	my @result = ();
	foreach my $p_lp (keys %plp) {
		foreach my $m_lp (keys %mlp) {
			push @result,$p_lp.$m_lp;
		}
	}

#	print "@result\n";

	return @result;

}

sub anchorFirstPair {#
	my ($ref_pwd,$ref_G) = @_;

	my @marker = keys %{$ref_G};
	my (@firstPair,@strongestLinkPair) = ();

	my $maxLOD = 0;
	my $slp_maxLOD = 0;

	for (my $i=0;$i<@marker-1 ;$i++) {
		for (my $j=$i+1;$j<@marker ;$j++) {
			if (defined $ref_pwd->{$marker[$i]}{$marker[$j]}) {

				my @lp = split /,/,$ref_pwd->{$marker[$i]}{$marker[$j]}{'lp'};

				if (@lp == 1) {

					if ($slp_maxLOD < $ref_pwd->{$marker[$i]}{$marker[$j]}{'lod'}) {

						$slp_maxLOD = $ref_pwd->{$marker[$i]}{$marker[$j]}{'lod'};

						@firstPair = ($marker[$i],$marker[$j]);
					}

				}

				if ($maxLOD < $ref_pwd->{$marker[$i]}{$marker[$j]}{'lod'}) {

					$maxLOD = $ref_pwd->{$marker[$i]}{$marker[$j]}{'lod'};

					@strongestLinkPair = ($marker[$i],$marker[$j]);

				}
			}	
		}
	}

	if (@firstPair == 0) {

		my @anchorLoci = grep {my ($p1,$p2,$m1,$m2) = $ref_G->{$_}{'type'}=~/(.)(.)x(.)(.)/; $p1 ne $p2 and $m1 ne $m2} keys %{$ref_G};

		if (@anchorLoci == 0) {

			@firstPair = @strongestLinkPair;
		}else{
			
			my $anchor_maxLOD = 0;

			foreach (@anchorLoci) {

				my ($max_marker) = sort {$ref_pwd->{$_}{$b}{'lod'} <=> $ref_pwd->{$_}{$a}{'lod'}} keys %{$ref_pwd->{$_}};

				if ($anchor_maxLOD < $ref_pwd->{$_}{$max_marker}{'lod'}) {

					$anchor_maxLOD = $ref_pwd->{$_}{$max_marker}{'lod'};

					@firstPair = ($_,$max_marker);
				}
			}
		}
	}

	return @firstPair;
}

sub chooseNextLoci {#
	my ($determined,$undetermined,$ref_pwd) = @_;

	my %sumLOD = ();
	foreach my $uLoci (keys %{$undetermined}) {
		
		my $sum = 0;

		foreach my $dLoci (@{$determined}) {

			next unless (exists $ref_pwd->{$uLoci}{$dLoci}{'lod'}) ;
			$sum += $ref_pwd->{$uLoci}{$dLoci}{'lod'};

		}

		$sumLOD{$sum}{$uLoci} = 1;
		
	}

	my ($max) = sort {$b <=> $a} keys %sumLOD;

	my @nextLociSet = keys %{$sumLOD{$max}};

	return $nextLociSet[0];
}


sub checkDataSet {#
	my ($pwdLociSet,$gLociSet) = @_;

	if (@{$gLociSet} == 0) {
		print "no loci found in genotype file!!!\n";
		exit;
	}

	my @overlap = @{intersection($pwdLociSet,$gLociSet)};
	my @missInfo = su($gLociSet,\@overlap);

	if (@missInfo != 0) {

		warn "The following Marekrs' pwd info dosen't exist in pwd file,please check\n";

		print join("\t",@missInfo),"\n";

		exit(-1);
	}
}

sub intersection {#
	my ($A,$B)=@_;
	my %uniqA=map {$_,1} @{$A};
	my %uniqB=map {$_,1} @{$B};
	my %merge=();
	my %overlap=();
	foreach  (keys %uniqA,keys %uniqB) {
		$merge{$_}++ && $overlap{$_}++;
	}
	my @result = keys %overlap;
	return \@result;
}

sub su {#
	my ($A,$B)=@_;
	my %uniqA=map {$_,1} @{$A};
	my %uniqB=map {$_,1} @{$B};
	my %merge=();
	my %overlap=();
	foreach  (keys %uniqA,keys %uniqB) {
		$merge{$_}++ && $overlap{$_}++;
	}
	my @result = grep {$merge{$_} == 1} keys %merge;
	return @result;
}

sub loadPWD {#
	my ($ref_pwd,$fIn) = @_;

	my %pwdMarker = ();

	open (IN,"$fIn") or die $!;
	while (<IN>) {

		chomp;
		next if (/^$/ or /^;/) ;

		my ($marker1,$marker2,$r,$mLOD,$lLOD,$fi,$integrity,$lp) =split;

		$ref_pwd->{$marker1}{$marker2}{'r'}=$r;
		$ref_pwd->{$marker2}{$marker1}{'r'}=$r;
		$ref_pwd->{$marker1}{$marker2}{'lod'}=$mLOD;
		$ref_pwd->{$marker2}{$marker1}{'lod'}=$mLOD;
		$ref_pwd->{$marker1}{$marker2}{'lp'}=$lp;
		$ref_pwd->{$marker2}{$marker1}{'lp'}=$lp;

		$pwdMarker{$marker1}=1;
		$pwdMarker{$marker2}=1;
		
	}
	
	close (IN) ;
	return keys %pwdMarker;
}

sub loadGtype {#
	my ($ref_G,$fIn) = @_;

	open (IN,"$fIn") or die $!;
#	my $head = <IN>;
#	my @head = split /\s+/,$head;
#	$nind = @head - 2;
	my $return=0;
	while (<IN>) {
		chomp;
		next if (/^$/ or /^;/ or /[Tt]ype/ or /^name/ or /^popt/ or /^nloc/ or /^nind/ or /^#/) ;

		my ($marker,$type,@indi) = split;
		my $phase;
		if (/\{..\}/) {
			($marker,$type,$phase,@indi)=split;
			$return=1;
		}
		$type=~s/<|>//g;
		$ref_G->{$marker}{'type'} = $type;
		$ref_G->{$marker}{'indi'} = \@indi;
		$ref_G->{$marker}{'phase'}= $phase;
		$nind = @indi unless (defined $nind) ;
	}
	close (IN) ;
	return $return;
}

sub sum {#
	my (@arr) = @_;

	my $sum = 0;

	$sum+=$_ for(@arr);

	return $sum;
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
			This program is used to determine linkage phase for CP population.

			The format of pwd file is like this:

			============================================================================================================
				;2012-03-27 16:20:41
				;4325 pairs out of 4950, 625 not informative
				;name = LG1

				LOCI_1	LOCI_2	R	Modified_LOD	Linkage_LOD		Fisher_Information	Possible_Linkage_Phase
				Marker0 Marker1 0.0000000000000000e+00  2.9447218304520607e+01  3.0102999566398115e+01  inf     1.0000  CC,CR
				Marker0 Marker2 0.0000000000000000e+00  2.9481967655774760e+01  3.0102999566398115e+01  inf     1.0000  RR,CR
			=============================================================================================================

			Above pwd  file can be obtained using programe:

				$Bin/calculatePWDforCP.pl


Usage:
  Options:
  -p	<file>	PWD file , forced
  -g	<file>	Genotype file, forced
  -k	<str>	Key of output file,forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -r	<float>	Maximum  threshold for r,optional,default [0.4]
  -LOD	<float>	Minimum thresshold for LOD score, optional, default [1.0]
  -h		Help

USAGE
	print $usage;
	exit;
}