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
my ($Map,$fmap,$fix_list,$fKey,$dOut,$log,$checkDepth,$checkMarker,$ratio,$miss,$imput_mode,$correct,$loc,$phaseFileType,$extend,$sus_ratio,$diff_ratio,$unsure_ratio,$dis_threshold,$window_type,$cycle_threshold);
GetOptions(
				"help|?" =>\&USAGE,
				"m:s"=>\$Map,
				"l:s"=>\$loc,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				"imput:s"=>\$imput_mode,

				"fix:s"=>\$fix_list,

				"window_type:s"=>\$window_type,
				"marker:i"=>\$extend,
				"map:f"=>\$dis_threshold,
				"fmap:s"=>\$fmap,
				"sus:f"=>\$sus_ratio,
				"diff_ratio"=>\$diff_ratio,
				"unsure_ratio"=>\$unsure_ratio,
				"cycle_threshold:s"=>\$cycle_threshold,
				) or &USAGE;
&USAGE unless ($Map and $loc and $fKey);
#------------------------------------------------------------------
# Global parameter  
#------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$fKey="$dOut/$fKey";

$window_type||="marker";
$sus_ratio||=0.06;
$dis_threshold||=20;

$imput_mode         =     defined $imput_mode ? $imput_mode : 1;
$extend             =     defined $extend ? $extend : 6;
$diff_ratio         =     defined $diff_ratio ? $diff_ratio : 0.85;
$unsure_ratio       =     defined $unsure_ratio ? $unsure_ratio : 0.5;
$cycle_threshold    =     defined $cycle_threshold ? $cycle_threshold : 3;

if ($window_type ne "marker" && $window_type ne "map" ) {
	print "-window_type $window_type is not both \"marker\" and \"map\", wrong value!\n";
	die;
}

#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------
my (%haplo,@marker,%indiToCheck,%depth,%genotype,%allel,%depthStat,%cross_type,%G,%linkagePhase)=();
my %fix = ();   ## fix marker 
my %suspicious;
my %sus_marker;
my %map_dis;
my %map;
my @map;

my %order=();
my %order_marker=();
my $order=0;
my %correct;
my $off_sum;

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------

######### Get loc #################
my (@LOChead)=();
open (L,"$loc") or die $!;
$/="\n";
while (<L>) {
	chomp;
	next if (/^$/ || /^;/) ;
	if (/^name/|| /^popt/|| /^nloc/|| /^nind/) {
		if (/^nind/) {
			(undef,$off_sum)=split /=/,$_;
		}
		push @LOChead,$_;
		next;
	}
	my ($marker,$cross_type,$linkPhase,@geno)=split;
	$cross_type=~s/<|>//g;
	$G{"original"}{$marker}=\@geno;
	$cross_type{$marker}=$cross_type;
	$linkPhase=~s/[\{\}]//g;
	$linkagePhase{$marker}=$linkPhase;
}
close (L) ;

#### get fix marker info 
if (defined $fix_list and -f $fix_list) {
	open In,$fix_list;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ or /^;/);
		my @loci = split /\s+/,$_;
		$fix{$_}=1 for(@loci);
	}
	close In;
}

#######  Get map file ##############
open (M,"$Map") or die $!;
while (<M>) {
	chomp;
	next if (/^$/ || /^group/ || /^;/) ;

	my ($loci) = split /\s+/,$_;

	if (exists $G{'original'}{$loci}) {
		for (my $i=0;$i<@{$G{'original'}{$loci}} ;$i++) {
			my $haplo = determineHaploSource($cross_type{$loci},$linkagePhase{$loci},$G{'original'}{$loci}->[$i]);
			my @unit = split //,$haplo; 
			push @{$haplo{$i}},[@unit];
		}
	}else{
		print "info of $loci in map file not found in loc file \n";
		next;
	}
	
	push @marker,$loci;
	$order_marker{$order} = $loci;
	$order{$loci} = $order;
	$order++;

}
close (M) ;

my $marker_sum=@marker;

########### Get map #######################
if($window_type eq "map"){
	open (M,$fmap) or die $!;
	$/="\n";
	while (<M>) {
		chomp;
		next if($_=~/^$/ || /^;/ || /unmap/ || /group/);
		next if($_!~/Marker/);
		my ($marker,$map_locus)=split /\s+/,$_;
		$map{$marker}=$map_locus;
	}
	close (M) ;

	for (my $i=0;$i<@marker;$i++) {
		foreach my $marker2 (keys %map) {

			next if($marker2 eq $marker[$i] || !exists $map{$marker[$i]});
			$map_dis{$marker[$i]}{$marker2}=abs($map{$marker[$i]}-$map{$marker2});
		}
	}

}


#-------------------------------------------------------------------
# Correct process 
#-------------------------------------------------------------------

############# correct halpoSource ################
open (INF,">$fKey.info") or die $!;
open (SUS,">$fKey.suspi.stat") or die $!;

####统计marker的可疑子代数
foreach my $indi (sort {$a <=> $b} keys %haplo) {
	for (my $chain_num=0;$chain_num<2;$chain_num++) {
		for (my $i=0;$i<@{$haplo{$indi}} ;$i++) {
			next if($haplo{$indi}->[$i]->[$chain_num] ne "1" && $haplo{$indi}->[$i]->[$chain_num] ne "2");
			my $before=$haplo{$indi}->[$i]->[$chain_num];
			my $after=$haplo{$indi}->[$i]->[$chain_num];
			for (my $j=$i-1;$j>0;$j--) {
				if ($haplo{$indi}->[$j]->[$chain_num] eq "1" || $haplo{$indi}->[$j]->[$chain_num] eq "2") {
					$before=$haplo{$indi}->[$j]->[$chain_num];
					last;
				}
			}
			for (my $j=$i+1;$j<@{$haplo{$indi}};$j++) {
				if ($haplo{$indi}->[$j]->[$chain_num] eq "1" || $haplo{$indi}->[$j]->[$chain_num] eq "2") {
					$after=$haplo{$indi}->[$j]->[$chain_num];
					last;
				}
			}
			if ($before eq $after && $haplo{$indi}->[$i]->[$chain_num] ne $before ) {
				$suspicious{$marker[$i]}{$indi}=1;
			}
		}
	}
}
foreach my $marker (keys %suspicious) {
	
	my $sus_off_sum=scalar keys %{$suspicious{$marker}};
	print SUS $marker,"\t",$sus_off_sum,"\n";
	if ($sus_off_sum>$off_sum*$sus_ratio) {
		$sus_marker{$marker}=$sus_off_sum;
	}
}
close (SUS) ;

print "Suspicious marker:\n";
foreach my $marker (sort {$a cmp $b} keys %sus_marker ) {
	print $marker,"\t";
}
print "\n";

######## correct ########


&correct($extend);
my $cycle=1;
while ($cycle++<$cycle_threshold ) {
	&correct($extend);
}


###### correct gene ################
foreach my $off_num (keys %correct) {
	foreach my $marker_num (keys %{$correct{$off_num}}) {
		my $gene=&haplo2gene($cross_type{$order_marker{$marker_num}},$linkagePhase{$order_marker{$marker_num}},$haplo{$off_num}->[$marker_num]->[0],$haplo{$off_num}->[$marker_num]->[1]);
		$gene = $G{'original'}{$order_marker{$marker_num}}->[$off_num] if ($G{'original'}{$order_marker{$marker_num}}->[$off_num] ne '--' && $gene eq '--') ;  ## don't miss hk to -- 
		$G{"correct"}{$order_marker{$marker_num}}{$off_num}=$gene;
	}
}

#-------------------------------------------------------------------
# Print
#-------------------------------------------------------------------

### corrected loc file 
open (C,">$fKey.correct.loc") or die $!;
print C ";",join(",",&GetTime()),"\n\n";
foreach my $head (@LOChead) {
	print C $head,"\n";
}
print C "\n";

foreach my $marker (@marker) {
	print C join("\t",($marker,"<$cross_type{$marker}>","{$linkagePhase{$marker}}")),"\t";
	for (my $i=0;$i<$off_sum;$i++) {
		if (defined $G{"correct"}{$marker}{$i}) {
			print C $G{"correct"}{$marker}{$i},"\t";
		}else{
			print C $G{"original"}{$marker}->[$i],"\t";
		}
	}
	print C "\n";
}
close (C) ;

### corrected log file 
open (LOG,">$fKey.correct.loc.log") or die $!;
foreach my $marker (keys %indiToCheck) {
	foreach my $indi (keys %{$indiToCheck{$marker}}) {
		print LOG $marker,"\t",$indi,"\t",$G{"original"}{$marker}->[$indi],"\t",$G{"correct"}{$marker}{$indi},"\n" if ($G{"original"}{$marker}->[$indi] ne $G{"correct"}{$marker}{$indi}) ;
	}
}
close (LOG) ;

### corrected phase file 
open (OUT,">","$fKey.correct.phase") or die $!;
for (my $i=0;$i<$marker_sum;$i++) {
	print OUT $order_marker{$i},"\t",$cross_type{$order_marker{$i}},"\t","{$linkagePhase{$order_marker{$i}}}","\t";
	for (my $j=0;$j<$off_sum;$j++) {
		print OUT $haplo{$j}->[$i]->[0],$haplo{$j}->[$i]->[1],"\t";
	}
	print OUT "\n";
}
close (OUT) ;

### corrected phase log file
open (LOG,">$fKey.correct.phase.log") or die $!;
foreach my $off_num (sort {$a <=> $b} keys %correct) {
	foreach my $marker_num (keys %{$correct{$off_num}}) {
		foreach my $chain_num (sort {$a <=> $b} keys %{$correct{$off_num}{$marker_num}}) {
			print LOG $off_num,"\t",$chain_num,"\t",$order_marker{$marker_num},"\t",$correct{$off_num}{$marker_num}{$chain_num}{"original"},"\t",$correct{$off_num}{$marker_num}{$chain_num}{"correct"},"\n";
		}
	}
}
close (LOG) ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub determineHaploSource {#
	my ($type,$linkPhase,$progenyGeno)=@_;

	return "--" if ($progenyGeno eq "--") ;

	my $haploSource='';;
	my (%parentAllel,@gameteCombination)=();

	my %haploIndex=(
		"0"=>"11",
		"1"=>"12",
		"2"=>"21",
		"3"=>"22",
	);

	my ($p1,$p2,$m1,$m2)= $type =~/(.)(.)x(.)(.)/ ;
	@{$parentAllel{"P"}}=($p1,$p2);
	@{$parentAllel{"M"}}=($m1,$m2);

	my ($PLinkPhase, $MLinkPhase) = split //,$linkPhase ;
	if ($PLinkPhase eq '1'){
		@{$parentAllel{"P"}} = reverse @{$parentAllel{"P"}};
	}
	if ($MLinkPhase eq '1'){
		@{$parentAllel{"M"}} = reverse @{$parentAllel{"M"}} ;
	}

	foreach my $pGamete (@{$parentAllel{"P"}}) {
		foreach my $mGamete (@{$parentAllel{"M"}}) {
			push @gameteCombination,$pGamete.$mGamete;
		}
	}

	for (my $j=0;$j<@gameteCombination ;$j++) {
		
		my @haplo=split //,$gameteCombination[$j];
		my $haplo=join("",sort {$a cmp $b} @haplo);

		if ($haplo eq $progenyGeno) {

			my ($allel1,$allel2) = $progenyGeno =~/(.)(.)/;

			if ($p1 eq $p2) {

				my ($lp) = $haploIndex{$j} =~/.(.)/;
				$haploSource = "0".$lp;
				last;
				
			}elsif ($m1 eq $m2) {

				my ($lp) = $haploIndex{$j} =~/(.)./;
				$haploSource = $lp."0";
				last;

			}elsif (join("",sort {$a cmp $b} ($p1,$p2)) eq join("",sort {$a cmp $b} ($m1,$m2)) and $allel1 ne $allel2) {
				
				$haploSource = "--";
				last;

			}else{
				
				$haploSource = $haploIndex{$j};
			}
		}
		
	}
	return $haploSource;

}

sub correct {#
	my $extend=$_[0];
	for (my $off_num=0;$off_num<$off_sum;$off_num++) {
		my @chain=();
		for (my $i=0;$i<$marker_sum;$i++) {
			my ($chain1,$chain2)=@{$haplo{$off_num}->[$i]};
			push @{$chain[0]},$chain1;
			push @{$chain[1]},$chain2;
		}

		for (my $chain_num=0;$chain_num<2;$chain_num++) {
			my $mark=0;
			print INF "deal with off:$off_num\t","chain:",$chain_num,"\n";

			####处理整条亲本来源完全一致以及只有一个marker不一致的情况
			my $total_sum_1=0;
			my $total_sum_2=0;
			for (my $i=0;$i<@{$chain[$chain_num]};$i++) {
				if ($chain[$chain_num][$i] eq "1") {
					$total_sum_1++;
				}elsif($chain[$chain_num][$i] eq "2"){
					$total_sum_2++;
				}
			}
			if ($total_sum_2<2 && $total_sum_1>10) {
				for (my $i=0;$i<@{$chain[$chain_num]};$i++) {
					next if($chain[$chain_num][$i] eq "1" || $chain[$chain_num][$i] eq "0");
					next if (exists $fix{$marker[$i]}) ; ## important 
					$haplo{$off_num}->[$i]->[$chain_num]='1';
					$correct{$off_num}{$i}{$chain_num}{"original"}=$chain[$chain_num][$i];
					$correct{$off_num}{$i}{$chain_num}{"correct"}='1';
					$indiToCheck{$order_marker{$i}}{$off_num}=1;
				}
				next;
			}elsif($total_sum_1<2 && $total_sum_2>10) {
				for (my $i=0;$i<@{$chain[$chain_num]};$i++) {
					next if($chain[$chain_num][$i] eq "2" || $chain[$chain_num][$i] eq "0");
					next if (exists $fix{$marker[$i]}) ; ## important 
					$haplo{$off_num}->[$i]->[$chain_num]='2';
					$correct{$off_num}{$i}{$chain_num}{"original"}=$chain[$chain_num][$i];
					$correct{$off_num}{$i}{$chain_num}{"correct"}='2';
					$indiToCheck{$order_marker{$i}}{$off_num}=1;
				}
				next;
			}

			####具体分析处理

			for (my $to_be_confirm=0;$to_be_confirm<$marker_sum;) {
				if(exists $sus_marker{$marker[$to_be_confirm]}){
					print INF "suspicous marker\t$marker[$to_be_confirm]\n";
					$to_be_confirm++;
					print INF "to_be_confirm:\t",$to_be_confirm,"\n";
					next;
				}

				if (exists $fix{$marker[$to_be_confirm]}) { ### add by macx 2012-10-24
					$to_be_confirm++;
					next;
				}

				my $start;
				my $stop;
				my %sum;
				my %downstream_sum;
				my %upstream_sum;

				if($window_type eq "map"){
					my @window_marker=();
					foreach my $marker (keys %{$map_dis{$order_marker{$to_be_confirm}}}) {
	
						if ($map_dis{$order_marker{$to_be_confirm}}{$marker}<=$dis_threshold) {
							push @window_marker,$marker;
	
							}
						}
						if (@window_marker<5) {

							$mark++;
							if ($mark>3) {
								print "The window by map distance is too small, please up the parameter:-map\n"; 
								die;
							}
							$to_be_confirm++;
							next;
						}
						my @window_marker_sorted=sort {$order{$a} <=> $order{$b}} (@window_marker);
						if ($order{$window_marker_sorted[0]}>$to_be_confirm) {
							$start=$to_be_confirm;
						}else{
							$start=$order{$window_marker_sorted[0]};
						}
						if ($order{$window_marker_sorted[$#window_marker_sorted]}<$to_be_confirm) {
							$stop=$to_be_confirm+1;
						}else{
							$stop=$order{$window_marker_sorted[$#window_marker_sorted]};
						}
	
					}elsif($window_type eq "marker"){
						if ($to_be_confirm-$extend<0) {
							$start=0;
						}else{
							$start=$to_be_confirm-$extend;
						}
						if ($to_be_confirm+$extend>$marker_sum) {
							$stop=$marker_sum;
						}else{
							$stop=$to_be_confirm+$extend;
						}
					}
					print INF "\n";
					print INF $order_marker{$to_be_confirm}," ",$to_be_confirm,"\t",$off_num,"\t",$chain_num,"\n";
					print INF "start\t",$start,"\t","stop\t",$stop,"\n";

					for (my $j=$start;$j<$stop;$j++) {
						next if(exists $sus_marker{$marker[$j]});
						push @{$sum{$chain[$chain_num][$j]}},$j;
						if ($j>$to_be_confirm) {
							push @{$downstream_sum{$chain[$chain_num][$j]}},$j;
						}elsif($j<$to_be_confirm){
							push @{$upstream_sum{$chain[$chain_num][$j]}},$j;
						}
					}
					foreach my $haplo (keys %sum) {
						print INF $haplo,":\n";
						for (my $i=0;$i<@{$sum{$haplo}};$i++) {
							print INF $sum{$haplo}->[$i],"\t";
						}
						print INF "\n";
					}
					print INF "upstream:\n";
					foreach my $haplo (keys %upstream_sum) {
						print INF $haplo,":\n";
						for (my $i=0;$i<@{$upstream_sum{$haplo}};$i++) {
							print INF $upstream_sum{$haplo}->[$i],"\t";
						}
						print INF "\n";
					}
					print INF "downstream:\n";
					foreach my $haplo (keys %downstream_sum) {
						print INF $haplo,":\n";
						for (my $i=0;$i<@{$downstream_sum{$haplo}};$i++) {
							print INF $downstream_sum{$haplo}->[$i],"\t";
						}
						print INF "\n";
					}


			
					if ((!defined $sum{'1'} && !defined $sum{'-'}) || (!defined $sum{'2'} && !defined $sum{'-'})) {
						print INF "(!defined \$sum{'1'} && !defined \$sum{'-'}) || (!defined \$sum{'2'} && !defined \$sum{'-'})","\t",$order_marker{$to_be_confirm}," ($to_be_confirm)\t",$off_num,"\t",$chain_num,"\n";
						$to_be_confirm=$stop;
						print INF "to_be_confirm\t",$to_be_confirm,"\n";
						next;
					}elsif((!defined $sum{'1'} && defined $sum{'-'}) || (!defined $sum{'2'} && defined $sum{'-'})){
						if ($imput_mode) {
							my $haplo;
							if (defined $sum{'1'}) {
								$haplo='1';
							}elsif(defined $sum{'2'}){
								$haplo='2';
							}else{
								print INF "!defined \$sum{'1'} && !defined \$sum{'2'}\t$order_marker{$to_be_confirm} ($to_be_confirm)\t$off_num\t$chain_num\n";
								$to_be_confirm=$stop;
								print INF "to_be_confirm\t",$to_be_confirm,"\n";
								next;
							}

							if ($chain[$chain_num][$to_be_confirm] eq '-') {
								if ($stop!=$marker_sum) {
									if (($stop-$to_be_confirm>5 && (!exists $downstream_sum{$haplo} || @{$downstream_sum{$haplo}}<0.1*($stop-$to_be_confirm)) ) || ($to_be_confirm-$start>5 && (!defined $upstream_sum{$haplo} || @{$upstream_sum{$haplo}}<0.1*($to_be_confirm-$start) ) ) ) {
										print INF "upstream and downstream\n";
										my ($num1,$num_,$min_num,$mid_num,$max_num);
										if (defined $downstream_sum{$haplo}) {
											$num1=$downstream_sum{$haplo}->[0];
										}else{
											$num1=$stop;
										}
										if ($imput_mode) {
											if (defined $downstream_sum{'-'}) {
												$num_=$downstream_sum{'-'}->[0];
											}else{
												$num_=$stop;
											}
											($min_num,$max_num)=sort {$a <=> $b} ($num1,$num_);
										}else{
											$min_num=$num1;
										}
										$to_be_confirm=$min_num;

										print INF "to_be_confirm\t",$to_be_confirm,"\n";
										next;

									}
							}
							$haplo{$off_num}->[$to_be_confirm]->[$chain_num]=$haplo;
							$correct{$off_num}{$to_be_confirm}{$chain_num}{"original"}='-';
							$correct{$off_num}{$to_be_confirm}{$chain_num}{"correct"}=$haplo;
							$indiToCheck{$order_marker{$to_be_confirm}}{$off_num}=1;
							my $s_=@{$sum{'-'}};
							my $s_haplo=@{$sum{$haplo}};


							print INF $order_marker{$to_be_confirm},"\t",$off_num,"\t",$chain_num,"\t",'-',"\t",$s_,"\t",$haplo,"\t",$s_haplo,"\n";
							for (my $k=0;$k<@{$chain[$chain_num]};$k++) {
								print INF "\t",$chain[$chain_num][$k];
								if ($k==$to_be_confirm) {
									print INF "*";
								}
								print INF "\n";
							}
						}

						my $sign=0;
						for (my $j=0;$j<@{$sum{'-'}};$j++) {
							if ($sum{'-'}->[$j]>$to_be_confirm) {
								$to_be_confirm=$sum{'-'}->[$j];
								$sign=1;
								last;
							}
						}
	
						if ($sign==0) {
							$to_be_confirm=$stop;
						}
						print INF "to_be_confirm\t",$to_be_confirm,"\n";
						next;
					}else{
						$to_be_confirm=$stop;
						print INF "NOT Imput!\n";
						print INF "to_be_confirm\t",$to_be_confirm,"\n";
						next;
					}
				}

				my $sum1=@{$sum{'1'}};
				my $sum2=@{$sum{'2'}};
				my $sum0;
				my $sum_;
				if (defined $sum{'0'}) {
					$sum0=@{$sum{'0'}};
				}else{
					$sum0=0;
				}
				if (defined $sum{'-'}) {
					$sum_=@{$sum{'-'}};
				}else{
					$sum_=0;
				}
				my ($less_sum,$more_sum)=sort {$a <=> $b} ($sum1,$sum2);
				if($sum0+$sum_>0.3*($stop-$start)){
					if ((abs($sum1-$sum2)<$more_sum*$diff_ratio) && $less_sum!=1 ) {
						print INF "abs(\$sum1-\$sum2) && \$less_sum!=1\t",abs($sum1-$sum2),"\t",$more_sum*$diff_ratio,"\t$order_marker{$to_be_confirm} ($to_be_confirm)\t$off_num\t$chain_num\n";
						$to_be_confirm++;
						print INF "to_be_confirm\t",$to_be_confirm,"\n";
						next;
					}elsif($less_sum==1 && $more_sum<2 ){
						print INF "\$less_sum==1 && \$more_sum<2\t",$less_sum,"\t",$more_sum*$diff_ratio,"\t$order_marker{$to_be_confirm} ($to_be_confirm)\t$off_num\t$chain_num\n";
						$to_be_confirm++;
						print INF "to_be_confirm\t",$to_be_confirm,"\n";
						next;
					}
				}else{
					if ((abs($sum1-$sum2)<$more_sum*0.8)) {
						print INF "(abs(\$sum1-\$sum2)<\$more_sum*0.7\t",abs($sum1-$sum2),"\t",$more_sum*$diff_ratio,"\t$order_marker{$to_be_confirm} ($to_be_confirm)\t$off_num\t$chain_num\n";
						$to_be_confirm++;
						print INF "to_be_confirm\t",$to_be_confirm,"\n";
						next;
					}
				}

				my ($less,$more)=sort {@{$sum{$a}} <=> @{$sum{$b}}} ('1','2');

				if ($chain[$chain_num][$to_be_confirm] eq $less ) {
					if ($stop!=$marker_sum) {
						if ((exists $downstream_sum{$less} && $stop-$to_be_confirm>5 && (!defined $downstream_sum{$more} || @{$downstream_sum{$more}}<0.15*($stop-$to_be_confirm)) ) || (!exists $downstream_sum{$less} && $stop-$to_be_confirm>15 && !defined $downstream_sum{$more}) || (exists $upstream_sum{$less} && $to_be_confirm-$start>5 && (!defined $upstream_sum{$more} || @{$upstream_sum{$more}}<0.15*($to_be_confirm-$start))) || (!exists $upstream_sum{$less} && $to_be_confirm-$start>15 && !defined $upstream_sum{$more}) ) {
							print INF "upstream and downstream\n";
							my ($num1,$num2,$num_,$min_num,$mid_num,$max_num);
							if (defined $downstream_sum{'1'}) {
								$num1=$downstream_sum{'1'}->[0];
							}else{
								$num1=$stop;
							}
							if (defined $downstream_sum{'2'}) {
								$num2=$downstream_sum{'2'}->[0];
							}else{
								$num2=$stop;
							}
							if ($imput_mode) {
								if (defined $downstream_sum{'-'}) {
									$num_=$downstream_sum{'-'}->[0];
								}else{
									$num_=$stop;
								}
								($min_num,$mid_num,$max_num)=sort {$a <=> $b} ($num1,$num2,$num_);
							}else{
								($min_num,$max_num)=sort {$a <=> $b} ($num1,$num2);
							}
							$to_be_confirm=$min_num;

							print INF "to_be_confirm\t",$to_be_confirm,"\n";
							next;
						}
					}

					$correct{$off_num}{$to_be_confirm}{$chain_num}{"original"}=$less;
					$correct{$off_num}{$to_be_confirm}{$chain_num}{"correct"}=$more;
					$haplo{$off_num}->[$to_be_confirm]->[$chain_num]=$more;
					$indiToCheck{$order_marker{$to_be_confirm}}{$off_num}=1;
					my $t=@{$sum{$less}};
					my $t2=@{$sum{$more}};
					print INF $order_marker{$to_be_confirm},"\t",$off_num,"\t",$chain_num,"\t",$less,"\t",$t,"\t",$more,"\t",$t2,"\n";
					for (my $k=0;$k<@{$chain[$chain_num]};$k++) {
						print INF "\t",$chain[$chain_num][$k];
						if ($k==$to_be_confirm) {
							print INF "*";
						}
						print INF "\n";
					}
				}

				if ($imput_mode) {
					if ($chain[$chain_num][$to_be_confirm] eq '-') {
						if ($stop!=$marker_sum) {
						if ((exists $downstream_sum{$less} && $stop-$to_be_confirm>5 && (!defined $downstream_sum{$more} || @{$downstream_sum{$more}}<0.15*($stop-$to_be_confirm)) ) || (!exists $downstream_sum{$less} && $stop-$to_be_confirm>15 && !defined $downstream_sum{$more}) || (exists $upstream_sum{$less} && $to_be_confirm-$start>5 && (!defined $upstream_sum{$more} || @{$upstream_sum{$more}}<0.15*($to_be_confirm-$start))) || (!exists $upstream_sum{$less} && $to_be_confirm-$start>15 && !defined $upstream_sum{$more}) ) {
								print INF "upstream and downstream\n";
								my ($num1,$num2,$num_,$min_num,$mid_num,$max_num);
								if (defined $downstream_sum{'1'}) {
									$num1=$downstream_sum{'1'}->[0];
								}else{
									$num1=$stop;
								}
								if (defined $downstream_sum{'2'}) {
									$num2=$downstream_sum{'2'}->[0];
								}else{
									$num2=$stop;
								}
								if ($imput_mode) {
									if (defined $downstream_sum{'-'}) {
										$num_=$downstream_sum{'-'}->[0];
									}else{
										$num_=$stop;
									}
									($min_num,$mid_num,$max_num)=sort {$a <=> $b} ($num1,$num2,$num_);
								}else{
									($min_num,$max_num)=sort {$a <=> $b} ($num1,$num2);
								}
								$to_be_confirm=$min_num;

								print INF "to_be_confirm\t",$to_be_confirm,"\n";
								next;
							}
						}
						$correct{$off_num}{$to_be_confirm}{$chain_num}{"original"}='-';
						$correct{$off_num}{$to_be_confirm}{$chain_num}{"correct"}=$more;
						$haplo{$off_num}->[$to_be_confirm]->[$chain_num]=$more;
						$indiToCheck{$order_marker{$to_be_confirm}}{$off_num}=1;
						my $t2=@{$sum{$more}};
						print INF $order_marker{$to_be_confirm},"\t",$off_num,"\t",$chain_num,"\t",$more,"\t",$t2,"\n";
						for (my $k=0;$k<@{$chain[$chain_num]};$k++) {
							print INF "\t",$chain[$chain_num][$k];
							if ($k==$to_be_confirm) {
								print INF "*";
							}
							print INF "\n";
						}
					}
				}

				my @to_off_num=();
				if ($imput_mode) {
					if (defined $sum{'-'}) {
						@to_off_num=sort {$a <=> $b} (@{$sum{$less}},@{$sum{'-'}});
					}else{
						@to_off_num=sort {$a <=> $b} @{$sum{$less}};
					}
				}else{
					@to_off_num=sort {$a <=> $b} @{$sum{$less}};
				}
				my $sign=0;
				for (my $j=0;$j<@to_off_num;$j++) {
					if ($to_off_num[$j]>$to_be_confirm) {
						$to_be_confirm=$to_off_num[$j];
						$sign=1;
						last;
					}
				}
				if ($sign==0) {
					$to_be_confirm=$stop;
				}

				print INF "to_be_confirm\t",$to_be_confirm,"\n";
				next;
			}
		}
	}
}

sub haplo2gene {#条件：单体来源判断是根据“连锁相0表示正，1表示反”来判断的
	my ($cross_type,$phase,@haplo)=@_;
	my $gene;
	my @gene;
	$cross_type=~s/x//g;
	my @cross_type = split //,$cross_type;
	my @phase = split //,$phase;
	for (my $i=0;$i<2;$i++) {
		if ($phase[$i] eq '1') {
			my $t=$cross_type[$i*2];
			$cross_type[$i*2] = $cross_type[$i*2+1];
			$cross_type[$i*2+1] = $t;
		}
		if ($haplo[$i] eq '1') {
			$gene[$i]=$cross_type[$i*2];
		}elsif($haplo[$i] eq '2'){
			$gene[$i]=$cross_type[$i*2+1];
		}else{
			if ($cross_type[$i*2] eq $cross_type[$i*2+1]) {
				$gene[$i] = $cross_type[$i*2];
			}else{
				$gene[$i]='-';
				$gene='--';
			}
		}
	}
	my @unit=sort {$a cmp $b} @gene;
	if (!defined $gene) {
		$gene=$unit[0].$unit[1];
	}
	return $gene;
}

sub extend_length {#
	my ($ref_chain,$marker_num)=@_;
	my $marker_sum=@{$ref_chain};
	my $length=0;
	my $point=$marker_num;
	for (my $i=$marker_num;$i>=0;$i--) {
		if ($ref_chain->[$i] eq '0' || $ref_chain->[$i] eq '-'|| $ref_chain->[$i] eq $ref_chain->[$marker_num]) {
			$length++;
		}else{
			last;
		}
	}
	for (my $i=$marker_num;$i<$marker_sum;$i++) {
		if ($ref_chain->[$i] eq '0' || $ref_chain->[$i] eq '-'|| $ref_chain->[$i] eq $ref_chain->[$marker_num]) {
			$length++;
		}else{
			last;
		}
	}
	return $length;
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
	-window_type参数越大，-diff_ratio比值越低，纠补越粗暴。
	
Usage:
  Options:
  -m		<file>	Map file, forced
  -l		<file>	Loc file with linkage phase info,forced
  -fix		<file>	Fix_list,hold marker list,not correct

  -k		<str>	Key of output file,forced
  -d		<str>	dir of output file,default ./

  -sus		<float>	suspi ratio,default 0.06
  -imput	<0,1>	1:imput,0:not imput,default 1

  -window_type	<str>	"marker" or "map",default "marker"
	marker:
	  -marker	<int>	extend marker sum,default 6
	map:
	  -fmap		<file>	map file,forced
	  -map		<float>	extend map distance,default 20


  -diff_ratio	<float>	diff ratio,default 0.85
  -cycle_threshold	<int>	cycle times,default 3
  -h		Help

USAGE
	print $usage;
	exit;
}
