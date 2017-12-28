#!/usr/bin/perl 
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
my ($fIn,$fOut,$group,$window,$step);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"g:s"=>\$group,
				"o:s"=>\$fOut,
				"w:s"=>\$window,
				"s:s"=>\$step,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $group);
my %group;
my %rgroup;
$window||=20000;
$step||=5000;
$rgroup{Total}=1;
if ($group){
	open In,$group;
	while (<In>){
		chomp;
		next if ($_ eq ""||/^$/ || /^#/);
		my ($id,$gr,undef)=split(/\s+/,$_);
		$group{$id}=$gr;
		$rgroup{$gr}=1;
	}
	close In;
}
open In,$fIn;
open Out,">$fOut";
my @Indi;
my %Geno;
my %Alle;
while (<In>){
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/){
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
		push @Indi,@indi;
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@Geno)=split(/\t/,$_);
		$pos=join("\t",$chr,$pos);
		my @format=split(/\:/,$format);
		my @All=split(",",join(",",$ref,$alt));
		for (my $i=0;$i<@Geno;$i++){
			my @geno=split(/\:/,$Geno[$i]);
			my $gid="Total";
			if (exists $group{$Indi[$i]}){
				$gid=$group{$Indi[$i]};
			}
			for (my $j=0;$j<@geno;$j++){
				if ($format[$j] eq "GT"){
					my @alle=split(/\//,$geno[$j]);
					$Alle{Total}{$pos}{$alle[0]}++;
					$Alle{Total}{$pos}{$alle[1]}++;
					$Alle{$gid}{$pos}{$alle[0]}++;
					$Alle{$gid}{$pos}{$alle[1]}++;
					$Geno{Total}{$pos}{$geno[$j]}++;
					$Geno{$gid}{$pos}{$geno[$j]}++;
				}	
			}
		}
	}
}
close In;
my $shanon=shanon_calc(\%Alle);
my $pic=pic_calc(\%Alle);
my $home=home_calc(\%Alle);
my @shanon=split(/\s+/,$shanon);
my @pic=split(/\s+/,$pic);
my @home=split(/\s+/,$home);
my @out;
my $n=0;
for (my $i=0;$i<@shanon;$i++) {
	push @{$out[$i]},$shanon[$i];
	push @{$out[$i]},$pic[$i];
	push @{$out[$i]},$home[$i];
}
my @Out;
for (my $i=0;$i<@out;$i++) {
	push @Out,join("\:",@{$out[$i]});
}
print Out join("\t",sort keys %rgroup),"\n";
print Out join("\t","shanon:pic:hete:homo",@Out),"\n";
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub home_calc{#Average
	my ($Geno)=@_;
	my @return;
	foreach my $pop (sort keys %$Geno){
		my $sumH=0;
		my $sumJ=0;
		my $nPos=0;
		foreach my $pos (sort keys %{$$Geno{$pop}}) {
			my @alle=sort{$$Geno{$pop}{$pos}{$b} <=> $$Geno{$pop}{$pos}{$a}} keys %{$$Geno{$pop}{$pos}};
			my @value=sort{$b<=>$a} values %{$$Geno{$pop}{$pos}};
			my $J=0;
			my $H=1;
			for (my $i=0;$i<@value;$i++){
				my $p=$value[$i]/sum(\@value);
				$J+=$p*$p;
			}
			$H=1-$J;
			$sumH+=$H;
			$sumJ+=$J;
			$nPos++;
		}
		my $AH=$sumH/$nPos;
		my $AJ=$sumH/$nPos;
		push @return,"$AH:$AJ";
	}
	return join("\t",@return);
}
sub pic_calc{
	my ($Geno)=@_;
	my @return;
	foreach my $pop(sort keys %$Geno){
		my $sumPIC;
		my $nPos=0;
		foreach my $pos (sort keys %{$$Geno{$pop}}) {
			my @alle=sort{$$Geno{$pop}{$pos}{$b} <=> $$Geno{$pop}{$pos}{$a}} keys %{$$Geno{$pop}{$pos}};
			my @value=sort{$b<=>$a} values %{$$Geno{$pop}{$pos}};
			my $PIC=1;
			for(my $i=0;$i<@value;$i++){
				for(my $j=1;$j<@value;$j++){
					my $p=$value[$i]/sum(\@value);
					my $q=$value[$j]/sum(\@value);
					$PIC-=$p*$p*$q*$q*2;
				}
				my $p=$value[$i]/sum(\@value);
				$PIC-=$p*$p;
			}
			$sumPIC+=$PIC;
			$nPos++;
		}
		push @return,$sumPIC/$nPos;
	}
	return join("\t",@return);
}
sub shanon_calc{
	my ($Geno)=@_;
	my @return;
	foreach my $pop(sort keys %$Geno){
		my $sumH=0;
		foreach my $pos (sort keys %{$$Geno{$pop}}) {
			my @alle=sort{$$Geno{$pop}{$pos}{$b} <=> $$Geno{$pop}{$pos}{$a}} keys %{$$Geno{$pop}{$pos}};
			my @value=sort{$b<=>$a} values %{$$Geno{$pop}{$pos}}; 
			my $H=0;
			for (my $i=0;$i<@value;$i++){
				$H+=-1* $value[$i]/sum(\@value) * log($value[$i]/sum(\@value));
			}
			$sumH+=$H;
		}
		push @return,$sumH;
	}
	return join("\t",@return);
}
sub maf_calc{
	my $Geno=@_;
	my @maf;
	foreach my $pop(sort keys %$Geno){
		my @alle=sort{$$Geno{$pop}{$b} <=> $$Geno{$pop}{$a}} keys %{$$Geno{$pop}};
		my @value=sort{$b<=>$a} values %{$$Geno{$pop}}; 
		my $info;
		if (scalar @alle ==0){
			$info=join(":","-","-");
		}else{
			$info=join(":",join(",",@alle),join(",",@value));
		}
		if(scalar @alle == 1){
			push @maf,$info.":0";
		}elsif(scalar @alle == 0){
			push @maf,$info.":-";
		}else{
			push @maf,$info.":".$alle[1]/sum(@value);
		}
	}
	return join("\t",@maf);
}
sub hw_calc{
	my ($Alle,$Geno)=@_;
	my @HW;
	my @return;
	foreach my $pop(sort keys %$Geno){
		my @alle=sort{$$Alle{$pop}{$b} <=> $$Alle{$pop}{$a}} keys %{$$Alle{$pop}};
		my @value=sort{$b<=>$a} values %{$$Alle{$pop}};
		my %HWstat;
		for (my $i=0;$i<@alle;$i++){#now only for two
			for (my $j=1;$j<@alle;$j++){
				my $p=$$Alle{$pop}{$alle[$i]}/sum(@value);
				my $q=$$Alle{$pop}{$alle[$j]}/sum(@value);
				my $geno=join("/",sort($alle[$i],$alle[$j]));
				$HWstat{ob}{$geno}=$$Geno{$pop}{$geno}/sum(values %{$$Geno{$pop}});
				$HWstat{ex}{$geno}=$p*$q*2;
			}
			my $geno=join("/",$alle[$i],$alle[$i]);
			$HWstat{ob}{$geno}=$$Geno{$pop}{$geno}/sum(values %{$$Geno{$pop}});
			my $p=$$Alle{$pop}{$alle[$i]}/sum(@value);
			$HWstat{ex}{$geno}=$p*$p;
		}
		my @order=sort {$HWstat{ob}{$a}<=>$HWstat{ob}{$b}}keys %{$HWstat{ob}};
		my @od;my @ex;
		for(my $i=0;$i<@order;$i++){
			push @od,$HWstat{ob}{$order[$i]};
			push @ex,$HWstat{ex}{$order[$i]};
		}#
		my $od=join(":",@od);
		my $ex=join(":",@ex);
		my $seg=Segregation($od,$ex,1);
		push @return,$seg;
	}
	return join("\t",@return);

}
sub sum{
	my ($a)=@_;
	my $sum=0;
	foreach my $n(@$a){
		$sum+=$n;
	}
	return $sum;
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

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang on development

Usage:
  Options:
	-i	<file>	input vcf file 
	-g	<file>	input group file
	-o	<file>	output file
	-w	<num>	window size default 20k
	-s	<step>	step size default 5k

USAGE
	print $usage;
	exit;
}

