#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$Key,$enz1);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
my $enz2 = "";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"e1:s"=>\$enz1,
	"e2:s"=>\$enz2,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
			) or &USAGE;
&USAGE unless ($fIn and $dOut and $enz1 and $Key);
my $mkdir=1;
$mkdir=(mkdir $dOut) if (!-d $dOut);
die "Error make dir $dOut" if($mkdir == 0);
die "Error input Genome $fIn!\n" if (!-f $fIn );
my %enzyme=(
	"EcoRI"=>"GAATTC\t1",
	"MseI"=>"TTAA\t1",
	"PstI"=>"CTGCAG\t1",
	"TaqaI"=>"TCGA\t1",
);
$enz2=$enz1 if ($enz2 eq "");
die "Error input enzyme name!$enz1" if (!exists $enzyme{$enz1});
die "Error input enzyme name!$enz2" if (!exists $enzyme{$enz2});

my %cut;
my ($cut1,$len1)=split(/\s+/,$enzyme{$enz1});
my ($cut2,$len2)=split(/\s+/,$enzyme{$enz2});
$cut{$cut1}{len}=$len1;
$cut{$cut2}{len}=$len2;
$cut{$cut1}{nam}=$enz1;
$cut{$cut2}{nam}=$enz2;
my @cut=keys %cut;
open In,$fIn;
open Detail,">$dOut/$Key.enzyme.detail";
#open Fa,">$dOut/$Key.enzyme.fa";
#open Fa1,">$dOut/$Key.enzyme.fa1";
#open Fa2,">$dOut/$Key.enzyme.fa2";
$/=">";
my %RAD;
my %GBS;
my %len;
my $RADtotal=0;
my $Totallen=0;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($id,@line)=split(/\n/,$_,);
	my $seq=join("",@line);
	my %pos;
	$Totallen+=length($seq);
	foreach my $cut (@cut) {
		while ($seq =~ m/($cut)/g) {
			my $pos=pos($seq)-length($cut)+$cut{$cut}{len};
			$pos{$pos}=$cut{$cut}{nam};
		}
	}
	my @pos=sort{$a<=>$b}keys %pos;
	 if(scalar @pos ==0){print  "error genome,$id,no enzyme site\n";next;}
	my $fragment=substr($seq,0,$pos[0]);
	my $len=$pos[0]-0;
	$RAD{$len}+=1;
	$RADtotal+=1;
	$len{$len}=1;
	if ($len > 200) {
		print Detail $id,"\t","0","\t",$pos[0],"\t",$pos[0]-0,"\t","N-$pos{$pos[0]}\n";
		#print Fa ">$id\_0\_$pos[0]\n".$fragment,"\n";
		my $read1=substr($fragment,0,150);
		my $read2=substr($fragment,-150,150);
		$read2=reverse($read2);
		$read2=~tr/ATGC/TACG/;
		#print Fa1 ">$id\_$pos[-1]\_".length($fragment)."\_1\n$read1\n";
		#print Fa2 ">$id\_$pos[-1]\_".length($fragment)."\_2\n$read2\n";
	}
	for (my $i=0;$i<@pos-1;$i++) {
		my $fragment=substr($seq,$pos[$i],$pos[$i+1]-$pos[$i]);
		$len=$pos[$i+1]-$pos[$i];
		next if ($len < 200);
		print Detail $id,"\t",$pos[$i],"\t",$pos[$i+1],"\t",$pos[$i+1]-$pos[$i],"\t","$pos{$pos[$i]}-$pos{$pos[$i+1]}","\n";
		#print Fa ">$id\_$pos[$i]\_$pos[$i+1]\n".$fragment,"\n";
		if ($enz1 eq $enz2) {
			$RAD{$len}+=2;
			$RADtotal+=2;
			$len{$len}=1;
		}else{
			if ($pos{$pos[$i]} ne $pos{$pos[$i+1]}) {
				$GBS{$len}++;
				$len{$len}=1;
			}
		}
		my $read1=substr($fragment,0,150);
		my $read2=substr($fragment,-150,150);
		$read2=reverse($read2);
		$read2=~tr/ATGC/TACG/;
		#print Fa1 ">$id\_$pos[-1]\_".length($fragment)."\_1\n$read1\n";
		#print Fa2 ">$id\_$pos[-1]\_".length($fragment)."\_2\n$read2\n";
	}
	$len=length($seq)-$pos[-1];
	$fragment=substr($seq,$pos[-1],length($seq)-$pos[-1]);
	if ($len > 200) {
		print Detail $id,"\t","0","\t",$pos[0],"\t",$pos[0]-0,"\t","N-$pos{$pos[0]}\n";
		#print Fa ">$id\_0\_$pos[0]\n".$fragment,"\n";
		my $read1=substr($fragment,0,150);
		my $read2=substr($fragment,-150,150);
		$read2=reverse($read2);
		$read2=~tr/ATGC/TACG/;
		#print Fa1 ">$id\_$pos[-1]\_".length($fragment)."\_1\n$read1\n";
		#print Fa2 ">$id\_$pos[-1]\_".length($fragment)."\_2\n$read2\n";
	}
	$RAD{$len}+=1;
	$len{$len}=1;
	$RADtotal+=1;

}
close In;
close Detail;
#close Fa;
#close Fa1;
#close Fa2;
open Stat,">$dOut/$Key.pre-design.stat";
my %sumRAD;
my %sumGBS;
my %covGBS;
my %covRAD;
open Out,">$dOut/$Key.enzyme.draw";
print Out "#lenth\tnumber\n";
foreach my $len (sort {$a<=>$b} keys %len) {
	if ($enz1 eq $enz2) {
		$RAD{$len}||=0;
		print Out $len,"\t",$RAD{$len},"\n";
	}else{
		$GBS{$len}||=0;
		print Out $len,"\t",$GBS{$len},"\n";
	}
	if ($enz1 ne $enz2) {
		$GBS{$len}||=0;
		if ($len >= 280 and $len <=330) {
			$sumGBS{"280-330"}+=$GBS{$len};
			if ($len <= 300) {
				$covGBS{"280-330"}+=$GBS{$len}*$len/$Totallen;
			}else{
				$covGBS{"280-330"}+=$GBS{$len}*300/$Totallen;
			}
		}elsif ($len >= 330 and $len <=380) {
			$sumGBS{"330-380"}+=$GBS{$len};
			$covGBS{"330-380"}+=$GBS{$len}*300/$Totallen;
		}elsif ($len >= 380 and $len <=430) {
			$sumGBS{"380-430"}+=$GBS{$len};
			$covGBS{"380-430"}+=$GBS{$len}*300/$Totallen;
		}elsif ($len >= 430 and $len <=480) {
			$sumGBS{"430-480"}+=$GBS{$len};
			$covGBS{"430-480"}+=$GBS{$len}*300/$Totallen;
		}elsif ($len >= 480 and $len <=530) {
			$sumGBS{"480-530"}+=$GBS{$len};
			$covGBS{"480-530"}+=$GBS{$len}*300/$Totallen;
		}elsif($len >= 530 and $len <=580) {
			$sumGBS{"530-580"}+=$GBS{$len};
			$covGBS{"530-580"}+=$GBS{$len}*300/$Totallen;
		}
	}else{
		$RAD{$len}||=0;
		if ($len >= 530) {
			$sumRAD{"280-330"}+=$RAD{$len};
			$sumRAD{"330-380"}+=$RAD{$len};
			$sumRAD{"380-430"}+=$RAD{$len};
			$sumRAD{"430-480"}+=$RAD{$len};
			$sumRAD{"480-530"}+=$RAD{$len};
			$sumRAD{"530-580"}+=$RAD{$len};

			$covRAD{"280-330"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"330-380"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"380-430"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"430-480"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"480-530"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"530-580"}+=$RAD{$len}*300/$Totallen;

		}elsif ($len >= 480) {
			$sumRAD{"280-330"}+=$RAD{$len};
			$sumRAD{"330-380"}+=$RAD{$len};
			$sumRAD{"380-430"}+=$RAD{$len};
			$sumRAD{"430-480"}+=$RAD{$len};
			$sumRAD{"480-530"}+=$RAD{$len};
			$covRAD{"280-330"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"330-380"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"380-430"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"430-480"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"480-530"}+=$RAD{$len}*300/$Totallen;

		}elsif ($len >= 430) {
			$sumRAD{"280-330"}+=$RAD{$len};
			$sumRAD{"330-380"}+=$RAD{$len};
			$sumRAD{"380-430"}+=$RAD{$len};
			$sumRAD{"430-480"}+=$RAD{$len};
			$covRAD{"280-330"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"330-380"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"380-430"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"430-480"}+=$RAD{$len}*300/$Totallen;

		}elsif ($len >= 380) {
			$sumRAD{"280-330"}+=$RAD{$len};
			$sumRAD{"330-380"}+=$RAD{$len};
			$sumRAD{"380-430"}+=$RAD{$len};
			$covRAD{"280-330"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"330-380"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"380-430"}+=$RAD{$len}*300/$Totallen;

		}elsif ($len >= 330) {
			$sumRAD{"280-330"}+=$RAD{$len};
			$sumRAD{"330-380"}+=$RAD{$len};
			$covRAD{"280-330"}+=$RAD{$len}*300/$Totallen;
			$covRAD{"330-380"}+=$RAD{$len}*300/$Totallen;

		}elsif ($len >= 280) {
			$sumRAD{"280-330"}+=$RAD{$len};
			$covRAD{"280-330"}+=$RAD{$len}*300/$Totallen;
		}
	}
}
close Out;

foreach my $value ((keys %covRAD,keys %covGBS)) {
	if ($enz1 eq $enz2){
		$covRAD{$value}=int($covRAD{$value}*10000)/10000;
	}else{
		$covGBS{$value}=int($covGBS{$value}*10000)/10000;
	}
}
print Stat "#Genome\ttype\tenzyme\t280-330\t330-380\t380-430\t430-480\t480-530\t530-580\n";
if ($enz1 eq $enz2){
	print Stat "$Key\tRAD-num\t"."$enz1+$enz2","\t",join("\t",$sumRAD{"280-330"},$sumRAD{"330-380"},$sumRAD{"380-430"},$sumRAD{"430-480"},$sumRAD{"480-530"},$sumRAD{"530-580"}),"\n";
	print Stat "$Key\tRAD-num\t"."$enz1+$enz2","\t",join("\t",$covRAD{"280-330"},$covRAD{"330-380"},$covRAD{"380-430"},$covRAD{"430-480"},$covRAD{"480-530"},$covRAD{"530-580"}),"\n";
}else{
	print Stat "$Key\tGBS-num\t"."$enz1+$enz2","\t",join("\t",$sumGBS{"280-330"},$sumGBS{"330-380"},$sumGBS{"380-430"},$sumGBS{"430-480"},$sumGBS{"480-530"},$sumGBS{"530-580"}),"\n";
	print Stat "$Key\tGBS-cov\t"."$enz1+$enz2","\t",join("\t",$covGBS{"280-330"},$covGBS{"330-380"},$covGBS{"380-430"},$covGBS{"430-480"},$covGBS{"480-530"},$covGBS{"530-580"}),"\n";
}
close Stat;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -e1	<str>	input enzyme1 name
  -e2	<str>	input enzyme2 name
  -o	<dir>	output dir name
  -k	<str>	output keys of filename
  -h         Help

USAGE
        print $usage;
        exit;
}
