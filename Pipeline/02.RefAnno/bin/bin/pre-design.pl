#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$ref,
	"o:s"=>\$fOut,
	) or &USAGE;
&USAGE unless ($ref and $fOut);
my %enzyme=(
	"EcoRI"=>"GAATTC",
	"MseI"=>"TTAA",
	"PstI"=>"CTGCAG",
	"TaqaI"=>"TCGA",
);
my @enzyme=keys %enzyme;
my %range=(
	"330-380"=>1,
	"380-430"=>1,
	"430-480"=>1,
	"480-530"=>1,
	"530-580"=>1,
);
open In,$ref;
$/=">";
my %RAD;
my %GBS;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,@seq)=split(/\n/,$_);
	my $seq=join("",@seq);
	for (my $i=0;$i<@enzyme;$i++) {
		my %pos;
		my $cut=$enzyme{$enzyme[$i]};
		while ($seq =~ m/($cut)/g) {
			my $pos=pos($seq)-length($cut);
			$pos{$pos}=1;
		}
		my @pos=sort {$a<=>$b} keys %pos;
		for (my $j=0;$j<@pos-1;$j++) {
			my $l=$pos[$j+1]-$pos[$j];
			$RAD{$enzyme[$i]}{"330-380"}++ if($l > 330);
			$RAD{$enzyme[$i]}{"380-430"}++ if($l > 380);
			$RAD{$enzyme[$i]}{"430-480"}++ if($l > 430);
			$RAD{$enzyme[$i]}{"480-530"}++ if($l > 480);
			$RAD{$enzyme[$i]}{"530-580"}++ if($l > 530);

		}
	}
	for (my $i=0;$i<@enzyme;$i++) {
		for (my $j=$i+1;$j<@enzyme;$j++) {
			my $cut1=$enzyme{$enzyme[$i]};
			my $cut2=$enzyme{$enzyme[$j]};
			my %pos;
			while ($seq =~ m/($cut1|$cut2)/g) {
				my $pos=pos($seq)-length($cut1);
				$pos{$pos}=1;
			}
			my @pos=sort {$a<=>$b} keys %pos;
			my $str=join("\t",$enzyme[$i],$enzyme[$j]);
			for (my $k=0;$k<@pos-1;$k++) {
				my $len=$pos[$k+1]-$pos[$k];
				$GBS{$str}{"330-380"}++ if ($len >= 330 and $len <=380);
				$GBS{$str}{"380-430"}++ if ($len >= 380 and $len <=430);
				$GBS{$str}{"430-480"}++ if ($len >= 430 and $len <=480);
				$GBS{$str}{"480-530"}++ if ($len >= 480 and $len <=530);
				$GBS{$str}{"530-580"}++ if ($len >= 530 and $len <=580);
			}
		}
	}
}
close In;
open Out,">$fOut";
my @eGBS=keys %GBS;
my @eRAD=keys %RAD;
print Out join("\t","#enzyme",join("\t",@eGBS),join("\t",@eRAD)),"\n";
foreach my $range (keys %range) {
	my @out;
	push @out,$range;
	foreach my $e (@eGBS) {
		print $GBS{$e}{$range};die;
		push @out,$GBS{$e}{$range};
	}
	foreach my $e (@eRAD) {
		push @out,$RAD{$e}{$range};
	}
	print Out join("\t",@out),"\n";
}
close Out;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	reformat genome,rename scaffold name at genome fa file and gff file
eg:
	perl $Script -i Genome.fa -g Genome.gff -k keyname -o dir 

Usage:
  Options:
  -i	<file>	input genome name,fasta format,
  -o	<str>	output file prefix

  -h         Help

USAGE
        print $usage;
        exit;
}
