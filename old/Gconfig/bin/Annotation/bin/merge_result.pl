#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$Gene);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"g:s"=>\$Gene,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $Gene);
$/=">";
open In,$Gene;
my @gene;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@line)=split(/\n/,$_);
	push @gene,$id;
}
close In;
$/="\n";
my @anno=glob("$fIn/*.anno");
my %info;
my %head;
foreach my $anno (@anno) {
	my $id=(split(/\./,basename($anno)))[1];
	open In,$anno;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/);
		if (/^#/) {
			my (undef,@Head)=split(/\t/,$_);
			$head{$id}=join("|",@Head);
		}else {
			my ($gene,@info)=split(/\t/,$_);
			$info{$gene}{$id}=join("|",@info);
		}
	}
	close In;
}
open Out,">$fOut";
print Out "#geneid\tchr\tpos1\tpos2\t",$head{nr},"\t",$head{swissprot},"\t",$head{kegg},"\t",$head{go},"\t",$head{eggnog},"\n";
foreach my $gene (@gene) {
	if (!exists $info{$gene}{nr}) {
		my @id=split(/\t/,$head{nr});
		$info{$gene}{nr}="--";
	}
	if (!exists $info{$gene}{swissprot}) {
		my @id=split(/\t/,$head{swissprot});
		$info{$gene}{swissprot}="--";
	}
	if (!exists $info{$gene}{kegg}) {
		my @id=split(/\t/,$head{kegg});
		$info{$gene}{kegg}="--";
	}
	if (!exists $info{$gene}{go}) {
		my @id=split(/\t/,$head{go});
		$info{$gene}{go}="--";
	}
	if (!exists $info{$gene}{eggnog}) {
		my @id=split(/\t/,$head{eggnog});
		$info{$gene}{eggnog}="--";
	}
	my @id=split(/\_/,$gene,4);
	print Out join("\t",$id[3],$id[0],$id[1],$id[2]),"\t",$info{$gene}{nr},"\t",$info{$gene}{swissprot},"\t",$info{$gene}{kegg},"\t",$info{$gene}{go},"\t",$info{$gene}{eggnog},"\n";
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
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
