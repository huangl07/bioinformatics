#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($snp,$indel,$anno,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"snp:s"=>\$snp,
				"indel:s"=>\$indel,
				"anno:s"=>\$anno,
				"out:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($snp and $indel and $anno and $fOut);
open In,$snp;
my @type;
my %info;
my %effect;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#GeneName/) {
		(undef,undef,undef,undef,@type)=split(/\s+/,$_);
		next;
	}
	my ($gname,$gid,$tid,$biotype,@num)=split(/\s+/,$_);
	$info{$gid}=join("\t",$gname,$gid,$tid,$biotype);
	$info{$gname}=join("\t",$gname,$gid,$tid,$biotype);
	$info{$tid}=join("\t",$gname,$gid,$tid,$biotype);
	my $info=join("\t",$gname,$gid,$tid,$biotype);
	$effect{$info}{high}+=$type[0];
	$effect{$info}{low}+=$type[1];
	$effect{$info}{middle}+=$type[1];
	$effect{$info}{modifer}+=$type[2];
}
close In;
open In,$indel;
my @type;
my %info;
my %effect;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#GeneName/) {
		(undef,undef,undef,undef,@type)=split(/\s+/,$_);
		next;
	}
	my ($gname,$gid,$tid,$biotype,@num)=split(/\s+/,$_);
	next if (exists $info{$gid} || exists $info{$gname} || exists $info{$tid});
	$info{$gid}=join("\t",$gname,$gid,$tid,$biotype);
	$info{$gname}=join("\t",$gname,$gid,$tid,$biotype);
	$info{$tid}=join("\t",$gname,$gid,$tid,$biotype);
	my $info=join("\t",$gname,$gid,$tid,$biotype);
	$effect{$info}{high}+=$type[0];
	$effect{$info}{low}+=$type[1];
	$effect{$info}{middle}+=$type[1];
	$effect{$info}{modifer}+=$type[2];
}
close In;
open In,$anno;
open Out,">$fOut.stat";
my %ko;
my %go;
my %eggnog;
my $Total=0;
my $Erich=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$NR,$UNI,$KEGG,$GO,$EGGNOG)=split(/\s+/,$_);
	my $info=$info{$id};
	print Out join("\t",$info{$id},$effect{$info}{high},$effect{$info}{middle},$effect{$info}{low},$effect{$info}{modifer},$NR,$UNIPROT,$KEGG,$GO,$EGGNOG),"\n";
	$Total++;
	$Erich++ if($effect{$info}{high}+$effect{$info}{middle} > 0);
	my ($ko,$kdetail)=split(/\|/,$KEGG);
	my ($GO,$gdetail)=split(/\|/,$GO);
	my ($EGGNOG,$edetail)=split(/\|/,$EGGNOG);
	my @ko=split(/\:/,$ko);
	my @kdetail=split(/\:/,$kdetail);
	for (my $i=0;$i<@ko;$i++) {
		$ko{$ko[$i]}{tnum}++;
		$ko{$ko[$i]}{enum}++ if($effect{$info}{high}+$effect{$info}{middle} > 0);
		$ko{$ko[$i]}{detail}=$kdetail[$i];
	}
	my @go=split(/\,/,$GO);
	my @gdetail=split(/\:/,$gdetail)
	for (my $i=0;$i<@go;$i++) {
		$go{$go[$i]}{tnum}++;
		$go{$go[$i]}{enum}++ if($effect{$info}{high}+$effect{$info}{middle} > 0);
		$go{$go[$i]}{detail}=$gdetail[$i];
	}
	my @egg=split(/\,/,$eggnog);
	my @edetail=split(/\:/,$eggnog)
	for (my $i=0;$i<@egg;$i++) {
		$egg{$egg[$i]}{tnum}++;
		$egg{$egg[$i]}{enum}++ if($effect{$info}{high}+$effect{$info}{middle} > 0);
		$egg{$egg[$i]}{detail}=$edetail[$i];
	}
}
close In;
close Out;
open Out,">$fOut.GO.enrich";
foreach my $go (sort keys %GO) {
	print Out $go{$go}{tnum},"\t",$go{$go}{enum},"\t",$Total,"\t",$Enrich,"\n";
}
close Out;
open Out,">$fOut.kegg.enrich";
foreach my $go (sort keys %ko) {
	print Out $ko{$go}{tnum},"\t",$ko{$go}{enum},"\t",$Total,"\t",$Enrich,"\n";
}

close Out;
open Out,">$fOut.eggnog.enrich";
foreach my $go (sort keys %egg) {
	print Out $egg{$go}{tnum},"\t",$egg{$go}{enum},"\t",$Total,"\t",$Enrich,"\n";
}
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-snp	<file>	input snp gene txt file 
	-indel	<file>	input indel gene txt file
	-anno	<file>	anno file
	-out	<file>	output file
USAGE
	print $usage;
	exit;
}
