#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp , $indel , $anno , $fOut );
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"snp:s"=>\$snp,
	"indel:s"=>\$indel,
	"anno:s"=>\$anno,
	"out:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($snp and $indel and $anno and $fOut );
open In,$snp;
my %info;
my %stat;
my %eff;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	my ($gname,$gid,$tid,$biotype,$high,$low,$moderate,$modifier,undef)=split(/\t/,$_);
	my $info=join("\t",$gname,$gid,$tid,$biotype);
	$info{$gname}=$info;
	$info{$gid}=$info;
	$info{$tid}=$info;
	$eff{$info}++ if($high+$moderate > 0);
	$stat{$info}{high}+=$high;
	$stat{$info}{low}+=$low;
	$stat{$info}{middle}+=$moderate;
	$stat{$info}{unknow}+=$modifier;
}
close In;
open In,$indel;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	my ($gname,$gid,$tid,$biotype,$high,$low,$moderate,$modifier,undef)=split(/\t/,$_);
	my $info=join("\t",$gname,$gid,$tid,$biotype);
	$info{$gname}=$info;
	$info{$gid}=$info;
	$info{$tid}=$info;
	$eff{$info}++ if($high+$moderate > 0);
	$stat{$info}{high}+=$high;
	$stat{$info}{low}+=$low;
	$stat{$info}{middle}+=$moderate;
	$stat{$info}{unknow}+=$modifier;
}
close In;
open In,$anno;
open Out,">$fOut.summary";
my %kdetail;
my %gdetail;
my %enrich;
my %edetail;
my %fun;
print Out join("\t","#Gene_name","Gene_id","Transcript_id","Bio_Type","Chr","Pos1","Pos2","High","Moderate","Low","Modifier","NR-ID","NR-ANNO","Uni-ID","Uni-ANNO","KEGG-ID","KEGG-ANNO","GO-ID","GO-ANNO","EggNOG-ID","EggNOG-ANNO"),"\n";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		next;	
	}else{
		my ($id,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$expre)=split(/\t/,$_);
		my @ids=split(/:/,$id);
		my $pos=join("\t",$ids[2],$ids[3],$ids[4]);
		$info{$ids[0]}||=join("\t",$ids[1],$ids[0],"--","--");
		my $info=$info{$ids[0]};
		$stat{$info}{high}||=0;
		$stat{$info}{low}||=0;
		$stat{$info}{middle}||=0;
		$stat{$info}{unknow}||=0;
		$fun{total}{nr}++ if($nrid ne "--");
		$fun{eff}{nr}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $nrid ne "--");
		$fun{total}{uni}++ if($uniid ne "--");
		$fun{eff}{uni}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $uniid ne "--");
		$fun{total}{kegg}++ if($koid ne "--");
		$fun{eff}{kegg}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $koid ne "--");
		$fun{total}{go}++ if($goid ne "--");
		$fun{eff}{go}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $goid ne "--");
		$fun{total}{eggnog}++ if($eid ne "--");
		$fun{eff}{eggnog}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $eid ne "--");
		print Out join("\t",$info,$pos,$stat{$info}{high},$stat{$info}{middle},$stat{$info}{low},$stat{$info}{unknow},$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$eid,$eanno),"\n";
		my @koid=split(/:/,$koid);
		my @kdetail=split(/:/,$koanno);
		for (my $i=0;$i<@koid;$i++) {
			next if ($koid[$i] eq "--");
			$enrich{$koid[$i]}{total}++;
			$enrich{$koid[$i]}{enrich}++ if($stat{$info}{high}+$stat{$info}{middle} > 0);
			$kdetail{$koid[$i]}=$kdetail[$i];
		}
		my @goid=split(/,/,$goid);
		my @gdetail=split(/:/,$goanno);
		for (my $i=0;$i<@goid;$i++) {
			next if ($goid[$i] eq "--");
			$enrich{$goid[$i]}{total}++;
			$enrich{$goid[$i]}{enrich}++ if($stat{$info}{high}+$stat{$info}{middle} > 0);
			$gdetail{$goid[$i]}=$gdetail[$i];
		}
		my ($eid,$eanno)=split(/:/,$egid);
		my @eid=split(/,/,$eid);
		my @edetail=split(/;/,$eanno);
		for (my $i=0;$i<@eid;$i++) {
			next if ($eid[$i] eq "--");
			$enrich{$eid[$i]}{total}++;
			$enrich{$eid[$i]}{enrich}++ if($stat{$info}{high}+$stat{$info}{middle} > 0);
			$edetail{$eid[$i]}=$edetail[$i];
		}
	}
}
close Out;
close In;
open Out,">$fOut.kegg.stat";
print Out join("\t","#koid","ko_detail","eff_variant","all_gene","total_eff","total_gene"),"\n";
foreach my $koid (sort keys %kdetail) {
	$enrich{$koid}{enrich}||=0;
	$enrich{$koid}{total}||=0;
	print Out join("\t",$koid,$kdetail{$koid},$enrich{$koid}{enrich},$enrich{$koid}{total},scalar keys %eff,scalar keys %stat),"\n";
}
close Out;
open Out,">$fOut.go.stat";
print Out join("\t","#goid","go_detail","eff_variant","all_gene","total_eff","total_gene"),"\n";
foreach my $goid (sort keys %gdetail) {
	$enrich{$goid}{enrich}||=0;
	$enrich{$goid}{total}||=0;
	print Out join("\t",$goid,$gdetail{$goid},$enrich{$goid}{enrich},$enrich{$goid}{total},scalar keys %eff,scalar keys %stat),"\n";
}
close Out;
open Out,">$fOut.eggnog.stat";
print Out join("\t","#eggnogid","eggnog_detail","eff_variant","all_gene","total_eff","total_gene"),"\n";
foreach my $eggnog (sort keys %edetail) {
	$enrich{$eggnog}{enrich}||=0;
	$enrich{$eggnog}{total}||=0;
	print Out join("\t",$eggnog,$edetail{$eggnog},$enrich{$eggnog}{enrich},$enrich{$eggnog}{total},scalar keys %eff,scalar keys %stat),"\n";
}
close Out;
open Out,">$fOut.stat.csv";
$fun{total}{nr}||=0;
$fun{eff}{nr}||=0;
$fun{total}{uni}||=0;
$fun{eff}{uni}||=0;
$fun{total}{kegg}||=0;
$fun{eff}{kegg}||=0;
$fun{total}{go}||=0;
$fun{eff}{go}||=0;
$fun{total}{eggnog}||=0;
$fun{eff}{eggnog}||=0;
print Out join("\t","#type","NR","Uniprot","KEGG","GO","EGGNOG"),"\n";
print Out join("\t","total",$fun{total}{nr},$fun{total}{uni},$fun{total}{kegg},$fun{total}{go},$fun{total}{eggnog}),"\n";
print Out join("\t","eff",$fun{eff}{nr},$fun{eff}{uni},$fun{eff}{kegg},$fun{eff}{go},$fun{eff}{eggnog}),"\n";
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -snp	<file>	input snp file name
  -indel	<file>	input indel file name
  -anno	<file>	input anno file name
  -out	<key>	output keys of file name
  -h         Help

USAGE
        print $usage;
        exit;
}
