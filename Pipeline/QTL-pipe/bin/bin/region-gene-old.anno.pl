#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gff , $anno , $fOut ,$select);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$select,
	"a:s"=>\$anno,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($select and $anno and $fOut );
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:

Usage:
  Options:
  -i	<file>  input select file name
  -a	<file>  input anno file name
  -o	<key>   output keys of file name
  -h         Help

USAGE
        print $usage;
        exit;
}

open In,$select;
my %region;
my $qtl;
while (<In>) {
	chomp;
	s/\"//g;
	next if ($_ eq ""||/^$/ ||/#/ ||/mark1/);
	my ($marker,$chr,$pos,$lod,$var,$pm1,$pm2,$start,$end,$mark1,$mark2)=split(/\s+/,$_);
	my ($chr1,$start1)=split(/\_/,$mark1);
	my ($chr2,$start2)=split(/\_/,$mark2);
	if ($chr1 ne $chr2) {
		die "error region!,$chr1\t$chr2\n";
	}else{
		#$number=(split(/\D+/,$chr))[-1];
		my $number=$chr;
		$region{$number}{$qtl}=join("\t",sort{$a<=>$b}($start1,$start2));
	}
}
close In;

open In,$anno;
my %kdetail;
my %gdetail;
my %enrich;
my %edetail;
my %info;
my %stat;
my %astat;
my $head;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		$head=$_;
		next;
	}
	#my ($Gene_name,$Gene_id,$Transcript_id,$Bio_Type,$chr,$Pos1,$Pos2,$High,$Moderate,$Low,$Modifier,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$express)=split(/\t/,$_);
	my ($Gene_name,$Gene_id,$Transcript_id,$Bio_Type,$High,$Moderate,$Low,$Modifier,$chr,$Pos1,$Pos2,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$express)=split(/\s+/,$_);
	my $regioned=0;
	my $chrid=(split(/\D+/,$chr))[-1];
		#foreach my $region (sort keys %{$region{$chr}}) {
	foreach my $number (sort{$a<=>$b} keys %region){
		my ($pos3,$pos4)=split(/\t/,$region{$number}{$qtl});
		if($chrid eq $number){
			if (($Pos1 > $pos3 && $Pos1 <$pos4) ||($Pos2 > $pos3 && $Pos2 < $pos4) || ($pos3 > $Pos1 && $pos3 < $Pos2) || ($pos4 > $Pos1 && $pos4 < $Pos2)) {
				#push @{$info{$chr}{$region}},$_;
				push @{$info{$number}{$qtl}},$_;
				$stat{$number}{$qtl}{total}++;
				$stat{$number}{$qtl}{totalnr}++ if($nrid ne "--");
				$stat{$number}{$qtl}{totaluni}++ if($uniid ne "--");
				$stat{$number}{$qtl}{totalkegg}++ if($koid ne "--");
				$stat{$number}{$qtl}{totalgo}++ if($goid ne "--");
				$stat{$number}{$qtl}{totaleggnog}++ if($egid ne "S");
				$stat{$number}{$qtl}{eff}++ if($High+$Moderate > 0);
				$stat{$number}{$qtl}{effnr}++ if($nrid ne "--" && $High+$Moderate > 0);
				$stat{$number}{$qtl}{effuni}++ if($uniid ne "--" && $High+$Moderate > 0 );
				$stat{$number}{$qtl}{effkegg}++ if($koid ne "--" && $High+$Moderate > 0);
				$stat{$number}{$qtl}{effgo}++ if($goid ne "--" && $High+$Moderate > 0);
				$stat{$number}{$qtl}{effeggnog}++ if($egid ne "S" && $High+$Moderate > 0);
				$regioned=1;
			}
		}
	}
	$astat{total}++;
	$astat{eff}++ if($regioned == 1);
	my @koid=split(/:/,$koid);
	my @kdetail=split(/:/,$koanno);
	for (my $i=0;$i<@koid;$i++) {
		next if ($koid[$i] eq "--");
		$enrich{$koid[$i]}{total}++;
		$enrich{$koid[$i]}{enrich}++ if($regioned ==1);
		$kdetail{$koid[$i]}=$kdetail[$i];
	}
	my @goid=split(/,/,$goid);
	my @gdetail=split(/:/,$goanno);
	for (my $i=0;$i<@goid;$i++) {
		next if ($goid[$i] eq "--");
		$enrich{$goid[$i]}{total}++;
		$enrich{$goid[$i]}{enrich}++ if($regioned ==1);
		$gdetail{$goid[$i]}=$gdetail[$i];
	}
	my ($eid,$eanno)=split(/:/,$egid);
	my @eid=split(/,/,$eid);
	my @edetail=split(/;/,$eanno);
	for (my $i=0;$i<@eid;$i++) {
		next if ($eid[$i] eq "--");
		$enrich{$eid[$i]}{total}++;
		$enrich{$eid[$i]}{enrich}++ if($regioned ==1);
		$edetail{$eid[$i]}=$edetail[$i];
	}
}
close In;
open Out,">$fOut.total";
print Out "#\@chr\tpos1\tpos2\ttotal\teff\ttotalnr\ttotaluni\ttotalkegg\ttotalgo\ttotaleggnog\teffnr\teffuni\teffkegg\teffgo\teffeggnog\n";
print Out "$head\n";
foreach my $number (sort{$a<=>$b} keys %region) {
	foreach my $qtl (sort keys %{$stat{$number}}) {
		$stat{$number}{$qtl}{totalnr}||=0;
		$stat{$number}{$qtl}{totaluni}||=0;
		$stat{$number}{$qtl}{totalkegg}||=0;
		$stat{$number}{$qtl}{totalgo}||=0;
		$stat{$number}{$qtl}{totaleggnog}||=0;
		$stat{$number}{$qtl}{effnr}||=0;
		$stat{$number}{$qtl}{effuni}||=0;
		$stat{$number}{$qtl}{effkegg}||=0;
		$stat{$number}{$qtl}{effgo}||=0;
		$stat{$number}{$qtl}{effeggnog}||=0;
		$stat{$number}{$qtl}{total}||=0;
		$stat{$number}{$qtl}{eff}||=0;
		#print Out join("\t","\@$chr",$region,$stat{$chr}{$region}{total},$stat{$chr}{$region}{eff},$stat{$chr}{$region}{totalnr},$stat{$chr}{$region}{totaluni},$stat{$chr}{$region}{totalkegg},$stat{$chr}{$region}{totalgo},$stat{$chr}{$region}{totaleggnog},$stat{$chr}{$region}{effnr},$stat{$chr}{$region}{effuni},$stat{$chr}{$region}{effkegg},$stat{$chr}{$region}{effgo},$stat{$chr}{$region}{effeggnog}),"\n";
		print Out join("\t","\@chr$number",$qtl,$stat{$number}{$qtl}{total},$stat{$number}{$qtl}{eff},$stat{$number}{$qtl}{totalnr},$stat{$number}{$qtl}{totaluni},$stat{$number}{$qtl}{totalkegg},$stat{$number}{$qtl}{totalgo},$stat{$number}{$qtl}{totaleggnog},$stat{$number}{$qtl}{effnr},$stat{$number}{$qtl}{effuni},$stat{$number}{$qtl}{effkegg},$stat{$number}{$qtl}{effgo},$stat{$number}{$qtl}{effeggnog}),"\n";
		next if (!defined $info{$number}{$qtl});
		print Out join("\n",@{$info{$number}{$qtl}}),"\n";
	}
}
close Out;
open Out,">$fOut.eff";
print Out "#\@chr\tpos1\tpos2\n";
print Out "$head\n";
foreach my $number (sort{$a<=>$b} keys %region) {
	foreach my $qtl (sort keys %{$region{$number}}) {
		print Out "\@chr$number\t$qtl\n";
		next if (!defined $info{$number}{$qtl});
		foreach my $line (@{$info{$number}{$qtl}}) {
			my @split=split(/\t/,$line);
			print Out $line,"\n" if($split[7]+$split[8] > 0);
		}
	}
}
close Out;

open Out,">$fOut.kegg.stat";
foreach my $koid (sort keys %kdetail) {
	$enrich{$koid}{enrich}||=0;
	$enrich{$koid}{total}||=0;
	$astat{totalkegg}||=0;
	$astat{effkegg}||=0;
	print Out join("\t",$koid,$kdetail{$koid},$enrich{$koid}{enrich},$enrich{$koid}{total},$astat{eff},$astat{total}),"\n";
}
close Out;
open Out,">$fOut.go.stat";
foreach my $goid (sort keys %gdetail) {
	$enrich{$goid}{enrich}||=0;
	$enrich{$goid}{total}||=0;
	$astat{totalgo}||=0;
	$astat{effgo}||=0;
	print Out join("\t",$goid,$gdetail{$goid},$enrich{$goid}{enrich},$enrich{$goid}{total},$astat{eff},$astat{total}),"\n";
}
close Out;
open Out,">$fOut.eggnog.stat";
foreach my $eggnog (sort keys %edetail) {
	$enrich{$eggnog}{enrich}||=0;
	$enrich{$eggnog}{total}||=0;
	$astat{totaleggnog}||=0;
	$astat{effeggnog}||=0;
	print Out join("\t",$eggnog,$edetail{$eggnog},$enrich{$eggnog}{enrich},$enrich{$eggnog}{total},$astat{eff},$astat{total}),"\n";
}
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
