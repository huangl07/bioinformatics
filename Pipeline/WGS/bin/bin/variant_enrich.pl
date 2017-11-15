#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$fAnno);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"s:s"=>\$fAnno,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $fAnno);
open In,$fIn;
my %genelist;
my %rich;
my @ID;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/#GeneName/) {
		(undef,@ID)=split(/\t/,$_);
	}
	my ($id,@info)=split(/\t/,$_);
	$genelist{$id}=join("\t",@info[8..$#info]);
}
close In;
open In,$fAnno;
open Out,">$fOut";
my %TotalG;
my %TotalK;
my %TotalE;
my %ListG;
my %ListK;
my %ListE;
print Out "#geneid\tchr\tstart\tend\t",join("\t",@ID[8..$#ID],"NR","UNIPROT","KEGG","GO","EGGNOG"),"\n";
while (<In>) {
	next if ($_ eq ""||/^$/||/^#/);
	my ($id,$chr,$start,$end,$nr,$uni,$KO,$GO,$EGGNOG)=split(/\t/,$_);
	print Out $id,"\t",$chr,"\t",$start,"\t",$end,"\t",$genelist{$id},"\t",join("\t",$nr,$uni,$KO,$GO,$EGGNOG),"\n";
	my ($koid,$KID,$ko,$kdes)=split(/\|/,$KO);
	print $KO;die;
	my @ko=split(/:/,$ko);
	my @kdes=split(/:/,$kdes);
	my ($goid,$godes)=split(/\|/,$GO);
	my @ogoid=split(/\:/,$goid);
	my @goid;
	for (my $i=0;$i<@ogoid;$i++) {
		push @goid,"GO:$ogoid[$i]" if($ogoid[$i] ne "GO");
	}
	my @godes=split(/:/,$godes);
	my ($eid,$edes)=split(/\|/,$EGGNOG);
	my @eid=split(/\:/,$eid);
	my @edes=split(/\:/,$edes);
	for (my $i=0;$i<@ko;$i++) {
		$TotalK{$ko[$i]}{num}++;
		$TotalK{$ko[$i]}{des}=$kdes[$i];
		if (!exists $genelist{$id}) {
			$ListK{$ko[$i]}{num}++;
			$ListK{$ko[$i]}{des}=$kdes[$i];
		}

	}
	for (my $i=0;$i<@goid;$i++) {
		$TotalG{$goid[$i]}{num}++;
		$TotalG{$goid[$i]}{des}=$godes[$i];
		if (!exists $genelist{$id}) {
			$ListG{$goid[$i]}{num}++;
			$ListG{$goid[$i]}{des}=$kdes[$i];
		}
	}
	for (my $i=0;$i<@eid;$i++) {
		$TotalE{$eid[$i]}{num}++;
		$TotalE{$eid[$i]}{des}=$edes[$i];
		if (!exists $genelist{$id}) {
			$ListE{$eid[$i]}{num}++;
			$ListE{$eid[$i]}{des}=$kdes[$i];
		}
	}

}
close In;
close Out;
open Out,">$fOut.kegg.enrich";
foreach my $ID (sort keys %TotalK) {
	$ListK{$ID}{num}||=0;
	print Out join("\t",$ID,$ListK{$ID}{num},scalar keys %genelist,$TotalK{$ID}{num})
}
close Out;
open Out,">$fOut.go.enrich";
foreach my $ID (sort keys %TotalK) {
	$ListK{$ID}{num}||=0;
	print Out join("\t",$ID,$ListK{$ID}{num},$TotalK{$ID}{num})
}
close Out;
open Out,">$fOut.eggnog.enrich";
foreach my $ID (sort keys %TotalK) {
	$ListK{$ID}{num}||=0;
	print Out join("\t",$ID,$ListK{$ID}{num},$TotalK{$ID}{num})
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
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
