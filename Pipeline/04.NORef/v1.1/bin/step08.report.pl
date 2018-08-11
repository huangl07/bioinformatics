#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcfdir,$dOut,$dShell,$statdir,$fastqc);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcfdir:s"=>\$vcfdir,
	"statdir:s"=>\$statdir,
	"fastqc:s"=>\$fastqc,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($vcfdir and $dOut and $dShell and $statdir and $fastqc);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
$dOut=ABSOLUTE_DIR($dOut);
$vcfdir=ABSOLUTE_DIR($vcfdir);
$statdir=ABSOLUTE_DIR($statdir);
$fastqc=ABSOLUTE_DIR($fastqc);
mkdir "$dOut/Table" if (!-d "$dOut/Table");
my @fqstat=glob("$fastqc/*.stat");
my %fqstat;
foreach my $fqstat (@fqstat) {
	open In,$fqstat;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ ||/#/);
		next if(!/Total/);
		my(undef,$sample,$readnum,$basenum,$a,$t,$g,$c,$n,$GC,$Q30,$Q20,undef)=split(/\s+/,$_);
		my $sampleID=(split(/\-/,$sample))[0];
		$fqstat{$sampleID}{readnum}+=$readnum;
		$fqstat{$sampleID}{basenum}+=$basenum;
		$fqstat{$sampleID}{GC}+=$GC*$basenum;
		$fqstat{$sampleID}{Q30}+=$Q30*$basenum;
	}
	close In;
}
open Out,">$dOut/Table/cleanData.xls";
print Out "#sampleID\tclean Reads\tclean base\tGC(%)\tQ30(%)\n";
foreach my $sample (sort keys %fqstat) {
	my @out;
	push @out,$sample;
	push @out,$fqstat{$sample}{readnum};
	push @out,$fqstat{$sample}{basenum};
	push @out,sprintf("%.2f",$fqstat{$sample}{GC}/$fqstat{$sample}{basenum});
	push @out,sprintf("%.2f",$fqstat{$sample}{Q30}/$fqstat{$sample}{basenum});
	print Out join("\t",@out),"\n";
}
close Out;
`ln -s $statdir/snp.stat $dOut/Table`;
`ln -s $statdir/tag.stat $dOut/Table`;
mkdir "$dOut/Result" if (!-d "$dOut/Result");
`ln -s $vcfdir/populations.snps.vcf $dOut/Result`;
`ln -s $vcfdir/populations.tag $dOut/Result`;
mkdir "$dOut/Figure" if (!-d "$dOut/Figure");
`ln -s $fastqc/fig/*.pdf $dOut/Figure`;
`ln -s $fastqc/fig/*.png $dOut/Figure`;
`ln -s $statdir/*.pdf $dOut/Figure`;
`ln -s $statdir/*.png $dOut/Figure`;



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
  -vcfdir	<dir>	input vcf file
  -statdir	<dir>	stat dir
  -fastqc	<dir>	fastqc dir
  -out	<dir>	split windows sh
  -dsh	<dir>	output work sh	
  -h         Help

USAGE
        print $usage;
        exit;
}
