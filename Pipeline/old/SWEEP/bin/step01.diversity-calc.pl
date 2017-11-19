#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$gro);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"gro:s"=>\$gro,
	"dep:s"=>\$dep,
			) or &USAGE;
&USAGE unless ($vcf and $out and $dsh );
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$vcf=ABSOLUTE_DIR($vcf);
$gro=ABSOLUTE_DIR($gro);
open In,$gro;
my %gro;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$pop)=split(/\s+/,$_);
	$gro{$id}=$pop;
}
close In;
my @group=sort keys %gro;
open SH,">$dsh/step01.diversity-calc.sh";
for (my $i=0;$i<@group;$i++) {
	for (my $j=1;$j<@group;$i++) {
		print SH "vcftools --vcf $vcf --weir-fst-pop $group[$i] --weir-fst-pop $group[$j] --fst-window-size 200000 --fst-window-step 10000 --out $out/$group[$i]-$group[$j]\n";
	}
	print SH "vcftools --vcf $vcf --keep $group[$i] --window-pi 2000000 --window-pi-step 10000 --out $out/$group[$i]\n";
	print SH "vcftools --vcf $vcf --keep $group[$i] --TajimaD 100000 --out $out/$group[$i]\n";
	print SH "cd $out/ && OmegaPlus-F -name $group[$i] -input $vcf  -grid 100 -threads 8 -minwin 1000 -maxwin 10000 -seed 12345 -sampleList $group[$i] && perl $Bin/bin/omega.pl -i $out/SweeD_Report.$group[i].sweep -o $out/Sweed.$group[$i].result\n";
	print SH "cd $out/ && SweeD-MPFR-P -name $group[$i] -input $vcf  -grid 100 -threads 8 -sampleList $group[$i] && perl $Bin/bin/sweed.pl -i $out/SweeD_Report.$group[i].sweep -o $out/Sweed.$group[$i].result\n";

}
print SH "vcftools --vcf $vcf --window-pi 2000000 --window-pi-step 10000 "
close SH;
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
  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -dsh	<dir>	output work shell
  -gro	<str>	group list
  -maf	<num>	maf filter default 0.05
  -mis	<num>	mis filter default 0.3
  -dep	<num>	dep filter default 2

  -h         Help

USAGE
        print $usage;
        exit;
}
