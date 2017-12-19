#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$pop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"pop:s"=>\$pop,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
			) or &USAGE;
&USAGE unless ($vcf and $out and $pop);
$vcf=ABSOLUTE_DIR($vcf);
$pop=ABSOLUTE_DIR($pop);

mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$dsh="$out/work_sh";
mkdir $dsh if (!-d $dsh);
$mis||=0.3;
$maf||=0.05;
$dep||=2;
$mis=1-$mis;
open SH,">$dsh/01.filtered.sh";
print SH "vcftools --vcf $vcf  --out $out/pop --max-missing $mis --maf $maf --minDP $dep --recode \n ";
close SH;
open SH,">$dsh/02.calculate.sh";
open SH2,">$dsh/03.draw-select.sh";
open In,$pop;
my %group;
my %filehand;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$gid)=split(/\s+/,$_);
	if (!exists $filehand{$gid}) {
		open $filehand{$gid},">$out/$gid.list";
	}
	print {$filehand{$gid}} "$id\n";
}
close In;
my @groid=sort keys %filehand;
for (my $i=0;$i<@groid;$i++) {
	print SH "vcftools --vcf $out/pop.recode.vcf --keep $out/$groid[$i].list --out $out/$groid[$i]  --window-pi 2000000 --window-pi-step 10000 \n";
	print SH "vcftools --vcf $out/pop.recode.vcf --keep $out/$groid[$i].list --out $out/$groid[$i] --TajimaD 10000 \n";
	print SH2 "Rscript $Bin/bin/pi-tajima.R --tajima $out/$groid[$i].Tajima.D --pi1 $out/$groid[$i].windowed.pi --pi2 $out/$groid[$i].windowed.pi --out $out/$groid[$i]\n";
	for (my $j=$i+1;$j<@groid;$j++) {
		print SH "vcftools --vcf $out/pop.recode.vcf --weir-fst-pop $out/$groid[$i].list --weir-fst-pop $out/$groid[$j].list --out $out/$groid[$i]-$groid[$j] --fst-window-size 2000000 --fst-window-step 10000 \n";
		print SH2 "Rscript $Bin/bin/fst-pi.R --fst $out/$groid[$i]-$groid[$j].windowed.weir.fst --pi1 $out/$groid[$i].windowed.pi --pi2 $out/$groid[$j].windowed.pi --out $out/$groid[$i]-$groid[$j] \n";
		print SH2 "Rscript $Bin/bin/manhattan.R --fst $out/$groid[$i]-$groid[$j].windowed.weir.fst --pi1 $out/$groid[$i].windowed.pi --pi2 $out/$groid[$j].windowed.pi --out $out/$groid[$i]-$groid[$j] --tajima1 $out/$groid[$i].Tajima.D --tajima2 $out/$groid[$j].Tajima.D \n";
	}
}
close SH;
close SH2;
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
  -pop	<str>	group list
  -maf	<num>	maf filter default 0.05
  -mis	<num>	mis filter default 0.3
  -dep	<num>	dep filter default 2

  -h         Help

USAGE
        print $usage;
        exit;
}
