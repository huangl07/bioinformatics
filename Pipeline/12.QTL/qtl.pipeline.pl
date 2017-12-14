#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$ann,$pop,$btl,$vcf,$dir,$key);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dir:s"=>\$dir,
	"key:s"=>\$key,
	"out:s"=>\$out,
	"ann:s"=>\$ann,
	"vcf:s"=>\$vcf,
	"pop:s"=>\$pop,
	"btl:s"=>\$btl,
	) or &USAGE;
&USAGE unless ($key and $out and $ann and $vcf and $out);
mkdir $out if (!-d $out);
$pop||="F2";
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$ann=ABSOLUTE_DIR($ann);
$dir=ABSOLUTE_DIR($dir);
my $mark="$dir/$key";
my $oout=$out;
mkdir "$out/work_flow" if (!-d "$out/work_flow");
mkdir "$out/work_sh" if (!-d "$out/work_sh");
my $dsh="$out/work_sh";
$out="$out/work_flow";
mkdir "$out/fig" if (!-d "$out/fig");
open SH,">$dsh/qtl1.sh";
if ($btl) {
	print SH "Rscript $Bin/bin/btl.R --input $mark --output $out/qtl --pop $pop&& ";
}else{
	print SH "Rscript $Bin/bin/qtl.R --input $mark --output $out/qtl --pop $pop && ";
}
close SH;
open SH2,">$dsh/qtl2.sh";
my @out=glob("$out/*.qtl.csv");
foreach my $trt (sort @out) {
	my $trtname=basename($trt);
	$trtname=~s/\.qtl\.csv//g;
	print SH "perl $Bin/bin/region-variant.pl -i $vcf -o $out/$trtname.qtl.vcf -r $trt && ";
	print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out/$trtname.gene -i $trt && ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input $out/$trtname.gene.kegg.stat --output  $out/$trtname.gene.kegg.stat --top 1 && ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input  $out/$trtname.gene.gene.go.stat --output  $out/$trtname.gene.go.stat --top 1&& ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input  $out/$trtname.gene.gene.eggnog.stat --output  $out/$trtname.gene.eggnog.stat --eggnog \n";
}
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
  -vcf	<file>	input vcf file
  -out	<dir>	output dir
  -mark	<str>	input key words of input qtl file(*.csv *.trt)
  -dir	<dir>	input qtl dir
  -ann	<file>	input ann file
  -pop	<str>	pop file
  -btl			binary trt or not
	
  -h         Help

USAGE
        print $usage;
        exit;
}
