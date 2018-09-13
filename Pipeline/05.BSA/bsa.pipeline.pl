#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$pid,$bid,$popt,$ann,$chr,$gff,);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"gff:s"=>\$gff,
	"chr:s"=>\$chr,
	"ann:s"=>\$ann,
	"pid:s"=>\$pid,
	"bid:s"=>\$bid,
	"popt:s"=>\$popt,
			) or &USAGE;
&USAGE unless ($vcf and $out and $bid and $ann);
mkdir $out if (!-d $out);
$popt||="F2";
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$ann=ABSOLUTE_DIR($ann);
my $oout=$out;
mkdir "$out/work_flow" if (!-d "$out/work_flow");
mkdir "$out/work_sh" if (!-d "$out/work_sh");
my $dsh="$out/work_sh";
$out="$out/work_flow";
mkdir "$out/fig" if (!-d "$out/fig");
open SH,">$dsh/bsa.sh";
my @bid=split(/\,/,$bid);
if (scalar @bid == 1) {
	print SH "perl $Bin/bin/mutmap.pl -vcf $vcf -out $out/index-calc.result -bid $bid ";
	if ($pid) {
		print SH "-pid $pid && ";
	}else{
		print SH "&& ";
	}
	print SH "Rscript $Bin/bin/slidingwin.R --infile $out/index-calc.result.index --outfile $out/sliding-win --col 3 --win 2000000 --step 10000 --method bp && ";
	print SH "Rscript $Bin/bin/manhattan.R --input $out/sliding-win.result --output $out/fig/bsa --col 1,3,4,5 && ";
	print SH "perl $Bin/bin/region-variant.pl -i $out/index-calc.result.index -o $out/region.threshold.variant -r $out/sliding-win.threshold.select && ";
	print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out/region.threshold.gene -i $out/sliding-win.threshold.select && ";
	print SH "perl $Bin/bin/region-vcf.pl -i $vcf -o $out/region.threshold.vcf -r $out/sliding-win.threshold.select && ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.kegg.stat --output $out/fig/region.threshold.gene.kegg.stat --top 1 && ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.go.stat --output $out/fig/region.threshold.gene.go.stat --top 1&& ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.eggnog.stat --output $out/fig/region.threshold.gene.eggnog.stat --eggnog && ";
	print SH "perl $Bin/bin/draw-circos.pl --windows 100000 --gwindows 10000 -vcf $vcf --index $out/sliding-win.result --column 4 --chrlist $chr --gff $gff\n";

}else{
	print SH "perl $Bin/bin/qtlseq.pl -vcf $vcf -out $out/index-calc.result -bid $bid -popt $popt ";
	if ($pid) {
		print SH "-pid $pid && ";
	}else{
		print SH "&& ";
	}
	print SH "Rscript $Bin/bin/slidingwin.R --infile $out/index-calc.result.index --outfile $out/sliding-win --col 5 --win 2000000 --step 10000 --method bp && ";
	print SH "Rscript $Bin/bin/manhattan.R --input $out/sliding-win.result --output $out/bsa --col 1,3,4,5,6,7 && ";
	print SH "perl $Bin/bin/region-variant.pl -i $out/index-calc.result.index -o $out/region.threshold.variant -r $out/sliding-win.threshold.select && ";
	print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out/region.threshold.gene -i $out/sliding-win.threshold.select && ";
	print SH "perl $Bin/bin/region-vcf.pl -i $vcf -o $out/region.threshold.vcf -r $out/sliding-win.threshold.select && ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.kegg.stat --output $out/region.threshold.gene.kegg.stat --top 1 && ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.go.stat --output $out/region.threshold.gene.go.stat --top 1&& ";
	print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.eggnog.stat --output $out/region.threshold.gene.eggnog.stat --eggnog&& ";
	print SH "perl $Bin/bin/draw-circos.pl --windows 100000 --gwindows 10000 -vcf $vcf --index $out/sliding-win.result --column 5 --chrlist $chr --gff $gff\n"
}

close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  --Resource mem=20G --CPU 1 $dsh/bsa.sh";
`$job`;
mkdir "$oout/result/" if (!-d "$oout/result");
my $resultdir="$oout/result";
mkdir "$resultdir/report" if (!-d "$oout/report");
mkdir "$resultdir/data/" if (!-d "$resultdir/data");
mkdir "$resultdir/fig/" if (!-d "$resultdir/fig");
`ln -s $out/index-calc.result.stat $resultdir/report/3-12.xls`;
`grep @ $out/region.threshold.vcf.total|sed 's/@//g'|less -S > $resultdir/report/3-13.xls`;
`perl $Bin/gatherdata.pl -i $out/region.threshold.gene.total -o $resultdir/report/`;
`ln -s $out/index-calc.result $resultdir/data`;
`ln -s $out/region*.total $resultdir/data`;
`ln -s $out/region*.eff $resultdir/data`;
`ln -s $out/fig/*.detail $resultdir/data`;
`ln -s $out/fig/*.pdf $resultdir/fig`;
`ln -s $out/fig/*.png $resultdir/fig`;

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
  -vcf	<file>	input vcf file
  -out	<dir>	output dir
  -pid	<id>	input p id
  -bid	<id>	input b id
  -ann	<file>	input ann file
  -h         Help

USAGE
        print $usage;
        exit;
}
