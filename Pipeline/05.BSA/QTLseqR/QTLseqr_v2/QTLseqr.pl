#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$ref,$ann,$highBulkid,$lowBulkid);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"ref:s"=>\$ref,
	"ann:s"=>\$ann,
	"HB:s"=>\$highBulkid,
	"LB:s"=>\$lowBulkid,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($vcf and $out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
mkdir "$out/work_flow" if (!-d "$out/work_flow");
mkdir "$out/work_sh" if (!-d "$out/work_sh");
mkdir "$out/result" if (!-d "$out/result");
my $resultdir = "$out/result";
mkdir "$resultdir/Result" if (!-d "$out/Result");
mkdir "$resultdir/Table/" if (!-d "$resultdir/Table");
mkdir "$resultdir/Figure/" if (!-d "$resultdir/Figure");
my $dsh = "$out/work_sh";
$out = "$out/work_flow";

open SH,">$dsh/01.vcf.convert.sh";
print SH "java -Xmx50G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T VariantsToTable -R $ref -V $vcf -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL -o $out/pop.vcf.table";
open SH,">$dsh/02.bsa.sh";
print SH "Rscript $Bin/bin/QTLseqr.R --input $out/pop.vcf.table --output $out --HB $highBulkid --LB $lowBulkid && grep -v sca $out/G.analysis.result|sed 's/chr//g' >$out/G.nosca.result && ";
print SH "perl $Bin/bin/region.pl -select $out/filter.index.result -out $out/qtl.csv && ";
print SH "perl $Bin/bin/variant.vcf.pl -G $out/G.analysis.result -vcf $vcf -o $out && ";
print SH "perl $Bin/bin/region-variant.pl -i $out/variant.region.vcf -o $out/region.threshold.variant -r $out/qtl.csv && ";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out/region.threshold.gene -i $out/qtl.csv && ";
print SH "perl $Bin/bin/region-vcf.pl -i $vcf -o $out/region.threshold.vcf -r $out/qtl.csv && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.kegg.stat --output $out/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.go.stat --output $out/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out/region.threshold.gene.eggnog.stat --output $out/region.threshold.gene.eggnog.stat --eggnog && ";

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=50G --CPU 2  $dsh/01.vcf.convert.sh";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=20G --CPU 2  $dsh/02.bsa.sh";
`$job`;

`cp $out/variant.region.stat $resultdir/Table/3-12.xls`;
`perl $Bin/gatherdata.pl -i $out/region.threshold.gene.total -o $resultdir/Table`;
`cut -d "," -f1,3-15 $out/qtl.csv >$resultdir/Result/qtl.region.xls`;
`grep "\@chr" $out/region.threshold.variant.total|sed 's/@//g' > $resultdir/Table/3-13.xls `;
`cp $out/region*.total $resultdir/Result`;
`cp $out/region*.eff $resultdir/Result`;
`cp $out/*.stat $resultdir/Result`;
`cp $out/*.pdf $resultdir/Figure`;
`cp $out/*.png $resultdir/Figure`;
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
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf, pop.final.vcf
	"ref:s"=>\$ref, ref.fa
	"HB:s"=>\$highBulkid, highBulkid
	"LB:s"=>\$lowBulkid, lowBulkid
	"out:s"=>\$out, output dir
	"ann:s"=>\$ann, pop.summary
USAGE
        print $usage;
        exit;
}
