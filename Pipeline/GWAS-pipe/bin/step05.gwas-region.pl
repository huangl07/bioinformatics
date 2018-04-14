#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($list1,$list2,$anno,$vcf,$out,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "list1:s"=>\$list1,
    "list2:s"=>\$list2,
    "anno:s"=>\$anno,
    "vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($vcf and $anno and $list1 and $list2 and $out and $dsh);
mkdir $out if (!-d $out);
mkdir "$out/Bonferroni1" if (!-d "$out/Bonferroni1");
mkdir "$out/Bonferroni2" if (!-d "$out/Bonferroni2");
mkdir $dsh if (!-d $dsh);
my $out1="$out/Bonferroni1";
my $out2="$out/Bonferroni2";
$out1=ABSOLUTE_DIR($out1);
$out2=ABSOLUTE_DIR($out2);

$dsh=ABSOLUTE_DIR($dsh);
$vcf=ABSOLUTE_DIR($vcf);
$list1=ABSOLUTE_DIR($list1);
$list2=ABSOLUTE_DIR($list2);
$anno=ABSOLUTE_DIR($anno);

open SH,">$dsh/step05.gwas-region.sh";
open IN,$list1;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($trit,$file)=split/\t/,$_;
    print SH "perl $Bin/bin/region-gene.pl -anno $anno  -out $out1/$trit  -select $file && Rscript $Bin/bin/eff-enrich.R --input $out1/$trit.kegg.stat --output $out1/$trit.kegg.stat --top 1 && Rscript $Bin/bin/eff-enrich.R --input $out1/$trit.go.stat --output $out1/$trit.go.stat --top 1 && Rscript $Bin/bin/eff-enrich.R --input $out1/$trit.eggnog.stat --output  $out1/$trit.eggnog.stat --eggnog && perl $Bin/bin/region-vcf.pl -i $vcf  -o $out1/$trit.vcf -r $file\n";
    }
close IN;
open IN,$list2;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($trit,$file)=split/\t/,$_;
    print SH "perl $Bin/bin/region-gene.pl -anno $anno  -out $out2/$trit  -select $file && Rscript $Bin/bin/eff-enrich.R --input $out2/$trit.kegg.stat --output $out2/$trit.kegg.stat --top 1 && Rscript $Bin/bin/eff-enrich.R --input $out2/$trit.go.stat --output $out2/$trit.go.stat --top 1 && Rscript $Bin/bin/eff-enrich.R --input $out2/$trit.eggnog.stat --output  $out2/$trit.eggnog.stat --eggnog && perl $Bin/bin/region-vcf.pl -i $vcf  -o $out2/$trit.vcf -r $file\n";
}
close IN;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step05.gwas-region.sh";
`$job`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	Gwas region stat
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -list1 <file>  input gwas Bonferroni1 list
  -list2 <file>  input gwas Bonferroni2 list
  -anno  <file>  input anno file
  -vcf  <file> input vcf file
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
