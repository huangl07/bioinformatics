#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($option,$trt,$tassel,$out,$dsh,$hmp);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "trt:s"=>\$trt,
    "tassel:s"=>\$tassel,
    "hmp:s"=>\$hmp,
	"out:s"=>\$out,
    "dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($trt and $tassel and $hmp and $out and $dsh);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$trt=ABSOLUTE_DIR($trt);
$tassel=ABSOLUTE_DIR($tassel);
$dsh=ABSOLUTE_DIR($dsh);
$hmp=ABSOLUTE_DIR($hmp);

open SH1,">$dsh/step02.gwas-gapit.sh";
open SH2,">$dsh/step02.gwas-mvp.sh";
open SH3,">$dsh/step02.gwas-tassel1.sh";
open SH4,">$dsh/step02.gwas-tassel2.sh";
open OUT,">$out/gwas.dir.list";

open IN,$trt;

while(<IN>){
    $_=~s/[\n\r]//g;
    my ($trit,$file)=split/\t/,$_;
    mkdir "$out/$trit" if (!-d "$out/$trit");
    print OUT "$trit\t$out/$trit\n";
    print SH1 "Rscript $Bin/bin/GWAS-GAPIT.R  --hapmap $hmp --trait $file --out $out/$trit/ \n";    
    print SH2 "Rscript $Bin/bin/MVP.single.R --hapmap $hmp --trait $file --output $out/$trit/ \n";
}
close IN;
open IN,$tassel;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($trit,$file)=split/\t/,$_;
    print SH3 "/mnt/ilustre/users/dna/.env/bin/tassel -Xmx30G -fork1 -h $hmp -KinshipPlugin -endPlugin -export $out/$trit/kinship && /mnt/ilustre/users/dna/.env/bin/tassel -Xmx30G -fork1 -h $hmp -PrincipalComponentsPlugin -covariance true -endPlugin -export $out/$trit/pca\n";
    print SH4 "/mnt/ilustre/users/dna/.env/bin/tassel -fork1 -h $hmp -fork2 -t $file -fork3 -r $out/$trit/pca1.txt -excludeLastTrait -fork4 -k $out/$trit/kinship.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel Optimum -export $out/$trit/$trit\_tassel\n";
    }
close IN;
close SH1;
close SH2;
close SH3;
close SH4;
close OUT;

my $job1="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=200G --CPU 8 $dsh/step02.gwas-gapit.sh";
#`$job1`;
my $job2="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=150G --CPU 8 $dsh/step02.gwas-mvp.sh";
#`$job2`;
my $job3="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=100G --CPU 8 $dsh/step02.gwas-tassel1.sh";
`$job1;$job2;$job3`;
my $job4="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=100G --CPU 8 $dsh/step02.gwas-tassel2.sh";
`$job4`;
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
	filter snp and retrive fa 
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -trt  <file>  input trit list
  -tassel   <file>  input tassel trit list
  -hmp  <file>    input hapmap file
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
