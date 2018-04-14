#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($variation,$bid,$pid,$cg,$out,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "pid:s"=>\$pid,
    "bid:s"=>\$bid,
	"variation:s"=>\$variation,
    "cg:s"=>\$cg,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($variation and $cg and $out);
mkdir $out if (!-d $out);
mkdir "$out/01.fastq_qc" if (!-d "$out/01.fastq_qc");
mkdir "$out/02.reference" if (!-d "$out/02.reference");
mkdir "$out/03.map_stat" if (!-d "$out/03.map_stat");
mkdir "$out/04.variant-stat" if (!-d "$out/04.variant-stat");
mkdir "$out/05.annovar" if (!-d "$out/05.annovar");
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$cg=ABSOLUTE_DIR($cg);
my $qc=ABSOLUTE_DIR("$out/01.fastq_qc");
my $ref=ABSOLUTE_DIR("$out/02.reference");
my $map=ABSOLUTE_DIR("$out/03.map_stat");
my $var=ABSOLUTE_DIR("$out/04.variant-stat");
my $anno=ABSOLUTE_DIR("$out/05.annovar");
$dsh=ABSOLUTE_DIR($dsh);
$variation=ABSOLUTE_DIR($variation);

######
######01.fastq_qc
my @atgc=glob("$variation/01.fastq-qc/*.atgc");
mkdir "$qc/atgc" if (!-d "$qc/atgc");
for(@atgc){
    system("cp $_ $qc/atgc/\n");
    }
mkdir "$qc/qual" if (!-d "$qc/qual");
my @qual=glob("$variation/01.fastq-qc/*.qual");
for(@qual){
    system("cp $_ $qc/qual/\n");
    }
mkdir "$qc/stat" if (!-d "$qc/stat");
system("cp $variation/14.report/Table/3-3.xls $qc/stat/QC.xls\n");
######
######02.reference
open OUT,">$out/info.list";
my @pid=split/\,/,$pid;
my @bid=split/\,/,$bid;

print OUT "PID\t",join("\t",@pid),"\n";
print OUT "BID\t",join("\t",@bid);
close OUT;
system("cp $cg $ref/ref.changelog \n");
system("cp $variation/02.ref-config/ref.chrlist $ref/ \n");
#system("cp $variation/02.ref-config/ref.fa $ref/ \n");
system("cp $variation/02.ref-config/ref.fa $ref && gzip -9 $ref/ref.fa \n");
system("cp $variation/02.ref-config/ref.log $ref/ \n");
######
######03.map_stat
mkdir "$map/coverage" if (!-d "$map/coverage");
mkdir "$map/depth" if (!-d "$map/depth");
mkdir "$map/insert" if (!-d "$map/insert");
mkdir "$map/result.stat" if (!-d "$map/result.stat");
my @coverage=glob("$variation/06.map-stat/*.depth.fordraw\n");
my @depth=glob("$variation/06.map-stat/*.coverage\n");
my @insert=glob("$variation/06.map-stat/*.insert\n");

for(@depth){my $newfile=basename($_);my $files=(split/\./,$newfile)[0];system("cp $_ $map/depth/$files.depth\n")};
for(@coverage){my $newfile=basename($_);my $files=(split/\./,$newfile)[0];system("cp $_ $map/coverage/$files.coverage\n")};
for(@insert){system("cp $_ $map/insert/\n")};
system("cp $variation/14.report/Table/Total.mapped.detail $map/result.stat/\n");
######
######04.variant-stat
mkdir "$var/snp" if (!-d "$var/snp");
mkdir "$var/indel" if (!-d "$var/indel");
system("cp $variation/13.variant-stat/snp.effects $var/snp/\n");
system("cp $variation/13.variant-stat/snp.region $var/snp/\n");
system("cp $variation/14.report/Table/3-6.xls $var/snp/snp.stat\n");
system("cp $variation/13.variant-stat/indel.effects $var/indel/\n");
system("cp $variation/13.variant-stat/indel.len $var/indel/\n");
system("cp $variation/13.variant-stat/indel.region $var/indel/\n");
system("cp $variation/14.report/Table/3-10.xls $var/indel/indel.stat \n");
######
######05.annovar
if(-e "$variation/10.annovar/pop.final.vcf.gz"){
    system("cp $variation/10.annovar/pop.final.vcf.gz $anno/ \n");
    }else{
        system("cp $variation/10.annovar/pop.final.vcf $anno && gzip -9 $anno/pop.final.vcf \n");
        }
open IN,"$variation/10.annovar/pop.summary";
open OUT,">$anno/pop.summary";
while(<IN>){
    $_=~s/[\n\r]//g;
    if(/^#/){
        print OUT "$_\n";
     }else{
            my @array=split/\t/,$_;
            $array[15]=~s/\:/\,/g;
			if ($array[19]!~/:/){
				$array[19]="$array[19]:$array[20]";
				$array[20]=$array[19];
			}
			print OUT join("\t",@array),"\n";
    }
}
close IN;
close OUT;
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
	
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -pid <string>  example:p1,p2
  -bid  <string>  example:b1,b2
  -variation	<dir>	input variation dir;
  -cg   <file>  change.log file
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
