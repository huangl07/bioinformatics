#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ulist,$clist,$dOut,$dShell,$slist,$popmap);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ulist:s"=>\$ulist,
	"clist:s"=>\$clist,
	"slist:s"=>\$slist,
	"gro:s"=>\$popmap,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($ulist and $clist and $slist and $dOut and $dShell and $popmap);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if(!-d $dShell);
$dOut=ABSOLUTE_DIR($dOut);
$dShell=ABSOLUTE_DIR($dShell);
$popmap=ABSOLUTE_DIR($popmap);
open In,$ulist;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$ustacks)=split(/\s+/,$_);
	`ln -s $ustacks* $dOut/`;
	my $ustack=basename($ustacks);
}
close In;
open In,$clist;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$cstacks)=split(/\s+/,$_);
	`ln -s $cstacks* $dOut/`;
}
close In;
open In,$slist;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$sstacks)=split(/\s+/,$_);
	`ln -s $sstacks* $dOut/`;
}
open SH,">$dShell/step06.genotype.sh";
print SH "/mnt/ilustre/users/dna/.env/stacks-2.1/bin/tsv2bam -P $dOut -M $popmap -t 36 && ";
print SH "/mnt/ilustre/users/dna/.env/stacks-2.1/bin/gstacks -P $dOut -M $popmap -t 36 && ";
print SH "/mnt/ilustre/users/dna/.env/stacks-2.1/bin/populations -P $dOut -M $popmap -t 36 -O $dOut --vcf && ";
print SH "perl $Bin/bin/tag-generate.pl -vcf $dOut/populations.snps.vcf -catalog $dOut/catalog.tags.tsv.gz -out $dOut/populations.tag && ";
print SH "perl $Bin/bin/vcf-convert.pl -i $dOut/populations.snps.vcf -o $dOut/pop.recode.vcf\n";
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=200G --CPU 16 --Nodes 1 $dShell/step06.genotype.sh";
`$job`;
print "$job\tdone!\n";

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
	"ulist:s"=>\$ulist,   input ulist
	"clist:s"=>\$clist,	  input clist
	"slist:s"=>\$slist,	  input slist
	"gro:s"=>\$popmap,	input group.list
	"out:s"=>\$dOut,   output dir
	"dsh:s"=>\$dShell,  work_sh dir
	 -h         Help

USAGE
        print $usage;
        exit;
}
