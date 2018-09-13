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
my $check_sample = `wc -l $slist/sstacks.list`;
chomp $check_sample;
$check_sample = (split(/\s+/,$check_sample))[0];
my $check_log = `ls $slist/\*.matches.tsv.gz|wc -l`;
chomp $check_log;
$check_log = (split(/\s+/,$check_log))[0];
if ($check_sample ne $check_log) {
	print "There is some wrong in step 5 ,please check!";die;
}
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
my $Slist = "$slist/sstacks.list";
open In,$Slist;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$sstacks)=split(/\s+/,$_);
	`ln -s $sstacks* $dOut/`;
}
open SH,">$dShell/step06.genotype.sh";
print SH "/mnt/ilustre/users/dna/.env/stacks-2.1/bin/tsv2bam -P $dOut -M $popmap -t 16 && ";
print SH "/mnt/ilustre/users/dna/.env/stacks-2.1/bin/gstacks -P $dOut -M $popmap -t 16 && ";
print SH "/mnt/ilustre/users/dna/.env/stacks-2.1/bin/populations -P $dOut -M $popmap -t 16 -O $dOut --vcf --fasta_loci && ";
print SH "perl $Bin/bin/tag-generate.pl -vcf $dOut/populations.snps.vcf -catalog $dOut/catalog.tags.tsv.gz -out $dOut/populations.tag && ";
print SH "perl $Bin/bin/vcf-convert.pl -i $dOut/populations.snps.vcf -o $dOut/pop.recode.vcf\n";
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=256G --CPU 16 --Nodes 1 $dShell/step06.genotype.sh";
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
	"ulist:s"=>\$ulist,   input ulist.list 
	"clist:s"=>\$clist,	  input clist.list  
	"slist:s"=>\$slist,	  input slist dir
	"gro:s"=>\$popmap,	input group.list
	"out:s"=>\$dOut,   output dir
	"dsh:s"=>\$dShell,  work_sh dir
	 -h         Help

USAGE
        print $usage;
        exit;
}
