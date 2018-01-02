#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($sv,$out,$bam,$ref);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"sv:s"=>\$sv,
	"out:s"=>\$out,
	"bam:s"=>\$bam,
	"ref:s"=>\$ref,
			) or &USAGE;
&USAGE unless ($sv and $out and $bam and $ref);
mkdir $out if (!-d $out);
my $dsh="$out/work_sh";
mkdir $dsh if (!-d $dsh);
$bam=ABSOLUTE_DIR($bam);
$dsh=ABSOLUTE_DIR($dsh);
$out=ABSOLUTE_DIR($out);
open SH,">$dsh/transgenic0.sh";
print SH "ln -s $ref $out/ref.fa && ";
print SH "makeblastdb -in $out/ref.fa -input_type fasta -dbtype nucl \n ";
close SH;

open In,$bam;
my %bam;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$bam)=split(/\s+/,$_);
	$bam{$sample}=$bam;
}
close In;
open In,$sv;
open SH,">$dsh/transgenic1.sh";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$sv)=split(/\s+/,$_);
	print SH "perl $Bin/bin/bed-creat.pl -i $sv -r $ref -o $out/$sample.bed && ";
	print SH "samtools view -L $out/$sample.bed -h $bam{$sample}|samtools fastq -N -1 $out/$sample.1.fastq -2 $out/$sample.2.fastq - && ";
	open Config,">$out/$sample.config";
	print Config "max_rd_len=150\n";
	print Config "[LIB]\n";
	print Config "avg_ins=400\n";
	print Config "reverse_seq=0\n";
	print Config "asm_flags=3\n";
	print Config "rank=1\n";
	print Config "q1=$out/$sample.1.fastq\n";
	print Config "q2=$out/$sample.2.fastq\n";
	close Config;
	print SH "SOAPdenovo-127mer all -s $out/$sample.config -K 63 -R -o $out/$sample -p 20 1> $out/$sample.log 2>$out/$sample.err && ";
	print SH "blastn -db $out/ref.fa -query $out/$sample.scafSeq -out $out/$sample.blast -evalue=10e-5 -outfmt 6\n";
}
close In;
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

  -sv	<file>	input sv anno files
  -out	<dir>	output dir
  -bam	<file>	bam file
  -ref	<file>	input reference files
  -h         Help

USAGE
        print $usage;
        exit;
}
