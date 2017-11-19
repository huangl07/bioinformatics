#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$fgff,$dOut,$dShell);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"gff:s"=>\$fgff,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($ref and $fgff and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$ref=ABSOLUTE_DIR($ref);
$fgff=ABSOLUTE_DIR($fgff);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
open SH,">$dShell/step00.genome-config.sh";
`ln -s $ref $dOut/ref.fa`;
`ln -s $fgff $dOut/ref.gff`;
print SH "cd $dOut && bwa index $dOut/ref.fa\n";
print SH "cd $dOut && java -jar /mnt/ilustre/users/dna/.env/bin/picard.jar CreateSequenceDictionary REFERENCE=$dOut/ref.fa && ";
print SH "perl $Bin/bin/chr.pl -i $dOut/ref.dict -o $dOut/ref.chrlist\n";
print SH "cd $dOut && samtools faidx $dOut/ref.fa\n";
print SH "cd $dOut && perl $Bin/bin/annovar/gffread.pl -i $dOut/ref.gff -o $dOut/ref.gtf && ";
print SH "/mnt/ilustre/users/dna/.env/bin/gtfToGenePred -allErrors -genePredExt $dOut/ref.gtf $dOut/ref\_refGene.txt && ";
print SH "$Bin/bin/annovar/retrieve_seq_from_fasta.pl -format refGene -seqfile $dOut/ref.fa -outfile $dOut/ref\_refGeneMrna.fa $dOut/ref\_refGene.txt \n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dShell/step00.genome-config.sh";
print $job;
`$job`;
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
	fq thanslate to ref format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -ref	<file>	input genome ref file
  -gff	<file>	input gff file
  -out	<dir>	output dir
  -dsh	<dir>	output worksh dir
  -h         Help

USAGE
        print $usage;
        exit;
}
