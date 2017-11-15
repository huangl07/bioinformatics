#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$gff,$out,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($ref and $gff and $out and $dsh);
mkdir $out if (!-d $out);
$ref=ABSOLUTE_DIR($ref);
$gff=ABSOLUTE_DIR($gff);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/02.ref-config.sh";
`ln -s $ref $out/ref.fa`;
`ln -s $gff $out/ref.gff`;
print SH "cd $out && bwa index $out/ref.fa\n";
print SH "cd $out && java -jar /mnt/ilustre/users/dna/.env/bin/picard.jar CreateSequenceDictionary REFERENCE=$out/ref.fa && ";
print SH "perl $Bin/bin/chr.pl -i $out/ref.dict -o $out/ref.chrlist\n";
print SH "cd $out && samtools faidx $out/ref.fa\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dsh/02.ref-config.sh";
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
