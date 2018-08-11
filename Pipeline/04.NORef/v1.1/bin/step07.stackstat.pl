#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$dOut,$dShell,$ulist,$fastqc);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"ulist:s"=>\$ulist,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($vcf and $dOut and $dShell and $ulist);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
$dOut=ABSOLUTE_DIR($dOut);
$vcf=ABSOLUTE_DIR($vcf);
$ulist=ABSOLUTE_DIR($ulist);
open SH,">$dShell/step07.stack-stat.sh";
print SH "perl $Bin/bin/snp.stat.pl -i $vcf -o $dOut/snp.stat -m $dOut/snp.matrix && ";
print SH "Rscript $Bin/bin/diff_matrix.R --i $dOut/snp.matrix --o $dOut/snp.diff\n";
print SH "perl $Bin/bin/variant_qual.pl -i $vcf -o1 $dOut/snp.depth &&";
print SH "Rscript $Bin/bin/snp_qual.R --dep $dOut/snp.depth --o $dOut/snp.qual \n";
print SH "perl $Bin/bin/merge.pl -i $ulist -o $dOut/tag.stat\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=20G --CPU 1 --Nodes 1 $dShell/step07.stack-stat.sh";
print "$job\n";
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
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input vcf file
  -sstacks	<file>	sstacks list
  -out	<dir>	split windows sh
  -dsh	<dir>	output work sh	
  -h         Help

USAGE
        print $usage;
        exit;
}
