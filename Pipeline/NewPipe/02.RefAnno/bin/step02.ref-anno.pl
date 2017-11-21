#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$out,$gff,$chr,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"out:s"=>\$out,
	"gff:s"=>\$gff,
	"chr:s"=>\$chr,
	"dsh:s"=>\$dsh,
	) or &USAGE;
&USAGE unless ($ref and $out and $dsh and $gff and $chr);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
open SH,">$dsh/step01.new-ref.sh";
print SH "perl $Bin/bin/GRename.pl -i $ref -g $gff -o $out/ref -f $chr && ";
print SH "perl $Bin/bin/getGeneFasta.pl -i $ref -o $out/ref.gene.fa -g $gff ^^ ";
print SH "perl $Bin/bin/pre-design.pl -i $ref -o $out/ref.predesign\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step01.new-ref.sh";
`$job`;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -i	<file>	input genome name,fasta format,
  -g	<file>	input genome gff file,
  -o	<str>	output file prefix
  -f	<file>	chromosome change file

  -h         Help

USAGE
        print $usage;
        exit;
}
