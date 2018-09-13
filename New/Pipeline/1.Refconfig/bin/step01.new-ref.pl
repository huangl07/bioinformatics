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
&USAGE unless ($ref and $out and $dsh and $gff);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step01.new-ref.sh";
print SH "gffread -o $out/ref.packaged.gff $gff && perl $Bin/bin/GRename.pl -i $ref -g $out/ref.packaged.gff -o $out/ref  ";
if ($chr) {
	print SH "-f $chr ";
}
print SH " && perl $Bin/bin/getGeneFasta.pl -i $out/ref.fa -o $out/ref.gene.fa -g $out/ref.gff && ";
print SH "perl $Bin/bin/pre-design.pl -i $ref -o $out/ref.predesign\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step01.new-ref.sh";
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

Usage:
  Options:
  -ref	<file>	input genome name,fasta format,
  -gff	<file>	input genome gff file,
  -out	<dir>	output data prefix
  -chr	<file>	chromosome change file
  -dsh	<dir>	output work sh dir

  -h         Help

USAGE
        print $usage;
        exit;
}
