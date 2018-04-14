#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$ref,$type,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"ref:s"=>\$ref,
	"type:s"=>\$type,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($vcf and $ref and $type and $out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$vcf=ABSOLUTE_DIR($vcf);
$ref=ABSOLUTE_DIR($ref);

my $tmp=time();
open SH,">$out/work_sh/primer.$tmp.sh";
print SH "perl $Bin/bin/01.misa.pl -input $vcf -output $out -type $type && ";
print SH "perl $Bin/bin/02.primer_p3in.pl -input $ref -output $out -type $type && ";
print SH "cd ~dna/Environment/biotools/primer3/src/ && ./primer3_core < $out.$type.p3in > $out.$type.p3out && ";
print SH "perl $Bin/bin/03.merge_p3out.pl $out.$type.p3out $out.$type.misa && ";
print SH "perl $Bin/bin/04.result.pl -input $out.$type.result -output $out -type $type ";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl --Resource mem=16G --CPU 1 $out/work_sh/primer.$tmp.sh";
`job`;

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
  -vcf	<file>	input variant vcf file
  -ref	<dir>	input ref genome file
  -out	<dir>	output file dir
  -type	<str>	marker type
  -h         Help

USAGE
        print $usage;
        exit;
}
