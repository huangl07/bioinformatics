#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($vcf,$dsh,$out,$popt,$pid,$mid);
GetOptions(
				"help|?" =>\&USAGE,
				"vcf:s"=>\$vcf,
				"popt:s"=>\$popt,
				"out:s"=>\$out,
				"dsh:s"=>\$dsh,
				"pid:s"=>\$pid,
				"mid:s"=>\$mid,
				) or &USAGE;
&USAGE unless ($vcf and $out and $dsh and $popt);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step01.vcf-convert.sh";
print SH "perl $Bin/bin/vcf2marker.pl -vcf $vcf -out $out/pop.primary.marker -PID $pid -MID $mid && ";
print SH "perl $Bin/bin/markerfilter.pl -input $out/pop.primary.marker -output $out/pop -popt $popt \n";
close SH;
my $job="perl  /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dsh/step01.vcf-convert.sh";
`$job`;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:

	-vcf	<file>	input vcf file
	-popt	<str>	population type
	-out	<dir>	output file dir
	-dsh	<dir>	output shell dir
	-pid	<str>	paternal id
	-mid	<str>	maternal id

USAGE
	print $usage;
	exit;
}
