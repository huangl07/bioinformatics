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
my ($mlod,$marker,$chrlist,$out,$dsh,$nchr,$ref,$popt);
GetOptions(
				"help|?" =>\&USAGE,
				"mlod:s"=>\$mlod,
				"marker:s"=>\$marker,
				"nchr:s"=>\$nchr,
				"popt:s"=>\$popt,
				"ref"=>\$ref,
				"out:s"=>\$out,
				"dsh:s"=>\$dsh,
				) or &USAGE;
&USAGE unless ($mlod and $out and $dsh);
$mlod=ABSOLUTE_DIR($mlod);
$marker=ABSOLUTE_DIR($marker);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step04.grouping.sh";
if ($ref) {
	print SH "perl $Bin/bin/linkage_by_ref.pl -i $mlod -2 $marker -o $out -k Total && ";
#	if ($popt eq "CP") {
#		print SH "perl $Bin/bin/splitbyLG-CP.pl -l $out/Total.lg -i $marker -d $out/ -t $popt";
#	}else{
#		print SH "perl $Bin/bin/splitbyLG-NOCP.pl -l $out/Total.lg -i $marker -d $out/ -t $popt";
#	}
}else{
	print SH "perl $Bin/bin/linkage_by_mlod.pl -i $mlod -k Total -d $out -n $nchr && ";
#	if ($popt eq "CP") {
#		print SH "perl $Bin/bin/splitbyLG-CP.pl -l $out/Total.lg -i $marker -d $out/ -t $popt";
#	}else{
#		print SH "perl $Bin/bin/splitbyLG-NOCP.pl -l $out/Total.lg -i $marker -d $out/ -t $popt";
#	}

}
close SH;
my $mem=`du $mlod`;
chomp $mem;
$mem=(split(/\s+/,$mem))[0];
$mem=int($mem/1000000)+3;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=$mem"."G --CPU 1 $dsh/step04.grouping.sh";
print $job;
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
	-mlod	<file>	input vcf list
	-marker	<file>	marker file
	-popt	<str>	population type
	-nchr	<num>	chr num
	-out	<out>	output dir
	-dsh	<dir>	output work shell file
USAGE
	print $usage;
	exit;
}
