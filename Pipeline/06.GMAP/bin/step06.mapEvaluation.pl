#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$out,$dsh,$popt,$mark,$ref);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"mark:s"=>\$mark,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"ref"=>\$ref,
	"popt:s"=>\$popt,
			) or &USAGE;
&USAGE unless ($dmap and $out and $dsh and $popt);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$mark=ABSOLUTE_DIR($mark);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$dmap=ABSOLUTE_DIR($dmap);
open SH1,">$dsh/step06.mapEvaluation1.sh";
open SH2,">$dsh/step06.mapEvaluation2.sh";
if ($popt ne "CP") {
	$popt=lc($popt);
	print SH1 "perl $Bin/bin/MapMergeNOCP.pl -dmap $dmap -o $out && ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.map -o $out/sexAver.mapstat \n ";
	print SH2 "Rscript $Bin/bin/drawmap.R --mark $out/total  --out $out/total --pop $popt \n";
	print SH2 "Rscript $Bin/bin/drawbinNOCP.R --mark $out/total.csv  --out $out/total.bin \n";
	if ($ref) {
		print SH2 "Rscript $Bin/bin/drawAligmentRalationMap.pl -m $out/total.map --out $out/ -k total.phy\n";
	}
}else{
	print SH1 "perl $Bin/bin/MapMergeCP.pl -dmap $dmap -o $out --mark $mark && ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.sexAver.map -o $out/sexAver.mapstat \n ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.male.map -o $out/male.mapstat \n ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.female.map -o $out/female.mapstat \n ";
	print SH2 "Rscript $Bin/bin/drawmap.R --mark $out/total  --out $out/total --pop cp \n";
	print SH2 "Rscript $Bin/bin/drawbinCP.R --mark $out/total.sexAver.phase  --out $out/total.sexAver.bin  \n";
	print SH2 "Rscript $Bin/bin/drawbinCP.R --mark $out/total.male.phase  --out $out/total.male.bin \n";
	print SH2 "Rscript $Bin/bin/drawbinCP.R --mark $out/total.female.phase  --out $out/total.female.bin \n";
	if ($ref) {
		print SH2 "Rscript $Bin/bin/drawAligmentRalationMap.pl -m $out/total.sexAver.map --out $out/ -k total.phy\n";
	}
}
close SH1;
close SH2;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step06.mapEvaluation1.sh";
	`$job`;
	$job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step06.mapEvaluation2.sh";
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -dmap	<dir>	input dmap dir
  -mark	<file>	input mark file
  -out	<dir>	output dir
  -dsh	<dir>	worksh dir
  -popt	<srt>	population type
  -h         Help

USAGE
        print $usage;
        exit;
}
