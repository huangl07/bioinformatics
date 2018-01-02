#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$out,$dsh,$pop,$mark,$ref);
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
	"popt:s"=>\$pop,
			) or &USAGE;
&USAGE unless ($dmap and $out and $dsh and $pop);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$mark=ABSOLUTE_DIR($mark);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$dmap=ABSOLUTE_DIR($dmap);
open SH1,">$dsh/step06.mapEvaluation1.sh";
open SH2,">$dsh/step06.mapEvaluation2.sh";
if ($pop ne "CP") {
	$pop=lc($pop);
	print SH1 "perl $Bin/bin/MapMergeNOCP.pl -dmap $dmap -o $out -adjust && ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.map -o $out/total.mapstat \n ";
	print SH2 "Rscript $Bin/bin/markerinfo.R --input $out/total.mapstat --output $out/fig/total \n";
	print SH2 "perl $Bin/bin/markerinfo.pl -map $out/total.map -input $out/total.marker --pop $pop -out $out/total.marker.info \n";
	print SH2 "Rscript $Bin/bin/drawmap.R --mark $out/total  --out $out/fig --pop $pop \n";
	print SH2 "Rscript $Bin/bin/drawbinNOCP.R --mark $out/total.csv  --out $out/fig/total.bin \n";
	if ($ref) {
		print SH2 "perl $Bin/bin/drawAligmentRalationMap.pl -m $out/total.map -o $out/fig/ -k total.phy\n";
	}
}else{
	print SH1 "perl $Bin/bin/MapMergeCP.pl -dmap $dmap -o $out --mark $mark -adjust && ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.sexAver.map -o $out/sexAver.mapstat \n ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.male.map -o $out/male.mapstat \n ";
	print SH2 "perl $Bin/bin/mapEstimate.pl -i $out/total.female.map -o $out/female.mapstat \n ";
	print SH2 "perl $Bin/bin/markerinfo.pl -map $out/total.sexAver.map -input $out/total.loc --pop $pop -out $out/total.sexAver.info \n";
	print SH2 "perl $Bin/bin/markerinfo.pl -map $out/total.male.map -input $out/total.loc --pop $pop -out $out/total.male.info \n";
	print SH2 "perl $Bin/bin/markerinfo.pl -map $out/total.female.map -input $out/total.loc --pop $pop -out $out/total.female.info \n";
	print SH2 "Rscript $Bin/bin/markerinfo.R --input $out/total.sexAver.mapstat --output $out/fig/total.sexAver \n";
	print SH2 "Rscript $Bin/bin/markerinfo.R --input $out/total.male.mapstat --output $out/fig/total.male \n ";
	print SH2 "Rscript $Bin/bin/markerinfo.R --input $out/total.female.mapstat --output $out/fig/total.female \n ";
	print SH2 "Rscript $Bin/bin/drawmap.R --mark $out/total  --out $out/fig --pop cp \n";
	print SH2 "Rscript $Bin/bin/drawbinCP.R --mark $out/total.sexAver.phase  --out $out/fig/total.sexAver.bin  \n";
	print SH2 "Rscript $Bin/bin/drawbinCP.R --mark $out/total.male.phase  --out $out/fig/total.male.bin \n";
	print SH2 "Rscript $Bin/bin/drawbinCP.R --mark $out/total.female.phase  --out $out/fig/total.female.bin \n";
	if ($ref) {
		print SH2 "Rscript $Bin/bin/drawAligmentRalationMap.pl -m $out/total.sexAver.map --out $out/fig/ -k total.phy\n";
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
