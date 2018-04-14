#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$list,$vcf,$marker,$out,$dsh,$pop,$type);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"list:s"=>\$list,
	"vcf:s"=>\$vcf,
	"marker:s"=>\$marker,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"popt:s"=>\$pop,
	"type:s"=>\$type,
			) or &USAGE;
&USAGE unless ($dmap and $list and $vcf and $marker and $out and $dsh and $pop);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:                 $Script
Description:
        fq thanslate to fa format
        eg:
        perl $Script -i -o -k -c

Usage:
  Options:
  -dmap <dir>   input dmap dir(10.mapEvalue)
  -list	<file>	input ref.chrlist
  -vcf	<file>	input pop.final.vcf
  -marker <dir>  input 01.vcf-convert dir
  -out  <dir>   output dir
  -dsh  <dir>   worksh dir
  -popt <srt>   population type
  -type	<str>	chr or sca
  -h         Help

USAGE
        print $usage;
        exit;
}
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
#$marker=ABSOLUTE_DIR($marker);
#$out=ABSOLUTE_DIR($out);
#$dsh=ABSOLUTE_DIR($dsh);
#$map=ABSOLUTE_DIR($map);
mkdir "$out/Figure" if (!-d "$out/Figure");
mkdir "$out/Result" if (!-d "$out/Result");
mkdir "$out/Table" if (!-d "$out/Table");
$type||="chr";
open SH1,">$dsh/step08.merge-result.sh";
if ($pop ne "CP") {
	#system("cp $dmap/total.marker.info $out/Result/ \n");
	open In,"$dmap/total.marker.info";
	open Out,">$out/Result/total.marker.info.stat";
	print Out"#MarkerID\tGroupID\tDis\ttype\tNind\tNmiss\tGeno\tNGeno\tSegretion\n";
	while(<In>){
		chomp;
		next if($_=~/^#/|| /^$/);
		my($mak,$group,$dis,$type,$nind,$nmiss,$geno,$ngeno,$pv,$seg)=split/\t/,$_;
		my $info=join("\t",$mak,$group,$dis,$type,$nind,$nmiss,$geno,$ngeno,$seg);
		print Out $info,"\n";
	}
	close In;
	close Out;
	system("cp $dmap/total.marker $out/Result/ \n");
	system("cp $dmap/fig/total.bin.* $out/Figure/ \n");
	#system("cp $dmap/fig/total.map.* $out/Figure/ \n");
	system("cp $dmap/fig/total.phy.* $out/Figure/ \n");
	system("cp $dmap/fig/total.lg.* $out/Figure/ \n");
	print SH1 "perl $Bin/bin/merge-marker.vcf.pl -vcf $vcf -mark $marker/pop.primary.marker -out $out/Table/3-12.xls\n";
	print SH1 "perl $Bin/bin/merge-marker_inlgs_type.pl -type $type -vcf $vcf -mark $marker/pop.filtered.marker -out $out/Table/3-13.xls\n";
	open In,"$dmap/total.mapstat";
	open Out,">$out/Table/3-14.xls";
	print Out "#LGID\tMarker Number\tTotal Distance\tAvarage Distance\tGap < 5cM(%)\tMax Gap\n";
	while (<In>){
		chomp;
		next if ($_=~/^#/|| /^$/);
		my($lg,$marknum,$uniqnum,$totaldis,$averdis,$gapnum,$maxgap)=split(/\t/,$_);
		my $gapratio=sprintf("%.2f",100*($uniqnum - $gapnum)/$uniqnum);
		print Out "$lg\t$uniqnum\t$totaldis\t$averdis\t$gapratio\t$maxgap\n";
	}
	close In;
	close Out;
	print SH1 "perl $Bin/bin/merge-wuliandyichuan.pl -list $list -mark $dmap/total.map -out $out/Table/3-15.xls \n";
}else{
	my @makinfo=glob("$dmap/total.*info.marker.stat");
	for(@makinfo){
		open In,"$_";
		if($_=~/total.male.info/){
			#open In,"$_";
			open Out,">$out/Result/total.male.marker.info.stat";
		}elsif($_=~/total.female.info/){
			#open In,"$_";
			open Out,">$out/Result/total.female.marker.info.stat";
		}elsif($_=~/total.sexAver.info/){
			open Out,">$out/Result/total.sexAver.marker.info.stat";
		}
		print Out "#MarkerID\tGroupID\tDis\ttype\tNind\tNmiss\tGeno\tNGeno\tSegretion\n";
		while(<In>){
			chomp;
			next if ($_=~/^#/|| /^$/);
			my($mak,$group,$dis,$type,$nind,$nmiss,$geno,$ngeno,$pv,$seg)=split/\t/,$_;
			my $info=join("\t",$mak,$group,$dis,$type,$nind,$nmiss,$geno,$ngeno,$seg);
			print Out $info,"\n";
		}
		close In;
		close Out;
	}
	system("cp $dmap/total.qtl $out/Result/ \n");
	system("cp $dmap/fig/total.female.phy.* $out/Figure/ \n");
	system("cp $dmap/fig/total.male.phy.* $out/Figure/ \n");
	system("cp $dmap/fig/total.sexAver.phy.* $out/Figure/ \n");
	system("cp $dmap/fig/total.female.lg.* $out/Figure/ \n");
	system("cp $dmap/fig/total.male.lg.* $out/Figure/ \n");
	system("cp $dmap/fig/total.sexAver.lg.* $out/Figure/ \n");
	system("cp $dmap/fig/total.female.bin.* $out/Figure/ \n");
	system("cp $dmap/fig/total.male.bin.* $out/Figure/ \n");
	system("cp $dmap/fig/total.sexAver.bin.* $out/Figure/ \n");
	
	print SH1 "perl $Bin/bin/merge-marker.vcf.pl -vcf $vcf -mark $marker/pop.primary.marker -out $out/Table/3-12.xls \n";
	print SH1 "perl $Bin/bin/merge-marker_inlgs_type.pl -type $type -vcf $vcf -mark $marker/pop.filtered.marker -out $out/Table/3-13.xls\n";
	my @mapstat=glob("$dmap/total.*mapstat");
	for(@mapstat){
		if($_=~/total.male.mapstat/){
			open In,"$_";
			open Out,">$out/Table/3-14.xls";
		}elsif($_=~/total.female.mapstat/){
			open In,"$_";
			open Out,">$out/Table/3-15.xls";
		}elsif($_=~/total.sexAver.mapstat/){
			open In,"$_";
			open Out,">$out/Table/3-16.xls";
		}
		print Out "#LGID\tMarker Number\tTotal Distance\tAvarage Distance\tGap < 5cM(%)\tMax Gap\n";
		while (<In>){
			chomp;
			next if ($_=~/^#/|| /^$/);
			my($lg,$marknum,$uniqnum,$totaldis,$averdis,$gapnum,$maxgap)=split(/\t/,$_);
			my $gapratio=sprintf("%.2f",100*($uniqnum - $gapnum)/$uniqnum);
			print Out "$lg\t$uniqnum\t$totaldis\t$averdis\t$gapratio\t$maxgap\n";
		}
		close In;
		close Out;
	}
	print SH1 "perl $Bin/bin/merge-wuliandyichuan.pl -list $list -mark $dmap/total.sexAver.map -out $out/Table/3-17.xls\n";
}
close SH1;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=12G --CPU 1 $dsh/step08.merge-result.sh";
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
