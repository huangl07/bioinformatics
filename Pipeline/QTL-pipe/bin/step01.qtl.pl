#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$pop,$btl,$trit,$dmap,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"trit:s"=>\$trit,
	"out:s"=>\$out,
	"popt:s"=>\$pop,
	"btl:s"=>\$btl,
	"dsh:s"=>\$dsh,
	) or &USAGE;
&USAGE unless ($out and $dmap and $trit and $pop);

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:
        fq thanslate to fa format
        eg:
        perl $Script -i -o -k -c
Usage:
  Options:
  -dmap	<dir>	input the gmap's result dir(*.csv or *.phase and *.map )
  -trit	<file>	input the qtl's trit list file
  -out	<dir>	output dir
  -popt	<str>	pop type
  -btl		binary trt or not
  -dsh	<dir>	worksh dir
        
  -h         Help

USAGE
        print $usage;
        exit;
}
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$trit=ABSOLUTE_DIR($trit);
#mkdir "$out/vcf" if (!-d "$out/vcf");
#mkdir "$out/gene" if (!-d "$out/gene");

#open SH,">$dsh/step01.qtl.sh";
if ($btl) {
	if ($pop=~/CP/){
		mkdir "$out/male" if (!-d "$out/male");
		mkdir "$out/female" if (!-d "$out/female");
		mkdir "$out/sexAver" if (!-d "$out/sexAver");
		open In,$trit;
		open SH1,">$dsh/step01.qtl.male.sh";
		open SH2,">$dsh/step01.qtl.female.sh";
		open SH3,">$dsh/step01.qtl.sexAver.sh";
		while(<In>){
			chomp;
			next if ($_ eq ""|| /^$/);
			#my ($tritname,$tritfile)=split(/\t/,$_);
			my $tritfile=$_;
			my @mark=glob("$dmap/total*.phase");
			for(@mark){
				if($_=~/total.male.phase/){
					print SH1 "Rscript $Bin/bin/btl-CP.R --map $dmap/total.male.map --loc $dmap/total.male.phase --trt $tritfile --out $out/male --num 1000 &&\n";
				}elsif($_=~/total.female.phase/){
					print SH2 "Rscript $Bin/bin/btl-CP.R --map $dmap/total.female.map  --loc $dmap/total.female.phase --trt $tritfile --out $out/female --num 1000 &&\n";
				}elsif($_=~/total.sexAver.phase/){
					print SH3 "Rscript $Bin/bin/btl-CP.R --map $dmap/total.sexAver.map  --loc $dmap/total.sexAver.phase --trt $tritfile --out $out/sexAver --num 1000 &&\n";
				}
			}
		}
		close In;
		close SH1;
		close SH2;
		close SH3;
		my $sh="cat $dsh/step01.qtl.male.sh $dsh/step01.qtl.female.sh $dsh/step01.qtl.sexAver.sh >$dsh/step01.sh";
		`$sh \n`; 
		my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=40G --CPU 1  --Maxjob 20 $dsh/step01.sh \n";
		`$job \n`;
		print "job done \n";
	}elsif($pop=~/F2/){
		open In,$trit;
		open SH,">$dsh/step01.qtl.sh";
		while(<In>){
			chomp;
			next if ($_ eq ""|| /^$/);
			my $tritfile=$_;
			print SH "Rscript $Bin/bin/btl-NOCP.R --mark $dmap/total.csv --trt $tritfile --pop bcsft --bc 2 --f 2 --out  $out --num 1000 && \n";
		}
		close In;
		close SH;
		my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=40G --CPU 1  --Maxjob 20 $dsh/step01.qtl.sh";
		`$job \n`;
	}
}else{
	if($pop=~/F2/){
		open In,$trit;
		open SH,">$dsh/step01.qtl.sh";
		while(<In>){
			chomp;
			next if ($_ eq ""|| /^$/);
			#my ($tritname,$tritfile)=split(/\t/,$_);
			my $tritfile=$_;
			print SH "Rscript $Bin/bin/qtl-NOCP.R --mark $dmap/total.csv --trt $tritfile --out $out --pop bcsft --num 1000  --bc 2 --f 2 && \n";
		}
		close In;
		close SH;
		my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=40G --CPU 1 --Maxjob 20 $dsh/step01.qtl.sh";
		`$job\n`;
		print "job done\n";
	}elsif($pop=~/CP/){
		mkdir "$out/male" if (!-d "$out/male");
		mkdir "$out/female" if (!-d "$out/female");
		mkdir "$out/sexAver" if (!-d "$out/sexAver");
		open In,$trit;
		open SH1,">$dsh/step01.qtl.male.sh";
		open SH2,">$dsh/step01.qtl.female.sh";
		open SH3,">$dsh/step01.qtl.sexAver.sh";
		while(<In>){
			chomp;
			next if ($_ eq ""|| /^$/);
			#my ($tritname,$tritfile)=split(/\t/,$_);
			my $tritfile=$_;
			my @mark=glob("$dmap/total*.phase");
			for(@mark){
				if($_=~/total.male.phase/){
					print SH1 "Rscript $Bin/bin/qtl-CP.R --map $dmap/total.male.map --loc $dmap/total.male.phase --trt $tritfile --out $out/male --num 1000 &&\n";
				}elsif($_=~/total.female.phase/){
					print SH2 "Rscript $Bin/bin/btl-CP.R --map $dmap/total.female.map  --loc $dmap/total.female.phase --trt $tritfile --out $out/female --num 1000 &&\n";
				}elsif($_=~/total.sexAver.phase/){
					print SH3 "Rscript $Bin/bin/btl-CP.R --map $dmap/total.sexAver.map  --loc $dmap/total.sexAver.phase --trt $tritfile --out $out/sexAver --num 1000 &&\n";
				}
			}
		}
		close In;
		close SH1;
		close SH2;
		close SH3; 
		my $sh="cat $dsh/step01.qtl.male.sh $dsh/step01.qtl.female.sh $dsh/step01.qtl.sexAver.sh >$dsh/step01.sh";
		`$sh\n`;
		my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=40G --CPU 1 --Maxjob 20 $dsh/step01.sh";
		`$job\n`;
		print "job done\n";
	}
}
#my @jobs=glob("$dsh/step01.*sh");
#for(@jobs){
#	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=40G --CPU 1  --Maxjob 20 $_ \n";
#	`$job\n`;
#	print "done\n";
#}
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
