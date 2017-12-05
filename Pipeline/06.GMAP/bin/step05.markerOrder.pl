#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$out,$dsh,$popt,$cycle,$lg);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"lg:s"=>\$lg,
	"gen:s"=>\$fIn,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"popt:s"=>\$popt,
	"cycle:s"=>\$cycle,
			) or &USAGE;
&USAGE unless ($fIn and $out and $dsh and $popt);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$fIn=ABSOLUTE_DIR($fIn);
if ($lg) {
	$lg=ABSOLUTE_DIR($lg);
	open SH,">$dsh/step05-0.split.sh";
	if ($popt eq "CP") {
		print SH "perl $Bin/bin/splitbyLG-CP.pl -l $lg -i $fIn -d $out/ -t $popt";
	}else{
		print SH "perl $Bin/bin/splitbyLG-NOCP.pl -l $lg -i $fIn -d $out/ -t $popt";
	}
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-0.split.sh";
	`$job`;
	$fIn=ABSOLUTE_DIR("$out/pri.marker.list");
}
$cycle||=1;
if ($popt eq "CP") {
	mkdir "$out/pwd" if (!-d "$out/pwd");
	open In,$fIn;
	open SH,">$dsh/step05-$cycle.1.calculatePWD.sh";
	my %lg;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($lg,$file)=split(/\s+/,$_);
		die $file if (!-f $file);
		$lg{$lg}=$file;
		my $head;
		my @Marker;
		my $read;
		open Marker,$file;
		while ($read=<Marker>) {
			chomp;
			next if ($read eq ""||$read =~ /^$/);
			if ($read =~ /^#/) {
				$head=$read;
			}else{
				push @Marker,$read;
			}
		}
		close Marker;
		my $only=200;
		my $split=int(scalar @Marker/$only)+1;
		for (my $i=0;$i<$split;$i++) {
			open Out,">$out/pwd/$lg.sub.$i.genotype";
			print Out $head;
			for (my $j= $only*$i;$j<@Marker;$j++) {
				print Out $Marker[$j],"\n";
			}
			close Out;
			if ($i == $split-1) {
				print SH "perl $Bin/calculatePWDforCP.pl -i $out/pwd/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/pwd/ \n";
			}else{
				print SH "perl $Bin/calculatePWDforCP.pl -i $out/pwd/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/pwd/ -s $only\n";
			}
		}
	}
	close In;
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-$cycle.1.calculatePWD.sh";
	`$job`;
	open SH,">$dsh/step05-$cycle.2.mapping.sh";
	open List,">$out/marker.list";
	foreach my $lg (sort keys %lg) {
		print SH "cat $out/pwd/$lg.sub.*.pwd > $out/pwd/$lg.pwd && cat $out/pwd/$lg.sub.*.pwd.detail > $out/pwd/$lg.pwd.detail && ";
		print SH "perl $Bin/bin/linkagePhase.pl -p $out/pwd/$lg.pwd -g $lg{$lg} -k $lg -d $out/ && ";
		print SH "perl $Bin/bin/extractPwdViaLP.pl -i  $out/pwd/$lg.pwd.detail -l $out/$lg.loc -k $lg -d $out && ";
		print SH "sgsMap -loc $out/$lg.loc -pwd $out/$lg.pwd -k $out/$lg &&";
		print SH "perl $Bin/bin/smooth-CP.pl -m $out/$lg.sexAver.map -l $out/$lg.loc -k $lg -d $out\n";
		print List $lg,"\t","$out/$lg.correct.loc\n";
	}
	close List;
	close SH;
	 $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-$cycle.2.mapping.sh";
	`$job`;
}else{
	open SH,">$dsh/step05-$cycle.markerOrder.sh";
	open List,">$out/marker.list";
	open In,$fIn;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($lg,$file)=split(/\s+/,$_);
		print SH "MSTMap $file $out/$lg.out &&";
		print SH "perl $Bin/bin/smooth-NOCP.pl -i $file -m $out/$lg.out -o $out/$lg.marker \n ";
		print List "$lg\t$out/$lg.marker\n";
	}
	close SH;
	close In;
	close List;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-$cycle.markerOrder.sh";
	print $job;
	`$job`;
}





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
  -gen	<file>	input gen file
  -out	<dir>	output dir
  -dsh	<dir>	worksh dir
  -popt	<srt>	population type
  -h         Help

USAGE
        print $usage;
        exit;
}
