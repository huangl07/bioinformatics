#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$out,$dsh,$popt,$cycle,$lg,$ref);
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
	"ref"=>\$ref,
	"cycle:s"=>\$cycle,
			) or &USAGE;
&USAGE unless ($fIn and $out and $dsh and $popt);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$fIn=ABSOLUTE_DIR($fIn);
mkdir "$out/pri-pwd" if (!-d "$out/pri-pwd");
mkdir "$out/pwd" if (!-d "$out/pwd");

$cycle||=1;
if ($lg) {
	$lg=ABSOLUTE_DIR($lg);
	open SH,">$dsh/step05-0-1.split.sh";
	if ($popt eq "CP") {
		print SH "perl $Bin/bin/splitbyLG-CP.pl -l $lg -i $fIn -d $out/ -t $popt ";
	}else{
		print SH "perl $Bin/bin/splitbyLG-NOCP.pl -l $lg -i $fIn -d $out/ -t $popt";
	}
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-0-1.split.sh";
	`$job`;
	$fIn=ABSOLUTE_DIR("$out/pri.marker.list");
	if ($ref) {
		if ($popt eq "CP") {
			open SH,">$dsh/step05-0-2.calculatePWD.sh";
			open In,$fIn;
			my %lg;
			while (<In>) {
				chomp;
				next if ($_ eq ""||/^$/);
				my ($lg,$file)=split(/\s+/,$_);
				die $file if (!-f $file);
				`ln -s $file $out/$lg.pri.loc`;
				$lg{$lg}=$file;
				my $head;
				my @Marker;
				my $read="";
				open Marker,$file;
				while ($read=<Marker>) {
					chomp;
					next if ($read eq ""||$read =~ /^$/ || /^;/);
					if ($read =~ /^#/ || $read =~/=/ ) {
						$head.=$read."\n";
					}else{
						push @Marker,$read;
					}
				}
				close Marker;
				my $only=200;
				my $split=int(scalar @Marker/$only)+1;
				for (my $i=0;$i<$split;$i++) {
					open Out,">$out/pri-pwd/$lg.sub.$i.genotype";
					print Out $head;
					for (my $j= $only*$i;$j<@Marker;$j++) {
						print Out $Marker[$j],"\n";
					}
					close Out;
					if ($i == $split-1) {
						print SH "perl $Bin/bin/calculatePWDforCP.pl -i $out/pri-pwd/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/pri-pwd/ \n";
					}else{
						print SH "perl $Bin/bin/calculatePWDforCP.pl -i $out/pri-pwd/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/pri-pwd/ -p $only\n";
					}
				}
			}
			close In;
			close SH;
			my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-0-2.calculatePWD.sh";
			`$job`;
			open SH,">$dsh/step05-0-3.ref.sh";
			open List,">$out/ref.marker.list";
			foreach my $lg (sort keys %lg) {
				print SH "cat $out/pri-pwd/$lg.sub.*.pwd|less|sort|uniq > $out/pri-pwd/$lg.pri.pwd && cat $out/pri-pwd/$lg.sub.*.pwd.detail|less|sort|uniq > $out/pri-pwd/$lg.pri.pwd.detail && ";
				my $loc="$out/$lg.pri.loc";
				print SH "perl $Bin/bin/linkagePhase.pl -p $out/pri-pwd/$lg.pri.pwd -g $lg{$lg} -k $lg.pri -d $out/ && ";
				print SH "perl $Bin/bin/extractPwdViaLP.pl -i  $out/pri-pwd/$lg.pri.pwd.detail -l $loc -k $lg.pri -d $out && ";
				print SH "perl $Bin/bin/smooth-CP.pl -m $out/$lg.pri.map -l $loc -k $lg.pri -d $out\n";
				print List $lg,"\t","$out/$lg.pri.correct.loc\n";
			}
			close List;
			close SH;
			$job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-0-3.ref.sh";
			`$job`;
			$fIn=ABSOLUTE_DIR("$out/ref.marker.list");
		}else{
			open SH,">$dsh/step05-0-2.ref.sh";
			open In,$fIn;
			open Out,">$out/ref.marker.list";
			while (<In>) {
				chomp;
				next if ($_ eq ""||/^$/);
				my ($lg,$file)=split(/\s+/,$_);
				print Out join("\t",$lg,"$out/$lg.cor.marker"),"\n";
				print SH "perl $Bin/bin/smooth-NOCP.pl -i $file -m $out/$lg.pri.map -o $out/$lg.cor.marker \n ";
			}
			close Out;
			close In;
			close SH;
			my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-0-2.ref.sh";
			`$job`;
			$fIn=ABSOLUTE_DIR("$out/ref.marker.list");
		}
	}
}
$cycle||=1;
if ($popt eq "CP") {
	open In,$fIn;
	open SH,">$dsh/step05-$cycle.1.calculatePWD.sh";
	my %lg;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($lg,$file)=split(/\s+/,$_);
		die $file if (!-f $file);
		`ln -s $file $out/$lg.pri.loc`;
		$lg{$lg}=$file;
		my $head;
		my @Marker;
		my $read="";
		open Marker,$file;
		while ($read=<Marker>) {
			chomp;
			next if ($read eq ""||$read =~ /^$/ || /^;/);
			if ($read =~ /^#/ || $read =~/=/ ) {
				$head.=$read."\n";
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
				print SH "perl $Bin/bin/calculatePWDforCP.pl -i $out/pwd/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/pwd/ \n";
			}else{
				print SH "perl $Bin/bin/calculatePWDforCP.pl -i $out/pwd/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/pwd/ -p $only\n";
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
		print SH "cat $out/pwd/$lg.sub.*.pwd|less|sort|uniq > $out/pwd/$lg.pwd && cat $out/pwd/$lg.sub.*.pwd.detail|less|sort|uniq > $out/pwd/$lg.pwd.detail && ";
		my $loc="$out/$lg.loc";
		if ($cycle == 1) {
			print SH "perl $Bin/bin/linkagePhase.pl -p $out/pwd/$lg.pwd -g $lg{$lg} -k $lg -d $out/ && ";
		}else{
			$loc="$out/$lg.pri.loc"
		}
		print SH "perl $Bin/bin/extractPwdViaLP.pl -i  $out/pwd/$lg.pwd.detail -l $loc -k $lg -d $out && ";
		print SH "sgsMap -loc $loc -pwd $out/$lg.pwd -k $out/$lg &&";
		print SH "perl $Bin/bin/smooth-CP.pl -m $out/$lg.sexAver.map -l $loc -k $lg -d $out\n";
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
		print SH "MSTmap $file $out/$lg.out &&";
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
  -lg	<file>	input lg file
  -gen	<file>	input gen file
  -out	<dir>	output dir
  -dsh	<dir>	worksh dir
  -popt	<srt>	population type
  -ref mapping by ref
  -h         Help

USAGE
        print $usage;
        exit;
}
