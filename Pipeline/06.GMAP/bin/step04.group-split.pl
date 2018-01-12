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
open SH,">$dsh/step04-1.grouping.sh";
my $lg="$out/Total.lg";
if ($ref) {
	print SH "perl $Bin/bin/linkage_by_ref.pl -i $mlod -2 $marker -o $out -k Total -add && ";
	print SH "perl $Bin/bin/splitbyLG-CP.pl -l $lg -i $marker -d $out/ -t $popt ";
}else{
	print SH "perl $Bin/bin/linkage_by_mlod.pl -i $mlod -k Total -d $out -n $nchr && ";
	print SH "perl $Bin/bin/splitbyLG-NOCP.pl -l $lg -i $marker -d $out/ -t $popt";
}
close SH;
my $mem=`du $mlod`;
chomp $mem;
$mem=(split(/\s+/,$mem))[0];
$mem=int($mem/1000000)+3;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=$mem"."G --CPU 1 $dsh/step04-1.grouping.sh";
print $job;
`$job`;
if ($ref) {
	if ($popt eq "CP") {
		open In,"$out/pri.marker.list";
		open SH,">$dsh/step04-2.calculatePWD.sh";
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
		close SH;
		close In;
		my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-0-2.calculatePWD.sh";
		`$job`;
		open SH,">$dsh/step04-3.ref.sh";
		open List,">$out/marker.list";
		foreach my $lg (sort keys %lg) {
			print SH "cat $out/pri-pwd/$lg.sub.*.pwd|less|sort|uniq > $out/pri-pwd/$lg.pri.pwd && cat $out/pri-pwd/$lg.sub.*.pwd.detail|less|sort|uniq > $out/pri-pwd/$lg.pri.pwd.detail && ";
			print SH "perl $Bin/bin/linkagePhase.pl -p $out/pri-pwd/$lg.pri.pwd -g $lg{$lg} -k $lg.pri -d $out/ && ";
			print SH "perl $Bin/bin/smooth-CP.pl -m $out/$lg.pri.map -l $out/$lg.pri.loc -k $lg.pri -d $out -diff_ratio 0.7 \n";
			print List $lg,"\t","$out/$lg.pri.correct.loc\n";
		}
		close List;
		close SH;
		$job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step04-3.ref.sh";
		`$job`;
	}else{
		open SH,">$dsh/step05-0-2.ref.sh";
		open In,"$out/pri.marker.list";
		open Out,">$out/marker.list";
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
	}
}else{
	`ln -s $out/pri.marker.list $out/marker.list`;
}
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
	-ref			split by ref or not
	-popt	<str>	population type
	-out	<out>	output dir
	-dsh	<dir>	output work shell file
USAGE
	print $usage;
	exit;
}
