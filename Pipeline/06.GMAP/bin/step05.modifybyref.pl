#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$only,$dsh,$cycle,$ref);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"lg:s"=>\$lg,
	"gen:s"=>\$gen,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dsh,
	"ref:s"=>\$ref,
			) or &USAGE;
&USAGE unless ($fIn and $dOut and $dsh);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
$cycle||=5;
open SH,">$dsh/step$cycle-0-2.split-calclodr.sh";
open In,$lg;
$/=">";
my %lg;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@info)=split(/\n/,$_);
	for (my $i=0;$i<@info;$i++) {
		my @marker=split(/\s+/,$info[$i]);
		foreach my $m (@marker) {
			$lg{$m}=$id;
		}
	}
}
close In;
$/="\n";
open In,$fIn;
my %lgpos;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,undef)=split(/\s+/,$_);
	my $lgID=$lg{$id};
	push @{$lgpos{$lgID}},$_;
}
close In;
foreach my $perl  (sort keys %lgpos) {
	my $only=200;
	my $split=int(scalar @{$lgpos{$lg}}/$only)+1;
	for (my $i=0;$i<$split;$i++) {
		open Out,">$out/sub/$lg.sub.$i.genotype";
		print Out $head;
		for (my $j= $only*$i;$j<@{$lgpos{$lg}};$j++) {
			print Out $lgpos{$lg}[$j],"\n";
		}
		close Out;
		if ($i == $split-1) {
			print SH "perl $Bin/bin/calculatePWDforCP.pl -i $out/sub/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/sub/ \n";
		}else{
			print SH "perl $Bin/bin/calculatePWDforCP.pl -i $out/sub/$lg.sub.$i.genotype -k $lg.sub.$i -d $out/sub/ -p $only\n";
		}
	}
}
close In;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-0-2.calculatePWD.sh";
`$job`;
close SH;
open SH,">$dsh/step03-2.merge-mlod.sh";
open List,">$out/marker.list";
foreach my $lg (sort keys %lg) {
	print SH "cat $out/sub/$lg.sub.*.pwd|less|sort|uniq > $out/$lg.pwd\n";
	print SH "cat $out/sub/$lg.sub.*.pwd.detail|less|sort|uniq > $out/$lg.pwd.detail \n";
	print List "$lg $out/$lg.pwd $out/$lg.pwd.detail\n";
	print SH "perl $Bin/bin/linkagePhase.pl -p $out/$lg.pwd -g $lg{$lg} -k $lg -d $out/ && ";
}
close SH;
close List;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dsh/step03-2.merge-mlod.sh";
print $job;
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
  -input	<file>	input file name
  -out	<dir>	output file dir
  -dsh	<dir>	output worksh dir
  -split	<num>	split number 
  -h         Help

USAGE
        print $usage;
        exit;
}
