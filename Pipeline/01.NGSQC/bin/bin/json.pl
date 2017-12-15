#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use JSON;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
my $json=new JSON;
my $js;
my $sample=(split(/\./,(split(/\//,$fIn))[-1]))[0];
open In,$fIn;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	$js.=$_."\n";
}
close In;
my $obj = $json->decode($js);

my @b1rawA=@{$obj->{"read1_before_filtering"}->{"content_curves"}->{"A"}};
my @b1rawT=@{$obj->{"read1_before_filtering"}->{"content_curves"}->{"T"}};
my @b1rawC=@{$obj->{"read1_before_filtering"}->{"content_curves"}->{"C"}};
my @b1rawG=@{$obj->{"read1_before_filtering"}->{"content_curves"}->{"G"}};
my @b1rawN=@{$obj->{"read1_before_filtering"}->{"content_curves"}->{"N"}};

my @b1cleanA=@{$obj->{"read1_after_filtering"}->{"content_curves"}->{"A"}};
my @b1cleanT=@{$obj->{"read1_after_filtering"}->{"content_curves"}->{"T"}};
my @b1cleanC=@{$obj->{"read1_after_filtering"}->{"content_curves"}->{"C"}};
my @b1cleanG=@{$obj->{"read1_after_filtering"}->{"content_curves"}->{"G"}};
my @b1cleanN=@{$obj->{"read1_after_filtering"}->{"content_curves"}->{"N"}};

my @b2rawA=@{$obj->{"read2_before_filtering"}->{"content_curves"}->{"A"}};
my @b2rawT=@{$obj->{"read2_before_filtering"}->{"content_curves"}->{"T"}};
my @b2rawC=@{$obj->{"read2_before_filtering"}->{"content_curves"}->{"C"}};
my @b2rawG=@{$obj->{"read2_before_filtering"}->{"content_curves"}->{"G"}};
my @b2rawN=@{$obj->{"read2_before_filtering"}->{"content_curves"}->{"N"}};

my @b2cleanA=@{$obj->{"read2_after_filtering"}->{"content_curves"}->{"A"}};
my @b2cleanT=@{$obj->{"read2_after_filtering"}->{"content_curves"}->{"T"}};
my @b2cleanC=@{$obj->{"read2_after_filtering"}->{"content_curves"}->{"C"}};
my @b2cleanG=@{$obj->{"read2_after_filtering"}->{"content_curves"}->{"G"}};
my @b2cleanN=@{$obj->{"read2_after_filtering"}->{"content_curves"}->{"N"}};


my @b1rawQ=@{$obj->{"read1_before_filtering"}->{"quality_curves"}->{"mean"}};
my @b2rawQ=@{$obj->{"read2_before_filtering"}->{"quality_curves"}->{"mean"}};


my @b1cleanQ=@{$obj->{"read1_after_filtering"}->{"quality_curves"}->{"mean"}};
my @b2cleanQ=@{$obj->{"read2_after_filtering"}->{"quality_curves"}->{"mean"}};


open Out,">$fOut.stat";
print Out join("\t","#sampleID","rawreads","rawdata","rawq20","rawq30","rawGC","ada","cleanreads","cleandata","cleanq20","cleanq30","cleanGC"),"\n";
my @out;
push @out,$sample;
push @out,$obj->{"summary"}->{"before_filtering"}->{"total_reads"};
push @out,$obj->{"summary"}->{"before_filtering"}->{"total_bases"};
push @out,$obj->{"summary"}->{"before_filtering"}->{"q20_rate"};
push @out,$obj->{"summary"}->{"before_filtering"}->{"q30_rate"};
push @out,average(\@b1rawC,\@b1rawG,\@b2rawC,\@b2rawG);
push @out,$obj->{"adapter_cutting"}->{"adapter_trimmed_reads"}/$obj->{"summary"}->{"before_filtering"}->{"total_reads"};
push @out,$obj->{"summary"}->{"after_filtering"}->{"total_reads"};
push @out,$obj->{"summary"}->{"after_filtering"}->{"total_bases"};
push @out,$obj->{"summary"}->{"after_filtering"}->{"q20_rate"};
push @out,$obj->{"summary"}->{"after_filtering"}->{"q30_rate"};
push @out,average(\@b1cleanG,\@b1cleanC,\@b2cleanG,\@b2cleanG);
print Out join("\t",@out),"\n";
close Out;
open Out,">$fOut.raw.atgcn\n";
for (my $i=0;$i<@b1rawA;$i++) {
	print Out join("\t",$i,$b1rawA[$i],$b1rawT[$i],$b1rawG[$i],$b1rawC[$i],$b1rawN[$i]),"\n";
}
for (my $i=0;$i<@b2rawA;$i++) {
	print Out join("\t",$i+scalar @b1rawA,$b2rawA[$i],$b2rawT[$i],$b2rawG[$i],$b2rawC[$i],$b2rawN[$i]),"\n";
}
close Out;
open Out,">$fOut.clean.atgcn\n";
for (my $i=0;$i<@b1cleanA;$i++) {
	print Out join("\t",$i,$b1cleanA[$i],$b1cleanT[$i],$b1cleanG[$i],$b1cleanC[$i],$b1cleanN[$i]),"\n";
}
for (my $i=0;$i<@b2cleanA;$i++) {
	print Out join("\t",$i+scalar @b1cleanA,$b2cleanA[$i],$b2cleanT[$i],$b2cleanG[$i],$b2cleanC[$i],$b2cleanN[$i]),"\n";
}
close Out;
open Out,">$fOut.raw.qual\n";
for (my $i=0;$i<@b1rawA;$i++) {
	print Out join("\t",$i,$b1rawQ[$i]),"\n";
}
for (my $i=0;$i<@b2rawA;$i++) {
	print Out join("\t",$i,$b2rawQ[$i]),"\n";
}
close Out;
open Out,">$fOut.clean.quan\n";
for (my $i=0;$i<@b1cleanA;$i++) {
	print Out join("\t",$i,$b1cleanQ[$i]),"\n";
}
for (my $i=0;$i<@b2cleanA;$i++) {
	print Out join("\t",$i,$b2cleanQ[$i]),"\n";
}
close Out;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub average{
	my ($a,$b,$c,$d)=@_;
	my $sum=0;
	my $n=0;
	for (my $i=0;$i<@{$a};$i++) {
		$sum+=$$a[$i];
		$sum+=$$b[$i];
		$n++;
	}
	for (my $i=0;$i < @{$c} ;$i++) {
		$sum+=$$c[$i];
		$sum+=$$d[$i];
		$n++;
	}
	return($sum/$n);
}
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

Usage:
  Options:
  -i	<file>	input dict name
  -qc	<file>	out qc file name
  -qual	<file>	out qual file name
  -base	<file>	out base file name
  -h         Help

USAGE
        print $usage;
        exit;
}
