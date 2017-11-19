#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$trt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"t:s"=>\$trt,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$trt;
my $Head;
my %Trt;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/SampleID/) {
		$Head=$_;
	}else{
		my ($id,@info)=split(/\s+/,$_);
		$Trt{$id}=join("\t",@info);
	}
}
close In;
open In,$fIn;
open Out,">$fOut.hapmap";
my %INDI;
my @Indi;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	s/\.trim//g;
	if (/rs/) {
		my ($snpid,$alleles,$chro,$pos,$strand,$assembly,$center,$protLSID,$assayLSID,$panel,$QCcode);
		($snpid,$alleles,$chro,$pos,$strand,$assembly,$center,$protLSID,$assayLSID,$panel,$QCcode,@Indi)=split(/\t/,$_);
		my @out;
		for (my $i=0;$i<@Indi;$i++) {
			next if (!exists $Trt{$Indi[$i]});
			$INDI{$Indi[$i]}=1; 
			push @out,$Indi[$i];
		}
		print Out join("\t",$snpid,$alleles,$chro,$pos,$strand,$assembly,$center,$protLSID,$assayLSID,$panel,$QCcode,@out),"\n";
	}else{
		my ($snpid,$alleles,$chro,$pos,$strand,$assembly,$center,$protLSID,$assayLSID,$panel,$QCcode,@indi)=split(/\t/,$_);
		my @out;
		for (my $i=0;$i<@indi;$i++) {
			next if (!exists $INDI{$Indi[$i]});
			push @out,$indi[$i];
		}
		print Out join("\t",$snpid,$alleles,$chro,$pos,$strand,$assembly,$center,$protLSID,$assayLSID,$panel,$QCcode,@out),"\n";
	}
}
close In;
close Out;
open Out,">$fOut.trt";
print Out $Head,"\n";
foreach my $id (sort keys %Trt) {
	next if (!exists $INDI{$id});
	print Out $id,"\t",$Trt{$id},"\n";
}
close Out;
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
