#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($Afile,$Bfile,$Cfile);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"inputA:s"=>\$Afile
	"inputB:s"=>\$Bfile,
	"outputC:s"=>\$Cfile
			) or &USAGE;
&USAGE unless ($Afilte and $Bfile and $Cfile);
open In,$Afile;
my %info;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my $id=(split(/\t/,$_))[0];
	my @ids=split(/\|/,$id);
	my $ids=join("|",$id[0],$id[1],$id[2],$id[3]);
	$info{$ids}++;
}
close In
$/=">"
open In,$Bfile;
print Out,">$Cfile";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@line)=split(/\n/,$_);
	my @ids=split(/\|/,$id);
	my $ids=join("|",$id[0],$id[1],$id[2],$id[3]);
	if (exists $info{$ids}) {
		print Out join(">$ids\n").join("\n",@line),"\n";
	}
}
close In;
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
  -A	<file>	input A file name
  -B	<file>	input B file name
  -C	<file>	output C file
  -h         Help

USAGE
        print $usage;
        exit;
}
