#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gtf,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$gtf,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($gtf and $fOut);
open In,$gtf;
open Out,">$fOut";
my %id;
while (<In>) {
	chomp;
	next if ($_ eq  ""||/^$/);
	my @info=split(/\t/,$_);
	my $geneid="";
	my $genename="";
	my $transid="";
	if ($info[-1]=~/gene_name "([^\"\']*)";/) {
		$genename=$1;
	}
	if($info[-1]=~/gene_id "([^\"\']*)";/) {
		$geneid=$1;
	}
	if ($info[-1]=~/transcript_id "([^\"\']*)";/) {
		$transid=$1;
	}
	my $id;
	if ( $genename ne "" && $geneid eq ""  &&  $transid eq "") {
		$id="$genename.1";
		print Out join("\t",@info[0..7],"transcript_id \"$id\"; gene_id \"$id\"; gene_name \"$id\""),"\n";
	}elsif ( $geneid ne "" && $genename eq "" &&  $transid eq "") {
		$id="$geneid.1";
		print Out join("\t",@info[0..7],"transcript_id \"$id\"; gene_id \"$id\"; gene_name \"$id\""),"\n";
	}elsif ( $transid ne "" &&  $genename eq "" &&  $geneid eq "") {
		$id="$transid.1";
		print Out join("\t",@info[0..7],"transcript_id \"$id\"; gene_id \"$id\"; gene_name \"$id\""),"\n";
	}elsif ( $transid ne "" &&  $genename ne "" &&  $geneid ne "") {
		my $key=join("\t",$transid,$genename,$geneid);
		$id{$key}=$genename;
		$id="$genename.".scalar keys %id;
		print Out join("\t",@info[0..7],"transcript_id \"$id\"; gene_id \"$id\"; gene_name \"$id\""),"\n";
	}else{
		print $_;die;
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
