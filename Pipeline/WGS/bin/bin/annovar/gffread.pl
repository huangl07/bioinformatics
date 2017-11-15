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
my $genename="";
my $oldchr="";
my $genepos=-1;
my %out;
while (<In>) {
	chomp;
	next if ($_ eq  ""||/^$/ ||/#/);
	my @info=split(/\t/,$_);
	if ($oldchr ne $info[0] || $info[4] > $genepos) {
		$genename="";
		$oldchr=$info[0];
	}
	if ($info[2] eq "gene" || $info[2] eq "mRNA" || $info[2] eq "pseudoRNA"){
		$genepos=$info[4];
		if($info[8] =~ /Genebank:([^\"\'\;]*)/){
			$genename=$1;
		}elsif($info[8]=~/Name=([^\"\'\;]*)/) {
			$genename=$1;
		}elsif ( $info[8]=~/gene=([^\"\'\;]*)/) {
			$genename=$1;
		}elsif ($info[8]=~/geneID=([^\"\'\;]*)/) {
			$genename=$1;
		}elsif ($info[8] =~ /gene_name=([^\"\'\;]*)/) {
			$genename=$1;
		}elsif ($info[8] =~/ID=([^\"\'\;]*)/) {
			$genename.=":$1";
		}
	}
	if (/exon/ || /CDS/) {
		if ($genename eq ""){
			print STDOUT "no genename\t$_\n";
			next;
		};
		print Out join("\t",@info[0..7],"transcript_id \"$genename\" gene_id \"$genename\" gene_name \"$genename\""),"\n";
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
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
