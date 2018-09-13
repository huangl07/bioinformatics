#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($nr,$kegg,$go,$eggnog,$uniprot,$refdir,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"nr:s"=>\$nr,
	"kegg:s"=>\$kegg,
	"go:s"=>\$go,
	"eggnog:s"=>\$eggnog,
	"uniprot:s"=>\$uniprot,
	"ref:s"=>\$refdir,
	"out:s"=>\$out,
	) or &USAGE;
&USAGE unless ($nr and $kegg and $go and $eggnog and $uniprot and $out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
if (defined $refdir) {
	`ln -s $refdir/ref.fa $out/ref.fa`;
	`ln -s $refdir/ref.gff $out/ref.gff`;
	`ln -s $refdir/ref.predesign $out/ref.pre-design`;
}
`paste $nr $uniprot $kegg $go $eggnog|cut -f 1,2,3,5,6,8,9,11,12,14,15 > $out/anno.summary`;

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

Usage:
  Options:

	-nr	<file>	nr anno result
	-kegg	<file>	kegg anno result
	-go	<file>	go anno result
	-eggnog	<file>	eggnog anno result
	-uniprot	<file>	uniprot anno result
	-refdir	<dir>	ref dir
	-out	<output>	output dir

  -h         Help

USAGE
        print $usage;
        exit;
}
