#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fa,$out,$dsh,$type);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fa:s"=>\$fa,
	"type:s"=>\$type,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	) or &USAGE;
&USAGE unless ($out and $dsh and $fa);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$type||="nuc";
$dsh=ABSOLUTE_DIR($dsh);
my $dGO="/mnt/ilustre/users/long.huang/DataBase/2018-8-27/GO/go_weekly-seqdb.fasta";###
my $blast="diamond";
if ($type eq "pro") {
	$blast=" blastp";
}else{
	$blast=" blastx";
}
open SH,">$dsh/step05.GO1.sh";
open In,$fa;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my $fname=basename($_);
	print SH "diamond9 ",$blast;###
	print SH " --db $dGO --query $_ --evalue 10e-10 --outfmt 5 --threads 8 --out $out/$fname.GO.blast\n";
}
close SH;
close In;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step05.GO1.sh --Resource mem=5G --CPU 8 --maxjob=20";
`$job`;
open SH,">$dsh/step05.GO2.sh";
print SH "perl $Bin/bin/GOanno.pl -i $fa -d $out/ -o $out/GO.anno\n";
close SH;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step05.GO2.sh";
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
Update:			update database file
Description:

Usage:
  Options:
  -fa	<file>	input gene fasta format,
  -type	<file>	input seq type
  -out	<dir>	output data prefix
  -dsh	<dir>	output work sh dir

  -h         Help

USAGE
        print $usage;
        exit;
}
