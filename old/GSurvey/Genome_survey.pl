#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fIn,
	"out:s"=>\$dOut,
			) or &USAGE;
&USAGE unless ($fIn and $dOut);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
my $dsh="$dOut/work_sh";
mkdir $dsh if(!-d $dsh);
open In,$fIn;
open JF,">$dOut/jellyfish.fqlist";
open GCE,">$dOut/gce.fqlist";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($id,$fq1,$fq2)=split(/\s+/,$_);
	if (!-f $fq1 || !-f $fq2) {
		die "$fq1\n$fq2";
	}
	if ($fq1 =~ /gz/) {
		print JF "gunzip -c $fq1\n";
	}else{
		print JF "$fq1\n";
	}
	if ($fq2 =~ /gz/) {
		print JF "gunzip -c $fq2\n";
	}else{
		print JF "$fq2\n"
	}
	print GCE "$fq1\n$fq2\n";
}
close In;
close JF;
close GCE;
open SH,">$dsh/step01.kmer-count.sh";
print SH "jellyfish count -m 17 -s 1000M -t 8 -g $dOut/jellyfish.fqlist -C -o $dOut/jellyfish.kmer.hash && ";
print SH "jellyfish histo -h 255 -o $dOut/jellyfish.freq.stat $dOut/jellyfish.kmer.hash && ";
print SH "Rscript $Bin/bin/genomescope.R  $dOut/jellyfish.freq.stat 17 150 $dOut/jellyfish\n";
print SH "kmer_freq_pread -k 17 -l $dOut/gce.fqlist -c -t 8 -p $dOut/gce 2>$dOut/gce.kmerfreq.log && ";
print SH "gce -f $dOut/gce.freq.stat > $dOut/gce.esimate && ";
print SH "Rscript $Bin/bin/genomescope.R $dOut/gce.freq.stat 17 150 $dOut/gce \n";
print SH "kmer_freq_hash -k 17 -l $dOut/gce.fqlist -c -t 8 -p $dOut/gce-hash 2>$dOut/gce.kmerfreq-hash.log \n"; 
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource --mem-100G --cpu 8 $dsh/step01.kmer-count.sh";
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
  -fqlist	<file>	input file name
  -out		<file>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}
