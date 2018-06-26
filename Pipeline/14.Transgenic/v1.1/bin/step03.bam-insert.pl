#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamlist,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bamlist:s"=>\$bamlist,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($bamlist and $out and $dsh);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$bamlist=ABSOLUTE_DIR($bamlist);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step03.bam-insert.sh";

open IN,$bamlist;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($id,$bam)=split/\t/,$_;
    print SH "perl $Bin/bin/retrive.insert.bam.v2.0.pl -id $id -bam $bam -out $out/ && perl $Bin/bin/retrive.insert.id.pl -bam $bam -id $id -bed $out/$id.bed -out $out/\n";
    }
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=50G $dsh/step03.bam-insert.sh";
`$job`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script 

Usage:
  Options:
  -bamlist  <file> bam list file  
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
