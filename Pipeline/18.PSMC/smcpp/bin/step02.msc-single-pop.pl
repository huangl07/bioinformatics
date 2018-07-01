#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$group,$chr,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "vcf:s"=>\$vcf,
	"group:s"=>\$group,
    "chr:s"=>\$chr,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($group);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step02.smc-single-pop.sh";
open IN,$group;
my %hash;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($sample,$pop)=split/\t/,$_;
    push @{$hash{$pop}},$sample;
    }
close IN;
foreach my $pop (keys %hash){
    #my $pops=join("\,",@{$hash{$pop}});
    for(my $i=1;$i<=$chr;$i++){
        my $chr="chr$i";
        my $outfile="$pop.$chr.smc.gz";
        print SH "smc++ vcf2smc $vcf $out/$outfile $chr $pop:",join("\,",@{$hash{$pop}}),"\n";
    }
}
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step02.smc-single-pop.sh";
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
	run smc
	eg:
	perl $Script -group -chr -out -dsh -h

Usage:
  Options:
  -vcf  <file>  the filter vcf file
  -group    <file>	input group list files#tab format eg:sample1    group;
  -chr  <number> input the number of chr
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
