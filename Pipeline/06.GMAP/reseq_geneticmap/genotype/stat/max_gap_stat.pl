#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
open (OUT,">$fOut") or die$!;
my @files =glob ("$fIn/*.map");
foreach $fIn (@files) {
	my ($num,$tra)=$fIn=~/(chr\d+).(\w+)\.map/i;
	open (IN,"$fIn") or die $!;
	my %words;
	my $max=0;
	my $i=1;
	my $j=1;
	while (<IN>) {
		chomp;
		next if (/^$/);
		if (/Marker|Block/i) {
		my @mar=split (/\t/,$_);
		$words{$i}=$mar[1];
		$i++;
		}
	}
	my @num=(sort {$a <=> $b} values %words);
	my $size=@num;
	while ($j<$size) {
		my $cha=sprintf "%.3f",$num[$j]-$num[$j-1];
		#my $cha=$num[$j]-$num[$j-1];
		#print "$num[$j]\t$num[$j-1]\t$cha\n";
		#print OUT "$cha \n";
		if ($cha>$max) {
			$max=$cha;
		}
		$j++;
	}

print OUT "$num\t$max\n";
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
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
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: /Wheat_Polymorphismmap_BMK131206-748/Analysis
Version: $version
Contact:Yuan ZhengWen <yuanzw\@biomarker.com.cn> 
Description:
Usage:
  Options:
  -i mapÎÄ¼þ
  -o     
 
  -h         Help

USAGE
	print $usage;
	exit;
}
