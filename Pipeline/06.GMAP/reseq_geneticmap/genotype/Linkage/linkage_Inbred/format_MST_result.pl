#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$locfile,$fKey,$pwdfile);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"l:s"=>\$locfile,
				"k:s"=>\$fKey,
				"o:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($fIn and $fKey and $locfile and $dOut );

mkdir ($dOut) if (!-d $dOut);
$fIn = Cwd::abs_path($fIn);
$locfile = Cwd::abs_path($locfile);
$dOut = Cwd::abs_path($dOut);

#
# read loc 
#
open (LOC ,$locfile) or die $!;
my $head ;
my %marker_info;
while (<LOC>) {
	chomp;
	if (/=/){
		$head .= "$_\n" ;
	}else{
		s/\bA\b/a/g;
		s/\bB\b/b/g;
		s/\b\bU/-/g;
		s/\bX\b/h/g;
		my ($marker,@genotype) = split ;
		next if (@genotype < 5);
		$marker_info{$marker} = \@genotype;
	}
}
close (LOC) ;

#
# output map 
#
open (MSTR ,$fIn) or die $! ;
open (NEWLOC,">$dOut/$fKey.order.loc") or die $!;
open (MAP,">$dOut/$fKey.map") or die $!;

print NEWLOC "$head\n" ;
while (<MSTR>) {
	chomp;
	next if (/^\s*$/ || /^;/) ;
	my ($marker, $length) = split ;
	if ($marker ne "group") {
		print MAP "$_\n" ;
	}else {
		print MAP "group\t0\n";
		next;
	}
	print NEWLOC $marker,"\t" , join("\t",@{$marker_info{$marker}}) , "\n" ;
}
close (MSTR) ;
close (MAP) ;
close (NEWLOC) ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage: 将MSTmap结果转化为画图所需的格式并提取对应顺序的loc文件
  Options:

  -help		USAGE
  -i		MSTmap map file, forced 
  -l		input loc file, forced
  -k		output file stem, forced 
  -o		output directory 

USAGE
	print $usage;
	exit;
}

