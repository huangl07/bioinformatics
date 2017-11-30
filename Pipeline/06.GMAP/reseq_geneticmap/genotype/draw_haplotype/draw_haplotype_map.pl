#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
use newPerlBase;
my %config=%{readconf("/home/maxl/script.cfg")};	
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($loc,$fOut,$genotypefile);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$loc,
				"genotype:s"=>\$genotypefile,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($loc and $fOut);

my $head ;
if (defined $genotypefile) {

	open (Gen ,$genotypefile) or die $!;
	while (<Gen>) {
		chomp;
		if (/MarkerID/i || m/BlockID/i){
			$head = $_ ;
			last ;
		}
	}
	close (Gen) ;
}

open (Loc ,$loc) or die $!;
open (Mat ,">$fOut.matrix") or die $!;
while (<Loc>) {
	s/\ba\b/0/g;
	s/\bb\b/1/g;
	s/\bh\b/3/g;
#	s/-/-/g;
	chomp;
	my @temp =split;
	next if (/=/ || /^\s*$/) ;
	if (/MarkerID/i || m/BlockID/i){
		$head = $_ ;
		print Mat $head,"\n" ;
		next ;
	}

	print Mat join("\t",@temp) ,"\n"; 
}

close (Mat) ;
close (Loc) ;

`$config{'matrix2png'} -data $fOut.matrix -rcd -size 10:10 -mincolor 10:70:243 -midcolor  5:108:20 -maxcolor 100:20:37  -missingcolor grey -startrow 2 >$fOut.haplo.png` ;

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
Program:
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage:
  Options:
  -help			USAGE
  -i			loc file, ordered by map, forced 
  -o			output file stem, forced 
  -genotype		genotype file, optional 
USAGE
	print $usage;
	exit;
}

