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
my ($genotype,$type,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$genotype,
				"t:s"=>\$type,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($genotype and $type and $fOut);


my $nind ;
my %mark_info;

#
# read in genotype
#
open ( Gen , $genotype ) or die $! ;
while (<Gen>) {
	next if (/^MarkerID/ || /^\s+$/) ;
	last if (/^locus/) ;

	s/\baa\b/a/g;
	s/\bbb\b/b/g;
	s/\bab\b/h/g;
	s/\-\-/\-/g;

	my ($marker,$type,@sam_type) = split ;
	$nind = @sam_type ;
	next if (@sam_type < 5 ) ;
	$mark_info{$marker}  = \@sam_type ;
}
close (Gen) ;

#
# output loc file 
#

open (Loc ,">",$fOut) or die $!;
my @marker = keys %mark_info ;
my $nloc = @marker;
my $name = basename ($fOut) ;
my @name = split/\./,$name;
$name=$name[1]; 
#$name =~ s/\.loc//ig ;
#$name=~/\.(?<chr>\w+)(?<num>\d+)\./;
#print "$+{chr}$+{num}";
print Loc "name = $name\n" ;
print Loc "popt = $type\n" ;
print Loc "nloc = $nloc\n" ;
print Loc "nind = $nind\n\n" ;

foreach my $marker (sort {($a=~/(\d+)/)[0] <=> ($b=~/(\d+)/)[0]} @marker ) {
	print Loc $marker,"\t", join("\t",@{$mark_info{$marker}}) ,"\n";
}
close (Loc) ;

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

Usage: 将近交群体genotype文件转换成JoinMap格式loc文件 

  Options:
  -i		genotype file, forced
  -t		population type, forced   
  -o		output file
  -h		help

USAGE
	print $usage;
	exit;
}

