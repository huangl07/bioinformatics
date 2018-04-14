#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$type);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use utf8;
my $version="1.0.0";
GetOptions(
        "help|?" =>\&USAGE,
        "input:s"=>\$fin,
        "output:s"=>\$fout,
		"type:s"=>\$type,
                        ) or &USAGE;
&USAGE unless ($fin );
open In,$fin;
open Out,">$fout.$type.report.xls";
print Out "ID\tchr number\tTotal number\tType\tRef\tAlt\tMarker size(bp)\t";
print Out "FORWARD PRIMER1 (5'-3')\tTemperature\tsize\tREVERSE PRIMER1 (5'-3')\tTemperature\tsize\tPRODUCT1 size (bp)\tstart (bp)\tend (bp)\t";
print Out "FORWARD PRIMER2 (5'-3')\tTemperature\tsize\tREVERSE PRIMER2 (5'-3')\tTemperature\tsize\tPRODUCT2 size (bp)\tstart (bp)\tend (bp)\t";
print Out "FORWARD PRIMER3 (5'-3')\tTemperature\tsize\tREVERSE PRIMER3 (5'-3')\tTemperature\tsize\tPRODUCT3 size (bp)\tstart (bp)\tend (bp)\n";
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/|| /#/);
	my ($id,$chr_nr,$num,$type,$ref,$alt,$size,$start,$end,$nstart,$nend,$lpo,$lto,$lso,$rpo,$rto,$rso,$os,$ostart,$oend,$lpt,$ltt,$lst,$rpt,$rtt,$rst,$ts,$tstart,$tend,$lpr,$ltr,$lsr,$rpr,$rtr,$rsr,$rs,$rstart,$rend)=split(/\s+/,$_);
	my $astart=$nstart + $ostart ;
	my $aend=$nstart + $oend ;
	my $bstart=$nstart + $tstart ;
	my $bend=$nstart + $tend ;
	my $cstart=$nstart + $rstart ;
	my $cend=$nstart + $rend ;
	print Out join("\t",$id,$chr_nr,$num,$type,$ref,$alt,$size,$lpo,$lto,$lso,$rpo,$rto,$rso,$os,$astart,$aend,$lpt,$ltt,$lst,$rpt,$rtt,$rst,$ts,$bstart,$bend,$lpr,$ltr,$lsr,$rpr,$rtr,$rsr,$rs,$cstart,$cend),"\n";
}
close In;
close Out;
#########################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#########################################################################
sub  USAGE{#
	my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:                 $Script
Description:
        fq thanslate to fa format
        eg:
        perl $Script -i -o -k -c

Usage:
  Options:
  -input	<file>  input file name
  -output	<file>  output base name
  -type	<str>	marker type
  -h         Help

USAGE
        print $usage;
        exit;
}

