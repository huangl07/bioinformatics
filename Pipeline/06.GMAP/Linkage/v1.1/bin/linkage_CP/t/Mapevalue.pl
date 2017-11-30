#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use lib "/share/nas1/huangl/perlModuals/Statistics-LineFit/lib/Statistics/";
use LineFit;
my $BEGIN_TIME=time();
my $version="1.0.0";
$Script=~s/\.pl//g;
my @Times=localtime();
my $year=$Times[5]+1990;
my $month=$Times[4]+1;
my $day=$Times[3];
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
my @map=glob("$fIn/*.map");
open Out,">>$fOut" or die $!;
foreach my $map (@map) {
	my $name=basename($map);
	$name=~s/\.map//g;
	my $loc=$map;
	$loc=~s/\.map/.loc/g;
	open In,$map;
	my %marker;
	my $max=0;
	while (<In>) {
		s/\r//g;
		chomp;
		next if ($_ eq "" || /^$/ || /^;/ || /group/  || /Position/ || /markers/ || /Printing/);
		my @temp=split(/\s+/,$_);
		if (scalar @temp == 0) {
			next;
		}
		if ($temp[0] eq "") {
			shift @temp;
		}
		$marker{$temp[0]}=$temp[1];
		$max=$temp[1] if ($max < $temp[1]);
	}
	my @y;
	my @x;

	foreach my $k (sort keys %marker) {
		push @y,$marker{$k};
		my $t=$k;
		$t=~s/Marker//g;
		push @x,$t;
	}
	my $lineFit = Statistics::LineFit->new();
	 $lineFit->setData (\@x, \@y);
	my ($intercept, $slope) = $lineFit->coefficients();
	my $rSquared = $lineFit->rSquared();
	print Out $map,"\t",$max,"\t",$rSquared,"\n";
	close In;
#    mkdir "tmp" if (!-d "tmp");				
#	`perl /share/nas1/shenl/bin/Joinmap/statHaploSourceNoise.pl -i $loc -M $map -o $name -od tmp`;
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath{
        my ($type,$input) = @_;

        my $return;

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chop $return;
                $return .="/".$file;
                chdir($pwd);
        }
        return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[$month:$day:$year:]
	Contact:Huang Long <huangl\@biomarker.com.cn>
	Options:
		-i	<file>	input file
		-o	<file>	output file
		-h	Help

USAGE
	print $usage;
	exit;
}
