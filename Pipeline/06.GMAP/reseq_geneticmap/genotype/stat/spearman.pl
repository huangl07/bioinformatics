#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use lib "/share/nas1/huangl/perlModuals/Statistics-RankCorrelation-0.1203/lib/Statistics/";
use RankCorrelation;
use File::Basename qw(basename dirname);
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
my ($fIn,$fOut,$fPosi);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"p:s"=>\$fPosi,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
my @map=glob("$fIn/*.map");
my %out;
open In,$fPosi;
my %Pos;
my %filter;
while (<In>) {
	chomp;
	s/\,/\t/g;
	next if ($_ eq "" ||/^$/);
	my ($id,$chr,$pos,undef)=split(/\s+/,$_);
	$Pos{$chr}{$id}=$pos;
	$filter{$id}=$chr;
}
close In;
open Out,">$fOut";
foreach my $map (@map) {
	open In,$map;
	my @marker1;
	my @marker2;
	my $max=0;
	my $chr=(split(/\./,basename($map)))[0];
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /group/ || /^;/) ;
		my ($ID,$pos)=split(/\s+/,$_);
		next if (!exists $filter{$ID});
		push @marker1,$pos;
		push @marker2,$Pos{$filter{$ID}}{$ID};
	}
	close In;
	my $c = Statistics::RankCorrelation->new( \@marker2, \@marker1,   );
	my $n=$c->spearman;
	print Out basename($map),"\t",$max,"\t",$n,"\n";
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
		-p	<file>	posi file
		-h	Help

USAGE
	print $usage;
	exit;
}
