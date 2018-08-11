#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($input,$output);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"output:s"=>\$output,
			) or &USAGE;
&USAGE unless ($input and $output);
open In,$input;
my %genewise_outs;
my $target_seq;
$/ = "genewise output";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if ( /Score\s+(\d+\.\d+)\s+bits/ ) {
        $genewise_outs{$_} = $1;
    }
    elsif ( /Target\s+Sequence\s+(.*?)\n/msg ) {
        $target_seq = $1;
    }

}
close In;
my $best=(sort {$genewise_outs{$b}<=>$genewise_outs{$a}} keys %genewise_outs)[0];
my @blocks = split( "//", $best );
my $frameshift = "NO";
my $stopsite   = "NO";
my $translation;
if ( $blocks[0] =~ /\!/ ) {
    $frameshift = "YES";
}

if ( $blocks[0] =~ /\*/ ) {
    $stopsite = "YES";
}
open Out,">$output";
print Out "$frameshift\t$stopsite\n";
close Out;
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
	
Usage:
  Options:
  -input	<file>	input dir
  -output	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
