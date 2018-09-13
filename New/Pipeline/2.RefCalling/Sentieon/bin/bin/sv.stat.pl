#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
			) or &USAGE;
&USAGE unless ($fIn );
$dOut||="./";
open In,$fIn;
$/="#";
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my($format,@info)=split(/\n/,$_);
	my @b; #=grep(".+",@info);
	foreach(@info){
		next if($_ eq "");
		my $a=$_;
		$a=~ s/\,/\t/g;
		push @b,$a;
	}
	if($format =~ "Variantss by type"){
		open OUT,">$dOut/sv.stat";
		print OUT join("\n",@b),"\n";
		close OUT;
	}
	if($format =~ "Count by effects"){
		open OUT,">$dOut/sv.effects";
		print OUT join("\n",@b),"\n";
		close OUT;
	}
	if($format =~ "Count by genomic region"){
		open OUT,">$dOut/sv.region";
		print OUT join("\n",@b),"\n";
		close OUT;
	}
}
close In;
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}
