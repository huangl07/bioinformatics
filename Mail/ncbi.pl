#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"output:s"=>\$output,
			) or &USAGE;
&USAGE unless ($input and $output);
	open In,"esearch -db pubmed -query \"$input\"\|efetch -format abstract|";
	open Out,">$output";
	$/="\n\n\n";
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		next if (/Article in Chinese/);
		my @line=split(/\n\n/,$_);
		my ($title,$info);
		my @mail;
		for (my $i=0;$i<@line;$i++) {
			if ($i == 0) {
				$info=$line[$i];
				$info=~s/\n//g;
				$info=~s/^\d+\.//g;
			}
			if ($i == 1) {
				$title=$line[$i];
				$title=~s/\n//g;
			}elsif ($line[$i] =~ /Author information/) {
				my @address=split(/\n\(/,$line[$i]);
				foreach my $add (@address) {
					if ($add=~/(\S+\@\S+)/) {
						my $mail=$1;
						$mail=~s/\.$//g;
						$mail=~s/\)//g;
						push @mail,$mail;
					}
				}
			}
		}
		if (!defined $title) {
			next;
		}
		print Out $title,"\t",$info,"\t",join("\,",@mail),"\n";
	}
	close In;
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
