#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$fOn,$fmb);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i|fasta:s"=>\$fIn,
	"m|model:s"=>\$fOn,
	"o|nex:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn && $fOn && $fOut);
open In,$fIn || die $!;
open Out,">$fOut.nex" || die $!;
my %matrix;
my $name_len=0;
my $len;
$/=">";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($name,@bases)=split /\n/,$_;
	my $bases=join("",@bases);
	$matrix{$name}=$bases;
	if(!defined $len){
		$len=length($bases);
	}else{
		if($len!=length($bases)){
		die "fasta must have same length\n";
		}
	}
	$name_len=length($name) if(length($name)>$name_len);
}
close In;
$/="\n";
my $ntax=scalar keys %matrix;
print Out "#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=$ntax;\n  TAXLABELS\n";
print Out "    $_\n" foreach (sort keys %matrix);
print Out "    ;\nEND;\nBEGIN data;\n  DIMENSIONS ntax=$ntax nchar=$len;\n  FORMAT\n    datatype=DNA\n    missing=?\n    gap=-\n    interleave\n    ;\n  MATRIX\n";
my $N=$name_len."s";
foreach my $key (sort keys %matrix){
	print Out "    ";
	printf Out "%-$N",$key;
	print Out " ","$matrix{$key}\n";
}
print Out "    ;\nEND;\n";
close Out;
my $file=basename("$fOut.nex");
open On,"$fOn" || die "$!\n";
open Out,">$fOut.run.nex" || die "$!\n";
print Out "begin mrbayes;\n\tset autoclose=yes nowarn=yes;\n\texecute $file;\n";
my %model=(
	"GTR"=>6,
	"HKY"=>2,
	"F81"=>1, 
	"mixed"=>"mixed",
);
my %rates=(
	"F81"=>"equal",
	"GTR"=>"gamma",
	0=>"adgamma",
	1=>"propinv",
	"GTR+I"=>"invgamma",
);
while(<On>){
	chomp;
	next if (/^$/|| $_ eq "");
	if(/Model:\s+(\w+)/){
		my $model=$1;
		if ($model=~/GTR+I/) {
			$model="GTR+I";
		}
		print Out "\tlset nst=$model{$model} rates=$rates{$model} ngammacat=4;\n";
	}
	if(/Frequencies:\s+(.*)/){
		my @freq=split(/\s+/,$1);
		print Out "\tprset statefreqpr=fixed(",join(",",@freq).");\n";
		last;
	}
}
print Out "\tmcmcp nchains=4 ngen=90000000 stoprule=yes stopval=0.01 samplefreq=1000 printfreq=1000;\n";
print Out "\tmcmc;\n";
print Out "\tsump;\n";
print Out "\tsumt;\n";
print Out "\tquit;\n";
print Out "end;\n";
close Out;
#######################################################################################
 print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub DIR{
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;
		$dir=`pwd`;
		chomp $dir;
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
sub USAGE {
        my $usage=<<"USAGE";
Contact:        qingmei.cui\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i pop.fasta -o
Usage:
  Options:
  -i	<file>	input fasta 
  -m	<file>	input modeltest output file
  -o	<file>	output nex 
  -b	<file>	mr.sh
  -h         Help

USAGE
        print $usage;
        exit;
}
