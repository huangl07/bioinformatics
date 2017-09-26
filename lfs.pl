#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($package,$prefix,$step,$profile);
GetOptions(
				"help|?" =>\&USAGE,
				"package:s"=>\$package,
				"prefix:s"=>\$prefix,
				"profile:s"=>\$profile,
				"step:s"=>\$step,
				) or &USAGE;
&USAGE unless ($package and $prefix);
mkdir $package;
mkdir $prefix;
$package=ABSOLUTE_DIR($package);
$prefix=ABSOLUTE_DIR($prefix);
$profile=ABSOLUTE_DIR($profile);
open Out,">$profile";
print Out "export publib=$prefix\n";
print Out "export packge=$package\n";
print Out "export PATH=$prefix/bin:/usr/bin\n";
close Out;
chdir $package;
open LOG,">lfs.log";
$step||=1;
if ($step == 1) {
	open SH,">step01.gcc.sh";
	print SH "wget -c http://mirrors.ustc.edu.cn/gnu/gcc/gcc-7.2.0/gcc-7.2.0.tar.xz\n";
	print SH "tar -xf gcc-7.2.0.tar.xz\n";
	print SH "cd gcc-7.2.0/ && ./contrib/download_prerequisites\n";
	print SH "./configure --prefix=$prefix --disable-multilib\n";
	print SH "make clean && make -j8 && make install\n";
	close SH;
	#`sh step01.gcc.sh`;
}
close LOG;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-package	<dir>	input package dir 
	-prefix 	<dir>	output prefix dir
	-profile	<file>	intput profile file
	-step	<num>	step control
		-1	gcc
		-2	bin
USAGE
	print $usage;
	exit;
}
