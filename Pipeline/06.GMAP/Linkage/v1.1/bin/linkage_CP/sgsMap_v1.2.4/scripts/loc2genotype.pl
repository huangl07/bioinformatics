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
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				) or &USAGE;
&USAGE unless ($fIn);

my @file=();
if (-f $fIn) {
	push @file,$fIn;
}else{
	push @file,glob("$fIn/*.loc");
}
if (! @file) {
	warn "no file $fIn or on files in $fIn direction.\n";
	exit (-1);
}


foreach my $file (@file) {
	
	my $dir = dirname($file);
	$fOut=basename($file);

	if ($fOut =~ /(\S+)\.loc$/) {
		$fOut = $1.".genotype";
	}else{
		$fOut.=".genotype";
	}
	
	$fOut = $dir."/".$fOut;
	print $fOut,"\n";

	my $header_done_flag=0;
	my %Header=();

	open (IN,"<",$file) or die $!;
	$/="\n";

	while (! $header_done_flag) {
		my $line=<IN>;
		chomp $line;
		$line=~s/\r//g;
		next if ($line =~ /^$/ || $line =~ /^;/ || $line =~ /^\s+;/) ;
		last if ($line =~ /locus numbers/) ;

		if ($line =~ /name\s+\=\s+(\S+)/) {
			$Header{"name"}=$1;
		}elsif($line =~ /popt\s+\=\s+(\S+)/){
			$Header{"popt"}=$1;
		}elsif($line =~ /nloc\s+\=\s+(\S+)/){
			$Header{"nloc"}=$1;
		}elsif($line =~ /nind\s+\=\s+(\S+)/){
			$Header{"nind"}=$1;
		}

		if (exists $Header{"name"} && exists $Header{"popt"} && exists $Header{"nloc"} && exists $Header{"nind"} ) {
			$header_done_flag=1;
		}
	}


	my $Marker_str="";
	$/="\n";
	while (my $line=<IN>) {
		chomp $line;
		$line=~s/\r//g;
		next if ($line =~ /^$/ || $line =~ /^;/ ) ;
		last if ($line =~ /locus numbers/ || $line =~ /individual names/) ;
		my ($real,$anno)=split(";",$line);
		$Marker_str.=" ".$real;
	}
	close (IN) ;


	open (OUT,">",$fOut) or die $!;

	print OUT "MarkerID\ttype";
	for (my $i=0;$i<$Header{"nind"} ;$i++) {
		print OUT "\t","aa";
	}
	print OUT "\n";

	$Marker_str=~ s/^\s+//;
	$Marker_str=~ s/\{\S+\}//g;
	$Marker_str=~ s/\(\S+\)//g;
	my @marker=split /\s+/,$Marker_str;

	my @marker_arr=();
	for (my $i=0;$i<@marker ;$i++) {
		if ($i>0 && $i % ($Header{"nind"}+2) == 0) {
			print OUT join("\t",@marker_arr),"\n";
			@marker_arr=();
		}
		push @marker_arr,$marker[$i];
	}

	if (@marker_arr) {
		print OUT join("\t",@marker_arr),"\n";
		@marker_arr=();
	}

	close (OUT) ;

}



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
Program: loc2genotype.pl (把JoinMap生成的loc文件转换为genotype文件)
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>
Description:
	把JoinMap的loc文件转换为genotype文件，loc文件可能是换行的，本程序也支持。
	loc文件中的nloc，nind的参数是正确的。否则转换会出错
	程序伪造了个体的编号，都是aa
	程序测试完全正确，放心使用。

Usage:
  Options:
  -i <file|dir>  loc file, loc format, forced
  -h         Help

USAGE
	print $usage;
	exit;
}
