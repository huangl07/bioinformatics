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
my ($fIn,$fKey);
GetOptions(
				"help|?" =>\&USAGE,
				"l:s"=>\$fIn,
				"k:s"=>\$fKey,
				) or &USAGE;
&USAGE unless ($fIn && $fKey);
#######################################################################################

# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------

my $header_done_flag=0;
my %Header=();

open (IN,"<",$fIn) or die $!;
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
	next if ($line =~ /^$/ || $line =~ /^;/) ;
	last if ($line =~ /locus numbers/ || $line =~ /individual names/) ;
	$line=~s/^\s+//;
	my ($real,$anno)=split(";",$line);
	$Marker_str.=" ".$real;
}
close (IN) ;
$Marker_str=~s/^\s+//;


my @marker=split /\s+/,$Marker_str;
#print "@marker","\n";exit;

my @SexAver_arr=();
my @Male_arr=();
my @Female_arr=();


my @marker_arr=();
for (my $i=0;$i<$Header{"nloc"} ;$i++) {
	
	my @this_marker=();
	my $cross_type ="";

	push @this_marker,shift @marker;

	while ($marker[0] =~ /\<\S+\>/ || $marker[0] =~ /\(\S+\)/ || $marker[0] =~ /\{\S+\}/) {

		if ($marker[0] =~ /\<\S+\>/) {
			($cross_type) = $marker[0] =~ /\<(\S+)\>/;
		}

		push @this_marker, shift @marker;
	}
	for (my $j=0;$j<$Header{"nind"} ;$j++) {
		push @this_marker, shift @marker;
	}

	push @SexAver_arr, join("\t",@this_marker);

	if ($cross_type eq "lmxll") {
		push @Male_arr,join("\t",@this_marker);

	}elsif($cross_type eq "nnxnp"){
		push @Female_arr,join("\t",@this_marker);

	}elsif($cross_type eq "abxcd"){
		my @male_arr = ();
		my @female_arr = ();

		foreach my $elem (@this_marker) {
			my $male_elem = $elem;
			my $female_elem = $elem;

			$male_elem=~s/\babxcd\b/lmxll/g;
			$male_elem=~s/\bac\b/ll/g;
			$male_elem=~s/\bbc\b/lm/g;
			$male_elem=~s/\bad\b/ll/g;
			$male_elem=~s/\bbd\b/lm/g;

			$female_elem=~s/\babxcd\b/nnxnp/g;
			$female_elem=~s/\bac\b/nn/g;
			$female_elem=~s/\bbc\b/nn/g;
			$female_elem=~s/\bad\b/np/g;
			$female_elem=~s/\bbd\b/np/g;

			push @male_arr, $male_elem;
			push @female_arr, $female_elem;

		}

		push @Male_arr,join("\t",@male_arr);
		push @Female_arr,join("\t",@female_arr);
			
	}elsif($cross_type eq "efxeg"){
		my @male_arr = ();
		my @female_arr = ();

		foreach my $elem (@this_marker) {
			my $male_elem = $elem;
			my $female_elem = $elem;

			$male_elem=~s/\befxeg\b/lmxll/g;
			$male_elem=~s/\bee\b/ll/g;
			$male_elem=~s/\bef\b/lm/g;
			$male_elem=~s/\beg\b/ll/g;
			$male_elem=~s/\bfg\b/lm/g;

			$female_elem=~s/\befxeg\b/nnxnp/g;
			$female_elem=~s/\bee\b/nn/g;
			$female_elem=~s/\bef\b/nn/g;
			$female_elem=~s/\beg\b/np/g;
			$female_elem=~s/\bfg\b/np/g;

			push @male_arr, $male_elem;
			push @female_arr, $female_elem;

		}

		push @Male_arr,join("\t",@male_arr);
		push @Female_arr,join("\t",@female_arr);
			
	}

}


# male and female marker's linkage phase need to be modified, such as lmxll from {00} to {0-}
for (my $i=0;$i<@Male_arr ; $i++) {

	if ($Male_arr[$i] =~ /{..}/) {
		$Male_arr[$i] =~ s/.}/\-}/;
	}
}

for (my $i=0;$i<@Female_arr ; $i++) {

	if ($Female_arr[$i] =~ /{..}/) {
		$Female_arr[$i] =~ s/{./{\-/;
	}
}


open (SEXAVER,">","$fKey.sexAver.loc") or die $!;
my $header_str="";
$header_str.="name = "."SA.".$Header{"name"}."\n";
$header_str.="popt = ".$Header{"popt"}."\n";
$header_str.="nloc = ".scalar @SexAver_arr."\n";
$header_str.="nind = ".$Header{"nind"}."\n";
print SEXAVER $header_str,"\n",join("\n",@SexAver_arr),"\n";
close (SEXAVER) ;


open (MALE,">","$fKey.male.loc") or die $!;
$header_str="";
$header_str.="name = "."M.".$Header{"name"}."\n";
$header_str.="popt = ".$Header{"popt"}."\n";
$header_str.="nloc = ".scalar @Male_arr."\n";
$header_str.="nind = ".$Header{"nind"}."\n";
print MALE $header_str,"\n",join("\n",@Male_arr),"\n";
close (MALE) ;



open (FEMALE,">","$fKey.female.loc") or die $!;
$header_str="";
$header_str.="name = "."F.".$Header{"name"}."\n";
$header_str.="popt = ".$Header{"popt"}."\n";
$header_str.="nloc = ".scalar @Female_arr."\n";
$header_str.="nind = ".$Header{"nind"}."\n";
print FEMALE $header_str,"\n",join("\n",@Female_arr),"\n";
close (FEMALE) ;







#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:generateSexSpeficLocFiles_v1.pl
Version: $version
Contact: Jiang Chuanbei <jiangchb\@biomarker.com.cn> <jiangxiao8\@qq.com>
Description:

	根据输入的loc文件，通过拟测交的方式生成雌性，雄性，中性loc文件。
	给定的loc文件可一个带有连锁相，也可以没有连锁相信息。
	
Usage:
  Options:
	-l	<input file>   	input loc files, only loc files are required,
	-k	<key>	    	key of output files, required,
	-h					Help

USAGE
	print $usage;
	exit;
}
