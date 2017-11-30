#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use lib "$Bin/Statistics-RankCorrelation-0.1205/lib/Statistics/";
use RankCorrelation;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

#modified by liangsh <2015.9.24>
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$locfile,$posifile,$msto,$fKey);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"l:s"=>\$locfile,
				"p:s"=>\$posifile,
				"m:s"=>\$msto,
				"k:s"=>\$fKey,
				"o:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($fIn and $fKey and $locfile and $dOut && $posifile && $msto);

mkdir ($dOut) if (!-d $dOut);
$fIn = Cwd::abs_path($fIn);
$locfile = Cwd::abs_path($locfile);
$posifile = Cwd::abs_path($posifile);
$msto = Cwd::abs_path($msto);
$dOut = Cwd::abs_path($dOut);

##########################################read posifile <liangsh>
## read posi
#
open (POS,"<$posifile") or die $!;
my %posi;
while (<POS>) {
	chomp;
	my @a = split /\s+/,$_;
	$posi{$a[0]} = $a[2];
}
close POS;
##########################################
#
# read loc 
#
open (LOC ,$locfile) or die $!;
my $head ;
my %marker_info;
while (<LOC>) {
	chomp;
	if (/=/){
		$head .= "$_\n" ;
	}else{
next if (/^$/);
		s/\bA\b/a/g;
		s/\bB\b/b/g;
		s/\bU\b/-/g;
		s/\bX\b/h/g;
		my ($marker,@genotype) = split ;
		next if (@genotype < 5);
		$marker_info{$marker} = \@genotype;
	}
}
close (LOC) ;

############################judge the map order by spearman cofficient	<liangsh>
open (IN,"<$fIn") or die $!;
my (@marker,@md,@group);
while (<IN>) {
	chomp;
	next if (/^\s*$/ || /^;/ || /group/) ;
	my ($marker,$gd) = split;
	push @marker,$marker;
	push @md,$gd;
	push @group,$marker,$gd;
}
close IN;
my @ref;
for (my $i=0; $i<@marker; $i++) {
	if (! $posi{$marker[$i]}){die "The posifile is wrong!\n";};
	push @ref,$posi{$marker[$i]};
}
my $c = Statistics::RankCorrelation->new(\@md, \@ref);
my $n = $c->spearman;
#print $n,"\n";
####################



#
# output map 
#

open (NEWLOC,">$dOut/$fKey.order.loc") or die $!;
open (MAP,">$dOut/$fKey.map") or die $!;
if ($n>0) {
	print MAP "group\t0\n";
	print NEWLOC "$head\n";
	for (my $i=0; $i<@group;$i+=2) {
		print MAP "$group[$i]\t$group[$i+1]\n";
		print NEWLOC $group[$i],"\t" , join("\t",@{$marker_info{$group[$i]}}) , "\n" ;
	}
}elsif ($n<0) {
	print MAP "group\t0\n";
	print NEWLOC "$head\n";
	printf MAP "%s\t%.3f\n",$group[-2],0.000;
	print NEWLOC $group[-2],"\t" , join("\t",@{$marker_info{$group[-2]}}) , "\n" ;
	my $sum;
	for (my $i=2; $i<@group; $i+=2) {
			$sum += ($group[-$i+1]-$group[-$i-1]);
			$sum =~ s/^\s+//;
			printf MAP "%s\t%.3f\n",$group[-$i-2],$sum;
			print NEWLOC $group[-$i-2],"\t" , join("\t",@{$marker_info{$group[-$i-2]}}) , "\n" ;
		}

}

close NEWLOC;
close MAP;


#############
#formate the order of mst.o file
##############
open (IN,"<$msto") or die $!;
open (OUT,">$dOut/$fKey.mst_order.o") or die $!;
my $flag = 0;
#my %
if ($n>0) {
	while (<IN>) {
		chomp;
		if (/==========/) {
			print OUT $_,"\n";
			$flag++;
			next;
		}
		if ($flag == 2) {
			print OUT $_,"\n";
		}else{next;}
	}
}

my @mar;
my %rf;
if ($n<0) {
	while (<IN>) {
		chomp;
		if (/==========/) {
			print OUT $_,"\n";
			$flag++;
			next;
		}
		next unless ($flag == 2);
		if (/^\s+/) {
			s/^\s+//;
			@mar = split /\s+/,$_;
		}
		my ($marker,@recom)=split;
		for (my $i=0; $i<@recom; $i++) {
			$rf{$marker}{$mar[$i]} = $recom[$i];
		}
	}
##
	my @ordermap;
	for (my $i=1; $i<=@marker; $i++) {
		push @ordermap,$marker[-$i];
	}
	my $string = join ("\t",@ordermap);
	print OUT "\t$string\n";
	foreach my $m (@ordermap) {
		#die "$m\n";
		print OUT $m,"\t";
		for (my $i=0; $i<@ordermap; $i++) {
			print OUT $rf{$m}{$ordermap[$i]},"\t";
		}
		print OUT "\n";
	}
	
}
close IN;
close OUT;









=pop

	open (MSTR ,$fIn) or die $! ;
	open (NEWLOC,">$dOut/$fKey.order.loc") or die $!;
	open (MAP,">$dOut/$fKey.map") or die $!;
	
	print NEWLOC " $head \n" ;
	while (<MSTR>) {
		chomp;
		next if (/^\s*$/ || /^;/) ;
		my ($marker, $length) = split ;
		if ($marker ne "group") {
			print MAP "$_\n" ;
	}else {
			print MAP "group\t0\n";
			next;
		}
		print NEWLOC $marker,"\t" , join("\t",@{$marker_info{$marker}}) , "\n" ;
	}
	close (MSTR) ;
	close (MAP) ;
	close (NEWLOC) ;



=cut 



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
Program: $0
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage: 将MSTmap结果转化为画图所需的格式并提取对应顺序的loc文件
  Options:

  -help		USAGE
  -i		MSTmap map file, forced 
  -l		input loc file, forced
  -p		input posi file forced
  -m		input mst.o file forced
  -k		output file stem, forced 
  -o		output directory 

USAGE
	print $usage;
	exit;
}

