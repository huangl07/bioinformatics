#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$fLG,$fKey,$type);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"l:s"=>\$fLG,
				
				"d:s"=>\$dOut,
				"k:s"=>\$fKey,
				
				"t:s"=>\$type,
				
				) or &USAGE;
&USAGE unless ($fIn and $dOut and $fLG and $fKey);

mkdir ($dOut) if (!-d $dOut) ;
$fIn  = Cwd::abs_path($fIn) ;
$dOut = Cwd::abs_path($dOut) ;
$fLG  = Cwd::abs_path($fLG);


#
# read lg file 
#
my %lg_marker;
$/="\>";

open ( LG , $fLG ) or die $!;
<LG>;
while (<LG>) {
	s/\>//g;
	next if (/^$/) ;

	my ($head ,@line) = split /\n/ , $_ ;
	my $lg = ( split /\s+/ , $head ) [0] ;
	
	my @markers = ();
	foreach my $line (@line) {
		my @tmp = split/\s+/,$line;
		push @markers, @tmp;
	}
	$lg_marker{$lg}= \@markers ;
}
close (LG) ;


#
# read genotype file 
#
my $head ;
my %marker_genotype;

open (GEN , "$fIn") or die $!;
$/ = "\n" ;
while (<GEN>) {
	chomp;
	next if (/^\s*$/ ) ;

	$head = $_ if ( /Type/ || /type/) ;
	my ($marker,@genotype) = split /\s+/ , $_ ;
	$marker_genotype{$marker} = join("\t",@genotype) ;
}
close (GEN) ;

#
# output 
#

foreach my $lg (sort {($a=~/(\d+)/)[0] <=> ($b=~/(\d+)/)[0]}  keys %lg_marker ) {
	
	open (LG,">","$dOut/$fKey.$lg.genotype") or die $!;
	
	print LG "$head\n";
	foreach my $marker (@{$lg_marker{$lg}}) {
		next if (!exists $marker_genotype{$marker});
		print LG $marker , "\t" ,$marker_genotype{$marker},"\n" ;
	}
	close (LG) ;

	if ( (-f "$Bin/genotype2loc.pl" ) && defined $type) {
		my $command = "perl $Bin/genotype2loc.pl -i $dOut/$fKey.$lg.genotype -o $dOut/$fKey.$lg.loc -t $type" ;
		`$command ` ;
	}
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
Program: $0
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage:  将genotype文件按连锁群分割
  Options:
  -help			USAGE,
  -i			genotype file， forced
  -l			linkage group  file, forced
  -k			output file stem, forced 
  -d			output directory
  -t			population type  
   
USAGE
	print $usage;
	exit;
}

