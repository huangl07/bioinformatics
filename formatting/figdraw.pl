#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fMap,$fPos,$fLoc,$fChr,$dOut,$Key);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fMap,
				"p:s"=>\$fPos,
				"l:s"=>\$fLoc,
				"c:s"=>\$fChr,
				"d:s"=>\$dOut,
				"k:s"=>\$Key,
				) or &USAGE;
&USAGE unless ($fMap and $dOut and $Key and $fPos and $fLoc and $fChr);
mkdir $dOut if (!-d $dOut);
mkdir "$dOut/tmp" if (!-d "$dOut/tmp");
my($groupnumber);
open (MAPIN,"$fMap") || die "Can't open $fMap\n" ;
open (POSIN,"$fPos") || die "Can't open $fPos\n" ;
open (LOCIN,"$fLoc") || die "Can't open $fLoc\n" ;
open MAPTMP,">$dOut/tmp/$Key.tmp.map" || die $!;
open POSTMP,">$dOut/tmp/$Key.tmp.pos" || die $!;
open LOCTMP,">$dOut/tmp/$Key.tmp.loc" || die $!;

##format .map file for map-plot.R and colinearity-plot.R
print MAPTMP "id\tcm\tlg\n";
while (<MAPIN>) {
	chomp;
	if (/^group (.*)/){
		$groupnumber=(split(/\s+/,$_))[-1];
	}else{
		print MAPTMP "$_\t$groupnumber\n";
	}
}

##format .filt.pos file for colinearity-plot.R
print POSTMP "id\tchr\tpos\n";
while (<POSIN>) {
	next if /^#MarkerID/;
	print POSTMP $_;
}
close MAPIN;
close POSIN;
close MAPTMP;
close POSTMP;

##combine $Key.tmp.map and $Key.R.final.loc file for lod-genetype.R
open (MAPTMPIN,"$dOut/tmp/$Key.tmp.map") || die "Can't open $dOut/tmp/$Key.tmp.map\n" ;
my %table;
while (<MAPTMPIN>){
	next if /^id/;
	chomp;
	my($id,$cm,$lg) = split /\t/;
	$table{$id}{'id'} = $id;
	$table{$id}{'cm'} = $cm;
	$table{$id}{'lg'} = $lg;
}
close MAPTMPIN;
while (<LOCIN>){
	chomp;
	if (/^#MarkerID/){
		my ($junk,@indi) = split /\t/;
		my $combindi = join(",",@indi);
		print LOCTMP "id,,,$combindi\n";
	}
	else{
		my ($id,@geno) = split /\t/;
		my $combgeno = join(",",@geno);
		$table{$id}{'geno'} = uc $combgeno;
	}
}
for (sort idorder keys %table){
	print LOCTMP "$table{$_}{'id'},$table{$_}{'lg'},$table{$_}{'cm'},$table{$_}{'geno'}\n";
}
close LOCIN;
close LOCTMP;

`Rscript $Bin/map-plot.R -m $dOut/tmp/$Key.tmp.map -k $Key -d $dOut`; ##draw fig 3-9
`Rscript $Bin/plotrh.R  -m $dOut/tmp/$Key.mappos -c $fChr -d $dOut -k $Key`;##draw fig 3-11
`Rscript $Bin/colinearity-plot.R -m $dOut/tmp/$Key.tmp.map -p $dOut/tmp/$Key.tmp.pos -d $dOut -k $Key`;##draw fig 3-14

my $job ="Rscript $Bin/plotrf.R -m $dOut/tmp/$Key.tmp.loc -d $dOut -k $Key &&\n";##draw fig 3-10 3-13
$job .= "Rscript $Bin/plotgeno.R -m $dOut/tmp/$Key.tmp.loc -d $dOut -k $Key &&\n";##draw fig 3-12
open SH, ">$dOut/lodg.sh" || die $!;
print SH $job;
close SH;
my $jobfil = "$dOut/lodg.sh";
print "$Bin/../qsub-sge.pl -resource mem=100G -ppn 16 $jobfil\n";
`$Bin/../qsub-sge.pl -resource mem=100G -ppn 16 $jobfil`;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub idorder{
	(my $da = $a) =~ s/\D+//g;
	(my $db = $b) =~ s/\D+//g;
	$da <=> $db;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: figdraw.pl for fig 3-9 3-10 3-12 3-13 3-14
Version: $version
Contact: Ao.li

Usage:
  Options:
	-i <file>		input final.map file, forced
	-p <file>		input filt.pos file, forced
	-l <file>		input R.final.loc file, forced
	-c <file>		input filt.chr file, forced
	-d <file>		output dir,forced
	-k <str>  		output keys of file name,forced
	-h     			Help

USAGE
	print $usage;
	exit;
}
