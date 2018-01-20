#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
my $BEGIN_TIME=time();
my $version="1.2.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fPos,$fXls,$dOut,$Key);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"p:s"=>\$fPos,
				"x:s"=>\$fXls,
				"d:s"=>\$dOut,
				"k:s"=>\$Key,
				) or &USAGE;
&USAGE unless ($fIn and $dOut and $Key and $fXls and $fPos);
mkdir $dOut if (!-d $dOut);
mkdir "$dOut/tmp" if (!-d "$dOut/tmp");
open (MAPIN,"$fIn") || die "Can't open $fIn\n" ;
open (POSIN,"$fPos") || die "Can't open $fPos\n" ;
open (XLSIN,"$fXls") || die "Can't open $fXls\n" ;
open CMBOUT,">$dOut/tmp/$Key.mappos";
open BINOUT,">$dOut/tmp/$Key.bincount";#table 3-15
open DISOUT,">$dOut/tmp/$Key.discount";#table 3-16
open MISSOUT,">$dOut/tmp/$Key.misscount";#table 3-17

print BINOUT "LG\tSNP markers\tBin markers\tAverage length(Mb)\tAverage Gene Numbers\n";
print DISOUT "LG ID\tbinNum\tTotal cM\tAver cM\tMax Gap\n";
print MISSOUT "Group ID\tMissing Precent(%)\n";
my (%bintable, %distable);#%bintable for 3-15, distable for 3-16/17
my ($groupnumber,$gap,$premarkercm);

##combine pos file and map file
my %table;
my $groupnumber2;
open (MAPTMPIN,"$fIn") || die "Can't open $fIn\n" ;
while (<MAPTMPIN>){
	chomp;
        if (/^group (.*)/){
                $groupnumber2 = $1;
        }else{
        	chomp;
	        my($id,$cm) = split /\t/;
      		$table{$id}{'id'} = $id;
        	$table{$id}{'cm'} = $cm;
        	$table{$id}{'lg'} = $groupnumber2;
	}
}
close MAPTMPIN;

while (<POSIN>){
	chomp;
	next if /^#MarkerID/;
	my ($id,$chr,$pos) = split /\t/;
	$table{$id}{'pos'} = $pos if exists $table{$id}{'id'};
}
close POSIN;

for (sort idorder keys %table){
        print CMBOUT "$table{$_}{'id'}\t$table{$_}{'lg'}\t$table{$_}{'cm'}\t$table{$_}{'pos'}\n";
}
close CMBOUT;

my %tmp3;
while (<MAPIN>) {##count binNum, this can be overtaken by triple dim Hash like %tmp1
	chomp;
	if (/^group (.*)/){
		$groupnumber = $1;
		$distable{$groupnumber}{'maxgap'}=0;##reset
		$bintable{$groupnumber}{'Binmarkers'}=0;		
		$distable{$groupnumber}{'binNum'}=0;
		$premarkercm =0;##reset
	}else{
		my($markerid, $markercm) = split /\t/;
		$bintable{$groupnumber}{'SNPmarkers'}++;
		if ($tmp3{$groupnumber}{$markercm}++){#get bin number, singleton counts 
			$bintable{$groupnumber}{'Binmarkers'}++;		
			$distable{$groupnumber}{'binNum'}++;
		} 
		$gap=$markercm-$premarkercm;
		$distable{$groupnumber}{'maxgap'}= $gap if $gap > $distable{$groupnumber}{'maxgap'};
		$premarkercm = $markercm;
	}
}

###count bin length for each cm
open CMBIN,"$dOut/tmp/$Key.mappos";
my (%tmp1, %tmp2);
my $dispos;
my (@sumpos, @sortpos);
while (<CMBIN>){
	chomp;
	my($id,$lg,$cm,$pos)=split /\t/;
	push @{$tmp1{$lg}{$cm}}, $pos;
}
close CMBIN;

for my $lg(keys %tmp1){
	for my $cm(keys %{$tmp1{$lg}}){
		@sortpos = sort {$a<=>$b} @{$tmp1{$lg}{$cm}};
		$dispos = $sortpos[-1] - $sortpos[0];
		#print "$sortpos[-1] - $sortpos[0] = $dispos \n";
		$tmp2{$lg} += $dispos;
	}
	print "group $lg has dis with $tmp2{$lg}\n";
	#print "$bintable{$lg}{'Binmarkers'}==\n";
	if ($bintable{$lg}{'Binmarkers'} != 0){
		$bintable{$lg}{'averLength'} = $tmp2{$lg}/(1000000*$bintable{$lg}{'Binmarkers'});#Mb
	}else{
		$bintable{$lg}{'averLength'} = 0;
	}
}

<XLSIN>;
while (<XLSIN>){
	chomp;
	my($LG,$nloc,$nind,$sum,$miss,$single,$noise,$missp,$singlep,$noisep,$gap5,$gap5p,$dist,$avgdist) = split /\t/;
	$distable{$LG}{'TotalcM'}=$dist;
	if($distable{$LG}{'binNum'}!=0){
		$distable{$LG}{'avercM'}=$dist/$distable{$LG}{'binNum'}
	}else{
		$distable{$LG}{'avercM'}=0;
	}
	$distable{$LG}{'MissingPercent'}=$missp;
}

for (sort keys %bintable){##table 3-15 output
	print BINOUT "$_\t$bintable{$_}{'SNPmarkers'}\t$bintable{$_}{'Binmarkers'}\t$bintable{$_}{'averLength'}\t--\n";
}

my($bintotal,$cmtotal,$avgcmtotal);
my $maxgaptotal=0;
for (sort keys %distable){##table 3-16 3-17 output
	printf DISOUT "%s\t%d\t%.3f\t%.6f\t%.2f\n",$_,$distable{$_}{'binNum'},$distable{$_}{'TotalcM'},$distable{$_}{'avercM'},$distable{$_}{'maxgap'};
	printf MISSOUT "%s\t%s\n",$_,$distable{$_}{'MissingPercent'};
	$bintotal+=$distable{$_}{'binNum'};
	$cmtotal+=$distable{$_}{'TotalcM'};
	$maxgaptotal = $distable{$_}{'maxgap'} if $distable{$_}{'maxgap'} > $maxgaptotal;
}
$avgcmtotal=$cmtotal/$bintotal;
printf DISOUT "Total\t%d\t%.3f\t%.6f\t%.2f\n",$bintotal,$cmtotal,$avgcmtotal,$maxgaptotal;

close XLSIN;
close MAPIN;
close BINOUT;
close DISOUT;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub idorder {
        (my $da = $a) =~ s/\D+//g;
        (my $db = $b) =~ s/\D+//g;
        $da <=> $db;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: binsum.pl for table 3-15 3-16 3-17
Version: $version
Contact: Ao.li

Usage:
  Options:
	-i <file>		input map file, forced
	-p <file>		input pos file, forced
	-x <file>		input \$Key.MSTmap.stat.xls file, forced
	-d <file>		output dir,forced
	-k <str>  		output keys of file name,forced
	-h     			Help

USAGE
	print $usage;
	exit;
}
