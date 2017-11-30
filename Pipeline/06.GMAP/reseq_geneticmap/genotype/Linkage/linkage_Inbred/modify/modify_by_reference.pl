#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
$Script=~s/\.pl//g;
my @Times=localtime();
my $year=$Times[5]+1990;
my $month=$Times[4]+1;
my $day=$Times[3];
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fPosi,$dOut,$fLg,$fGeno,$Key,$sus,$D,$M,$Population,$diff_ratio);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fPosi,
				"l:s"=>\$fLg,
				"g:s"=>\$fGeno,
				"p:s"=>\$Population,
				"k:s"=>\$Key,
				"o:s"=>\$dOut,
				"sus:s"=>\$sus,
				"diff_ratio:s"=>\$diff_ratio,
				"D:s"=>\$D,
				"M:s"=>\$M,
				) or &USAGE;
&USAGE unless ($fPosi and $dOut and $fLg);
$diff_ratio||=0.85;
$Population||="CP";
$D=0.95;
$M=0.85;
$sus=0.06;
mkdir $dOut if (!-d $dOut);
$dOut=AbsolutePath("dir",$dOut);
open In,$fGeno;
my $Head;
my %Geno;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	if (/Type/ || /type/) {
		$Head = $_;
	}else{
		my ($id,$str)=split(/\s+/,$_);
		$Geno{$id}=$_;
	}
}
close In;

open In,$fPosi;
my %Pos;
my %info;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($ID,$CHR,$START,$END)=split(/\s+/,$_);
	next if (!exists $Geno{$ID});
	$END = $START if (!defined $END);
	my $POS=($START+$END)/2;
	if (exists $Pos{$CHR}{$POS}) {
		print $ID,"\n",$Pos{$CHR}{$POS},"\n";
		die;
	}
	$Pos{$CHR}{$POS}=$ID;
	$info{$ID}=$CHR;
}
close In;

open In,$fLg;
$/=">";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($id,$info)=split(/\n/,$_);
	$id=(split(/\s+/,$id))[0];
	my @ID=split(/\s+/,$info);
	my %exi;
	for (my $i=0;$i<@ID;$i++) {
		$exi{$ID[$i]}=1;
	}
	my $chr=$info{$ID[0]};
	mkdir "$dOut/primary/" if (!-d "$dOut/primary/");
	open Out,">$dOut/primary/$Key.$chr.genotype";
	open Map,">$dOut/primary/$Key.$chr.map";
	print Map "group 1\n";
	print Out $Head,"\n";
	foreach my $pos (sort {$a<=>$b} keys %{$Pos{$chr}}) {
		my $ID=$Pos{$chr}{$pos};
		next if (!exists $exi{$ID});
		print Out $Geno{$ID},"\n";
		print Map $ID,"\t",$pos,"\n";
	}
	close Out;
	close Map;
}
close In;
$/="\n";
mkdir "$dOut/smooth" if (!-d "$dOut/smooth");
if ($Population eq "CP") {
	my @genotype=glob("$dOut/primary/*.genotype");
	open SH,">$dOut/$Key.smooth.sh";
	foreach my $g (@genotype) {
		my $map=$g;
		$map=~s/\.genotype/\.map/g;
		my $key=basename($g);
		print SH "perl $Bin/bin/calculatePWDviaQsub.pl -i $g -d $dOut/pwd -k $key && perl $Bin/bin/linkagePhase.pl -p $dOut/pwd/$key.pwd -g $g -k $key -d $dOut/phase && perl $Bin/bin/smooth_CP.pl -m $map -l $dOut/phase/$key.loc -k $key -d $dOut/smooth -cycle_threshold 3  -sus $sus\n";
	}
	close SH;
}else{
	my @genotype=glob("$dOut/primary/*.genotype");
	open SH,">$dOut/$Key.smooth.sh";
	foreach my $g (@genotype) {
		my $map=$g;
		my @map=split(/\./,$map);
		$map[-1]="map";
		$map=join(".",@map);
		my $key=basename($g);
		print SH "perl $Bin/bin/smooth_Inbred.pl  -m $map -i $g -k $key -o $dOut/smooth -D $D -miss $M &&\n";
	}
	close SH;
}
#my $host=`hostname`;
#if ($host =~ /cluster/) {
	`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $dOut/$Key.smooth.sh -reqsub`;
#}else{
#	`ssh cluster -Y qsub-sge.pl $dOut/$Key.smooth.sh -reqsub`;
#}
if ($Population eq "CP") {
	my @loc=glob("$dOut/smooth/*.correct.loc");
	open Out,">$dOut/$Key.final.genotype";
	print Out "$Head\n";
	foreach my $l (@loc) {
		#sprint $l,"\n\n\n";
		open In,$l;
		while (<In>) {
			next if ($_ eq "" || /^$/ || /\=/ || /^;/);
			my @a=split(/\s+/,$_);
			$a[1]=~s/<//g;
			$a[1]=~s/>//g;
			print Out $a[0],"\t",$a[1],"\t",join("\t",@a[3..$#a]),"\n";
		}
		close In;
	}
	close Out;
}else{
	my @genotype=glob("$dOut/smooth/*.genotype");
	open Out,">$dOut/$Key.final.genotype";
	print Out "$Head\n";
	foreach my $l (@genotype) {
		open In,$l;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /^$/ || /\=/ || /type/ || /Type/);
			print Out $_,"\n";
		}
		close In;
	}
	close Out;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath{
        my ($type,$input) = @_;

        my $return;

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chop $return;
                $return .="/".$file;
                chdir($pwd);
        }
        return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[$month:$day:$year:]
	Contact:Huang Long <huangl\@biomarker.com.cn>
	Options:

		-i	<file>	input Pos file [ID	CHR	START	END]
		-l	<file>	linkage grouping file
		-g	<file>	genotyping file
		-p	<str>	population of Genotype default [CP];
		-o	<dir>	output dir
		-k	<str>	output keys of filename
		######################CP###########################

		-sus	<float>	error Rate default 0.06
		-diff_ratio	<floar>	default 0.085
		######################Inbreed######################

		-D	<float>  threshold values to change a loci, defaut [0.95]
		-M	<float>  minimim miss ratio, [0.85]

		-h	Help

USAGE
	print $usage;
	exit;
}
