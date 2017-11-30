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
my @Original_ARGV=@ARGV;

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fKey,$dOut,$reference,$modify,$fPosi,$onlyFirst,$start,$end,$stepSize,$minGroup,$maxGroup,$nChro,$minMarkerNum,$stepminGroup,$stepmaxGroup);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$dOut,
				"k:s"=>\$fKey,
				"c:s"=>\$fPosi,
				"n:s"=>\$nChro,
				"reference"=>\$reference,
				"modify:s"=>\$modify,
				"p:i"=>\$onlyFirst,
				"m:i"=>\$minMarkerNum,
				"b:s"=>\$start,
				"e:s"=>\$end,
				"s:s"=>\$stepSize,
				"minGroup:s"=>\$minGroup,
				"maxGroup:s"=>\$maxGroup,
				"stepminGroup:s"=>\$stepminGroup,
				"stepmaxGroup:s"=>\$stepmaxGroup,
				) or &USAGE;
&USAGE unless ($fIn and $fKey and $nChro);
$dOut||="./";
$fIn = AbsolutePath("file",$fIn);
mkdir $dOut if (!-d $dOut);
$dOut = AbsolutePath("dir",$dOut);
$fPosi = AbsolutePath("file",$fPosi) if (defined $fPosi);
$minMarkerNum||=400;
$onlyFirst||=200;
$modify ||=1;
$start||=3;        
$end||=20;         
$stepSize||=1;     
$minGroup||=20;    
$maxGroup||=1000;
$stepminGroup||=10;
$stepmaxGroup||=100;
open Log,">$dOut/$fKey.log";
print Log "############################################\n";
print Log "perl $Bin/$Script.pl"."\t".join("\t",@Original_ARGV)."\n";
print Log "############################################\n";
my $Step =0;
if ($Step == 0) {#### calculate mLOD for Linkage
	my $job = "perl $Bin/bin/calculate_mLOD_via_qsub.pl -i $fIn -d $dOut -k $fKey -p $onlyFirst -m $minMarkerNum 1";
	print Log "$job\n";
	`$job`;
	print Log "Done!\n";
	$Step ++;
}
if ($Step == 1) {
	mkdir "$dOut/tmp/" if (!-d "$dOut/tmp");
	my $job;
	if ($reference) {
		mkdir "$dOut/tmp/ref" if (!-d "$dOut/tmp/ref");
		$job="perl $Bin/bin/linkage_by_ref.pl -i $fIn -o $dOut/tmp/ref -k $fKey -2 $fPosi";
	}else{
		open SH,">$dOut/tmp/$fKey.linkgrouping.sh";
		while ($minGroup < $maxGroup) {
			print SH "perl $Bin/bin/linkage_by_mlod.pl -i $dOut/$fKey.mLOD -k $fKey\_$minGroup\_$maxGroup -d $dOut/tmp/$minGroup\_$maxGroup -n $nChro -b $start -e $end -s $stepSize -minGroup $minGroup -maxGroup $maxGroup &&\n";
			$minGroup += $stepminGroup;
			$maxGroup -= $stepmaxGroup;
		}
		close SH;
#		my $host=`hostname`;
#		if ($host =~ /cluster/) {
			$job = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $dOut/tmp/$fKey.linkgrouping.sh";
#		}else{
#			$job = "ssh cluster -Y qsub-sge.pl $dOut/tmp/$fKey.linkgrouping.sh";
#		}
	}
	print Log "$job\n";
	`$job`;
	print Log "Done!\n";
	$Step ++;


}
if ($Step == 2) {
	if (defined $reference && $modify == 1) {#有染色体水平的数据的校正
		mkdir "$dOut/tmp" if (!-d "$dOut/tmp");
		my @Group=glob("$dOut/tmp/ref/*.genotype");
		open SH,">$dOut/tmp/$fKey.modify.sh";
		foreach my $Gm (@Group) {
			my $maxGroup=`wc -l $Gm`;
			chomp $maxGroup;
			$maxGroup --;
			my $minGroup = int($maxGroup * 0.7);
			my $LG=basename($Gm);
			print SH "perl $Bin/bin/linkage_by_mlod.pl -g $Gm -i $dOut/$fKey.mLOD -k $LG -d $dOut/tmp/$LG -n 1 -b 4 -e $end -s $stepSize -minGroup $minGroup -maxGroup $maxGroup &&\n";
		}
		close SH;
		my $host=`hostname`;
		my $job;
#		if ($host =~ /cluster/) {
			$job = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc 20 --interval 10 --resource vf=5G  $dOut/tmp/$fKey.modify.sh";
#		}else{
#			$job = "ssh cluster -Y qsub-sge.pl $dOut/tmp/$fKey.modify.sh";
#		}
		` $job `;

		my @LG=glob("$dOut/tmp/*/*.*.lg");
		open Out,">$dOut/$fKey.final.lg";
		foreach my $LG (@LG) {
			my $LGID=(split(/\./,basename($LG)))[1];
			open In,$LG;
			$/=">";
			while (<In>) {
				chomp;
				next if ($_ eq "" || /^$/);
				my ($ID,$info)=split(/\n/,$_);
				my ($id,$line)=split(/\s+/,$ID,2);
				print Out ">$LGID\t$line\n$info\n";
			}
			close In;
			$/="\n";
		}
		close Out;
	}elsif (defined $fPosi && -f $fPosi) {
		mkdir "$dOut/all_scheme" if (!-d "$dOut/all_scheme");
		`ln -s $dOut/tmp/*/all_scheme/*.lg $dOut/all_scheme`;
		my @lg = glob("$dOut/all_scheme/*.lg");
		my %info;	
		my %sum;
		foreach my $lg (@lg) {
			open In,$lg;
			$/=">";
			while (<In>){
				chomp;
				next if ($_ eq "" || /^$/);
				my ($id,$info)=split(/\n/,$_);
				my $lgID=(split(/\s+/,$id))[0];
				my @info=split(/\s+/,$info);
				$sum{$lg}+=scalar @info;
				for (my $i=0;$i<@info;$i++) {
					$info{$lg}{$info[$i]}=$lgID;
				}
			}
			$/="\n";
			close In;
		}
		open In,$fPosi;
		my %posi;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /^$/);
			my ($marker,$chr,$start,$end)=split(/\s+/,$_);
			$posi{$chr}{$marker}=1;
		}
		close In;
		my %lgerror;
		my %stat;
		my %info_test;
		foreach my $lg (sort keys %info) {
			my %stat;
			foreach my $chr (sort keys %posi) {
				foreach my $marker (sort keys %{$posi{$chr}}) {
					next if (!exists $info{$lg}{$marker});
					$stat{$info{$lg}{$marker}}++;
					$info_test{$lg}{$info{$lg}{$marker}}{$chr}++;
				}
				next if (scalar keys %stat <= 1);
				my ($a,$b)=sort{$a<=>$b} values %stat;
				$lgerror{$lg}+=$a;
			}
		}
		open Out,">$dOut/lg_pick.check";
		foreach my $lg (sort keys %info_test) {
			print Out ">$lg\t$sum{$lg}\n";
			foreach my $lgID (sort keys %{$info_test{$lg}}) {
				print Out $lgID,"\t";
				print Out join("\t",sort{$info_test{$lg}{$lgID}{$a}<=>$info_test{$lg}{$lgID}{$b}} keys %{$info_test{$lg}{$lgID}});
				print Out "\t",join("\t",sort {$a<=>$b} values %{$info_test{$lg}{$lgID}}),"\n";
			}
		}
		close Out;
		if (scalar keys %lgerror == 0) {
			my $lg_pick=(sort{$sum{$a}<=>$sum{$b}}keys %sum)[0];
			`cp $lg_pick $dOut/$fKey.lg`;
		}else{
			my $lg_pick=(sort { $lgerror{$a} <=> $lgerror{$b}} keys %lgerror)[0];
			`cp $lg_pick $dOut/$fKey.lg`;
		}
	}else{
		mkdir "$dOut/all_scheme" if (!-d "$dOut/all_scheme");
		`ln -s $dOut/tmp/*/all_scheme/*.lg $dOut/all_scheme`;
		my @lg=glob("$dOut/all_scheme/*.lg");
		my %lg;								  
		foreach my $lg (@lg) {

			open In,$lg;
			$/=">";
			while (<In>){
				chomp;
				next if ($_ eq "" || /^$/);
				my ($id,$info)=split(/\n/,$_);
				my $lgID=(split(/\s+/,$id))[0];
				$lg{$lg}+=(split(/\s+/,$id))[6];
			}
			$/="\n";
			close In;
		}
		my $lg_pick=(sort{$lg{$a}<=>$lg{$b}} keys %lg)[0];
		`cp $lg_pick $dOut/$fKey.lg`;
	}
}

close Log;













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
	Version:$version	
	Contact:Huang Long <huangl\@biomarker.com.cn>
	###############################################
	
	Note: if you use reference to grouping,please make sure your reference is on CHROMOSOME level!!!
	
	###############################################
	Options:
		-i	<file>	input Genotype file	[forced]
		-o	<dir>	output dir	[forced]
		-k	<str>	output keys of filename	[forced]
		-c	<file>	input posi file	[optinal]

		################reference###########################
		-reference
		
			-modify	<num>	modify genome result default 1
		
		################no reference########################

		-n	<int>   Species\' chromosome number,forced

		################grouping parameter#############################

		-p	<int>	Integer skip number marker at split file, default 200
		-m	<int>	Minimum threshold of markers\' number for qsub, optional, default [400]

		-b	<int>   Start lod,optional,default [3]
		-e	<int>   End lod,optional,default [20]
		-s	<int>   Stepsize of lod,optional,default [1]

		-minGroup     <int>   Minimum threshold of group size,optional,default [20]
		-maxGroup     <int>   Maximum threshold of group size,optional,default [1000]
		
		-stepminGroup	<int>	step of minGroup default 10
		-stepmaxGroup	<int>	step of maxGroup default 100

		-h	Help

USAGE
	print $usage;
	exit;
}
