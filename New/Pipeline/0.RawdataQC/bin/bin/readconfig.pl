#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use utf8;
use Spreadsheet::ParseExcel;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$fqdir,$RunID);
GetOptions(
				"help|?" =>\&USAGE,
				"config:s"=>\$fIn,
				"outdir:s"=>\$dOut,
				"runID:s"=>\$RunID,
				) or &USAGE;
&USAGE unless ($fIn and $dOut and $RunID);
my %Barcode=(
	"10"=>"GGCTAC",
	"16"=>"ACCTCT",
	"18"=>"ATAATC",
	"46"=>"TGCTTA",
	"52"=>"TTCCAC",
	"53"=>"CTCC",
	"61"=>"AGGC",
	"62"=>"GATC",
	"63"=>"TCAC",
	"66"=>"TCACC",
	"78"=>"CATCT",
	"85"=>"TCGTT",
	"96"=>"ACGTGTT",
	"100"=>"CGCTGAT",
	"101"=>"CGGTAGA",
	"104"=>"TAGCGGA",
	"108"=>"ACGACTAC",
	"111"=>"TGCAAGGA",
	"114"=>"CCGGATAT",
	"119"=>"CCATGGGT",
	"121"=>"CGTGTGGT",
	"122"=>"GCTGTGGA",
	"123"=>"GGATTGGT",
	"124"=>"GTGAGGGT",
);
my %MSE=(
	"M1"=>"AACG",
	"M3"=>"TTACA",
	"M6"=>"CCTCAG",
	"M9"=>"GAAGATCC",
);
my %TAG=(
	"T1"=>"AGAACATA",
	"T5"=>"TTCAATG",
	"T7"=>"GAACTA",
	"T12"=>"CCTA",
);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir "$dOut/barcode" if (!-d "$dOut/barcode");
#open RAD,">$dOut/RAD.config";
my %barcode;
my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->Parse($fIn);
if ( !defined $workbook ) {
    die $parser->error(), ".\n";
}
my %config;
my %laneID;
my %out;
my @id;
for my $worksheet ( $workbook->worksheets() ) {
	my $name=$worksheet->get_name();
	if ( $name eq "文库表" || $name eq "Sheet1") {
		my ( $row_min, $row_max ) = $worksheet->row_range();
		my ( $col_min, $col_max ) = $worksheet->col_range();
		my $LaneCol;
		my $ProjCol;
		my $LibCol;
		my $TypeCol;
		my $SampleCol;
		my $NeadCol;
		my $DataCol;
		my $MajorCol;
		for (my $i=0;$i<=$col_max;$i++) {
			my $cell=$worksheet->get_cell(0,$i);
			last if (!defined $cell || $cell->value() eq "");
			if ($cell->value() eq "lane") {
				$LaneCol=$i;
			}elsif ($cell->value() eq "合同编号") {
				$ProjCol=$i;
			}elsif ($cell->value() eq "文库编号") {
				$LibCol=$i;
			}elsif ($cell->value() eq "样品名称") {
				$SampleCol=$i;
			}elsif ($cell->value() eq "合同数据量") {
				$NeadCol=$i;
			}elsif ($cell->value() eq "文库类型") {
				$TypeCol=$i;
			}elsif ($cell->value() eq "数据名") {
				$DataCol=$i;
			}elsif ($cell->value() eq "美吉编号") {
				$MajorCol=$i;
			}
		}
		for (my $i=1;$i<=$row_max;$i++) {
				my $cell=$worksheet->get_cell($i,$LaneCol);
				last if (!defined $cell || $cell->value() eq "");
				my $laneID=$cell->value();
				$cell=$worksheet->get_cell($i,$ProjCol);
				if (!defined $cell || $cell->value() eq "" || $cell->value() !~ /\D+\d+/){
					print $name,"\t",$i,"\t","ProjCol\n";
					next;
				}
				my $projectID=$cell->value();
				$cell=$worksheet->get_cell($i,$LibCol);
				if (!defined $cell || $cell->value() eq ""){
					print $name,"\t",$i,"\t","LibCol\n";
					next;
				}
				my $libID=$cell->value();
				$cell=$worksheet->get_cell($i,$SampleCol);
				if (!defined $cell || $cell->value() eq ""){
					print $name,"\t",$i,"\t","$SampleCol\n";
					next;
				}
				my $sampleID=$cell->value();
				$cell=$worksheet->get_cell($i,$MajorCol);
				my $majorID;
				if (!defined $cell || $cell->value() eq ""){
					$majorID="-";
				}else{
					$majorID=$cell->value();
					push @id,$majorID;
				}
				$cell=$worksheet->get_cell($i,$TypeCol);
				if (!defined $cell || $cell->value() eq ""){
					print $name,"\t",$i,"\t","TypeCol\n";
					next;
				}
				my $libtype="PE";
				if ($cell->value() =~ /RAD/) {
					$libtype="RAD";
				}elsif ($cell->value() =~ /GBS/) {
					$libtype="GBS";
				}
				$cell=$worksheet->get_cell($i,$NeadCol);
				if (!defined $cell || $cell->value() eq ""){
					print $name,"\t",$i,"\t","NeadCol\n";
					next;
				}
				push @{$laneID{$libID}},$laneID;
				my $sampleIDneed=$cell->value();
				my $dataID=$libID;
				if ($libtype =~ /PE/) {
					$config{join("\t",$laneID,$projectID,$libID,$libtype,$majorID,$sampleID,$sampleIDneed,"-","-")}=1;
				}
		}
	}elsif ($name eq "样品表"|| $name eq "Sheet2") {
		#上机位置	文库类型	合同编号	客户姓名	文库编号	美吉编号	样品名称	文库大小	数据量	数据单位	酶切浓度	样品体积	补水体积	Barcode	备注	"RawData
		my ( $row_min, $row_max ) = $worksheet->row_range();
		my ( $col_min, $col_max ) = $worksheet->col_range();
		my $LaneCol;
		my $ProjCol;
		my $LibCol;
		my $TypeCol;
		my $SampleCol;
		my $NeadCol;
		my $BarcodeCol;
		my $MajorCol;
		for (my $i=0;$i<=$col_max;$i++) {
			my $cell=$worksheet->get_cell(0,$i);
			last if (!defined $cell || $cell->value() eq "");
			if ($cell->value() eq "合同编号") {
				$ProjCol=$i;
			}elsif ($cell->value() eq "文库编号") {
				$LibCol=$i;
			}elsif ($cell->value() eq "样品名称") {
				$SampleCol=$i;
			}elsif ($cell->value() eq "数据量") {
				$NeadCol=$i;
			}elsif ($cell->value() eq "文库类型") {
				$TypeCol=$i;
			}elsif ($cell->value() eq "Barcode") {
				$BarcodeCol=$i;
			}elsif ($cell->value() eq "美吉编号") {
				$MajorCol=$i;
			}
		}
		for (my $i=1;$i<=$row_max;$i++) {
			my $cell=$worksheet->get_cell($i,$TypeCol);
			if (!defined $cell || $cell->value() eq ""){
				print $name,"\t",$i,"\t","TypeCol\n";
				next;
			}
			my $libtype;
			if ($cell->value() =~ /RAD/) {
				$libtype="RAD";
			}elsif ($cell->value() =~ /GBS/) {
				$libtype="GBS";
			}
			$cell=$worksheet->get_cell($i,$ProjCol);
			if (!defined $cell || $cell->value() eq "" || $cell->value() !~ /\D+\d+/){
				print $name,"\t",$i,"\t","ProjCol\n";
				next;
			}
			my $projectID=$cell->value();
			$cell=$worksheet->get_cell($i,$LibCol);
			if (!defined $cell || $cell->value() eq ""){
				print $name,"\t",$i,"\t","LibCol\n";
				next;
			}
			my $libID=$cell->value();
			$cell=$worksheet->get_cell($i,$SampleCol);
			if (!defined $cell || $cell->value() eq ""){
				print $name,"\t",$i,"\t","$SampleCol\n";
				next;
			}
			my $sampleID=$cell->value();
			my $majorID;
			$cell=$worksheet->get_cell($i,$MajorCol);
			if (!defined $cell || $cell->value() eq ""){
				$majorID="-";
			}else{
				$majorID=$cell->value();
				push @id,$majorID;
			}
			$cell=$worksheet->get_cell($i,$BarcodeCol);
			if (!defined $cell || $cell->value() eq ""){
				print $name,"\t",$i,"\t","BarcodeCol\n";
				next;
			}
			my $index=$cell->value();
			$cell=$worksheet->get_cell($i,$NeadCol);
			if (!defined $cell || $cell->value() eq ""){
				print $name,"\t",$i,"\t","NeadCol\n";
				next;
			}
			my $need=$cell->value();
			#if (!exists $barcode{$libID}) {
			#	open $barcode{$libID},">$dOut/barcode/$libID.barcode.txt";
			#}
			my $enzyme;
			if ($libtype =~ "GBS") {
				my $barcode=$index;
				$barcode=~/(T\d+)/;
				my $e1=$1;
				$barcode=~/(M\d+)/;
				my $e2=$1;
				$enzyme="TaqaI\tMseI";
				$barcode{$libID}{join("\t",$TAG{$e1},$MSE{$e2},$sampleID)}=1;;
				#print {$barcode{$libID}} $TAG{$e1},"\t",$MSE{$e2},"\t",$sampleID,"\n";
			}else{
				$index=~/(\d+)/;
				my $barcode=$1;
				$index=~/(\w+)/;
				$enzyme=$1;
				$enzyme="$enzyme\t$enzyme";
				#print $index,"\t",$barcode,"\t",$Barcode{$barcode};die;
				$barcode{$libID}{join("\t",$Barcode{$barcode},$sampleID)}=1;
				#print {$barcode{$libID}} $Barcode{$barcode},"\t",$sampleID,"\n";
			}
			next if (!exists $laneID{$libID});
			foreach my $laneID (@{$laneID{$libID}}) {
				$config{join("\t",$laneID,$projectID,$libID,$libtype,$majorID,$sampleID,$need,$enzyme,"$dOut/barcode/$libID.barcode.txt")}=1;
			}
			
		}
	}

}
foreach my $libID (sort keys %barcode) {
	next if (!exists $laneID{$libID});
	open Out,">$dOut/barcode/$libID.barcode.txt";
	print Out join("\n",keys %{$barcode{$libID}});
	close Out;
}
open Config,">$dOut/$RunID.config";
print Config "#RunID\tLaneID\tProjectID\tLibID\tLibType\tSampleID\tSampleNeed\tEnzyme\tBarcode\n";
foreach my $config (keys %config) {
	print Config "$RunID\t$config","\n";
}
close Config;

open OUT,">$dOut/sample.list";
print OUT join("\n",@id);
close OUT;
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
	-config	<file>	input run xls file [forced]
	-outdir	<dir>	output dir,[forced]
	-fqdir	<dir>	"/mnt/ilustre/upload/hiseq/hiseq4000/20170816nXten/"
	-run	<num>	runID [forced]
USAGE
	print $usage;
	exit;
}
