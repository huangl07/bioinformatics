#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output,$annotation,$level);
use Data::Dumper;
use GO::OntologyProvider::OboParser;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"obofile:s"=>\$input,
	"annotation:s"=>\$annotation,
	"output:s"=>\$output,
	"level:s"=>\$level,
			) or &USAGE;
&USAGE unless ($input and $output);
$level||=2;
my $process   = GO::OntologyProvider::OboParser->new(ontologyFile => $input,
                                                     aspect       => 'P');
my $component = GO::OntologyProvider::OboParser->new(ontologyFile => $input,
                                                     aspect       => 'C');
my $function  = GO::OntologyProvider::OboParser->new(ontologyFile => $input,
                                                     aspect       => 'F');
my @pnodes = $process->allNodes;
my @cnodes = $component->allNodes;
my @fnodes = $function->allNodes;
my %terms;
foreach my $nodes (@pnodes) {
	$terms{$nodes->goid}="biological_process";
}
foreach my $nodes (@cnodes) {
	$terms{$nodes->goid}="cellular_component";
}
foreach my $nodes (@fnodes) {
	$terms{$nodes->goid}="molecular_function";
}
$/="\n";
open In,$annotation;
my %stat;
my %idterm;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^#/);
	my ($geneId,$goID,$term)=split(/\t/,$_);
	next if ($goID eq "--");
	my @goid=split(/\,/,$goID);
	foreach my $goid (@goid) {
		my $node;
		next if (!exists $terms{$goid} || $goid eq "GO:0003674" || $goid eq "GO:0005575" ||$goid eq "GO:0008150");
		if ($terms{$goid} eq "biological_process") {
			$node= $process->nodeFromId($goid);
		}elsif ($terms{$goid} eq "molecular_function") {
			$node= $function->nodeFromId($goid);
		}elsif ($terms{$goid} eq "cellular_component") {
			$node= $component->nodeFromId($goid);
		}
		my @pathsToRoot = $node->pathsToRoot;
		next if ($node->meanLengthOfPathsToRoot < $level+1);
		foreach my $path (@pathsToRoot) {
			my $goid1=$path->[$level]->goid;
			$stat{$goid1}{$geneId}=1;
			$idterm{$goid1}=$path->[$level]->term;
		}
	}
}
close In;
open Out,">$output";
print Out "goid\tgoterm\tis_a\tgenes\n";
foreach my $goid (sort {$terms{$a} cmp $terms{$b}} keys %stat) {
	print Out $goid,"\t",$idterm{$goid},"\t",$terms{$goid},"\t",join(",",keys %{$stat{$goid}}),"\n";
}
close Out;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
Usage:
  Options:
  -input	<file>	input reference file name
  -output	<file>	input gff file name
  -h         Help

USAGE
        print $usage;
        exit;
}
