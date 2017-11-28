#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($pop,$out,$dsh,$maf,$mis,$dep,$gro,$sam);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"pop:s"=>\$pop,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"sam:s"=>\$sam,
			) or &USAGE;
&USAGE unless ($pop and $out and $dsh );
$pop=ABSOLUTE_DIR($pop);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$sam=ABSOLUTE_DIR($sam);
open SH,">$dsh/step01.structure1.sh";
open In,$pop;
open List,">$out/cv.list";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	next if (!/structure/);
	my (undef,$id,$bed)=split(/\s+/,$_);
	for (my $i=2;$i<=20;$i++) {
		print SH "cd $out && admixture $bed $i --cv -j8 > $out/$id.$i.log && paste $sam $out/$id.$i.Q > $out/$id.$i.xls && ";
		print SH "Rscript $Bin/bin/structure.R --infile $out/$id.$i.xls --outfile $out/$id.$i\n";
		print List "$id\t$i\t$out/$id.$i.log\t$out/$id.$i.xls\n";
	}
}
close List;
close In;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step01.structure1.sh";
`$job`;
open SH,">$dsh/step01.structure2.sh\n";
print SH "perl $Bin/bin/CVerror.pl -i $out/cv.list -o $out";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step01.structure2.sh";
`$job`;

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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -pop	<file>	input pop list 
  -out	<dir>	output dir
  -dsh	<dir>	output work shell
  -h         Help

USAGE
        print $usage;
        exit;
}
