#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$pop,$dIn,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"in:s"=>\$dIn,
	"out:s"=>\$out,
	"pop:s"=>\$pop,
	"dsh:s"=>\$dsh,
) or &USAGE;
&USAGE unless ($out and $dIn);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:
        fq thanslate to fa format
        eg:
        perl $Script -i -o -k -c
Usage:
  Options:
  -in	<dir>	input dir
  -out	<dir>	output dir
  -pop	<str>	pop type
  -dsh	<dir>	worksh
        
  -h         Help

USAGE
        print $usage;
        exit;
}

mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
#$pop||="F2";
#$out=ABSOLUTE_DIR($out);
#$vcf=ABSOLUTE_DIR($vcf);
#$ann=ABSOLUTE_DIR($ann);
open SH,">$dsh/step03.merge.sh";
if($pop=~/CP/){
	mkdir "$out/male" if (!-d "$out/male");
	mkdir "$out/female" if (!-d "$out/female");
	mkdir "$out/sexAver" if (!-d "$out/sexAver");
	mkdir "$out/male/Result" if(!-d "$out/male/Result");
	mkdir "$out/male/Figure" if (!-d "$out/male/Figure");
	mkdir "$out/male/Table" if (!-d "$out/male/Table");
	mkdir "$out/female/Result" if(!-d "$out/female/Result");
	mkdir "$out/female/Figure" if(!-d "$out/female/Figure");
	mkdir "$out/female/Table" if (!-d "$out/female/Table");
	mkdir "$out/sexAver/Result" if(!-d "$out/sexAver/Result");
	mkdir "$out/sexAver/Figure" if(!-d "$out/sexAver/Figure");
	mkdir "$out/sexAver/Table" if (!-d "$out/sexAver/Table");
	print SH "perl $Bin/bin/merge.qtl.pl -i $dIn/01.qtl/male/ -o $out/male/3-18.xls \n";
	print SH "perl $Bin/bin/merge.vcf.pl -i $dIn/02.anno/male/ -o $out/male/3-19.xls \n";
	print SH "perl $Bin/bin/merge.gene.pl -i $dIn/02.anno/male/ -o $out/male/3-20.xls \n";
	print SH "perl $Bin/bin/merge.qtl.pl -i $dIn/01.qtl/female/ -o $out/female/3-18.xls \n";
	print SH "perl $Bin/bin/merge.vcf.pl -i $dIn/02.anno/female/ -o $out/female/3-19.xls \n";
	print SH "perl $Bin/bin/merge.gene.pl -i $dIn/02.anno/female/ -o $out/female/3-20.xls \n";
	print SH "perl $Bin/bin/merge.qtl.pl -i $dIn/01.qtl/sexAver/ -o $out/sexAver/3-18.xls \n";
	print SH "perl $Bin/bin/merge.vcf.pl -i $dIn/02.anno/sexAver/ -o $out/sexAver/3-19.xls \n";
	print SH "perl $Bin/bin/merge.gene.pl -i $dIn/02.anno/sexAver/ -o $out/sexAver/3-20.xls \n";
	system("ln -s $dIn/01.qtl/male/*.qtl.p* $out/male/Figure/ \n");
	system("ln -s $dIn/01.qtl/male/*.scan.p* $out/male/Figure/ \n");
	system("ln -s $dIn/02.anno/male/*.png $out/male/Figure/ \n");
	system("ln -s $dIn/02.anno/male/*.pdf $out/male/Figure/ \n");
	system("ln -s $dIn/01.qtl/male/*.qtl.csv $out/male/Result/ \n");
	system("ln -s $dIn/01.qtl/male/*.scan.csv $out/male/Result/ \n");
	system("ln -s $dIn/02.anno/male/*.total $out/male/Result/ \n");
	system("ln -s $dIn/02.anno/male/*.eff $out/male/Result/ \n");
	system("ln -s $dIn/01.qtl/female/*.qtl.p* $out/female/Figure/ \n");
	system("ln -s $dIn/01.qtl/female/*.scan.p* $out/female/Figure/ \n");
	system("ln -s $dIn/02.anno/female/*.png $out/female/Figure/ \n");
	system("ln -s $dIn/02.anno/female/*.pdf $out/female/Figure/ \n");
	system("ln -s $dIn/01.qtl/female/*.qtl.csv $out/female/Result/ \n");
	system("ln -s $dIn/01.qtl/female/*.scan.csv $out/male/Result/ \n");
	system("ln -s $dIn/02.anno/female/*.total $out/female/Result/ \n");
	system("ln -s $dIn/02.anno/female/*.eff $out/female/Result/ \n");
	system("ln -s $dIn/01.qtl/sexAver/*.scan.p* $out/sexAver/Figure/ \n");
	system("ln -s $dIn/01.qtl/sexAver/*.scan.p* $out/sexAver/Figure/ \n");
	system("ln -s $dIn/02.anno/sexAver/*.png $out/sexAver/Figure/ \n");
	system("ln -s $dIn/02.anno/sexAver/*.pdf $out/sexAver/Figure/ \n");
	system("ln -s $dIn/01.qtl/sexAver/*.scan.csv $out/sexAver/Result/ \n");
	system("ln -s $dIn/01.qtl/sexAver/*.qtl.csv $out/sexAver/Result/ \n");
	system("ln -s $dIn/02.anno/sexAver/*.total $out/sexAver/Result/ \n");
	system("ln -s $dIn/02.anno/sexAver/*.eff $out/sexAver/Result/ \n");
}else{
	mkdir "$out/Result" if (!-d "$out/Result");
	mkdir "$out/Figure" if (!-d "$out/Figure");
	mkdir "$out/Table" if (!-d "$out/Table");
	mkdir "$out/Result/qtl" if (!-d "$out/Result/qtl");
	mkdir "$out/Result/vcf" if (!-d "$out/Result/vcf");
	mkdir "$out/Result/gene" if (!-d "$out/Result/gene");
	mkdir "$out/Result/anno" if (!-d "$out/Result/anno");
	mkdir "$out/Figure/qtl" if (!-d "$out/Figure/qtl");
	mkdir "$out/Figure/anno" if (!-d "$out/Figure/anno");
	print SH "perl $Bin/bin/merge.qtl.pl -i $dIn/01.qtl/ -o $out/Table/3-15.xls \n";
	print SH "perl $Bin/bin/merge.vcf.pl -i $dIn/02.anno/ -o $out/Table/3-16.xls \n";
	print SH "perl $Bin/bin/merge.gene.pl -i $dIn/02.anno/ -o $out/Table/3-17.xls \n";
	system("ln -s $dIn/01.qtl/*.scan.p* $out/Figure/qtl/ \n");
	system("ln -s $dIn/01.qtl/*.qtl.p* $out/Figure/qtl/ \n");
	system("ln -s $dIn/02.anno/*.png $out/Figure/anno/ \n");
	system("ln -s $dIn/02.anno/*.pdf $out/Figure/anno/ \n");
	system("ln -s $dIn/01.qtl/*.qtl.csv $out/Result/qtl/ \n");
	system("ln -s $dIn/01.qtl/*.scan.csv $out/Result/qtl/ \n");
	system("ln -s $dIn/02.anno/*.vcf.total $out/Result/vcf");
	system("ln -s $dIn/02.anno/*.vcf.eff $out/Result/vcf");
	system("ln -s $dIn/02.anno/*.gene.total $out/Result/gene/");
	system("ln -s $dIn/02.anno/*.gene.eff $out/Result/gene/");
	system("ln -s $dIn/02.anno/*.stat $out/Result/anno/"); 
}
my $job="sh $dsh/step03.merge.sh \n";	
`$job`;
print "job done\n";

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#$pavfile=&ABSOLUTE_DIR($pavfile);
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

