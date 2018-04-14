#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$ann,$pop,$btl,$vcf,$qtl,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"out:s"=>\$out,
	"ann:s"=>\$ann,
	"vcf:s"=>\$vcf,
	"pop:s"=>\$pop,
	"qtl:s"=>\$qtl,
	"dsh:s"=>\$dsh,
) or &USAGE;
&USAGE unless ($out and $ann and $vcf and $qtl);
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
  -out	<dir>	output dir
  -ann	<file>	input ann file
  -vcf	<file>  input vcf file
  -pop	<str>	pop type
  -qtl	<dir>	qtl dir(scan,qtl.csv)
  -dsh	<dir>	worksh
        
  -h         Help

USAGE
        print $usage;
        exit;
}

mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$pop||="F2";
#$out=ABSOLUTE_DIR($out);
#$vcf=ABSOLUTE_DIR($vcf);
#$ann=ABSOLUTE_DIR($ann);
mkdir "$out/vcf" if (!-d "$out/vcf");
mkdir "$out/gene" if (!-d "$out/gene");
mkdir "$out/stat" if (!-d "$out/stat");
my(@out,$trt,$job);
if($pop=~/CP/){
	mkdir "$out/male" if (!-d "$out/male");
	mkdir "$out/female" if (!-d "$out/female");
	mkdir "$out/sexAver" if (!-d "$out/sexAver");
	open SH1,">$dsh/step02.anno.male.sh";
	open SH2,">$dsh/step02.anno.female.sh";
	open SH3,">$dsh/step02.anno.sexAver.sh";
	@out=glob("$qtl/male/*.qtl.csv");
	for(@out){
		$trt=(split(/\//,$_))[-1];
		$trt=~s/\.qtl\.csv//g;
		print SH1 "perl $Bin/bin/region-vcf.pl -i $vcf -o $out/male/$trt.vcf -r $_ && ";
		print SH1 "perl $Bin/bin/region-gene-new.anno.pl -a $ann -o $out/male/$trt.gene -i $_ && ";
		print SH1 "Rscript $Bin/bin/eff-enrich.R --input $out/male/$trt.gene.kegg.stat --output $out/male/$trt.gene.kegg.stat  --top 1 && ";
		print SH1 "Rscript $Bin/bin/eff-enrich.R --input $out/male/$trt.gene.gene.go.stat --output $out/male/$trt.gene.go.stat --top 1&& ";
		print SH1 "Rscript $Bin/bin/eff-enrich.R --input $out/male/$trt.gene.gene.eggnog.stat --output $out/male/$trt.gene.eggnog.stat --eggnog \n";
	}
	@out=glob("$qtl/female/*.qtl.csv");
	for (@out){
		$trt=(split(/\//,$_))[-1];
		$trt=~s/\.qtl\.csv//g;
		print SH2 "perl $Bin/bin/region-vcf.pl -i $vcf -o $out/female/$trt.vcf -r $_ && ";
		print SH2 "perl $Bin/bin/region-gene-new.anno.pl -a $ann -o $out/female/$trt.gene -i $_ && ";
		print SH2 "Rscript $Bin/bin/eff-enrich.R --input $out/female/$trt.gene.kegg.stat --output $out/female/$trt.gene.kegg.stat  --top 1 && ";
		print SH2 "Rscript $Bin/bin/eff-enrich.R --input $out/female/$trt.gene.go.stat --output $out/female/$trt.gene.go.stat --top 1&& ";
		print SH2 "Rscript $Bin/bin/eff-enrich.R --input $out/female/$trt.gene.eggnog.stat --output $out/female/$trt.gene.eggnog.stat --eggnog \n";
	}
	@out=glob("$qtl/sexAver/*.qtl.csv");
	for (@out){
		$trt=(split(/\//,$_))[-1];
		$trt=~s/\.qtl\.csv//g;
		print SH3 "perl $Bin/bin/region-vcf.pl -i $vcf -o $out/sexAver/$trt.vcf -r $_ && ";
		print SH3 "perl $Bin/bin/region-gene-new.anno.pl -a $ann -o $out/sexAver/$trt.gene -i $_ && ";
		print SH3 "Rscript $Bin/bin/eff-enrich.R --input $out/sexAver/$trt.gene.kegg.stat --output $out/sexAver/$trt.gene.kegg.stat  --top 1 && ";
		print SH3 "Rscript $Bin/bin/eff-enrich.R --input $out/sexAver/$trt.gene.go.stat --output $out/sexAver/$trt.gene.go.stat --top 1 && ";
		print SH3 "Rscript $Bin/bin/eff-enrich.R --input $out/sexAver/$trt.gene.eggnog.stat --output $out/sexAver/$trt.gene.eggnog.stat --eggnog \n";
	}
	close SH1;
	close SH2;
	close SH3;
	my $newsh=`cat $dsh/step02.anno.sexAver.sh $dsh/step02.anno.male.sh $dsh/step02.anno.female.sh`;
	$job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=12G  --CPU 1 --Maxjob 30 $newsh \n";
}else{
	#mkdir "$out/vcf" if (!-d "$out/vcf");
	#mkdir "$out/gene" if (!-d "$out/gene");
	#mkdir "$out/stat" if (!-d "$out/stat");
	open SH,">$dsh/step02.anno.sh";
	my @out=glob("$qtl/*.qtl.csv");
	for(@out){
		my $trt=(split(/\//,$_))[-1];
		$trt=~s/\.qtl\.csv//g;
		#print SH "perl $Bin/bin/region-variant.pl -i $vcf -o $out/variant/$trtname.qtl.vcf -r $trt && ";
		print SH "perl $Bin/bin/region-vcf.pl -i $vcf -o $out/$trt.vcf -r $_ \n ";
		print SH "perl $Bin/bin/region-gene-new.anno.pl -a $ann -o $out/$trt.gene -i $_ && ";
		print SH "Rscript $Bin/bin/eff-enrich.R --input $out/$trt.gene.kegg.stat --output  $out/$trt.gene.kegg.stat --top 1 && ";
		print SH "Rscript $Bin/bin/eff-enrich.R --input  $out/$trt.gene.gene.go.stat --output  $out/$trt.gene.go.stat --top 1&& ";
		print SH "Rscript $Bin/bin/eff-enrich.R --input  $out/$trt.gene.gene.eggnog.stat --output  $out/$trt.gene.eggnog.stat --eggnog \n";
	}
	close SH;
	$job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=12G  --CPU 1 --Maxjob 30 $dsh/step02.anno.sh \n";
}

#$job.="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=12G  --CPU 1 --Maxjob 30 $dsh/step02.anno.sh \n";

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

