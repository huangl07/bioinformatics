#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$group,$env,$enum,$dir,$dsh,$pnum);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
    "group:s"=>\$group,
    "pnum:s"=>\$pnum,
    "env:s"=>\$env,
    "enum:s"=>\$enum,
	"dir:s"=>\$dir,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($vcf and $dir and $pnum and $enum and $dsh);
mkdir $dir if (!-d $dir);
mkdir $dsh if (!-d $dsh);
$dir=ABSOLUTE_DIR($dir);
$dsh=ABSOLUTE_DIR($dsh);
$vcf=ABSOLUTE_DIR($vcf);
$group=ABSOLUTE_DIR($group);

open OUT,">$dir/VCF_BAYENV.spid";
if($group){
    print OUT "PARSER_FORMAT=VCF\n";
    print OUT "VCF_PARSER_QUAL_QUESTION=\n";
    print OUT "VCF_PARSER_POP_FILE_QUESTION=$group\n";
    print OUT "VCF_PARSER_PLOIDY_QUESTION=\n";
    print OUT "VCF_PARSER_POP_QUESTION=TRUE\n";
    print OUT "VCF_PARSER_GTQUAL_QUESTION=\n";
    print OUT "VCF_PARSER_MONOMORPHIC_QUESTION=Faulse\n";
    print OUT "VCF_PARSER_IND_QUESTION=\n";
    print OUT "VCF_PARSER_REGION_QUESTION=\n";
    print OUT "VCF_PARSER_READ_QUESTION=\n";
    print OUT "VCF_PARSER_PL_QUESTION=\n";
    print OUT "VCF_PARSER_EXC_MISSING_LOCI_QUESTION=\n";
    print OUT "WRITER_FORMAT=BAYENV\n";
    print OUT "BAYENV_WRITER_SAMPLE_FILE_QUESTION=\n";
    print OUT "BAYENV_WRITER_HALF_MISSING_QUESTION=\n";
    print OUT "BAYENV_WRITER_WRITE_INFO_FILE_QUESTION=ture\n";
    print OUT "BAYENV_WRITER_WRITE_SAMPLE_FILE_QUESTION=\n";
    print OUT "BAYENV_WRITER_INFO_FILE_QUESTION=SNP\n";
    }
close OUT;
`java -jar /mnt/ilustre/users/dna/.env/bin/PGDSpider2-cli.jar -inputfile $vcf -inputformat VCF -outputfile $dir/snpsfile.table -outputformat BAYENV -spid $dir/VCF_BAYENV.spid`;
open SH,">$dsh/bayenv2.sh";
my $maxnum=$pnum+1;
print SH "cd $dir && bayenv2 -i snpsfile.table -p $pnum -k 100000 -r 63479 > matrix.table && tail -$maxnum $dir/matrix.table | head -$pnum > matrix.$pnum && cd $dir && sh $Bin/bin/run.bayenv.sh snpsfile.table $env matrix.$pnum $pnum 10000 $enum && perl $Bin/bin/bayenv.result.pl -bf $dir/pop.bf -out $dir \n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/bayenv2.sh";
`$job`;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input vcf files;must be filtered
  -group <file> input group list;split by \\t
  -env <file> env file
  -enum <number> env number
  -pnum <number> pop number
  -dir	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
