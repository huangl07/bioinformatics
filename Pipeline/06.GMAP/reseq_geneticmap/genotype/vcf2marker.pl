#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"i=s","o=s","P=s","M=s","Pd=s","Md=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || !defined($opts{P}) || !defined($opts{M}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:
		-i           vcf file              <infile>     must be given
		-o           result file           <outfile>    must be given
		-P           paternal sample ID    <string>     must be given
		-M           maternal sample ID    <string>     must be given
		-Pd          paternal min depth    [int]        optional [4]
		-Md          maternal min depth    [int]        optional [4]
		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
## get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;
my $paternal_ID = $opts{P} ;
my $maternal_ID = $opts{M} ;
my $paternal_depth = defined $opts{Pd} ? $opts{Pd} : 4 ;
my $maternal_depth = defined $opts{Md} ? $opts{Md} : 4 ;

## reading vcf file and 
open (IN,"$infile") || die "Can't open $infile,$!\n" ;
open (OUT,">$outfile") || die "Can't creat $outfile,$!\n" ;
my @culumns = ();
my %hcul2sam = ();
my ($p_cul, $m_cul) ;
my $total_marker = 0 ;
my $marker_num = 0 ;
while (<IN>){
	chomp ;
	next if (m/^\##/ || m/^\s*$/) ;
	if (m/^\#CHROM/){
		print OUT "#MarkerID\tChr\tpos\ttype\tP\tM" ;
		my @a = split ;
		my @samples = @a[9..$#a] ;
		for (my $i=0; $i<@samples; $i++){
			if ($samples[$i] eq $paternal_ID){
				$p_cul = $i ;
			}
			elsif ($samples[$i] eq $maternal_ID){
				$m_cul = $i ;
			}
			else{
				$hcul2sam{$i} = $samples[$i] ;
				print OUT "\t$samples[$i]" ;
				push @culumns, $i ;
			}
		}
		print OUT "\n" ;
		next ;
	}
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genotypes) = split ;
	$total_marker ++ ;
	next if (length($ref) != 1 || length($alt) != 1) ;
	# select homozygosis SNP and different between parents
	my ($flag, $p_type, $p_depth, $m_type, $m_depth) = &judge_parents_genotype($genotypes[$p_cul], $genotypes[$m_cul]);
	next if ($flag == 0);
	next if ($p_depth eq '.' || $m_depth eq '.') ;
	next if ($p_depth < $paternal_depth || $m_depth < $maternal_depth) ;
	$marker_num++ ;
	# out put marker file
	print OUT "Marker$marker_num\t$chr\t$pos\taaxbb\ta,$p_depth\tb,$m_depth" ;
	for my $cul (@culumns){
		my ($flag, $s_type1, $s_depth1, $s_type2, $s_depth2) = &judge_sample_genotype($genotypes[$cul], $p_type, $m_type);
		if ($flag == 0){
			print OUT "\t-" ;
		}
		elsif ($flag == 1){
			print OUT "\t$s_type1,$s_depth1" ;
		}
		elsif ($flag == 2){
			print OUT "\t$s_type1,$s_depth1,$s_type2,$s_depth2" ;
		}
	}
	print OUT "\n" ;
}
close(IN);
close(OUT);
open (ST, ">$outfile.stat") || die "Can't creat $outfile.stat, $!\n" ;
print ST "Total_marker:\t$total_marker\n" ;
print ST "Remain_marker:\t$marker_num\n" ;
close(ST);
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir [$in]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	my $start_time = &show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		my $end_time = &show_log("Error: command fail: $cmd");
		&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
		exit(1);
	}
	my $end_time = &show_log("done.");
	&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
	return ;
}

## qsub
sub qsub()
{
	my ($shfile, $maxproc) = @_ ;
	if (`hostname` =~ /cluster/){
		my $cmd = "qsub-sge.pl --maxproc $maxproc --reqsub $shfile --independent" ;
		&run_or_die($cmd);
	}
	else{
		my $cmd = "ssh cluster -Y qsub-sge.pl --maxproc $maxproc --reqsub $shfile --independent" ;
		&run_or_die($cmd);
	}

	return ;
}

## qsub_mult($shfile, $max_proc, $job_num)
sub qsub_mult()
{
	my ($shfile, $max_proc, $job_num) = @_ ;
	if ($job_num > 500){
		my @shfiles = &cut_shfile($shfile);
		for my $file (@shfiles){
			&qsub($file, $max_proc);
		}
	}
	else{
		&qsub($shfile, $max_proc) ;
	}
}

#my @shfiles = &cut_shfile($shfile);
sub cut_shfile()
{
	my ($file) = @_ ;
	my @files = ();
	my $num = 0 ;
	open (IN, $file) || die "Can't open $file, $!\n" ;
	(my $outprefix = $file) =~ s/.sh$// ;
	while (<IN>){
		chomp ;
		if ($num % 500 == 0){
			close(OUT);
			my $outfile = "$outprefix.sub_".(int($num/500)+1).".sh" ;
			push @files, $outfile ;
			open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
		}
		print OUT $_, "\n" ;
		$num ++ ;
	}
	close(IN);
	close(OUT);

	return @files ;
}

#my ($flag, $p_type, $p_depth, $m_type, $m_depth) = &judge_parents_genotype($genotypes[$hparents_cul{P}], $genotypes[$hparents_cul{M}]);
sub judge_parents_genotype()
{
	my ($p_genotype, $m_genotype) = @_ ;
	my ($flag, $p_type, $m_type) ;

	my ($p_geno, $p_depth) = (split /\:/, $p_genotype)[0,2] ;
	my ($m_geno, $m_depth) = (split /\:/, $m_genotype)[0,2] ;

	if ($p_geno eq '0/0' && $m_geno eq '1/1'){
		$flag = 1 ;
		$p_type = 0 ;
		$m_type = 1 ;
	}
	elsif ($p_geno eq '1/1' && $m_geno eq '0/0'){
		$flag = 1 ;
		$p_type = 1 ;
		$m_type = 0 ;
	}
	else{
		$flag = 0 ;
	}

	return ($flag, $p_type, $p_depth, $m_type, $m_depth) ;
}

#my ($flag, $s_type1, $s_depth1, $s_type2, $s_depth2) = &judge_sample_genotype($genotypes[$cul], $p_type, $m_type);
sub judge_sample_genotype()
{
	my ($sample_genotype, $p_type, $m_type) = @_ ;
	my ($flag, $s_type1, $s_depth1, $s_type2, $s_depth2) ;

	my ($s_genotype, $allele_depth, $total_depth) = (split /\:/, $sample_genotype)[0,1,2] ;
	if ($s_genotype eq '0/0'){
		$flag = 1 ;
		$s_depth1 = $total_depth ;
		if ($total_depth eq '.'){
			$s_depth1 = 0 ;
		}
		if ($p_type == 0){
			$s_type1 = 'a' ;
		}
		else{
			$s_type1 = 'b' ;
		}
	}
	elsif ($s_genotype eq '1/1'){
		$flag = 1 ;
		$s_depth1 = $total_depth ;
		if ($total_depth eq '.'){
			$s_depth1 = 0 ;
		}
		if ($p_type == 0){
			$s_type1 = 'b' ;
		}
		else{
			$s_type1 = 'a' ;
		}
	}
	elsif ($s_genotype eq '0/1'){
		$flag = 2 ;
		$s_type1 = 'a' ;
		$s_type2 = 'b' ;
		my ($alle1_depth, $alle2_depth) = split /\,/, $allele_depth ;
		if ($p_type == 0){
			$s_depth1 = $alle1_depth ;
			$s_depth2 = $alle2_depth ;
		}
		else{
			$s_depth1 = $alle2_depth ;
			$s_depth2 = $alle1_depth ;
		}
	}
	else {
		$flag = 0 ;
	}
	return ($flag, $s_type1, $s_depth1, $s_type2, $s_depth2) ;
}

