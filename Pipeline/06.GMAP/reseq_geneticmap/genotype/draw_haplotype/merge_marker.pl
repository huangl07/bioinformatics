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
GetOptions(\%opts,"i=s","o=s","n=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || !defined($opts{n}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:
		-i           marker of block file          <infile>     must be given
		-o           result file                   <outfile>    must be given
		-n           marker number to merge        <int>        must be given
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
my $marker_num = $opts{n} ;

## reading marker file and merge
open (IN,"$infile") || die "Can't open $infile,$!\n" ;
open (OUT,">$outfile") || die "Can't creat $outfile,$!\n" ;
my $num = 0 ;
my @marker_info = ();
while (<IN>){
	chomp ;
	next if (m/^\#/ || m/^\s*$/);
	my ($markerID, $type, @genotypes) = split ;
	next if ($markerID =~ m/ID$/i) ;
	push @marker_info, [$markerID, $type, @genotypes] ;
	$num ++ ;
	if ($num % $marker_num == 0){
		my ($markerID, $type, @genotypes) = &merge_marker(\@marker_info);
		print OUT "$markerID\t$type\t", join("\t",@genotypes), "\n" ;
		@marker_info = ();
	}
}
if ($num % $marker_num != 0){
	my ($markerID, $type, @genotypes) = &merge_marker(\@marker_info);
	print OUT "$markerID\t$type\t", join("\t",@genotypes), "\n" ;
	@marker_info = ();
}
close(IN);
close(OUT);
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

#my ($markerID, $type, @genotypes) = &merge_marker(\@marker_info);
sub merge_marker()
{
	my ($amarker_info) = @_ ;
	my ($final_markerID, $final_type, @final_genotypes) ;

	my %hgenotype = ();
	for (my $i=0; $i<@{$amarker_info}; $i++){
		my ($markerID, $type, @genotypes) = @{$amarker_info->[$i]} ;
		$final_markerID = $markerID if (!defined $final_markerID);
		$final_type = $type if (!defined $final_type) ;
		for (my $j=0; $j<@genotypes; $j++){
			$hgenotype{$j}{$genotypes[$j]}++ ;
		}
	}

	for my $i (sort {$a<=>$b} keys %hgenotype){
		my ($genotype) = sort {$hgenotype{$i}{$b} <=> $hgenotype{$i}{$a}} keys %{$hgenotype{$i}} ;
		$final_genotypes[$i] = $genotype ;
	}

	return ($final_markerID, $final_type, @final_genotypes);
}


