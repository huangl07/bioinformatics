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
GetOptions(\%opts,"g=s","od=s","p=s","max=s","f=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{g}) || !defined($opts{od}) || !defined($opts{p}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		This script is write to separate genotype file into chromosome genotype files
		Version: $ver

	Usage:
		-g           genotype file               <infile>     must be given
		-od          result file dir             <outdir>     must be given
		-p           prefix of outfile           <string>     must be given
		-max         max outfile number          [int]        optional [500]
		-f           whether swap the type       [int]        optional [1] [1: yes, 0: no]
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
my $gtpfile = $opts{g} ;
my $outdir = $opts{od} ;
mkdir $outdir if (!-d $outdir);
my $prefix = $opts{p} ;
my $max_file_num = defined $opts{max} ? $opts{max} : 500 ;
my $swap_flag = defined $opts{f} ? $opts{f} : 1 ;

## reading genotype file and separate to chr
&separate_genotype_file($gtpfile, $outdir, $prefix);


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

#&separate_genotype_file($gtpfile, $outdir, $prefix);
sub separate_genotype_file()
{
	my ($gtpfile, $outdir, $prefix) = @_ ;

	open (IN, $gtpfile) || die "Can't open $gtpfile, $!\n" ;
	my %hfile = ();
	my $file_num = 0 ;
	my @head ;
	while (<IN>){
		chomp ;
		next if (m/^\s*$/);
		my ($markerID, $chr, $pos, $type, $P, $M, @sample_types) = split ;
		#next if ($markerID =~ /MarkerID/i || $markerID =~ m/BlockID/i);
		if ($markerID =~ /MarkerID/i || $markerID =~ m/BlockID/i){
			@head = ($markerID, $chr, $pos, $type, $P, $M, @sample_types) ;
			next ;
		}
		if (!defined $hfile{$chr}){
			$file_num++ ;
			next if ($file_num > $max_file_num);
			my $outfile = "$outdir/$prefix.$chr.genotype" ;
			open ($hfile{$chr}, ">$outfile") || die "Can't creat $outfile, $!\n" ;
			print {$hfile{$chr}} "$head[0]\t$head[3]\t", join("\t", @head[6..$#head]),"\n" ;
		}
		if (defined $hfile{$chr}){
			if ($swap_flag == 1){
				for (my $i=0; $i<@sample_types; $i++){
					$sample_types[$i] = &swap_type($sample_types[$i]);
				}
			}
			print {$hfile{$chr}} $markerID, "\t", $type, "\t", join("\t", @sample_types), "\n" ;
		}
	}
	close(IN);

	for my $chr (keys %hfile){
		close($hfile{$chr}) ;
	}

	return ;
}

#my $type = &swap_type($sample_types[$i]);
sub swap_type()
{
	my ($raw_type) = @_ ;

	if ($raw_type =~ /-/){
		return ('--');
	}
	elsif ($raw_type =~ /\,/){
		my ($t1, $d1, $t2, $d2) = split /\,/, $raw_type ;
		if (!defined $t2){
			return ($t1.$t1);
		}
		else{
			return ($t1.$t2);
		}
	}
	else{
		return ($raw_type);
	}
}

