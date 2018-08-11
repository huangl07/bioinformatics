#!/usr/bin/env perl
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $force;
GetOptions(
	"force"=>\$force,
) or &USAGE;

my $who=`whoami`;
chomp $who;
my $qstat=`qstat -fu $who|sed 's/^\\s*//g'`;
my @line=split("=",$qstat);
for (my $j=0;$j<@line;$j++){
	my @line2=split(/\n+/,$line[$j]);
	if (scalar @line2 > 1){
		$line[$j]=join("",@line2[0..$#line2-1])."\n".$line2[-1];
	}
}
my $jobs=join("=",@line);
my @jobs=split(/\n/,$jobs);
my $shell=shift;
&USAGE if (!defined $shell||$shell eq "");
if ($shell eq "all"){
	print "all jobs from $who will all be delete! please check\n\n";
	open Out,">delete.sh";
	foreach my $job(@jobs){
		if ($job=~/Job Id: (\d+)\.majorbio\.cluster\.com/){
			print Out "delete $1\n";
		}
	}
	close Out;
}else{
	my $jobid;
	$shell=PATH($shell);
	foreach my $job(@jobs){
		if ($job=~/Job Id: (\d+)\.majorbio\.cluster\.com/){
			$jobid=$1;	
		}
		if($job=~/Error_Path = login02\.local:([^\s+]*)/ ){
			my $path=$1;
			if ($path=~/$shell\.\d+\.qsub\/sub\d+/){
				print Out "delete $jobid\n";
			}
		}
		
	}
}
print "run your nead :\n sh delete.sh";
if ($force){
	`sh delete.sh`;
}
sub PATH{
        my ($in)=@_;
        my $return="";
		my $cur_dir=`pwd`;chomp$cur_dir;
        if(-f $in)
        {
                my $dir=dirname($in);
                my $file=basename($in);
                chdir $dir;$dir=`pwd`;chomp $dir;
                $return="$dir/$file";
        }
        elsif(-d $in)
        {
                chdir $in;$return=`pwd`;chomp $return;
        }
        else
        {
                warn "Warning just for file and dir\n";
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
	This program will generate a delete.sh for qdel some jobs
	if work.sh is all it will generate a shell by usrid [carefully]
	if work.sh provide it will generate a delete.sh for qdel jobs which qsub-sge.pl throwed by the work.sh
	eg:
	perl $Script [work.sh]
Usage:
  Options:
	--force	<str>	delete all jobs not only generate a delete sh
	--h         Help

USAGE
        print $usage;
        exit;
}
