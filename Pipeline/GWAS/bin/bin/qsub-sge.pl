#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
my ($Queue,$CPU,$Nodes,$maxjob,$Resource);
GetOptions(
	"Queue:s"=>\$Queue,
	"CPU:s"=>\$CPU,
	"Nodes:s"=>\$Nodes,
	"Maxjob:s"=>\$maxjob,
	"Resource:s"=>\$Resource,
) or &USAGE;
$Queue ||= "dna";
$maxjob ||= 20;
$Resource ||= "mem=3G";
$CPU ||=1;
$Nodes||=1;
my $hostname=`hostname`;
if ($hostname !~/local/) {
	print  "error cluster! must head cluster or local cluster!";
}
my $worksh = shift;
&USAGE if (!defined $worksh||$worksh eq "");
$worksh=PATH($worksh);
my $id=$$;
my $workdir=$worksh.".".$id."."."qsub";
mkdir $workdir if (!-d $workdir);
my $basename=basename($worksh);
die "error input file $worksh" if (!-f $worksh); 
$Resource=~/\=(.+)G/;
$Resource="mem=$1G";
open In,$worksh||die "error open $worksh!";
my %worksh;
my $job_mark="00001";
open Log,">$worksh.$id.log";
my $qsub;
if ($hostname =~/local/) {
	$qsub="qsub -l nodes=$Nodes:ppn=$CPU -l $Resource -q $Queue";
}else{
	$qsub="ssh -Y login02 && qsub -l nodes=$Nodes:ppn=$CPU -l $Resource -q $Queue";
}
print Log $qsub;
my @Shell;
#my @SHfile;
my %filehand;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my $id=$basename."_".$job_mark;
	my @line=split(/\s+/,$_);
	my $l=join(" ",@line);
	$l=~s/&&$//g;
	open Out,">$workdir/$id.sh";
	print Out $l." && echo This-work-is-complete! && touch $workdir/$id.sh.check\n";
	close Out;
	my $job="$workdir/$id.sh";
	push @Shell,$job;
	#push @SHfile,"$workdir/$id.sh";
	$job_mark++;
}
close In;
chdir($workdir);
my $Cycle=0;
my %job;
my %jobid;
my $nshell=scalar @Shell;
while (scalar @Shell!=0) {
	my $check=Checkjob(\%job);
	foreach my $job (@Shell) {
		while (scalar keys %job == $maxjob) {
			sleep(60);
			$check=Checkjob(\%job);
		}
		if (scalar keys %job < $maxjob) {
			my $jobs="$qsub $job";
			my $return=`$jobs `;
			while ($return!~/^\d+/) {
				sleep(5);
				$return=`$jobs`;
			}
			$return=(split(/\./,$return))[0];
			open $filehand{$return},">$job.m$return"; 
			print {$filehand{$return}} "#jobid\tmemused\tvmemused\truntime\tcputime\thostname\n";
			$job{$return}=$job;
			$jobid{$job}=$return;
		}
	}
	print Log "qsub done! now check!";
	while (scalar keys %job != 0) {
		$check=Checkjob(\%job);
		print Log "cycle:$Cycle:there were ",scalar keys %job," still running! Please wait!\n";	
		sleep(60);
	}
	my @newshell;
	foreach my $job (@Shell) {
		if (!-f "$job.check") {
			push @newshell,$job;
		}
	}
	@Shell=();
	push @Shell,@newshell;
	%job=();
	print Log "cycle:$Cycle: there were ",scalar @Shell," not done! reqsub\n";
	$Cycle++;
	last if ($Cycle > 5 || $nshell == scalar @Shell);
}

print Log scalar @Shell,"still not done! poor luck !\n";
print Log "$worksh \n qsub done£¡\nTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
close Log;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
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
sub Checkjob{
	my ($qsub)=@_;
	my $job=0;
	foreach my $id (sort keys %$qsub) {
		my $check=system("qstat -f $id 1>$workdir/check.log 2>&1");
		open LOG,"$workdir/check.log";
		my @line=<LOG>;
		close LOG;
		my $line=join("\n",@line);
		if ($line=~/Unknown Job Id/) {
			delete $$qsub{$id};
			next;
		}
		$job++;
		my $mem=0;
		my $vmem=0;
		my $host="--";
		my $time=0;
		my $ctime=0;
		foreach my $line (@line) {
			$mem=(split(/\s+/,$line))[-1] if ($line =~ /resources_used.mem/);
			$vmem=(split(/\s+/,$line))[-1] if ($line=~/resources_used.vmem/);
			$host=(split(/\s+/,$line))[-1] if ($line=~/exec_host/);
			$time=(split(/\s+/,$line))[-1] if ($line=~/resources_used.walltime/);
			$ctime=(split(/\s+/,$line))[-1] if ($line=~/resources_used.cput/);
		}
		print {$filehand{$id}} "$id\t$mem\t$vmem\t".$time."\t",$ctime,"\t",$host,"\n"; 
		sleep(2);
	}
	return $job;
}
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description: 
	This program throw the jobs and control them running on linux SGE system. 
	It reads jobs from an input shell file. One line is the smallest unit of a single job.
	It will create a path at work.sh's dir and do qsub by one line.
	It will check all jobs doing and stderr output is empty.
	eg:
	perl $Script --Queue dna --Resource mem=3G --CPU 8  --Nodes 1 --Maxjob 20 ./work.sh
Usage:
  Options:
	--Queue	<str>	set the queue to use,default "dna"
	--Resource <str>	set the memeory used for one job,default "mem=10G"
	--CPU <num>	set the cpu number to used for one job,default 8
	--Nodes <num>	set the procesor to submit for one jobs,default 1      
	--Maxjob <num>   set the maximum number of jobs in a qsub line, default 20
	--h         Help

USAGE
        print $usage;
        exit;
}
