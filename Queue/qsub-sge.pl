#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use MIME::Lite;

my $version="1.0.0";
my ($Queue,$CPU,$Nodes,$maxjob,$Resource,$runtime);
GetOptions(
	"Queue:s"=>\$Queue,
	"CPU:s"=>\$CPU,
	"Nodes:s"=>\$Nodes,
	"Maxjob:s"=>\$maxjob,
	"Resource:s"=>\$Resource,
	"Runtime:s"=>\$runtime,
) or &USAGE;
$Queue ||= "dna";
$maxjob ||= 20;
$Resource ||= "mem=3G";
$runtime||="walltime=240:00:00";
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
	$qsub="qsub -l nodes=$Nodes:ppn=$CPU -l $Resource -l $runtime -q $Queue";
}else{
	$qsub="ssh -Y login02 && qsub -l nodes=$Nodes:ppn=$CPU -l $Resource -l $runtime -q $Queue";
}
print Log $qsub;
my @Shell;
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
	$job_mark++;
}
close In;
chdir($workdir);
my %job;
my %queue,
my $nqsubed=0;
my %Shell;
while ($nqsubed ne scalar @Shell ) {
	my $start=$nqsubed;
	my $end=$nqsubed+$maxjob-scalar keys %job;
	$end=scalar @Shell if ($end > scalar @Shell);
	for (my $i=$start;$i<$end;$i++) {
		my $job=$Shell[$i];
		my $jobs="$qsub $job";
		my $return=`$jobs `;
		while ($return!~/^\d+/) {
			sleep(10);
			$return=`$jobs`;
		}
		$return=(split(/\./,$return))[0];
		open $filehand{$return},">$job.m$return"; 
		print {$filehand{$return}} "#jobid\tmemused\tvmemused\truntime\tcputime\thostname\n";
		$Shell{$job}{e}="$job.e$return";
		$Shell{$job}{o}="$job.o$return";
		$job{$return}=$job;
		$queue{$return}=time();
	}
	$nqsubed=$end;
	my $check=Checkjob(\%job,\%queue,$qsub,\%filehand);
	while (scalar keys %job == $maxjob) {
		sleep(300);
		$check=Checkjob(\%job,\%queue,$qsub,\%filehand);
	}
}
print Log "qsub done! now check!";
while (scalar keys %job != 0) {
	my $check=Checkjob(\%job,\%queue,$qsub,\%filehand);
	print Log "there were ",scalar keys %job," still running! Please wait!\n";	
	sleep(60);
}
print Log scalar @Shell,"still not done! poor luck !\n";
print Log "$worksh \n qsub done£¡\nTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
close Log;
open Out,">$worksh.reqsub.sh";
foreach my $Shell (@Shell) {
	if (!-f "$Shell.check" || !-f $Shell{$Shell}{e} || !-f $Shell{$Shell}{o}) {
		open In,$Shell;
		while (<In>) {
			chomp;
			next if ($_ eq ""||/^$/);
			my @line=split(/\&\&/,$_);
			print Out join("&&",@line[0..$#line-2]),"\n";
		}
		close In;
	}
}
close Out;
my $nline=`wc -l $worksh.reqsub.sh`;
$nline=chomp($nline);
my $who=`whoami`;
chomp($who);
if ($nline > 0) {
	SendAttachMail("poor luck!you must check $worksh.reqsub.sh",'',"$who\@majorbio.com");
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub SendAttachMail
{
    my $mail_content = shift;
    my $mail_attach = shift;
	my $mail_to=shift;
	my $path=dirname($mail_attach);
	my $file=basename($mail_attach);
	my $mail_from="long.huang\@majorbio.com";
	my $mail_cc="long.huang\@majorbio.com";
	my $mail_subject="poor luck,qsub error";
    my $msg=MIME::Lite->new(
                    From=>$mail_from,
                    To=>$mail_to,
                    Cc=>$mail_cc,
                    Subject=>$mail_subject,
                    Type=>'TEXT',
                    Data=>$mail_content,);
	if ($mail_attach != "") {
		    $msg->attach(
            Type=>'AUTO',
            Path=>$path,
            Filename=>$file,);
	}
	$msg->send('smtp', "smtp.majorbio.com", AuthUser=>"dna\@majorbio.com", AuthPass=>"majorbio-genome-2017" );
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
sub Checkjob{
	my ($qsub,$queue,$reqsub,$filehand)=@_;
	my $job=0;
	my $whoami=`whoami`;
	my $check=`qstat -fu $whoami`;
	my @line=split(/\n/,$check);
	my $jobid="";
	my %checked;
	my $mem=0;
	my $vmem=0;
	my $host="--";
	my $time=0;
	my $ctime=0;
	my $job_stat="R";
	foreach my $l (@line) {
		if ($l=~/Job Id/) {
			if ($jobid ne "") {
				print {$$filehand{$jobid}} "$jobid\t$mem\t$vmem\t".$time."\t",$ctime,"\t",$host,"\n"; 
			}
			my $id=$l;
			$id=~s/\D//g;
			if (exists $$qsub{$id}) {
				$jobid=$id;
				$checked{$id}=1;
				$job++;
			}else{
				$jobid="";
			}
		}elsif ($jobid ne "") {
			my $id=$jobid;
			$mem=(split(/\s+/,$l))[-1] if ($l =~ /resources_used.mem/);
			$vmem=(split(/\s+/,$l))[-1] if ($l=~/resources_used.vmem/);
			$host=(split(/\s+/,$l))[-1] if ($l=~/exec_host/);
			$time=(split(/\s+/,$l))[-1] if ($l=~/resources_used.walltime/);
			$ctime=(split(/\s+/,$l))[-1] if ($l=~/resources_used.cput/);
			$job_stat=(split(/\s+/,$l))[-1] if ($l =~ /job_state/);
			if ($job_stat eq "Q") {
				if (!exists $$queue{$id}) {
					$$queue{$id}=time();
				}else{
					$time=time()-$$queue{$id};
					if ($time > 7200) {
						`qdel $id`;
						my $jobsh=$$qsub{$id};
						my $jobs="$reqsub $jobsh";
						my $return=`$jobs `;
						while ($return!~/^\d+/) {
							sleep(5);
							$return=`$jobs`;
						}
						$return=(split(/\./,$return))[0];
						delete $$qsub{$id};
						delete $$queue{$id};
						close $$filehand{$id};
						open $$filehand{$return},">$jobsh.m$return";
						print {$$filehand{$return}} "#jobid\tmemused\tvmemused\truntime\tcputime\thostname\n";
						$$qsub{$return}=$jobsh;
						$$queue{$return}=time();
					}
				}
			}
		}else{
			next;
		}
	}
	foreach my $id (sort keys %$qsub) {
		if (!exists $checked{$id}) {
			delete $$qsub{$id};
			delete $$queue{$id}
		}
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
	--Runtime	<str> set the maximum time for one job running default "walltime=240:00:00"
	--h         Help

USAGE
        print $usage;
        exit;
}
