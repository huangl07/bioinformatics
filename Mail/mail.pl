#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($maillist,$Subject,$Attach,$Content);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"Mail:s"=>\$maillist,
	"Subject:s"=>\$Subject,
	"Attach:s"=>\$Attach,
	"Content:s"=>\$Content,
			) or &USAGE;
&USAGE unless ($maillist and $Subject and $Content);
 my $mail_content = `less -S $Content`;
 my $mail_from="long.huang\@majorbio.com";
 my $mail_subject=$Subject;
 my $mail_attach=$Attach;
 my $path=dirname($Attach);
 my $file=basename($Attach);
open In,$maillist;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my $mail_to=$_;
    my $msg=MIME::Lite->new(
                    From=>$mail_from,
                    To=>$mail_to,
                    Subject=>$mail_subject,
                    Type=>'TEXT',
                    Data=>$mail_content,);
	if ($mail_attach != "") {
		    $msg->attach(
            Type=>'AUTO',
            Path=>$path,
            Filename=>$file,);
	}
	$msg->send('smtp', "smtp.majorbio.com", AuthUser=>"long.huang\@majorbio.com", AuthPass=>"cwn.711hl" );

}
close In;
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
