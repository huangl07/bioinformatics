use strict;
use POSIX qw(log10);

open IN,"pop.bf";
while(<IN>){
    $_=~s/[\n\r]//g;
    my (undef,$test1,$test2,@others)=split/\t/,$_;
    my $log_10 = log10($test2);
    if($log_10>1.5){
        print $log_10,"\n";
        }
    }
