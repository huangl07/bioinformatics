#!/usr/bin/perl -w
use FindBin qw($Bin $Script);
BEGIN{ #add absolute path of required moduals to array \@INC
	push @INC,"$Bin/draw/blib/";
	push @INC,"$Bin/Graph-0.94/lib/";
}
use Statistics::Regression;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";



my $weight = 2;
my $mapFunction = 'Kosambi';
my %mapFunction=(
	"Kosambi"=>\&Kosambi,
	"Haldane"=>\&Haldane,
	);



####################################record  .. RECORD a novel method for ordering loci on a genetic linkage map.
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$fpwd);
GetOptions(
				"help|?" =>\&USAGE,
				"loc:s"=>\$fIn,
				"pwd:s"=>\$fpwd,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);

$fIn = &ABSOLUTE_DIR($fIn);

open (LOG,">$fOut.record.log") or die $!;

my $cur_time;

######################
##########get data
#####################

my %num;my $genotype;my $rate;my %bin;
my $recom;
($genotype ,$rate ,$recom )= recominition_num( $fIn , $fpwd , \%num ,\%bin);



#####################
########order 
####################

$cur_time = &GetTime;
print LOG "\nORDER:\n\n$cur_time\n";
my @order ;
my $Count = &order(\%num,\@order );

$cur_time = &GetTime;
print LOG $cur_time,"\n";

######################
#######ripple
#######################
$cur_time = &GetTime;

print LOG "\nRipple:\n $cur_time\n";
print LOG join ("\n" , @order);

my $ripple_num = ripple(\@order,\%num);



####################
##########add bins
####################


my %temp;
#print join ("\n",@order) ;
foreach my $marker (keys %bin) {
	push @{$temp{$bin{$marker}}} , $marker ;
}

my @order_temp;
for (my $i= 0 ; $i < @order ;$i ++) {
	push @order_temp , sort {$a cmp $b} @{$temp{$order[$i]}};
}
@order = @order_temp ;
#print  "\n\n",join ("\n",@order) ;


######print OUT 

print LOG "Calculate map distance \n";

$cur_time = &GetTime;
print LOG $cur_time,"\n";


###
my $head;
open (LOC , $fIn) or die $!;
	while (<LOC>) {
		chomp;
		if (/\s=\s/) {
			$head .= "$_\n";
		}
	}
close (LOC) ;


open (OUT,">$fOut.order.loc") or die $!;
open (Map,">$fOut.map") or die $!;

print OUT $head,"\n";
map {print OUT $_ , "\t" , join ("\t",@{$genotype->{$_}}),"\n"} @order ;
print Map "group 0\n" ;
#print Map $order[0] ,"\t" ,"0\n";

my $mapdis ;

$mapdis =  weighted_linear_regression( $recom ,\@order) ;

my @dis ;
$dis[0] = 0 ;
for (my $ i = 0 ; $i < @{$mapdis} ; $i++) {
	$dis[$i+1] += abs($_) foreach (@{$mapdis}[0..$i]);
}

for (my $i = 0 ; $i< @order ;$i++) {

		print Map $order[$i],"\t", sprintf("%.3f", $dis[$i]),"\n"  ;
}


close (OUT) ;
close (Map) ;

close (LOG) ;


#################
##############subs
#################



sub recominition_num {#\loc \pwd \%rate \%bin
	my ($loc,$pwd,$num,$bin) = @_ ;
	my %bin_temp ;
	my %fil ;my %pos;my $indx_num;
	print LOG "Check bin:\n" ;
	open (LOC, $loc) or die $!;
	while (<LOC>) {
		next if (/=/ || /^\s*$/) ;
		my ($marker,@genotype) = split /\s+/ , $_ ;
		my $bin_test = join ("" , @genotype ) ;
		if (defined $bin_temp{$bin_test}) {
			print LOG $bin_temp{$bin_test} ,"\t" ,$marker,"\n\n";
			$bin->{$marker} = $bin_temp{$bin_test};
		}else{
			$bin->{$marker} = $marker ;
			$bin_temp{$bin_test} = $marker ;
		}
		
		$fil{$marker} = \@genotype ;
		$indx_num = @genotype;
		for (my $i=0; $i<@genotype ;$i++) {
			$pos{$marker}{$i} = 1 if ($genotype[$i] eq "-") ;
		}
	}
	close (LOC) ;

	my %recom ;

	my %rate;
	open (PWD,$pwd) or die $!;
	while (<PWD>) {
		my $count = 0 ;
		next if (/^;/ || /^\s*$/ || /=/) ;
		my ($marker1,$marker2,@rate) = split /\s+/ , $_ ;
		next if (@rate < 2) ;
		$recom{$marker1}{$marker2}{"r"} = $rate[0] ;
		$recom{$marker2}{$marker1}{"r"} = $rate[0] ;
		$recom{$marker1}{$marker2}{"lod"} = $rate[1] ;
		$recom{$marker2}{$marker1}{"lod"} = $rate[1] ;
		
		$marker1 = $bin->{$marker1};
		$marker2 = $bin->{$marker2};
		next if ($marker1 eq $marker2) ;
		$rate{$marker1}{$marker2} = $rate[0];
		$rate{$marker2}{$marker1} = $rate[0];
		my @pos1 = sort {$a <=> $b} keys %{$pos{$marker1}};
		my @pos2 = sort {$a <=> $b} keys %{$pos{$marker2}};
		foreach my $pos (@pos1) {
			if (defined $pos{$marker2}{$pos}) {
				$count ++;
			}
		}
		$count = $indx_num - $count ;
		$num->{$marker1}->{$marker2} = $count * $rate[0] ;
		$num->{$marker2}->{$marker1} = $count * $rate[0] ;
	}
	close (PWD) ;
	return (\%fil,\%rate ,\%recom);

}



##############


sub order {#\%num,\@order
	my ($num ,$order) = @_ ;
	my $Count;
	my @marker = sort {$a cmp $b} keys %{$num}; 
	my %marker ;
	$marker{$_} = 1 foreach (@marker) ;
	my $tt=0;
	my $ini_maker1;my $ini_maker2;


#	my ($r1 ,$r2);
#	while ($tt==0) {
#		($r1 ,$r2) = (int rand (@marker) , int rand (@marker) );
#		next if ($r1 == $r2) ;
#		$ini_maker1 = $marker[$r1] ; $ini_maker2 = $marker [$r2] ;
#		if ($num -> {$ini_maker1}->{$ini_maker2} < 0.45 * @marker ) {
#			$tt = 1;
#		}
#	}


	my $temp = 10000;
	foreach my $marker1 (keys %{$num}) {
		foreach my $marker2 (keys %{$num->{$marker1}}) {
			if ($num{$marker1}{$marker2} < $temp ) {
				$temp = $num{$marker1}{$marker2} ;
				$ini_maker1 = (sort {$a cmp $b} ($marker1 ,$marker2 ))[0]; 
				$ini_maker2 = (sort {$a cmp $b} ($marker1 ,$marker2 ))[1];
			}
		}
	}

	print LOG "initial marker : $ini_maker1\t$ini_maker2 \n " ;
	
	$order -> [0] = $ini_maker1;
	$order -> [1] = $ini_maker2;

	delete $marker{$ini_maker1};
	delete $marker{$ini_maker2};
	
	foreach my $marker ( sort {$a cmp $b}  keys %marker) {

		print LOG  "\n$marker \t insert :\n" ;
		$Count = &insert ($order , $num , $marker );
		print LOG join ("\t",@{$order}),"\n";
		my $test = ripple ($order,$num);
		print LOG  "ripple_num:\t",$test,"\n",join ("\t",@{$order}),"\n\n" if ( $test != 0 );
	}
	return $Count;
}

###################

sub insert {#\@order \%num $marker
	my ( $order , $num , $marker) = @_ ;
	my $Count = 10000000;my @temp_order;

	for (my $i = 0; $i <= @{$order} ;$i++) {
		my @test =();
		@test =( $marker , @{$order} ) if ($i == 0);
		@test =( @{$order} , $marker ) if ($i == @{$order}) ;

		if ($i != 0 && $i != @{$order} ){
			my $length = @{$order} - 1;
			@test = ( @{$order}[0..$i-1] , $marker , @{$order}[$i..$length] ) ;
		}

		my $count ;
		for (my $j=0; ;$j++) {
			last if (!defined $test[$j+1]) ;
			$count += $num -> {$test[$j]} -> {$test[$j+1]};
		}
		if ($count < $Count){
			@temp_order = @test ;
			$Count = $count;
		}
	}
	@{$order} = @temp_order ;

	return $Count;
}



###############




sub ripple {#
	my ($order,$count) = @_;
	my $count_num = 0;
	RIPPLE:for (my $i=2;$i<@order-1 ;$i++) {
			for (my $j=0;$j<@order-$i+1 ;$j++) {	
				last RIPPLE if ( $count_num > 2000 ) ;
				my $diff;
				if ($j == 0) {
					$diff = $count->{$order->[$j]}{$order->[$j+$i]} - $count->{$order->[$j+$i-1]}{$order->[$j+$i]} ;
					$diff = 0 if (abs($diff) < 10e-3) ;
					if ($diff < 0) {
						splice(@{$order},$j,$i,reverse @{$order}[$j..$j+$i-1]);
						$count_num ++ ; 
						goto RIPPLE;
					}
				}elsif($j == @order-$i){
					$diff = $count->{$order->[$j+$i-1]}{$order->[$j-1]} - $count->{$order->[$j-1]}{$order->[$j]} ;
					$diff = 0 if (abs($diff) < 10e-3 ) ;
					if ( $diff < 0) {
						splice(@{$order},$j,$i,reverse @{$order}[$j..$j+$i-1]);
						$count_num ++ ;
						goto RIPPLE;
					}
				}else{
					$diff = $count->{$order->[$j]}{$order->[$j+$i]} + $count->{$order->[$j-1]}{$order->[$i+$j-1]} - $count->{$order->[$j-1]}{$order->[$j]} - $count->{$order->[$i+$j-1]}{$order->[$i+$j]} ;
					$diff = 0 if (abs($diff) < 10e-3) ;
					if ( $diff < 0 ) {
						splice(@{$order},$j,$i,reverse @{$order}[$j..$j+$i-1]);
#						print LOG $count->{$order->[$j]}{$order->[$j+$i]} + $count->{$order->[$j-1]}{$order->[$i+$j-1]} - $count->{$order->[$j-1]}{$order->[$j]} - $count->{$order->[$i+$j-1]}{$order->[$i+$j]} ,"\t";
#						print LOG join("\t",@{$order}[$j..$j+$i-1]),"\n" ;
						$count_num ++ ; 
						goto RIPPLE
						
					}
				}
				
			}
	}
	return $count_num ;
}

########

sub weighted_linear_regression {
	
	my ($pairWise,$markerList)=@_;
	
	my $variablesNum=@{$markerList}-1;
	my @variablename;
	for ( my $i=0;$i<@{$markerList}-1;$i++) {
		push @variablename,$i;
		}
	my $obj=Statistics::Regression->new($variablesNum,\@variablename);
	
	for (my $i=0;$i<@{$markerList}-1;$i++) {
			
			my @coefficientVector=();
			my $const=0;
			
			for (my $j=0;$j<@{$markerList}-1;$j++) {
					
					my ($min,$max)=sort {$a <=> $b} ($i,$j);
					my $sum=0;
					
					for(my $h=0;$h<=$min;$h++){
						for (my $k=$max+1;$k<@{$markerList};$k++) {
							
							next if (!defined $pairWise->{$markerList->[$h]}{$markerList->[$k]}) ;
							$sum+=$pairWise->{$markerList->[$h]}{$markerList->[$k]}{"lod"}**$weight;
							
						}
					}
					
					push @coefficientVector,$sum;
					
					if ($j<=$i) {
						
						for (my $k=$i+1;$k<@{$markerList};$k++) {
							
							next if (!defined $pairWise->{$markerList->[$j]}{$markerList->[$k]}) ;
							$const+=$pairWise->{$markerList->[$j]}{$markerList->[$k]}{"lod"}**$weight*(&{$mapFunction{$mapFunction}}($pairWise->{$markerList->[$j]}{$markerList->[$k]}{"r"}));
							
						}

					}
			}

#			print join("\t",@coefficientVector,$const),"\n";

			$obj->include($const,\@coefficientVector);
	}

#	print "==================================================\n";
	
	my @result=$obj->theta;
	return \@result;
}




########
sub Haldane {#
	my $r=shift;
	if ($r>=0.5) {
		$r = 0.4999;
	}
	my $result=-100*(1/2)*log(1-2*$r);
	return $result;
}

sub Kosambi {#
	my $r=shift;
	if ($r>=0.5) {
		$r = 0.4999;
	}
	my $result=100*(1/4)*log((1+2*$r)/(1-2*$r));
	return $result;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Readfasta {
	my $infile=shift;
	my $seqId_seq=shift;

	my $c=0;
	open IN, $infile || die $!;
	my $seqId;
	while (<IN>) {
		if (/^>(\S+)/) {
			$seqId=$1;
			$c++;
		}
		else {
			$_=~s/\s//g;
			$seqId_seq->{$seqId}.=$_;
		}
	}
	close IN;
	return $c;
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
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
Program: record  ..
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage:
  Options
  -help		USAGE
  -loc		loc file
  -pwd		pwd file
  -o		outfiles . fOut.map  fOut.order.loc
USAGE
	print $usage;
	exit;
}





