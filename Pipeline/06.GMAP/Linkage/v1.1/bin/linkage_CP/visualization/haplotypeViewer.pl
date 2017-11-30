#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#######写程序、修改程序一定要注明时间，人物，目的，对程序添加一些注释、描述信息########
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fMap,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"m:s"=>\$fMap,
				) or &USAGE;
&USAGE unless ($fIn  and $fMap and $fOut);


#2012-1-16 更改判断单体来源函数的骨架，同时纠正不能判断单一亲本连锁相缺失情况下的ef基因型来源的bug

my @Marker;

	#read map file
	open(MAP,$fMap);
	while (<MAP>) {
		chomp;
		if (/^M/) {
			my($marker,undef)=split /\s+/,$_;
			push @Marker,$marker;
		}
	}
	close(MAP);

	#loc2loc
	&Loc2loc($fIn);

	#judgment source
	&JudgmentScurce($fIn);

	#repair 
	&Repair($fIn);

	my %Repair;
	my $marker1;my $marker2;my @tqNumber;my @reNumber;
	my $tq_num=0;
	open(TQ1,"$fOut.tq");
	my $head1=<TQ1>;
	while (<TQ1>) {
		chomp;
		($marker1,@tqNumber)=split /\s+/,$_;
		@{$Repair{$marker1}{"TQ"}}=@tqNumber;
	}
	close(TQ1);

	open(TQRE,">$fOut.repairDis.tq");# add at 2011-11-17 by wus
	print TQRE $head1;
	open(RE,"$fOut.repair.tq");
	my $head2=<RE>;
	while (<RE>) {
		chomp;
		($marker2,@reNumber)=split /\s+/,$_;
		if (exists $Repair{$marker2}) {
			print TQRE $marker2,"\t";
			for (my $i=0;$i<@reNumber ;$i+=3) {
				if ($Repair{$marker2}{"TQ"}[$i] eq 0) {
					if ($reNumber[$i] eq "-") {
						print TQRE "0\t";
					}else{
						print TQRE $reNumber[$i],"\t";
					}
				}else{
					print TQRE $reNumber[$i],"\t";
				}

				if ($Repair{$marker2}{"TQ"}[$i+1] eq 0) {
					if ($reNumber[$i+1] eq "-") {
						print TQRE "0\t";
					}else{
						print TQRE $reNumber[$i+1],"\t";
					}
				}else{
					print TQRE $reNumber[$i+1],"\t";
				}
				
				print TQRE "0\t";
			}
		}
		print TQRE "\n";
	}
	close(RE);
	close(TQRE);
	# add at 2011-11-17 by wus
	`$Bin/matrix2png -data $fOut.tq -rcd -size 8:8 -mincolor white  -midcolor green -maxcolor 0:40:120 > $fOut.original.tq.png`;
	`$Bin/matrix2png -data $fOut.repair.tq -rcd -size 8:8 -mincolor white  -midcolor green -maxcolor 0:40:120 > $fOut.repair.tq.png`;
	`$Bin/matrix2png -data $fOut.repairDis.tq -rcd -size 8:8 -mincolor white  -midcolor green -maxcolor 0:40:120 > $fOut.repairDis.tq.png`;

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
sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;

        if ($type eq 'dir')
        {
			my $pwd = `pwd`;
			chomp $pwd;
			chdir($input);
			$return = `pwd`;
			chomp $return;
			chdir($pwd);
        }
        elsif($type eq 'file')
        {
			my $pwd = `pwd`;
			chomp $pwd;

			my $dir=dirname($input);
			my $file=basename($input);
			chdir($dir);
			$return = `pwd`;
			chomp $return;
			$return .="\/".$file;
			chdir($pwd);
        }
        return $return;
}

sub Repair{
	my ($tqfile)=@_;
	my @Markers;
	my $indiNum;
	my %chr;
	my %chrRepair;

	open (IN,"$fOut.tq") or die $!;
	my $head=<IN>;
	while (<IN>) {
		chomp;
		my ($marker,@number)=split /\s+/,$_;
		$indiNum=@number;
		push @Markers,$marker;
		for (my $i=0;$i<@number;$i++) {
			push @{$chr{$i}},$number[$i];
		}
	}
	close (IN) ;

	my @arr1=@{$chr{0}};
	my @arr=&cutArr(@arr1);

	for (my $i=0;$i<$indiNum;$i++) {
		@{$chrRepair{$i}}=&cutArr(@{$chr{$i}});
	}

	open (OUT,">","$fOut.repair.tq") or die $!;
	print OUT $head;
	my $pos=0;
	foreach my $marker (@Markers) {
		print OUT $marker,"\t";
		for (my $i=0;$i<$indiNum;$i++) {
			print OUT ${$chrRepair{$i}}[$pos],"\t";
		}
		print OUT "\n";
		$pos++;
	}
	close(OUT);
}

sub cutArr {#
	my @array=@_;
	my @ordarray=@array;
	for (my $i=0;$i<@array;$i++) {
		if ($array[$i] eq 0) {
			$array[$i]="-"; 
		}
	}
	my $left=$array[0];
	my $part=1;
	my %arrPart;
	my %flag;
	my @repair;
	for (my $i=0;$i<@array;$i++) {
		if ($array[$i] eq $left) {
			push @{$arrPart{$part}},$array[$i];
			$flag{$part}=$array[$i];
		}
		else {
			$part++;
			$left=$array[$i];
			push @{$arrPart{$part}},$array[$i];
			$flag{$part}=$array[$i];
		}
	}
	if ($part==1) {
		if (${$arrPart{1}}[0] eq "-") {
			return @ordarray;
			next;
		}
		return @{$arrPart{1}};
		next;
	}
	for (my $i=1;$i<=$part;$i++) {
		if ($flag{$i} eq "-") {
			$flag{$i-1}=exists $flag{$i-1} ? $flag{$i-1} : $flag{$i+1};
			$flag{$i+1}=exists $flag{$i+1} ? $flag{$i+1} : $flag{$i-1};
			if ($flag{$i-1} eq $flag{$i+1}) {
				push @repair,$flag{$i-1} foreach (@{$arrPart{$i}}) ;
			}
			else {
				push @repair,@{$arrPart{$i}};
			}
		}
		else {
			push @repair,@{$arrPart{$i}};
		}
	}
	return @repair;
}


sub Loc2loc {
	my ($fIn)=@_;
	my $header_done_flag=0;
	my %Header=();

	open (IN,$fIn) or die $!;
	$/="\n";
	my $header_str="";
	while (! $header_done_flag) {
		my $line=<IN>;
		$line=~s/\r//g;
		$header_str.=$line;
		chomp $line;
		next if ($line =~ /^$/ || $line =~ /^;/ || $line =~ /^\s+;/) ;
		last if ($line =~ /locus numbers/) ;

		if ($line =~ /name\s+\=\s+(\S+)/) {
			$Header{"name"}=$1;
		}elsif($line =~ /popt\s+\=\s+(\S+)/){
			$Header{"popt"}=$1;
		}elsif($line =~ /nloc\s+\=\s+(\S+)/){
			$Header{"nloc"}=$1;
		}elsif($line =~ /nind\s+\=\s+(\S+)/){
			$Header{"nind"}=$1;
		}

		if (exists $Header{"name"} && exists $Header{"popt"} && exists $Header{"nloc"} && exists $Header{"nind"} ) {
			$header_done_flag=1;
		}
	}

	my $Marker_str="";
	$/="\n";
	while (my $line=<IN>) {
		$line=~s/\r//g;
		chomp $line;
		next if ($line =~ /^$/ || $line =~ /^;/ ) ;
		last if ($line =~ /locus numbers/ || $line =~ /individual names/) ;
		my ($real,$anno)=split(";",$line);
		$Marker_str.=" ".$real;
	}
	close (IN) ;

	open (OUT,">$fOut.phaseloc") or die $!;
	print OUT $header_str,"\n";
	my $curr_offerset=0;
	$Marker_str=~ s/^\s+//;
	my @marker=split /\s+/,$Marker_str;
	my @marker_arr=();
	for (my $i=0;$i<@marker ;$i++) {
		if ($i>0 && $i % ($Header{"nind"}+3+$curr_offerset) == 0) {
			print OUT join("\t",@marker_arr),"\n";
			@marker_arr=();
		}
		push @marker_arr,$marker[$i];
	}

	if (@marker_arr) {
		print OUT join("\t",@marker_arr),"\n";
		@marker_arr=();
	}
	close (OUT) ;
}


#判断单体来源
#判断子代基因型来自亲本的第几个基因型，-为缺失，0为无法判断，1为亲本第一位，2为亲本第二位
sub JudgmentScurce {
	my ($fIn)=@_;
	open(OUT,">$fOut.phase");
	open(TQ,">$fOut.tq");
	my $mark=0;
	my %file;#保存格式转换好的loc文件内容
	open(IN,"$fOut.phaseloc");
	while (<IN>) {
		chomp;
		s/\r//g;
		next if(/^$/ || /""/ || /^;/ || /^name/ || /^popt/ || /^nloc/ || /^nind/);
		my ($marker,@line)=split/\s+/,$_;
		$file{$marker}{"line"}=$_;
	}
	close(IN);


	for (my $i=0;$i<@Marker ;$i++) {
		my $marker=$Marker[$i];
		my @temp;
		if (!exists $file{$marker}) {
			print STDOUT "not exists $marker in loc file \n";
			exit;
		}else{
			@temp=split /\s+/,$file{$marker}{"line"};
			if ($mark == 0) {
				print TQ "ID\t";
				for (my $i=3;$i<@temp ;$i++) {
					print TQ $i-2,"\t",$i-2,"\t",$i-2,"\t";
				}
				print TQ "\n";
				$mark++;
			}
		}
		$temp[1] = "<".$temp[1].">" if ($temp[1] !~/<\w+>/) ; ## added by macx 2012-5-24

		my @parent=split("",$temp[1]);# 保存亲本基因型 lmxll
		my @phase=split("",$temp[2]);# 保存连锁相 {0-}
		print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t";
		print TQ $temp[0],"\t";

		for (my $i=3;$i<@temp ;$i++) {
			my @Ind=split("",$temp[$i]);#循环子代基因型
			my $f1;my $f2;my $m1;my $m2;#亲本第一二位基因型
			my $off1;my $off2;#子代来自哪个亲本
			my %hash;

			#根据连锁相判断亲本第一、二位的基因型
			if ($phase[1] eq "0") {
				$f1=$parent[1];
				$f2=$parent[2];
				$hash{"F"}{$parent[1]}++;
				$hash{"F"}{$parent[2]}++;
			}elsif($phase[1] eq "1"){
				$f1=$parent[2];
				$f2=$parent[1];
				$hash{"F"}{$parent[1]}++;
				$hash{"F"}{$parent[2]}++;
			}else{
				$f1=$f2="-";
				$hash{"F"}{$parent[1]}++;
				$hash{"F"}{$parent[2]}++;
			}
			if ($phase[2] eq "0") {
				$m1=$parent[4];
				$m2=$parent[5];
				$hash{"M"}{$parent[4]}++;
				$hash{"M"}{$parent[5]}++;
			}elsif($phase[2] eq "1"){
				$m1=$parent[5];
				$m2=$parent[4];
				$hash{"M"}{$parent[4]}++;
				$hash{"M"}{$parent[5]}++;
			}else{
				$m1=$m2="-";
				$hash{"M"}{$parent[4]}++;
				$hash{"M"}{$parent[5]}++;
			}
			
			#子代基因型缺失
			if ($Ind[0] eq "-" && $Ind[1] eq "-") {
				$off1="-";
				$off2="-";
				print OUT "--\t";
				print TQ "-\t-\t0\t";
			}
			elsif($parent[1] eq $parent[4] && $parent[2] eq $parent[5] && $parent[1] ne $parent[2]){
			#基因型为hkxhk，只能判断hh、kk类型，无法区分hk类型
				if ($Ind[0] ne $Ind[1]) {#hk类型
					$off1=0;#not sure
					$off2=0;#not sure
				}else{#hh或kk类型
					if ($f1 eq "-" && $f2 eq "-") {
						$off1=0;
						if ($m1 eq "-" && $m2 eq "-") {
							$off2=0;
						}else{
							if($Ind[1] eq $m1){#第一位应来自父本，第二位应来自母本，没有亲本位置不确定情况
								$off2=1;
							}elsif($Ind[1] eq $m2){
								$off2=2;
							}else{
								$off2=0;
							}
						}
					}
					if ($m1 eq "-" && $m2 eq "-") {
						$off2=0;
						if ($f1 ne "-" && $f2 ne "-") {
							if ($Ind[0] eq $f1) {
								$off1=1;
							}elsif($Ind[0] eq $f2){
								$off1=2;
							}else{
								$off1=0;
							}
						}
					}
					if ($f1 ne "-" && $f2 ne "-" && $m1 ne "-" && $m2 ne "-") {
						if ($Ind[0] eq $f1) {#判断第一个基因型
							$off1=1;
						}elsif($Ind[0] eq $f2){
							$off1=2;
						}else{
							$off1=0;
						}
						if ($Ind[1] eq $m1) {#判断第二个基因型
							$off2=1;
						}elsif($Ind[1] eq $m2){
							$off2=2;
						}else{
							$off2=0;
						}
					}
				}
				print OUT $off1,$off2,"\t";
				print TQ $off1,"\t",$off2,"\t0\t";
			}
			elsif($parent[1] ne $parent[2]  && $parent[1] ne $parent[4] && $parent[1] ne $parent[5]){
			#基因型为abxcd，子代类型应只有ac、ad、bc、bd，不能有aa、bb、cc、dd、ab、cd类型
				my $geno=$Ind[0].$Ind[1];
				if ($geno eq "ab" || $geno eq "cd" || $geno eq "aa" || $geno eq "bb" || $geno eq "cc" || $geno eq "dd") {
					print STDOUT "the genotype $marker have $geno is error\n";
					$off1=0;
					$off2=0;
				}else{
					if ($f1 eq "-" && $f2 eq "-") {
						$off1=0;
						if ($m1 eq "-" && $m2 eq "-") {
							$off2=0;
						}else{
							if($Ind[1] eq $m1){#第一位应来自父本，第二位应来自母本，没有亲本位置不确定情况
								$off2=1;
							}elsif($Ind[1] eq $m2){
								$off2=2;
							}else{
								$off2=0;
							}
						}
					}
					if ($m1 eq "-" && $m2 eq "-") {
						$off2=0;
						if ($f1 ne "-" && $f2 ne "-") {
							if ($Ind[0] eq $f1) {
								$off1=1;
							}elsif($Ind[0] eq $f2){
								$off1=2;
							}else{
								$off1=0;
							}
						}
					}
					if ($f1 ne "-" && $f2 ne "-" && $m1 ne "-" && $m2 ne "-") {
						if ($Ind[0] eq $f1) {#判断第一个基因型
							$off1=1;
						}elsif($Ind[0] eq $f2){
							$off1=2;
						}else{
							$off1=0;
						}
						if ($Ind[1] eq $m1) {#判断第二个基因型
							$off2=1;
						}elsif($Ind[1] eq $m2){
							$off2=2;
						}else{
							$off2=0;
						}
					}
				}
				print OUT $off1,$off2,"\t";
				print TQ $off1,"\t",$off2,"\t0\t";
			}elsif($parent[1] eq $parent[4] && $parent[2] ne $parent[5] && $parent[1] ne $parent[2] && $parent[4] ne $parent[5]){
				#基因型为efxeg，只应有ee、ef、eg、fg类型，不应有ff、gg类型
				my $geno=$Ind[0].$Ind[1];
				if ($geno eq "ff" || $geno eq "gg") {
					print STDOUT "the genotype $marker have $geno is error\n";
					$off1=0;
					$off2=0;
				}else{
					if ($f1 eq "-" && $f2 eq "-") {
						$off1=0;
						if ($m1 eq "-" && $m2 eq "-") {
							$off2=0;
						}else{
							#第一、二位不确定来自父本还是母本，有亲本位置不确定情况
							if (exists $hash{"F"}{$Ind[0]} && exists $hash{"M"}{$Ind[1]}) {
								if ($Ind[1] eq $m1) {
									$off2=1;
								}elsif($Ind[1] eq $m2){
									$off2=2;
								}else{
									$off2=0;
								}
							}elsif(exists $hash{"F"}{$Ind[1]} && exists $hash{"M"}{$Ind[0]}){
								if ($Ind[0] eq $m1) {
									$off2=1;
								}elsif($Ind[0] eq $m2){
									$off2=2;
								}else{
									$off2=0;
								}
							}else{
								$off2=0;
							}
						}
					}
					if ($m1 eq "-" && $m2 eq "-") {
						$off2=0;
						if ($f1 ne "-" && $f2 ne "-") {
							if (exists $hash{"F"}{$Ind[0]} && exists $hash{"M"}{$Ind[1]}) {
								if ($Ind[0] eq $f1) {
									$off1=1;
								}elsif($Ind[0] eq $f2){
									$off1=2;
								}else{
									$off1=0;
								}
							}elsif(exists $hash{"F"}{$Ind[1]} && exists $hash{"M"}{$Ind[0]}){
								if ($Ind[1] eq $f1) {
									$off1=1;
								}elsif($Ind[1] eq $f2){
									$off1=2;
								}else{
									$off1=0;
								}
							}else{
								$off1=0;
							}
						}
					}
					if ($f1 ne "-" && $f2 ne "-" && $m1 ne "-" && $m2 ne "-") {
						if ($Ind[0] eq $f1 && $Ind[1] eq $m1) {
							$off1=1;
							$off2=1;
						}elsif($Ind[0] eq $f2 && $Ind[1] eq $m1){
							$off1=2;
							$off2=1;
						}elsif($Ind[0] eq $f1 && $Ind[1] eq $m2){
							$off1=1;
							$off2=2;
						}elsif($Ind[0] eq $f2 && $Ind[1] eq $m2){
							$off1=2;
							$off2=2;
						}
						elsif($Ind[1] eq $f1 && $Ind[0] eq $m1) {
							$off1=1;
							$off2=1;
						}elsif($Ind[1] eq $f2 && $Ind[0] eq $m1){
							$off1=2;
							$off2=1;
						}elsif($Ind[1] eq $f1 && $Ind[0] eq $m2){
							$off1=1;
							$off2=2;
						}elsif($Ind[1] eq $f2 && $Ind[0] eq $m2){
							$off1=2;
							$off2=2;
						}else{
							$off1=0;
							$off2=0;
						}
					}
				}
				print OUT $off1,$off2,"\t";
				print TQ $off1,"\t",$off2,"\t0\t";
			}elsif($parent[1] ne $parent[2] && $parent[1] eq $parent[4] && $parent[1] eq $parent[5] ){#lmxll
				my $geno=$Ind[0].$Ind[1];
				if ($geno eq "mm") {
					print STDOUT "the genotype $marker have $geno is error\n";
					$off1=0;
					$off2=0;
				}else{
					if ($f1 eq "-" && $f2 eq "-") {
						$off1=0;
						if ($m1 eq "-" && $m2 eq "-") {
							$off2=0;
						}else{
							if ($Ind[0] eq $f1 && $Ind[1] eq $m1) {
								$off2=1;
							}elsif($Ind[0] eq $f2 && $Ind[1] eq $m2){
								$off2=2;
							}elsif($Ind[1] eq $f1 && $Ind[0] eq $m1){
								$off2=1;
							}elsif($Ind[1] eq $f2 && $Ind[0] eq $m2){
								$off2=2;
							}else{
								$off2=0;
							}
 					}
					}
					if ($m1 eq "-" && $m2 eq "-") {
						$off2=0;
						if ($f1 ne "-" && $f2 ne "-") {
							if ($Ind[0] eq $f1 && exists $hash{"M"}{$Ind[1]}) {
								$off1=1;
							}elsif($Ind[0] eq $f2 && exists $hash{"M"}{$Ind[1]}){
								$off1=2;
							}elsif($Ind[1] eq $f1 && exists $hash{"M"}{$Ind[0]}){
								$off1=1;
							}elsif($Ind[1] eq $f2 && exists $hash{"M"}{$Ind[0]}){
								$off1=2;
							}else{
								$off1=0;
							}
						}
					}
					if ($f1 ne "-" && $f2 ne "-" && $m1 ne "-" && $m2 ne "-") {
						if ($Ind[0] eq $Ind[1]) {
							if ($f1 ne $f2) {
								if ($Ind[0] eq $f1) {
									$off1=1;
								}elsif($Ind[0] eq $f2){
									$off1=2;
								}
							}else{
								$off1=0;#not sure
							}
							if ($m1 ne $m2) {
								if ($Ind[1] eq $m1) {
									$off2=1;
								}elsif($Ind[1] eq $m2){
									$off2=2;
								}
							}else{
								$off2=0;#not sure
							}
						}else{
							if (exists $hash{"F"}{$Ind[0]} && exists $hash{"M"}{$Ind[1]}) {
								if ($f1 eq $f2) {
									$off1=0;#not sure
								}else{
									if ($Ind[0] eq $f1) {
										$off1=1;
									}elsif($Ind[0] eq $f2){
										$off1=2;
									}
								}
								if ($m1 eq $m2) {
									$off2=0;#not sure
								}else{
									if ($Ind[1] eq $m1) {
										$off2=1;
									}elsif($Ind[1] eq $m2){
										$off2=2;
									}
								}
							}elsif(exists $hash{"F"}{$Ind[1]} && exists $hash{"M"}{$Ind[0]}){
								if ($f1 eq $f2) {
									$off1=0;#not sure
								}else{
									if ($Ind[1] eq $f1) {
										$off1=1;
									}elsif($Ind[1] eq $f2){
										$off1=2;
									}
								}
								if ($m1 eq $m2) {
									$off2=0;#not sure
								}else{
									if ($Ind[0] eq $m1) {
										$off2=1;
									}elsif($Ind[0] eq $m2){
										$off2=2;
									}
								}
							}
						}
					}
				}
				print OUT $off1,$off2,"\t";
				print TQ $off1,"\t",$off2,"\t0\t";
			}elsif($parent[1] eq $parent[2] && $parent[1] eq $parent[4] && $parent[1] ne $parent[5]){#nnxnp
				my $geno=$Ind[0].$Ind[1];
				if ($geno eq "pp") {
					print STDOUT "the genotype $marker have $geno is error\n";
					$off1=0;
					$off2=0;
				}else{
					if ($f1 eq "-" && $f2 eq "-") {
						$off1=0;
						if ($m1 eq "-" && $m2 eq "-") {
							$off2=0;
						}else{
							if ($Ind[0] eq $m1 && exists $hash{"F"}{$Ind[1]}) {
								$off2=1;
							}elsif($Ind[0] eq $m2 && exists $hash{"F"}{$Ind[1]}){
								$off2=2;
							}elsif($Ind[1] eq $m1 && exists $hash{"M"}{$Ind[0]}){
								$off2=1;
							}elsif($Ind[1] eq $m2 && exists $hash{"M"}{$Ind[0]}){
								$off2=2;
							}else{
								$off2=0;
							}
 						}
					}
					if ($m1 eq "-" && $m2 eq "-") {
						$off2=0;
						if ($f1 ne "-" && $f2 ne "-") {
							if ($Ind[0] eq $f1 && $Ind[1] eq $m1) {
								$off1=1;
							}elsif($Ind[0] eq $f2 && $Ind[1] eq $m2){
								$off1=2;
							}elsif($Ind[1] eq $f1 && $Ind[0] eq $m1){
								$off1=1;
							}elsif($Ind[1] eq $f2 && $Ind[0] eq $m2){
								$off1=2;
							}else{
								$off1=0;
							}
						}
					}
					if ($f1 ne "-" && $f2 ne "-" && $m1 ne "-" && $m2 ne "-") {
						if ($Ind[0] eq $Ind[1]) {
							if ($f1 ne $f2) {
								if ($Ind[0] eq $f1) {
									$off1=1;
								}elsif($Ind[0] eq $f2){
									$off1=2;
								}
							}else{
								$off1=0;#not sure
							}
							if ($m1 ne $m2) {
								if ($Ind[1] eq $m1) {
									$off2=1;
								}elsif($Ind[1] eq $m2){
									$off2=2;
								}
							}else{
								$off2=0;#not sure
							}
						}else{
							if (exists $hash{"F"}{$Ind[0]} && exists $hash{"M"}{$Ind[1]}) {
								if ($f1 eq $f2) {
									$off1=0;#not sure
								}else{
									if ($Ind[0] eq $f1) {
										$off1=1;
									}elsif($Ind[0] eq $f2){
										$off1=2;
									}
								}
								if ($m1 eq $m2) {
									$off2=0;#not sure
								}else{
									if ($Ind[1] eq $m1) {
										$off2=1;
									}elsif($Ind[1] eq $m2){
										$off2=2;
									}
								}
							}elsif(exists $hash{"F"}{$Ind[1]} && exists $hash{"M"}{$Ind[0]}){
								if ($f1 eq $f2) {
									$off1=0;#not sure
								}else{
									if ($Ind[1] eq $f1) {
										$off1=1;
									}elsif($Ind[1] eq $f2){
										$off1=2;
									}
								}
								if ($m1 eq $m2) {
									$off2=0;#not sure
								}else{
									if ($Ind[0] eq $m1) {
										$off2=1;
									}elsif($Ind[0] eq $m2){
										$off2=2;
									}
								}
							}
						}
					}
				}
				print OUT $off1,$off2,"\t";
				print TQ $off1,"\t",$off2,"\t0\t";
			}
		}
		print OUT "\n";
		print TQ "\n";
	}
	close(IN);
	close(OUT);
	close(TQ);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:IndividualLinkagePhaseNumber.pl
call:NULL
Author: wushuang <wus\@biomarker.com.cn> 
Version: $version
Date: 2012-1-16 更改函数骨架
Description:
	根据亲本的连锁相和基因型，统计子代个体的基因型的亲本类型，loc文件要包含连锁相信息
	-：缺失 gray bar
	0: 无法判断 white bar 
	1：子代基因型来自亲本第一个等位基因型 green bar
	2：子代基因型来自亲本第二个等位基因型 blue bar
Usage:
  Options:
  -i <file>   Input loc file, forced
  -m <file>   Input map file, forced
  -o <key>    Output file key 
              *.phase          转换的连锁相文件
              *.tq             画图的数据
              *.tq.png         图片
              *.original.*     原始单体来源
              *.repairDis.*    区分缺失和无法判断
              *.repair.*       不区分缺失和无法判断
  -h         Help

USAGE
	print $usage;
	exit;
}
