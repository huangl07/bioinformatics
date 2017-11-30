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

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn1,$fIn2,$fIn3,$outdir,$fOut,$fOut1,$Reverse,$switch);
GetOptions(
				"help|?" =>\&USAGE,
				"i1:s"=>\$fIn1,
				"i2:s"=>\$fIn2,
				"i3:s"=>\$fIn3,
				"o:s"=>\$fOut,
				"o1:s"=>\$fOut1,
				"d:s"=>\$outdir,
				"r:s"=>\$Reverse,
				"s:s"=>\$switch,
				) or &USAGE;
$switch||="no";
&USAGE unless ($switch eq "run");

#------------------------------------------------read in files
my %all;
if (!defined $fIn1 and !defined $fIn2 and !defined $fIn3) {#without input files, search current diretory for *.male.map *.female.map *.sexAver.map
	my @files;my %exame;my @names;my @useful;
	@files=glob ("./*.map");
	foreach my $file (@files) {
		my ($key)=$file=~/\.\/(.*)\.\w+\.map/;
		push @names,$key;
		$exame{$file}=1;
	}
	my %saw;
	@names=grep(!$saw{$_}++,@names);#delete the same key of files
	for (my $i=0;$i<@names ;$i++) {#judge if a given key,for example AAA,has corresponding files,exactly like AAA.male.map AAA.female.map AAA.sexAver.map 
		if (defined $exame{"./".$names[$i].".male.map"} && defined $exame{"./$names[$i].female.map"}  && defined $exame{"./$names[$i].sexAver.map"}){
			push @useful,$names[$i] ;
		}else{
			print STDOUT "Can't draw LG$names[$i] for lacking correspongding male information!\n" if (!defined $exame{"./".$names[$i].".male.map"}) ;
			print STDOUT "Can't draw LG$names[$i] for lacking correspongding female information!\n" if (!defined $exame{"./".$names[$i].".female.map"}) ;
			print STDOUT "Can't draw LG$names[$i] for lacking correspongding sexAver information!\n" if (!defined $exame{"./".$names[$i].".sexAver.map"}) ;
		}
	}
FOR:foreach my $useful (@useful) {#get the information and judge whether the male,female,sexAver three files are related 
		$/="\n";

		$fIn1="./$useful.sexAver.map";
		my $linesCount1="linesCount1";
		`wc  "$fIn1"  >$linesCount1 `;
		open (L1,$linesCount1) or die $!;
		if (<L1>=~/^0/) {
			print "$fIn1 is empty please chech it and corresponding files!\n";
			next FOR;
		}
		open (IN1,$fIn1) or die $!;
		while (<IN1>) {
			chomp;
			next if /^\s/ || /^;/;  
			$all{$useful}{"7"}{$1}=$2 if ( /^(\S+)\s+(\S+)/);
		} 
		close(IN1);

		$fIn2="./$useful.female.map";
		my $linesCount2="linesCount2";
		`wc  "$fIn2"  >$linesCount2 `;
		open (L2,$linesCount2) or die $!;											 
		if (<L2>=~/^0/) {
			print "$fIn2 is empty please chech it and corresponding files!\n";
			delete 	$all{$useful};
			next FOR;
		}
		open (IN2,$fIn2) or die $!;
WHILE1:	while (<IN2>) {
			chomp;
			next if /^\s/ || /^;/;
			if (/^(\S+)\s+(\S+)/){
				if (!exists $all{$useful}{"7"}{$1}) {
					print STDOUT "the $fIn2 file contains markers like $1 that can't been foud in sexAver file, please enter y to continue or n to stop this png\n";
					my $option=<STDIN>;
					chomp ($option);
					if ($option eq  "y"){
						$all{$useful}{"5"}{$1}=$2;
						my @restlines=<IN2>;
						for (my $i=0;$i<@restlines ;$i++) {
							chomp ($restlines[$i]);					
							next if $restlines[$i]=~/^\s/ || /^;/;
							$restlines[$i]=~/^(\S+)\s+(\S+)/;
							$all{$useful}{"5"}{$1}=$2;
						}											
					}elsif($option eq "n"){
						delete $all{$useful};
						next FOR;
					}else{
						print "please enter an appropriate  option as asked\n";
						redo WHILE1;
					}

				}else{
				  $all{$useful}{"5"}{$1}=$2;
				}
			}
		}
		close(IN2);
		$fIn3="./$useful.male.map";
		my $linesCount3="linesCount3";
		`wc  "$fIn3"  >$linesCount3 `;
		open (L3,$linesCount3) or die $!;
		if (<L3>=~/^0/) {
			print "$fIn3 is empty please chech it and corresponding files!\n";
			delete 	$all{$useful};
			next FOR;
		}
		open (IN3,$fIn3) or die $!;
WHILE2:	while (<IN3>) {
			chomp;
			next if /^\s/ || /^;/;	
			if (/^(\S+)\s+(\S+)/){
				if (!exists $all{$useful}{"7"}{$1}) {
					print STDOUT "the $fIn3 file contains markers like $1 that can't been foud in sexAver file, please enter y to continue or n to stop this png\n";
					my $option=<STDIN>;
					chomp ($option);
					if ($option eq  'y'){
						$all{$useful}{"9"}{$1}=$2;
						my @restlines=<IN3>;
						for (my $i=0;$i<@restlines ;$i++) {
							chomp ($restlines[$i]);
							next if $restlines[$i]=~/^\s/ || /^;/;
							$restlines[$i]=~/^(\S+)\s+(\S+)/;
							$all{$useful}{"9"}{$1}=$2;
						}	
					}elsif($option eq 'n'){
						delete 	 $all{$useful};
						next FOR;
					}else{
						print "please enter an appropriate  option as asked\n";
						redo WHILE2;
					}

				}else{
				  $all{$useful}{"9"}{$1}=$2;
				}
			}
		} 
		close(IN3);

	}
}else{#with input files,get the information and judge whether the groups with the same number are related,judge if the file is empty 
	&USAGE unless ($switch eq "run" and $fIn1 and $fIn2 and $fIn3);
	print "get input files: $fIn1\t$fIn2\t$fIn3\n";
	my ($file1_group,$file2_group,$file3_group); 
	$/="\n";
	my $linesCount2="linesCount2";
		`wc  "$fIn2"  >$linesCount2 `;
		open (L2,$linesCount2) or die $!;
		if (<L2>=~/^0/) {
			print "$fIn2 is empty please chech it and corresponding files!\n";
			`rm linesCount?` ;
			die;
		}
	open (IN2,$fIn2) or die $!;
	while (<IN2>) {
		chomp;
		next if /^\s/ || /^;/;  
		$file2_group=$1 if /^group\s+(\d+)/;
		$all{$file2_group}{"7"}{$1}=$2 if ( /^(\S+)\s+(\S+)/);
	} 
	close(IN2);
	
	open (IN1,$fIn1) or die $!;
	my $linesCount1="linesCount1";
		`wc  "$fIn1"  >$linesCount1 `;
		open (L1,$linesCount1) or die $!;
		if (<L1>=~/^0/) {
			print "$fIn1 is empty please chech it and corresponding files!\n";
			`rm linesCount?` ;
			die;
		}
F1:	while (<IN1>) {
		chomp;
		next if /^\s/ || /^;/;
		$file1_group=$1 if /^group\s+(\d+)/;
		if  ( /^(\S+)\s+(\S+)/){
			if (exists $all{$file1_group}{"7"}{$1}) {
				$all{$file1_group}{"5"}{$1}=$2;
			}else{
				print "LG$file1_group of $fIn1 contains marker:$1 can't been found in corresponding group in sexAver file,please enter y to continue or n to stop this png\n";
				my $option=<STDIN>;
				chomp ($option);
				if ($option eq 'y') {
					$all{$file1_group}{"5"}{$1}=$2;
					while (<IN1>){
						chomp;
						next if (/^;/ || /^$/); 
						/^(\S+)\s+(\S+)/;
						$all{$file1_group}{"5"}{$1}=$2;
						if (/^group\s+(\d+)/){
							$file1_group=$1;
							last;
						}
					}
				}elsif($option eq 'n'){					
						while (<IN1>){
							if (/^group\s+(\d+)/){
								$file1_group=$1;
								last;
							}
						}
					delete $all{$file1_group};
					next;
				}else{
					print "please enter an appropriate option as asked\n";
					redo F1;
				}
			}
		}
	}
	close(IN1);

		my $linesCount3="linesCount3";
		`wc  "$fIn3"  >$linesCount3 `;
		open (L3,$linesCount3) or die $!;
		if (<L3>=~/^0/) {
			print "$fIn3 is empty please chech it and corresponding files!\n";
			`rm linesCount?` ;
			die;
		}
	open (IN3,$fIn3) or die $!;
F2:	while (<IN3>) {
		chomp;
		next if /^\s/ || /^;/;
		$file3_group=$1 if /^group\s+(\d+)/;
		if  (/^(\S+)\s+(\S+)/){
			if (exists $all{$file3_group}{"7"}{$1}) {
				$all{$file3_group}{"9"}{$1}=$2;
			}else{
				print "$file3_group of $fIn3 contains marker:$1 can't been found in corresponding group in sexAver file,please enter y to continue or n to stop this png\n";
				my $option=<STDIN>;
				chomp ($option);
				if ($option eq 'y') {
					$all{$file3_group}{"9"}{$1}=$2;
					while (<IN3>){
						chomp;
						next if (/^;/ || /^$/); 
						/^(\S+)\s+(\S+)/;
						$all{$file3_group}{"9"}{$1}=$2;
						if (/^group\s+(\d+)/){
							$file3_group=$1;
							last;
						}
					}
				}elsif($option eq 'n'){
					while (<IN3>){
						if (/^group\s+(\d+)/){
							$file3_group=$1;
							last;
						}
					}
						delete $all{$file3_group};
						next;
				}else{
					print "please enter an appropriate option as asked\n";
					redo F2;
				}


			}
		}
	} 
	close(IN3);

	
}
print "read data done!!!\n";
#print Dumper %all;die;
#print Dumper %all;die;
$Reverse||=0;
$outdir||="synteny";
`mkdir $outdir` unless (-d $outdir); 
my $workdir=`pwd`;chomp $workdir;

#print Dumper %all;die;
# ------------------------------------------------------------------
# Draw SVG parameter
# ------------------------------------------------------------------
my $LeftMargin=50;
my $RightMargin=50;
my $TopMargin=50;
my $BottomMargin=50;

my $TextColor="black";
my $TextSize="12";
my $TextWidth=20;
my $TextSizeX="14";

my $TitleHeight=10;
my $titleHeight=10;
my $Title_title_space=20;
my $title_chro_space=10;

my $MarkerTxtSpace=70;
my $MarkerLineSpace=20;
my $Marker_LineSpace=40;
my $CmTxtSpace=60;
#my $CmLineSpace=30;

my $LineSize=1;
my $LineColor="#6B8E23";
my $ShortLineColor="black";

my $space_betw_chro=180;

my $ChroWidth=35;
my $ChroColor="#999900";
my $factor=10;
my $MiddleColor="rgb(247,251,209)";

#-------------------------------------------------得到三个文件中相应group的最大遗传距离，并存到一个数组中    
my %max_chro_length;
foreach my $groups (sort keys %all) { 
	my $max=0;
	foreach my $whichK (sort {$a <=> $b} keys %{$all{$groups}}) {
		my @marker=(sort {$all{$groups}{$whichK}{$a} <=> $all{$groups}{$whichK}{$b}} keys  %{$all{$groups}{$whichK}} );
		 $max=$all{$groups}{$whichK}{$marker[-1]} ;
		 push @{$max_chro_length{$groups}},$max;
	}
	unshift  @{$max_chro_length{$groups}},"0" if (!defined $all{$groups}{'5'});
	push   @{$max_chro_length{$groups}},"0" if (!defined $all{$groups}{'9'});
} 

my ($pos,$gong)=&get_colinearity_data(\%all);
my %pos=%$pos;my %gong=%$gong;
foreach my $groups (sort keys %all) {
	print STDOUT "LG$groups is drawing......\n";
	my $png = defined $fOut ? "$fOut.$groups":$groups;  
	open (SVG,">","$outdir/$png.FAM.svg") or die $!;
	my @tem=@{$max_chro_length{$groups}};
	my @temp = sort {$a <=> $b} @tem;
#	print @tem,"\n";die;
	my $max_length=$temp[-1];
	$factor= (keys %{$all{$groups}{'7'}})*$TextSizeX >=$max_length*$factor ? ((keys %{$all{$groups}{'7'}})*$TextSizeX)/$max_length:$factor;
	$factor= (keys %{$all{$groups}{'5'}})*$TextSizeX >=$tem[0]*$factor ? ((keys %{$all{$groups}{'5'}})*$TextSizeX)/$tem[0]:$factor;
	$factor= (keys %{$all{$groups}{'9'}})*$TextSizeX >=$tem[2]*$factor ? ((keys %{$all{$groups}{'9'}})*$TextSizeX)/$tem[2]:$factor;
	my $PaperWidth=$LeftMargin+$RightMargin+($MarkerTxtSpace*3+$MarkerLineSpace*6+$ChroWidth*3+$space_betw_chro*2)+$CmTxtSpace*3+$Marker_LineSpace*6;
	my $PaperHeight=$TopMargin+$BottomMargin+$max_length*$factor+$ChroWidth+$TitleHeight+$titleHeight+$Title_title_space+$title_chro_space;
	print SVG &svg_paper($PaperWidth,$PaperHeight);
	print SVG "<defs>\n<linearGradient id=\"linked_group\"  x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n<stop offset=\"0%\" style=\"stop-color:$ChroColor;stop-opacity:1\"/>\n<stop offset=\"50%\" style=\"stop-color:$MiddleColor;stop-opacity:1\"/>\n<stop offset=\"100%\" style=\"stop-color:$ChroColor;stop-opacity:1\"/>\n</linearGradient>\n</defs>","\n";
	print SVG &svg_txt($PaperWidth/2,$TopMargin+$TitleHeight,$TextSize,$TextColor,"LG $png");
	my $x1=$LeftMargin+$MarkerLineSpace+$MarkerTxtSpace+$Marker_LineSpace;
	my $x2=$LeftMargin+$MarkerLineSpace*3+$MarkerTxtSpace*2+$space_betw_chro+$ChroWidth+$CmTxtSpace+$Marker_LineSpace*3;
	my $x3=$LeftMargin+$MarkerLineSpace*5+$MarkerTxtSpace*3+$space_betw_chro*2+$ChroWidth*2+$CmTxtSpace*2+$Marker_LineSpace*5;

	
	print SVG &svg_txt($x1+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$Title_title_space,$TextSize,$TextColor,"$png.female");
	print SVG &svg_circle($x1+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ChroWidth/2,$ChroColor,"linked_group"),"\n";
	print SVG &svg_circle($x1+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+($tem[0])*$factor+$Title_title_space+$title_chro_space,$ChroWidth/2,$ChroColor,"linked_group"),"\n"; 
	print SVG &svg_rect($x1,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ChroWidth,$tem[0]*$factor,$ChroColor,"linked_group") ;
	
	print SVG &svg_txt($x2+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$Title_title_space,$TextSize,"$TextColor","$png.sexAver");
	print SVG &svg_circle($x2+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ChroWidth/2,$ChroColor,"linked_group"),"\n";
	print SVG &svg_circle($x2+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+($tem[1])*$factor+$Title_title_space+$title_chro_space,$ChroWidth/2,$ChroColor,"linked_group"),"\n"; 
	print SVG &svg_rect($x2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ChroWidth,$tem[1]*$factor,$ChroColor,"linked_group") ;
	
	print SVG &svg_txt($x3+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$Title_title_space,$TextSize,"$TextColor","$png.male");
	print SVG &svg_circle($x3+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ChroWidth/2,$ChroColor,"linked_group"),"\n";
	print SVG &svg_circle($x3+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+($tem[2])*$factor+$Title_title_space+$title_chro_space,$ChroWidth/2,$ChroColor,"linked_group"),"\n"; 
	print SVG &svg_rect($x3,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ChroWidth,$tem[2]*$factor,$ChroColor,"linked_group") ;
	

	foreach my $whichK (sort {$a <=> $b} keys %{$all{$groups}}) {
	
		my @Marker=keys %{$all{$groups}{$whichK}};
		my $k=0;
		foreach my $marker(sort {$all{$groups}{$whichK}{$a} <=> $all{$groups}{$whichK}{$b}} keys  %{$all{$groups}{$whichK}} ){
			my $y=$all{$groups}{$whichK}{$marker}*$factor;
			my $Y=$all{$groups}{$whichK}{$marker};

			if ($whichK eq "5") {
				my $length=$tem[0]/(scalar @Marker)*$factor;
				my $MY=$length*$k;
				print SVG &svg_txt($LeftMargin,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$TextSize,"$TextColor","$marker");
				print SVG &svg_line($x1-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x1,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
				print SVG &svg_line($x1-$MarkerLineSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x1-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);

				print SVG &svg_txt($x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$TextSize,"$TextColor","$Y Cm");
				print SVG &svg_line($x1+$ChroWidth,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x1+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
				print SVG &svg_line($x1+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
				$pos{$groups}{"5"}{$marker}=$MY;
				$k++;
			}elsif($whichK eq "7") {
				my $length;
				$length=$tem[1]/(scalar @Marker)*$factor;
				my $MY=$length*$k;
				if ($pos{$groups}{"7k"} eq "reverse") {
					my $bottomY=$tem[1]*$factor+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space;
					print SVG &svg_txt($x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x2-$MarkerLineSpace,$bottomY-$y,$x2,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2-$MarkerLineSpace-$Marker_LineSpace,$bottomY-$MY,$x2-$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$Y Cm");
					print SVG &svg_line($x2+$ChroWidth,$bottomY-$y,$x2+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$ShortLineColor,$LineSize);
					$pos{$groups}{"7Ynew"}{$marker}=$bottomY-$MY;
				}else{
					print SVG &svg_txt($x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x2-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x2,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2-$MarkerLineSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x2-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$TextSize,"$TextColor","$Y Cm");
					print SVG &svg_line($x2+$ChroWidth,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x2+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
					$pos{$groups}{"7ynoreverse"}{$marker}=$MY;
				}
			$k++;
			}elsif ($whichK eq "9"){
				my $length;my $chrolength;
				$length=$tem[2]/(scalar @Marker)*$factor;
				$chrolength=$tem[2];
				my $MY=$length*$k;
				if ($pos{$groups}{"9k"} eq "reverse") {
					my $bottomY=$chrolength*$factor+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space;
					print SVG &svg_txt($x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x3-$MarkerLineSpace,$bottomY-$y,$x3,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3-$MarkerLineSpace-$Marker_LineSpace,$bottomY-$MY,$x3-$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$Y Cm");
					print SVG &svg_line($x3+$ChroWidth,$bottomY-$y,$x3+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$ShortLineColor,$LineSize);
					$pos{$groups}{"9Ynew"}{$marker}=$bottomY-$MY;
				}else{
					print SVG &svg_txt($x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x3-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x3,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3-$MarkerLineSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x3-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$TextSize,"$TextColor","$Y Cm");
					print SVG &svg_line($x3+$ChroWidth,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x3+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$ShortLineColor,$LineSize);
					$pos{$groups}{"9ynoreverse"}{$marker}=$MY;
				}
				$k++;
			}
		}
	}
	my @marker1=keys %{$all{$groups}{"5"}};
	my @marker2=keys %{$all{$groups}{"7"}};
	foreach my $marker1 (@marker1) {
		if (defined $pos{$groups}{"7Ynew"}{$marker1}) {
			print SVG &svg_line($x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"5"}{$marker1}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"7Ynew"}{$marker1},$LineColor,$LineSize);
		}elsif (defined $all{$groups}{"7"}{$marker1}) {
			print SVG &svg_line($x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"5"}{$marker1}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"7ynoreverse"}{$marker1}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$LineColor,$LineSize);
		}
	}
	if ($pos{$groups}{"7k"} eq "reverse") {
		foreach my $marker2 (@marker2) {
			if (defined $pos{$groups}{"9Ynew"}{$marker2}) {
				print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"7Ynew"}{$marker2},$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"9Ynew"}{$marker2},$LineColor,$LineSize);
			}elsif (defined $all{$groups}{"9"}{$marker2}) {
				print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"7Ynew"}{$marker2},$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"9ynoreverse"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$LineColor,$LineSize);
			}
		}
	}else{
		foreach my $marker2 (@marker2) {
			if (defined $pos{$groups}{"9Ynew"}{$marker2}) {
				print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"7ynoreverse"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"9Ynew"}{$marker2},$LineColor,$LineSize);
			}elsif (defined $all{$groups}{"9"}{$marker2}) {
				 print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"7ynoreverse"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"9ynoreverse"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space,$LineColor,$LineSize);
			}
		}
	}
	print SVG &svg_end ();
	close (SVG);
	chdir $outdir;
	`perl $Bin/svg2xxx_release/svg2xxx  "$png.FAM.svg" -t -w $PaperWidth -h $PaperHeight `;
	chdir $workdir;
}
foreach my $groups (sort keys %gong) {
	my $png = defined $fOut ? $fOut:$groups;  
	open (STAT,">$outdir/$png.FAM.synteny.stat.xls") or die $!;
	print STAT "LG$png\n";
	print STAT"${$gong{$groups}}[0]\/${$gong{$groups}}[2]\t${$gong{$groups}}[1]\/${$gong{$groups}}[3]\n";
}

`rm linesCount?`;


















#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub get_colinearity_data {#  共线性数量，公共marker数量的计算及染色体是否要翻转的判断
	my %copyall=%{$_[0]};
	my %position;my %gongline;
	foreach my $groups (sort keys %copyall) {
		if ($Reverse==0) {
			my ($gongline1,$gongline2,$a1,$a2);
			my (@order1,@order11,@order2,@order22);
			my $k=0;my $i=0;my $j=0;my $p=0;
			my @marker2=(sort {$copyall{$groups}{"7"}{$a} <=> $copyall{$groups}{"7"}{$b}} keys  %{$copyall{$groups}{"7"}} );
			my @Rmarker2=reverse @marker2;
			my @marker1=(sort {$copyall{$groups}{"5"}{$a} <=> $copyall{$groups}{"5"}{$b}} keys  %{$copyall{$groups}{"5"}} );
			my @marker3=(sort {$copyall{$groups}{"9"}{$a} <=> $copyall{$groups}{"9"}{$b}} keys  %{$copyall{$groups}{"9"}} );
			my @Rmarker3=reverse @marker3;
			foreach my $marker2 (@marker2) {
				 $k++;
				$position{$groups}{"7"}{$marker2}=$k;
				
			}
			foreach my $Rmarker2 (@Rmarker2) {
				 $j++;
				$position{$groups}{"77"}{$Rmarker2}=$j;
				
			}
			foreach my $marker3 (@marker3) {
				 $p++; 
				$position{$groups}{"9"}{$marker3}=$p;
			
			}
			foreach my $Rmarker3 (@Rmarker3) {
				 $i++;
				$position{$groups}{"99"}{$Rmarker3}=$i;
			
			}
			foreach my $marker1 (@marker1) {
				if (defined $position{$groups}{"7"}{$marker1}) {
				push  @order1,$position{$groups}{"7"}{$marker1};
				}
				if (defined $position{$groups}{"77"}{$marker1}) {
				push @order11,$position{$groups}{"77"}{$marker1};
				}
			}
			$a1=@order1;
			if (scalar @order1==0 || scalar @order11==0){
				$position{$groups}{"7k"}="not reverse";
			}else{
				my @result1=PickMaxGroup(@order1);
				my @result11=PickMaxGroup(@order11);

		#print scalar @result1,"\n",scalar @result11,"\n";
				if (scalar @result1 <  scalar @result11) {
					$gongline1=scalar @result11;
					$position{$groups}{"7k"}="reverse";
					@marker2=@Rmarker2;
				}else{
					$position{$groups}{"7k"}="not reverse";
					$gongline1=scalar @result1;
				}	
			}
			foreach my $marker2 (@marker2) {
				if (defined $position{$groups}{"9"}{$marker2}) {
				push  @order2,$position{$groups}{"9"}{$marker2};
				}
				if (defined $position{$groups}{"99"}{$marker2}) {
				push  @order22,$position{$groups}{"99"}{$marker2};
				}
			}
			$a2=@order2;
			if (scalar @order2==0 || scalar @order22==0 ) {
				$position{$groups}{"9k"}="not reverse";
			}else{
				my  @result2=PickMaxGroup(@order2);
				my  @result22=PickMaxGroup(@order22);
				if (scalar @result2 <  scalar @result22) {
					$gongline2=scalar @result22;
					$position{$groups}{"9k"}="reverse";
				}else{
					$position{$groups}{"9k"}="not reverse";
					$gongline2=scalar @result2;
				}	
			}
			$gongline1= !defined $gongline1 ? 0:$gongline1; 
			$gongline2= !defined $gongline2 ? 0:$gongline2;
			$a1= !defined $a1 ? 0:$a1;
			$a2= !defined $a2 ? 0:$a2;
			push @{$gongline{$groups}},$gongline1;
			push @{$gongline{$groups}},$gongline2;
			push @{$gongline{$groups}},$a1;
			push @{$gongline{$groups}},$a2;
		}else{
		$position{$groups}{"9k"}="not reverse";
		$position{$groups}{"7k"}="not reverse";
		}
		
	}
	return (\%position,\%gongline);
}

sub PickMaxGroup
{
    my(@Pos)=@_;
    
    my @dp;#[$pos][0]-total [$pos][1]-prePos
    my $MaxEnd=0;
    my $MaxScore=0;
    for(my $i=0;$i<@Pos;$i++)
    {
        my $max_total=1;
        my $pre=-1;
        for(my $j=$i-1;$j>=0;$j--)
        {
            if($Pos[$i]>=$Pos[$j])
            {
                if($dp[$j][0]>=$max_total)
                {
                    $max_total=1+$dp[$j][0];
                    $pre=$j;
                    #last;
                }
            }
        }
        $dp[$i][0]=$max_total;
        $dp[$i][1]=$pre;
        if($MaxScore<$max_total)
        {
            $MaxScore=$max_total;
            $MaxEnd=$i;
        }
    }
    
    #find max
    my @result;
    for(my $i=$MaxEnd;$i>=0;)
    {
        push @result,$i;
        $i=$dp[$i][1];
    }
    return reverse @result;
}

sub svg_paper (){#&svg_paper(width,height,[color])
	my $svg_drawer = "litc"."@"."biomarker\.com\.cn";
	chomp $svg_drawer;
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\">\n";
	$line.="<Drawer>$svg_drawer</Drawer>\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub svg_end (){#
	return "</svg>\n";
}

sub svg_txt (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_line (){#&svg_line(x1,y1,x2,y2,color,width,[opacity])
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[6]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" opacity=\"$svg_x[6]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_polyline (){#colorfill,colorstroke,width,\@point
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_width=shift(@svg_x);
	my $svg_points=join(" ",@{$svg_x[-1]});
	my $line="<polyline fill=\"$svg_color\" stroke=\"$svg_color2\" stroke-width=\"$svg_width\" points=\"$svg_points\"/>\n";

	#print "$line\n";
	return $line;

	#<polyline points="0,0 0,20 20,20 20,40 40,40 40,60" style="fill:white;stroke:red;stroke-width:2"/>
}

sub svg_rect () {#&svg_rect(x,y,width,height,color,Id,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[6]) {
		$svg_x[6]=1;
	}
	my $line="<rect style=\"fill:url(#$svg_x[5])\" x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[6]\"/>\n";
	return $line;
}

sub svg_rect1 () {#&svg_rect(x,y,width,height,color,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[6]) {
		$svg_x[6]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[6]\"/>\n";
	return $line;
}
sub svg_polygon () {#colorfill,colorstroke,coloropacity,point1,point2,...
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_points=join(" ",@svg_x);
	my $line="<polygon fill=\"$svg_color\" stroke=\"$svg_color2\" opacity=\"$svg_trans\" points=\"$svg_points\"/>\n";
	return $line;
}

sub svg_ellipse () {#&svg_ellipse(cx,cy,rx,ry,colorfill,colorstroke,width,[coloropacity])
	my @svg_x=@_;
	my $line= "<ellipse cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" rx=\"$svg_x[2]\" ry=\"$svg_x[3]\" fill=\"$svg_x[4]\" stroke=\"$svg_x[5]\" stroke-width=\"$svg_x[6]\"/>\n";
	if (defined $svg_x[7]) {
		$line="<ellipse cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" rx=\"$svg_x[2]\" ry=\"$svg_x[3]\" fill=\"$svg_x[4]\" stroke=\"$svg_x[5]\" stroke-width=\"$svg_x[6]\" opacity=\"$svg_x[7]\"/>\n";
	}
	return $line;
}

sub svg_circle () {#&svg_circle(cx,cy,r,color)
	my @svg_x=@_;
	my $line="<circle style=\"fill:url(#$svg_x[4])\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" r=\"$svg_x[2]\" stroke=\"$svg_x[3]\" stroke-width=\"0\" fill=\"$svg_x[3]\"/>";
	return $line;
}


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:draw_ThreeGroup.pl
Version: version1.1
Contact: Xin HuaiGen <xinhg\@biomarker.com.cn>
Description:you can use it like:
	        perl draw_ThreeGroup.pl  -s run
	        or:
	        perl draw_ThreeGroup.pl  -i1 *.female.map -i2 *.sexAver.map -i3 *.male.map  -o xxx  -s run

Usage:map.pl <可用于画三个图谱间的共线性关系，从左往右依次是雌中雄,运行时不输入文件，会在当前路径下找形如*.female.map、*.male.map、*.sexAver.map的文件，并绘图>
  Options:
  -i1 <file>  *.map file1,format,with less markers,usually the *.female.map
  -i2 <file>  *.map file2,format,whth most markers,usually the *.sexAver.map
  -i3 <file>  *.map file3,format,whit less markers,usually the *.male.map
  -o <string>  key of output file,if not given,output file would be named after group number 
  -d <string>  output dir, default,"synteny"
  -r <number>  optional,0 means judge whether reverse,any other numbers mean don\'t reverse straightly ,default,0
  -s <control> must be given as 'run' or it won't run
  -h         Help

USAGE
	print $usage;
	exit;
}														
