#!/usr/bin/perl -w
# Author: Thomas Thiel
# Program name: primer3_in.pl
# Description: creates a PRIMER3 input file based on SSR search results

open (IN,"<$ARGV[0]") || die ("\nError: Couldn't open misa.pl results file (*.misa) !\n\n");

my $filename = $ARGV[0];
$filename =~ s/\.misa//;
open (SRC,"<$filename") || die ("\nError: Couldn't open source file containing original FASTA sequences !\n\n");
open (OUT,">$filename.p3in");

undef $/;
$in = <IN>;
study $in;

$/= ">";

my $count;
while (<SRC>)
  {
  next unless (my ($id,$seq) = /(.*?)\n(.*)/s);
  $id=(split(/\s+/,$id))[0];
  $seq =~ s/[\d\s>]//g;#remove digits, spaces, line breaks,...
  while ($in =~ /$id\t(\d+)\t\S+\t\S+\t(\d+)\t(\d+)/g)
    {
    my ($ssr_nr,$size,$start) = ($1,$2,$3);
    $count++;
	my $right;
	my $left=$start;
	if ($start < 300) {
		$left=0;
	}
	if ($left+300 > length($seq)) {
		$left=length($seq);
	}
	my $newseq=substr($seq,$right,$left-$right);
    print OUT "PRIMER_SEQUENCE_ID=$id"."_$start\nSEQUENCE=$newseq\n";
	print Out "PRIMER_TASK=generic"
    print OUT "PRIMER_PRODUCT_SIZE_RANGE=100-280\n";
    print OUT "TARGET=",$right,",",$size+6,"\n";
	print Out "SEQUENCE_INTERNAL_EXCLUDED_REGION=$right,$size"
	print Out "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,$right-$start,",$right+$size,"$left-$start-$size,\n";
    print OUT "PRIMER_MAX_END_STABILITY=250\n=\n"
    };
  };
print "\n$count records created.\n";