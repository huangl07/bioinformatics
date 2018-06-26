#!/bin/bash
#just a small bash script to calculate BFs for all SNPs from SNPFILE
#please copy this script into the same directory as bayenv and execute it there
#please see the Bayenv2 manual for details about usage
#make this script executable (chmod +x calc_bf.sh)
#Usage: ./calc_bf.sh <Name of your SNPSFILE> <Name of your ENVFILE> <Name of your MATFILE> <Nuber of populations> <Number of MCMC iterations> <Number of environmental factors>

SNPFILE=$1
ENVFILE=$2
MATFILE=$3
POPNUM=$4
ITNUM=$5
ENVNUM=$6

split -a 10 -d -l 2 $SNPFILE snp_batch
for f in $(ls snp_batch*)
    do
        bayenv2 -i $f -m $MATFILE -e $ENVFILE -p $POPNUM -k $ITNUM -n $ENVNUM -t -r 429 -o pop
        done

        rm -f snp_batch*
