#!/usr/bin/env Rscript
times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'kmer','k',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript base_distrbution.r --i  --o  
	
Usage:
	--input     input histo file
	--output    outdir
	--kmer	 kmer length
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input) ) { print_usage(spec) }
if ( is.null(opt$output) ) { print_usage(spec) }
if ( is.null(opt$kmer) ) { print_usage(spec) }
if(dir.exists(opt$output)){dir.create(opt$output)} 
library(findGSE)

findGSE(histo=opt$input,sizek=opt$kmer, outdir=opt$output)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)

