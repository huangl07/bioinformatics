#!/usr/bin/env Rscript
# load library
times<-Sys.time()
library(ape)
library('getopt');
options(bitmapType='cairo')

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'group', 'g', 1 , "character",
	'outgroup','u',1,"character",
	'raxml', 'raxml' , 1 , "character"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
Options: 
	--help		NULL 		get this help
	--infile 	character 	the input file [forced]
	--outfile 	character 	the filename for output graph [forced]
	--raxml	character	the raxml software for draw
	--outgroup	character	the outgroup sample split by ,
	\n")
	q(status=1);
}

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
#if ( is.null(opt$group) )	{ print_usage(spec) }
if ( is.null(opt$group) ){print("none colour tree!");}
library(ggtree)
if (is.null(opt$raxml)){
	raxml<-read.tree(file=opt$infile)
}else{
	raxml<-read.raxml(file=opt$infile)
}

if (!is.null(opt$outgroup)){
	raxml<-root(raxml,which(raxml$tip.label %in% as.list(strsplit(opt$outgroup, split=","))[[1]]))
}
write.tree(raxml,file="opt$outfile")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()