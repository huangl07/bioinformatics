library('getopt');

spec = matrix(c(
	'hapmap','h',0,'character',
	'trait','t',0,'character',
	'output','o',0,'character',
	'help','c',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript indel_len.r --i  --o  
	
Usage:
	--hapmap	the input hapmap file
	--trait	the trait file 
	--output	the output dir
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$hapmap)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }

library(multtest)
library(gplots)
library("grid")
library(LDheatmap)
library(genetics)
library(ape)
library(Matrix)
library(EMMREML)
library(compiler) 
library("scatterplot3d")

source("/mnt/ilustre/users/long.huang/Software/GAPIT/emma.txt");
source("/mnt/ilustre/users/long.huang/Software/GAPIT/gapit_functions.txt")

myY <- read.table(opt$trait, head = TRUE)
myG <- read.table(opt$hapmap, head = FALSE)
setwd(opt$output);
myGAPIT <- GAPIT(Y=myY,G=myG,PCA.total=3,Geno.View.output=FALSE,)
