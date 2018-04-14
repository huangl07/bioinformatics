library('getopt');
options(bitmapType='cairo')

spec = matrix(c(
	'hapmap','h',0,'character',
	'trait','t',0,'character',
	'output','o',0,'character',
	'CV','c',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--hapmap	the input hapmap file
	--trait	the trait file 
	--output	the output dir
	--CV	input CV file
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$hapmap)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

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

source("/mnt/ilustre/users/dna/Environment/biotools/GAPIT/emma.txt");
source("/mnt/ilustre/users/dna/Environment/biotools/GAPIT/gapit_functions.txt")

myY <- read.table(opt$trait, head = TRUE,na.strings="NaN")
myG <- read.table(opt$hapmap, head = FALSE)
if(!dir.exists(opt$output)){dir.create(opt$output)}
if (is.null(opt$CV)){
	setwd(opt$output)
	myGAPIT <- GAPIT(Y=myY,G=myG,PCA.total=3,Prior="memory")
}else{
	myCV <- read.table(opt$CV, head = FALSE)
	names(myCV)[1]="Taxa";
	for(i in 2:length(colnames(myCV))){names(myCV)[i]=paste("Q",i-1,sep="");}
	setwd(opt$output)
	myGAPIT<-GAPIT(Y=myY,G=myG,CV=myCV,Prior="memory")
}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
