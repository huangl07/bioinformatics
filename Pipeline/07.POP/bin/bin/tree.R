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
	--group		character	the group file for draw
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

if (!is.null(opt$group)){
	gro=read.table(opt$group,header=FALSE);
	col=unique(gro$V2);
	cls<-NULL
	for (i in 1:length(col)){
		cls<-c(cls,list(as.character(gro$V1[gro$V2==col[i]])));
	}
	names(cls)=paste("gro",1:length(col),sep="")
	print(col)
	raxml <- groupOTU(raxml,cls)
}
g<-NULL;
if (!is.null(opt$outgroup)){
	raxml<-root(raxml,which(raxml$tip.label %in% as.list(strsplit(opt$outgroup, split=","))[[1]]))
}




g<-ggtree(raxml,size=.2,layout="rectangular")+geom_tiplab(size=1,align=TRUE,linesize=.2)
if (!is.null(opt$raxml)){
g<-g+geom_text2(size=1,color="black",aes(subset=!isTip, label=bootstrap,color="black"))
}
if(!is.null(opt$group)){
g<-g+aes(color=group)+scale_color_manual(values=c(rainbow(length(col)+1)))
}

pdf(paste(opt$outfile,".rectangular.tree.pdf",sep=""))
print(g)
dev.off()

png(paste(opt$outfile,".rectangular.tree.png",sep=""))
print(g)
dev.off()
g<-NULL;

g<-ggtree(raxml,size=.2,layout="circular")+geom_tiplab2(size=2,align=TRUE,linesize=.2,aes(angle=angle))
if (!is.null(opt$raxml)){
g<-g+geom_text2(size=2,color="black",aes(subset=!isTip, label=bootstrap,color="black"))
}
if(!is.null(opt$group)){
g<-g+aes(color=group)+scale_color_manual(values=c(rainbow(length(col)+1)))
}

pdf(paste(opt$outfile,".circular.tree.pdf",sep=""))
print(g)

dev.off()
png(paste(opt$outfile,".circular.tree.png",sep=""))
print(g)

dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()