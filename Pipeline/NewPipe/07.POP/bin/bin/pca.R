#!/usr/bin/env Rscript 
times<-Sys.time()

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
	'varfile', 'v', 1 , "character"
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
	--varfile	character	the val file
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
if ( is.null(opt$varfile) )	{ print_usage(spec) }


pca.file<-read.table(opt$infile,header=TRUE)
if (!is.null(opt$group)){
	pop.list<-read.table(opt$group,header=TRUE)
	names(pop.list)=c("id","popid")
	pop.id<-as.vector(unique(pop.list$popid));
	color.list<-rainbow(length(pop.id));
	for (i in 1:length(pop.id)){
		pop.list$colour[pop.list$popid==pop.id[i]]=color.list[i]
	}
}else{
	pop.id<-c(1);
	color.list<-rainbow(1);
}
pdf(paste(opt$outfile,".pc1vspc2.pdf",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC2,col=pop.list$colour,xlab="PC1",ylab="PC2",main="PC1 vs PC2")
	legend("bottom",col=color.list,legend=pop.id,bty=T,pch=1)
}else{
	plot(pca.file$PC1,pca.file$PC2,xlab="PC1",ylab="PC2",main="PC1 vs PC2")
}
dev.off();
png(paste(opt$outfile,".pc1vspc2.png",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC2,col=pop.list$colour,xlab="PC1",ylab="PC2",main="PC1 vs PC2")
	legend("bottom",col=color.list,legend=pop.id,bty=T,pch=1)
}else{
	plot(pca.file$PC1,pca.file$PC2,xlab="PC1",ylab="PC2",main="PC1 vs PC2")
}
dev.off();
pdf(paste(opt$outfile,".pc1vspc3.pdf",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC3,col=pop.list$colour,xlab="PC1",ylab="PC2",main="PC1 vs PC2")
	legend("bottom",col=color.list,legend=pop.id,bty=T,pch=1)
}else{
	plot(pca.file$PC1,pca.file$PC3,xlab="PC1",ylab="PC2",main="PC1 vs PC3")
}
dev.off();
png(paste(opt$outfile,".pc1vspc3.png",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC1,pca.file$PC3,col=pop.list$colour,xlab="PC1",ylab="PC3",main="PC1 vs PC3")
	legend("bottom",col=color.list,legend=pop.id,bty=T,pch=1)
}else{
	plot(pca.file$PC1,pca.file$PC3,xlab="PC1",ylab="PC2",main="PC1 vs PC3")
}
dev.off();
pdf(paste(opt$outfile,".pc2vspc3.pdf",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC2,pca.file$PC3,col=pop.list$colour,xlab="PC3",ylab="PC3",main="PC2 vs PC3")
	legend("bottom",col=color.list,legend=pop.id,bty=T,pch=1)
}else{
	plot(pca.file$PC2,pca.file$PC3,xlab="PC1",ylab="PC2",main="PC2 vs PC3")
}
png(paste(opt$outfile,".pc2vspc3.png",sep=""))
if (length(pop.id) >1){
	plot(pca.file$PC2,pca.file$PC3,col=pop.list$colour,xlab="PC2",ylab="PC3",main="PC2 vs PC3")
	legend("bottom",col=color.list,legend=pop.id,bty=T,pch=1)
}else{
	plot(pca.file$PC2,pca.file$PC3,xlab="PC1",ylab="PC2",main="PC2 vs PC3")
}
dev.off();

value<-read.table(opt$varfile,header=FALSE);
sum<-sum(value$V1)
value$V1<-value$V1/sum*100
pdf(paste(opt$outfile,".val.pdf",sep=""),width=800,height=800);
valmax=length(value$V1);
barplot(value$V1,col="blue",xlab="PCAs",ylab="Variance%",beside=FALSE,names.arg=c(1:valmax),border=TRUE,ylim=c(0,100),main="Variance of 
PCAs")
dev.off()
png(paste(opt$outfile,".val.png",sep=""),width=800,height=800);
barplot(value$V1,col="blue",xlab="PCAs",ylab="Variance%",beside=FALSE,names.arg=c(1:valmax),border=TRUE,ylim=c(0,100),main="Variance of 
PCAs")
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
