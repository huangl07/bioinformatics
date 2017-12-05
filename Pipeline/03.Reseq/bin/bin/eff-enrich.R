#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--input	the input  file
	--output	the out file 
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

data<-read.table(opt$input,sep="\t",comment.char="^",head=TRUE)
names(data)<-c("id","description","k","M","n","N");
data<-na.omit(data)
pvalue<-phyper(data$k,data$M,data$N-data$M,data$n,lower.tail=FALSE);
qvalue<-p.adjust(pvalue,method="fdr")
outd<-data.frame(id=data$id,des=data$description,eff=data$k,total=data$n,pvalue=1-pvalue,qvalue=1-qvalue)
outd<-outd[order(outd$eff,decreasing=T),]
write.table(file=paste(opt$output,"detail",sep="."),outd,row.name=FALSE);
order<-order(outd$qvalue)[0:20]
draw<-data.frame(id=data$id,des=data$description,eff=data$k/data$n,total=data$M/data$N,pvalue=1-pvalue,qvalue=1-qvalue)
ymax=max(max(draw$eff[order]),max(draw$total[order]))

pdf(paste(opt$output,"pdf",sep="."))
barplot(rbind(draw$eff[order],draw$total[order]),col=c("red","blue"),beside=TRUE,names.arg=outd$id[order],las=2,ylim=c(0,ymax/0.8),legend=c("eff","total"),main=("Enrichment Analysis"))
box()
dev.off()
png(paste(opt$output,"png",sep="."))
barplot(rbind(draw$eff[order],draw$total[order]),col=c("red","blue"),beside=TRUE,names.arg=outd$id[order],las=2,ylim=c(0,ymax/0.8),legend=c("eff","total"),main=("Enrichment Analysis"))
box()
dev.off()



escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
