#!/usr/bin/env Rscript
times<-Sys.time();
library('getopt');
options(bitmapType='cairo');

spec = matrix(c(
	'i','a',0,'character',
	'o','b',0,'character',
	's','s',0,'character',
	'help','c',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript base_distrbution.r --i  --o  
	
Usage:
	--i     base_distrbution
	--o	     base_distrbution picture name
	--s	single
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$i) ) { print_usage(spec) }
if ( is.null(opt$o) ) { print_usage(spec) }

xmax =2000;
a=read.table(opt$i)
b=cumsum(as.numeric(a[,2]));
sums=sum(as.numeric(a[,2]))
b=b/sums;
c=a[,2]/sums;
lmax=max(c);
library(plotrix)
xmean=a$V1[match(lmax,c)]
if (xmean < 10){xmax=20;}else{xmax=100}
if (is.null(opt$s)){
	pdf(paste(opt$o,".pdf",sep=""))

	twoord.plot(a[,1],c,a[,1],b,xlim=c(0,xmax),lylim=c(0,lmax*1.2),rylim=c(0,1.1),lcol="red3",rcol="blue",xlab="Sequencing depth",ylab="Percent of base",rylab="Percent of cumulative base",type=c("l","l"))
	legend("right",legend=c("Percent of base","Percent of cumulative base"),col=c("red3","blue"),lty=1,bty="n",text.col=c("red3","blue"),cex=0.5)
	dev.off()
	png(paste(opt$o,".png",sep=""))
	twoord.plot(a[,1],c,a[,1],b,xlim=c(0,xmax),lylim=c(0,lmax*1.2),rylim=c(0,1.1),lcol="red3",rcol="blue",xlab="Sequencing depth",ylab="Percent of base",rylab="Percent of cumulative base",type=c("l","l"))
	legend("right",legend=c("Percent of base","Percent of cumulative base"),col=c("red3","blue"),lty=1,bty="n",text.col=c("red3","blue"),cex=0.5)
	dev.off()
}else{
	pdf(paste(opt$o,".pdf",sep=""))
	plot(a[,1],c,col="blue",type="l",xlab="Sequencing depth",ylab="Percent of base",xlim=c(0,xmax))
	dev.off();
	png(paste(opt$o,".png",sep=""))
	plot(a[,1],c,col="blue",type="l",xlab="Sequencing depth",ylab="Percent of base",xlim=c(0,xmax))
	dev.off();
}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)

