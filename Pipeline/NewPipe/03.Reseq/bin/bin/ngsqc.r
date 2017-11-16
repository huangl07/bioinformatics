#!/usr/bin/env Rscript
times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'base','b',0,'character',
	'qual','q',0,'character',
	'out','o',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript QC.r --base  --qual  --key  --od
	
Usage:
	--base	base distribution file
	--qual	base quality file
	--key	sample ID
	--od	output dir
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$base) ) { print_usage(spec) }
if ( is.null(opt$qual) ) { print_usage(spec) }
if ( is.null(opt$out) ) { print_usage(spec) }

read1<-read.table(opt$base,header=FALSE)
read2<-read.table(opt$qual,header=FALSE)

read1<-read1[complete.cases(read1),]
read2<-read2[complete.cases(read2),]

num<-(length(read2[,1])/2)
for (i in 1:(length(read2[,1])-1)) {if ((as.numeric(read2[i+1,1])-as.numeric(read2[i,1]))!=1) {num<-i}}

quan<-read2[,ncol(read2)]
quan<-10^(quan/10*-1)
max_quan=ceiling(max(quan)/0.001)*0.001
Aper<-read1[,2]
Cper<-read1[,5]
Gper<-read1[,4]
Tper<-read1[,3]
Nper<-read1[,6]
max_per=ceiling(mean(c(Aper,Cper,Gper,Tper,Nper))/0.25)*0.25*2
axix_at<-c(1,50,100,num,num+50,num+100,length(read2[,1]))
axix_txt<-c(1,50,100,num,50,100,length(read2[,1])/2)


pdf(file=paste(opt$od,"/",opt$out,".qual.pdf",sep=""))
barplot(quan*100,col='springgreen2',space=0,ylab='Error rate(%)',border= NA,ylim=c(0,max_quan*100),xlab='Reads position(bp)',main='Quality distribution')
axis(1,labels=axix_txt,at=axix_at)
abline(v=num+0.5, lty=2,col='darkgray')
box()
dev.off()
png(file=paste(opt$od,"/",opt$out,".qual.png",sep=""))
barplot(quan*100,col='springgreen2',space=0,ylab='Error rate(%)',border= NA,ylim=c(0,max_quan*100),xlab='Reads position(bp)',main='Quality distribution')
axis(1,labels=axix_txt,at=axix_at)
abline(v=num+0.5, lty=2,col='darkgray')
box()
dev.off()

pdf(file=paste(opt$od,"/",opt$out,".base.pdf",sep=""))
plot(Aper,col='red',type = 'l',xlab='Reads position(bp)',ylab='Percent(%)',ylim=c(0,max_per),xaxt="n",lty=1,lwd=1.5,main="Base distribution")
axis(1,labels=axix_txt,at=axix_at)
abline(v=num+0.5, lty=2,col='darkgray')
lines(Tper,col='orange',lty=1,lwd=1.5)
lines(Gper,col='green',lty=1,lwd=1.5)
lines(Cper,col='blue',lty=1,lwd=1.5)
lines(Nper,col='black',lty=1,lwd=1.5)
legend("topright",c("A","T","G","C","N"),lty=c(1,1,1,1,1),col=c("red","orange","green","blue","black"))
dev.off()

png(file=paste(opt$od,"/",opt$out,".base.png",sep=""))
plot(Aper,col='red',type = 'l',xlab='Reads position(bp)',ylab='Percent(%)',ylim=c(0,max_per),xaxt="n",lty=1,lwd=1.5,main="Base distribution")
axis(1,labels=axix_txt,at=axix_at)
abline(v=num+0.5, lty=2,col='darkgray')
lines(Tper,col='orange',lty=1,lwd=1.5)
lines(Gper,col='green',lty=1,lwd=1.5)
lines(Cper,col='blue',lty=1,lwd=1.5)
lines(Nper,col='black',lty=1,lwd=1.5)
legend("topright",c("A","T","G","C","N"),lty=c(1,1,1,1,1),col=c("red","orange","green","blue","black"))
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
