#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}
if ( !is.null(opt$help)) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outfile)){ print_usage(spec) }

times<-Sys.time()
ld<-read.table(opt$infile,head=TRUE,na.strings=c("nan","-nan"))
distance=ld$V1
R2=ld$V2
n=length(R2)
HW.st<-c(C=0.01)
HW.nonlinear<-nls(R2~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
newLD<-data.frame(distance,fpoints)
maxld<-max(newLD$fpoints)
decay05<-newLD$distance[which.min(abs(newLD$fpoints-maxld/2))]
decay01<-newLD$distance[which.min(abs(newLD$fpoints-0.1))]
newLD<-newLD[order(newLD$distance),]
pdf(paste(opt$outfile,"pdf",sep="."))
plot(distance,R2,type="p",col="blue",main="LD decay distribution",ylab="R square",xlab="Distance",sub=paste("decay05:",decay05,"decay01:",decay01,sep=" "));
lines(newLD,col="red",lwd=2)
dev.off()

png(paste(opt$outfile,"png",sep="."))
plot(distance,R2,type="p",col="blue",main="LD decay distribution",ylab="R square",xlab="Distance",sub=paste("decay05:",decay05,"decay01:",decay01,sep=" "));
lines(newLD,col="red",lwd=2)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
