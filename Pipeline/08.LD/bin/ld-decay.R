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
ldfile<-read.table(opt$infile,head=TRUE)
files=ldfile$file
popid=ldfile$popid
col<-rainbow(5)
pop.id<-paste("pop",c(1:5))
for (i in 1:5){
	ld=read.table(file=as.character(files[i,1]),head=TRUE,comment.char=":");
	distance=ld$X.Dist
	distance=ld$X.Dist
	R2=ld$Mean_r.2
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
	newLD$distance<-newLD$distance/1000
	if (i != 1){
		lines(newLD,col=col[i],lwd=2)
	}else{
		plot(newLD,type="l",lwd=2,ylim=c(0,1),col=col[i],main="LD decay",ylab="R^2",xlab="distance(kb)")
	}
}
legend("topright",col=col,legend=pop.id,pch=1,cex=1)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
