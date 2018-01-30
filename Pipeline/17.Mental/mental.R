library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'dif','d',0,'character',
	'env','e',0,'character',
	'out','o',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--dif	the dif matrix file
	--env	the env file 
	--out	the output file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$dif)){ print_usage(spec) }
if ( is.null(opt$env)){ print_usage(spec) }
if ( is.null(opt$out)){ print_usage(spec) }

times<-Sys.time()
library(ape)
library(ade4)
library(vegan)
d1<-read.table(opt$dif,head=FALSE);
d2<-read.table(opt$env,head=TRUE);
d1<-d1[order(d1[,1],decreasing=T),];
d2<-d2[order(d2[,1],decreasing=T),];
if(!dir.exists(opt$out)){dir.create(opt$out)}
setwd(opt$out);
d2.names<-colnames(d2)
data<-NULL;
for (i in 2:length(d2.names)){
	matrix2<-as.matrix(dist(d2[[i]]));
	mantel<-mantel.test(d1[,-1], matrix2, graph = FALSE)
	vegan<-mantel(xdis = d1[-1], ydis = matrix2)
	if (!is.null(data)){
		data<-rbind(data,data.frame(names=d2.names[i],z.stat=mantel$z.stat,p=mantel$p,vegan<-vegan$statistic));
	}else{
		data<-data.frame(names=d2.names[i],z.stat=mantel$z.stat,p=mantel$p,vegan<-vegan$statistic)
	}
	pdf(paste(d2.names[i],".pdf",sep=""));
	plot(unlist(d1[,-1]),as.vector(unlist(matrix2)),main="Mantel test",ylab=d2.names[i],xlab="genetic distance");
	d<-data.frame(x=unlist(d1[,-1]),y=as.vector(unlist(matrix2)))
	d<-d[order(d[,1],decreasing=T),]
	lm.sol<-lm(d$y~d$x)
	abline(lm.sol,col="blue")
	dev.off();
	png(paste(d2.names[i],".png",sep=""));
	plot(unlist(d1[,-1]),as.vector(unlist(matrix2)),main="Mantel test",ylab=d2.names[i],xlab="genetic distance");
	lm.sol<-lm(d$y~d$x)
	abline(lm.sol,col="blue")
	dev.off();
}
write.table(file="mantel.test.xls",data,row.names=FALSE);
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
