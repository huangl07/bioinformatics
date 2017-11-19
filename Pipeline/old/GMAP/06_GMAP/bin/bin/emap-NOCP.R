library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'mark','m',1,'character',
	'out','o',1,'character',
	'pop','p',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript emap.R --mark --out --pop
	
Usage:
	--mark	map file
	--out	out dir
	--pop	pop type
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) ) { print_usage(spec) }
if ( is.null(opt$pop) ) { print_usage(spec) }

if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}

d<-read.cross(file=opt$mark,format="csvr",crosstype=opt$pop)
setwd(opt$out)
d<-jittermap(d)
d<-est.rf(d)
summary<-summaryMap(d)
write.table(file="summary.table",summary)

pdf("total.map.pdf")
plotMap(d,alternate.chrid=TRUE,shift=TRUE)
dev.off()

pdf("total.rf.pdf")
plotRF(d,col.scheme="redblue",alternate.chrid=TRUE,mark.diagonal=TRUE)
dev.off()
chrname<-chrnames(d);
for(i in chrname){
	pdf(paste(i,".rf.pdf",sep=""))
	plotRF(d,chr=i,col.scheme="redblue",alternate.chrid=TRUE)
	dev.off()
}
pdf("total.geno.pdf");
plotGeno(d,include.xo=TRUE,horizontal=TRUE,min.sep=2)
dev.off()
for(i in chrname){
	pdf(paste(i,".geno.pdf",sep=""))
	plotGeno(d,chr=i,include.xo=TRUE,horizontal=TRUE,min.sep=2)
	dev.off()
}
pdf("total.miss.pdf");
par(mfrow=c(1,2), las=1)
plot(ntyped(d), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(d, "mar"), ylab="No. typed individuals",
main="No. genotypes by marker")
dev.off()
pdf("total.dup.pdf")
cg <- comparegeno(d)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
dev.off()

escaptime=Sys.time()-times;
print("Done!\n")
print(escaptime)