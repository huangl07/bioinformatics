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
pdf("total.lg.pdf");
plot(d)
dev.off()
pdf("total.lg.png");
plot(d)
dev.off()

pdf("total.rf.pdf")
plotRF(d,col.scheme="redblue",alternate.chrid=TRUE,mark.diagonal=TRUE,what="rf")
dev.off()
png("total.rf.png")
plotRF(d,col.scheme="redblue",alternate.chrid=TRUE,mark.diagonal=TRUE,what="rf")
dev.off()

chrname<-chrnames(d);
for(i in chrname){
	pdf(paste(i,".rf.pdf",sep=""))
	plotRF(d,chr=i,col.scheme="redblue",alternate.chrid=TRUE,what="rf")
	dev.off()
	png(paste(i,".rf.png",sep=""))
	plotRF(d,chr=i,col.scheme="redblue",alternate.chrid=TRUE,what="rf")
	dev.off()
}
png("total.maker.png")
profileMark(d, stat.type = c("seg.dist", "miss", "dxo"),layout = c(1,3),type="l",cex=0.7,byChr=TRUE)
dev.off()
pdf("total.marker.pdf")
profileMark(d, stat.type = c("seg.dist", "miss", "dxo"),layout = c(1,3),type="l",cex=0.7,byChr=TRUE)
dev.off()

escaptime=Sys.time()-times;
print("Done!\n")
print(escaptime)
