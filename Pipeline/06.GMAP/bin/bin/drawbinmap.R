times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'mark','m',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
	'num','n',1,'character',
	'pop','p',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript Rqtl-NOCP.r --mark  --out --num --pop
	
Usage:
	--mark	map file
	--out	out dir
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

d<-read.table(opt$mark,head=TRUE,sep=",")
d<-na.omit(d)
cols=c("green", "blue", "red","grey40")
colnames(d)[2:3]=c("chr","pos")
chr=unique(d$chr);
for(i in chr){
	pos=d[d$chr==i]
	d$start[d$chr==i]=(d$pos[c(1,1:(length(pos)-1))]+d$pos[c(1:length(pos))])/2
	d$end[d$chr==i]=(d$pos[c(2:(length(pos)),length(pos))]+d$pos[c(1:length(pos))])/2
}
chr.len <- tapply(d$end, d$chr, max)
x.brks <- cumsum(chr.len) - chr.len/2
chr.cum.len <- c(0, cumsum(chr.len)[-length(chr.len)])
names(chr.cum.len) <- names(chr.len)
d$start <- d$start + chr.cum.len[d$chr]
d$end <- d$end + chr.cum.len[d$chr]
names(cols) <- c("A", "B","H", "-")

pdf(paste(opt$out,"pdf",sep="."));
plot(5, type = "n", xlim = c(0, max(d$end)), ylim = c(0,ncol(d) - 4), xaxt = "n")
axis(side = 1, at = x.brks, labels = unique(d$chr))
for (i in 5:ncol(d)) {
    rect(xleft = d$start, ybottom = i - 5, xright = d$end,ytop = i - 4, col = cols[d[, i]], border = NA)
}
abline(v = c(0, cumsum(chr.len)), col = "grey80", lwd = 0.5)
dev.off()

png(paste(opt$out,"png",sep="."));
plot(5, type = "n", xlim = c(0, max(d$end)), ylim = c(0,ncol(d) - 4), xaxt = "n")
axis(side = 1, at = x.brks, labels = unique(d$chr))
for (i in 5:ncol(d)) {
    rect(xleft = d$start, ybottom = i - 5, xright = d$end,ytop = i - 4, col = cols[d[, i]], border = NA)
}
abline(v = c(0, cumsum(chr.len)), col = "grey80", lwd = 0.5)
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
