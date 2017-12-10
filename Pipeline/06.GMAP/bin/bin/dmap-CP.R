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
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) ) { print_usage(spec) }
if ( is.null(opt$pop) ) { print_usage(spec) }

if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}
library('qtl');
library('ASMap');
opt$pop=tolower(opt$pop)
if(opt$pop == "cp"){
	d<-read.cross(genfile=paste(opt$mark,"loc",sep="."),phefile=paste(opt$mark,"trt",sep="."),mapfile=paste(opt$mark,"sexAver.map",sep="."),format="mapqtl")
	d1<-read.cross(genfile=paste(opt$mark,"loc",sep="."),phefile=paste(opt$mark,"trt",sep="."),mapfile=paste(opt$mark,"male.map",sep="."),format="mapqtl")
	d2<-read.cross(genfile=paste(opt$mark,"loc",sep="."),phefile=paste(opt$mark,"trt",sep="."),mapfile=paste(opt$mark,"female.map",sep="."),format="mapqtl")
	setwd(opt$out)
	d<-jittermap(d)
	d<-est.rf(d)
	d1<-jittermap(d1)
	d1<-est.rf(d1)
	d2<-jittermap(d2)
	d2<-est.rf(d2)
	pdf("total.sexAver.pdf");
	plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	png("total.sexAver.png");
	plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	pdf("total.male.pdf");
	plotMap(d1,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	png("total.male.png");
	plotMap(d1,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	pdf("total.female.pdf");
	plotMap(d2,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	png("total.female.png");
	plotMap(d2,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	pdf("heatMap.sexAver.pdf")
	heatMap(d,lmax=50)
	dev.off()
	png("heatMap.sexAver.png")
	heatMap(d,lmax=50)
	dev.off()
	pdf("heatMap.male.pdf")
	heatMap(d1,lmax=50)
	dev.off()
	png("heatMap.male.png")
	heatMap(d1,lmax=50)
	dev.off()
	pdf("heatMap.female.pdf")
	heatMap(d2,lmax=50)
	dev.off()
	png("heatMap.female.png")
	heatMap(d2,lmax=50)
	dev.off()
	for(i in chrname){
		pdf(paste(i,".heatMap.sexAver.pdf",sep=""))
		heatMap(d,chr=i,lmax=50)
		dev.off()
		png(paste(i,".heatMap.sexAver.png",sep=""))
		heatMap(d,chr=i,lmax=50)
		dev.off()
	}
	for(i in chrname){
		pdf(paste(i,".heatMap.male.pdf",sep=""))
		heatMap(d1,chr=i,lmax=50)
		dev.off()
		png(paste(i,".heatMap.male.png",sep=""))
		heatMap(d1,chr=i,lmax=50)
		dev.off()
	}
	for(i in chrname){
		pdf(paste(i,".heatMap.female.pdf",sep=""))
		heatMap(d2,chr=i,lmax=50)
		dev.off()
		png(paste(i,".heatMap.female.png",sep=""))
		heatMap(d2,chr=i,lmax=50)
		dev.off()
	}
}else{
	d<-read.cross(file=paste(opt$mark,"csv",sep="."),format="csvsr",crosstype=opt$pop)
	setwd(opt$out)
	d<-jittermap(d)
	d<-est.rf(d)
	pdf("total.lg.pdf");
	plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	png("total.lg.png");
	plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	dev.off()
	pdf("total.rf.pdf")
	heatMap(d,lmax=50)
	dev.off()
	png("total.rf.png")
	heatMap(d,lmax=50)
	dev.off()
	for(i in chrname){
		pdf(paste(i,".rf.pdf",sep=""))
		heatMap(d,chr=i,lmax=50)
		dev.off()
		png(paste(i,".rf.png",sep=""))
		heatMap(d,chr=i,lmax=50)
		dev.off()
	}
}


escaptime=Sys.time()-times;
print("Done!\n")
print(escaptime)
