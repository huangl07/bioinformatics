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
	--trt	trt file
	--pop	pop type
	--out	out dir
	--num	pm number
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) ) { print_usage(spec) }
if ( is.null(opt$trt) ) { print_usage(spec) }
if ( is.null(opt$pop) ) { print_usage(spec) }
if ( is.null(opt$num) ) { opt$num=1000; }else{opt$num=as.numeric(opt$num)}
if ( is.null(opt$out) ) { opt$out="./";}

d<-read.cross(file=opt$mark,phefile=opt$trt,format="csvsr",crosstype=opt$pop,na.string="NaN")
if(!dir.exists(opt$out)){dir.create(opt$out)}
setwd(opt$out);
d<-jittermap(d)
d<-sim.geno(d)
d<-calc.genoprob(d)
phe.name<-colnames(d$pheno);
ncol=ceiling(sqrt(length(phe.name)));
nrow=ncol;
pdf(paste("pheno.pdf",sep="/"),width=30*ncol,height=40*nrow)
par(mfrow=c(ncol,nrow))
for (i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype"){next;}
	plotPheno(d,pheno.col=phe.name[i])
}
dev.off()
png(paste("pheno.png",sep="/"))
par(mfrow=c(ncol,nrow))
for (i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype"){next;}
	plotPheno(d,pheno.col=phe.name[i])
}
dev.off()
qtls<-matrix()
print(length(phe.name))
for(i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype"){next;}
	eff<-effectscan(d,pheno.col=phe.name[i],draw=FALSE);
	scan<-cim(d,pheno.col=phe.name[i]);
	scan.pm<-cim(d,pheno.col=phe.name[i],n.perm=opt$num);
	markerid<-find.marker(d,chr=eff$chr,pos=eff$pos)
	outd<-data.frame(markerid=markerid,chr=scan$chr,pos=scan$pos,lod=scan$lod,eff=eff$a);
	write.table(file=paste(phe.name[i],".scan.csv",sep=""),sep="\t",outd,row.names=FALSE)
	write.table(file=paste(phe.name[i],".pm.csv",sep=""),sep="\t",scan.pm);
	scan.result<-summary(scan, perms=scan.pm, pvalues=TRUE)
	if(min(scan.result$pval) >0.1){
		scan.result<-summary(scan,format="tabByCol",threshold=3,drop=1)
		pm.result<-c(3,2.5)
		legend=pm.result
	}else{	
		pm.result<-summary(scan.pm,alpha=c(0.01,0.05))
		scan.result<-summary(scan,format="tabByCol",perms=scan.pm,alpha=0.1,drop=1)
		legend=paste(rownames(pm.result),round(pm.result,2))
	}
	pdf(file=paste(phe.name[i],".scan.pdf",sep=""))
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()
	png(file=paste(phe.name[i],".scan.png",sep=""))
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()

	qtlname=paste(phe.name[i],c(1:length(scan.result$lod$chr)))
	qtl<-makeqtl(d,chr=scan.result$lod$chr,pos=scan.result$lod$pos,qtl.name=qtlname)
	fitqtl<-fitqtl(cross=d,qtl=qtl,get.est=TRUE,pheno.col=i)
	markerid<-find.marker(d,chr=qtl$chr,pos=qtl$pos)
	var<-fitqtl$result.drop[,"%var"]
	if (length(qtl) == 1){var<-fitqtl$result.full["Model","%var"]}
	data<-data.frame(marker=markerid,chr=scan.result$lod$chr,pos=scan.result$lod$pos,lod=scan.result$lod$lod,var=var,pm1=pm.result[1],pm2=pm.result[2])
	for(j in 1:length(qtlname)){
		insert<-bayesint(scan,chr=qtl$chr[j],expandtomarkers=FALSE,prob=0.99)
		data$start[j]=min(insert$pos);
		data$end[j]=max(insert$pos);
		data$mark1[j]=find.marker(d,chr=qtl$chr[j],data$start[j])
		data$mark2[j]=find.marker(d,chr=qtl$chr[j],data$end[j])
	}
	write.table(file=paste(phe.name[i],".qtl.csv",sep=""),sep="\t",data,row.names=FALSE)
	pdf(paste(phe.name[i],".qtl.pdf",sep=""))
	plot(qtl)
	dev.off()
	pdf(paste(phe.name[i],".PXG.pdf",sep=""))
	plotPXG(d,data$marker,pheno.col=phe.name[i])
	dev.off()
	png(paste(phe.name[i],".qtl.png",sep=""))
	plot(qtl)
	dev.off()
	png(paste(phe.name[i],".PXG.png",sep=""))
	plotPXG(d,data$marker,pheno.col=phe.name[i])
	dev.off()

}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)