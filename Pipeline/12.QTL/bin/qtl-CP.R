library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'map','m',1,'character',
	'loc','l',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
	'num','n',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript Rqtl_CP.r --map --loc --trt --out --num
	
Usage:
	--map	map file
	--loc	loc file
	--trt	trt file
	--out	out dir
	--num	pm number
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$map) ) { print_usage(spec) }
if ( is.null(opt$loc) ) { print_usage(spec) }
if ( is.null(opt$trt) ) { print_usage(spec) }
if ( is.null(opt$num) ) { opt$num=1000; }
if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}


d<-read.cross(mapfile=opt$map,genfile=opt$loc,phefile=opt$trt,format="mapqtl",crosstype="4way")
setwd(opt$out);
d<-jittermap(d)
d<-sim.geno(d)
d<-calc.genoprob(d)
phe.name<-colnames(d$pheno)
nrow=4;
ncol=ceiling(length(phe.name)/4)
if(ncol==1){nrow=length(phe.name)}
ncol=ceiling(sqrt(length(phe.name)));
nrow=ncol;
pdf("pheno.pdf",width=30*ncol,height=40*nrow)
par(mfrow=c(ncol,nrow))
for (i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype"){next;}
	plotPheno(d,pheno.col=phe.name[i])
}
dev.off()
png("pheno.png")
par(mfrow=c(ncol,nrow))
for (i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype"){next;}
	plotPheno(d,pheno.col=phe.name[i])
}
dev.off()
for(i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype" | phe.name[i]=="sampleID"){next;}
	print(paste(opt$method,"trait",phe.name[i],sep="\t"))
	scan<-scanone(d,pheno.col=phe.name[i]);
	scan.pm<-scanone(d,pheno.col=phe.name[i],n.perm=opt$num);
	markerid<-find.marker(d,chr=scan$chr,pos=scan$pos)
	outd<-data.frame(markerid=markerid,chr=scan$chr,pos=scan$pos,lod=scan$lod);
	write.table(file=paste(phe.name[i],".scan.csv",sep=""),sep="\t",outd,row.names=FALSE)
	write.table(file=paste(phe.name[i],".pm.csv",sep=""),sep="\t",scan.pm);
	scan.result<-summary(scan, perms=scan.pm, pvalues=TRUE)
	if(min(scan.result$pval) >0.1){
		scan.result<-summary(scan,format="tabByCol",threshold=3,drop=1)
		if(length(rownames(scan.result$lod)) < 1){
			scan.result<-summary(scan,format="tabByCol",threshold=2.5,drop=1)
		}
		pm.result<-c(3,2.5)
		legend=pm.result
	}else{	
		scan.result<-summary(scan,format="tabByCol",perms=scan.pm,alpha=0.1,drop=1)
		if(length(rownames(scan.result$lod)) < 1){
			scan.result<-summary(scan,format="tabByCol",alpha=0.05,drop=1)
		}

		pm.result<-summary(scan.pm,alpha=c(0.01,0.05))
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
	if(length(rownames(scan.result$lod)) < 1){
		next;
	}
	qtlname=paste(phe.name[i],c(1:length(rownames(scan.result$lod))))
	qtl<-makeqtl(d,chr=scan.result$lod$chr,pos=scan.result$lod$pos,qtl.name=qtlname)
	fitqtl<-fitqtl(cross=d,qtl=qtl,get.est=TRUE,pheno.col=i)
	markerid<-find.marker(d,chr=qtl$chr,pos=qtl$pos)
	var<-fitqtl$result.drop[,"%var"]
	if (length(qtl$name) == 1){var<-fitqtl$result.full["Model","%var"]}
	data<-data.frame(marker=markerid,chr=scan.result$lod$chr,pos=scan.result$lod$pos,lod=scan.result$lod$lod,var=var,pm1=pm.result[1],pm2=pm.result[2])
	for(j in 1:length(qtlname)){
		insert<-bayesint(scan,chr=qtl$chr[j],expandtomarkers=FALSE,prob=0.99)
		data$start[j]=min(insert$pos);
		data$end[j]=max(insert$pos);
		data$mark1[j]=find.marker(d,chr=qtl$chr[j],data$start[j])
		data$mark2[j]=find.marker(d,chr=qtl$chr[j],data$end[j])
		pdf(paste(phe.name[i],".",qtl$name[j],".PXG.pdf",sep=""))
		plotPXG(d,data$marker[j],pheno.col=phe.name[i])
		dev.off()
		png(paste(phe.name[i],".",qtl$name[i],".PXG.png",sep=""))
		plotPXG(d,data$marker[j],pheno.col=phe.name[i])
		dev.off()
	}
	write.table(file=paste(phe.name[i],".qtl.csv",sep=""),sep="\t",data,row.names=FALSE)
	pdf(paste(phe.name[i],".qtl.pdf",sep=""))
	plot(qtl)
	dev.off()
	png(paste(phe.name[i],".qtl.png",sep=""))
	plot(qtl)
	dev.off()
}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)