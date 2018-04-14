library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'mark','m',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
	'num','n',1,'character',
	'pop','p',1,'character',
	'bc','b',1,'character',
	'f','f',1,'character',
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
	--bc	bc gen for bcsft
	--f		f gen for bcsft
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
if(opt$pop =="bcsft" & is.null(opt$bc) & is.null(opt$f)){print_usage(spec)}
if(opt$pop =="bcsft") {	
	d<-read.cross(file=opt$mark,phefile=opt$trt,format="csvsr",crosstype=opt$pop,na.string="NaN",BC.gen=opt$bc,F.gen=opt$f)
}else{
	d<-read.cross(file=opt$mark,phefile=opt$trt,format="csvsr",crosstype=opt$pop,na.string="NaN")
}
if(!dir.exists(opt$out)){dir.create(opt$out)}
setwd(opt$out);
d<-jittermap(d)
d<-sim.geno(d)
d<-calc.genoprob(d)
phe.name<-colnames(d$pheno);
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
print(length(phe.name))
for(i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype"){next;}
	eff<-effectscan(d,pheno.col=phe.name[i],draw=FALSE);
	scan<-scanone(d,pheno.col=i,model="binary");
	scan.pm<-scanone(d,pheno.col=i,model="binary",n.perm=1000);
	markerid<-find.marker(d,chr=eff$chr,pos=eff$pos)
	outd<-data.frame(markerid=markerid,chr=scan$chr,pos=scan$pos,lod=scan$lod,eff=eff$a);
	write.table(file=paste(phe.name[i],".scan.csv",sep=""),sep="\t",outd,row.names=FALSE)
	write.table(file=paste(phe.name[i],".pm.csv",sep=""),sep="\t",scan.pm);
	scan.result<-summary(scan, perms=scan.pm, pvalues=TRUE)
	pm.result<-summary(scan.pm,alpha=c(0.01,0.05,0.1))
	legend=paste(rownames(pm.result),round(pm.result,2))
	threshold=pm.result[1,1]
	scan.result<-summary(scan,format="tabByCol",threshold=pm.result[1,1],drop=1)
	if(length(scan.result$lod$chr) < 1){
		scan.result<-summary(scan,format="tabByCol",threshold=pm.result[2,1],drop=1)
		threshold=pm.result[2,1]
		if(length(scan.result$lod$chr) <1){
			scan.result<-summary(scan,format="tabByCol",threshold=pm.result[3,1],drop=1)
			threshold=pm.result[3,1]
			if(length(scan.result$lod$chr) <1){
				pm.result<-c(3,2.5)
				scan.result<-summary(scan,format="tabByCol",threshold=3,drop=1);
				threshold=3
				if(length(scan.result$lod$chr) <1){
					threshold=2.5
					scan.result<-summary(scan,format="tabByCol",threshold=2.5,drop=1);
				}
				legend=pm.result
			}else{
				pm.result<-summary(scan.pm,alpha=c(0.05,0.1))
				legend=paste(rownames(pm.result),round(pm.result,2))
			}
		}else{
			pm.result<-summary(scan.pm,alpha=c(0.01,0.05))
			legend=paste(rownames(pm.result),round(pm.result,2))
		}
	}else{
		pm.result<-summary(scan.pm,alpha=c(0.01,0.05));
		legend=paste(rownames(pm.result),round(pm.result,2))
	}	pdf(file=paste(phe.name[i],".scan.pdf",sep=""))
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()
	png(file=paste(phe.name[i],".scan.png",sep=""))
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()
	qdata<-NULL
	n=0;
	for (j in chr){
		subd=which(outd$chr==j & outd$lod > theshold[1])
		if(length(subd) < 1){next;}
		start=1000;
		end=-1;
		for(k in c(2:length(subd))){
			if(subd[k]-subd[k-1] < 2){
				if(subd[k-1] < start){start=subd[k-1]}
				if(subd[k] > end){end=subd[k]}
			}else{	
				if (!is.null(qdata)){
					qdata<-rbind(qdata,data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end]))
				}else{
					qdata<-data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end])
				}
				n=n+1
				start=subd[k]
				end=subd[k]
			}
		}
		if(start != 1000){
			n=n+1;
			if (!is.null(qdata)){
				qdata<-rbind(qdata,data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end]))
			}else{
				qdata<-data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end])
			}
		}
	}
	qtlname=paste(phe.name[i],c(1:length(qdata$n)))
	qtl<-makeqtl(d,chr=qdata$chr,pos=qdata$pos,qtl.name=qtlname)
	if(opt$pop =="bcsft") {	
		fitqtl<-fitqtl(cross=d,qtl=qtl,pheno.col=i)
	}else{
		fitqtl<-fitqtl(cross=d,qtl=qtl,pheno.col=i,get.est=TRUE)
	}
	markerid<-find.marker(d,chr=qtl$chr,pos=qtl$pos)
	var<-fitqtl$result.drop[,"%var"]
	if (length(qtl$name) == 1){var<-fitqtl$result.full["Model","%var"]}
	data<-data.frame(marker=markerid,chr=qdata$chr,pos=qdata$pos,lod=qdata$lod,var=var,pm1=pm.result[1],pm2=pm.result[2],start=outd$pos[start],end=outd$pos[end])
	for(j in 1:length(qtlname)){
		data$mark1[j]=find.marker(d,chr=qtl$chr[j],data$start[j])
		data$mark2[j]=find.marker(d,chr=qtl$chr[j],data$end[j])
		pdf(paste(phe.name[i],".",qtl$name[j],".PXG.pdf",sep=""))
		plotPXG(d,data$marker[j],pheno.col=i)
		dev.off()
		png(paste(phe.name[i],".",qtl$name[i],".PXG.png",sep=""))
		plotPXG(d,data$marker[j],pheno.col=i)
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
