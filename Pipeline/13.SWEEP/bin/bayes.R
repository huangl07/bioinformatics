#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'help','h',0,"logical",
	'infile','i',1,"character",
	'outfile','o',1,"character"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
Options: 
	--help		NULL 		get this help
	--infile 	character 	bayscan file
	--outfile 	character 	the filename for output graph [forced]
	")
	q(status=1);
}
if(is.null(opt$infile)){print_usage(spec)}
if(is.null(opt$outfile)){print_usage(spec)}

res=read.table(opt$infile);
FDR=0.05;
pos=0.35;
colfstat=5
colq=colfstat-2
size=1
highlight=NULL;name_highlighted=F;add_text=T;
highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
ok_outliers=TRUE
if (sum(res[,colq]<=FDR)==0)
	ok_outliers=FALSE;

res[res[,colq]<=0.0001,colq]=0.0001

pdf(paste(opt$outfile,"bayes.pdf",sep="."));
plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
if (name_highlighted) {
 	if (length(highlight_rows)>0) {
 		#text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
 	}
}else {
	points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
	# add names of loci over p and vertical line
	if (ok_outliers & add_text) {
		#text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
	}
}
lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
dev.off();
png(paste(opt$outfile,"bayes.png",sep="."));
plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
if (name_highlighted) {
 	if (length(highlight_rows)>0) {
 		text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
 	}
}else {
	points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
	# add names of loci over p and vertical line
	if (ok_outliers & add_text) {
		#text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
	}
}
lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
dev.off()

write.table(file=paste(opt$outfile,"outliers",sep="."),list("outliers"=outliers,"nb_outliers"=length(outliers)))





escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
