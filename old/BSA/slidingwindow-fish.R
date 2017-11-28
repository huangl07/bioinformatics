times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'col','c',0,'character',
	'win','w',0,'character',
	'step','s',0,'character',
	'method','m',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--infile	the input hapmap file
	--outfile	the trait file 
	--col	the col of chr pos index
	--win	the window size
	--step	the step size
	--method	default bp or num
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outfile)){ print_usage(spec) }
if ( is.null(opt$win)){opt$win=2000000;}
if ( is.null(opt$step)){opt$step=opt$win/40}
if ( is.null(opt$method)){opt$method="bp"}
library(qqman)

data<-read.table(opt$infile,head=FALSE);
collist <- unlist(strsplit(opt$col, split=","))
if(length(collist) ==3){
	chr<-data[[as.numeric(collist[1])]]
	pos<-data[[as.numeric(collist[2])]]
	index<-data[[as.numeric(collist[3])]]
}else{
	chr<-data[[as.numeric(collist[1])]]
	pos<-data[[as.numeric(collist[2])]]
	index1<-data[[as.numeric(collist[3])]]
	index2<-data[[as.numeric(collist[4])]]
	delta<-data[[as.numeric(collist[5])]]
}
twin<-function(pos,p1,p2){
	a<-pos[pos > p1 & pos < p2]
	if(is.null(a)){
		return(0);
	}
	return(length(a))
}
mwin<-function(pos,index,p1,p2){
	a<-index[pos > p1 & pos < p2]
	if(is.null(a) || length(a)<10){
		return(NULL);
	}
	return(mean(a))
}
print("sliding window")
chrname=unique(chr)
info<-NULL
for (i in 1:(length(chrname))){
	print(chrname[i]);
	chrpos=pos[which(chr==chrname[i])]
	backpos=chrpos;
	if (opt$method=="num"){chrpos=c(1:length(chrpos))}
	if(length(collist) ==3){
		chrindex=index[which(chr==chrname[i])]
	}else{
		chrindex1=index1[which(chr==chrname[i])]
		chrindex2=index2[which(chr==chrname[i])]
		chrdelta=delta[which(chr==chrname[i])]
	}
	chrlen=max(chrpos);
	win=ceiling(chrlen/as.numeric(opt$step));
	ss=c(1:win)
	pos1=ss*as.numeric(opt$step)-as.numeric(opt$win)/2;
	pos2=ss*as.numeric(opt$step)+as.numeric(opt$win)/2;
	if (opt$method=="num"){
		pos1[pos1 <= 0]=1;
		pos2[pos2 >= chrlen]=chrlen;
		pos1=backpos[pos1]
		pos2=backpos[pos2]
	}else{
		pos1[pos1 < 0]=0;
		pos2[pos2 > chrlen]=chrlen;
	}
	x=data.frame(pos1,pos2);
	if(length(collist) ==3){
		wmean=apply(x,MARGIN=1,function(x,y,z,a) mwin(backpos,chrindex,x[1],x[2]));
		twin=apply(x,MARGIN=1,function(x,y,z) twin(backpos,x[1],x[2]));
		info<-rbind(info,data.frame(chr=chrname[i],pos1=pos1,pos2=pos2,index=wmean,total=twin))
		info<-na.omit(info);
	}else{
		wmean1=apply(x,MARGIN=1,function(x,y,z,a) mwin(backpos,chrindex1,x[1],x[2]));
		wmean2=apply(x,MARGIN=1,function(x,y,z,a) mwin(backpos,chrindex2,x[1],x[2]));
		delta=apply(x,MARGIN=1,function(x,y,z,a) mwin(backpos,chrdelta,x[1],x[2]));
		twin=apply(x,MARGIN=1,function(x,y,z) twin(backpos,x[1],x[2]));
		print(twin)
		info<-rbind(info,data.frame(chr=chrname[i],pos1=pos1,pos2=pos2,index1=wmean1,index2=wmean2,delta=wmean3,total=twin))
		info<-na.omit(info);
	}
}
if(length(collist) ==3){
	write.table(file=paste(opt$outfile,".sliding",sep=""),info)
}else{
	write.table(file=paste(opt$outfile,".sliding",sep=""),info)
}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()