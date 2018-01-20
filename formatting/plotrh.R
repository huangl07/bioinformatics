#! /mnt/ilustre/app/pub/R/bin/Rscript
library("ggplot2")
library(getopt)
# opt<-list()
# opt$dir<-"D:/RStudio/projects/fig311/"
# opt$map<-"D:/RStudio/projects/fig311/wangjun254.mappos"
# opt$chr<-"D:/RStudio/projects/fig311/wangjun254.filt.chr"
# opt$fKey<-"fig11"

opt = getopt(matrix(c(
  'map','m',1,'character',
  'chr','c',2,'character',
  'dir','d',3,'character',
  'fKey','k',4,'character',
  'help','h',0,'logical'
),byrow=TRUE, ncol=4));
usage<-function(){
  cat("This script is used to draw recombination hotspot
Usage   Rscript plotrh.R for fig 3-11 [options]
Options:        
	-m, --map	mappos file with marker codes, forced
	-c, --chr	filt.chr file with lg chr and chr length, forced
	-d, --dir	outputdir, forced
	-k, --fKey	output file stem, forced
	-h, --help	print display this help and exit
")
  q(status=1);
}
if (!is.null(opt$help) ) { usage() }
if (is.null(opt$map) ) { usage() }
if (is.null(opt$chr) ) { usage() }
if (is.null(opt$dir) ) { usage() }
if (is.null(opt$fKey)) { usage() }

##Get Rscript file dir##
getProgramName<-function(arguments){
  args <- commandArgs(trailingOnly = FALSE)
  sub("--file=", "", args[grep("--file=", args)])
}
flname <- getProgramName()
bindir <- dirname(flname)

# Multiple plot function
# 
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow =T)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


stepsize=500000#
#k=2600000
dat1<-read.table(opt$map,col.names = c("id","lg","cm","pos"))
dat2<-read.table(opt$chr,header=T,col.names = c("lg","chr","length"))

totalcmmb=vector()
figlist<-list()
#tmptab3<-data.frame()
for (i in dat2$lg){
  tmptab1<-dat1[dat1$lg==i,]
  j=1
  k=0
  tmptab2<-data.frame()
  while ((k+stepsize)<=dat2$length[i]){
    tmpvec<-tmptab1[(tmptab1$pos>=k)&(tmptab1$pos<=(k+stepsize)),3]
    if (length(tmpvec)==0){tmpvec<-0}
    tmptab2[j,1]<-(max(tmpvec)-min(tmpvec))*1000000/stepsize
    tmptab2[j,2]<-k/1000000
    tmptab2[j,3]<-dat2[i,2]
    k=k+stepsize
    j=j+1
  }
#  tmptab3<-rbind(tmptab3,tmptab2)
  totalcmmb<-union(totalcmmb,tmptab2[,1])
  colnames(tmptab2)=c("cmmb","mpos","chr")
  p = ggplot(tmptab2, aes(x = mpos, y = cmmb))+
    geom_bar(stat="identity",fill = "black",width=0.5)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.border = element_blank(),panel.grid.minor = element_blank())+
    theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+
    labs(fill="",x=paste(dat2[i,2],"(Mb)",sep=""),y="cM/Mb")
  figlist<-c(figlist,list(p))
}
#colnames(tmptab3)=c("cmmb","mpos","chr")

colnum<-0
if (length(dat2$lg)==1){colnum=1
}else if (length(dat2$lg)<=4){colnum=2
}else if (length(dat2$lg)<=9){colnum=3
}else if (length(dat2$lg)<=16){colnum=4
}else {colnum=5
}

dir.create(opt$dir,showWarnings=F,recursive =T)
nam1=paste(opt$dir,"/figure_3_11_",opt$fKey,".pdf",sep="")
pdf(nam1,width=16,height=9)
multiplot(plotlist = figlist,cols=colnum)
dev.off()
nam1p=paste(opt$dir,"/figure_3_11_",opt$fKey,".png",sep="")
png(nam1p,width=16,height=9,units="in", res=300)
multiplot(plotlist = figlist,cols=colnum)
dev.off()



