#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript bin/cog_anno_plot.r in.stat out.png")
	print("1) in.stat: the stat file for COG anno")
	print("2) out.png: the filename for output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 ) {
	print(args)
	usage()
	stop("the length of args != 2")
}



# load library
library(ggplot2)
library(grid)


# get args
args <-commandArgs(TRUE)


# reading data
data <- read.delim(args[1], header=TRUE, sep="\t")
head(data)


# plot
df <- data.frame(Frequency=data[,3], group=data[,1])
data$ration=round(100*as.vector(data[,3])/sum(as.vector(data[,3])),digits=2)
labels <- paste(data[,1],": " ,data[,2]," [",data[,3],"~",data$ration,"%","]",sep="")
p <- ggplot(data=df, aes(x=group, y=Frequency)) + geom_bar(aes(fill=group), stat="identity")
p <- p + scale_fill_discrete(name="", breaks=sort(data[,1]), labels=sort(labels))
p <- p + theme(legend.key.size=unit(0.5, "cm")) 
p <- p + labs(x="Function Class", title="COG Function Classification of Consensus Sequence")


# output plot
png(filename=args[2], height = 3000, width = 5000, res = 500, units = "px")
print(p)
dev.off()




