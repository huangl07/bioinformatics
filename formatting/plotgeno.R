#! /mnt/ilustre/app/pub/R/bin/Rscript

library(getopt)
opt = getopt(matrix(c(
  'map','m',1,'character',
  'dir','d',2,'character',
  'fKey','k',3,'character',
  'help','h',0,'logical'
),byrow=TRUE, ncol=4));
usage<-function(){
cat("This script is used to draw genotype figure
Usage   Rscript plotgeno.R for fig 3-12[options]
Options:
	-m, --map	map file with marker codes, forced
	-d, --dir	outputdir, forced
	-k, --fKey	output file stem, forced
	-h, --help	print display this help and exit
")
  q(status=1);
}
if (!is.null(opt$help) ) { usage() }
if (is.null(opt$map) ) { usage() }
if (is.null(opt$dir) ) { usage() }
if (is.null(opt$fKey)) { usage() }

##Get Rscript file dir##
getProgramName<-function(arguments){
  args <- commandArgs(trailingOnly = FALSE)
  sub("--file=", "", args[grep("--file=", args)])
}
flname <- getProgramName()
bindir <- dirname(flname)

##install.packages("qtl")
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if (!is.installed("qtl")){
  rqtlpack = paste(bindir,"/qtl_1.39-5.tar.gz",sep="")
  install.packages(rqtlpack, repos=NULL)}
library(qtl)

##plotgeno
plotgeno<-function (x, chr, ind, include.xo = TRUE, horizontal = TRUE, 
    cutoff = 4, min.sep = 2, cex = 1.2, color=c("deepskyblue", "#7CFC00","gold","orange","black","red"),pch=4,...) 
{
    cross <- x
    if (!any(class(cross) == "cross")) 
        stop("Input should have class \"cross\".")
    if (missing(chr)) 
        chr <- names(cross$geno)[1]
    cross <- subset(cross, chr = chr)
    if (nchr(cross) > 1) 
        cross <- subset(cross, chr = names(cross$geno)[1])
    if (!missing(ind)) {
        if (is.null(getid(cross))) 
            cross$pheno$id <- 1:nind(cross)
        if (!is.logical(ind)) 
            ind <- unique(ind)
        cross <- subset(cross, ind = ind)
    }
    id <- getid(cross)
    if (is.null(id)) 
        id <- 1:nind(cross)
    use.id <- TRUE
    type <- class(cross)[1]
    old.las <- par("las")
    on.exit(par(las = old.las))
    par(las = 1)
    if (!("errorlod" %in% names(cross$geno[[1]]))) {
        warning("First running calc.errorlod.")
        cross <- calc.errorlod(cross, error.prob = 0.01)
    }
    errors <- matrix(0, ncol = ncol(cross$geno[[1]]$data), nrow = nrow(cross$geno[[1]]$data))
    dimnames(errors) <- dimnames(cross$geno[[1]]$data)
    top <- top.errorlod(cross, names(cross$geno)[1], cutoff, 
        FALSE)
    if (length(top) > 0) 
        for (i in 1:nrow(top)) errors[match(top[i, 2], id), as.character(top[i, 
            3])] <- 1
    map <- cross$geno[[1]]$map
    if (is.matrix(map)) 
        map <- map[1, ]
    L <- diff(range(map))
    min.d <- L * min.sep/100
    d <- diff(map)
    d[d < min.d] <- min.d
    map <- cumsum(c(0, d))
    cross$geno[[1]]$map <- map
    n.ind <- nrow(errors)
    data <- cross$geno[[1]]$data
    chrtype <- class(cross$geno[[1]])
    if (chrtype == "X" && (type == "f2" || type == "bc")) 
        data <- reviseXdata(type, sexpgm = getsex(cross), geno = data, 
            cross.attr = attributes(cross), force = TRUE)
    if (include.xo) {
        if (type != "4way") {
            xoloc <- locateXO(cross)
            xoloc <- data.frame(ind = rep(1:length(xoloc), sapply(xoloc, 
                length)), loc = unlist(xoloc), stringsAsFactors = TRUE)
        }
        else {
            mcross <- dcross <- cross
            class(mcross)[1] <- class(dcross)[1] <- "bc"
            mcross$geno[[1]]$data[!is.na(data) & data == 1 | 
                data == 3 | data == 5] <- 1
            mcross$geno[[1]]$data[!is.na(data) & data == 2 | 
                data == 4 | data == 6] <- 2
            mcross$geno[[1]]$data[!is.na(data) & data == 7 | 
                data == 8 | data == 9 | data == 10] <- NA
            dcross$geno[[1]]$data[!is.na(data) & data == 1 | 
                data == 2 | data == 7] <- 1
            dcross$geno[[1]]$data[!is.na(data) & data == 3 | 
                data == 4 | data == 8] <- 2
            dcross$geno[[1]]$data[!is.na(data) & data == 5 | 
                data == 6 | data == 9 | data == 10] <- NA
            mxoloc <- locateXO(mcross)
            mxoloc <- data.frame(ind = rep(1:length(mxoloc), 
                sapply(mxoloc, length)), loc = unlist(mxoloc), 
                stringsAsFactors = TRUE)
            dxoloc <- locateXO(dcross)
            dxoloc <- data.frame(ind = rep(1:length(dxoloc), 
                sapply(dxoloc, length)), loc = unlist(dxoloc), 
                stringsAsFactors = TRUE)
        }
    }
    args <- list(...)
    if ("main" %in% names(args)) 
        themain <- args$main
    else themain <- paste("Lg", names(cross$geno)[1])
    if ("xlim" %in% names(args)) 
        thexlim <- args$xlim
    else thexlim <- NULL
    if ("ylim" %in% names(args)) 
        theylim <- args$ylim
    else theylim <- NULL
    if (type == "4way") {
        jit <- 0.15
        mdata <- data
        ddata <- data
        mdata[!is.na(data) & (data == 1 | data == 3 | data == 
            5)] <- 1
        mdata[!is.na(data) & (data == 2 | data == 4 | data == 
            6)] <- 2
        mdata[!is.na(data) & (data == 7 | data == 8)] <- NA
        ddata[!is.na(data) & (data == 1 | data == 2 | data == 
            7)] <- 1
        ddata[!is.na(data) & (data == 3 | data == 4 | data == 
            8)] <- 2
        ddata[!is.na(data) & (data == 5 | data == 6)] <- NA
        if (horizontal) {
            if (is.null(thexlim)) 
                thexlim <- c(0, max(map))
            if (is.null(theylim)) 
                theylim <- c(n.ind + 1, 0)
            plot(0, 0, type = "n", xlab = "Marker", ylab = "Individual", 
                main = themain, ylim = theylim, xlim = thexlim, 
                yaxt = "n", yaxs = "i",xaxt="n")
            segments(0, 1:n.ind - jit, max(map), 1:n.ind - jit)
            segments(0, 1:n.ind + jit, max(map), 1:n.ind + jit)
            if (use.id) 
                axis(side = 2, at = 1:n.ind, labels = id,cex.axis=cex-0.4)
            else axis(side = 2, at = 1:n.ind,cex.axis=cex-0.4)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 1] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind - jit, pch = 21, col = "black", bg = color[1], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 2] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind - jit, pch = 21, col = "black", bg = color[3], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 9] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind - jit, pch = 21, col = "black", bg = color[4], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 10] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind - jit, pch = 21, col = "black", bg = color[5], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 1] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind + jit, pch = 21, col = "black", bg = color[1], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 2] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind + jit, pch = 21, col = "black", bg = color[3], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 9] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind + jit, pch = 21, col = "black", bg = color[4], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 10] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind + jit, pch = 21, col = "black", bg = color[5], 
                cex = cex)
            u <- par("usr")
            segments(map, u[3], map, u[3] - 1/2)
            segments(map, u[4], map, u[4] + 1/2)
            if (any(errors != 0)) {
                ind <- rep(1:n.ind, length(map))
                ind[errors != 1] <- NA
                points(x, ind - jit, pch = 0, col = color[6], 
                  cex = cex + 0.4, lwd = 2)
                points(x, ind + jit, pch = 0, col = color[6], 
                  cex = cex + 0.4, lwd = 2)
            }
            if (include.xo) {
                points(mxoloc$loc, mxoloc$ind - jit, pch = pch, 
                  col = "blue", lwd = 2,cex = cex)
                points(dxoloc$loc, dxoloc$ind + jit, pch = pch, 
                  col = "blue", lwd = 2,cex = cex)
            }
        }
        else {
            if (is.null(theylim)) 
                theylim <- c(max(map), 0)
            if (is.null(thexlim)) 
                thexlim <- c(0, n.ind + 1)
            plot(0, 0, type = "n", ylab = "Marker", xlab = "Individual", 
                main = themain, xlim = thexlim, ylim = theylim, 
                xaxt = "n", xaxs = "i",yaxt="n")
            segments(1:n.ind - jit, 0, 1:n.ind - jit, max(map))
            segments(1:n.ind + jit, 0, 1:n.ind + jit, max(map))
            if (use.id) 
                axis(side = 1, at = 1:n.ind, labels = id,cex.axis=cex-0.4)
            else axis(side = 1, at = 1:n.ind,cex.axis=cex-0.4)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 1] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind - jit, y, pch = 21, col = "black", bg = color[1], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 2] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind - jit, y, pch = 21, col = "black", bg = color[3], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 9] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind - jit, y, pch = 21, col = "black", bg = color[4], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(mdata)] <- NA
            ind <- tind
            ind[!is.na(mdata) & mdata != 10] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind - jit, y, pch = 21, col = "black", bg = color[5], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 1] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind + jit, y, pch = 21, col = "black", bg = color[1], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 2] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind + jit, y, pch = 21, col = "black", bg = color[3], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 9] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind + jit, y, pch = 21, col = "black", bg = color[4], 
                cex = cex)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(ddata)] <- NA
            ind <- tind
            ind[!is.na(ddata) & ddata != 10] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind + jit, y, pch = 21, col = "black", bg = color[5], 
                cex = cex)
            u <- par("usr")
            segments(u[1], map, (u[1] + 1)/2, map)
            segments(u[2], map, (n.ind + u[2])/2, map)
            if (any(errors != 0)) {
                ind <- rep(1:n.ind, length(map))
                ind[errors != 1] <- NA
                points(ind - jit, y, pch = 0, col = color[6], 
                  cex = cex + 0.4, lwd = 2)
                points(ind + jit, y, pch = 0, col = color[6], 
                  cex = cex + 0.4, lwd = 2)
            }
            if (include.xo) {
                points(mxoloc$ind - jit, mxoloc$loc, pch = pch, 
                  col = "blue", lwd = 2,cex = cex)
                points(dxoloc$ind + jit, dxoloc$loc, pch = pch, 
                  col = "blue", lwd = 2,cex = cex)
            }
        }
    }
    else {
        if (horizontal) {
            if (is.null(thexlim)) 
                thexlim <- c(0, max(map))
            if (is.null(theylim)) 
                theylim <- c(n.ind + 0.5, 0.5)
            plot(0, 0, type = "n", xlab = "Marker", ylab = "Individual", 
                main = themain, ylim = theylim, xlim = thexlim, 
                yaxt = "n",xaxt="n")
            segments(0, 1:n.ind, max(map), 1:n.ind)
            if (use.id) 
                axis(side = 2, at = 1:n.ind, labels = id,cex.axis=cex-0.4)
            else axis(side = 2,cex.axis=cex-0.4)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(data)] <- NA
            ind <- tind
            ind[!is.na(data) & data != 1] <- NA
            x <- rep(map, rep(n.ind, length(map)))
            points(x, ind, pch = 21, col = "black", bg = color[1], 
                cex = cex)
            ind <- tind
            ind[!is.na(data) & data != 2] <- NA
            if (type == "f2" || (type == "bc" && chrtype == "X")) 
                points(x, ind, pch = 21, col = "black", bg = color[2], 
                  cex = cex)
            else points(x, ind, pch = 21, col = "black", bg = color[3], 
                cex = cex)
            if (type == "f2" || (type == "bc" && chrtype == "X")) {
                ind <- tind
                ind[!is.na(data) & data != 3] <- NA
                points(x, ind, pch = 21, col = "black", bg = color[3], 
                  cex = cex)
            }
            if (type == "f2") {
                ind <- tind
                ind[!is.na(data) & data != 4] <- NA
                points(x, ind, pch = 21, col = "black", bg = color[4], 
                  cex = cex)
                ind <- tind
                ind[!is.na(data) & data != 5] <- NA
                points(x, ind, pch = 21, col = "black", bg = color[5], 
                  cex = cex)
            }
            u <- par("usr")
            segments(map, u[3], map, u[3] - 1/2)
            segments(map, u[4], map, u[4] + 1/2)
            if (any(errors != 0)) {
                ind <- rep(1:n.ind, length(map))
                ind[errors != 1] <- NA
                points(x, ind, pch = 0, col = color[6], cex = cex + 
                  0.4, lwd = 2)
            }
            if (include.xo) 
                points(xoloc$loc, xoloc$ind, pch = pch, col = "blue", 
                  lwd = 2,cex = cex)
        }
        else {
            if (is.null(theylim)) 
                theylim <- c(max(map), 0)
            if (is.null(thexlim)) 
                thexlim <- c(0.5, n.ind + 0.5)
            plot(0, 0, type = "n", ylab = "Marker", xlab = "Individual", 
                main = themain, xlim = thexlim, ylim = theylim, 
                xaxt = "n",yaxt="n")
            segments(1:n.ind, 0, 1:n.ind, max(map))
            if (use.id) 
                axis(side = 1, at = 1:n.ind, labels = id,cex.axis=cex-0.4)
            else axis(side = 1,cex.axis=cex-0.4)
            tind <- rep(1:n.ind, length(map))
            tind[is.na(data)] <- NA
            ind <- tind
            ind[!is.na(data) & data != 1] <- NA
            y <- rep(map, rep(n.ind, length(map)))
            points(ind, y, pch = 21, col = "black", bg = color[1], 
                cex = cex)
            ind <- tind
            ind[!is.na(data) & data != 2] <- NA
            if (type == "f2" || (type == "bc" && chrtype == "X")) 
                points(ind, y, pch = 21, col = "black", bg = color[2], 
                  cex = cex)
            else points(ind, y, pch = 21, col = "black", bg = color[3], 
                cex = cex)
            if (type == "f2" || (type == "bc" && chrtype == "X")) {
                ind <- tind
                ind[!is.na(data) & data != 3] <- NA
                points(ind, y, pch = 21, col = "black", bg = color[3], 
                  cex = cex)
            }
            if (type == "f2") {
                ind <- tind
                ind[!is.na(data) & data != 4] <- NA
                points(ind, y, pch = 21, col = "black", bg = color[4], 
                  cex = cex)
                ind <- tind
                ind[!is.na(data) & data != 5] <- NA
                points(ind, y, pch = 21, col = "black", bg = color[5], 
                  cex = cex)
            }
            u <- par("usr")
            segments(u[1], map, (u[1] + 1)/2, map)
            segments(u[2], map, (n.ind + u[2])/2, map)
            if (any(errors != 0)) {
                ind <- rep(1:n.ind, length(map))
                ind[errors != 1] <- NA
                points(ind, y, pch = 0, col = color[6], cex = cex + 
                  0.4, lwd = 2)
            }
            if (include.xo) 
                points(xoloc$ind, xoloc$loc, pch = pch, col = "blue", 
                  lwd = 2,cex = cex)
        }
    }
    invisible()
}

g<-read.cross("csvr", "", opt$map)
s<-summary(g)
nm<-t(data.frame(s$n.mar))

##draw fig 3-12, this may take tens of hours
g1<-calc.errorlod(g,map.function="kosambi")
#pos<-seq(1,length(g1$pheno[,1]))

nam3d=paste(opt$dir,"/figure_3_12_",opt$fKey,"/pdf",sep="")
dir.create(nam3d,showWarnings=F,recursive =T)
nam3pd=paste(opt$dir,"/figure_3_12_",opt$fKey,"/png",sep="")
dir.create(nam3pd,showWarnings=F,recursive =T)

for(i in 1:length(names(g$geno))){
  if (s$n.ind/10 < 10) {
    h<-10
  }
  else {
    h<-s$n.ind/10
  }
  if (nm[,i]/5 < h*0.8){
    w<-h*0.8
  }
  else{
    w<-nm[,i]/5
  }
  #name<-paste("Lg",names(g$geno)[i],"geno.pdf",sep="")
  nam3=paste(nam3d,"/figure_3_12_",opt$fKey,"_",names(g$geno)[i],".pdf",sep="")
  pdf(nam3,height=h,width=w)
  plotgeno(g1,chr=names(g$geno)[i],min.sep=3,cex=0.8)
  dev.off
  nam3p=paste(nam3pd,"/figure_3_12_",opt$fKey,"_",names(g$geno)[i],".png",sep="")
  png(nam3p,height=h,width=w,units = "in", res=300)
  plotgeno(g1,chr=names(g$geno)[i],min.sep=3,cex=0.8)
  dev.off()
}
#axis(side=2, at=pos,las=1,labels=g1$pheno[,1],cex.axis=0.5)
