library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'vcf','v',0,'character',
	'out','o',0,'character',
	'pop','p',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--vcf	the input file
	--o	the out file \
	--pop	the col number for draw
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$vcf)){ print_usage(spec) }
if ( is.null(opt$out)){ print_usage(spec) }
if ( is.null(opt$pop)){ print_usage(spec) }
times<-Sys.time()
library("SNPRelate")
vcf.infile<-opt$vcf
metadata.infile<- opt$pop
transformNegativeFSTtozero<- TRUE
keepLowerTriangularMatrixOnly<- TRUE

#########
# Parsing VCF file 
########

snpgdsVCF2GDS(vcf.infile, "ccm.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm.gds")

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))


metadata=read.table(file = metadata.infile,header = T,sep = "\t",stringsAsFactors = F)
metadata=metadata[order(factor(metadata$samples, levels = sample.id)),]
pop_code=metadata$pop
poplevels=levels(as.factor(pop_code))

#####################################################
################## Pairwise FST  ####################
#####################################################



sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# pairwise populations matrix creation
res= outer(X= poplevels , Y= poplevels, 
           FUN = function(h,k){
             paste(h,k,sep = "/")
           }
)
colnames(res)=poplevels
rownames(res)=poplevels

as.data.frame(res)


# pairwise population matrix FST calculation

for(i in poplevels) {
  for(j in poplevels) {
    popelem= unlist(strsplit(res[i,j],"/"))
    
    #takes selection of samples an population to use for each pair
    flag<- pop_code %in% c(popelem[1],popelem[2])
    samp.sel<- sample.id[flag]
    pop.sel<- pop_code[flag]
    result<-NULL;
    if (popelem[1]==popelem[2]){result$Fst="0"}else{
      result = snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), 
                         autosome.only=FALSE, method="W&C84")
    }
    res[i,j]=as.character(result$Fst)
  }
}

########
# final edits and prints
########

# cell transformation from characters to numerics

res=data.frame(apply(res, 2, function(x) as.numeric(as.character(x))))
rownames(res)=poplevels

# trasforming negative FST into zero
# Negative Fst are technical artifact of the computation (see Roesti el al. 2012) and will be automatically replaced with zero inside this function.
if (transformNegativeFSTtozero==TRUE){
  for(i in poplevels) {
    for(j in poplevels) {
      if (res[i,j]<0){res[i,j]<-0} 
    }
  }
}

# keep only the lower triangular matrix
# (set the upper triangular to zero)
if(keepLowerTriangularMatrixOnly == TRUE){
  res[upper.tri(res)]<-0
}


########
# PRINT pairwise FST matrix to text file
########
outfile <- paste(opt$out,".matrix", sep="")
write.table(x = res, file = outfile, sep = "\t", dec = ".", 
            quote = F, row.names = T,col.names = NA)




###############################################################
############## FSTPAIRWISE TABLE PLOT #########################
# taken from Arlequin's pairFstMatrix.r (Author: Heidi Lischer)
# with very few edits
###############################################################


numericMatrix=res

# preliminar functions
#----Mirror matrix (left-right)----
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

#----Rotate matrix 270 clockworks----
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}

Matrix <- rotate270.matrix(numericMatrix)



ColorRamp <- colorRampPalette(c("white", "steelblue1", "blue3"))

outfileGraphic <- paste(opt$out,".png", sep="")
# outfileGraphic <- paste(outfile, "pairFstMatrix ", timeAttr, ".pdf", sep="")

#save graphic
png(outfileGraphic, width=1300, height=1300, res=144)
  # pdf(outfileGraphic, width = 10, height = 10)  
  
  smallplot <- c(0.874, 0.9, 0.18, 0.83)
  bigplot <- c(0.13, 0.85, 0.14, 0.87)
  
  old.par <- par(no.readonly = TRUE)
  
  # draw legend --------------------------------
  par(plt = smallplot)
  
  # get legend values
  Min <- min(Matrix, na.rm=TRUE)
  Max <- max(Matrix, na.rm=TRUE)
  binwidth <- (Max - Min) / 64
  y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
  z <- matrix(y, nrow = 1, ncol = length(y))
  
  image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)
  
  # adjust axis if only one value exists
  if(Min == Max){
    axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
  } else {
    axis(side=4, las = 2, cex.axis=0.8)
  }
  
  box()
  mtext(text=expression(bold(F[ST])), side=4, line=2.5, cex=1.1)
  
  
  #draw main graphic ---------------------------
  a <- ncol(numericMatrix)
  b <- nrow(numericMatrix)
  
  x <- c(1:a)
  y <- c(1:b)
  
  par(new = TRUE, plt = bigplot)
  
  image(x,y,as.matrix(Matrix), col=ColorRamp(64),
        main=expression(bold(Matrix~of~pairwise~F[ST])), xlab="",
        ylab="", axes=FALSE)
  box()
  
  #add labels
  Labels=poplevels
  if(is.null(Labels)){
    axis(1, at = c(1:a))
    axis(2, at = c(1:b), labels=c(b:1))
    mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,
          font=2)
    mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,
          font=2)
  } else{
    axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=0.75,
         las=2)
    axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=0.75,
         las=2)
  }
  
  
  par(old.par)  #reset graphic parameters

dev.off()

outfileGraphic <- paste(opt$out,".pdf", sep="")
# outfileGraphic <- paste(outfile, "pairFstMatrix ", timeAttr, ".pdf", sep="")

#save graphic
pdf(outfileGraphic, width=1300, height=1300)
  # pdf(outfileGraphic, width = 10, height = 10)  
  
  smallplot <- c(0.874, 0.9, 0.18, 0.83)
  bigplot <- c(0.13, 0.85, 0.14, 0.87)
  
  old.par <- par(no.readonly = TRUE)
  
  # draw legend --------------------------------
  par(plt = smallplot)
  
  # get legend values
  Min <- min(Matrix, na.rm=TRUE)
  Max <- max(Matrix, na.rm=TRUE)
  binwidth <- (Max - Min) / 64
  y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
  z <- matrix(y, nrow = 1, ncol = length(y))
  
  image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)
  
  # adjust axis if only one value exists
  if(Min == Max){
    axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
  } else {
    axis(side=4, las = 2, cex.axis=0.8)
  }
  
  box()
  mtext(text=expression(bold(F[ST])), side=4, line=2.5, cex=1.1)
  
  
  #draw main graphic ---------------------------
  a <- ncol(numericMatrix)
  b <- nrow(numericMatrix)
  
  x <- c(1:a)
  y <- c(1:b)
  
  par(new = TRUE, plt = bigplot)
  
  image(x,y,as.matrix(Matrix), col=ColorRamp(64),
        main=expression(bold(Matrix~of~pairwise~F[ST])), xlab="",
        ylab="", axes=FALSE)
  box()
  
  #add labels
  Labels=poplevels
  if(is.null(Labels)){
    axis(1, at = c(1:a))
    axis(2, at = c(1:b), labels=c(b:1))
    mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,
          font=2)
    mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,
          font=2)
  } else{
    axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=0.75,
         las=2)
    axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=0.75,
         las=2)
  }
  
  
  par(old.par)  #reset graphic parameters

dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
