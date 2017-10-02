library(qtl)
mapthis<-read.cross("csv",file="mapthis.csv",estimate.map=FALSE)
pdf("mapthis.miss.pdf")
par(mfrow=c(1,2), las=1, cex=0.8)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual",type="h")
plot(ntyped(mapthis,"mar"), ylab="No. typed individuals", main="No. genotypes by individual",type="h")
dev.off()
#miss
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < nind(mapthis)*0.7])
mapthis <- drop.markers(mapthis, todrop)
#seg
gt <- geno.table(mapthis)
todrop <- rownames(gt[gt$P.value < 0.05,])
mapthis <- drop.markers(mapthis, todrop)
#group
mapthis <- est.rf(mapthis)
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
#plot.RF
pdf("mapthis.rf.pdf")
plotRF(mapthis)
dev.off()
toswitch <- markernames(mapthis, chr=c(5, 7:11))
mapthis <- switchAlleles(mapthis, toswitch)
toswitch <- markernames(mapthis,chr=c(4,6))
mapthis <-switchAlleles(mapthis,toswitch)
mapthis <- est.rf(mapthis)
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
mapthis <- orderMarkers(mapthis)


