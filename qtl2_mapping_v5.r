## Install and load libraries (library location C:\ProgramData\R\win-library\3.6).

# install.packages("qtl")
# install.packages("ASMap")
# install.packages("LinkageMapView")
# install.packages("qtl2", repos="https://rqtl.org/qtl2cran")
# install.packages("corrplot")
# remotes::install_github("tavareshugo/qtl2helper")

##### Genetic linkage map construction (QTL and ASMap) #############################################################################

library(qtl)
library(ASMap)
library(LinkageMapView)
library(ggplot2)

## Reading the genotype and phenotype data.
Cross <- read.cross("csv", ".", "BDpop6.csv", crosstype="dh", na.strings=c("-","u"), genotypes=c("a","b"), alleles=c("A","B"))
summary(Cross)

## Plot missing genotypes.
par(mfrow=c(1,1), las=1)
plotMissing(Cross)

## Find and remove individuals with large number of missing genotypes (>=25).
sg <- statGen(Cross, bychr=FALSE, stat.type="miss", id="ID")
sg$miss
sub1 <- subset(Cross, ind = sg$miss<25)
plotMissing(sub1)

## Identify individuals with likely identical genotypes (clones).
gc <- genClones(sub1, tol=0.95, id="ID")
gc$cgd

## If some of the combinations of identical individuals are unlikely to be clones, they should be removed from the table prior to merging clones into a consensus genotype.
cgd <- gc$cgd[-c(3),]  # Add row numbers of individuals that are unlikely to be clones, to remove these from the table.
sub2 <- fixClones(sub1, cgd, consensus=TRUE, id="ID")

## If all combinations of identical individuals in the table are likely to be clones, the genotypes of the clones are merged into a consensus genotype.
# sub2 <- fixClones(sub1, gc$cgd, consensus=TRUE, id="ID")

## Check for segregation distortion.
profileMark(sub2, stat.type=c("seg.dist","prop","miss"), crit.val="bonf", layout=c(1,4), type="l")
profileMark(sub2, chr="1", stat.type=c("prop","miss"), display.markers=TRUE, mark.line=TRUE)
profileMark(sub2, chr="2", stat.type=c("prop","miss"), display.markers=TRUE, mark.line=TRUE)
profileMark(sub2, chr="4", stat.type=c("prop","miss"), display.markers=TRUE, mark.line=TRUE)
profileMark(sub2, chr="5", stat.type=c("prop","miss"), display.markers=TRUE, mark.line=TRUE)

## Remove strongly distorted markers.
sm1 <- statMark(sub2, stat.type="marker")$marker$AA
dm1 <- markernames(sub2)[(sm1>0.85) | (sm1<0.15)]
sub3 <- drop.markers(sub2, dm1)

## Remove markers with large number of missing genotypes (>50%).
sm2 <- statMark(sub3, stat.type="marker")$marker$missing
dm2 <- markernames(sub3)[(sm2>0.2)]
sub4 <- drop.markers(sub3, dm2)
plotMissing(sub4)

## Pull out markers before genetic linkage map construction, markers can be pushed back in a later stage.
# pp <- pp.init(seg.thresh=0.05, seg.ratio=NULL, miss.thresh=0.1)  # Defining the thresholds using pp.init.
pp <- pp.init(seg.thresh="bonf", seg.ratio=NULL, miss.thresh=0.2)  # Defining the thresholds using pp.init.
sub4 <- pullCross(sub4, type="co.located")
# sub4 <- pullCross(sub4, type="seg.distortion", pars=pp)
# sub4 <- pullCross(sub4, type="missing", pars=pp)

## Print an overview of pulled markers and the summary
sub4$co.located$table
# sub4$seg.distortion$table
# sub4$missing$table

## Summarize the dataset used for constructing a genetic linkage map.
summary(sub4)

## Generate a P-value graph to determine threshold for marker clustering.
pValue(dist=seq(5,50, by=5), pop.size=100:200, map.function="kosambi", LOD=FALSE)

## Construct the genetic linkage map:
## 1) By respecting the inputted marker order.
# sub5 <- mstmap(sub4, id="ID", bychr=TRUE, anchor=TRUE, dist.fun="kosambi", objective.fun="COUNT", p.value=1e-6, detectBadData=TRUE, trace=TRUE)
## 2) By ordering markers within linkage groups.
sub5 <- mstmap(sub4, id="ID", bychr=TRUE, dist.fun="kosambi", objective.fun="COUNT", p.value=1e-2, detectBadData=TRUE, trace=TRUE)
## 3) By bulking all data and reconstructing entire linkage map.
# sub5 <- mstmap(sub4, id="ID", bychr=FALSE, dist.fun="kosambi", objective.fun="COUNT", p.value=1e-6, detectBadData=TRUE, trace=TRUE)
chrlen(sub5)

## Heatmap of the estimated pairwise recombination fractions and LOD linkage between markers.
heatMap(sub5, lmax = 20)

## Profile genotype statistics.
pg <- profileGen(sub5, bychr=FALSE, stat.type=c("xo","dxo","miss"),
                 id="ID", xo.lambda=13, layout=c(1,3), lty=2)

## Remove all individuals with an above expected number of recombinations.
sub6 <- subsetCross(sub5, ind=!pg$xo.lambda)

## Re-construct the genetic linkage map.
sub7 <- mstmap(sub6, id="ID", bychr=TRUE, dist.fun="kosambi", objective.fun="COUNT", p.value=1e-4, detectBadData=TRUE, trace=TRUE)
chrlen(sub7)
nmar(sub7)

## Profile marker and interval statistics.
profileMark(sub7, stat.type=c("seg.dist","prop","dxo","recomb"), layout=c(1,5), type="l")

## Pushing back pulled out markers.
pp <- pp.init(seg.thresh = 0.05, seg.ratio = NULL, miss.thresh=0.4, max.rf=0.25, min.lod=3)  # Defining the thresholds using pp.init.
sub8 <- pushCross(sub7, type="co.located")
# sub8 <- pushCross(sub8, type="seg.distortion", pars=pp)
# sub8 <- pushCross(sub8, type="missing", pars=pp)
pull.map(sub8)

## Break up a linkage group.
# slist <- list("2"="2.m.5")
# sub8 <- breakCross(sub8, split=slist)

## Merge linkage groups.
# mlist <- list("2"=c("2.1","2.2","2.3"))
mlist <- list("2"=c("2.1","2.2","2.3"),
              "3"=c("3.1","3.2"),
              "4"=c("4.1","4.2","4.3"),
              "7"=c("7.1","7.2"),
              "10"=c("10.1.1","10.1.2","10.2"),
              "12"=c("12.1","12.2"),
              "13"=c("13.1","13.2"))
sub8 <- mergeCross(sub8, merge=mlist)

## Move a marker from one linkage group to another linkage group.
# pull.map(sub8, chr=c("9","12"))
# sub8 <- movemarker(sub8, "Ch9_16209", 9, 0.0)

## Re-construct the genetic linkage map by only reordering markers within linkage groups.
sub9 <- mstmap(sub8, id="ID", bychr=TRUE, trace=TRUE, p.value=2)
chrlen(sub9)

## Heatmap of the estimated pairwise recombination fractions and LOD linkage between markers.
heatMap(sub9, lmax = 20)

## Profile marker and interval statistics.
# profileMark(sub9, stat.type="marker", display.markers=TRUE, layout=c(1,5), type="l")
# profileMark(sub9, stat.type="interval", map.function="kosambi", display.markers=TRUE, layout=c(1,5), type="l")
profileMark(sub9, stat.type=c("seg.dist","dxo","mrf","lod"), use.dist=FALSE, display.markers=TRUE,
            mark.line=TRUE, layout=c(1,4), type="l")

## Inspect the genetic linkage map to check if all linkage groups are in the correct orientation.
pull.map(sub9, as.table=TRUE)

## If necessary flip linkage group(s).
sub9 <- flip.order(sub9, chr=c(1,3:12))

## Compare the impact of changing the order of markers within a linkage group and change the order.
# pull.map(sub9, chr=2, as.table=TRUE)
# compareorder(sub9, chr=2, c(1:3,6,4:5,7:12), map.function="kosambi")
# sub9 <- switch.order(sub9, chr=2, c(1:3,6,4:5,7:12), map.function="kosambi")

## Extract mapping distances and plot the final genetic linkage map.
pull.map(sub9, as.table=TRUE)
plotMap(sub9, chr=c(1:13), show.marker.names=TRUE)

## Pull out all individuals used to generate the genetic linkage map.
getid(sub9)

## Pull out all markers used to generate the genetic linkage map.
markernames(sub9)

## Print the final genetic linkage map to file.
outfile <- "BD - Final genetic linkage map.pdf"
lmv.linkage.plot(sub9, outfile, mapthese=c(1:13), dupnbr=TRUE)
lmv.linkage.plot(sub9, outfile, mapthese=c(1:13),
                 lgperrow=5,
                 pdf.width=12, pdf.height=10, pdf.pointsize=10, dupnbr=TRUE,
                 cex.main=2, col.main="black")

## Print a density map to file.
outfile <- "BD - Final density map.pdf"
lmv.linkage.plot(sub9, outfile, denmap=TRUE)

## Save the genetic linkage map in an R readable format (RDS)
saveRDS(sub9, file = "Final genetic linkage map.RDS")

##### QTL analysis for binary scored phenotypes (QTL2) #############################################################################

library(qtl)
library(qtl2)
library(qtl2helper)

## Read the RDS file generated above
# sub9 <- readRDS("Final genetic linkage map.RDS")

## Histogram of phenotypes.
par(mfrow=c(3,2), las=1)
for(i in 2:6)
  plotPheno(sub9, pheno.col=i)

## Convert QTL dataset into QTL2 format (if continued after generating the genetic linkage map above).
sub10 <- convert2cross2(sub9)

## Add phenotype data to sub10
# pheno <- read.table("Phenotypes.csv", sep=",", header=T, check.names=FALSE)
# sub11 <- add_pheno(sub10, pheno, idcol="ID", retain_all=TRUE)

## Summarize the data
summary(sub10)

## Insert pseudo markers into the map.
map <- insert_pseudomarkers(sub10$gmap, step=1)

## Calculate genotype probabilities.
probs <- calc_genoprob(sub10, map, error_prob=0.002, map_function="kosambi")

## Genome scan by Haley-Knott regression using the binary model (phenotypes with 0 and 1 scores).
out_bin <- scan1(probs, sub10$pheno, model="binary")

## Permutation test to determine the LOD significance threshold.
operm_bin <- scan1perm(probs, sub10$pheno, model="binary", n_perm=1000)

summary(operm_bin, alpha=c(0.2, 0.1, 0.05))

## Set LOD significance threshold.
thr <- summary(operm_bin, alpha=0.05)

## Plot LOD curves for five phenotypes.
pheno_names(sub10)

svg("agaricus_mapping.svg",width=10,height=15)
par(mfrow=c(4,1), las=1, mar=c(4, 3, 1, 1),cex=1.5)
ymx <- maxlod(out_bin)
plot(out_bin, map, lodcolumn=2, col="black",ylim=c(0,ymx*1.02),bgcolor="#FAFAFA",altbgcolor="#F0F0F0")
add_threshold(map, summary(operm_bin, alpha=0.05), col="grey70", lty=2,lwd=3)
title("A) Set 1: BA compatibility",xlab=NA)
plot(out_bin, map, lodcolumn=3, col="black",ylim=c(0,ymx*1.02),bgcolor="#FAFAFA",altbgcolor="#F0F0F0")
add_threshold(map, summary(operm_bin, alpha=0.05), col="grey70", lty=2,lwd=3)
title("B) Set 2: DA compatibility",font=2)
plot(out_bin, map, lodcolumn=4, col="black",ylim=c(0,ymx*1.02),bgcolor="#FAFAFA",altbgcolor="#F0F0F0")
add_threshold(map, summary(operm_bin, alpha=0.05), col="grey70", lty=2,lwd=3)
title("C) Set 3: DC compatibility",font=3)
plot(out_bin, map, lodcolumn=5, col="black",ylim=c(0,ymx*1.02),bgcolor="#FAFAFA",altbgcolor="#F0F0F0")
add_threshold(map, summary(operm_bin, alpha=0.05), col="grey70", lty=2,lwd=3)
title("D) Set 4: BC compatibility",font=4,outer=F)
dev.off()



## Function to calculate the percentage of phenotype variation explained (PVE) for a QTL.
getPVE <- function(LOD, N) {
  100 * (1 - 10^((-2*LOD)/N))
}

## Find significant QTL peaks and add the PVE value and peak marker (incl pseudomarkers) for each peak.
peaks <- find_peaks(out_bin, map, threshold=summary(operm_bin), peakdrop=1, prob=0.95, expand2markers=TRUE)
peaks$pve <- getPVE(peaks$lod, n_ind_pheno(sub10))
peaks$peak_marker_ps <- find_marker(map, chr=peaks$chr, pos=peaks$pos)
peaks$peak_marker <- find_marker(sub10$gmap, chr=peaks$chr, pos=peaks$pos)
peaks$lo_ci_marker <- find_marker(sub10$gmap, chr=peaks$chr, pos=peaks$ci_lo)
peaks$hi_ci_marker <- find_marker(sub10$gmap, chr=peaks$chr, pos=peaks$ci_hi)
peaks

## Estimate and plot QTL effect along a chromosome (i.e. of phenotype "compBC" on chromosome 7).
coeff_bin <- scan1coef(probs[,"7"], sub10$pheno[,"compBC"], model="binary")
par(mfrow=c(1,1), mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred")
plot(coeff_bin, map["7"], columns=1:2, col=col)
last_coef <- unclass(coeff_bin)[nrow(coeff_bin),]
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

## Plot phenotypes against genotypes at a single putative QTL position
par(mfcol=c(1,4))
g <- maxmarg(probs, map, chr=7, pos=87.00, return_char=TRUE)
plot_pxg(g, sub10$pheno[, "compBA"], ylab="Phenotype", main="compBA QTL on chromosome 7",sort=F)
plot_pxg(g, sub10$pheno[, "compBC"], ylab="Phenotype", main="compBC QTL on chromosome 7",sort=F)
plot_pxg(g, sub10$pheno[, "compDA"], ylab="Phenotype", main="compDA QTL on chromosome 7",sort=F)
plot_pxg(g, sub10$pheno[, "compDC"], ylab="Phenotype", main="compDC QTL on chromosome 7",sort=F)

par(mfcol=c(1,4))
g <- maxmarg(probs, map, chr=6, pos=50.00, return_char=TRUE)
plot_pxg(g, sub10$pheno[, "compBA"], ylab="Phenotype", main="compBA QTL on chromosome 6",sort=F)
plot_pxg(g, sub10$pheno[, "compBC"], ylab="Phenotype", main="compBC QTL on chromosome 6",sort=F)
plot_pxg(g, sub10$pheno[, "compDA"], ylab="Phenotype", main="compDA QTL on chromosome 6",sort=F)
plot_pxg(g, sub10$pheno[, "compDC"], ylab="Phenotype", main="compDC QTL on chromosome 6",sort=F)


