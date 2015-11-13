library(reshape)
library(ggplot2)
library(preprocessCore)
library(pca3d)
library(DESeq2)
library(anota)
library(edgeR)
library(gage)
library(gageData)
library(XLConnect)
library(pathview)
library(gplots)
library(RColorBrewer)
library(anotaBatch)
library(GOstats)
library(cluster)
library(sciplot)
library(GO.db)

# Make two annotation sets for the later GO analyses HOW TO CHOOSE?

# 1
data(go.sets.hs)
data(egSymb)
go.gs.sym <- lapply(go.sets.hs, eg2sym)


# 2 Initially, dataset 1 was used for GAGE analysis but org.Hs.eg$BP is used in GOstats so we try to make a GAGE compatible datset from org.Hs.eg
GOdf <- as.data.frame(org.Hs.egGO2EG)
BPs <- unique(subset(GOdf, GOdf[,"Ontology"]=="BP")[,"go_id"])
MFs <- unique(subset(GOdf, GOdf[,"Ontology"]=="MF")[,"go_id"])
CCs <- unique(subset(GOdf, GOdf[,"Ontology"]=="CC")[,"go_id"])

db          <- as.list(org.Hs.egGO2EG)
db_uniq     <- lapply(db, function(x) unique(x))
db_uniq_sym <- lapply(db_uniq, eg2sym)
db_uniq_sym_BP <- db_uniq_sym[names(db_uniq_sym) %in% BPs]
db_uniq_sym_MF <- db_uniq_sym[names(db_uniq_sym) %in% MFs]
db_uniq_sym_CC <- db_uniq_sym[names(db_uniq_sym) %in% CCs]

goTerms            <- as.list(GOTERM)
GOids_BP              <- names(db_uniq_sym_BP)
termNames_BP          <- sapply(GOids_BP, function(x) Term(goTerms[[x]]))
names(db_uniq_sym_BP) <- termNames_BP

goTerms            <- as.list(GOTERM)
GOids_MF              <- names(db_uniq_sym_MF)
termNames_MF          <- sapply(GOids_MF, function(x) Term(goTerms[[x]]))
names(db_uniq_sym_MF) <- termNames_MF

goTerms            <- as.list(GOTERM)
GOids_CC              <- names(db_uniq_sym_CC)
termNames_CC          <- sapply(GOids_CC, function(x) Term(goTerms[[x]]))
names(db_uniq_sym_CC) <- termNames_CC

# User defined functions
SplitData <- function(data, anot, RNAorigin = "RNA", polysomal = "Polysomal", total = "Total"){
  if (RNAorigin %in% colnames(anot)) {
  } else {
    stop("There is no RNA column!")
  }
  if (identical(colnames(data), rownames(anot))) {
    sampleNames <- dimnames(data)[[2]]
    
    ind1 <-  anot[, RNAorigin] == polysomal
    dataP <- data[, ind1]
    anotP <- anot[ind1, ]
    rowNamesP <- sampleNames[ind1]
    
    ind2 <-  anot[, RNAorigin] == total
    dataC <- data[, ind2]
    anotC <- anot[ind2, ]
    rowNamesC<- sampleNames[ind2]
    
    rownames(anotP) <- rowNamesP
    colnames(anotP) <- colnames(anot)
    rownames(anotC) <- rowNamesC
    colnames(anotC) <- colnames(anot)
    
    colnames(dataP) <- rowNamesP
    rownames(dataP) <- rownames(data)
    colnames(dataC) <- rowNamesC
    rownames(dataC) <- rownames(data)
  } else {
    stop("Rownames of the anotation and column names of the data are not identical!")
  }
  newList        <- list(dataP, anotP, dataC, anotC)
  names(newList) <- c("dataP", "anotP", "dataC", "anotC")
  return(newList)
}

CalculateFDRList <- function(x, method = "BH") {
  FDR <- p.adjust(x$Pvalue, method = method)
  x$FDR<-FDR
  return(x)
}

# Read in output of the rpkmforgenes.py script from present working directory
rawData <- read.table("NKcellSamples_rpkms_filtered.txt", header = TRUE)

# Select only the count data from the samples, removing among other things the RPKM
dataCount           <- rawData[, c(50:96)]
rownames(dataCount) <- rawData[, 1]
dataCount           <- as.matrix(dataCount)

# Remove rows with any zeros to get a good dataset
dataNoZero <- dataCount[apply(dataCount == 0, 1, sum) <= 2, ]

# Remove genes merged by underscore (?) since the genemodel does not allow to seperate them CHECK IF NECESSARY
noMergeGenes      <- !rownames(dataNoZero) %in% grep("_", rownames(dataNoZero), value = TRUE)
dataNoZeroNoMerge <- dataNoZero[noMergeGenes, ]

# Give correct names to data samples
sampleNames                 <- gsub("\\.1", "", colnames(dataCount))
colnames(dataNoZeroNoMerge) <- sampleNames
data                        <- dataNoZeroNoMerge 

# Read in the anotation file
anot <- as.matrix(read.table("Anot.csv", header = TRUE, sep=","))

# Change order of the samples to correct order, there was discrepancy between the two supplied anotation files, PCAs revealed the most likely correct order
anot <- anot[c(1:17, 20, 21, 18, 19, 22:47), ]

# Perform some renaming
rownames(anot)                                                   <- colnames(data)
anot[anot  ==  "48 hours activation "]                           <- "A"
anot[anot  ==  "48 hours activation + 24h cytokine deprivation"] <- "S"

# Remove samples before further analysis, missing poly data for one donor, seemed to be a duplicate sample in other donor, remove the starved samples from one and the activated from the other
Samples_to_remove <- c(" 1", "28", "29", "18", "19", "42", "43")
keepSamplesInd    <- !anot[, "Samp_ID"] %in% Samples_to_remove
dataRed           <- data[, keepSamplesInd]
anotRed           <- anot[keepSamplesInd, ]

# Take rlog of data
dataRed_log           <- rlog(dataRed)
rownames(dataRed_log) <- rownames(dataRed)

# Split data accoring to anotation file, see function 
splitDat <- SplitData(dataRed_log, anotRed)

dataP <- splitDat[[1]]
dataC <- splitDat[[3]]
anotP <- splitDat[[2]]
anotC <- splitDat[[4]]

# Reduce data and anot to the starved samples
indStarvedC <- anotC[, "Treatment"] == "S"
indStarvedP <- anotP[, "Treatment"] == "S"
anotC_starved <- anotC[indStarvedC, ]
anotP_starved <- anotP[indStarvedP, ]

dataC_starved <- dataC[,indStarvedC]
dataP_starved <- dataP[,indStarvedP]

# Reduce data and anot now to the activated samples
indActivatedC <- anotC[, "Treatment"] == "A"
indActivatedP <- anotP[, "Treatment"] == "A"
anotC_activated <- anotC[indActivatedC, ]
anotP_activated <- anotP[indActivatedP, ]

dataC_activated <- dataC[, indActivatedC]
dataP_activated <- dataP[, indActivatedP]

# Create the phenoVec and batchVec needed for the anota analysis. Take care of sample order!!
anotaPhenoVec_starved   <- as.vector(anotC_starved[, "Cytokine"])
anotaPhenoVec_activated <- as.vector(anotC_activated[, "Cytokine"])
batchVec_starved        <- as.vector(anotC_starved[, "Donor"])
batchVec_activated      <- as.vector(anotC_activated[, "Donor"])

# Choose contrast to be analyzed (IL15 vs IL2)
contrasts <- as.matrix(c(1, -1), ncol = 1, nrow = 2)

# Perform  anota analysis for the starved samples first: do starved and activated seperately because there seems to be an interaction effect betweeen cytoine and treatment (anota cannot handle them yet)
anotaQcOut <- anotaPerformQc(dataT    = dataC_starved, 
                             dataP    = dataP_starved, 
                             phenoVec = anotaPhenoVec_starved
                             )

anotaResidOut <- anotaResidOutlierTest(anotaQcObj = anotaQcOut)

starvTranslation <- anotaGetSigGenesBatch(dataT      = dataC_starved, 
                                          dataP      = dataP_starved, 
                                          contrasts  = contrasts, 
                                          phenoVec   = anotaPhenoVec_starved, 
                                          anotaQcObj = anotaQcOut, 
                                          batchVec   = batchVec_starved
                                          )

starvCytosolic <- anotaGetSigGenesBatch(dataP      = dataC_starved, 
                                        phenoVec   = anotaPhenoVec_starved,
                                        contrasts  = contrasts,  
                                        anotaQcObj = anotaQcOut, 
                                        batchVec   = batchVec_starved
                                        )

starvPolysomal <- anotaGetSigGenesBatch(dataP      = dataP_starved, 
                                        phenoVec   = anotaPhenoVec_starved,
                                        contrasts  = contrasts,  
                                        anotaQcObj = anotaQcOut, 
                                        batchVec   = batchVec_starved
                                        )

# Perform anota analysis for the activated samples
anotaQcOut <- anotaPerformQc(dataT    = dataC_activated, 
                             dataP    = dataP_activated, 
                             phenoVec = anotaPhenoVec_activated
                             )

anotaResidOut <- anotaResidOutlierTest(anotaQcObj = anotaQcOut)

actTranslation <- anotaGetSigGenesBatch(dataT      = dataC_activated, 
                                        dataP      = dataP_activated,
                                        contrasts  = contrasts,  
                                        phenoVec   = anotaPhenoVec_activated, 
                                        anotaQcObj = anotaQcOut, 
                                        batchVec   = batchVec_activated
                                        )

actCytosolic <- anotaGetSigGenesBatch(dataP      = dataC_activated, 
                                      phenoVec   = anotaPhenoVec_activated,
                                      contrasts  = contrasts,  
                                      anotaQcObj = anotaQcOut, 
                                      batchVec   = batchVec_activated
                                      )

actPolysomal <- anotaGetSigGenesBatch(dataP      = dataP_activated, 
                                      phenoVec   = anotaPhenoVec_activated, 
                                      anotaQcObj = anotaQcOut, 
                                      batchVec   = batchVec_activated, 
                                      contrasts  = contrasts
                                      )

# Parameters for the anotaPlotSigGenes function
mySlopeP     <- 0.01
mySelDeltaPT <- log2(1.25)
mySelDeltaP  <- log2(1.25)
myRvmPAdj    <- 0.25
myMinEff     <- log2(1.25)
myMinEffNoT  <- log2(1.5)
myMinSlope   <- -0.5
myMaxSlope   <- 1.5
makePlots    <- TRUE

starvedGenes_translation <- anotaPlotSigGenes(selContr    = 1, 
                                              anotaSigObj = starvTranslation, 
                                              maxRvmPAdj  = myRvmPAdj, 
                                              slopeP      = mySlopeP, 
                                              selDeltaPT  = mySelDeltaPT, 
                                              selDeltaP   = mySelDeltaP, 
                                              sortBy      = "Eff", 
                                              minEff      = myMinEff, 
                                              performPlot = makePlots, 
                                              minSlope    = myMinSlope, 
                                              maxSlope    = myMaxSlope, 
                                              fileName    = "starvedGenes_Translation_no4.pdf"
                                              )

starvedGenes_cytosolic <- starvCytosolic$apvStatsRvm[[1]][abs(starvCytosolic$apvStatsRvm[[1]][, "apvEff"]) > myMinEffNoT & starvCytosolic$apvStatsRvm[[1]][, "apvRvmPAdj"] < myRvmPAdj, ]

starvedGenes_polysomal <- starvPolysomal$apvStatsRvm[[1]][abs(starvPolysomal$apvStatsRvm[[1]][, "apvEff"]) > myMinEffNoT & starvPolysomal$apvStatsRvm[[1]][, "apvRvmPAdj"] < myRvmPAdj, ]

activatedGenes_translation <- anotaPlotSigGenes(selContr    = 1, 
                                                anotaSigObj = actTranslation, 
                                                maxRvmPAdj  = myRvmPAdj, 
                                                slopeP      = mySlopeP, 
                                                selDeltaPT  = mySelDeltaPT, 
                                                selDeltaP   = mySelDeltaP, 
                                                sortBy      = "Eff", 
                                                minEff      = myMinEff, 
                                                performPlot = makePlots, 
                                                minSlope    = myMinSlope, 
                                                maxSlope    = myMaxSlope, 
                                                fileName    = "activatedGenes_Translation_no4.pdf"
                                                )

activatedGenes_cytosolic <- actCytosolic$apvStatsRvm[[1]][abs(actCytosolic$apvStatsRvm[[1]][, "apvEff"]) > myMinEffNoT & actCytosolic$apvStatsRvm[[1]][, "apvRvmPAdj"] < myRvmPAdj, ]

activatedGenes_polysomal <- actPolysomal$apvStatsRvm[[1]][abs(actPolysomal$apvStatsRvm[[1]][, "apvEff"]) > myMinEffNoT & actPolysomal$apvStatsRvm[[1]][, "apvRvmPAdj"] < myRvmPAdj, ]

# Save image after anota analysis
#save.image("NoZero_fiveDonors_log21")

# Make a scatterplot of cytosolic effect versus polysomal effect for the starved samples
effPoly    <- starvPolysomal$apvStatsRvm[[1]][, "apvEff"]
effCyto    <- starvCytosolic$apvStatsRvm[[1]][, "apvEff"]
scatterMat <- cbind(effCyto, effPoly)

# Select the significant genes to be plotted in different color
sigTrans   <- rownames(starvedGenes_translation$selectedRvmData)
sigPoly    <- rownames(starvedGenes_polysomal)

# Make smoothscatter without significant genes specified
pdf("scatterplot_log21_transGenes.pdf", useDingbats = FALSE)
smoothScatter(scatterMat, 
              nrpoints = 0, 
              xlab     = "Cytosolic RNA", 
              ylab     = "Polysome-associated RNA", 
              xlim     = c(-4, 4), 
              ylim     = c(-4, 4)
)
abline(h = 0)
abline(v = 0)

# Add points to scatterplot
transGenes <- intersect(sigTrans, sigPoly)
points(scatterMat[sigPoly, ], cex = 0.5, col = "green", pch = 16)
points(scatterMat[transGenes, ], cex = 0.5, col = "red", pch = 16)
points(scatterMat[c("BIRC5", "TOP2A", "RACGAP1", "CKS2", "NUSAP1"),, drop = F], cex = 1, col = "blue", pch = 16)
dev.off()

# Calculate the Rsquared
summary(lm(effCyto ~ effPoly))$r.squared

##########################
##########################
##########################


# Make a plot of the density of the Padj values for the starved samples

# Create data frame with adjusted p values
pvalDF <- data.frame(cyto        = starvCytosolic$apvStatsRvm[[1]][, "apvRvmPAdj"],
                     poly        = starvPolysomal$apvStatsRvm[[1]][, "apvRvmPAdj"],
                     translation = starvTranslation$apvStatsRvm[[1]][, "apvRvmPAdj"]
)

pdf(file = "pvalDensities_starved.pdf")
par(mar = c(5, 3, 2, 2) + 0.1)
plot(density(pvalDF$cyto, bw = 0.05), 
     col  = "blue",
     lwd  = 2, 
     main = "5 Donors pVal_adj Distribution"
     )
lines(density(pvalDF$poly, bw = 0.05), 
      col = "green", 
      lwd = 2
       )
lines(density(pvalDF$translation, bw = 0.05), 
      col = "red", 
      lwd = 2
      )
legend("topright", 
       legend = c("Cytosolic RNA", "Polysome-associated RNA", "Translation (Anota)"), 
       fill   = c("blue", "green", " red"), 
       cex    = 1.5
       )
dev.off()

#########################
#########################
#########################


# Plot the positive control genes

geneSet <- c("IL2RA", "NCAM1")

# Polysome-associated RNA
pdf("Polysome-associated RNA.pdf")

# Below for loops are to make sure that the effect of IL15 per gene is also per donor, normal subtraction of two matrices subsetted with the anotation matrix would result in different donor order for IL2 and IL15 subgroup
diffMat_act   <- c()
for(i in unique(anotP_activated[, "Donor"])) {
  tmp15       <- actPolysomal$inputData$dataP[geneSet, anotP_activated[, "Cytokine"] == "IL15" & anotP_activated[, "Donor"] == i]
  tmp2        <- actPolysomal$inputData$dataP[geneSet, anotP_activated[, "Cytokine"] == "IL2" & anotP_activated[, "Donor"] == i]
  tmpDiff     <- tmp15-tmp2
  diffMat_act <- cbind(diffMat_act, tmpDiff)
}

colnames(diffMat_act) <- unique(anotP_activated[, "Donor"])

diffMat_starv <- c()
for(i in unique(anotP_starved[, "Donor"])) {
  tmp15         <- starvPolysomal$inputData$dataP[geneSet, anotP_starved[, "Cytokine"] == "IL15" & anotP_starved[, "Donor"] == i]
  tmp2          <- starvPolysomal$inputData$dataP[geneSet, anotP_starved[, "Cytokine"] == "IL2" & anotP_starved[, "Donor"] == i]
  tmpDiff       <- tmp15-tmp2
  diffMat_starv <- cbind(diffMat_starv, tmpDiff)
}

colnames(diffMat_starv) <- unique(anotP_starved[, "Donor"])

diffMat_act_t   <- t(diffMat_act) 
diffMat_starv_t <- t(diffMat_starv) 
diffMat         <- cbind(diffMat_act_t, diffMat_starv_t)
colnames(diffMat) <- c("IL-2RA_act", "CD56_act", "IL-2RA_depr", "CD56_depr")
ggplotDat       <- melt(diffMat)
ggplotDat$X2    <- factor(ggplotDat$X2, levels=c("IL-2RA_act","CD56_act","IL-2RA_depr","CD56_depr"))
bargraph.CI(ggplotDat$X2, response=ggplotDat$value)
dev.off()

# Cytosolic RNA
pdf("Cytosol RNA.pdf")

# Below for loops are to make sure that the effect of IL15 per gene is also per donor, normal subtraction of two matrices subsetted with the anotation matrix would result in different donor order for IL2 and IL15 subgroup
diffMat_act   <- c()
for(i in unique(anotC_activated[, "Donor"])) {
  tmp15       <- actCytosolic$inputData$dataP[geneSet, anotC_activated[, "Cytokine"] == "IL15" & anotC_activated[, "Donor"] == i]
  tmp2        <- actCytosolic$inputData$dataP[geneSet, anotC_activated[, "Cytokine"] == "IL2" & anotC_activated[, "Donor"] == i]
  tmpDiff     <- tmp15-tmp2
  diffMat_act <- cbind(diffMat_act, tmpDiff)
}

colnames(diffMat_act) <- unique(anotC_activated[, "Donor"])

diffMat_starv <- c()
for(i in unique(anotC_starved[, "Donor"])) {
  tmp15         <- starvCytosolic$inputData$dataP[geneSet, anotC_starved[, "Cytokine"] == "IL15" & anotC_starved[, "Donor"] == i]
  tmp2          <- starvCytosolic$inputData$dataP[geneSet, anotC_starved[, "Cytokine"] == "IL2" & anotC_starved[, "Donor"] == i]
  tmpDiff       <- tmp15-tmp2
  diffMat_starv <- cbind(diffMat_starv, tmpDiff)
}

colnames(diffMat_starv) <- unique(anotC_starved[, "Donor"])

diffMat_act_t   <- t(diffMat_act) 
diffMat_starv_t <- t(diffMat_starv) 
diffMat         <- cbind(diffMat_act_t, diffMat_starv_t)
colnames(diffMat) <- c("IL-2RA_act", "CD56_act", "IL-2RA_depr", "CD56_depr")
ggplotDat       <- melt(diffMat)
ggplotDat$X2    <- factor(ggplotDat$X2, levels=c("IL-2RA_act","CD56_act","IL-2RA_depr","CD56_depr"))
bargraph.CI(ggplotDat$X2, response=ggplotDat$value)


dev.off()

########################
########################
########################

# Make a polysome and cytosolic heatmap of the polysomally significant starved genes
geneSel <- rownames(starvedGenes_polysomal)

# Poly data for heatmap
diffMat_starv_poly   <- c()
for(i in unique(anotP_starved[, "Donor"])) {
  tmp15              <- starvPolysomal$inputData$dataP[geneSel, anotP_starved[, "Cytokine"] == "IL15" & anotP_starved[, "Donor"] == i]
  tmp2               <- starvPolysomal$inputData$dataP[geneSel, anotP_starved[, "Cytokine"] == "IL2" & anotP_starved[, "Donor"] == i]
  tmpDiff            <- tmp15-tmp2
  diffMat_starv_poly <- cbind(diffMat_starv_poly, tmpDiff)
}
colnames(diffMat_starv_poly) <- unique(anotP_starved[, "Donor"])

# Cyto data for heatmap
diffMat_starv_cyto   <- c()
for(i in unique(anotP_starved[, "Donor"])){
  tmp15              <- starvCytosolic$inputData$dataP[geneSel, anotC_starved[, "Cytokine"] == "IL15" & anotP_starved[, "Donor"] == i]
  tmp2               <- starvCytosolic$inputData$dataP[geneSel, anotC_starved[, "Cytokine"] == "IL2" & anotP_starved[, "Donor"] == i]
  tmpDiff            <- tmp15-tmp2
  diffMat_starv_cyto <- cbind(diffMat_starv_cyto, tmpDiff)
}
colnames(diffMat_starv_cyto) <- unique(anotP_starved[, "Donor"])

mat <- cbind(diffMat_starv_poly, diffMat_starv_cyto)

transGenesInd <- rownames(mat) %in% rownames(starvedGenes_translation$selectedRvmData)
transGenesInd <- replace(transGenesInd, transGenesInd == FALSE, "white")
transGenesInd <- replace(transGenesInd, transGenesInd == TRUE, "black")

distMat     <- dist(mat)
distMatTree <- hclust(distMat)

pdf(file = "polysomalStarved_heatmap_log21.pdf")
lwid = c(1.5, 4)
lhei = c(0.8, 4)
heatmap.2(mat, 
          col           = "bluered", 
          breaks        = c(seq(-2, -0.5, length = 10), seq(-0.49, 0.5, length = 10), seq(0.51, 2, length = 10)), 
          dendrogram    = "row", 
          trace         = "none", 
          density.info  = "none",
          Rowv          = as.dendrogram(distMatTree),
          Colv          = NULL, 
          cexCol        = 2, 
          labRow        = "",
          RowSideColors = transGenesInd,
          lwid          = lwid,
          lhei          = lhei
)
dev.off()

##########################
##########################
##########################

##gage on Fold Changes
geneFC_cyto  <- starvCytosolic$apvStatsRvm[[1]][, "apvEff"]
geneFC_poly  <- starvPolysomal$apvStatsRvm[[1]][, "apvEff"]
geneFC_trans <- starvTranslation$apvStatsRvm[[1]][, "apvEff"]

dataSets <- c("geneFC_cyto", "geneFC_poly","geneFC_trans")

goTerms               <- list()
goTermsEss_up         <- list() 
goTermsEss_down       <- list() 
goTermsEss_genes_up   <- list()
goTermsEss_genes_down <- list()

gsets <- db_uniq_sym_BP
for (i in dataSets) {
  gageOut                    <- gage(get(i), 
                                     gsets = gsets , 
                                     ref   = NULL, 
                                     samp  = NULL
                                     )
  goTerms[[i]]               <- gageOut
  gageOut_greater_ess        <- esset.grp(gageOut$greater, 
                                          get(i), 
                                          gsets  = gsets , 
                                          ref    = NULL, 
                                          samp   = NULL, 
                                          use.q  = TRUE, 
                                          cutoff = 0.1
                                          )
  gageOut_less_ess           <- esset.grp(gageOut$less, 
                                          get(i), 
                                          gsets  = gsets , 
                                          ref    = NULL, 
                                          samp   = NULL, 
                                          use.q  = TRUE, 
                                          cutoff = 0.1
                                          )
  goTermsEss_up[[i]]         <- gageOut_greater_ess$essentialSets
  goTermsEss_down[[i]]       <- gageOut_less_ess$essentialSets
  goTermsEss_genes_up[[i]]   <- gageOut_greater_ess$coreGeneSets
  goTermsEss_genes_down[[i]] <- gageOut_less_ess$coreGeneSets
}

# Below we perform gage analysis for the cytosolic and polysome associated RNA

cutoff <- 0.05
cutoffTrans <- 0.05
# Choose number of GOs to plot
geneNo <- "Inf"

psu     <- head(goTerms$geneFC_poly$greater, geneNo)
tsu     <- head(goTerms$geneFC_trans$greater, geneNo)
csu     <- head(goTerms$geneFC_cyto$greater, geneNo)
psu_red <- subset(psu,psu[, "q.val"] < cutoff)
tsu_red <- subset(tsu,tsu[, "q.val"] < cutoffTrans)
csu_red <- subset(csu,csu[, "q.val"] < cutoff)

psd     <- head(goTerms$geneFC_poly$less, geneNo)
tsd     <- head(goTerms$geneFC_trans$less, geneNo)
csd     <- head(goTerms$geneFC_cyto$less, geneNo)
psd_red <- subset(psd,psd[, "q.val"] < cutoff)
tsd_red <- subset(tsd,tsd[, "q.val"] < cutoffTrans)
csd_red <- subset(csd,csd[, "q.val"] < cutoff)


allSigProc <- unique(c(rownames(psu_red), rownames(psd_red),
                       rownames(tsu_red), rownames(tsd_red),
                       rownames(csu_red), rownames(csd_red)))

allSigProcMat <- matrix(ncol=6, nrow=length(allSigProc), data=1)
rownames(allSigProcMat) <- allSigProc
colnames(allSigProcMat) <- c("Cytosolic down","Polysome-associated down","Translation down",
                             "Cytosolic up","Polysome-associated up", "Translation up")

allSigProcMat[allSigProc, "Polysome-associated down"] <- psd[allSigProc, "q.val"]
allSigProcMat[allSigProc, "Cytosolic down"] <- csd[allSigProc, "q.val"]
allSigProcMat[allSigProc, "Translation down"] <- tsd[allSigProc, "q.val"]

allSigProcMat[allSigProc, "Polysome-associated up"] <- psu[allSigProc, "q.val"]
allSigProcMat[allSigProc, "Cytosolic up"] <- csu[allSigProc, "q.val"]
allSigProcMat[allSigProc, "Translation up"] <- tsu[allSigProc, "q.val"]
allSigProcMat <- -log10(allSigProcMat)


goMat               <- as.matrix(allSigProcMat)
distGoMat           <- dist(goMat)
#distgoMatTree       <- hclust(distGoMat, method = "average")
col                 <- colorRampPalette(c("white", "blue", "red"))(256)
pdf(paste("Heatmap_Starved_genes_GOenrich_BP_", "top", geneNo, ".pdf", sep=""))
go_heatmap <- heatmap.2(goMat,
                        col          = col,
                        #breaks       = c(seq(0,7,length=257)),
                        margins      = c(5, 30),
                        keysize      = 1,
                        dendrogram   = "row", 
                        trace        = "none", 
                        density.info = "none",
                        labCol       = colnames(goMat),
                        Colv         = NULL,
                        cexCol       = 1,
                        cexRow       = 0.6,
                        key.xlab     = "-log10FDR",
                        main         = paste("top", geneNo, sep="_")
)
dev.off()

# After getting the gage analysis for the complete poly and cyto genelist, here we look into some specific GO terms

term1  <- "mitochondrial inner membrane"
term2  <- "Cellular Respiration"
term3  <- "antigen processing and presentation of peptide antigen via MHC class I"
term4  <- "Signaling Receptor Activity"
term5  <- "type I interferon signaling pathway"
term6  <- "tumor necrosis factor-mediated signaling pathway"
term7  <- "mitotic cell cycle"
term8  <- "programmed cell death"
term9  <- "canonical glycolysis"
term10 <- "translation"
term11 <- "regulation of immune response"

# Select all genes associated with a GO term
term    <- term10
geneSel <- unique(unlist(gsets[grep(paste(".*", term, ".*", sep=""), names(gsets), perl = TRUE, value = TRUE, ignore.case = TRUE)]))

polySigGenes <- rownames(starvedGenes_polysomal)
geneSelSig   <- intersect(geneSel, polySigGenes)

polyData15_S <- starvPolysomal$inputData$dataP[rownames(starvPolysomal$apvStatsRvm[[1]]) %in% geneSelSig, anotP_starved[, "Cytokine"] == "IL15"]
polyData2_S  <- starvPolysomal$inputData$dataP[rownames(starvPolysomal$apvStatsRvm[[1]]) %in% geneSelSig, anotP_starved[, "Cytokine"] == "IL2"]

# Reorder polyData2_S
polyData2_S_reorder <- polyData2_S[, c(5, 1, 3, 4, 2)] 

# Make plotting matrix
mat           <- polyData15_S-polyData2_S_reorder
colnames(mat) <- c("D1_IL15", "D2_IL15", "D3_IL15", "D4_IL15", "D5_IL15")

transGenesIndTerm <- rownames(mat) %in% transGenes
transGenesIndTerm <- replace(transGenesIndTerm, transGenesIndTerm == FALSE, "white")
transGenesIndTerm <- replace(transGenesIndTerm, transGenesIndTerm == TRUE, "black")

# Construct distance matrix and clustering
distMat     <- dist(mat)
distMatTree <- hclust(distMat)

pdf(file = paste(term, "_log21.pdf", sep = ""))
par(mar = c(5, 3, 2, 2) + 0.1)
lwid = c(1.5, 4)
lhei = c(0.8, 4)
plot <- heatmap.2(mat, 
                  col          = bluered,
                  breaks       = seq(-2, 2, length = 50),
                  dendrogram   = "row",
                  Colv         = NULL,
                  trace        = "none",  
                  density.info = "none",
                  Rowv         = as.dendrogram(distMatTree), 
                  cexCol       = 0.8,
                  cexRow       = 0.7,
                  main         = term,
                  lwid         = lwid,
                  lhei         = lhei,
                  RowSideColors = transGenesIndTerm,
)
dev.off()


#########################
#########################
#########################

uniGenes_starved_red <- rownames(starvedGenes_polysomal)
geneSel              <- uniGenes_starved_red

polyMeans_A_IL15  <- rowMeans(actPolysomal$inputData$dataP[rownames(actPolysomal$inputData$dataP) %in% geneSel, anotP_activated[, "Cytokine"] == "IL15"])
polyMeans_A_IL2   <- rowMeans(actPolysomal$inputData$dataP[rownames(actPolysomal$inputData$dataP) %in% geneSel, anotP_activated[, "Cytokine"] == "IL2"])
polyMeans_S_IL15  <- rowMeans(starvPolysomal$inputData$dataP[rownames(starvPolysomal$inputData$dataP) %in% geneSel, anotP_starved[, "Cytokine"] == "IL15"])
polyMeans_S_IL2   <- rowMeans(starvPolysomal$inputData$dataP[rownames(starvPolysomal$inputData$dataP) %in% geneSel, anotP_starved[, "Cytokine"] == "IL2"])
meanMat           <- cbind(polyMeans_A_IL15, polyMeans_A_IL2, polyMeans_S_IL15, polyMeans_S_IL2)
colnames(meanMat) <- c("IL15_activated", "IL2_activated", "IL15_starved", "IL2_starved")
meanMat_comb      <- meanMat[, c(1, 3, 2, 4)]
meanMat_comb_norm <- sweep(meanMat_comb, 1, rowMeans(meanMat_comb), "-")

# Exploratory analysis for number of clusters in per condition gene expression

# Kmeans clustering
nclusts   <- 10
silWidths <- vector(mode = "character", length = nclusts-1)
clustList <- list()
data      <- meanMat_comb_norm

for(i in 2:nclusts) {
  clustList[[i]] <- kmeans(x = data, centers = i, nstart = 10, iter.max = 25)
}
pdf("silhouettePlots.pdf")
for(i in 2:nclusts) {
  dissE <- daisy(data)
  dE2   <- dissE^2
  sk    <- silhouette(clustList[[i]]$cl, dE2)
  plot(sk)
}
dev.off()

for(i in 2:nclusts) {
  dissE          <- daisy(data)
  dE2            <- dissE^2
  sk             <- silhouette(clustList[[i]]$cl, dE2)
  silWidths[i-1] <- summary(sk)["avg.width"]
}
pdf("Clust_vs_silWidth_kmeans.pdf")
plot(2:nclusts, silWidths, type = "o")
dev.off()

# According to the evaluation of the plots above select here the correct number of clusters using the clustList (to make sure you use the same k means operation)
clustMat <- clustList[[4]]
m.kmeans <- cbind(data, clustMat$cluster)
transBar <- rownames(m.kmeans) %in% transGenes
transBar <- replace(transBar, transBar == TRUE, "black")
transBar <- replace(transBar, transBar == FALSE, "white")
transCol <- paste(m.kmeans[,5], transBar, sep="_")
m.kmeans <- data.frame(m.kmeans, transCol)
o        <- order(m.kmeans[, dim(m.kmeans)[2]])
m.kmeans <- m.kmeans[o, ]


clusts <- as.character(m.kmeans[, 5])
clusts <- replace(clusts, clusts == 1, "green")
clusts <- replace(clusts, clusts == 2, "blue")
clusts <- replace(clusts, clusts == 3, "brown")
clusts <- replace(clusts, clusts == 4, "orange")

plotDat <- m.kmeans[,1:dim(data)[2]]
plotDat <- as.matrix(plotDat)

transBar <- transBar[o]

pdf("Heatmap_AllCond_groupBar.pdf")
heatmap.2(plotDat, 
          col           = bluered,
          breaks        = c(seq(-1.5, 1.5, length = 30)),
          keysize       = 1,
          trace         = "none",  
          density.info  = "none", 
          Colv          = NULL,
          Rowv          = NULL,
          cexCol        = 1.4, 
          labRow        = "",
          srtCol        = 45,
          key.xlab      = "Log FC",
          RowSideColors = clusts
)
dev.off()

polyMeans <- matrix(nrow = length(unique(clusts)), ncol = 4)
polySDs   <- matrix(nrow = length(unique(clusts)), ncol = 4)
ind       <- 0

colnames(polyMeans) <- colnames(polySDs) <- colnames(data)
rownames(polyMeans) <- rownames(polySDs) <- unique(clusts)

for(color in unique(clusts)) {
  ind <- ind + 1
  
  tmpMean <- apply(m.kmeans[, 1:dim(data)[2]][clusts == color, ], 2, mean)
  tmpSD   <- apply(m.kmeans[, 1:dim(data)[2]][clusts == color, ], 2, sd)
  
  polyMeans[ind, ] <- tmpMean
  polySDs[ind, ]   <- tmpSD
}
pdf(file="groupDescription.pdf", useDingbats = FALSE)
par(mfcol = c(2, 2), las = 2)
for(i in 1:dim(polyMeans)[1]) {
  plot(polyMeans[i, ], 
       col  = rownames(polyMeans)[i], 
       type = "p", 
       pch  = 20, 
       ylim = c(-1, 1), 
       xlim = c(0.4, 4.5), 
       xaxt = "n", 
       main = paste(rownames(polyMeans)[i],  "group", sep=""), 
       xlab = "", 
       ylab = "Mean Normalized Expression", 
       cex  = 3
  )
  segments(1, polyMeans[i, 1] - polySDs[i, 1], 1, polyMeans[i, 1] + polySDs[i, 1], col = rownames(polyMeans)[i])
  segments(2, polyMeans[i, 2] - polySDs[i, 2], 2, polyMeans[i, 2] + polySDs[i, 2], col = rownames(polyMeans)[i])
  segments(3, polyMeans[i, 3] - polySDs[i, 3], 3, polyMeans[i, 3] + polySDs[i, 3], col = rownames(polyMeans)[i])
  segments(4, polyMeans[i, 4] - polySDs[i, 4], 4, polyMeans[i, 4] + polySDs[i, 4], col = rownames(polyMeans)[i])
  axis(at = c(1, 2, 3, 4), labels = NULL, side = 1)
}

dev.off()

backgroundNames <- select(org.Hs.eg.db,
                          keys    = rownames(actTranslation$apvStatsRvm[[1]]),
                          columns = c("ENTREZID"),
                          keytype = "SYMBOL"
)
background <- as.vector(na.omit(backgroundNames$ENTREZID))

# Add the ones that were duplicated in the select function (basically all NAs), you can find and add them manually
xx         <- as.list(org.Hs.egALIAS2EG)
dups       <- backgroundNames[duplicated(backgroundNames[,"ENTREZID"]),][,"SYMBOL"]
dupsID     <- as.vector(unlist(xx[dups]))
background <- c(background, dupsID)

allGroupsGO <- list()

# Remove or keep certain subgroups
clusts <- as.character(m.kmeans[, 5])
colors <- unique(clusts)



for(color in colors) {
  subGroupName   <- color
  subGroup       <- rownames(subset(m.kmeans, m.kmeans[, 5] == color))
  subGroupLength <- length(subGroup)
  subGroupList   <- select(org.Hs.eg.db,
                           keys    = subGroup,
                           columns = c("ENTREZID"),
                           keytype = "SYMBOL"
  )
  subGroupIDs    <- as.vector(na.omit(subGroupList$ENTREZID))
  params         <- new("GOHyperGParams", 
                        geneIds         = subGroupIDs,
                        universeGeneIds = background,
                        annotation      = "org.Hs.eg.db",
                        ontology        = "BP",
                        pvalueCutoff    = 1,
                        conditional     = TRUE, 
                        testDirection   = "over"
  )
  hgOver                      <- hyperGTest(params)
  tmp                         <- summary(hgOver)
  tmp$Length                  <- subGroupLength
  allGroupsGO[[subGroupName]] <- tmp
}

allGroupsGO <- lapply(allGroupsGO, function(x) { row.names(x) <- x[,"Term"]; x})

# Calculate the FDR of the resulting list
allGroupsGO <- lapply(allGroupsGO, CalculateFDRList)

# Filter the GOstats output
# Count per GO term more than 5% of "module" size
allGroupsGO_red <- lapply(allGroupsGO, function(x) x[x$Count > (x$Length)/20, ])
# GO term size less than 500 genes
allGroupsGO_red <- lapply(allGroupsGO, function(x) x[x$Size < 500, ])
##OddsRatio bigger than 2
allGroupsGO_red <- lapply(allGroupsGO_red, function(x) x[x$OddsRatio > 1.5, ])
# Select certain cutoff of FDR
allGroupsGO_red <- lapply(allGroupsGO_red, function(x) x[x$FDR < 0.1, ])

allSigGroupGO              <- unique(unlist(lapply(allGroupsGO_red, function(x) x[,"Term"])))
allSigGroupGOMat           <- matrix(ncol=4, nrow=length(allSigGroupGO), data=1)
rownames(allSigGroupGOMat) <- allSigGroupGO

allSigGroupGOMat[allSigGroupGO, 1] <- allGroupsGO[[1]][allSigGroupGO, "FDR"]
allSigGroupGOMat[allSigGroupGO, 2] <- allGroupsGO[[2]][allSigGroupGO, "FDR"]
allSigGroupGOMat[allSigGroupGO, 3] <- allGroupsGO[[3]][allSigGroupGO, "FDR"]
allSigGroupGOMat[allSigGroupGO, 4] <- allGroupsGO[[4]][allSigGroupGO, "FDR"]

allSigGroupGOMat <- -log10(allSigGroupGOMat)
allSigGroupGOMat[is.na(allSigGroupGOMat)] <- 0

goMat               <- as.matrix(allSigGroupGOMat)
distgoMat           <- dist(goMat)
distgoMatTree       <- hclust(distgoMat, method="complete")
col                 <- colorRampPalette(c("white", "blue", "red"))(256)
rowInd		  <- distgoMatTree$order
# rowInd_mod          <- rowInd[c(1:20, 52:59,21:51)]
pdf(file="GO_enrichment_noZero.pdf")
go_heatmap <- heatmap.2(goMat,
                        col           = col,
                        breaks        = c(seq(0, 4, length = 257)),
                        keysize       = 1,
                        margins       = c(4, 30),
                        dendrogram    = "none", 
                        trace         = "none", 
                        density.info  = "none",
                        labCol        = "",
                        Colv          = NULL,
                        Rowv          = rowInd,
                        cexCol        = 1,
                        cexRow        = 0.9,
                        ColSideColors = colors,
                        key.xlab      = "-log10FDR",
                        main          = "GO enrichment",
                        
)
dev.off()

#############################
#############################
#############################

# Do the Fisher exact test per group
allGenesSize <- dim(starvedGenes_polysomal)[1]

groupSize_1_trans <- sum(names(clustList[[4]]$cluster[clustList[[4]]$cluster == 1]) %in% transGenes)
groupSize_2_trans <- sum(names(clustList[[4]]$cluster[clustList[[4]]$cluster == 2]) %in% transGenes)
groupSize_3_trans <- sum(names(clustList[[4]]$cluster[clustList[[4]]$cluster == 3]) %in% transGenes)
groupSize_4_trans <- sum(names(clustList[[4]]$cluster[clustList[[4]]$cluster == 4]) %in% transGenes)

groupSize_1 <- allGroupsGO_red[[1]]$Length[1]
groupSize_2 <- allGroupsGO_red[[2]]$Length[1]
groupSize_3 <- allGroupsGO_red[[3]]$Length[1]
groupSize_4 <- allGroupsGO_red[[4]]$Length[1]

transGenesSize <- length(transGenes)



FE_1 <- fisher.test(matrix(c(groupSize_1_trans, 
                           transGenesSize - groupSize_1_trans,
                           groupSize_1 - groupSize_1_trans,
                           allGenesSize - groupSize_1_trans - (transGenesSize - groupSize_1_trans) - (groupSize_1 - groupSize_1_trans)),                            nrow=2, ncol=2), alternative ="t")

FE_2 <- fisher.test(matrix(c(groupSize_2_trans, 
                             transGenesSize - groupSize_2_trans,
                             groupSize_2 - groupSize_2_trans,
                             allGenesSize - groupSize_2_trans - (transGenesSize - groupSize_2_trans) - (groupSize_2 - groupSize_2_trans)),                            nrow=2, ncol=2), alternative ="t")


FE_3 <- fisher.test(matrix(c(groupSize_3_trans, 
                             transGenesSize - groupSize_3_trans,
                             groupSize_3 - groupSize_3_trans,
                             allGenesSize - groupSize_3_trans - (transGenesSize - groupSize_3_trans) - (groupSize_3 - groupSize_3_trans)),                            nrow=2, ncol=2), alternative ="t")


FE_4 <- fisher.test(matrix(c(groupSize_4_trans, 
                             transGenesSize - groupSize_4_trans,
                             groupSize_4 - groupSize_4_trans,
                             allGenesSize - groupSize_4_trans - (transGenesSize - groupSize_4_trans) - (groupSize_4 - groupSize_4_trans)),                            nrow=2, ncol=2), alternative ="t")

#################################
#################################
#################################

# Normalize per donor
normDataRed_log <- dataRed_log
for(i in unique(anotRed[, "Donor"])) {
  meanDonor_data                            <- rowMeans(dataRed_log[, anotRed[, "Donor"] == i])
  normDataRed_log[,anotRed[, "Donor"] == i] <- sweep(dataRed_log[, anotRed[, "Donor"] == i], 1, meanDonor_data)
}

normDataC <- dataC
for(i in unique(anotRed[, "Donor"])) {
  meanDonor_dataC                   <- rowMeans(dataC[, anotC[,"Donor"] == i])
  normDataC[,anotC[, "Donor"] == i] <- sweep(dataC[,anotC[, "Donor"] == i], 1, meanDonor_dataC)
}

normDataP <- dataP
for(i in unique(anotRed[, "Donor"])) {
  meanDonor_dataP                     <- rowMeans(dataP[,anotP[, "Donor"] == i])
  normDataP[,anotP[, "Donor"]  ==  i] <- sweep(dataP[,anotP[, "Donor"]  ==  i], 1, meanDonor_dataP)
}

pca    <- prcomp(t(normDataRed_log))
pcaTot <- prcomp(t(dataC))
pcaPol <- prcomp(t(dataP))

anot    <- anotRed
anotTot <- anotRed[anotRed[, "RNA"]  ==  "Total",]
anotPol <- anotRed[anotRed[, "RNA"]  ==  "Polysomal",]

group    <- paste(anot[, "Cytokine"], anot[, "Treatment"], sep = "_")
groupTot <- paste(anotTot[, "Cytokine"], anotTot[, "Treatment"], sep = "_")
groupPol <- paste(anotPol[, "Cytokine"], anotPol[, "Treatment"], sep = "_")

pdf("PCA_notNorm.pdf", useDingbats = FALSE)
pca2d(pca, group        = group, 
      show.centroids    = TRUE, 
      show.labels       = TRUE, 
      show.group.labels = TRUE
)
dev.off()
pdf("PCA_total_notNorm.pdf", useDingbats = FALSE)
pca2d(pcaTot, 
      group             = groupTot, 
      show.centroids    = TRUE,
      show.labels       = TRUE, 
      show.group.labels = TRUE
)
dev.off()
pdf("PCA_polysomal_notNorm.pdf", useDingbats = FALSE)
pca2d(pcaPol, 
      group             = groupPol, 
      show.centroids    = TRUE, 
      show.labels       = TRUE, 
      show.group.labels = TRUE
)
dev.off()

##########################################
##########################################
##########################################


pdf("QC_Heatmap.pdf", width = 8, height = 8, pointsize = 1/100)
par(cex = 2)
for(i in 1:dim(anotRed)[2]) {
  tmpData           <- normDataRed_log
  colnames(tmpData) <- anotRed[, i]
  tmpCor            <- cor(tmpData, method = "spearman")
  heatmap(tmpCor, main = colnames(anotRed)[i], symm=TRUE)
}
dev.off()
