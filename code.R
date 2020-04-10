

workingDir <- getwd()
dir.create("data")
dir.create("results")
dataDir <- file.path(workingDir, "data/")
resultsDir <- file.path(workingDir, "results/")


library(GEOquery)

my_gse <- getGEO(GEO = "GSE114626", destdir = dataDir)

getGEOSuppFiles(GEO = "GSE114626", makeDirectory = F, baseDir = dataDir)

BiocManager::install("hugene20sttranscriptcluster.db")


library(pd.hugene.2.0.st)
library(oligo)

celFiles <- list.celfiles("./data", full.names = T) #get list of .cel files
celFiles


library(Biobase)
my.targets <- read.AnnotatedDataFrame(file.path("./data", "targets.csv"), 
                                      header = T, row.names = 1,
                                      sep = ";")
my.targets


rawData <- read.celfiles(celFiles, phenoData = my.targets) #reads .cel files
rawData #ExpressionSet: both .cel and target in one object
pData(rawData)
my.targets@data$ShortName -> rownames(pData(rawData)) #change names of the samples
colnames(rawData) <- rownames((pData(rawData)))
head(rawData)


library(arrayQualityMetrics) #AQM
arrayQualityMetrics(rawData, outdir = file.path("./results", "QCDir.Raw")) 

library(ggplot2) #PCA plot
library(ggrepel)
plotPCA3 <- function(datos, labels, factor, title, scale, 
                     colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos), scale = scale)
  #plot adjustments
  dataDF <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100, 1)
  #main plot
  p1 <- ggplot(dataDF, aes(x=PC1, y =PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5, max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  #avoidind labels superposition
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels), segment.size = 0.25, size = size) +
    labs(x = c(paste("PC1", loads[1], "%")), y = c(paste("PC2", loads[2], "%"))) +
    ggtitle(paste("Principal Component Analysis for: ", title, sep = " ")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(values = colores)
}

plotPCA3(exprs(rawData), labels = my.targets@data$ShortName, factor = my.targets@data$CellType,
         title = "Raw data grouped by cell type", scale = F, size = 3,
         colores = c("red", "blue", "green", "yellow", "magenta"))

plotPCA3(exprs(rawData), labels = my.targets@data$ShortName, factor = my.targets@data$Agent,
         title = "Raw data grouped by agent", scale = F, size = 3,
         colores = c("red", "blue", "green", "yellow"))

boxplot(rawData, cex.axis = 0.5, las=2, which="all", #boxplot necesario arg which
        col = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)),
        main="Boxplot for arrays intensity: Raw Data")

plot(hclust(dist(t(exprs(rawData)))))

eset_rma <- rma(rawData) #RMA normalization (standard)

#quality control of normalized data

arrayQualityMetrics(eset_rma, outdir = file.path("./results", "QCDir.Norm"))

plotPCA3(exprs(eset_rma), labels = my.targets@data$ShortName, factor = my.targets@data$CellType,
         title = "Normalized data grouped by cell type", scale = F, size = 3,
         colores = c("red", "blue", "green", "yellow", "magenta"))

plotPCA3(exprs(eset_rma), labels = my.targets@data$ShortName, factor = my.targets@data$Agent,
         title = "Normalized data grouped by agent", scale = F, size = 3,
         colores = c("red", "blue", "green", "yellow"))

boxplot(eset_rma, cex.axis = 0.5, las=2, which="all", #boxplot necesario arg which
        col = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)),
        main="Boxplot for arrays intensity: Normalized Data")

plot(hclust(dist(t(exprs(eset_rma)))))


#NO funciona y no necesario
library(affyio)
get.celfile.dates(celFiles)

pData(eset_rma)
pData(my.targets)

targets <- read.csv2("./data/targets.csv", header = T, sep = ";")

pData(targets)

library(pvca)
pData(eset_rma)
pData(eset_rma) <- targets
#select the threshold
pct_threshold <- 0.6
#select the factor to analyze
batch.factors <- c("Agent")
#run the analysis
pvcaObj <- pvcaBatchAssess(eset_rma, batch.factors, pct_threshold)

pvcaObj$dat

#plot
bp <- barplot(pvcaObj$dat, xlab = "Effects",
              ylab = "Weighed average proportion variance",
              ylim = c(0, 1.1), col = c("mediumorchid"), las = 2,
              main = "PVCA estimation")
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.75, las = 2)
values <- round(pvcaObj$dat, 3)
text(bp, pvcaObj$dat, labels = values, pos=3, cex = 0.7)

#OK from here

#detecting most variable genes

sds <- apply(exprs(eset_rma), 1, sd) #plot of standard deviations
sds0 <- sort(sds)
plot(1:length(sds0), sds0, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))

#filtering least variable genes

library(genefilter)
library(hugene20sttranscriptcluster.db)
annotation(eset_rma) <- "hugene20sttranscriptcluster.db"
filtered <- nsFilter(eset_rma,
                     require.entrez = T, remove.dupEntrez = T,
                     var.filter = T, var.func = IQR, var.cutoff = 0.75,
                     filterByQuantile = T, feature.exclude = "^AFFX")
print(filtered$filter.log)

eset_filtered <- filtered$eset
eset_filtered
eset_rma

plot(hclust(dist(t(exprs(eset_filtered)))))

plotPCA3(exprs(eset_filtered), labels = my.targets@data$ShortName, factor = my.targets@data$Agent,
         title = "Filtered data grouped by agent", scale = F, size = 3,
         colores = c("red", "blue", "green", "yellow"))

boxplot(eset_filtered, cex.axis = 0.5, las=2, which="all", #boxplot necesario arg which
        col = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)),
        main="Boxplot for arrays intensity: Filtered Data")

#saving normalized and filterd data

write.csv(exprs(eset_rma), file = "./results/normalized.Data.csv")
write.csv(exprs(eset_filtered), file = "./results/normalized.Filtered.Data.csv")
save(eset_rma, eset_filtered, file = "./results/normalized.Data.Rda")

#defining the experimental setup: the design matrix
#linear models for microarrays method (limma)

if(!exists("eset_filtered")) load(file = "./results/normalized.Data.Rda")
library(limma)
designMat <- model.matrix(~0+Agent, pData(eset_filtered)) #design matrix
colnames(designMat) <- c("DMSO", "Meta", "Ortho", "Para")
print(designMat)

cont.matrix <- makeContrasts(PvsO = Para-Ortho,
                             PvsM = Para-Meta,
                             PvsD = Para-DMSO,
                             levels = designMat)
print(cont.matrix)

#model estimation and gene selection

library(limma)
fit <- lmFit(eset_filtered, designMat)
fit.main <- contrasts.fit(fit, cont.matrix)
fit.main <- eBayes(fit.main)
class(fit.main)


#obtaining lists of differentially expressed genes

topTab_PvsO <- topTable(fit.main, number = nrow(fit.main),
                               coef="PvsO", adjust="fdr")
head(topTab_PvsO)

topTab_PvsM <- topTable(fit.main, number = nrow(fit.main),
                             coef="PvsM", adjust="fdr")
head(topTab_PvsM)

topTab_PvsD <- topTable(fit.main, number = nrow(fit.main),
                       coef="PvsD", adjust="fdr")
head(topTab_PvsD)

#gene annotation

annotatedTopTable <- function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab <- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}

topAnnotated_PvsO <- annotatedTopTable(topTab_PvsO,
                                            anotPackage = "hugene20sttranscriptcluster.db")
write.csv(topAnnotated_PvsO, file="./results/topAnnotated_PvsO.csv")

head(topAnnotated_PvsO, 5)

topAnnotated_PvsM <- annotatedTopTable(topTab_PvsM,
                                              anotPackage = "hugene20sttranscriptcluster.db")
write.csv(topAnnotated_PvsM, file="./results/topAnnotated_PvsM.csv")

head(topAnnotated_PvsM, 5)

topAnnotated_PvsD <- annotatedTopTable(topTab_PvsD,
                                      anotPackage = "hugene20sttranscriptcluster.db")
write.csv(topAnnotated_PvsD, file="./results/topAnnotated_PvsD.csv")

head(topAnnotated_PvsD, 5)

#visualizing differential expression

library(hugene20sttranscriptcluster.db)
geneSymbols <- select(hugene20sttranscriptcluster.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS <- geneSymbols$SYMBOL
volcanoplot(fit.main, coef = 1, highlight = 4, names = SYMBOLS,
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))

#multiple comparisons

library(limma)
res <- decideTests(fit.main, method = "separate", adjust.method = "fdr", p.value = 0.1, lfc = 1)
sum.res.rows <- apply(abs(res),1,sum)
res.selected <- res[sum.res.rows!=0,]
print(summary(res))

vennDiagram(res.selected[,1:3], cex=0.8)
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC >1")

#heatmaps

probesInHeatmap <- rownames(res.selected)
HMdata <- exprs(eset_filtered)[rownames(exprs(eset_filtered)) %in% probesInHeatmap,]
geneSymbols <- select(hugene20sttranscriptcluster.db, rownames(HMdata), c("SYMBOL"))
SYMBOLS <- geneSymbols$SYMBOL
rownames(HMdata) <- SYMBOLS
write.csv(HMdata, file=file.path("./results/data4Heatmap.csv"))

my_palette <- colorRampPalette(c("blue", "red"))(n =299)
library(gplots)
dev.off()
heatmap.2(HMdata, Rowv = F, Colv = F,
          main = "Differentially expressed genes \n FDR < 0.1, logFC >=1.0",
          scale = "row", col = my_palette, sepcolor = "white",
          sepwidth = c(0.05,0.05), cexRow = 0.5, cexCol = 0.9,
          key = F, density.info = "histogram",
          ColSideColors = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)),
          tracecol = NULL, dendrogram = "none", srtCol = 30)

heatmap.2(HMdata, Rowv = T, Colv = T,
          main = "Differentially expressed genes \n FDR < 0.1, logFC >=1",
          scale = "row", col = my_palette, sepcolor = "white",
          sepwidth = c(0.05,0.05), cexRow = 0.5, cexCol = 0.9,
          key = F, density.info = "histogram",
          ColSideColors = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)),
          tracecol = NULL, dendrogram = "both", srtCol = 30)

#biological significance of resultsv (reactomePA)
library(org.Hs.eg.db)
listOfTables <- list(PvsO = topTab_PvsO,
                     PvsM = topTab_PvsM,
                     PvsD = topTab_PvsD)
listOfSelected <- list()
for(i in 1:length(listOfTables)){
  #select the toptable
  topTab <- listOfTables[[i]]
  #select the genes to be icluded in the analysis
  whichGenes <- topTab["adj.P.Val"]<0.15
  selectedIDs <- rownames(topTab)[whichGenes]
  #convert the IDD to Entrez
  EntrezIDs <- select(hugene20sttranscriptcluster.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}
sapply(listOfSelected, length)

mapped_genes2GO <- mappedkeys(org.Hs.egGO)
mapped_genes2KEGG <- mappedkeys(org.Hs.egPATH)
mapped_genes <- union(mapped_genes2GO, mapped_genes2KEGG)

library(ReactomePA)


listOfData <- listOfSelected[1:3]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

for (i in 1:length(listOfData)) {
  genesIn <- listOfData[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichPathway(gene = genesIn,
                                 pvalueCutoff = 0.05,
                                 readable = T,
                                 pAdjustMethod = "BH",
                                 organism = "human",
                                 universe = universe)
  cat("#############################")
  cat("\nComparison: ", comparison, "\n")
  print(head(enrich.result))
  
  if(length(rownames(enrich.result@result))!=0) {
    write.csv(as.data.frame(enrich.result),
              file=paste0("./results/", "ReactomePA.Results.", comparison, ".csv"),
              row.names = F)
    pdf(file=paste0("./results/", "ReactomePABarplot.", comparison, ".pdf"))
    print(barplot(enrich.result, showCategory=15, font.size=4,
                  title=paste0("Reactome Pathway Analysis for ", comparison, ". Barplot")))
    dev.off()
    
    pdf(file=paste0("./results/", "ReactomePAcnetplot.", comparison, ".pdf"))
    print(cnetplot(enrich.result, categorySize="geneNum", showCategory = 5,
                   vertex.label.cex = 0.2))
    dev.off()
  }
}

cnetplot(enrich.result, categorySize = "geneNum", showCategory = 5,
         vertex.label.cex = 0.2)











