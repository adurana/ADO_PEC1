
#Preparación del entorno

setwd("~/Desktop/Master UOC/Semestre 3/Asignaturas/157-Analisis de datos omicos/PEC1/ADO_PEC1")
workingDir <- getwd()
dir.create("data")
dir.create("results")
dataDir <- file.path(workingDir, "data/")
resultsDir <- file.path(workingDir, "results/")

#Obtención de datos

knitr::include_graphics("figuras/BusquedaGEO.png")

library(GEOquery)

my_gse <- getGEO(GEO = "GSE114626", destdir = dataDir)
getGEOSuppFiles(GEO = "GSE114626", makeDirectory = F, baseDir = dataDir)

targets <- read.csv2("./data/targets.csv", header = TRUE, sep = ";") 
knitr::kable(
  targets, booktabs = TRUE,
  caption = 'Contenido del fichero *targets.csv* para el estudio GSE114626')

#Lectura de datos

library(Biobase)
library(oligo)
library(pd.hugene.2.0.st)

celFiles <- list.celfiles("./data", full.names = T) #get list of .cel files
my.targets <- read.AnnotatedDataFrame(file.path("./data", "targets.csv"), 
                                      header = T, row.names = 1, sep = ";")
rawData <- read.celfiles(celFiles, phenoData = my.targets) #reads .cel files
pData(rawData)
my.targets@data$ShortName -> rownames(pData(rawData)) #change names of the samples
colnames(rawData) <- rownames((pData(rawData)))
rawData #ExpressionSet: both .cel and target in one object

#Control de calidad de datos crudos

library(arrayQualityMetrics)

arrayQualityMetrics(rawData, outdir = file.path("./results", "QCDir.Raw"))

knitr::include_graphics("figuras/indexRaw.png")

library(ggplot2)
library(ggrepel)

boxplot(rawData, cex.axis = 0.5, las=2, which="all", main="",
        col = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)))

plot(hclust(dist(t(exprs(rawData)))))

plotPCA <- function(datos, labels, factor, scale, 
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
  #avoiding labels superposition
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels), segment.size = 0.25, size = size) +
    labs(x = c(paste("PC1", loads[1], "%")), y = c(paste("PC2", loads[2], "%"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(values = colores)
}

plotPCA(exprs(rawData), labels = my.targets@data$ShortName, factor = my.targets@data$Agent,
        scale = F, size = 3, colores = c("red", "blue", "green", "yellow"))

#Normalización de datos

eset_rma <- rma(rawData)

#Control de calidad de datos normalizados

arrayQualityMetrics(eset_rma, outdir = file.path("./results", "QCDir.Norm"))

knitr::include_graphics("figuras/indexNorm.png")

boxplot(eset_rma, cex.axis = 0.5, las=2, which="all",
        col = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)))

plot(hclust(dist(t(exprs(eset_rma)))))

plotPCA(exprs(eset_rma), labels = my.targets@data$ShortName, factor = my.targets@data$Agent,
        scale = F, size = 3, colores = c("red", "blue", "green", "yellow"))

#Detección de los genes más variables

sds <- apply(exprs(eset_rma), 1, sd)
sds0 <- sort(sds)
plot(1:length(sds0), sds0, xlab="Genes desde el menos al más variable", ylab="Desviación estándar")
abline(v=length(sds)*c(0.9,0.95))

#Filtrado de genes

library(genefilter)
library(hugene20sttranscriptcluster.db)

annotation(eset_rma) <- "hugene20sttranscriptcluster.db"
filtered <- nsFilter(eset_rma,
                     require.entrez = T, remove.dupEntrez = T,
                     var.filter = T, var.func = IQR, var.cutoff = 0.75,
                     filterByQuantile = T, feature.exclude = "^AFFX")

eset_filtered <- filtered$eset

#Control de calidad de datos filtrados

arrayQualityMetrics(eset_filtered, outdir = file.path("./results", "QCDir.Filt"))

knitr::include_graphics("figuras/indexFilt.png")

boxplot(eset_filtered, cex.axis = 0.5, las=2, which="all",
        col = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)))

plot(hclust(dist(t(exprs(eset_filtered)))))

plotPCA(exprs(eset_filtered), labels = my.targets@data$ShortName, factor = my.targets@data$Agent,
        scale = F, size = 3, colores = c("red", "blue", "green", "yellow"))

write.csv(exprs(eset_rma), file = "./results/normalized.Data.csv")
write.csv(exprs(eset_filtered), file = "./results/filtered.Data.csv")
save(eset_rma, eset_filtered, file = "./results/normalized.filtered.Data.Rda")

#Diseño experimental

library(limma)

designMat <- model.matrix(~0+Agent, pData(eset_filtered))
colnames(designMat) <- c("DMSO", "Meta", "Ortho", "Para")

designMat

cont.matrix <- makeContrasts(PvsO = Para-Ortho,
                             PvsM = Para-Meta,
                             PvsD = Para-DMSO,
                             levels = designMat)

cont.matrix

#Modelización y selección de genes

fit <- lmFit(eset_filtered, designMat)
fit.main <- contrasts.fit(fit, cont.matrix)
fit.main <- eBayes(fit.main)

#Listado de genes diferencialmente expresados

topTab_PvsO <- topTable(fit.main, number = nrow(fit.main),
                        coef="PvsO", adjust="fdr")

topTab_PvsM <- topTable(fit.main, number = nrow(fit.main),
                        coef="PvsM", adjust="fdr")

topTab_PvsD <- topTable(fit.main, number = nrow(fit.main),
                        coef="PvsD", adjust="fdr")

head(topTab_PvsO)

head(topTab_PvsM)

head(topTab_PvsD)


#Anotación de genes

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

topAnnotated_PvsM <- annotatedTopTable(topTab_PvsM,
                                       anotPackage = "hugene20sttranscriptcluster.db")

topAnnotated_PvsD <- annotatedTopTable(topTab_PvsD,
                                       anotPackage = "hugene20sttranscriptcluster.db")

head(topAnnotated_PvsO)
head(topAnnotated_PvsM)
head(topAnnotated_PvsD)

write.csv(topAnnotated_PvsO, file="./results/topAnnotated_PvsO.csv")
write.csv(topAnnotated_PvsM, file="./results/topAnnotated_PvsM.csv")
write.csv(topAnnotated_PvsD, file="./results/topAnnotated_PvsD.csv")

#Visualización de la expresión diferencial

library(hugene20sttranscriptcluster.db)

geneSymbols <- select(hugene20sttranscriptcluster.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS <- geneSymbols$SYMBOL

volcanoplot(fit.main, coef = 1, highlight = 4, names = SYMBOLS)
abline(v=c(-1,1))

volcanoplot(fit.main, coef = 2, highlight = 4, names = SYMBOLS)
abline(v=c(-1,1))

volcanoplot(fit.main, coef = 3, highlight = 4, names = SYMBOLS)
abline(v=c(-1,1))

#Comparaciones múltiples

library(gplots)

res <- decideTests(fit.main, method = "separate", adjust.method = "fdr", p.value = 0.1, lfc = 1)
sum.res.rows <- apply(abs(res),1,sum)
res.selected <- res[sum.res.rows!=0,]
print(summary(res))

vennDiagram(res.selected[,1:3], cex=0.8)

probesInHeatmap <- rownames(res.selected)
HMdata <- exprs(eset_filtered)[rownames(exprs(eset_filtered)) %in% probesInHeatmap,]
geneSymbols <- select(hugene20sttranscriptcluster.db, rownames(HMdata), c("SYMBOL"))
SYMBOLS <- geneSymbols$SYMBOL
rownames(HMdata) <- SYMBOLS
my_palette <- colorRampPalette(c("blue", "red"))(n =299)

write.csv(HMdata, file=file.path("./results/data4Heatmap.csv"))

heatmap.2(HMdata, Rowv = T, Colv = T,
          scale = "row", col = my_palette, sepcolor = "white",
          sepwidth = c(0.05,0.05), cexRow = 0.5, cexCol = 0.9,
          key = F, density.info = "histogram",
          ColSideColors = c(rep("red", 4), rep("blue", 4), rep("green", 4), rep("yellow", 4), rep("magenta", 4)),
          tracecol = NULL, dendrogram = "both", srtCol = 30)

#Significación biológica de los resultados

library(org.Hs.eg.db)
library(ReactomePA)

listOfTables <- list(PvsO = topTab_PvsO,
                     PvsM = topTab_PvsM,
                     PvsD = topTab_PvsD)

listOfSelected <- list()

for(i in 1:length(listOfTables)){
  topTab <- listOfTables[[i]]
  whichGenes <- topTab["adj.P.Val"]<0.15
  selectedIDs <- rownames(topTab)[whichGenes]
  EntrezIDs <- select(hugene20sttranscriptcluster.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}

sapply(listOfSelected, length)

mapped_genes2GO <- mappedkeys(org.Hs.egGO)
mapped_genes2KEGG <- mappedkeys(org.Hs.egPATH)
mapped_genes <- union(mapped_genes2GO, mapped_genes2KEGG)

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
}

if(length(rownames(enrich.result@result))!=0) {
  write.csv(as.data.frame(enrich.result),
            file=paste0("./results/", "ReactomePA.Results.", comparison, ".csv"),
            row.names = F)
  pdf(file=paste0("./results/", "ReactomePABarplot.", comparison, ".pdf"))
  pdf(file=paste0("./results/", "ReactomePAcnetplot.", comparison, ".pdf"))  
}

Tab.react.PvsO <- read.csv2(file.path("./results/ReactomePA.Results.PvsO.csv"), 
                            sep = ",", header = TRUE, row.names = 1)

Tab.react.PvsO <- Tab.react.PvsO[1:5, c(1,4,5)]
knitr::kable(Tab.react.PvsO, booktabs = TRUE,
             caption = "Inicio del listado de Reactome para la comparación PvsO")

Tab.react.PvsM <- read.csv2(file.path("./results/ReactomePA.Results.PvsM.csv"), 
                            sep = ",", header = TRUE, row.names = 1)

Tab.react.PvsM <- Tab.react.PvsM[1:5, c(1,4,5)]
knitr::kable(Tab.react.PvsM, booktabs = TRUE,
             caption = "Inicio del listado de Reactome para la comparación PvsM")

Tab.react.PvsD <- read.csv2(file.path("./results/ReactomePA.Results.PvsD.csv"), 
                            sep = ",", header = TRUE, row.names = 1)

Tab.react.PvsD <- Tab.react.PvsD[1:5, c(1,4,5)]
knitr::kable(Tab.react.PvsD, booktabs = TRUE,
             caption = "Inicio del listado de Reactome para la comparación PvsD")

comparison <- comparisonsNames[1]

barplot(enrich.result, showCategory=10, font.size=10)

cnetplot(enrich.result, categorySize="geneNum", showCategory = 3,
         vertex.label.cex = 0.2)

#Resumen de resultados

listOfFiles <- dir("./results/") 
knitr::kable(
  listOfFiles, booktabs = TRUE,
  caption = "Listado de ficheros generados durante el análisis",
  col.names="ListadoFicheros"
)


