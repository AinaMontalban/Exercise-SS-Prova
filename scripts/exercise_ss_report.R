## ----setup, include=FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#install.packages("prettydoc")


## ----message=FALSE, warning=FALSE---------------------------------------------------------------
workingDir <- getwd()
data_dir <- file.path(workingDir, "data")
results_dir <- file.path(workingDir, "results")


## ----message=FALSE, warning=FALSE---------------------------------------------------------------
# load packages
library(oligo)
library(Biobase)
library(gplots)
library(ggplot2)
library(ggrepel)
library(affy)
library(GEOquery)
library(hgu133plus2cdf)
library(ggdendro)


## ----message=FALSE, warning=FALSE---------------------------------------------------------------
CELfiles_fullName <- list.celfiles(file.path(data_dir), full.names=TRUE) # save the full path


## ----message=FALSE, warning=FALSE, results='hide'-----------------------------------------------
raw_data <- read.celfiles(CELfiles_fullName)


## -----------------------------------------------------------------------------------------------
name_sample <- sampleNames(raw_data)
pData(raw_data)$name_sample <- name_sample
sampleNames <- sub(".*_", "", name_sample)
sampleNames <- sub(".CEL.gz$", "", name_sample)
sampleNames(raw_data) <- sampleNames


## ----------------------------------------------------------------------------------------------
 gse <- getGEO(filename = file.path(data_dir, "GSE23117_series_matrix.txt.gz"),
              GSEMatrix = FALSE)

## -----------------------------------------------------------------------------------------------
pheno_data <- pData(phenoData(gse))


## -----------------------------------------------------------------------------------------------
targets <- pheno_data[, c("geo_accession", "disease status:ch1")]
targets$group <- c(rep(1, 2), rep(0,4), rep(1,9))
targets$shortName <- c(rep("advanced_SS", 2), rep("control", 4), 
                       rep("control_SS", 1), rep("early_SS", 5), 
                       rep("moderate_SS", 3))


## -----------------------------------------------------------------------------------------------
# write a CSV file
write.csv(targets, file.path(data_dir, "targets.csv"), row.names = FALSE)


## -----------------------------------------------------------------------------------------------
columnDesc <-  data.frame(labelDescription= c("geo_accession",
                                              "disease status:ch1", "group", 
                                              "shortName"))
myAnnotDF <- new("AnnotatedDataFrame", data=targets, varMetadata= columnDesc)


## -----------------------------------------------------------------------------------------------
affy_object <- read.affybatch(filenames = CELfiles_fullName, 
                              phenoData = myAnnotDF)
show(affy_object)


## ----message=FALSE, warning=FALSE---------------------------------------------------------------
exp_norm <- affy::rma(affy_object)


## -----------------------------------------------------------------------------------------------
head(exprs(exp_norm))[,1:5]


## -----------------------------------------------------------------------------------------------
dist_eucl_clust <- hclust(dist(t(exprs(exp_norm))),
                               method="average")
dend <- as.dendrogram(dist_eucl_clust)


## -----------------------------------------------------------------------------------------------
dend_data <- dendro_data(dend, type = "rectangle")
lbs <- label(dend_data)$label
lst <- c()
for (ele in lbs){
    gName <- targets[targets$geo_accession==ele,"shortName"]
    lst <- c(lst, gName)
}
dend_data$labels[, "label"] <- lst


## -----------------------------------------------------------------------------------------------
p1 <- ggplot(segment(dend_data)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +  scale_y_reverse(expand = c(0.2, 0.2))

p1 <- p1 + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(),
     axis.text.y=element_blank(), axis.title.y=element_blank(),
     panel.background=element_rect(fill="white")) + labs(y="Height")

p1 <- p1 + geom_text(data = label(dend_data), 
              aes(x = x, y = y, label = label, color=label), 
              size = 3, vjust=0.5, hjust=0, show.legend = FALSE) 
p1

## ----include=FALSE------------------------------------------------------------------------------
pdf(file.path(results_dir, "images/dendrogram.pdf"))
print(p1)
dev.off()


## ----include=FALSE------------------------------------------------------------------------------
#showtext::showtext_auto()


## -----------------------------------------------------------------------------------------------
X <- t(exprs(exp_norm))
head(X)[1:4, 1:5]


## -----------------------------------------------------------------------------------------------
pca <- prcomp(X, scale. = FALSE)
Groups <- as.factor(targets$group)
shortName <- targets$shortName


## -----------------------------------------------------------------------------------------------
loads <- round(pca$sdev^2/sum(pca$sdev^2)*100,1)
pca_df <- data.frame(pca$x)


## -----------------------------------------------------------------------------------------------
p2 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Groups, label=shortName)) + 
       geom_point(aes(color=Groups), size=3) +
        scale_colour_discrete(name="Group", labels=c("Control", "SS"))

p2 <- p2 + theme_classic() + labs(x=c(paste("PC1 (",loads[1],"%)")),
                          y=c(paste("PC2 (",loads[2],"%)"))) 

p2 <- p2 + geom_text_repel(show.legend = FALSE)
p2


## ----include=FALSE------------------------------------------------------------------------------
pdf(file.path(results_dir, "images/pca.pdf"))
print(p2)
dev.off()

