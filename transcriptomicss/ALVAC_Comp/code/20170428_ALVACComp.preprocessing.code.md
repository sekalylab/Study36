---
title: ALVACComp gene-expression analysis
author: Slim Fourati
date: "02 March, 2020"
output: github_documents
---

Loading require packages

```r
suppressPackageStartupMessages(library(package = "knitr"))
suppressPackageStartupMessages(library(package = "parallel"))
suppressPackageStartupMessages(library(package = "Biobase"))
suppressPackageStartupMessages(library(package = "limma"))
suppressPackageStartupMessages(library(package = "impute"))
suppressPackageStartupMessages(library(package = "org.Hs.eg.db"))
suppressPackageStartupMessages(library(package = "pheatmap"))
suppressPackageStartupMessages(library(package = "grid"))
suppressPackageStartupMessages(library(package = "gtable"))
suppressPackageStartupMessages(library(package = "tidyverse"))
```

Set default options/variables

```r
workDir <- dirname(getwd())
opts_chunk$set(tidy = FALSE, fig.path = "../figure/")
options(stringsAsFactors  = FALSE,
        width             = 80,
        mc.cores          = detectCores() - 1,
        readr.num_columns = 0)
```

Read non-normalized matrix

```r
rawFile <- file.path(workDir,
		     "input",
		     "GA_illumina_expression.alvaccomp.matrix_non_norm.csv")
rawMat <- read_csv(file = rawFile, progress = FALSE)
```

Read arrays annotation

```r
arraysAnnotFile <- file.path(workDir,
			     "input",
			     "GA_illumina_expression.alvaccomp.metadata.csv")
arraysAnnotation <- read_csv(file = arraysAnnotFile, progress = FALSE)
# remove unused phenotypic information
arraysAnnotation <- select(arraysAnnotation,
                           -title,
                           -`source name`,
                           -organism,
                           -molecule,
                           -label,
                           -description,
                           -platform)
# remove prefix 'characteristics` of column names
setNames(arraysAnnotation, nm = gsub(pattern = "^[^:]+: (.+)$",
                                     replacement = "\\1",
                                     names(arraysAnnotation))) ->
  arraysAnnotation
```

Read features annotation

```r
featuresAnnotFile <- file.path(workDir,
			       "input/Illumina_HumanHT12_V4.rheMac3.chip")
featuresAnnotation <- read_tsv(file = featuresAnnotFile, progress = FALSE) %>%
    as.data.frame()
rownames(featuresAnnotation) <- featuresAnnotation$IlmnID
```

Create non-normalized ExpressionSet

```r
# format raw matrix
rNames <- rawMat$"ID_REF"
rawMat <- rawMat[, -grep(pattern = "ID_REF|Detection Pval",
                         colnames(rawMat))]
rawMat <- as.matrix(rawMat)
rownames(rawMat) <- rNames
# format phenodata
arraysAnnotation <- as.data.frame(arraysAnnotation)
rownames(arraysAnnotation) <- arraysAnnotation$"Sample name"
arraysAnnotation <- arraysAnnotation[colnames(rawMat), ]
# format feature annotation
featuresAnnotation <- as.data.frame(featuresAnnotation)
featuresAnnotation <- featuresAnnotation[rownames(rawMat), ]
# create ExpressionSet
esetRaw <- ExpressionSet(assayData   = rawMat,
                         phenoData   = AnnotatedDataFrame(arraysAnnotation),
                         featureData = AnnotatedDataFrame(featuresAnnotation))
# save raw ExpressionSet
save(esetRaw, file = file.path(workDir, "output/alvaccomp.esetRaw.RData"))
```

Quality control: kernel densities plot

```r
intensities <- as.vector(exprs(esetRaw))
# log2 transform (replace intensities < 1 by 1 to prevent -Inf)
intensities <- log2(pmax(intensities, 1))
arrays <- rep(sampleNames(esetRaw), each = nrow(esetRaw))
plotDF <- data.frame(intensities, arrays)
ggplot(data    = plotDF,
       mapping = aes(x     = intensities,
		     color = factor(arrays))) +
  scale_color_discrete(name = "scanning note") +
  stat_density(geom = "path", position = "identity") +
  labs(x     = "log2 intensity",
       title = "density plot [log2 raw intensities]") +
  theme_bw() +
  theme(legend.position = "none")
```

![plot of chunk density-plot](../figure/density-plot-1.png)

Quality control: multidimentional scaling plot

```r
# create tempory ExpressionSet and normalize the data
esetTemp <- esetRaw
rawMat <- exprs(esetTemp)
normMat <- normalizeBetweenArrays(rawMat, method = "quantile")
normMat <- log2(normMat)
# generate a multidimensional scaling plot
distMat <- dist(t(normMat))
mds <- cmdscale(d = distMat, k = ncol(normMat) - 1, eig = TRUE)
mdsStDev <- apply(mds$points, MARGIN = 2, FUN = sd)
mdsPercVar <- round((mdsStDev^2)/sum(mdsStDev^2) * 100)
tecRep <- grep(pattern = "_rep1", esetRaw$"Sample name", value = TRUE) %>%
  gsub(pattern = "_rep1", replacement = "")
dupArrays <- gsub(pattern = "_rep1", replacement = "", esetRaw$"Sample name")
dupArrays[!dupArrays %in% tecRep] <- NA
plotDF <- data.frame(V1               = mds$points[, 1],
                     V2               = mds$points[, 2],
                     `technical rep.` = dupArrays,
                     check.names      = FALSE)
ggplot(data = plotDF,
       mapping = aes(x = V1, y = V2, color = `technical rep.`)) +
  geom_point(size = 3) +
    labs(x     = paste0("1st dimension (", mdsPercVar[1], "%)"),
	 y     = paste0("2nd dimension (", mdsPercVar[2], "%)"),
	 title = "Multidimentional scaling plot") +
    theme_bw() +
    theme(legend.key = element_blank())
```

![plot of chunk mds-plot](../figure/mds-plot-1.png)

Normalizing raw expression

```r
eset <- esetRaw
# order esetRaw by idat file name and features by ProbeID
eset <- eset[order(as.numeric(fData(eset)$ProbeID)),
             order(eset$"idat file")]
# impute missing intensities (intensities = 0)
rawMat <- exprs(eset)
rawMat[rawMat == 0] <- NA
suppressWarnings(capture.output(rawMat <- impute.knn(data = rawMat)$data,
                                file = "/dev/null"))
exprs(eset) <- rawMat
# quantile normalized and log2 transform expression
normMat <- normalizeBetweenArrays(exprs(eset), method = "quantile")
# surrogate remplacement
normMat[normMat < 2^0.1] <- 2^0.1
normMat <- log2(normMat)
exprs(eset) <- normMat
# merge techinical replicates
exprsMat <- exprs(eset) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(`Sample name`, value, -rowname) %>%
  mutate(`Sample name` = gsub(pattern     = "_rep1",
                              replacement = "", 
                              `Sample name`)) %>%
  group_by(`Sample name`, rowname) %>%
  summarise(mu = mean(value)) %>%
  spread(`Sample name`, mu) %>%
  as.data.frame()
rownames(exprsMat) <- exprsMat$rowname
exprsMat$rowname <- NULL
eset <- eset[, colnames(exprsMat)]
exprs(eset) <- as.matrix(exprsMat)[featureNames(eset), ] 
# save normalized ExpressionSet
save(eset, file = file.path(workDir, "output/alvaccomp.eset.RData"))
```

quality control: gender check

```r
symbol2chr <- merge(as.data.frame(org.Hs.egSYMBOL),
                    as.data.frame(org.Hs.egCHR),
                    by = "gene_id")
pb2symbol <- strsplit(fData(eset)$SYMBOL, split = " /// ") %>%
  setNames(fData(eset)$IlmnID) %>%
  stack() %>%
  mutate(ind = as.vector(ind))
pb2chr <- merge(symbol2chr, pb2symbol, by.x = "symbol", by.y = "values")

pbY <- filter(pb2chr, chromosome %in% "Y") %>%
  .$ind

plotDF <- data.frame(mu = colMeans(exprs(eset)[pbY, ])) %>%
  rownames_to_column() %>%
  merge(pData(eset), by.x = "rowname", by.y = "Sample name")

ggplot(data = plotDF, mapping = aes(x = donor, y = mu)) +
  geom_jitter() +
  labs(y = "Average expression of chr.Y probes") +
  theme_bw()
```

![plot of chunk boxplot-chry](../figure/boxplot-chry-1.png)

Exploratory analysis: heatmap based on top 100 most varying transcripts

```r
bluered <- colorRampPalette(colors = c("blue", "white", "red"))
varList <- apply(exprs(eset), MARGIN = 1, FUN = var)
esetTemp <- eset[order(varList, decreasing = TRUE)[1:100], ]
esetTemp$donor <- factor(esetTemp$donor)
exploHeat <- pheatmap(mat            = exprs(esetTemp),
                      color          = bluered(100),
                      scale          = "row",
                      treeheight_row = 0,
                      annotation_col = pData(esetTemp)[,
                          c("donor",
                            "vaccination time",
                            "time after immunization")],
                      show_colnames  = FALSE,
                      show_rownames  = FALSE,
                      main           = paste0("top 100 varying probes",
                          "[Euclidean distance, complete linkage]"),
                      silent         = TRUE)
colorName <- textGrob("z-score", x = 0.5, y = 1.05, gp = gpar(fontface = "bold"))
exploHeat$gtable <- gtable_add_grob(exploHeat$gtable,
                                    colorName,
                                    t    = 3,
                                    l    = 5,
                                    b    = 5,
                                    clip = "off",
                                    name = "colorName")
grid.draw(exploHeat$gtable)
```

![plot of chunk exploratory-heatmap](../figure/exploratory-heatmap-1.png)

```r
# test association with time after immunization
tab <- table(cutree(exploHeat$tree_col, k = 2),
             esetTemp$"time after immunization")
fit <- fisher.test(tab)
print(tab)
```

```
##    
##     16 hours 24 hours 4 weeks 48 hours 72 hours pre
##   1        3        3       6       24       24   0
##   2       21       21       0        0        0   6
```

```r
print(fit)
```

```
## 
## 	Fisher's Exact Test for Count Data
## 
## data:  tab
## p-value < 2.2e-16
## alternative hypothesis: two.sided
```

```r
# test for clustering of samples by participants
participantList <- esetTemp$donor[exploHeat$tree_col$order]
y <- sum(participantList[-length(participantList)] == participantList[-1])
# derive p-value using a 1000 fold permutation test
B <- 1000
yhat <- mclapply(1:B, function(seed) {
  set.seed(seed = seed)
  participantList <- sample(esetTemp$donor)
  return(value = sum(participantList[-length(participantList)] ==
             participantList[-1]))
})
print(paste0("permutation test: p<=", max(mean(unlist(yhat) >= y), 1/B)))
```

```
## [1] "permutation test: p<=0.001"
```

Create a prevax-substracted ExpressionSet

```r
# identify complete pair of ENV-DMSO stimulated samples
esetPre <- eset[, eset$"vaccination time" %in% "pre-vaccination"]
esetBaselined <- eset[, eset$"vaccination time" != "pre-vaccination"]
exprs(esetBaselined) <- exprs(esetBaselined) - 
  exprs(esetPre[, match(esetBaselined$donor, table = esetPre$donor)])
# save baselined expression
save(esetBaselined, file = file.path(workDir, "output/alvaccomp.esetBaselined.RData"))
```







print session info

```r
sessionInfo()
```

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-apple-darwin19.2.0 (64-bit)
## Running under: macOS Catalina 10.15.3
## 
## Matrix products: default
## BLAS/LAPACK: /usr/local/Cellar/openblas/0.3.7/lib/libopenblasp-r0.3.7.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] forcats_0.4.0        stringr_1.4.0        dplyr_0.8.4         
##  [4] purrr_0.3.3          readr_1.3.1          tidyr_1.0.2         
##  [7] tibble_2.1.3         ggplot2_3.2.1        tidyverse_1.3.0     
## [10] gtable_0.3.0         pheatmap_1.0.12      org.Hs.eg.db_3.10.0 
## [13] AnnotationDbi_1.48.0 IRanges_2.20.1       S4Vectors_0.24.1    
## [16] impute_1.60.0        limma_3.42.0         Biobase_2.46.0      
## [19] BiocGenerics_0.32.0  knitr_1.27          
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.3         lubridate_1.7.4    lattice_0.20-38    assertthat_0.2.1  
##  [5] digest_0.6.23      R6_2.4.1           cellranger_1.1.0   backports_1.1.5   
##  [9] reprex_0.3.0       RSQLite_2.2.0      evaluate_0.14      highr_0.8         
## [13] httr_1.4.1         pillar_1.4.3       rlang_0.4.4        lazyeval_0.2.2    
## [17] readxl_1.3.1       rstudioapi_0.10    blob_1.2.1         labeling_0.3      
## [21] bit_1.1-15.1       munsell_0.5.0      broom_0.5.4        compiler_3.6.2    
## [25] modelr_0.1.5       xfun_0.12          pkgconfig_2.0.3    tidyselect_1.0.0  
## [29] fansi_0.4.1        withr_2.1.2        crayon_1.3.4       dbplyr_1.4.2      
## [33] nlme_3.1-143       jsonlite_1.6.1     lifecycle_0.1.0    DBI_1.1.0         
## [37] magrittr_1.5       scales_1.1.0       cli_2.0.1          stringi_1.4.5     
## [41] farver_2.0.3       fs_1.3.1           xml2_1.2.2         ellipsis_0.3.0    
## [45] generics_0.0.2     vctrs_0.2.2        RColorBrewer_1.1-2 tools_3.6.2       
## [49] bit64_0.9-7        glue_1.3.1         hms_0.5.3          colorspace_1.4-1  
## [53] rvest_0.3.5        memoise_1.1.0      haven_2.2.0
```

