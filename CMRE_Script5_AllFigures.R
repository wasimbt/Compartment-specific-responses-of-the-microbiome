#Copyright (c) 2023 University of Bern, Bern, Switzerland. 
#Author: Wasimuddin, Klaus Schlaeppi, Alban Ramette

# library(phyloseq); packageVersion("dada2") # 
# [1] ‘1.24.0’
# library(microbiome); packageVersion("microbiome") # 
# [1] ‘1.18.0’
# library(RColorBrewer); packageVersion("RColorBrewer") # 
# [1] ‘1.1.3’
# library(ggpubr); packageVersion("ggpubr") # 
# [1] ‘0.4.0’
# library(dplyr); packageVersion("dplyr") # 
# [1] ‘1.0.9’
# library(SpiecEasi); packageVersion("SpiecEasi") #
# [1] ‘1.1.2’
# library(network); packageVersion("network")
# [1] ‘1.17.2’
# library(intergraph); packageVersion("intergraph")
# [1] ‘2.0.2’
# library(ggnet); packageVersion("ggnet")
# [1] ‘0.1.0’
# library(igraph); packageVersion("igraph")
# [1] ‘1.3.2’
# library(sna); packageVersion("sna")
# [1] ‘2.7’
# library(poweRlaw); packageVersion("poweRlaw")
# [1] ‘0.70.6’
# library(tidyr); packageVersion("tidyr")
# ## [1] '0.8.1'
# library(ggplot2); packageVersion("ggplot2") 
# ## [1] '3.1.1'
# library(phyloseq); packageVersion("phyloseq") 
# ## [1] ‘1.26.1’
# library("data.table"); packageVersion("data.table")
# ## [1] ‘1.11.6’
# library(devtools); packageVersion("devtools")
# ## [1] ‘2.4.3’
# library("scales"); packageVersion("scales")
# ##[1] ‘1.2.0’


######################################################################################################
###########################----------------All plot scripts-------------------########################

######################################################################################################
###########################----------------Figure1 scripts--------------------########################
---
  title: "Figure 1 sample types"
author: "Klaus Schlaeppi"
date: "`r Sys.Date()`"
output:
  pdf_document:
  
  html_document:
  
fig_caption: yes
fig_height: 20
fig_width: 15
toc: yes
toc_float: yes
editor_options: 
  chunk_output_type: console
---
  
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=T)
```

```{r,include=FALSE}
rm(list=ls())
```


```{r libraries, echo=F, message=F, warning=F}
library(ggplot2); packageVersion("ggplot2")   # for creating the final graph at the end of the pipeline
library(phyloseq); packageVersion("phyloseq") # necessary for all analysis
library(vegan); packageVersion("vegan")       # necessary for distance calcs
library(gplots); packageVersion("gplots")     # heatmap.2()
library("VennDiagram"); packageVersion("VennDiagram") # 5-way venn diagramm
```


# data upload

```{r, include=FALSE, echo=F, message=F, warning=F}

# load data
CMRE9_r8100 <- readRDS("input/CMRE9_r8100.RDS") 
levels(sample_data(CMRE9_r8100)$Sampletype) <- c("Mouse", "Root", "Sediment", "Soil", "Water")
sample_data(CMRE9_r8100)$Sampletype <- factor(sample_data(CMRE9_r8100)$Sampletype, 
                                              levels = c("Sediment", "Water", "Soil", "Root", "Mouse")) # ordering samples
sample_data(CMRE9_r8100)$Treatment <- factor(sample_data(CMRE9_r8100)$Treatment, 
                                             levels = c("Ctr", "As", "Bx", "Tb" )) # ordering samples
```

# heatmap filtered data (1k abundant ASVs)  

## data normalization
# 1.) express as %RA.   
# 2.) filter >0.01% %RA.  
# 3.) presence/absence.   


```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}

### Data normalization
######################
## express as relative abundance, RA
dat <- t(otu_table(CMRE9_r8100))
dat.ra <- t( t(dat) / colSums(dat) ) * 100
# dim(dat.ra)

## filter >0.01% %RA
dat.ra.fi <- dat.ra[rowMeans(dat.ra) >= 0.01,] 
# dim(dat.ra.fi)

## write back into phyloseq object
CMRE9_r8100_ra_fi <- CMRE9_r8100
otu_table(CMRE9_r8100_ra_fi) <- dat.ra.fi
# CMRE9_r8100_ra_fi

## presence/absence
CMRE9_r8100_ra_fi_pa <- CMRE9_r8100_ra_fi
otu_table(CMRE9_r8100_ra_fi_pa)[otu_table(CMRE9_r8100_ra_fi_pa) > 0] <- 1
# otu_table(CMRE9_r8100_ra_fi_pa)[1:10,1:10]
CMRE9_r8100_ra_fi_pa

```

\newpage 

## heatmap
# 1) Jaccard distance (presence/absence of ASVs). 
# 2) clustering method "average". 

### default clustering    

```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}

### plotting with heatmap.2() function 
######################
CMRE9_r8100_ra_fi_pa <- t(CMRE9_r8100_ra_fi_pa) # requires rows as samples, columns as ASVs

### default clustering
# pdf("Fig_1B_heatmap_default_clustering.pdf", width=25/cm(1), height=15/cm(1), pointsize=5, fonts="Helvetica")
for_dendros <- heatmap.2(otu_table(CMRE9_r8100_ra_fi_pa), 
                         scale="none", trace="none", density.info="none",
                         labRow=NULL, cexCol=0.9,
                         main="default clustering", col=c("white", "black"),
                         margins = c(8,8), lhei=c(2,6),
                         hclust=function(x)hclust(d=x, method="average"),
                         dist=function(x){vegdist(x, method="jaccard")} )
# dev.off()
```

\newpage 

### manual clustering

```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}
## prepare a dendrogramm for the samples
distances_samples <- vegdist(otu_table(CMRE9_r8100_ra_fi_pa), method="jaccard")
dendro_samples <- hclust(distances_samples, method="average")
# plot(dendro_samples)

## prepare a dendrogramm for the ASVs
distance_ASVs <- vegdist(t(otu_table(CMRE9_r8100_ra_fi_pa)), method="jaccard")
dendro_ASVs <- hclust(distance_ASVs, method="average")
# plot(dendro_ASVs)

heatmap.2(otu_table(CMRE9_r8100_ra_fi_pa), 
          scale="none", trace="none", density.info="none",
          labRow=NULL, cexCol=0.9,
          main="manual clustering", col=c("white", "black"),
          margins = c(8,8), lhei=c(2,6),
          hclust=function(x)hclust(d=x, method="average"),
          dist=function(x){vegdist(x, method="jaccard")},
          Colv=as.dendrogram(dendro_ASVs),
          Rowv=as.dendrogram(dendro_samples) )

```

\newpage 

### sorting sample clustering  
(branches: Water, Sediment, Soil, Root to Mouse)  

```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}

# samples are in rows

### extract their order in the dendrogram 
# for_dendros$rowInd
sample_order <- sample_data(CMRE9_r8100_ra_fi_pa)[for_dendros$rowInd, 6] # take sample IDs and sort them as they are in the dendrogram
# # check sample_order
# head(sample_order)
# sample_order[260:270,]
# sample_order[c(1:5,260:266),]
# test <- sample_order[sample_order$Sampletype=="Mouse",]
# test[c(1:5,260:266),]

## new sorting for samples: Water, Sediment, Soil, Root to Mouse (but keep within-leaf order of the dendrogram)
sample_order_sorted <- rbind(sample_order[sample_order$Sampletype=="Water",],
                             sample_order[sample_order$Sampletype=="Sediment",],
                             sample_order[sample_order$Sampletype=="Soil",],
                             sample_order[sample_order$Sampletype=="Root",],
                             sample_order[sample_order$Sampletype=="Mouse",] )
# head(sample_order_sorted)

## apply new sample order to previously defined dengdrogramm dendro_samples
dendro_samples_reordered <- reorder(dendro_samples, wts=order(match(rownames(sample_order_sorted),                              rownames(otu_table(CMRE9_r8100_ra_fi_pa)))))
# plot(dendro_samples_reordered)

### heatmap with newly sorted dendrogram
heatmap.2(otu_table(CMRE9_r8100_ra_fi_pa), 
          scale="none", trace="none", density.info="none",
          labRow=NULL, cexCol=0.9,
          main="sample clustering sorted", col=c("white", "black"),
          margins = c(8,8), lhei=c(2,6),
          hclust=function(x)hclust(d=x, method="average"),
          dist=function(x){vegdist(x, method="jaccard")},
          Rowv=as.dendrogram(dendro_samples_reordered) )

```

\newpage 

### sorting ASV clustering  
(branches: blocks of Water-, Sediment-, Soil-, Root- to Mouse-ASVs)  

```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}

# ASVs are in columns

### extract their order in their dendrogram 
# for_dendros$rowInd
ASV_table_ordered <- otu_table(CMRE9_r8100_ra_fi_pa)[,for_dendros$colInd] # take colnames of ASV IDs and sort them as they are in the dendrogram
# dim(ASV_table_ordered)

### identifying ASV memberships in the different clusters...
## define a temporary design object
design <- sample_data(CMRE9_r8100_ra_fi_pa)[,6] 
# head(design)

### ...based on their summed presence/absence in their groups by sample types
## the identification of ASV memberships per clusters is done manually - eye-balling
# range of mouse-ASVs
ASV_table_ordered_Mouse <- ASV_table_ordered[rownames(design[design$Sampletype=="Mouse"]),]
# dim(ASV_table_ordered_Mouse)
ASV_table_ordered_Mouse_SUM <- colSums(ASV_table_ordered_Mouse)
# plot(ASV_table_ordered_Mouse_SUM)
# plot(ASV_table_ordered_Mouse_SUM[1:150])
# plot(ASV_table_ordered_Mouse_SUM[1:126])
# plot(ASV_table_ordered_Mouse_SUM[124:150])
Mouse_range <- 1:127

# range of water-ASVs
ASV_table_ordered_water <- ASV_table_ordered[rownames(design[design$Sampletype=="Water"]),]
# dim(ASV_table_ordered_water)
ASV_table_ordered_water_SUM <- colSums(ASV_table_ordered_water)
# plot(ASV_table_ordered_water_SUM)
# plot(ASV_table_ordered_water_SUM[100:300])
# plot(ASV_table_ordered_water_SUM[120:220])
# plot(ASV_table_ordered_water_SUM[120:200])
Water_range <- 128:198

# range of water/sediment and sediment-ASVs
ASV_table_ordered_Sedi <- ASV_table_ordered[rownames(design[design$Sampletype=="Sediment"]),]
# dim(ASV_table_ordered_Sedi)
ASV_table_ordered_Sedi_SUM <- colSums(ASV_table_ordered_Sedi)
# plot(ASV_table_ordered_Sedi_SUM[120:580])
# plot(ASV_table_ordered_Sedi_SUM[300:580]) #search break point before p80
# plot(ASV_table_ordered_Sedi_SUM[320:340])
WaterSedi_range <- 199:320
# plot(ASV_table_ordered_Sedi_SUM[500:600])
# plot(ASV_table_ordered_Sedi_SUM[553:580])
Sedi_range <- 321:553

# range of root- and root/soil-ASVs
ASV_table_ordered_root <- ASV_table_ordered[rownames(design[design$Sampletype=="Root"]),]
# dim(ASV_table_ordered_root)
ASV_table_ordered_root_SUM <- colSums(ASV_table_ordered_root)
# plot(ASV_table_ordered_root_SUM)
# plot(ASV_table_ordered_root_SUM[500:800])
# plot(ASV_table_ordered_root_SUM[732:750])
Root_range <- 554:732
# plot(ASV_table_ordered_root_SUM[700:900])
# plot(ASV_table_ordered_root_SUM[880:900])
# plot(ASV_table_ordered_root_SUM[896:900])
SoilRoot_range <- 733:896

# rest are range of soil-ASVs
Soil_range <- 897:1174

### combine ranges (use rev to flip the tree within the leaves)
ASV_order_sorted <- c(rev(Water_range), rev(WaterSedi_range), rev(Sedi_range), rev(Soil_range), rev(SoilRoot_range), rev(Root_range), rev(Mouse_range))
# get ASV-IDs from position indexes
ASV_ID_order_sorted <- colnames(ASV_table_ordered)[ASV_order_sorted]

## apply new ASV order to previously defined dendrogramm dendro_ASVs
dendro_ASVs_reordered <- reorder(dendro_ASVs, wts=order(match(ASV_ID_order_sorted, colnames(otu_table(CMRE9_r8100_ra_fi_pa)))))
# plot(dendro_ASVs_reordered)

### heatmap with newly sorted dendrogram for ASVs
heatmap.2(otu_table(CMRE9_r8100_ra_fi_pa), 
          scale="none", trace="none", density.info="none",
          labRow=NULL, cexCol=0.9,
          main="both clusterings sorted", col=c("white", "black"),
          margins = c(8,8), lhei=c(2,6),
          hclust=function(x)hclust(d=x, method="average"),
          dist=function(x){vegdist(x, method="jaccard")},
          Rowv=as.dendrogram(dendro_samples_reordered),
          Colv=as.dendrogram(dendro_ASVs_reordered) )

```

\newpage 

# Fig 1B   

```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}

### prepare color vector for the different sample types
design$cols <- design$Sampletype
levels(design$cols) <- c("slategray4","steelblue3","salmon4","darkolivegreen4","tan3")
# head(design)

### heatmap with newly defined color scheme
heatmap.2(otu_table(CMRE9_r8100_ra_fi_pa), 
          scale="none", trace="none", density.info="none",
          labRow=NULL, cexCol=0.9,
          main="both clusterings sorted, color", col=c("white", "black"),
          margins = c(8,8), lhei=c(2,6),
          hclust=function(x)hclust(d=x, method="average"),
          dist=function(x){vegdist(x, method="jaccard")},
          Rowv=as.dendrogram(dendro_samples_reordered),
          Colv=as.dendrogram(dendro_ASVs_reordered),
          RowSideColors=as.vector(design$cols) )


pdf("Fig_1B_heatmap_sorted_clustering.pdf", width=25/cm(1), height=15/cm(1), pointsize=5, fonts="Helvetica")
heatmap.2(otu_table(CMRE9_r8100_ra_fi_pa), 
          scale="none", trace="none", density.info="none",
          labRow=NULL, cexCol=0.9,
          main="both clusterings sorted, color", col=c("white", "black"),
          margins = c(8,8), lhei=c(2,6),
          hclust=function(x)hclust(d=x, method="average"),
          dist=function(x){vegdist(x, method="jaccard")},
          Rowv=as.dendrogram(dendro_samples_reordered),
          Colv=as.dendrogram(dendro_ASVs_reordered),
          RowSideColors=as.vector(design$cols) )

dev.off()
```

\newpage 

# Fig S5, heatmap all data (all 22k ASVs)  

```{r, include=FALSE, echo=F, message=F, warning=F}
CMRE9_r8100 # same as above

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 22457 taxa and 1054 samples ]
# sample_data() Sample Data:       [ 1054 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 22457 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 22457 tips and 22456 internal nodes ]
```

## data normalization
# 1.) express as %RA.   
# 2.) NO filtering >0.01% %RA.  
# 3.) presence/absence.   

```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}
### Data normalization
######################
# ## express as relative abundance, RA
# dat <- t(otu_table(CMRE9_r8100))
# dat.ra <- t( t(dat) / colSums(dat) ) * 100
# # dim(dat.ra)

## write back into phyloseq object
CMRE9_r8100_ra <- CMRE9_r8100
otu_table(CMRE9_r8100_ra) <- dat.ra
# CMRE9_r8100_ra_fi

## presence/absence
CMRE9_r8100_ra_pa <- CMRE9_r8100_ra
otu_table(CMRE9_r8100_ra_pa)[otu_table(CMRE9_r8100_ra_pa) > 0] <- 1
# otu_table(CMRE9_r8100_ra_fi_pa)[1:10,1:10]
# CMRE9_r8100_ra_pa
```

## heatmap
# 1) Jaccard distance (presence/absence of ASVs). 
# 2) clustering method "average". 

### default clustering    

```{r, echo=T, message=F, warning=F, fig.width=7, fig.height=5}

### plotting with heatmap.2() function 
######################
CMRE9_r8100_ra_pa <- t(CMRE9_r8100_ra_pa) # requires rows as samples, columns as ASVs

```

\newpage 

# Fig 1C, overlap analysis  
(detection of all 22k ASVs in the 5 compartments)  

```{r, echo=FALSE, message=F, warning=F, fig.width=5, fig.height=5}

# CMRE9_r8100 # same as above
# we subset the full object for each compartment, to determine the number of ASVs in each compartment 

### Subsets by compartment
Water_CMRE9_r8100 <- prune_samples(CMRE9_r8100@sam_data$Sampletype=="Water", CMRE9_r8100)
Sediment_CMRE9_r8100 <- prune_samples(CMRE9_r8100@sam_data$Sampletype=="Sediment", CMRE9_r8100)
Soil_CMRE9_r8100 <- prune_samples(CMRE9_r8100@sam_data$Sampletype=="Soil", CMRE9_r8100)
Root_CMRE9_r8100 <- prune_samples(CMRE9_r8100@sam_data$Sampletype=="Root", CMRE9_r8100)
Mouse_CMRE9_r8100 <- prune_samples(CMRE9_r8100@sam_data$Sampletype=="Mouse", CMRE9_r8100)

### per compartment, keep only ASVs with >= 1 counds (remove ASVs = 0)
Water_CMRE9_r8100 <- prune_taxa(taxa_sums(Water_CMRE9_r8100) > 0, Water_CMRE9_r8100) 
Sediment_CMRE9_r8100 <- prune_taxa(taxa_sums(Sediment_CMRE9_r8100) > 0, Sediment_CMRE9_r8100) 
Soil_CMRE9_r8100 <- prune_taxa(taxa_sums(Soil_CMRE9_r8100) > 0, Soil_CMRE9_r8100) 
Root_CMRE9_r8100 <- prune_taxa(taxa_sums(Root_CMRE9_r8100) > 0, Root_CMRE9_r8100) 
Mouse_CMRE9_r8100 <- prune_taxa(taxa_sums(Mouse_CMRE9_r8100) > 0, Mouse_CMRE9_r8100) 

### Get ASV-IDs of each compartment
Water_ASVs <- colnames(otu_table(Water_CMRE9_r8100))
length(Water_ASVs)
Sediment_ASVs <- colnames(otu_table(Sediment_CMRE9_r8100))
length(Sediment_ASVs)
Soil_ASVs <- colnames(otu_table(Soil_CMRE9_r8100))
length(Soil_ASVs)
Root_ASVs <- colnames(otu_table(Root_CMRE9_r8100))
length(Root_ASVs)
Mouse_ASVs <- colnames(otu_table(Mouse_CMRE9_r8100))
length(Mouse_ASVs)

length( unique( c(Water_ASVs, Sediment_ASVs, Soil_ASVs, Root_ASVs, Mouse_ASVs) ) )
# 22457 unique ASVs


### Get the intersections

## 2-way intersections
shared_Mouse_Root     <- intersect(Mouse_ASVs, Root_ASVs)
shared_Mouse_Sediment <- intersect(Mouse_ASVs, Sediment_ASVs)
shared_Mouse_Soil     <- intersect(Mouse_ASVs, Soil_ASVs)
shared_Mouse_Water    <- intersect(Mouse_ASVs, Water_ASVs)
shared_Root_Sediment  <- intersect(Root_ASVs, Sediment_ASVs)
shared_Root_Soil      <- intersect(Root_ASVs, Soil_ASVs)
shared_Root_Water     <- intersect(Root_ASVs, Water_ASVs)
shared_Sediment_Soil  <- intersect(Sediment_ASVs, Soil_ASVs)
shared_Sediment_Water <- intersect(Sediment_ASVs, Water_ASVs)
shared_Soil_Water     <- intersect(Soil_ASVs, Water_ASVs)

## 3-way intersections
shared_Mouse_Root_Sediment  <- intersect(intersect(Mouse_ASVs, Root_ASVs), Sediment_ASVs)
shared_Mouse_Root_Soil      <- intersect(intersect(Mouse_ASVs, Root_ASVs), Soil_ASVs)
shared_Mouse_Root_Water     <- intersect(intersect(Mouse_ASVs, Root_ASVs), Water_ASVs)
shared_Mouse_Sediment_Soil  <- intersect(intersect(Mouse_ASVs, Sediment_ASVs), Soil_ASVs)
shared_Mouse_Sediment_Water <- intersect(intersect(Mouse_ASVs, Sediment_ASVs), Water_ASVs)
shared_Mouse_Soil_Water     <- intersect(intersect(Mouse_ASVs, Soil_ASVs), Water_ASVs)
shared_Root_Sediment_Soil   <- intersect(intersect(Root_ASVs, Sediment_ASVs), Soil_ASVs)
shared_Root_Sediment_Water  <- intersect(intersect(Root_ASVs, Sediment_ASVs), Water_ASVs)
shared_Root_Soil_Water      <- intersect(intersect(Root_ASVs, Soil_ASVs), Water_ASVs)
shared_Sediment_Soil_Water  <- intersect(intersect(Sediment_ASVs, Soil_ASVs), Water_ASVs)

## 4-way intersections
shared_Mouse_Root_Sediment_Soil   <- Reduce(intersect, list(Mouse_ASVs, Root_ASVs, Sediment_ASVs, Soil_ASVs))
shared_Mouse_Root_Sediment_Water  <- Reduce(intersect, list(Mouse_ASVs, Root_ASVs, Sediment_ASVs, Water_ASVs))
shared_Mouse_Root_Soil_Water      <- Reduce(intersect, list(Mouse_ASVs, Root_ASVs, Soil_ASVs, Water_ASVs))
shared_Mouse_Sediment_Soil_Water  <- Reduce(intersect, list(Mouse_ASVs, Sediment_ASVs, Soil_ASVs, Water_ASVs))
shared_Root_Sediment_Soil_Water   <- Reduce(intersect, list(Root_ASVs, Sediment_ASVs, Soil_ASVs, Water_ASVs))

## 5-way intersection
shared_Mouse_Root_Sediment_Soil_Water <- Reduce(intersect, list(Mouse_ASVs, Root_ASVs, Sediment_ASVs, Soil_ASVs, Water_ASVs))


### Draw the 5-way venn diagram
# library("VennDiagram")
# ?draw.quintuple.venn
# area1 = Water
# area2 = Sediment
# area3 = Soil
# area4 = Root
# area5 = Mouse
# n12 = Water and Sediment
# n13 = ...

grid.newpage()
draw.quintuple.venn(area1=length(Water_ASVs), area2=length(Sediment_ASVs), area3=length(Soil_ASVs),
                    area4=length(Root_ASVs), area5=length(Mouse_ASVs),
                    n12=length(shared_Sediment_Water), n13=length(shared_Soil_Water), 
                    n14=length(shared_Root_Water), n15=length(shared_Mouse_Water), 
                    n23=length(shared_Sediment_Soil), n24=length(shared_Root_Sediment), n25=length(shared_Mouse_Sediment),
                    n34=length(shared_Root_Soil), n35=length(shared_Mouse_Soil), n45=length(shared_Mouse_Root),
                    n123=length(shared_Sediment_Soil_Water), n124=length(shared_Root_Sediment_Water), 
                    n125=length(shared_Mouse_Sediment_Water), n134=length(shared_Root_Soil_Water), 
                    n135=length(shared_Mouse_Soil_Water), n145=length(shared_Mouse_Root_Water), 
                    n234=length(shared_Root_Sediment_Soil), n235=length(shared_Mouse_Sediment_Soil), 
                    n245=length(shared_Mouse_Root_Sediment), n345=length(shared_Mouse_Root_Soil), 
                    n1234=length(shared_Root_Sediment_Soil_Water), n1235=length(shared_Mouse_Sediment_Soil_Water),
                    n1245=length(shared_Mouse_Root_Sediment_Water), n1345=length(shared_Mouse_Root_Soil_Water), 
                    n2345=length(shared_Mouse_Root_Sediment_Soil), n12345=length(shared_Mouse_Root_Sediment_Soil_Water),
                    cat.fontface=rep("plain", 5), cat.fontfamily=rep("Helvetica", 5), 
                    lty=rep("blank", 5), fill=c("steelblue3","slategray4","salmon4","darkolivegreen4","tan3"),
                    alpha=rep(0.5, 5), cex=rep(1.5), fontface=rep("bold", 31), fontfamily=rep("Helvetica", 31), 
                    cat.cex=rep(2.0))



pdf("Fig_1C_5-way_venn.pdf", width=15/cm(1), height=15/cm(1), pointsize=5, fonts="Helvetica")
grid.newpage()
draw.quintuple.venn(area1=length(Water_ASVs), area2=length(Sediment_ASVs), area3=length(Soil_ASVs),
                    area4=length(Root_ASVs), area5=length(Mouse_ASVs),
                    n12=length(shared_Sediment_Water), n13=length(shared_Soil_Water), 
                    n14=length(shared_Root_Water), n15=length(shared_Mouse_Water), 
                    n23=length(shared_Sediment_Soil), n24=length(shared_Root_Sediment), n25=length(shared_Mouse_Sediment),
                    n34=length(shared_Root_Soil), n35=length(shared_Mouse_Soil), n45=length(shared_Mouse_Root),
                    n123=length(shared_Sediment_Soil_Water), n124=length(shared_Root_Sediment_Water), 
                    n125=length(shared_Mouse_Sediment_Water), n134=length(shared_Root_Soil_Water), 
                    n135=length(shared_Mouse_Soil_Water), n145=length(shared_Mouse_Root_Water), 
                    n234=length(shared_Root_Sediment_Soil), n235=length(shared_Mouse_Sediment_Soil), 
                    n245=length(shared_Mouse_Root_Sediment), n345=length(shared_Mouse_Root_Soil), 
                    n1234=length(shared_Root_Sediment_Soil_Water), n1235=length(shared_Mouse_Sediment_Soil_Water),
                    n1245=length(shared_Mouse_Root_Sediment_Water), n1345=length(shared_Mouse_Root_Soil_Water), 
                    n2345=length(shared_Mouse_Root_Sediment_Soil), n12345=length(shared_Mouse_Root_Sediment_Soil_Water),
                    cat.fontface=rep("plain", 5), cat.fontfamily=rep("Helvetica", 5), 
                    lty=rep("blank", 5), fill=c("steelblue3","slategray4","salmon4","darkolivegreen4","tan3"),
                    alpha=rep(0.5, 5), cex=rep(2), fontface=rep("bold", 31), fontfamily=rep("Helvetica", 31), 
                    cat.cex=rep(2.0))
dev.off()
```


######################################################################################################
###########################----------------Figure2 scripts--------------------########################

############### Draw the boxplot for alpha diversity
Treatment <- (CMRE5_alpha_r8100$Treatment)
Timepoint <- (CMRE5_alpha_r8100$Timepoint)
Concentration <- (CMRE5_alpha_r8100$Concentration)
Lib <- (CMRE5_alpha_r8100$Lib)
CMRE5_alpha_r8100$Sampletype <- gsub("feces","Mouse",CMRE5_alpha_r8100$Sampletype)
CMRE5_alpha_r8100$Sampletype <- gsub("roots","Root",CMRE5_alpha_r8100$Sampletype)
CMRE5_alpha_r8100$Sampletype <- gsub("sedi","Sediment",CMRE5_alpha_r8100$Sampletype)
CMRE5_alpha_r8100$Sampletype <- gsub("water","Water",CMRE5_alpha_r8100$Sampletype)
CMRE5_alpha_r8100$Sampletype <- gsub("soil","Soil",CMRE5_alpha_r8100$Sampletype)
Sampletype <- (CMRE5_alpha_r8100$Sampletype)

Shannon <- (CMRE5_alpha_r8100$Shannon)

############################---Observed final box plot----##############################
#### for treatement and compartment 

CMRE9_alpha_r8100$Sampletype<- factor(CMRE9_alpha_r8100$Sampletype, levels = c("Water", "Sediment", "Soil", "Root", "Mouse"))
CMRE9_alpha_r8100$Treatment<- factor(CMRE9_alpha_r8100$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

######-- Fig. 2A --####

Fig2A_data <- CMRE9_alpha_r8100 %>% 
  select(Treatment,Sampletype,Observed,
         Concentration,Lib)

Fig2A_data_CTR <- Fig2A_data %>% filter(Treatment=="Ctr")

Fig2A_data_Obs_means <- Fig2A_data_CTR %>% 
  group_by(Sampletype) %>% 
  summarise(MeanObserved=mean(Observed),
            QuantDobserved=quantile(Observed,probs = 0.95))


Fig2A_data_anova <- aov(Observed ~ Treatment * Sampletype + Lib, data = Fig2A_data)

Fig2A_tukey <- TukeyHSD(Fig2A_data_anova)# Tukey's test
#The use of letters to indicate significant differences in pairwise comparisons is called compact letter display, and can simplify the visualisation and discussion of significant differences among means. We are going to use the multcompLetters4 function from the multcompView package. The arguments are the object from an aov function and the object from the TukeyHSD function.
# compact letter display
Fig2A_cld <-as.data.frame.list(multcompView::multcompLetters4(Fig2A_data_anova, Fig2A_tukey)[[4]])

Fig2A_cld_Ctr <- Fig2A_cld %>% 
  filter(grepl("Ctr", rownames(Fig2A_cld))) #select only the controls
rownames(Fig2A_cld_Ctr) <- gsub("Ctr:","",rownames(Fig2A_cld_Ctr) )

Fig2A_cld_Tbl <- data.frame(Sampletype =rownames(Fig2A_cld_Ctr),
                            L=Fig2A_cld_Ctr$Letters)

Fig2A_data_Obs_means_letters <- merge(Fig2A_data_Obs_means,Fig2A_cld_Tbl)

plot_Fig2A <- ggbarplot(Fig2A_data %>% filter(Treatment=="Ctr"), 
                        x = "Sampletype", 
                        y = "Observed", 
                        fill = "grey50",
                        width=0.2,
                        add = c("mean_se", "jitter"),
                        legend = "none",
                        ylim=c(0,1200)
)      +
  theme_bw() +
  geom_text(data = Fig2A_data_Obs_means_letters, 
            aes(x = Sampletype, y = QuantDobserved+50, label = L), 
            size = 5, vjust=-0.8, hjust =0.5)    +
  ylab("Observed ASV")+  xlab("")  + theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x = element_text(size=14, face = "bold"), axis.text.y = element_text(size=12, face = "bold"))  + theme(plot.margin = unit(c(1,1,-1.5,0), "lines"))

plot_Fig2A

###########----Fig. 2B --############################################
#removing the mean for each group
Fig2B_data_noC <- CMRE9_alpha_r8100 %>% filter(Treatment!="Ctr") %>% select(Sampletype,Observed,Treatment)
#removing the means
Fig2B_data_noC_diff <- merge(Fig2B_data_noC,Fig2A_data_Obs_means[,1:2])
Fig2B_data_noC_diff$Diff <- Fig2B_data_noC_diff$Observed-Fig2B_data_noC_diff$MeanObserved

#with letters
Fig2B_cld <- Fig2A_cld %>% filter(!grepl("Ctr", rownames(Fig2A_cld))) #select all, but the controls

#Compact letter display to indicate significant differences
Fig2B_cld_Tbl <- data.frame(Sampletype =rownames(Fig2B_cld),L=Fig2B_cld$Letters)
Fig2B_cld_Tbl1 <- separate(Fig2B_cld_Tbl,col="Sampletype",sep = ":",into=c("Treatment","Sampletype"))
Fig2B_data_Obs_means <- aggregate(Observed ~  Sampletype*Treatment, CMRE9_alpha_r8100, mean)
Fig2B_data_Obs_sd <- aggregate(Observed ~  Sampletype*Treatment, CMRE9_alpha_r8100, sd)
colnames(Fig2B_data_Obs_sd)[3] <- "sd"
#merge
Fig2B_data_mergeX.SD <- merge(Fig2B_data_Obs_means,Fig2B_data_Obs_sd,
                              by=c("Sampletype","Treatment"))
Fig2B_data_mergeX.SD.L <- merge(Fig2B_data_mergeX.SD,
                                Fig2B_cld_Tbl1,
                                by=c("Sampletype","Treatment"))

#merging with control group mean 
Fig2B_data_mergeX.SD.L.C <- merge(Fig2B_data_mergeX.SD.L,Fig2A_data_Obs_means[,1:2],by="Sampletype")

rename(Fig2B_data_mergeX.SD.L.C,MeanControl= MeanObserved )
Fig2B_data_mergeX.SD.L.C <- Fig2B_data_mergeX.SD.L.C %>% mutate(Diff=Observed-MeanObserved)
#removing the control
Fig2B_data_mergeX.SD.L.C_noC <- Fig2B_data_mergeX.SD.L.C %>% filter(Treatment!="Ctr") 

### rect  
Fig2B_data_mergeX.SD.L.C_noC_Sort <- Fig2B_data_mergeX.SD.L.C_noC %>% arrange(Sampletype,Treatment) %>% 
  mutate(MeanDiff=MeanObserved+Diff)

Fig2B_data_mergeX.SD.L.C_noC_Sort$N <- 1:nrow(Fig2B_data_mergeX.SD.L.C_noC_Sort)
Fig2B_data_mergeX.SD.L.C_noC_Sort

# AS"#E69F00", 
# BX"#009E73",  
# Tb"#D55E00" 
plot_Fig2B_prep <- ggplot(Fig2B_data_mergeX.SD.L.C_noC_Sort,
                          aes(xmin=N, 
                              xmax=N+0.8, #N+0.5
                              ymin=MeanObserved, 
                              ymax=MeanDiff, 
                              fill=Treatment )) + geom_rect() +
  scale_fill_manual(values=c("#E69F00", "#009E73", "#D55E00"))+
  labs(x= "", y = "Difference in Observed ASV \n to ctr mean") + theme_set(theme_bw()) +
  theme(  panel.grid.minor = element_blank(),
          # axis.text.x = element_blank(),
          # axis.ticks = element_blank(),
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=12), 
          axis.text.x  = element_text(colour="black", vjust=0.5, size=14), 
          axis.text.y  = element_text(colour="black", vjust=0.5, size=14))+
  ylim(0, 1200)

# plot_Fig2B_prep

NewY=Fig2B_data_mergeX.SD.L.C_noC_Sort$MeanDiff+ #adjusting the high of the labels
  c(0,0,0,
    0,0,0,
    50,60,40,
    0,0,0,
    -20,-20,-20)
plot_Fig2B <- plot_Fig2B_prep +
  geom_text(aes(x =N+0.25,y=NewY,label=L),
            position = position_dodge(width = 1), 
            angle = 0,size = 5, vjust=-2) +
  scale_x_continuous(breaks  = c(2+3*(0:4)), #
                     labels = c("Water","Sediment","Soil","Root","Mouse")
  ) +
  geom_segment(aes(x=1,xend=3+0.8, 
                   y=MeanObserved[1],
                   yend=MeanObserved[3])) +
  geom_segment(aes(x=4,xend=6+0.8,
                   y=MeanObserved[4],
                   yend=MeanObserved[4])) +
  geom_segment(aes(x=7,xend=9+0.8,
                   y=MeanObserved[7],
                   yend=MeanObserved[7])) +
  geom_segment(aes(x=10,xend=12+0.8,
                   y=MeanObserved[10],
                   yend=MeanObserved[10])) +
  geom_segment(aes(x=13,xend=15+0.8,
                   y=MeanObserved[13],
                   yend=MeanObserved[13])) +
  theme_bw() + theme(panel.grid.minor = element_blank())  + theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=14), axis.text.x = element_text(size=14, face = "bold"), axis.text.y = element_text(size=12, face = "bold"))  + theme(legend.position = "none") + theme(plot.margin = unit(c(1,1,-1.5,0), "lines"))


plot_Fig2B

# Combination of plots
require(patchwork)
C <- plot_Fig2A + plot_Fig2B
C
ggplot2::ggsave("Fig.2_AB_CMRE9_r8100_alphaDiv_observed_Tukey_Boxplot_centered_letters_bars_withoutlegand.tiff", 
                width = 10.30, height = 5.30, dpi=300)

############################---Supplementary Figure- Shannon final box plot --------##############################
#### for treatement and compartment 

CMRE5_alpha_r8100$Sampletype<- factor(CMRE5_alpha_r8100$Sampletype, levels = c("Water", "Sediment", "Soil", "Root", "Mouse"))
CMRE5_alpha_r8100$Treatment<- factor(CMRE5_alpha_r8100$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

######-- Fig. 2A --####

Fig2A_data <- CMRE5_alpha_r8100 %>% 
  select(Treatment,Sampletype,Shannon,
         Concentration,Lib)

Fig2A_data_CTR <- Fig2A_data %>% filter(Treatment=="Ctr")

Fig2A_data_Obs_means <- Fig2A_data_CTR %>% 
  group_by(Sampletype) %>% 
  summarise(MeanShannon=mean(Shannon),
            QuantDShannon=quantile(Shannon,probs = 0.95))


Fig2A_data_anova <- aov(Shannon ~ Treatment * Sampletype + Lib, data = Fig2A_data)

Fig2A_tukey <- TukeyHSD(Fig2A_data_anova)# Tukey's test
#The use of letters to indicate significant differences in pairwise comparisons is called compact letter display, and can simplify the visualisation and discussion of significant differences among means. We are going to use the multcompLetters4 function from the multcompView package. The arguments are the object from an aov function and the object from the TukeyHSD function.
# compact letter display
Fig2A_cld <-as.data.frame.list(multcompView::multcompLetters4(Fig2A_data_anova, Fig2A_tukey)[[4]])

Fig2A_cld_Ctr <- Fig2A_cld %>% 
  filter(grepl("Ctr", rownames(Fig2A_cld))) #select only the controls
rownames(Fig2A_cld_Ctr) <- gsub("Ctr:","",rownames(Fig2A_cld_Ctr) )

Fig2A_cld_Tbl <- data.frame(Sampletype =rownames(Fig2A_cld_Ctr),
                            L=Fig2A_cld_Ctr$Letters)

Fig2A_data_Obs_means_letters <- merge(Fig2A_data_Obs_means,Fig2A_cld_Tbl)

plot_Fig2A <- ggbarplot(Fig2A_data %>% filter(Treatment=="Ctr"), 
                        x = "Sampletype", 
                        y = "Shannon", 
                        fill = "grey50",
                        width=0.2,
                        add = c("mean_se", "jitter"),
                        legend = "none",
                        ylim=c(0,8)
)      +
  theme_bw() +
  geom_text(data = Fig2A_data_Obs_means_letters, 
            aes(x = Sampletype, y = QuantDShannon, label = L), 
            size = 5,  vjust=-0.8, hjust =0.5)    +
  ylab("Shannon")+  xlab("") + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=14, face = "bold"), axis.text.y = element_text(size=12, face = "bold"))  + theme(plot.margin = unit(c(1,1,-1.5,0), "lines"))

plot_Fig2A

###########----Fig. 2B --############################################
#removing the mean for each group
Fig2B_data_noC <- CMRE5_alpha_r8100 %>% filter(Treatment!="Ctr") %>% select(Sampletype,Shannon,Treatment)
#removing the means
Fig2B_data_noC_diff <- merge(Fig2B_data_noC,Fig2A_data_Obs_means[,1:2])
Fig2B_data_noC_diff$Diff <- Fig2B_data_noC_diff$Shannon-Fig2B_data_noC_diff$MeanShannon

#with letters
Fig2B_cld <- Fig2A_cld %>% filter(!grepl("Ctr", rownames(Fig2A_cld))) #select all, but the controls

#Compact letter display to indicate significant differences
Fig2B_cld_Tbl <- data.frame(Sampletype =rownames(Fig2B_cld),L=Fig2B_cld$Letters)
Fig2B_cld_Tbl1 <- separate(Fig2B_cld_Tbl,col="Sampletype",sep = ":",into=c("Treatment","Sampletype"))
Fig2B_data_Obs_means <- aggregate(Shannon ~  Sampletype*Treatment, CMRE5_alpha_r8100, mean)
Fig2B_data_Obs_sd <- aggregate(Shannon ~  Sampletype*Treatment, CMRE5_alpha_r8100, sd)
colnames(Fig2B_data_Obs_sd)[3] <- "sd"
#merge
Fig2B_data_mergeX.SD <- merge(Fig2B_data_Obs_means,Fig2B_data_Obs_sd,
                              by=c("Sampletype","Treatment"))
Fig2B_data_mergeX.SD.L <- merge(Fig2B_data_mergeX.SD,
                                Fig2B_cld_Tbl1,
                                by=c("Sampletype","Treatment"))

#merging with control group mean 
Fig2B_data_mergeX.SD.L.C <- merge(Fig2B_data_mergeX.SD.L,Fig2A_data_Obs_means[,1:2],by="Sampletype")

rename(Fig2B_data_mergeX.SD.L.C,MeanControl= MeanShannon )
Fig2B_data_mergeX.SD.L.C <- Fig2B_data_mergeX.SD.L.C %>% mutate(Diff=Shannon-MeanShannon)
#removing the control
Fig2B_data_mergeX.SD.L.C_noC <- Fig2B_data_mergeX.SD.L.C %>% filter(Treatment!="Ctr") 

### rect  
Fig2B_data_mergeX.SD.L.C_noC_Sort <- Fig2B_data_mergeX.SD.L.C_noC %>% arrange(Sampletype,Treatment) %>% 
  mutate(MeanDiff=MeanShannon+Diff)

Fig2B_data_mergeX.SD.L.C_noC_Sort$N <- 1:nrow(Fig2B_data_mergeX.SD.L.C_noC_Sort)
Fig2B_data_mergeX.SD.L.C_noC_Sort

# AS"#E69F00", 
# BX"#009E73",  
# Tb"#D55E00" 
plot_Fig2B_prep <- ggplot(Fig2B_data_mergeX.SD.L.C_noC_Sort,
                          aes(xmin=N, 
                              xmax=N+0.8, #N+0.5
                              ymin=MeanShannon, 
                              ymax=MeanDiff, 
                              fill=Treatment )) + geom_rect() +
  scale_fill_manual(values=c("#E69F00", "#009E73", "#D55E00"))+
  labs(x= "", y = "Difference in Shannon to ctr mean") + theme_set(theme_bw()) +
  theme(  panel.grid.minor = element_blank(),
          # axis.text.x = element_blank(),
          # axis.ticks = element_blank(),
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=12), 
          axis.text.x  = element_text(colour="black", vjust=0.5, size=14), 
          axis.text.y  = element_text(colour="black", vjust=0.5, size=14))+
  ylim(0, 8)

# plot_Fig2B_prep

NewY=Fig2B_data_mergeX.SD.L.C_noC_Sort$MeanDiff+ #adjusting the high of the labels
  c(0,0,-0.05,
    0,0,0,
    0.1,0.10,0.1,
    -0.1,-0.1,-0.1,
    -0.1,-0.1,-0.1
  )
plot_Fig2B <- plot_Fig2B_prep +
  geom_text(aes(x =N+0.25,y=NewY,label=L),
            position = position_dodge(width = 1), 
            angle = 0,size = 5, vjust=-2) +
  scale_x_continuous(breaks  = c(2+3*(0:4)), #
                     labels = c("Water","Sediment","Soil","Root","Mouse")
  ) +
  geom_segment(aes(x=1,xend=3+0.8, 
                   y=MeanShannon[1],
                   yend=MeanShannon[3])) +
  geom_segment(aes(x=4,xend=6+0.8,
                   y=MeanShannon[4],
                   yend=MeanShannon[4])) +
  geom_segment(aes(x=7,xend=9+0.8,
                   y=MeanShannon[7],
                   yend=MeanShannon[7])) +
  geom_segment(aes(x=10,xend=12+0.8,
                   y=MeanShannon[10],
                   yend=MeanShannon[10])) +
  geom_segment(aes(x=13,xend=15+0.8,
                   y=MeanShannon[13],
                   yend=MeanShannon[13])) +
  theme_bw() +theme(panel.grid.minor = element_blank()) + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=14, face = "bold"), axis.text.x = element_text(size=14, face = "bold"), axis.text.y = element_text(size=12, face = "bold")) + theme(plot.margin = unit(c(1,1,-1.5,0), "lines"))


plot_Fig2B

# Combination of plots
require(patchwork)
plot_Fig2A + plot_Fig2B

ggplot2::ggsave("Fig.2_Suppl_Shannon.tiff", 
                width = 10.80, height = 6.30, dpi=300)


######################################################################################################
###########----Fig. 2C --############################################
###########----CAP plot Water
CMRE5_r8100_plot_bray_cap_Water = plot_ordination(CMRE5_r8100_Water1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Water$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Water = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")


###########----CAP plot Sediment
CMRE5_r8100_plot_bray_cap_Sediment = plot_ordination(CMRE5_r8100_Sediment1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Sediment$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Sediment = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")


###########----CAP plot Soil
CMRE5_r8100_plot_bray_cap_Soil = plot_ordination(CMRE5_r8100_Soil1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Soil$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Soil = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")


###########----CAP plot Root
CMRE5_r8100_plot_bray_cap_Root = plot_ordination(CMRE5_r8100_Root1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Root$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Root = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")


###########----CAP plot Mouse
CMRE5_r8100_plot_bray_cap_mouse = plot_ordination(CMRE5_r8100_mouse1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_mouse$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_mouse = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")


###############---Combined figure for all beta cap compartments-------##############

A <- ggarrange(CMRE5_r8100_plot_bray_cap_Water, CMRE5_r8100_plot_bray_cap_Sediment, ncol = 3) 
B <- ggarrange(CMRE5_r8100_plot_bray_cap_Soil, CMRE5_r8100_plot_bray_cap_Root, CMRE5_r8100_plot_bray_cap_mouse, ncol = 3)
ggarrange (A, B, ncol = 1)

###############---Final figure 2-------##############

ggarrange(C, A, B, ncol = 1, heights = c(5, 3.6, 3.6), nrow = 3, align = "h") #

ggsave("Fig.2_3.tiff", width = 10.30, height = 9.30, dpi=600)

########################################################################################################################


######################################################################################################
###########################----------------Figure3 scripts--------------------########################

##################--Water--As---######

CMRE5_High_Water_CtrAs1_res1 <- as.data.frame(CMRE5_High_Water_CtrAs1_res)
CMRE5_High_Water_CtrAs1_res1$log2FoldChange1 = ifelse(CMRE5_High_Water_CtrAs1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Water_CtrAs1_res1$padj1 = ifelse(CMRE5_High_Water_CtrAs1_res1$padj > 0.05,"high","low")

CMRE5_High_Water_CtrAs1_res1$color_col <- paste(CMRE5_High_Water_CtrAs1_res1$log2FoldChange1, "-", CMRE5_High_Water_CtrAs1_res1$padj1)

CMRE5_High_Water_CtrAs1_res2 <- CMRE5_High_Water_CtrAs1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Water_CtrAs_MAplot3 <- ggplot(CMRE5_High_Water_CtrAs1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") +  scale_y_continuous(limits=c(-23,20), breaks=c(-20, -15, -10, -5,  0, 5, 10, 15, 20), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(0, 14), breaks=c(0, 2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Water--Bx---######

CMRE5_High_Water_CtrBx1_res1 <- as.data.frame(CMRE5_High_Water_CtrBx1_res)
CMRE5_High_Water_CtrBx1_res1$log2FoldChange1 = ifelse(CMRE5_High_Water_CtrBx1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Water_CtrBx1_res1$padj1 = ifelse(CMRE5_High_Water_CtrBx1_res1$padj > 0.05,"high","low")

CMRE5_High_Water_CtrBx1_res1$color_col <- paste(CMRE5_High_Water_CtrBx1_res1$log2FoldChange1, "-", CMRE5_High_Water_CtrBx1_res1$padj1)

CMRE5_High_Water_CtrBx1_res2 <- CMRE5_High_Water_CtrBx1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Water_CtrBx_MAplot3 <- ggplot(CMRE5_High_Water_CtrBx1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(0, 14), breaks=c(0, 2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Water--Tb---######

CMRE5_High_Water_CtrTb1_res1 <- as.data.frame(CMRE5_High_Water_CtrTb1_res)
CMRE5_High_Water_CtrTb1_res1$log2FoldChange1 = ifelse(CMRE5_High_Water_CtrTb1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Water_CtrTb1_res1$padj1 = ifelse(CMRE5_High_Water_CtrTb1_res1$padj > 0.05,"high","low")

CMRE5_High_Water_CtrTb1_res1$color_col <- paste(CMRE5_High_Water_CtrTb1_res1$log2FoldChange1, "-", CMRE5_High_Water_CtrTb1_res1$padj1)

CMRE5_High_Water_CtrTb1_res2 <- CMRE5_High_Water_CtrTb1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Water_CtrTb_MAplot3 <- ggplot(CMRE5_High_Water_CtrTb1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(0, 14), breaks=c(0, 2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Sediment--As---######

CMRE5_High_Sediment_CtrAs1_res1 <- as.data.frame(CMRE5_High_Sediment_CtrAs1_res)
CMRE5_High_Sediment_CtrAs1_res1$log2FoldChange1 = ifelse(CMRE5_High_Sediment_CtrAs1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Sediment_CtrAs1_res1$padj1 = ifelse(CMRE5_High_Sediment_CtrAs1_res1$padj > 0.05,"high","low")

CMRE5_High_Sediment_CtrAs1_res1$color_col <- paste(CMRE5_High_Sediment_CtrAs1_res1$log2FoldChange1, "-", CMRE5_High_Sediment_CtrAs1_res1$padj1)

CMRE5_High_Sediment_CtrAs1_res2 <- CMRE5_High_Sediment_CtrAs1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Sediment_CtrAs_MAplot3 <- ggplot(CMRE5_High_Sediment_CtrAs1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-20,20), breaks=c(-20, -15, -10, -5,  0, 5, 10, 15, 20 ), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(-2, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Sediment--Bx---######

CMRE5_High_Sediment_CtrBx1_res1 <- as.data.frame(CMRE5_High_Sediment_CtrBx1_res)
CMRE5_High_Sediment_CtrBx1_res1$log2FoldChange1 = ifelse(CMRE5_High_Sediment_CtrBx1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Sediment_CtrBx1_res1$padj1 = ifelse(CMRE5_High_Sediment_CtrBx1_res1$padj > 0.05,"high","low")

CMRE5_High_Sediment_CtrBx1_res1$color_col <- paste(CMRE5_High_Sediment_CtrBx1_res1$log2FoldChange1, "-", CMRE5_High_Sediment_CtrBx1_res1$padj1)
CMRE5_High_Sediment_CtrBx1_res2 <- CMRE5_High_Sediment_CtrBx1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Sediment_CtrBx_MAplot3 <- ggplot(CMRE5_High_Sediment_CtrBx1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(0, 11), breaks=c( 0, 2, 4, 6, 8, 10), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0)  + geom_hline(yintercept=0) 

##################--Sediment--Tb---######

CMRE5_High_Sediment_CtrTb1_res1 <- as.data.frame(CMRE5_High_Sediment_CtrTb1_res)
CMRE5_High_Sediment_CtrTb1_res1$log2FoldChange1 = ifelse(CMRE5_High_Sediment_CtrTb1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Sediment_CtrTb1_res1$padj1 = ifelse(CMRE5_High_Sediment_CtrTb1_res1$padj > 0.05,"high","low")

CMRE5_High_Sediment_CtrTb1_res1$color_col <- paste(CMRE5_High_Sediment_CtrTb1_res1$log2FoldChange1, "-", CMRE5_High_Sediment_CtrTb1_res1$padj1)
CMRE5_High_Sediment_CtrTb1_res2 <- CMRE5_High_Sediment_CtrTb1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Sediment_CtrTb_MAplot3 <- ggplot(CMRE5_High_Sediment_CtrTb1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(-2, 11), breaks=c(-2,  0, 2, 4, 6, 8, 10), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Soil--As---######

CMRE5_High_Soil_CtrAs1_res1 <- as.data.frame(CMRE5_High_Soil_CtrAs1_res)
CMRE5_High_Soil_CtrAs1_res1$log2FoldChange1 = ifelse(CMRE5_High_Soil_CtrAs1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Soil_CtrAs1_res1$padj1 = ifelse(CMRE5_High_Soil_CtrAs1_res1$padj > 0.05,"high","low")
CMRE5_High_Soil_CtrAs1_res1$color_col <- paste(CMRE5_High_Soil_CtrAs1_res1$log2FoldChange1, "-", CMRE5_High_Soil_CtrAs1_res1$padj1)

CMRE5_High_Soil_CtrAs1_res2 <- CMRE5_High_Soil_CtrAs1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Soil_CtrAs_MAplot3 <- ggplot(CMRE5_High_Soil_CtrAs1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(2, 14), breaks=c(2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Soil--Bx---######

CMRE5_High_Soil_CtrBx1_res1 <- as.data.frame(CMRE5_High_Soil_CtrBx1_res)
CMRE5_High_Soil_CtrBx1_res1$log2FoldChange1 = ifelse(CMRE5_High_Soil_CtrBx1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Soil_CtrBx1_res1$padj1 = ifelse(CMRE5_High_Soil_CtrBx1_res1$padj > 0.05,"high","low")
CMRE5_High_Soil_CtrBx1_res1$color_col <- paste(CMRE5_High_Soil_CtrBx1_res1$log2FoldChange1, "-", CMRE5_High_Soil_CtrBx1_res1$padj1)

CMRE5_High_Soil_CtrBx1_res2 <- CMRE5_High_Soil_CtrBx1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Soil_CtrBx_MAplot3 <- ggplot(CMRE5_High_Soil_CtrBx1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(2, 14), breaks=c(2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Soil--Tb---######

CMRE5_High_Soil_CtrTb1_res1 <- as.data.frame(CMRE5_High_Soil_CtrTb1_res)
CMRE5_High_Soil_CtrTb1_res1$log2FoldChange1 = ifelse(CMRE5_High_Soil_CtrTb1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Soil_CtrTb1_res1$padj1 = ifelse(CMRE5_High_Soil_CtrTb1_res1$padj > 0.05,"high","low")

CMRE5_High_Soil_CtrTb1_res1$color_col <- paste(CMRE5_High_Soil_CtrTb1_res1$log2FoldChange1, "-", CMRE5_High_Soil_CtrTb1_res1$padj1)
CMRE5_High_Soil_CtrTb1_res2 <- CMRE5_High_Soil_CtrTb1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Soil_CtrTb_MAplot3 <- ggplot(CMRE5_High_Soil_CtrTb1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(2, 14), breaks=c(2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Root--As---######

CMRE5_High_Root_CtrAs1_res1 <- as.data.frame(CMRE5_High_Root_CtrAs1_res)
CMRE5_High_Root_CtrAs1_res1$log2FoldChange1 = ifelse(CMRE5_High_Root_CtrAs1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Root_CtrAs1_res1$padj1 = ifelse(CMRE5_High_Root_CtrAs1_res1$padj > 0.05,"high","low")

CMRE5_High_Root_CtrAs1_res1$color_col <- paste(CMRE5_High_Root_CtrAs1_res1$log2FoldChange1, "-", CMRE5_High_Root_CtrAs1_res1$padj1)
CMRE5_High_Root_CtrAs1_res2 <- CMRE5_High_Root_CtrAs1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Root_CtrAs_MAplot3 <- ggplot(CMRE5_High_Root_CtrAs1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(-2, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Root--Bx---######

CMRE5_High_Root_CtrBx1_res1 <- as.data.frame(CMRE5_High_Root_CtrBx1_res)
CMRE5_High_Root_CtrBx1_res1$log2FoldChange1 = ifelse(CMRE5_High_Root_CtrBx1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Root_CtrBx1_res1$padj1 = ifelse(CMRE5_High_Root_CtrBx1_res1$padj > 0.05,"high","low")

CMRE5_High_Root_CtrBx1_res1$color_col <- paste(CMRE5_High_Root_CtrBx1_res1$log2FoldChange1, "-", CMRE5_High_Root_CtrBx1_res1$padj1)

CMRE5_High_Root_CtrBx1_res2 <- CMRE5_High_Root_CtrBx1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Root_CtrBx_MAplot3 <- ggplot(CMRE5_High_Root_CtrBx1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(-2, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Root--Tb---######

CMRE5_High_Root_CtrTb1_res1 <- as.data.frame(CMRE5_High_Root_CtrTb1_res)
CMRE5_High_Root_CtrTb1_res1$log2FoldChange1 = ifelse(CMRE5_High_Root_CtrTb1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Root_CtrTb1_res1$padj1 = ifelse(CMRE5_High_Root_CtrTb1_res1$padj > 0.05,"high","low")

CMRE5_High_Root_CtrTb1_res1$color_col <- paste(CMRE5_High_Root_CtrTb1_res1$log2FoldChange1, "-", CMRE5_High_Root_CtrTb1_res1$padj1)

CMRE5_High_Root_CtrTb1_res2 <- CMRE5_High_Root_CtrTb1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Root_CtrTb_MAplot3 <- ggplot(CMRE5_High_Root_CtrTb1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") + scale_y_continuous(limits=c(-6,6), breaks=c(-6, -4, -2,  0, 2,  4, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(-2, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 


##################--Mouse--As---######

CMRE5_High_Mouse_CtrAs1_res1 <- as.data.frame(CMRE5_High_Mouse_CtrAs1_res)
CMRE5_High_Mouse_CtrAs1_res1$log2FoldChange1 = ifelse(CMRE5_High_Mouse_CtrAs1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Mouse_CtrAs1_res1$padj1 = ifelse(CMRE5_High_Mouse_CtrAs1_res1$padj > 0.05,"high","low")
CMRE5_High_Mouse_CtrAs1_res1$color_col <- paste(CMRE5_High_Mouse_CtrAs1_res1$log2FoldChange1, "-", CMRE5_High_Mouse_CtrAs1_res1$padj1)

CMRE5_High_Mouse_CtrAs1_res2 <- CMRE5_High_Mouse_CtrAs1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Mouse_CtrAs_MAplot3 <- ggplot(CMRE5_High_Mouse_CtrAs1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5,  hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") +  scale_y_continuous(limits=c(-12, 12.5), breaks=c( -10, -5, 0, 5, 10), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(-2, 14), breaks=c(-2,  0, 2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1))+ geom_hline(yintercept=0) 

##################--Mouse--Bx---######

CMRE5_High_Mouse_CtrBx1_res1 <- as.data.frame(CMRE5_High_Mouse_CtrBx1_res)
CMRE5_High_Mouse_CtrBx1_res1$log2FoldChange1 = ifelse(CMRE5_High_Mouse_CtrBx1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Mouse_CtrBx1_res1$padj1 = ifelse(CMRE5_High_Mouse_CtrBx1_res1$padj > 0.05,"high","low")
CMRE5_High_Mouse_CtrBx1_res1$color_col <- paste(CMRE5_High_Mouse_CtrBx1_res1$log2FoldChange1, "-", CMRE5_High_Mouse_CtrBx1_res1$padj1)

CMRE5_High_Mouse_CtrBx1_res2 <- CMRE5_High_Mouse_CtrBx1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Mouse_CtrBx_MAplot3<- ggplot(CMRE5_High_Mouse_CtrBx1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none") +   scale_y_continuous(limits=c(-20,25), breaks=c(-20, -10,  0,  10, 20), labels = scales::number_format(accuracy = 0.1))  +   scale_x_continuous(limits=c(-2, 14), breaks=c(-2,  0, 2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##################--Mouse--Tb---######

CMRE5_High_Mouse_CtrTb1_res1 <- as.data.frame(CMRE5_High_Mouse_CtrTb1_res)
CMRE5_High_Mouse_CtrTb1_res1$log2FoldChange1 = ifelse(CMRE5_High_Mouse_CtrTb1_res1$log2FoldChange > 0.00,"pos","neg")
CMRE5_High_Mouse_CtrTb1_res1$padj1 = ifelse(CMRE5_High_Mouse_CtrTb1_res1$padj > 0.05,"high","low")

CMRE5_High_Mouse_CtrTb1_res1$color_col <- paste(CMRE5_High_Mouse_CtrTb1_res1$log2FoldChange1, "-", CMRE5_High_Mouse_CtrTb1_res1$padj1)
CMRE5_High_Mouse_CtrTb1_res2 <- CMRE5_High_Mouse_CtrTb1_res1 %>% mutate(color_col = recode(color_col, 'neg - low' = "C", 'pos - low'= "D", 'neg - high'= "A", 'pos - high'= "B"))

New1_ASV_relative_abundance_Mouse_CtrTb_MAplot3 <- ggplot(CMRE5_High_Mouse_CtrTb1_res2 %>% arrange(color_col), aes(x = log2(baseMean), y = log2FoldChange, color = color_col, alpha=0.7))  + geom_jitter(size=3.0,  height = 0.2) + scale_color_manual(values = c( "darkgray",  "darkgray" , "#B31B21", "#1465AC" ))  + theme(axis.title.y=element_blank()) + theme_bw() + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "Mean Relative Abundance (log2)", y = "Fold Change (log2)") + theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, hjust = 0.5, size=11), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=11)) + labs(color='log2FC') + theme(legend.position = "none")  +   scale_y_continuous(limits=c(-6, 7.5), breaks=c(-6, -3,  0,  3, 6), labels = scales::number_format(accuracy = 0.1)) +   scale_x_continuous(limits=c(-2, 14), breaks=c(-2,  0, 2, 4, 6, 8, 10, 12, 14), labels = scales::number_format(accuracy = 0.1)) + geom_hline(yintercept=0) 

##########################################################################################################################
#####################----Combined plot for all above----##################

figure3 <- ggarrange(New1_ASV_relative_abundance_Water_CtrAs_MAplot3, 
                    New1_ASV_relative_abundance_Sediment_CtrAs_MAplot3, 
                    New1_ASV_relative_abundance_Soil_CtrAs_MAplot3, 
                    New1_ASV_relative_abundance_Root_CtrAs_MAplot3, 
                    New1_ASV_relative_abundance_Mouse_CtrAs_MAplot3, 
                    New1_ASV_relative_abundance_Water_CtrBx_MAplot3, 
                    New1_ASV_relative_abundance_Sediment_CtrBx_MAplot3, 
                    New1_ASV_relative_abundance_Soil_CtrBx_MAplot3, 
                    New1_ASV_relative_abundance_Root_CtrBx_MAplot3, 
                    New1_ASV_relative_abundance_Mouse_CtrBx_MAplot3,
                    New1_ASV_relative_abundance_Water_CtrTb_MAplot3, 
                    New1_ASV_relative_abundance_Sediment_CtrTb_MAplot3, 
                    New1_ASV_relative_abundance_Soil_CtrTb_MAplot3, 
                    New1_ASV_relative_abundance_Root_CtrTb_MAplot3, 
                    New1_ASV_relative_abundance_Mouse_CtrTb_MAplot3,
                    ncol = 5, nrow = 3)

########################################################################################################


######################################################################################################
###########################----------------Figure4 scripts--------------------########################

##########-Mouse Ctr

CMRE5_Mouseb_Ctrc_otu_net1 <- symBeta(getOptBeta(CMRE5_Mouseb_Ctrc_otu_net))
colnames(CMRE5_Mouseb_Ctrc_otu_net1) <- rownames(CMRE5_Mouseb_Ctrc_otu_net1) <- colnames(CMRE5_Mouseb_Ctrc_otu)
Mouse.ig <- adj2igraph(getRefit(CMRE5_Mouseb_Ctrc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Mouseb_Ctrc)))
vsize <- rowMeans(clr(CMRE5_Mouseb_Ctrc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Mouse.ig)

otu.ids=colnames(CMRE5_Mouseb_Ctrc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Ctrc_otu_net)))

edges=E(Mouse.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Mouse.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Mouse.ig)$color=edge.colors
Mouse.net <- asNetwork(Mouse.Ctr.Mouse.ig)

net.hs <- hub_score(Mouse.Ctr.Mouse.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Mouse.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Mouse.net %v% "hubscore", levels = c("low","high"))) 

M_Ctr <- ggnet2(Mouse.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))

##########-Mouse As  

CMRE5_Mouseb_Asc_otu_net1 <- symBeta(getOptBeta(CMRE5_Mouseb_Asc_otu_net))
colnames(CMRE5_Mouseb_Asc_otu_net1) <- rownames(CMRE5_Mouseb_Asc_otu_net1) <- colnames(CMRE5_Mouseb_Asc_otu)

Mouse.ig <- adj2igraph(getRefit(CMRE5_Mouseb_Asc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Mouseb_Asc)))
vsize <- rowMeans(clr(CMRE5_Mouseb_Asc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Mouse.ig)

otu.ids=colnames(CMRE5_Mouseb_Asc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Asc_otu_net)))

edges=E(Mouse.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Mouse.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Mouse.ig)$color=edge.colors
Mouse.net <- asNetwork(Mouse.As.Mouse.ig)

net.hs <- hub_score(Mouse.As.Mouse.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Mouse.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Mouse.net %v% "hubscore", levels = c("low","high"))) 

M_As <- ggnet2(Mouse.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))

##########-Mouse Bx

CMRE5_Mouseb_Bxc_otu_net1 <- symBeta(getOptBeta(CMRE5_Mouseb_Bxc_otu_net))
colnames(CMRE5_Mouseb_Bxc_otu_net1) <- rownames(CMRE5_Mouseb_Bxc_otu_net1) <- colnames(CMRE5_Mouseb_Bxc_otu)  
Mouse.ig <- adj2igraph(getRefit(CMRE5_Mouseb_Bxc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Mouseb_Bxc)))
vsize <- rowMeans(clr(CMRE5_Mouseb_Bxc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Mouse.ig)

otu.ids=colnames(CMRE5_Mouseb_Bxc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Bxc_otu_net)))

edges=E(Mouse.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Mouse.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Mouse.ig)$color=edge.colors
Mouse.net <- asNetwork(Mouse.Bx.Mouse.ig)

net.hs <- hub_score(Mouse.Bx.Mouse.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Mouse.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Mouse.net %v% "hubscore", levels = c("low","high"))) 

M_Bx <- ggnet2(Mouse.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Mouse Tb

CMRE5_Mouseb_Tbc_otu_net1 <- symBeta(getOptBeta(CMRE5_Mouseb_Tbc_otu_net))
colnames(CMRE5_Mouseb_Tbc_otu_net1) <- rownames(CMRE5_Mouseb_Tbc_otu_net1) <- colnames(CMRE5_Mouseb_Tbc_otu)
Mouse.ig <- adj2igraph(getRefit(CMRE5_Mouseb_Tbc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Mouseb_Tbc)))
vsize <- rowMeans(clr(CMRE5_Mouseb_Tbc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Mouse.ig)

otu.ids=colnames(CMRE5_Mouseb_Tbc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Tbc_otu_net)))

edges=E(Mouse.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Mouse.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Mouse.ig)$color=edge.colors
Mouse.net <- asNetwork(Mouse.Tb.Mouse.ig)

net.hs <- hub_score(Mouse.Tb.Mouse.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Mouse.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Mouse.net %v% "hubscore", levels = c("low","high"))) 

M_Tb <- ggnet2(Mouse.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Soil Ctr

CMRE5_Soilb_Ctrc_otu_net1 <- symBeta(getOptBeta(CMRE5_Soilb_Ctrc_otu_net))
colnames(CMRE5_Soilb_Ctrc_otu_net1) <- rownames(CMRE5_Soilb_Ctrc_otu_net1) <- colnames(CMRE5_Soilb_Ctrc_otu)
Soil.ig <- adj2igraph(getRefit(CMRE5_Soilb_Ctrc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Soilb_Ctrc)))
vsize <- rowMeans(clr(CMRE5_Soilb_Ctrc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Soil.ig)

otu.ids=colnames(CMRE5_Soilb_Ctrc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Ctrc_otu_net)))

edges=E(Soil.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Soil.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Soil.ig)$color=edge.colors
Soil.net <- asNetwork(Soil.Ctr.Soil.ig)

net.hs <- hub_score(Soil.Ctr.Soil.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Soil.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Soil.net %v% "hubscore", levels = c("low","high"))) 

S_Ctr <- ggnet2(Soil.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Soil As

CMRE5_Soilb_Asc_otu_net1 <- symBeta(getOptBeta(CMRE5_Soilb_Asc_otu_net))
colnames(CMRE5_Soilb_Asc_otu_net1) <- rownames(CMRE5_Soilb_Asc_otu_net1) <- colnames(CMRE5_Soilb_Asc_otu)
Soil.ig <- adj2igraph(getRefit(CMRE5_Soilb_Asc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Soilb_Asc)))
vsize <- rowMeans(clr(CMRE5_Soilb_Asc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Soil.ig)

otu.ids=colnames(CMRE5_Soilb_Asc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Asc_otu_net)))

edges=E(Soil.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Soil.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Soil.ig)$color=edge.colors
Soil.net <- asNetwork(Soil.As.Soil.ig)

net.hs <- hub_score(Soil.As.Soil.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Soil.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Soil.net %v% "hubscore", levels = c("low","high"))) 

S_As <- ggnet2(Soil.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Soil Bx

CMRE5_Soilb_Bxc_otu_net1 <- symBeta(getOptBeta(CMRE5_Soilb_Bxc_otu_net))
colnames(CMRE5_Soilb_Bxc_otu_net1) <- rownames(CMRE5_Soilb_Bxc_otu_net1) <- colnames(CMRE5_Soilb_Bxc_otu)
Soil.ig <- adj2igraph(getRefit(CMRE5_Soilb_Bxc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Soilb_Bxc)))
vsize <- rowMeans(clr(CMRE5_Soilb_Bxc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Soil.ig)

otu.ids=colnames(CMRE5_Soilb_Bxc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Bxc_otu_net)))

edges=E(Soil.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Soil.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Soil.ig)$color=edge.colors
Soil.net <- asNetwork(Soil.Bx.Soil.ig)

net.hs <- hub_score(Soil.Bx.Soil.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Soil.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Soil.net %v% "hubscore", levels = c("low","high"))) 

S_Bx <- ggnet2(Soil.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Soil Tb

CMRE5_Soilb_Tbc_otu_net1 <- symBeta(getOptBeta(CMRE5_Soilb_Tbc_otu_net))
colnames(CMRE5_Soilb_Tbc_otu_net1) <- rownames(CMRE5_Soilb_Tbc_otu_net1) <- colnames(CMRE5_Soilb_Tbc_otu)
Soil.ig <- adj2igraph(getRefit(CMRE5_Soilb_Tbc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Soilb_Tbc)))
vsize <- rowMeans(clr(CMRE5_Soilb_Tbc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Soil.ig)

otu.ids=colnames(CMRE5_Soilb_Tbc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Tbc_otu_net)))

edges=E(Soil.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Soil.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Soil.ig)$color=edge.colors
Soil.net <- asNetwork(Soil.Tb.Soil.ig)

net.hs <- hub_score(Soil.Tb.Soil.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Soil.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Soil.net %v% "hubscore", levels = c("low","high"))) 

S_Tb <- ggnet2(Soil.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Root Ctr

CMRE5_Rootb_Ctrc_otu_net1 <- symBeta(getOptBeta(CMRE5_Rootb_Ctrc_otu_net))
colnames(CMRE5_Rootb_Ctrc_otu_net1) <- rownames(CMRE5_Rootb_Ctrc_otu_net1) <- colnames(CMRE5_Rootb_Ctrc_otu)
Root.ig <- adj2igraph(getRefit(CMRE5_Rootb_Ctrc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Rootb_Ctrc)))
vsize <- rowMeans(clr(CMRE5_Rootb_Ctrc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Root.ig)

otu.ids=colnames(CMRE5_Rootb_Ctrc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Ctrc_otu_net)))

edges=E(Root.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Root.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Root.ig)$color=edge.colors
Root.net <- asNetwork(Root.Ctr.Root.ig)

net.hs <- hub_score(Root.Ctr.Root.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Root.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Root.net %v% "hubscore", levels = c("low","high"))) 

R_Ctr <- ggnet2(Root.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Root As

CMRE5_Rootb_Asc_otu_net1 <- symBeta(getOptBeta(CMRE5_Rootb_Asc_otu_net))
colnames(CMRE5_Rootb_Asc_otu_net1) <- rownames(CMRE5_Rootb_Asc_otu_net1) <- colnames(CMRE5_Rootb_Asc_otu)
Root.ig <- adj2igraph(getRefit(CMRE5_Rootb_Asc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Rootb_Asc)))
vsize <- rowMeans(clr(CMRE5_Rootb_Asc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Root.ig)

otu.ids=colnames(CMRE5_Rootb_Asc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Asc_otu_net)))

edges=E(Root.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Root.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Root.ig)$color=edge.colors
Root.net <- asNetwork(Root.As.Root.ig)

net.hs <- hub_score(Root.As.Root.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Root.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Root.net %v% "hubscore", levels = c("low","high"))) 

R_As <- ggnet2(Root.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Root Bx

CMRE5_Rootb_Bxc_otu_net1 <- symBeta(getOptBeta(CMRE5_Rootb_Bxc_otu_net))
colnames(CMRE5_Rootb_Bxc_otu_net1) <- rownames(CMRE5_Rootb_Bxc_otu_net1) <- colnames(CMRE5_Rootb_Bxc_otu)
Root.ig <- adj2igraph(getRefit(CMRE5_Rootb_Bxc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Rootb_Bxc)))
vsize <- rowMeans(clr(CMRE5_Rootb_Bxc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Root.ig)

otu.ids=colnames(CMRE5_Rootb_Bxc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Bxc_otu_net)))

edges=E(Root.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Root.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Root.ig)$color=edge.colors
Root.net <- asNetwork(Root.Bx.Root.ig)

net.hs <- hub_score(Root.Bx.Root.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Root.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Root.net %v% "hubscore", levels = c("low","high"))) 

R_Bx <- ggnet2(Root.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Root Tb

CMRE5_Rootb_Tbc_otu_net1 <- symBeta(getOptBeta(CMRE5_Rootb_Tbc_otu_net))
colnames(CMRE5_Rootb_Tbc_otu_net1) <- rownames(CMRE5_Rootb_Tbc_otu_net1) <- colnames(CMRE5_Rootb_Tbc_otu)
Root.ig <- adj2igraph(getRefit(CMRE5_Rootb_Tbc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Rootb_Tbc)))
vsize <- rowMeans(clr(CMRE5_Rootb_Tbc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Root.ig)

otu.ids=colnames(CMRE5_Rootb_Tbc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Tbc_otu_net)))

edges=E(Root.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Root.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Root.ig)$color=edge.colors
Root.net <- asNetwork(Root.Tb.Root.ig)

net.hs <- hub_score(Root.Tb.Root.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Root.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Root.net %v% "hubscore", levels = c("low","high"))) 

R_Tb <- ggnet2(Root.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Water Ctr

CMRE5_Waterb_Ctrc_otu_net1 <- symBeta(getOptBeta(CMRE5_Waterb_Ctrc_otu_net))
colnames(CMRE5_Waterb_Ctrc_otu_net1) <- rownames(CMRE5_Waterb_Ctrc_otu_net1) <- colnames(CMRE5_Waterb_Ctrc_otu)
Water.ig <- adj2igraph(getRefit(CMRE5_Waterb_Ctrc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Waterb_Ctrc)))
vsize <- rowMeans(clr(CMRE5_Waterb_Ctrc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Water.ig)

otu.ids=colnames(CMRE5_Waterb_Ctrc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Ctrc_otu_net)))

edges=E(Water.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Water.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Water.ig)$color=edge.colors
Water.net <- asNetwork(Water.Ctr.Water.ig)

net.hs <- hub_score(Water.Ctr.Water.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Water.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Water.net %v% "hubscore", levels = c("low","high"))) 

W_Ctr <- ggnet2(Water.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Water As 

CMRE5_Waterb_Asc_otu_net1 <- symBeta(getOptBeta(CMRE5_Waterb_Asc_otu_net))
colnames(CMRE5_Waterb_Asc_otu_net1) <- rownames(CMRE5_Waterb_Asc_otu_net1) <- colnames(CMRE5_Waterb_Asc_otu)
Water.ig <- adj2igraph(getRefit(CMRE5_Waterb_Asc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Waterb_Asc)))
vsize <- rowMeans(clr(CMRE5_Waterb_Asc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Water.ig)

otu.ids=colnames(CMRE5_Waterb_Asc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Asc_otu_net)))

edges=E(Water.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Water.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Water.ig)$color=edge.colors
Water.net <- asNetwork(Water.As.Water.ig)

net.hs <- hub_score(Water.As.Water.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Water.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Water.net %v% "hubscore", levels = c("low","high"))) 

W_As <- ggnet2(Water.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Water Bx

CMRE5_Waterb_Bxc_otu_net1 <- symBeta(getOptBeta(CMRE5_Waterb_Bxc_otu_net))
colnames(CMRE5_Waterb_Bxc_otu_net1) <- rownames(CMRE5_Waterb_Bxc_otu_net1) <- colnames(CMRE5_Waterb_Bxc_otu)
Water.ig <- adj2igraph(getRefit(CMRE5_Waterb_Bxc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Waterb_Bxc)))
vsize <- rowMeans(clr(CMRE5_Waterb_Bxc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Water.ig)

otu.ids=colnames(CMRE5_Waterb_Bxc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Bxc_otu_net)))

edges=E(Water.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Water.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Water.ig)$color=edge.colors
Water.net <- asNetwork(Water.Bx.Water.ig)

net.hs <- hub_score(Water.Bx.Water.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Water.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Water.net %v% "hubscore", levels = c("low","high"))) 

W_Bx <- ggnet2(Water.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Water Tb

CMRE5_Waterb_Tbc_otu_net1 <- symBeta(getOptBeta(CMRE5_Waterb_Tbc_otu_net))
colnames(CMRE5_Waterb_Tbc_otu_net1) <- rownames(CMRE5_Waterb_Tbc_otu_net1) <- colnames(CMRE5_Waterb_Tbc_otu)
Water.ig <- adj2igraph(getRefit(CMRE5_Waterb_Tbc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Waterb_Tbc)))
vsize <- rowMeans(clr(CMRE5_Waterb_Tbc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Water.ig)

otu.ids=colnames(CMRE5_Waterb_Tbc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Tbc_otu_net)))

edges=E(Water.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Water.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Water.ig)$color=edge.colors
Water.net <- asNetwork(Water.Tb.Water.ig)

net.hs <- hub_score(Water.Tb.Water.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Water.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Water.net %v% "hubscore", levels = c("low","high"))) 

W_Tb <- ggnet2(Water.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Sediment Ctr

CMRE5_Sedimentb_Ctrc_otu_net1 <- symBeta(getOptBeta(CMRE5_Sedimentb_Ctrc_otu_net))
colnames(CMRE5_Sedimentb_Ctrc_otu_net1) <- rownames(CMRE5_Sedimentb_Ctrc_otu_net1) <- colnames(CMRE5_Sedimentb_Ctrc_otu)
Sediment.ig <- adj2igraph(getRefit(CMRE5_Sedimentb_Ctrc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Sedimentb_Ctrc)))
vsize <- rowMeans(clr(CMRE5_Sedimentb_Ctrc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Sediment.ig)

otu.ids=colnames(CMRE5_Sedimentb_Ctrc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Ctrc_otu_net)))

edges=E(Sediment.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Sediment.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Sediment.ig)$color=edge.colors
Sediment.net <- asNetwork(Sediment.Ctr.Sediment.ig)

net.hs <- hub_score(Sediment.Ctr.Sediment.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Sediment.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Sediment.net %v% "hubscore", levels = c("low","high"))) 

Se_Ctr <- ggnet2(Sediment.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Sediment As

CMRE5_Sedimentb_Asc_otu_net1 <- symBeta(getOptBeta(CMRE5_Sedimentb_Asc_otu_net))
colnames(CMRE5_Sedimentb_Asc_otu_net1) <- rownames(CMRE5_Sedimentb_Asc_otu_net1) <- colnames(CMRE5_Sedimentb_Asc_otu)
Sediment.ig <- adj2igraph(getRefit(CMRE5_Sedimentb_Asc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Sedimentb_Asc)))
vsize <- rowMeans(clr(CMRE5_Sedimentb_Asc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Sediment.ig)

otu.ids=colnames(CMRE5_Sedimentb_Asc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Asc_otu_net)))

edges=E(Sediment.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Sediment.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Sediment.ig)$color=edge.colors
Sediment.net <- asNetwork(Sediment.As.Sediment.ig)

net.hs <- hub_score(Sediment.As.Sediment.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Sediment.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Sediment.net %v% "hubscore", levels = c("low","high"))) 

Se_As <- ggnet2(Sediment.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Sediment Bx

CMRE5_Sedimentb_Bxc_otu_net1 <- symBeta(getOptBeta(CMRE5_Sedimentb_Bxc_otu_net))
colnames(CMRE5_Sedimentb_Bxc_otu_net1) <- rownames(CMRE5_Sedimentb_Bxc_otu_net1) <- colnames(CMRE5_Sedimentb_Bxc_otu)
Sediment.ig <- adj2igraph(getRefit(CMRE5_Sedimentb_Bxc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Sedimentb_Bxc)))
vsize <- rowMeans(clr(CMRE5_Sedimentb_Bxc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Sediment.ig)

otu.ids=colnames(CMRE5_Sedimentb_Bxc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Bxc_otu_net)))

edges=E(Sediment.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Sediment.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Sediment.ig)$color=edge.colors
Sediment.net <- asNetwork(Sediment.Bx.Sediment.ig)

net.hs <- hub_score(Sediment.Bx.Sediment.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Sediment.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Sediment.net %v% "hubscore", levels = c("low","high"))) 

Se_Bx <- ggnet2(Sediment.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))


##########-Sediment Tb

CMRE5_Sedimentb_Tbc_otu_net1 <- symBeta(getOptBeta(CMRE5_Sedimentb_Tbc_otu_net))
colnames(CMRE5_Sedimentb_Tbc_otu_net1) <- rownames(CMRE5_Sedimentb_Tbc_otu_net1) <- colnames(CMRE5_Sedimentb_Tbc_otu)
Sediment.ig <- adj2igraph(getRefit(CMRE5_Sedimentb_Tbc_otu_net),  rmEmptyNodes=T, vertex.attr=list(name=taxa_names(CMRE5_Sedimentb_Tbc)))
vsize <- rowMeans(clr(CMRE5_Sedimentb_Tbc_otu, 1))+6
am.coord <- layout.fruchterman.reingold(Sediment.ig)

otu.ids=colnames(CMRE5_Sedimentb_Tbc_otu_net1)
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Tbc_otu_net)))

edges=E(Sediment.ig)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(Sediment.ig,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"dodgerblue")
  }else if(beta<0){
    edge.colors=append(edge.colors,"lightsalmon")
  }
}

E(Sediment.ig)$color=edge.colors
Sediment.net <- asNetwork(Sediment.Tb.Sediment.ig)

net.hs <- hub_score(Sediment.Tb.Sediment.ig, weights=NA)$vector # 
net.hs1 <- ifelse(net.hs < 0.7, "low", "high") # 
Sediment.net %v% "hubscore" = ifelse(net.hs <0.7,"low","high") # 
factor_color <- sort(factor(Sediment.net %v% "hubscore", levels = c("low","high"))) 

Se_Tb <- ggnet2(Sediment.net, label = FALSE,  edge.color = "color", edge.size=0.5 , node.color = factor_color, palette = c("low" = "gray60", "high" = "dodgerblue4"), size= factor_color, size.palette = c("low" = 5, "high" = 15)) + guides(color = FALSE, size = FALSE) + theme(panel.border = element_rect(colour = "gray50", fill=NA, size=1))

####################################################
######------------Figure 4B------------############# 

############## Mouse Degree distributions ######################

Mouse.Ctr.feces.ig= readRDS("Mouse.Ctr.feces.ig.RDS")
Mouse.As.feces.ig= readRDS("Mouse.As.feces.ig.RDS")
Mouse.Bx.feces.ig= readRDS("Mouse.Bx.feces.ig.RDS")
Mouse.Tb.feces.ig= readRDS("Mouse.Tb.feces.ig.RDS")

Mouse.Ctr_all_deg <- as.data.frame(igraph::degree(Mouse.Ctr.feces.ig, mode="all")) 
Mouse.As_all_deg <- as.data.frame(igraph::degree(Mouse.As.feces.ig, mode="all"))
Mouse.Bx_all_deg <- as.data.frame(igraph::degree(Mouse.Bx.feces.ig, mode="all"))
Mouse.Tb_all_deg <- as.data.frame(igraph::degree(Mouse.Tb.feces.ig, mode="all"))

names(Mouse.Ctr_all_deg)[1]<-paste("Degree")
names(Mouse.As_all_deg)[1]<-paste("Degree")
names(Mouse.Bx_all_deg)[1]<-paste("Degree")
names(Mouse.Tb_all_deg)[1]<-paste("Degree")

Mouse.Ctr_all_deg$Treatment <- rep("Ctr",nrow(Mouse.Ctr_all_deg))
Mouse.As_all_deg$Treatment <- rep("As",nrow(Mouse.As_all_deg))
Mouse.Bx_all_deg$Treatment <- rep("Bx",nrow(Mouse.Bx_all_deg))
Mouse.Tb_all_deg$Treatment <- rep("Tb",nrow(Mouse.Tb_all_deg))

Mouse.Ctr_all_deg <- dplyr::as_tibble(Mouse.Ctr_all_deg, rownames = "ASV_ID")
Mouse.As_all_deg <- dplyr::as_tibble(Mouse.As_all_deg, rownames = "ASV_ID")
Mouse.Bx_all_deg <- dplyr::as_tibble(Mouse.Bx_all_deg, rownames = "ASV_ID")
Mouse.Tb_all_deg <- dplyr::as_tibble(Mouse.Tb_all_deg, rownames = "ASV_ID")

Mouse_all_deg <- rbind(Mouse.Ctr_all_deg, Mouse.As_all_deg, Mouse.Bx_all_deg, Mouse.Tb_all_deg)

Mouse_all_deg$Treatment<- factor(Mouse_all_deg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

Mouse_p1 <- ggplot(Mouse_all_deg, aes(x=Treatment,  y= Degree), colour = Treatment)  + geom_jitter(size=2.0, width = 0.25, aes(colour = Treatment)) + theme_set(theme_bw()) + theme(panel.border = element_rect(colour = "black", size=1), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=34, face="bold"), axis.text.y = element_text(size=24, face="bold"))  + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73",  "#D55E00")) + theme(legend.position = "none")  

############## Root Degree distributions ######################

Root.Ctr.Root.ig= readRDS("Root.Ctr.Root.ig.RDS")
Root.As.Root.ig= readRDS("Root.As.Root.ig.RDS")
Root.Bx.Root.ig= readRDS("Root.Bx.Root.ig.RDS")
Root.Tb.Root.ig= readRDS("Root.Tb.Root.ig.RDS")

Root.Ctr_all_deg <- as.data.frame(igraph::degree(Root.Ctr.Root.ig, mode="all")) 
Root.As_all_deg <- as.data.frame(igraph::degree(Root.As.Root.ig, mode="all"))
Root.Bx_all_deg <- as.data.frame(igraph::degree(Root.Bx.Root.ig, mode="all"))
Root.Tb_all_deg <- as.data.frame(igraph::degree(Root.Tb.Root.ig, mode="all"))

names(Root.Ctr_all_deg)[1]<-paste("Degree")
names(Root.As_all_deg)[1]<-paste("Degree")
names(Root.Bx_all_deg)[1]<-paste("Degree")
names(Root.Tb_all_deg)[1]<-paste("Degree")

Root.Ctr_all_deg$Treatment <- rep("Ctr",nrow(Root.Ctr_all_deg))
Root.As_all_deg$Treatment <- rep("As",nrow(Root.As_all_deg))
Root.Bx_all_deg$Treatment <- rep("Bx",nrow(Root.Bx_all_deg))
Root.Tb_all_deg$Treatment <- rep("Tb",nrow(Root.Tb_all_deg))

Root.Ctr_all_deg <- dplyr::as_tibble(Root.Ctr_all_deg, rownames = "ASV_ID")
Root.As_all_deg <- dplyr::as_tibble(Root.As_all_deg, rownames = "ASV_ID")
Root.Bx_all_deg <- dplyr::as_tibble(Root.Bx_all_deg, rownames = "ASV_ID")
Root.Tb_all_deg <- dplyr::as_tibble(Root.Tb_all_deg, rownames = "ASV_ID")

Root_all_deg <- rbind(Root.Ctr_all_deg, Root.As_all_deg, Root.Bx_all_deg, Root.Tb_all_deg)

Root_all_deg$Treatment<- factor(Root_all_deg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

Root_p1 <- ggplot(Root_all_deg, aes(x=Treatment,  y= Degree), colour = Treatment)  + geom_jitter(size=2.0, width = 0.25, aes(colour = Treatment)) + theme_set(theme_bw()) + theme(panel.border = element_rect(colour = "black", size=1), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=34, face="bold"), axis.text.y = element_text(size=24, face="bold"))  + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73",  "#D55E00")) + theme(legend.position = "none")

############## Soil Degree distributions ######################

Soil.Ctr.Soil.ig= readRDS("Soil.Ctr.Soil.ig.RDS")
Soil.As.Soil.ig= readRDS("Soil.As.Soil.ig.RDS")
Soil.Bx.Soil.ig= readRDS("Soil.Bx.Soil.ig.RDS")
Soil.Tb.Soil.ig= readRDS("Soil.Tb.Soil.ig.RDS")

Soil.Ctr_all_deg <- as.data.frame(igraph::degree(Soil.Ctr.Soil.ig, mode="all")) 
Soil.As_all_deg <- as.data.frame(igraph::degree(Soil.As.Soil.ig, mode="all"))
Soil.Bx_all_deg <- as.data.frame(igraph::degree(Soil.Bx.Soil.ig, mode="all"))
Soil.Tb_all_deg <- as.data.frame(igraph::degree(Soil.Tb.Soil.ig, mode="all"))

names(Soil.Ctr_all_deg)[1]<-paste("Degree")
names(Soil.As_all_deg)[1]<-paste("Degree")
names(Soil.Bx_all_deg)[1]<-paste("Degree")
names(Soil.Tb_all_deg)[1]<-paste("Degree")

Soil.Ctr_all_deg$Treatment <- rep("Ctr",nrow(Soil.Ctr_all_deg))
Soil.As_all_deg$Treatment <- rep("As",nrow(Soil.As_all_deg))
Soil.Bx_all_deg$Treatment <- rep("Bx",nrow(Soil.Bx_all_deg))
Soil.Tb_all_deg$Treatment <- rep("Tb",nrow(Soil.Tb_all_deg))

Soil.Ctr_all_deg <- dplyr::as_tibble(Soil.Ctr_all_deg, rownames = "ASV_ID")
Soil.As_all_deg <- dplyr::as_tibble(Soil.As_all_deg, rownames = "ASV_ID")
Soil.Bx_all_deg <- dplyr::as_tibble(Soil.Bx_all_deg, rownames = "ASV_ID")
Soil.Tb_all_deg <- dplyr::as_tibble(Soil.Tb_all_deg, rownames = "ASV_ID")

Soil_all_deg <- rbind(Soil.Ctr_all_deg, Soil.As_all_deg, Soil.Bx_all_deg, Soil.Tb_all_deg)
Soil_all_deg$Treatment<- factor(Soil_all_deg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

Soil_p1 <- ggplot(Soil_all_deg, aes(x=Treatment,  y= Degree), colour = Treatment)  + geom_jitter(size=2.0, width = 0.25, aes(colour = Treatment)) + theme_set(theme_bw())+ theme(panel.border = element_rect(colour = "black", size=1), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=34, face="bold"), axis.text.y = element_text(size=24, face="bold"))  + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73",  "#D55E00")) + theme(legend.position = "none") + scale_y_continuous(breaks = seq(0, 70, by=10), limits=c(0,70))

############## Water Degree distributions ######################

Water.Ctr.Water.ig= readRDS("Water.Ctr.Water.ig.RDS")
Water.As.Water.ig= readRDS("Water.As.Water.ig.RDS")
Water.Bx.Water.ig= readRDS("Water.Bx.Water.ig.RDS")
Water.Tb.Water.ig= readRDS("Water.Tb.Water.ig.RDS")

Water.Ctr_all_deg <- as.data.frame(igraph::degree(Water.Ctr.Water.ig, mode="all")) 
Water.As_all_deg <- as.data.frame(igraph::degree(Water.As.Water.ig, mode="all"))
Water.Bx_all_deg <- as.data.frame(igraph::degree(Water.Bx.Water.ig, mode="all"))
Water.Tb_all_deg <- as.data.frame(igraph::degree(Water.Tb.Water.ig, mode="all"))

names(Water.Ctr_all_deg)[1]<-paste("Degree")
names(Water.As_all_deg)[1]<-paste("Degree")
names(Water.Bx_all_deg)[1]<-paste("Degree")
names(Water.Tb_all_deg)[1]<-paste("Degree")

Water.Ctr_all_deg$Treatment <- rep("Ctr",nrow(Water.Ctr_all_deg))
Water.As_all_deg$Treatment <- rep("As",nrow(Water.As_all_deg))
Water.Bx_all_deg$Treatment <- rep("Bx",nrow(Water.Bx_all_deg))
Water.Tb_all_deg$Treatment <- rep("Tb",nrow(Water.Tb_all_deg))

Water.Ctr_all_deg <- dplyr::as_tibble(Water.Ctr_all_deg, rownames = "ASV_ID")
Water.As_all_deg <- dplyr::as_tibble(Water.As_all_deg, rownames = "ASV_ID")
Water.Bx_all_deg <- dplyr::as_tibble(Water.Bx_all_deg, rownames = "ASV_ID")
Water.Tb_all_deg <- dplyr::as_tibble(Water.Tb_all_deg, rownames = "ASV_ID")

Water_all_deg <- rbind(Water.Ctr_all_deg, Water.As_all_deg, Water.Bx_all_deg, Water.Tb_all_deg)
Water_all_deg$Treatment<- factor(Water_all_deg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

Water_p1 <- ggplot(Water_all_deg, aes(x=Treatment,  y= Degree), colour = Treatment)  + geom_jitter(size=2.0, width = 0.25, aes(colour = Treatment)) + theme_set(theme_bw()) + theme(panel.border = element_rect(colour = "black",  size=1), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=34, face="bold"), axis.text.y = element_text(size=24, face="bold"))  + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73",  "#D55E00")) + theme(legend.position = "none") + scale_y_continuous(breaks = seq(0, 35, by=10), limits=c(0,35))

############## Sediment Degree distributions ######################

Sediment.Ctr.Sediment.ig= readRDS("Sediment.Ctr.Sediment.ig.RDS")
Sediment.As.Sediment.ig= readRDS("Sediment.As.Sediment.ig.RDS")
Sediment.Bx.Sediment.ig= readRDS("Sediment.Bx.Sediment.ig.RDS")
Sediment.Tb.Sediment.ig= readRDS("Sediment.Tb.Sediment.ig.RDS")

Sediment.Ctr_all_deg <- as.data.frame(igraph::degree(Sediment.Ctr.Sediment.ig, mode="all"))
Sediment.As_all_deg <- as.data.frame(igraph::degree(Sediment.As.Sediment.ig, mode="all"))
Sediment.Bx_all_deg <- as.data.frame(igraph::degree(Sediment.Bx.Sediment.ig, mode="all"))
Sediment.Tb_all_deg <- as.data.frame(igraph::degree(Sediment.Tb.Sediment.ig, mode="all"))

names(Sediment.Ctr_all_deg)[1]<-paste("Degree")
names(Sediment.As_all_deg)[1]<-paste("Degree")
names(Sediment.Bx_all_deg)[1]<-paste("Degree")
names(Sediment.Tb_all_deg)[1]<-paste("Degree")

Sediment.Ctr_all_deg$Treatment <- rep("Ctr",nrow(Sediment.Ctr_all_deg))
Sediment.As_all_deg$Treatment <- rep("As",nrow(Sediment.As_all_deg))
Sediment.Bx_all_deg$Treatment <- rep("Bx",nrow(Sediment.Bx_all_deg))
Sediment.Tb_all_deg$Treatment <- rep("Tb",nrow(Sediment.Tb_all_deg))

Sediment.Ctr_all_deg <- dplyr::as_tibble(Sediment.Ctr_all_deg, rownames = "ASV_ID")
Sediment.As_all_deg <- dplyr::as_tibble(Sediment.As_all_deg, rownames = "ASV_ID")
Sediment.Bx_all_deg <- dplyr::as_tibble(Sediment.Bx_all_deg, rownames = "ASV_ID")
Sediment.Tb_all_deg <- dplyr::as_tibble(Sediment.Tb_all_deg, rownames = "ASV_ID")

Sediment_all_deg <- rbind(Sediment.Ctr_all_deg, Sediment.As_all_deg, Sediment.Bx_all_deg, Sediment.Tb_all_deg)
Sediment_all_deg$Treatment<- factor(Sediment_all_deg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

Sediment_p1 <- ggplot(Sediment_all_deg, aes(x=Treatment,  y= Degree), colour = Treatment)  + geom_jitter(size=2.0, width = 0.25, aes(colour = Treatment)) + theme_set(theme_bw()) + theme(panel.border = element_rect(colour = "black", size=1), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size=34, face="bold"), axis.text.y = element_text(size=24, face="bold"))  + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73",  "#D55E00")) + theme(legend.position = "none")

############################################################
######combined plot ####################

Figure_4 <- ggarrange(W_Ctr, Se_Ctr, S_Ctr, R_Ctr, M_Ctr, W_As, Se_As, S_As, R_As,  M_As,  W_Bx, Se_bx, S_Bx, R_Bx, M_Bx, W_Tb, Se_Tb, S_Tb, R_Tb, M_Tb, Water_p1, Sediment_p1, Soil_p1, Root_p1, Mouse_p1, ncol = 5, nrow = 5)


##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure 4---------###############################################

###################--NMDS plot

set.seed(101)# set seed for reproducibility
CMRE9_r8100_t_br.ord <- ordinate(CMRE9_r8100_t , "NMDS", "bray") 
# Before making plot, order variables
sample_data(CMRE9_r8100_t)$Sampletype<- factor(sample_data(CMRE9_r8100_t)$Sampletype, levels = c("water", "sedi", "soil", "roots", "feces"))
CMRE9_r8100_plot_bray_beta_diversity = plot_ordination(CMRE9_r8100_t, CMRE9_r8100_t_br.ord, color="Sampletype")  + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=21), axis.title.y = element_text(size=21), axis.text.x = element_text(size=19), axis.text.y = element_text(size=19)) + geom_point(size=2.50) + theme(legend.text=element_text(size=18), legend.title=element_text(size=20)) + scale_color_manual(labels = c("Water", "Sediment", "Soil", "Root", "Mouse"), values=c("steelblue3", "slategray4",  "salmon4", "darkolivegreen4", "tan3"))

##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure 6---------###############################################

###################--Beta Diversity-Box plots 

CMRE9_r8100_t  = readRDS("CMRE9_r8100_t.RDS")

CMRE9_r8100_t_feces <- subset_samples(CMRE9_r8100_t, Sampletype%in%c("feces"))
CMRE9_r8100_t_roots <- subset_samples(CMRE9_r8100_t, Sampletype%in%c("roots"))
CMRE9_r8100_t_soil <- subset_samples(CMRE9_r8100_t, Sampletype%in%c("soil"))
CMRE9_r8100_t_water <- subset_samples(CMRE9_r8100_t, Sampletype%in%c("water"))
CMRE9_r8100_t_sedi <- subset_samples(CMRE9_r8100_t, Sampletype%in%c("sedi"))

###############################--------- Bray-Curtis------------------####################

#######--For mouse -----#######

CMRE9_r8100_t_feces_As <- subset_samples(CMRE9_r8100_t_feces, Treatment%in%c("As"))
CMRE9_r8100_t_feces_Bx <- subset_samples(CMRE9_r8100_t_feces, Treatment%in%c("Bx"))
CMRE9_r8100_t_feces_Tb <- subset_samples(CMRE9_r8100_t_feces, Treatment%in%c("Tb"))
CMRE9_r8100_t_feces_Ctr <- subset_samples(CMRE9_r8100_t_feces, Treatment%in%c("Ctr"))

CMRE9_r8100_t_feces_As_br <- distance(CMRE9_r8100_t_feces_As, method = "bray") #bray-curtis 
CMRE9_r8100_t_feces_Bx_br <- distance(CMRE9_r8100_t_feces_Bx, method = "bray")
CMRE9_r8100_t_feces_Tb_br <- distance(CMRE9_r8100_t_feces_Tb, method = "bray")
CMRE9_r8100_t_feces_Ctr_br <- distance(CMRE9_r8100_t_feces_Ctr, method = "bray")

CMRE9_r8100_t_feces_As_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_feces_As_br))) ###converting distance matrix to list 
CMRE9_r8100_feces_As_br_pairlist <- CMRE9_r8100_t_feces_As_br_pairlist[as.numeric(CMRE9_r8100_t_feces_As_br_pairlist$col) > as.numeric(CMRE9_r8100_t_feces_As_br_pairlist$row),]
CMRE9_r8100_feces_As_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE9_r8100_feces_As_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE9_r8100_feces_As_br_pairlist$Treatment <- rep("As",nrow(CMRE9_r8100_feces_As_br_pairlist)) #make new column Sampletype and fill all Mouse

##

CMRE9_r8100_t_feces_Bx_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_feces_Bx_br))) ###converting distance matrix to list 
CMRE9_r8100_feces_Bx_br_pairlist <- CMRE9_r8100_t_feces_Bx_br_pairlist[as.numeric(CMRE9_r8100_t_feces_Bx_br_pairlist$col) > as.numeric(CMRE9_r8100_t_feces_Bx_br_pairlist$row),]
CMRE9_r8100_feces_Bx_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE9_r8100_feces_Bx_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE9_r8100_feces_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE9_r8100_feces_Bx_br_pairlist)) #make new column Sampletype and fill all Mouse

##

CMRE9_r8100_t_feces_Tb_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_feces_Tb_br))) ###converting distance matrix to list 
CMRE9_r8100_feces_Tb_br_pairlist <- CMRE9_r8100_t_feces_Tb_br_pairlist[as.numeric(CMRE9_r8100_t_feces_Tb_br_pairlist$col) > as.numeric(CMRE9_r8100_t_feces_Tb_br_pairlist$row),]
CMRE9_r8100_feces_Tb_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE9_r8100_feces_Tb_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE9_r8100_feces_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE9_r8100_feces_Tb_br_pairlist)) #make new column Sampletype and fill all Mouse

##

CMRE9_r8100_t_feces_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_feces_Ctr_br))) ###converting distance matrix to list 
CMRE9_r8100_feces_Ctr_br_pairlist <- CMRE9_r8100_t_feces_Ctr_br_pairlist[as.numeric(CMRE9_r8100_t_feces_Ctr_br_pairlist$col) > as.numeric(CMRE9_r8100_t_feces_Ctr_br_pairlist$row),]
CMRE9_r8100_feces_Ctr_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE9_r8100_feces_Ctr_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE9_r8100_feces_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE9_r8100_feces_Ctr_br_pairlist)) #make new column Sampletype and fill all Mouse

##################now merge all treatment csv files

CMRE9_r8100_feces_pairlist <- rbind(CMRE9_r8100_feces_Ctr_br_pairlist, CMRE9_r8100_feces_As_br_pairlist, CMRE9_r8100_feces_Bx_br_pairlist, CMRE9_r8100_feces_Tb_br_pairlist)

######################################

#######--For roots -----####### 

CMRE9_r8100_t_roots_As <- subset_samples(CMRE9_r8100_t_roots, Treatment%in%c("As"))
CMRE9_r8100_t_roots_Bx <- subset_samples(CMRE9_r8100_t_roots, Treatment%in%c("Bx"))
CMRE9_r8100_t_roots_Tb <- subset_samples(CMRE9_r8100_t_roots, Treatment%in%c("Tb"))
CMRE9_r8100_t_roots_Ctr <- subset_samples(CMRE9_r8100_t_roots, Treatment%in%c("Ctr"))

CMRE9_r8100_t_roots_As_br <- distance(CMRE9_r8100_t_roots_As, method = "bray") #bray-curtis 
CMRE9_r8100_t_roots_Bx_br <- distance(CMRE9_r8100_t_roots_Bx, method = "bray")
CMRE9_r8100_t_roots_Tb_br <- distance(CMRE9_r8100_t_roots_Tb, method = "bray")
CMRE9_r8100_t_roots_Ctr_br <- distance(CMRE9_r8100_t_roots_Ctr, method = "bray")

CMRE9_r8100_t_roots_As_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_roots_As_br))) ###converting distance matrix to list 
CMRE9_r8100_roots_As_br_pairlist <- CMRE9_r8100_t_roots_As_br_pairlist[as.numeric(CMRE9_r8100_t_roots_As_br_pairlist$col) > as.numeric(CMRE9_r8100_t_roots_As_br_pairlist$row),]
CMRE9_r8100_roots_As_br_pairlist$Sampletype <- rep("Root",nrow(CMRE9_r8100_roots_As_br_pairlist)) #make new column Sampletype and fill all Root
CMRE9_r8100_roots_As_br_pairlist$Treatment <- rep("As",nrow(CMRE9_r8100_roots_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE9_r8100_t_roots_Bx_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_roots_Bx_br))) ###converting distance matrix to list 
CMRE9_r8100_roots_Bx_br_pairlist <- CMRE9_r8100_t_roots_Bx_br_pairlist[as.numeric(CMRE9_r8100_t_roots_Bx_br_pairlist$col) > as.numeric(CMRE9_r8100_t_roots_Bx_br_pairlist$row),]
CMRE9_r8100_roots_Bx_br_pairlist$Sampletype <- rep("Root",nrow(CMRE9_r8100_roots_Bx_br_pairlist)) #make new column Sampletype and fill all Root
CMRE9_r8100_roots_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE9_r8100_roots_Bx_br_pairlist)) #make new column Sampletype and fill all Bx

##

CMRE9_r8100_t_roots_Tb_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_roots_Tb_br))) ###converting distance matrix to list 
CMRE9_r8100_roots_Tb_br_pairlist <- CMRE9_r8100_t_roots_Tb_br_pairlist[as.numeric(CMRE9_r8100_t_roots_Tb_br_pairlist$col) > as.numeric(CMRE9_r8100_t_roots_Tb_br_pairlist$row),]
CMRE9_r8100_roots_Tb_br_pairlist$Sampletype <- rep("Root",nrow(CMRE9_r8100_roots_Tb_br_pairlist)) #make new column Sampletype and fill all Root
CMRE9_r8100_roots_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE9_r8100_roots_Tb_br_pairlist)) #make new column Sampletype and fill all Tb

##

CMRE9_r8100_t_roots_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_roots_Ctr_br))) ###converting distance matrix to list 
CMRE9_r8100_roots_Ctr_br_pairlist <- CMRE9_r8100_t_roots_Ctr_br_pairlist[as.numeric(CMRE9_r8100_t_roots_Ctr_br_pairlist$col) > as.numeric(CMRE9_r8100_t_roots_Ctr_br_pairlist$row),]
CMRE9_r8100_roots_Ctr_br_pairlist$Sampletype <- rep("Root",nrow(CMRE9_r8100_roots_Ctr_br_pairlist)) #make new column Sampletype and fill all Root
CMRE9_r8100_roots_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE9_r8100_roots_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE9_r8100_roots_pairlist <- rbind(CMRE9_r8100_roots_Ctr_br_pairlist, CMRE9_r8100_roots_As_br_pairlist, CMRE9_r8100_roots_Bx_br_pairlist, CMRE9_r8100_roots_Tb_br_pairlist)

######################################

#######--For soil -----#######

CMRE9_r8100_t_soil_As <- subset_samples(CMRE9_r8100_t_soil, Treatment%in%c("As"))
CMRE9_r8100_t_soil_Bx <- subset_samples(CMRE9_r8100_t_soil, Treatment%in%c("Bx"))
CMRE9_r8100_t_soil_Tb <- subset_samples(CMRE9_r8100_t_soil, Treatment%in%c("Tb"))
CMRE9_r8100_t_soil_Ctr <- subset_samples(CMRE9_r8100_t_soil, Treatment%in%c("Ctr"))

CMRE9_r8100_t_soil_As_br <- distance(CMRE9_r8100_t_soil_As, method = "bray") #bray-curtis 
CMRE9_r8100_t_soil_Bx_br <- distance(CMRE9_r8100_t_soil_Bx, method = "bray")
CMRE9_r8100_t_soil_Tb_br <- distance(CMRE9_r8100_t_soil_Tb, method = "bray")
CMRE9_r8100_t_soil_Ctr_br <- distance(CMRE9_r8100_t_soil_Ctr, method = "bray")

CMRE9_r8100_t_soil_As_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_soil_As_br))) ###converting distance matrix to list 
CMRE9_r8100_soil_As_br_pairlist <- CMRE9_r8100_t_soil_As_br_pairlist[as.numeric(CMRE9_r8100_t_soil_As_br_pairlist$col) > as.numeric(CMRE9_r8100_t_soil_As_br_pairlist$row),]
CMRE9_r8100_soil_As_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE9_r8100_soil_As_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE9_r8100_soil_As_br_pairlist$Treatment <- rep("As",nrow(CMRE9_r8100_soil_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE9_r8100_t_soil_Bx_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_soil_Bx_br))) ###converting distance matrix to list 
CMRE9_r8100_soil_Bx_br_pairlist <- CMRE9_r8100_t_soil_Bx_br_pairlist[as.numeric(CMRE9_r8100_t_soil_Bx_br_pairlist$col) > as.numeric(CMRE9_r8100_t_soil_Bx_br_pairlist$row),]
CMRE9_r8100_soil_Bx_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE9_r8100_soil_Bx_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE9_r8100_soil_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE9_r8100_soil_Bx_br_pairlist)) #make new column Sampletype and fill all Bx
write.csv(CMRE9_r8100_soil_Bx_br_pairlist, file = "CMRE9_r8100_soil_Bx_br_pairlist.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save as CSV file
CMRE9_r8100_soil_Bx_br_pairlist = read.csv(file = "CMRE9_r8100_soil_Bx_br_pairlist.CSV", header = TRUE)

##

CMRE9_r8100_t_soil_Tb_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_soil_Tb_br))) ###converting distance matrix to list 
CMRE9_r8100_soil_Tb_br_pairlist <- CMRE9_r8100_t_soil_Tb_br_pairlist[as.numeric(CMRE9_r8100_t_soil_Tb_br_pairlist$col) > as.numeric(CMRE9_r8100_t_soil_Tb_br_pairlist$row),]
CMRE9_r8100_soil_Tb_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE9_r8100_soil_Tb_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE9_r8100_soil_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE9_r8100_soil_Tb_br_pairlist)) #make new column Sampletype and fill all Tb
write.csv(CMRE9_r8100_soil_Tb_br_pairlist, file = "CMRE9_r8100_soil_Tb_br_pairlist.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save as CSV file
CMRE9_r8100_soil_Tb_br_pairlist = read.csv(file = "CMRE9_r8100_soil_Tb_br_pairlist.CSV", header = TRUE)

##

CMRE9_r8100_t_soil_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_soil_Ctr_br))) ###converting distance matrix to list 
CMRE9_r8100_soil_Ctr_br_pairlist <- CMRE9_r8100_t_soil_Ctr_br_pairlist[as.numeric(CMRE9_r8100_t_soil_Ctr_br_pairlist$col) > as.numeric(CMRE9_r8100_t_soil_Ctr_br_pairlist$row),]
CMRE9_r8100_soil_Ctr_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE9_r8100_soil_Ctr_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE9_r8100_soil_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE9_r8100_soil_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE9_r8100_soil_pairlist <- rbind(CMRE9_r8100_soil_Ctr_br_pairlist, CMRE9_r8100_soil_As_br_pairlist, CMRE9_r8100_soil_Bx_br_pairlist, CMRE9_r8100_soil_Tb_br_pairlist)

######################################

#######--For water -----#######

CMRE9_r8100_t_water_As <- subset_samples(CMRE9_r8100_t_water, Treatment%in%c("As"))
CMRE9_r8100_t_water_Bx <- subset_samples(CMRE9_r8100_t_water, Treatment%in%c("Bx"))
CMRE9_r8100_t_water_Tb <- subset_samples(CMRE9_r8100_t_water, Treatment%in%c("Tb"))
CMRE9_r8100_t_water_Ctr <- subset_samples(CMRE9_r8100_t_water, Treatment%in%c("Ctr"))

CMRE9_r8100_t_water_As_br <- distance(CMRE9_r8100_t_water_As, method = "bray") #bray-curtis 
CMRE9_r8100_t_water_Bx_br <- distance(CMRE9_r8100_t_water_Bx, method = "bray")
CMRE9_r8100_t_water_Tb_br <- distance(CMRE9_r8100_t_water_Tb, method = "bray")
CMRE9_r8100_t_water_Ctr_br <- distance(CMRE9_r8100_t_water_Ctr, method = "bray")

CMRE9_r8100_t_water_As_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_water_As_br))) ###converting distance matrix to list 
CMRE9_r8100_water_As_br_pairlist <- CMRE9_r8100_t_water_As_br_pairlist[as.numeric(CMRE9_r8100_t_water_As_br_pairlist$col) > as.numeric(CMRE9_r8100_t_water_As_br_pairlist$row),]
CMRE9_r8100_water_As_br_pairlist$Sampletype <- rep("Water",nrow(CMRE9_r8100_water_As_br_pairlist)) #make new column Sampletype and fill all Water
CMRE9_r8100_water_As_br_pairlist$Treatment <- rep("As",nrow(CMRE9_r8100_water_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE9_r8100_t_water_Bx_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_water_Bx_br))) ###converting distance matrix to list 
CMRE9_r8100_water_Bx_br_pairlist <- CMRE9_r8100_t_water_Bx_br_pairlist[as.numeric(CMRE9_r8100_t_water_Bx_br_pairlist$col) > as.numeric(CMRE9_r8100_t_water_Bx_br_pairlist$row),]
CMRE9_r8100_water_Bx_br_pairlist$Sampletype <- rep("Water",nrow(CMRE9_r8100_water_Bx_br_pairlist)) #make new column Sampletype and fill all Water
CMRE9_r8100_water_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE9_r8100_water_Bx_br_pairlist)) #make new column Sampletype and fill all Bx

##

CMRE9_r8100_t_water_Tb_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_water_Tb_br))) ###converting distance matrix to list 
CMRE9_r8100_water_Tb_br_pairlist <- CMRE9_r8100_t_water_Tb_br_pairlist[as.numeric(CMRE9_r8100_t_water_Tb_br_pairlist$col) > as.numeric(CMRE9_r8100_t_water_Tb_br_pairlist$row),]
CMRE9_r8100_water_Tb_br_pairlist$Sampletype <- rep("Water",nrow(CMRE9_r8100_water_Tb_br_pairlist)) #make new column Sampletype and fill all Water
CMRE9_r8100_water_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE9_r8100_water_Tb_br_pairlist)) #make new column Sampletype and fill all Tb

##

CMRE9_r8100_t_water_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_water_Ctr_br))) ###converting distance matrix to list 
CMRE9_r8100_water_Ctr_br_pairlist <- CMRE9_r8100_t_water_Ctr_br_pairlist[as.numeric(CMRE9_r8100_t_water_Ctr_br_pairlist$col) > as.numeric(CMRE9_r8100_t_water_Ctr_br_pairlist$row),]
CMRE9_r8100_water_Ctr_br_pairlist$Sampletype <- rep("Water",nrow(CMRE9_r8100_water_Ctr_br_pairlist)) #make new column Sampletype and fill all Water
CMRE9_r8100_water_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE9_r8100_water_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE9_r8100_water_pairlist <- rbind(CMRE9_r8100_water_Ctr_br_pairlist, CMRE9_r8100_water_As_br_pairlist, CMRE9_r8100_water_Bx_br_pairlist, CMRE9_r8100_water_Tb_br_pairlist)

######################################

#######--For sediment -----#######

CMRE9_r8100_t_sedi_As <- subset_samples(CMRE9_r8100_t_sedi, Treatment%in%c("As"))
CMRE9_r8100_t_sedi_Bx <- subset_samples(CMRE9_r8100_t_sedi, Treatment%in%c("Bx"))
CMRE9_r8100_t_sedi_Tb <- subset_samples(CMRE9_r8100_t_sedi, Treatment%in%c("Tb"))
CMRE9_r8100_t_sedi_Ctr <- subset_samples(CMRE9_r8100_t_sedi, Treatment%in%c("Ctr"))

CMRE9_r8100_t_sedi_As_br <- distance(CMRE9_r8100_t_sedi_As, method = "bray") #bray-curtis 
CMRE9_r8100_t_sedi_Bx_br <- distance(CMRE9_r8100_t_sedi_Bx, method = "bray")
CMRE9_r8100_t_sedi_Tb_br <- distance(CMRE9_r8100_t_sedi_Tb, method = "bray")
CMRE9_r8100_t_sedi_Ctr_br <- distance(CMRE9_r8100_t_sedi_Ctr, method = "bray")

CMRE9_r8100_t_sedi_As_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_sedi_As_br))) ###converting distance matrix to list 
CMRE9_r8100_sedi_As_br_pairlist <- CMRE9_r8100_t_sedi_As_br_pairlist[as.numeric(CMRE9_r8100_t_sedi_As_br_pairlist$col) > as.numeric(CMRE9_r8100_t_sedi_As_br_pairlist$row),]
CMRE9_r8100_sedi_As_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE9_r8100_sedi_As_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE9_r8100_sedi_As_br_pairlist$Treatment <- rep("As",nrow(CMRE9_r8100_sedi_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE9_r8100_t_sedi_Bx_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_sedi_Bx_br))) ###converting distance matrix to list 
CMRE9_r8100_sedi_Bx_br_pairlist <- CMRE9_r8100_t_sedi_Bx_br_pairlist[as.numeric(CMRE9_r8100_t_sedi_Bx_br_pairlist$col) > as.numeric(CMRE9_r8100_t_sedi_Bx_br_pairlist$row),]
CMRE9_r8100_sedi_Bx_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE9_r8100_sedi_Bx_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE9_r8100_sedi_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE9_r8100_sedi_Bx_br_pairlist)) #make new column Sampletype and fill all Bx

##

CMRE9_r8100_t_sedi_Tb_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_sedi_Tb_br))) ###converting distance matrix to list 
CMRE9_r8100_sedi_Tb_br_pairlist <- CMRE9_r8100_t_sedi_Tb_br_pairlist[as.numeric(CMRE9_r8100_t_sedi_Tb_br_pairlist$col) > as.numeric(CMRE9_r8100_t_sedi_Tb_br_pairlist$row),]
CMRE9_r8100_sedi_Tb_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE9_r8100_sedi_Tb_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE9_r8100_sedi_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE9_r8100_sedi_Tb_br_pairlist)) #make new column Sampletype and fill all Tb

##

CMRE9_r8100_t_sedi_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE9_r8100_t_sedi_Ctr_br))) ###converting distance matrix to list 
CMRE9_r8100_sedi_Ctr_br_pairlist <- CMRE9_r8100_t_sedi_Ctr_br_pairlist[as.numeric(CMRE9_r8100_t_sedi_Ctr_br_pairlist$col) > as.numeric(CMRE9_r8100_t_sedi_Ctr_br_pairlist$row),]
CMRE9_r8100_sedi_Ctr_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE9_r8100_sedi_Ctr_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE9_r8100_sedi_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE9_r8100_sedi_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE9_r8100_sedi_pairlist <- rbind(CMRE9_r8100_sedi_Ctr_br_pairlist, CMRE9_r8100_sedi_As_br_pairlist, CMRE9_r8100_sedi_Bx_br_pairlist, CMRE9_r8100_sedi_Tb_br_pairlist)

######################################

##################-------now merge all Sampletype csv files-------###############

CMRE9_r8100_pairlist <- rbind(CMRE9_r8100_feces_pairlist, CMRE9_r8100_roots_pairlist, CMRE9_r8100_soil_pairlist, CMRE9_r8100_water_pairlist, CMRE9_r8100_sedi_pairlist)

######################################---Plot Bray------######################################

CMRE9_r8100_pairlist$Sampletype<- factor(CMRE9_r8100_pairlist$Sampletype, levels = c("Water",  "Sediment", "Soil", "Root",  "Mouse"))
CMRE9_r8100_pairlist$Treatment<- factor(CMRE9_r8100_pairlist$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE9_r8100_pairlist_bray_boxplot <- ggplot(CMRE9_r8100_pairlist, aes(Sampletype, value, fill = CMRE9_r8100_pairlist$Treatment)) + geom_boxplot() + scale_fill_manual(values=c("gray50","#E69F00", "#009E73", "#D55E00")) + labs(x= " ", y = "Bray-Curtis") + theme_set(theme_bw()) + theme(legend.position="none", axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x  = element_text(colour="black", vjust=0.5, size=16), axis.text.y  = element_text(colour="black", vjust=0.5, size=16)) + scale_y_continuous(breaks = seq(0.1, 0.9, by=0.15), limits=c(0.1,0.9)) 
CMRE9_r8100_pairlist_bray_boxplot1 <- CMRE9_r8100_pairlist_bray_boxplot + scale_x_discrete(labels=c("Ctr As Bx Tb \nWater", "Ctr As Bx Tb \nSediment", "Ctr As Bx Tb \nSoil", "Ctr As Bx Tb  \nPlant", "Ctr As Bx Tb \nAnimal"))
CMRE9_r8100_pairlist_bray_boxplot2 <- CMRE9_r8100_pairlist_bray_boxplot1 + stat_summary(fun.y=mean, geom="point", colour="black", shape=18, size=2, position=position_dodge(width=0.75))  

##############################################################################################

##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure 9---------###############################################

###################--Network degree distribution

########## Degree Distribution Figure (supplementary) for article for all compartments #####################


#### Mouse Degree distributions 

Mouse.Ctr.feces.ig= readRDS("Mouse.Ctr.feces.ig.RDS")
Mouse.As.feces.ig= readRDS("Mouse.As.feces.ig.RDS")
Mouse.Bx.feces.ig= readRDS("Mouse.Bx.feces.ig.RDS")
Mouse.Tb.feces.ig= readRDS("Mouse.Tb.feces.ig.RDS")

dd.Mouse.Ctr.feces.ig <- degree.distribution(Mouse.Ctr.feces.ig)
dd.Mouse.As.feces.ig <- degree.distribution(Mouse.As.feces.ig)
dd.Mouse.Bx.feces.ig <- degree.distribution(Mouse.Bx.feces.ig)
dd.Mouse.Tb.feces.ig <- degree.distribution(Mouse.Tb.feces.ig)

tiff(filename="Mouse_Degree_Ctr-As-Bx-Tb6.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,3.9,1,1))
plot(dd.Mouse.Ctr.feces.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.2, cex.axis=1.5, xlim=c(0, 72), ylim=c(0, 0.22), xlab="Degree", ylab="Frequency")
lines(dd.Mouse.Ctr.feces.ig, lwd = 2, col = "gray50")
lines(dd.Mouse.As.feces.ig, lwd = 2, col = "#E69F00")
lines(dd.Mouse.Bx.feces.ig, lwd = 2, col = "#009E73")
lines(dd.Mouse.Tb.feces.ig, lwd = 2, col = "#D55E00")
dev.off()

tiff(filename="Mouse_Degree_Ctr-As-Bx-Tb5_Zoom.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,4.1,1,1))
plot(dd.Mouse.Ctr.feces.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.5, cex.axis=1.8, xlim=c(10, 33), ylim=c(0, 0.10), xlab="Degree", ylab="Frequency")
lines(dd.Mouse.Ctr.feces.ig, lwd = 4, col = "gray50")
lines(dd.Mouse.As.feces.ig, lwd = 4, col = "#E69F00")
lines(dd.Mouse.Bx.feces.ig, lwd = 4, col = "#009E73")
lines(dd.Mouse.Tb.feces.ig, lwd = 4, col = "#D55E00")
dev.off()

#### Soil Degree distributions 
CMRE8_Soilb_Ctrc_Soil_ig= readRDS("CMRE8_Soilb_Ctrc_Soil_ig.RDS")
CMRE8_Soilb_Asc_Soil_ig= readRDS("CMRE8_Soilb_Asc_Soil_ig.RDS")
CMRE8_Soilb_Bxc_Soil_ig= readRDS("CMRE8_Soilb_Bxc_Soil_ig.RDS")
CMRE8_Soilb_Tbc_Soil_ig= readRDS("CMRE8_Soilb_Tbc_Soil_ig.RDS")

dd.Soil.Ctr.Soil.ig <- degree.distribution(CMRE8_Soilb_Ctrc_Soil_ig)
dd.Soil.As.Soil.ig <- degree.distribution(CMRE8_Soilb_Asc_Soil_ig)
dd.Soil.Bx.Soil.ig <- degree.distribution(CMRE8_Soilb_Bxc_Soil_ig)
dd.Soil.Tb.Soil.ig <- degree.distribution(CMRE8_Soilb_Tbc_Soil_ig)

tiff(filename="Soil_Degree_Ctr-As-Bx-Tb5.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,3.9,1,1))
plot(dd.Soil.Ctr.Soil.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.2, cex.axis=1.5, xlim=c(0, 72), ylim=c(0, 0.22), xlab="Degree", ylab="Frequency")
lines(dd.Soil.Ctr.Soil.ig, lwd = 2, col = "gray50")
lines(dd.Soil.As.Soil.ig, lwd = 2, col = "#E69F00")
lines(dd.Soil.Bx.Soil.ig, lwd = 2, col = "#009E73")
lines(dd.Soil.Tb.Soil.ig, lwd = 2, col = "#D55E00")
dev.off()

tiff(filename="Soil_Degree_Ctr-As-Bx-Tb5_Zoom.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,4.1,1,1))
plot(dd.Soil.Ctr.Soil.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.5, cex.axis=1.8, xlim=c(50, 72), ylim=c(0, 0.10), xlab="Degree", ylab="Frequency")
lines(dd.Soil.Ctr.Soil.ig, lwd = 4, col = "gray50")
lines(dd.Soil.As.Soil.ig, lwd = 4, col = "#E69F00")
lines(dd.Soil.Bx.Soil.ig, lwd = 4, col = "#009E73")
lines(dd.Soil.Tb.Soil.ig, lwd = 4, col = "#D55E00")
dev.off()

#### Root Degree distributions 
CMRE8_Rootb_Ctrc_Root_ig= readRDS("CMRE8_Rootb_Ctrc_Root_ig.RDS")
CMRE8_Rootb_Asc_Root_ig= readRDS("CMRE8_Rootb_Asc_Root_ig.RDS")
CMRE8_Rootb_Bxc_Root_ig= readRDS("CMRE8_Rootb_Bxc_Root_ig.RDS")
CMRE8_Rootb_Tbc_Root_ig= readRDS("CMRE8_Rootb_Tbc_Root_ig.RDS")

dd.Root.Ctr.Root.ig <- degree.distribution(CMRE8_Rootb_Ctrc_Root_ig)
dd.Root.As.Root.ig <- degree.distribution(CMRE8_Rootb_Asc_Root_ig)
dd.Root.Bx.Root.ig <- degree.distribution(CMRE8_Rootb_Bxc_Root_ig)
dd.Root.Tb.Root.ig <- degree.distribution(CMRE8_Rootb_Tbc_Root_ig)

tiff(filename="Root_Degree_Ctr-As-Bx-Tb5.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,3.9,1,1))
plot(dd.Root.Ctr.Root.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.2, cex.axis=1.5, xlim=c(0, 72), ylim=c(0, 0.22), xlab="Degree", ylab="Frequency")
lines(dd.Root.Ctr.Root.ig, lwd = 2, col = "gray50")
lines(dd.Root.As.Root.ig, lwd = 2, col = "#E69F00")
lines(dd.Root.Bx.Root.ig, lwd = 2, col = "#009E73")
lines(dd.Root.Tb.Root.ig, lwd = 2, col = "#D55E00")
dev.off()

tiff(filename="Root_Degree_Ctr-As-Bx-Tb5-zoom.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,4.1,1,1))
plot(dd.Root.Ctr.Root.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.5, cex.axis=1.8, xlim=c(20, 37), ylim=c(0, 0.10), xlab="Degree", ylab="Frequency")
lines(dd.Root.Ctr.Root.ig, lwd = 4, col = "gray50")
lines(dd.Root.As.Root.ig, lwd = 4, col = "#E69F00")
lines(dd.Root.Bx.Root.ig, lwd = 4, col = "#009E73")
lines(dd.Root.Tb.Root.ig, lwd = 4, col = "#D55E00")
dev.off()

#### Sediment Degree distributions 
CMRE9_Sedimentb_Ctrc_Sediment_ig= readRDS("CMRE9_Sedimentb_Ctrc_Sediment_ig.RDS")
CMRE9_Sedimentb_Asc_Sediment_ig= readRDS("CMRE9_Sedimentb_Asc_Sediment_ig.RDS")
CMRE9_Sedimentb_Bxc_Sediment_ig= readRDS("CMRE9_Sedimentb_Bxc_Sediment_ig.RDS")
CMRE9_Sedimentb_Tbc_Sediment_ig= readRDS("CMRE9_Sedimentb_Tbc_Sediment_ig.RDS")

dd.Sediment.Ctr.Sediment.ig <- degree.distribution(CMRE9_Sedimentb_Ctrc_Sediment_ig)
dd.Sediment.As.Sediment.ig <- degree.distribution(CMRE9_Sedimentb_Asc_Sediment_ig)
dd.Sediment.Bx.Sediment.ig <- degree.distribution(CMRE9_Sedimentb_Bxc_Sediment_ig)
dd.Sediment.Tb.Sediment.ig <- degree.distribution(CMRE9_Sedimentb_Tbc_Sediment_ig)

tiff(filename="Sediment_Degree_Ctr-As-Bx-Tb6.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,3.9,1,1))
plot(dd.Sediment.Ctr.Sediment.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.2, cex.axis=1.5, xlim=c(0, 72), ylim=c(0, 0.22), xlab="Degree", ylab="Frequency")
lines(dd.Sediment.Ctr.Sediment.ig, lwd = 2, col = "gray50")
lines(dd.Sediment.As.Sediment.ig, lwd = 2, col = "#E69F00")
lines(dd.Sediment.Bx.Sediment.ig, lwd = 2, col = "#009E73")
lines(dd.Sediment.Tb.Sediment.ig, lwd = 2, col = "#D55E00")
dev.off()

tiff(filename="Sediment_Degree_Ctr-As-Bx-Tb5_Zoom.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,4.1,1,1))
plot(dd.Sediment.Ctr.Sediment.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.5, cex.axis=1.8, xlim=c(35, 53), ylim=c(0, 0.10), xlab="Degree", ylab="Frequency")
lines(dd.Sediment.Ctr.Sediment.ig, lwd = 4, col = "gray50")
lines(dd.Sediment.As.Sediment.ig, lwd = 4, col = "#E69F00")
lines(dd.Sediment.Bx.Sediment.ig, lwd = 4, col = "#009E73")
lines(dd.Sediment.Tb.Sediment.ig, lwd = 4, col = "#D55E00")
dev.off()

#### Water Degree distributions 
CMRE9_Waterb_Ctrc_Water_ig= readRDS("CMRE9_Waterb_Ctrc_Water_ig.RDS")
CMRE9_Waterb_Asc_Water_ig= readRDS("CMRE9_Waterb_Asc_Water_ig.RDS")
CMRE9_Waterb_Bxc_Water_ig= readRDS("CMRE9_Waterb_Bxc_Water_ig.RDS")
CMRE9_Waterb_Tbc_Water_ig= readRDS("CMRE9_Waterb_Tbc_Water_ig.RDS")

dd.Water.Ctr.Water.ig <- degree.distribution(CMRE9_Waterb_Ctrc_Water_ig)
dd.Water.As.Water.ig <- degree.distribution(CMRE9_Waterb_Asc_Water_ig)
dd.Water.Bx.Water.ig <- degree.distribution(CMRE9_Waterb_Bxc_Water_ig)
dd.Water.Tb.Water.ig <- degree.distribution(CMRE9_Waterb_Tbc_Water_ig)

tiff(filename="Water_Degree_Ctr-As-Bx-Tb6.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,3.9,1,1))
plot(dd.Water.Ctr.Water.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.2, cex.axis=1.5, xlim=c(0, 72), ylim=c(0, 0.22), xlab="Degree", ylab="Frequency")
lines(dd.Water.Ctr.Water.ig, lwd = 2, col = "gray50")
lines(dd.Water.As.Water.ig, lwd = 2, col = "#E69F00")
lines(dd.Water.Bx.Water.ig, lwd = 2, col = "#009E73")
lines(dd.Water.Tb.Water.ig, lwd = 2, col = "#D55E00")
dev.off()

tiff(filename="Water_Degree_Ctr-As-Bx-Tb5_zoom.tif", width = 5, height = 4.5, units = 'in', res = 300, bg="white")
par(mar=c(3.6,4.1,1,1))
plot(dd.Water.Ctr.Water.ig, mgp = c(2.4, 0.8, 0), pch = 19, font=2, cex = 0.1, cex.lab=2.5, cex.axis=1.8, xlim=c(25, 38), ylim=c(0, 0.12), xlab="Degree", ylab="Frequency")
lines(dd.Water.Ctr.Water.ig, lwd = 4, col = "gray50")
lines(dd.Water.As.Water.ig, lwd = 4, col = "#E69F00")
lines(dd.Water.Bx.Water.ig, lwd = 4, col = "#009E73")
lines(dd.Water.Tb.Water.ig, lwd = 4, col = "#D55E00")
dev.off()

##############################################################################################

##############################################################################################################################





