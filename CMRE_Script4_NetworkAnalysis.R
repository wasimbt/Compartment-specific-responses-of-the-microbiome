#Copyright (c) 2023 University of Bern, Bern, Switzerland. 
#Author: Wasimuddin (wasim.bt@gmail.com; https://www.researchgate.net/profile/Wasim-Uddin)

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
# library(ggplot2); packageVersion("ggplot2")
# [1] ‘3.3.6’

########################################################################################################
########################################################################################################
########################################################################################################

######------------Mouse------------############

CMRE5_Mouse <- subset_samples(CMRE5, Sampletype == "feces")
CMRE5_Mousea = prune_taxa(taxa_sums(CMRE5_Mouse) > 0, CMRE5_Mouse) 
CMRE5_Mouseb <- prune_taxa(taxa_sums(CMRE5_Mousea) > 100, CMRE5_Mousea)
##################################################
######------------Mouse-Ctr------------############
CMRE5_Mouseb_Ctr <- subset_samples(CMRE5_Mouseb, Treatment == "Ctr" & Timepoint %in%c("1", "7"))
CMRE5_Mouseb_Ctra = prune_taxa(taxa_sums(CMRE5_Mouseb_Ctr) > 0, CMRE5_Mouseb_Ctr) #
##########Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(CMRE5_Mouseb_Ctra),
               MARGIN = ifelse(taxa_are_rows(CMRE5_Mouseb_Ctra), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(CMRE5_Mouseb_Ctra),
                    tax_table(CMRE5_Mouseb_Ctra))
prevalenceThreshold = 0.15 * nsamples(CMRE5_Mouseb_Ctra) #
Mouseb_Ctr_keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
CMRE5_Mouseb_Ctrb = prune_taxa(Mouseb_Ctr_keepTaxa, CMRE5_Mouseb_Ctra)

##########Add taxonomic classification to OTU ID
CMRE5_Mouseb_Ctrc <- format_to_besthit1(CMRE5_Mouseb_Ctrb)# run this function instead (see end of this script page)
CMRE5_Mouseb_Ctrc_otu <- t(otu_table(CMRE5_Mouseb_Ctrc)@.Data) #extract the otu table from phyloseq object
CMRE5_Mouseb_Ctrc_tax <- as.data.frame(tax_table(CMRE5_Mouseb_Ctrc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Mouseb_Ctrc_otu_net <- spiec.easi(CMRE5_Mouseb_Ctrc_otu, method='mb', lambda.min.ratio=1e-3, nlambda=50, pulsar.params=list(rep.num=99)) # increaseing permutations can decrease the number of edges (removes unstable edges)

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Ctrc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
Mouse.ig= readRDS("CMRE5_Mouseb_Ctrc_Mouse_ig.RDS")

nt_all_deg <- igraph::degree(Mouse.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Mouse.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Mouse.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Mouse.ig, "local", vids = V(Mouse.ig))
names(nt_all_transitivity)<- V(Mouse.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Mouse-As------------#############

CMRE5_Mouseb_As <- subset_samples(CMRE5_Mouseb, Treatment == "As" & Timepoint %in%c("1", "7"))
CMRE5_Mouseb_Asa = prune_taxa(taxa_sums(CMRE5_Mouseb_As) > 0, CMRE5_Mouseb_As) #remove these OTUs from class object

CMRE5_Mouseb_Asb = prune_taxa(Mouseb_Ctr_keepTaxa, CMRE5_Mouseb_Asa)
CMRE5_Mouseb_Asc <- format_to_besthit1(CMRE5_Mouseb_Asb)# run this function instead (see end of this script page)

CMRE5_Mouseb_Asc_otu <- t(otu_table(CMRE5_Mouseb_Asc)@.Data) #extract the otu table from phyloseq object
CMRE5_Mouseb_Asc_tax <- as.data.frame(tax_table(CMRE5_Mouseb_Asc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Mouseb_Asc_otu_net <- spiec.easi(CMRE5_Mouseb_Asc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Asc_otu_net)))

positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
Mouse.ig= readRDS("CMRE5_Mouseb_Asc_Mouse_ig.RDS")

nt_all_deg <- igraph::degree(Mouse.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Mouse.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Mouse.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Mouse.ig, "local", vids = V(Mouse.ig))
names(nt_all_transitivity)<- V(Mouse.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0
## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Mouse-Bx------------#############

CMRE5_Mouseb_Bx <- subset_samples(CMRE5_Mouseb, Treatment == "Bx" & Timepoint %in%c("1", "7"))
CMRE5_Mouseb_Bxa = prune_taxa(taxa_sums(CMRE5_Mouseb_Bx) > 0, CMRE5_Mouseb_Bx) #
CMRE5_Mouseb_Bxb = prune_taxa(Mouseb_Ctr_keepTaxa, CMRE5_Mouseb_Bxa)

CMRE5_Mouseb_Bxc <- format_to_besthit1(CMRE5_Mouseb_Bxb)# 
CMRE5_Mouseb_Bxc_otu <- t(otu_table(CMRE5_Mouseb_Bxc)@.Data) #
CMRE5_Mouseb_Bxc_tax <- as.data.frame(tax_table(CMRE5_Mouseb_Bxc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Mouseb_Bxc_otu_net <- spiec.easi(CMRE5_Mouseb_Bxc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) # 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Bxc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Mouse.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Mouse.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Mouse.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Mouse.ig, "local", vids = V(Mouse.ig))
names(nt_all_transitivity)<- V(Mouse.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Mouse-Tb------------#############

CMRE5_Mouseb_Tb <- subset_samples(CMRE5_Mouseb, Treatment == "Tb" & Timepoint %in%c("1", "7"))
CMRE5_Mouseb_Tba = prune_taxa(taxa_sums(CMRE5_Mouseb_Tb) > 0, CMRE5_Mouseb_Tb) #

CMRE5_Mouseb_Tbb = prune_taxa(Mouseb_Ctr_keepTaxa, CMRE5_Mouseb_Tba)
CMRE5_Mouseb_Tbc <- format_to_besthit1(CMRE5_Mouseb_Tbb)# 

CMRE5_Mouseb_Tbc_otu <- t(otu_table(CMRE5_Mouseb_Tbc)@.Data) #
CMRE5_Mouseb_Tbc_tax <- as.data.frame(tax_table(CMRE5_Mouseb_Tbc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Mouseb_Tbc_otu_net <- spiec.easi(CMRE5_Mouseb_Tbc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Mouseb_Tbc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Mouse.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Mouse.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Mouse.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Mouse.ig, "local", vids = V(Mouse.ig))
names(nt_all_transitivity)<- V(Mouse.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##########Perform Kolmogorov-Smirnov test

Mouse_Ctr <- read.csv("Network_properties_Mouse_Ctr.csv")
Mouse_As <- read.csv("Network_properties_Mouse_As.csv")
Mouse_Bx <- read.csv("Network_properties_Mouse_Bx.csv")
Mouse_Tb <- read.csv("Network_properties_Mouse_Tb.csv")

ks.test(Mouse_Ctr$nt_boot_degree, Mouse_As$nt_boot_degree)
ks.test(Mouse_Ctr$nt_boot_degree, Mouse_Bx$nt_boot_degree)
ks.test(Mouse_Ctr$nt_boot_degree, Mouse_Tb$nt_boot_degree)
ks.test(Mouse_As$nt_boot_degree, Mouse_Bx$nt_boot_degree)
ks.test(Mouse_As$nt_boot_degree, Mouse_Tb$nt_boot_degree)
ks.test(Mouse_Bx$nt_boot_degree, Mouse_Tb$nt_boot_degree)
ks.test(Mouse_Ctr$nt_boot_betweenness, Mouse_As$nt_boot_betweenness)
ks.test(Mouse_Ctr$nt_boot_betweenness, Mouse_Bx$nt_boot_betweenness)
ks.test(Mouse_Ctr$nt_boot_betweenness, Mouse_Tb$nt_boot_betweenness)
ks.test(Mouse_As$nt_boot_betweenness, Mouse_Bx$nt_boot_betweenness)
ks.test(Mouse_As$nt_boot_betweenness, Mouse_Tb$nt_boot_betweenness)
ks.test(Mouse_Bx$nt_boot_betweenness, Mouse_Tb$nt_boot_betweenness)
ks.test(Mouse_Ctr$nt_boot_closeness, Mouse_As$nt_boot_closeness)
ks.test(Mouse_Ctr$nt_boot_closeness, Mouse_Bx$nt_boot_closeness)
ks.test(Mouse_Ctr$nt_boot_closeness, Mouse_Tb$nt_boot_closeness)
ks.test(Mouse_As$nt_boot_closeness, Mouse_Bx$nt_boot_closeness)
ks.test(Mouse_As$nt_boot_closeness, Mouse_Tb$nt_boot_closeness)
ks.test(Mouse_Bx$nt_boot_closeness, Mouse_Tb$nt_boot_closeness)
ks.test(Mouse_Ctr$nt_boot_transitivity, Mouse_As$nt_boot_transitivity)
ks.test(Mouse_Ctr$nt_boot_transitivity, Mouse_Bx$nt_boot_transitivity)
ks.test(Mouse_Ctr$nt_boot_transitivity, Mouse_Tb$nt_boot_transitivity)
ks.test(Mouse_As$nt_boot_transitivity, Mouse_Bx$nt_boot_transitivity)
ks.test(Mouse_As$nt_boot_transitivity, Mouse_Tb$nt_boot_transitivity)
ks.test(Mouse_Bx$nt_boot_transitivity, Mouse_Tb$nt_boot_transitivity)

##########Hub Scores
Mouse_Ctr_net.hs <- hub_score(Mouse.Ctr.Mouse.ig)$vector
Mouse_Ctr_net.hs.sort <- sort(Mouse_Ctr_net.hs, decreasing = TRUE)

Mouse_As_net.hs <- hub_score(Mouse.As.Mouse.ig)$vector
Mouse_As_net.hs.sort <- sort(Mouse_As_net.hs, decreasing = TRUE)

Mouse_Bx_net.hs <- hub_score(Mouse.Bx.Mouse.ig)$vector
Mouse_Bx_net.hs.sort <- sort(Mouse_Bx_net.hs, decreasing = TRUE)

Mouse_Tb_net.hs <- hub_score(Mouse.Tb.Mouse.ig)$vector
Mouse_Tb_net.hs.sort <- sort(Mouse_Tb_net.hs, decreasing = TRUE)

##########power law testing Mouse-Ctr
Mouse.ig= readRDS("CMRE5_Mouseb_Ctrc_Mouse_ig.RDS")

Mouse.net <- asNetwork(Mouse.ig)
Mouse.ctr.spiec.deg <- degree(Mouse.net)

Mouse.ctr.spiec.deg_PL = displ$new(Mouse.ctr.spiec.deg)# power law distribution
Mouse.ctr.est = estimate_xmin(Mouse.ctr.spiec.deg_PL)
Mouse.ctr.spiec.deg_PL$xmin <- Mouse.ctr.est$xmin
Mouse.ctr.spiec.deg_PL$pars <- Mouse.ctr.est$pars

set.seed(101)
Mouse.ctr.bs <- bootstrap_p(Mouse.ctr.spiec.deg_PL, threads=4, no_of_sims = 1000)
Mouse.ctr.spiec.deg_N <- dislnorm$new(Mouse.ctr.spiec.deg)# Normal distribution
Mouse.ctr.spiec.deg_N$xmin <- Mouse.ctr.est$xmin
Mouse.ctr.spiec.deg_N$pars <- estimate_pars(Mouse.ctr.spiec.deg_N)
Mouse.ctr.bs <- bootstrap_p(Mouse.ctr.spiec.deg_N, threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.ctr.spiec.deg_PL, Mouse.ctr.spiec.deg_N)
Mouse.ctr.spiec.deg_E <- disexp$new(Mouse.ctr.spiec.deg)#  Exponential distribution
Mouse.ctr.spiec.deg_E$xmin <- Mouse.ctr.est$xmin
Mouse.ctr.spiec.deg_E$pars <- estimate_pars(Mouse.ctr.spiec.deg_E)

Mouse.ctr.bs <- bootstrap_p(Mouse.ctr.spiec.deg_E, threads=4, no_of_sims = 1000)
Mouse.comp <- compare_distributions(Mouse.ctr.spiec.deg_PL, Mouse.ctr.spiec.deg_E)
Mouse.ctr.spiec.deg_P <- dispois$new(Mouse.ctr.spiec.deg)# Poisson distribution
Mouse.ctr.spiec.deg_P$xmin <- Mouse.ctr.est$xmin
Mouse.ctr.spiec.deg_P$pars <- estimate_pars(Mouse.ctr.spiec.deg_P)
Mouse.ctr.bs <- bootstrap_p(Mouse.ctr.spiec.deg_P, threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.ctr.spiec.deg_PL, Mouse.ctr.spiec.deg_P)
Mouse.compNE <- compare_distributions(Mouse.ctr.spiec.deg_N, Mouse.ctr.spiec.deg_E)
Mouse.compEP <- compare_distributions(Mouse.ctr.spiec.deg_E, Mouse.ctr.spiec.deg_P)
Mouse.compNP <- compare_distributions(Mouse.ctr.spiec.deg_N, Mouse.ctr.spiec.deg_P)

##########power law testing Mouse-As
Mouse.ig= readRDS("CMRE5_Mouseb_Asc_Mouse_ig.RDS")

Mouse.net <- asNetwork(Mouse.ig)
Mouse.As.spiec.deg <- degree(Mouse.net)

Mouse.As.spiec.deg_PL = displ$new(Mouse.As.spiec.deg)# power law distribution
Mouse.As.est = estimate_xmin(Mouse.As.spiec.deg_PL)
Mouse.As.spiec.deg_PL$setXmin(Mouse.As.est)
Mouse.As.spiec.deg_PL$xmin <- Mouse.As.est$xmin
Mouse.As.spiec.deg_PL$pars <- Mouse.As.est$pars

set.seed(101)
bs <- bootstrap_p(Mouse.As.spiec.deg_PL, threads=4, no_of_sims = 1000)

Mouse.As.spiec.deg_N <- dislnorm$new(Mouse.As.spiec.deg)# Normal distribution
Mouse.As.spiec.deg_N$xmin <- Mouse.As.est$xmin
Mouse.As.spiec.deg_N$pars <- estimate_pars(Mouse.As.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Mouse.As.spiec.deg_N, threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.As.spiec.deg_PL, Mouse.As.spiec.deg_N)
Mouse.As.spiec.deg_E <- disexp$new(Mouse.As.spiec.deg)#  Exponential distribution
Mouse.As.spiec.deg_E$xmin <- Mouse.As.est$xmin
Mouse.As.spiec.deg_E$pars <- estimate_pars(Mouse.As.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Mouse.As.spiec.deg_E, threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.As.spiec.deg_PL, Mouse.As.spiec.deg_E)
Mouse.As.spiec.deg_P <- dispois$new(Mouse.As.spiec.deg)# Poisson distribution
Mouse.As.spiec.deg_P$xmin <- Mouse.As.est$xmin
Mouse.As.spiec.deg_P$pars <- estimate_pars(Mouse.As.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Mouse.As.spiec.deg_P, threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.As.spiec.deg_PL, Mouse.As.spiec.deg_P)
Mouse.compNE <- compare_distributions(Mouse.As.spiec.deg_N, Mouse.As.spiec.deg_E)
Mouse.compEP <- compare_distributions(Mouse.As.spiec.deg_E, Mouse.As.spiec.deg_P)
Mouse.compNP <- compare_distributions(Mouse.As.spiec.deg_N, Mouse.As.spiec.deg_P)

##########power law testing Mouse-Bx
Mouse.ig= readRDS("CMRE5_Mouseb_Bxc_Mouse_ig.RDS")
Mouse.net <- asNetwork(Mouse.ig)
Mouse.Bx.spiec.deg <- degree(Mouse.net)

Mouse.Bx.spiec.deg_PL = displ$new(Mouse.Bx.spiec.deg)# power law distribution
Mouse.Bx.est = estimate_xmin(Mouse.Bx.spiec.deg_PL)
Mouse.Bx.spiec.deg_PL$xmin <- Mouse.Bx.est$xmin
Mouse.Bx.spiec.deg_PL$pars <- Mouse.Bx.est$pars

set.seed(101)
bs <- bootstrap_p(Mouse.Bx.spiec.deg_PL,threads=4, no_of_sims = 1000)
Mouse.Bx.spiec.deg_N <- dislnorm$new(Mouse.Bx.spiec.deg)# Normal distribution
Mouse.Bx.spiec.deg_N$xmin <- Mouse.Bx.est$xmin
Mouse.Bx.spiec.deg_N$pars <- estimate_pars(Mouse.Bx.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Mouse.Bx.spiec.deg_N,threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.Bx.spiec.deg_PL, Mouse.Bx.spiec.deg_N)
Mouse.Bx.spiec.deg_E <- disexp$new(Mouse.Bx.spiec.deg)#  Exponential distribution
Mouse.Bx.spiec.deg_E$xmin <- Mouse.Bx.est$xmin
Mouse.Bx.spiec.deg_E$pars <- estimate_pars(Mouse.Bx.spiec.deg_E)

set.seed(101)
bs <- bootstrap_p(Mouse.Bx.spiec.deg_E,threads=4, no_of_sims = 1000)
Mouse.comp <- compare_distributions(Mouse.Bx.spiec.deg_PL, Mouse.Bx.spiec.deg_E)
Mouse.Bx.spiec.deg_P <- dispois$new(Mouse.Bx.spiec.deg)# Poisson distribution
Mouse.Bx.spiec.deg_P$xmin <- Mouse.Bx.est$xmin
Mouse.Bx.spiec.deg_P$pars <- estimate_pars(Mouse.Bx.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Mouse.Bx.spiec.deg_P,threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.Bx.spiec.deg_PL, Mouse.Bx.spiec.deg_P)
Mouse.compNE <- compare_distributions(Mouse.Bx.spiec.deg_N, Mouse.Bx.spiec.deg_E)
Mouse.compEP <- compare_distributions(Mouse.Bx.spiec.deg_E, Mouse.Bx.spiec.deg_P)
Mouse.compNP <- compare_distributions(Mouse.Bx.spiec.deg_N, Mouse.Bx.spiec.deg_P)

##########power law testing Mouse-Tb
Mouse.ig= readRDS("CMRE5_Mouseb_Tbc_Mouse_ig.RDS")

Mouse.net <- asNetwork(Mouse.ig)
Mouse.Tb.spiec.deg <- degree(Mouse.net)

Mouse.Tb.spiec.deg_PL = displ$new(Mouse.Tb.spiec.deg)# power law distribution
est = estimate_xmin(Mouse.Tb.spiec.deg_PL)
Mouse.Tb.spiec.deg_PL$xmin <- est$xmin
Mouse.Tb.spiec.deg_PL$pars <- est$pars

set.seed(101)
bs <- bootstrap_p(Mouse.Tb.spiec.deg_PL,threads=4, no_of_sims = 1000)
Mouse.Tb.spiec.deg_N <- dislnorm$new(Mouse.Tb.spiec.deg)# Normal distribution
Mouse.Tb.spiec.deg_N$xmin <- est$xmin
Mouse.Tb.spiec.deg_N$pars <- estimate_pars(Mouse.Tb.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Mouse.Tb.spiec.deg_N,threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.Tb.spiec.deg_PL, Mouse.Tb.spiec.deg_N)
Mouse.Tb.spiec.deg_E <- disexp$new(Mouse.Tb.spiec.deg)#  Exponential distribution
Mouse.Tb.spiec.deg_E$xmin <- est$xmin
Mouse.Tb.spiec.deg_E$pars <- estimate_pars(Mouse.Tb.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Mouse.Tb.spiec.deg_E,threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.Tb.spiec.deg_PL, Mouse.Tb.spiec.deg_E)
Mouse.Tb.spiec.deg_P <- dispois$new(Mouse.Tb.spiec.deg)# Poisson distribution
Mouse.Tb.spiec.deg_P$xmin <- est$xmin
Mouse.Tb.spiec.deg_P$pars <- estimate_pars(Mouse.Tb.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Mouse.Tb.spiec.deg_P,threads=4, no_of_sims = 1000)

Mouse.comp <- compare_distributions(Mouse.Tb.spiec.deg_PL, Mouse.Tb.spiec.deg_P)
Mouse.compNE <- compare_distributions(Mouse.Tb.spiec.deg_N, Mouse.Tb.spiec.deg_E)
Mouse.compEP <- compare_distributions(Mouse.Tb.spiec.deg_E, Mouse.Tb.spiec.deg_P)
Mouse.compNP <- compare_distributions(Mouse.Tb.spiec.deg_N, Mouse.Tb.spiec.deg_P)


########################################################################################################
########################################################################################################
########################################################################################################

######------------Soil------------############

CMRE5_Soil <- subset_samples(CMRE5, Sampletype == "soil")
CMRE5_Soila = prune_taxa(taxa_sums(CMRE5_Soil) > 0, CMRE5_Soil) 
CMRE5_Soilb <- prune_taxa(taxa_sums(CMRE5_Soila) > 100, CMRE5_Soila)
##################################################
######------------Soil-Ctr------------############
CMRE5_Soilb_Ctr <- subset_samples(CMRE5_Soilb, Treatment == "Ctr" & Timepoint %in%c("1", "7"))
CMRE5_Soilb_Ctra = prune_taxa(taxa_sums(CMRE5_Soilb_Ctr) > 0, CMRE5_Soilb_Ctr) #
##########Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(CMRE5_Soilb_Ctra),
               MARGIN = ifelse(taxa_are_rows(CMRE5_Soilb_Ctra), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(CMRE5_Soilb_Ctra),
                    tax_table(CMRE5_Soilb_Ctra))
prevalenceThreshold = 0.15 * nsamples(CMRE5_Soilb_Ctra) #
Soilb_Ctr_keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
CMRE5_Soilb_Ctrb = prune_taxa(Soilb_Ctr_keepTaxa, CMRE5_Soilb_Ctra)

##########Add taxonomic classification to OTU ID
CMRE5_Soilb_Ctrc <- format_to_besthit1(CMRE5_Soilb_Ctrb)# run this function instead (see end of this script page)
CMRE5_Soilb_Ctrc_otu <- t(otu_table(CMRE5_Soilb_Ctrc)@.Data) #extract the otu table from phyloseq object
CMRE5_Soilb_Ctrc_tax <- as.data.frame(tax_table(CMRE5_Soilb_Ctrc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Soilb_Ctrc_otu_net <- spiec.easi(CMRE5_Soilb_Ctrc_otu, method='mb', lambda.min.ratio=1e-3, nlambda=50, pulsar.params=list(rep.num=99)) # increaseing permutations can decrease the number of edges (removes unstable edges)

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Ctrc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
Soil.ig= readRDS("CMRE5_Soilb_Ctrc_Soil_ig.RDS")

nt_all_deg <- igraph::degree(Soil.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Soil.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Soil.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Soil.ig, "local", vids = V(Soil.ig))
names(nt_all_transitivity)<- V(Soil.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Soil-As------------#############

CMRE5_Soilb_As <- subset_samples(CMRE5_Soilb, Treatment == "As" & Timepoint %in%c("1", "7"))
CMRE5_Soilb_Asa = prune_taxa(taxa_sums(CMRE5_Soilb_As) > 0, CMRE5_Soilb_As) #remove these OTUs from class object

CMRE5_Soilb_Asb = prune_taxa(Soilb_Ctr_keepTaxa, CMRE5_Soilb_Asa)
CMRE5_Soilb_Asc <- format_to_besthit1(CMRE5_Soilb_Asb)# run this function instead (see end of this script page)

CMRE5_Soilb_Asc_otu <- t(otu_table(CMRE5_Soilb_Asc)@.Data) #extract the otu table from phyloseq object
CMRE5_Soilb_Asc_tax <- as.data.frame(tax_table(CMRE5_Soilb_Asc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Soilb_Asc_otu_net <- spiec.easi(CMRE5_Soilb_Asc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Asc_otu_net)))

positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
Soil.ig= readRDS("CMRE5_Soilb_Asc_Soil_ig.RDS")

nt_all_deg <- igraph::degree(Soil.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Soil.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Soil.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Soil.ig, "local", vids = V(Soil.ig))
names(nt_all_transitivity)<- V(Soil.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0
## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Soil-Bx------------#############

CMRE5_Soilb_Bx <- subset_samples(CMRE5_Soilb, Treatment == "Bx" & Timepoint %in%c("1", "7"))
CMRE5_Soilb_Bxa = prune_taxa(taxa_sums(CMRE5_Soilb_Bx) > 0, CMRE5_Soilb_Bx) #
CMRE5_Soilb_Bxb = prune_taxa(Soilb_Ctr_keepTaxa, CMRE5_Soilb_Bxa)

CMRE5_Soilb_Bxc <- format_to_besthit1(CMRE5_Soilb_Bxb)# 
CMRE5_Soilb_Bxc_otu <- t(otu_table(CMRE5_Soilb_Bxc)@.Data) #
CMRE5_Soilb_Bxc_tax <- as.data.frame(tax_table(CMRE5_Soilb_Bxc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Soilb_Bxc_otu_net <- spiec.easi(CMRE5_Soilb_Bxc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) # 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Bxc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Soil.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Soil.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Soil.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Soil.ig, "local", vids = V(Soil.ig))
names(nt_all_transitivity)<- V(Soil.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Soil-Tb------------#############

CMRE5_Soilb_Tb <- subset_samples(CMRE5_Soilb, Treatment == "Tb" & Timepoint %in%c("1", "7"))
CMRE5_Soilb_Tba = prune_taxa(taxa_sums(CMRE5_Soilb_Tb) > 0, CMRE5_Soilb_Tb) #

CMRE5_Soilb_Tbb = prune_taxa(Soilb_Ctr_keepTaxa, CMRE5_Soilb_Tba)
CMRE5_Soilb_Tbc <- format_to_besthit1(CMRE5_Soilb_Tbb)# 

CMRE5_Soilb_Tbc_otu <- t(otu_table(CMRE5_Soilb_Tbc)@.Data) #
CMRE5_Soilb_Tbc_tax <- as.data.frame(tax_table(CMRE5_Soilb_Tbc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Soilb_Tbc_otu_net <- spiec.easi(CMRE5_Soilb_Tbc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Soilb_Tbc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Soil.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Soil.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Soil.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Soil.ig, "local", vids = V(Soil.ig))
names(nt_all_transitivity)<- V(Soil.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##########Perform Kolmogorov-Smirnov test

Soil_Ctr <- read.csv("Network_properties_Soil_Ctr.csv")
Soil_As <- read.csv("Network_properties_Soil_As.csv")
Soil_Bx <- read.csv("Network_properties_Soil_Bx.csv")
Soil_Tb <- read.csv("Network_properties_Soil_Tb.csv")

ks.test(Soil_Ctr$nt_boot_degree, Soil_As$nt_boot_degree)
ks.test(Soil_Ctr$nt_boot_degree, Soil_Bx$nt_boot_degree)
ks.test(Soil_Ctr$nt_boot_degree, Soil_Tb$nt_boot_degree)
ks.test(Soil_As$nt_boot_degree, Soil_Bx$nt_boot_degree)
ks.test(Soil_As$nt_boot_degree, Soil_Tb$nt_boot_degree)
ks.test(Soil_Bx$nt_boot_degree, Soil_Tb$nt_boot_degree)
ks.test(Soil_Ctr$nt_boot_betweenness, Soil_As$nt_boot_betweenness)
ks.test(Soil_Ctr$nt_boot_betweenness, Soil_Bx$nt_boot_betweenness)
ks.test(Soil_Ctr$nt_boot_betweenness, Soil_Tb$nt_boot_betweenness)
ks.test(Soil_As$nt_boot_betweenness, Soil_Bx$nt_boot_betweenness)
ks.test(Soil_As$nt_boot_betweenness, Soil_Tb$nt_boot_betweenness)
ks.test(Soil_Bx$nt_boot_betweenness, Soil_Tb$nt_boot_betweenness)
ks.test(Soil_Ctr$nt_boot_closeness, Soil_As$nt_boot_closeness)
ks.test(Soil_Ctr$nt_boot_closeness, Soil_Bx$nt_boot_closeness)
ks.test(Soil_Ctr$nt_boot_closeness, Soil_Tb$nt_boot_closeness)
ks.test(Soil_As$nt_boot_closeness, Soil_Bx$nt_boot_closeness)
ks.test(Soil_As$nt_boot_closeness, Soil_Tb$nt_boot_closeness)
ks.test(Soil_Bx$nt_boot_closeness, Soil_Tb$nt_boot_closeness)
ks.test(Soil_Ctr$nt_boot_transitivity, Soil_As$nt_boot_transitivity)
ks.test(Soil_Ctr$nt_boot_transitivity, Soil_Bx$nt_boot_transitivity)
ks.test(Soil_Ctr$nt_boot_transitivity, Soil_Tb$nt_boot_transitivity)
ks.test(Soil_As$nt_boot_transitivity, Soil_Bx$nt_boot_transitivity)
ks.test(Soil_As$nt_boot_transitivity, Soil_Tb$nt_boot_transitivity)
ks.test(Soil_Bx$nt_boot_transitivity, Soil_Tb$nt_boot_transitivity)

##########Hub Scores
Soil_Ctr_net.hs <- hub_score(Soil.Ctr.Soil.ig)$vector
Soil_Ctr_net.hs.sort <- sort(Soil_Ctr_net.hs, decreasing = TRUE)

Soil_As_net.hs <- hub_score(Soil.As.Soil.ig)$vector
Soil_As_net.hs.sort <- sort(Soil_As_net.hs, decreasing = TRUE)

Soil_Bx_net.hs <- hub_score(Soil.Bx.Soil.ig)$vector
Soil_Bx_net.hs.sort <- sort(Soil_Bx_net.hs, decreasing = TRUE)

Soil_Tb_net.hs <- hub_score(Soil.Tb.Soil.ig)$vector
Soil_Tb_net.hs.sort <- sort(Soil_Tb_net.hs, decreasing = TRUE)

##########power law testing Soil-Ctr
Soil.ig= readRDS("CMRE5_Soilb_Ctrc_Soil_ig.RDS")

Soil.net <- asNetwork(Soil.ig)
Soil.ctr.spiec.deg <- degree(Soil.net)

Soil.ctr.spiec.deg_PL = displ$new(Soil.ctr.spiec.deg)# power law distribution
Soil.ctr.est = estimate_xmin(Soil.ctr.spiec.deg_PL)
Soil.ctr.spiec.deg_PL$xmin <- Soil.ctr.est$xmin
Soil.ctr.spiec.deg_PL$pars <- Soil.ctr.est$pars

set.seed(101)
Soil.ctr.bs <- bootstrap_p(Soil.ctr.spiec.deg_PL, threads=4, no_of_sims = 1000)
Soil.ctr.spiec.deg_N <- dislnorm$new(Soil.ctr.spiec.deg)# Normal distribution
Soil.ctr.spiec.deg_N$xmin <- Soil.ctr.est$xmin
Soil.ctr.spiec.deg_N$pars <- estimate_pars(Soil.ctr.spiec.deg_N)
Soil.ctr.bs <- bootstrap_p(Soil.ctr.spiec.deg_N, threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.ctr.spiec.deg_PL, Soil.ctr.spiec.deg_N)
Soil.ctr.spiec.deg_E <- disexp$new(Soil.ctr.spiec.deg)#  Exponential distribution
Soil.ctr.spiec.deg_E$xmin <- Soil.ctr.est$xmin
Soil.ctr.spiec.deg_E$pars <- estimate_pars(Soil.ctr.spiec.deg_E)

Soil.ctr.bs <- bootstrap_p(Soil.ctr.spiec.deg_E, threads=4, no_of_sims = 1000)
Soil.comp <- compare_distributions(Soil.ctr.spiec.deg_PL, Soil.ctr.spiec.deg_E)
Soil.ctr.spiec.deg_P <- dispois$new(Soil.ctr.spiec.deg)# Poisson distribution

Soil.ctr.spiec.deg_P$xmin <- Soil.ctr.est$xmin
Soil.ctr.spiec.deg_P$pars <- estimate_pars(Soil.ctr.spiec.deg_P)
Soil.ctr.bs <- bootstrap_p(Soil.ctr.spiec.deg_P, threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.ctr.spiec.deg_PL, Soil.ctr.spiec.deg_P)
Soil.compNE <- compare_distributions(Soil.ctr.spiec.deg_N, Soil.ctr.spiec.deg_E)
Soil.compEP <- compare_distributions(Soil.ctr.spiec.deg_E, Soil.ctr.spiec.deg_P)
Soil.compNP <- compare_distributions(Soil.ctr.spiec.deg_N, Soil.ctr.spiec.deg_P)


##########power law testing Soil-As
Soil.ig= readRDS("CMRE5_Soilb_Asc_Soil_ig.RDS")

Soil.net <- asNetwork(Soil.ig)
Soil.As.spiec.deg <- degree(Soil.net)

Soil.As.spiec.deg_PL = displ$new(Soil.As.spiec.deg)# power law distribution
Soil.As.est = estimate_xmin(Soil.As.spiec.deg_PL)
Soil.As.spiec.deg_PL$setXmin(Soil.As.est)
Soil.As.spiec.deg_PL$xmin <- Soil.As.est$xmin
Soil.As.spiec.deg_PL$pars <- Soil.As.est$pars

set.seed(101)
bs <- bootstrap_p(Soil.As.spiec.deg_PL, threads=4, no_of_sims = 1000)

Soil.As.spiec.deg_N <- dislnorm$new(Soil.As.spiec.deg)# Normal distribution
Soil.As.spiec.deg_N$xmin <- Soil.As.est$xmin
Soil.As.spiec.deg_N$pars <- estimate_pars(Soil.As.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Soil.As.spiec.deg_N, threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.As.spiec.deg_PL, Soil.As.spiec.deg_N)
Soil.As.spiec.deg_E <- disexp$new(Soil.As.spiec.deg)#  Exponential distribution
Soil.As.spiec.deg_E$xmin <- Soil.As.est$xmin
Soil.As.spiec.deg_E$pars <- estimate_pars(Soil.As.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Soil.As.spiec.deg_E, threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.As.spiec.deg_PL, Soil.As.spiec.deg_E)
Soil.As.spiec.deg_P <- dispois$new(Soil.As.spiec.deg)# Poisson distribution
Soil.As.spiec.deg_P$xmin <- Soil.As.est$xmin
Soil.As.spiec.deg_P$pars <- estimate_pars(Soil.As.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Soil.As.spiec.deg_P, threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.As.spiec.deg_PL, Soil.As.spiec.deg_P)
Soil.compNE <- compare_distributions(Soil.As.spiec.deg_N, Soil.As.spiec.deg_E)
Soil.compEP <- compare_distributions(Soil.As.spiec.deg_E, Soil.As.spiec.deg_P)
Soil.compNP <- compare_distributions(Soil.As.spiec.deg_N, Soil.As.spiec.deg_P)

##########power law testing Soil-Bx
Soil.ig= readRDS("CMRE5_Soilb_Bxc_Soil_ig.RDS")
Soil.net <- asNetwork(Soil.ig)
Soil.Bx.spiec.deg <- degree(Soil.net)

Soil.Bx.spiec.deg_PL = displ$new(Soil.Bx.spiec.deg)# power law distribution
Soil.Bx.est = estimate_xmin(Soil.Bx.spiec.deg_PL)
Soil.Bx.spiec.deg_PL$xmin <- Soil.Bx.est$xmin
Soil.Bx.spiec.deg_PL$pars <- Soil.Bx.est$pars

set.seed(101)
bs <- bootstrap_p(Soil.Bx.spiec.deg_PL,threads=4, no_of_sims = 1000)
Soil.Bx.spiec.deg_N <- dislnorm$new(Soil.Bx.spiec.deg)# Normal distribution
Soil.Bx.spiec.deg_N$xmin <- Soil.Bx.est$xmin
Soil.Bx.spiec.deg_N$pars <- estimate_pars(Soil.Bx.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Soil.Bx.spiec.deg_N,threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.Bx.spiec.deg_PL, Soil.Bx.spiec.deg_N)
Soil.Bx.spiec.deg_E <- disexp$new(Soil.Bx.spiec.deg)#  Exponential distribution
Soil.Bx.spiec.deg_E$xmin <- Soil.Bx.est$xmin
Soil.Bx.spiec.deg_E$pars <- estimate_pars(Soil.Bx.spiec.deg_E)

set.seed(101)
bs <- bootstrap_p(Soil.Bx.spiec.deg_E,threads=4, no_of_sims = 1000)
Soil.comp <- compare_distributions(Soil.Bx.spiec.deg_PL, Soil.Bx.spiec.deg_E)
Soil.Bx.spiec.deg_P <- dispois$new(Soil.Bx.spiec.deg)# Poisson distribution
Soil.Bx.spiec.deg_P$xmin <- Soil.Bx.est$xmin
Soil.Bx.spiec.deg_P$pars <- estimate_pars(Soil.Bx.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Soil.Bx.spiec.deg_P,threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.Bx.spiec.deg_PL, Soil.Bx.spiec.deg_P)
Soil.compNE <- compare_distributions(Soil.Bx.spiec.deg_N, Soil.Bx.spiec.deg_E)
Soil.compEP <- compare_distributions(Soil.Bx.spiec.deg_E, Soil.Bx.spiec.deg_P)
Soil.compNP <- compare_distributions(Soil.Bx.spiec.deg_N, Soil.Bx.spiec.deg_P)

##########power law testing Soil-Tb
Soil.ig= readRDS("CMRE5_Soilb_Tbc_Soil_ig.RDS")

Soil.net <- asNetwork(Soil.ig)
Soil.Tb.spiec.deg <- degree(Soil.net)

Soil.Tb.spiec.deg_PL = displ$new(Soil.Tb.spiec.deg)# power law distribution
est = estimate_xmin(Soil.Tb.spiec.deg_PL)
Soil.Tb.spiec.deg_PL$xmin <- est$xmin
Soil.Tb.spiec.deg_PL$pars <- est$pars

set.seed(101)
bs <- bootstrap_p(Soil.Tb.spiec.deg_PL,threads=4, no_of_sims = 1000)
Soil.Tb.spiec.deg_N <- dislnorm$new(Soil.Tb.spiec.deg)# Normal distribution
Soil.Tb.spiec.deg_N$xmin <- est$xmin
Soil.Tb.spiec.deg_N$pars <- estimate_pars(Soil.Tb.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Soil.Tb.spiec.deg_N,threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.Tb.spiec.deg_PL, Soil.Tb.spiec.deg_N)
Soil.Tb.spiec.deg_E <- disexp$new(Soil.Tb.spiec.deg)#  Exponential distribution
Soil.Tb.spiec.deg_E$xmin <- est$xmin
Soil.Tb.spiec.deg_E$pars <- estimate_pars(Soil.Tb.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Soil.Tb.spiec.deg_E,threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.Tb.spiec.deg_PL, Soil.Tb.spiec.deg_E)
Soil.Tb.spiec.deg_P <- dispois$new(Soil.Tb.spiec.deg)# Poisson distribution
Soil.Tb.spiec.deg_P$xmin <- est$xmin
Soil.Tb.spiec.deg_P$pars <- estimate_pars(Soil.Tb.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Soil.Tb.spiec.deg_P,threads=4, no_of_sims = 1000)

Soil.comp <- compare_distributions(Soil.Tb.spiec.deg_PL, Soil.Tb.spiec.deg_P)
Soil.compNE <- compare_distributions(Soil.Tb.spiec.deg_N, Soil.Tb.spiec.deg_E)
Soil.compEP <- compare_distributions(Soil.Tb.spiec.deg_E, Soil.Tb.spiec.deg_P)
Soil.compNP <- compare_distributions(Soil.Tb.spiec.deg_N, Soil.Tb.spiec.deg_P)

########################################################################################################
########################################################################################################
########################################################################################################

######------------Root------------############

CMRE5_Root <- subset_samples(CMRE5, Sampletype == "roots")
CMRE5_Roota = prune_taxa(taxa_sums(CMRE5_Root) > 0, CMRE5_Root) 
CMRE5_Rootb <- prune_taxa(taxa_sums(CMRE5_Roota) > 100, CMRE5_Roota)
##################################################
######------------Root-Ctr------------############
CMRE5_Rootb_Ctr <- subset_samples(CMRE5_Rootb, Treatment == "Ctr" & Timepoint %in%c("1", "7"))
CMRE5_Rootb_Ctra = prune_taxa(taxa_sums(CMRE5_Rootb_Ctr) > 0, CMRE5_Rootb_Ctr) #
##########Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(CMRE5_Rootb_Ctra),
               MARGIN = ifelse(taxa_are_rows(CMRE5_Rootb_Ctra), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(CMRE5_Rootb_Ctra),
                    tax_table(CMRE5_Rootb_Ctra))
prevalenceThreshold = 0.15 * nsamples(CMRE5_Rootb_Ctra) #
Rootb_Ctr_keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
CMRE5_Rootb_Ctrb = prune_taxa(Rootb_Ctr_keepTaxa, CMRE5_Rootb_Ctra)

##########Add taxonomic classification to OTU ID
CMRE5_Rootb_Ctrc <- format_to_besthit1(CMRE5_Rootb_Ctrb)# run this function instead (see end of this script page)
CMRE5_Rootb_Ctrc_otu <- t(otu_table(CMRE5_Rootb_Ctrc)@.Data) #extract the otu table from phyloseq object
CMRE5_Rootb_Ctrc_tax <- as.data.frame(tax_table(CMRE5_Rootb_Ctrc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Rootb_Ctrc_otu_net <- spiec.easi(CMRE5_Rootb_Ctrc_otu, method='mb', lambda.min.ratio=1e-3, nlambda=50, pulsar.params=list(rep.num=99)) # increaseing permutations can decrease the number of edges (removes unstable edges)

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Ctrc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
Root.ig= readRDS("CMRE5_Rootb_Ctrc_Root_ig.RDS")

nt_all_deg <- igraph::degree(Root.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Root.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Root.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Root.ig, "local", vids = V(Root.ig))
names(nt_all_transitivity)<- V(Root.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Root-As------------#############

CMRE5_Rootb_As <- subset_samples(CMRE5_Rootb, Treatment == "As" & Timepoint %in%c("1", "7"))
CMRE5_Rootb_Asa = prune_taxa(taxa_sums(CMRE5_Rootb_As) > 0, CMRE5_Rootb_As) #remove these OTUs from class object

CMRE5_Rootb_Asb = prune_taxa(Rootb_Ctr_keepTaxa, CMRE5_Rootb_Asa)
CMRE5_Rootb_Asc <- format_to_besthit1(CMRE5_Rootb_Asb)# run this function instead (see end of this script page)

CMRE5_Rootb_Asc_otu <- t(otu_table(CMRE5_Rootb_Asc)@.Data) #extract the otu table from phyloseq object
CMRE5_Rootb_Asc_tax <- as.data.frame(tax_table(CMRE5_Rootb_Asc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Rootb_Asc_otu_net <- spiec.easi(CMRE5_Rootb_Asc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Asc_otu_net)))

positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters

Root.ig= readRDS("CMRE5_Rootb_Asc_Root_ig.RDS")

nt_all_deg <- igraph::degree(Root.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Root.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Root.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Root.ig, "local", vids = V(Root.ig))
names(nt_all_transitivity)<- V(Root.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0
## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Root-Bx------------#############

CMRE5_Rootb_Bx <- subset_samples(CMRE5_Rootb, Treatment == "Bx" & Timepoint %in%c("1", "7"))
CMRE5_Rootb_Bxa = prune_taxa(taxa_sums(CMRE5_Rootb_Bx) > 0, CMRE5_Rootb_Bx) #
CMRE5_Rootb_Bxb = prune_taxa(Rootb_Ctr_keepTaxa, CMRE5_Rootb_Bxa)

CMRE5_Rootb_Bxc <- format_to_besthit1(CMRE5_Rootb_Bxb)# 
CMRE5_Rootb_Bxc_otu <- t(otu_table(CMRE5_Rootb_Bxc)@.Data) #
CMRE5_Rootb_Bxc_tax <- as.data.frame(tax_table(CMRE5_Rootb_Bxc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Rootb_Bxc_otu_net <- spiec.easi(CMRE5_Rootb_Bxc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) # 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Bxc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Root.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Root.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Root.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Root.ig, "local", vids = V(Root.ig))
names(nt_all_transitivity)<- V(Root.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Root-Tb------------#############

CMRE5_Rootb_Tb <- subset_samples(CMRE5_Rootb, Treatment == "Tb" & Timepoint %in%c("1", "7"))
CMRE5_Rootb_Tba = prune_taxa(taxa_sums(CMRE5_Rootb_Tb) > 0, CMRE5_Rootb_Tb) #

CMRE5_Rootb_Tbb = prune_taxa(Rootb_Ctr_keepTaxa, CMRE5_Rootb_Tba)
CMRE5_Rootb_Tbc <- format_to_besthit1(CMRE5_Rootb_Tbb)# 

CMRE5_Rootb_Tbc_otu <- t(otu_table(CMRE5_Rootb_Tbc)@.Data) #
CMRE5_Rootb_Tbc_tax <- as.data.frame(tax_table(CMRE5_Rootb_Tbc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Rootb_Tbc_otu_net <- spiec.easi(CMRE5_Rootb_Tbc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Rootb_Tbc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Root.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Root.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Root.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Root.ig, "local", vids = V(Root.ig))
names(nt_all_transitivity)<- V(Root.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##########Perform Kolmogorov-Smirnov test

Root_Ctr <- read.csv("Network_properties_Root_Ctr.csv")
Root_As <- read.csv("Network_properties_Root_As.csv")
Root_Bx <- read.csv("Network_properties_Root_Bx.csv")
Root_Tb <- read.csv("Network_properties_Root_Tb.csv")

ks.test(Root_Ctr$nt_boot_degree, Root_As$nt_boot_degree)
ks.test(Root_Ctr$nt_boot_degree, Root_Bx$nt_boot_degree)
ks.test(Root_Ctr$nt_boot_degree, Root_Tb$nt_boot_degree)
ks.test(Root_As$nt_boot_degree, Root_Bx$nt_boot_degree)
ks.test(Root_As$nt_boot_degree, Root_Tb$nt_boot_degree)
ks.test(Root_Bx$nt_boot_degree, Root_Tb$nt_boot_degree)
ks.test(Root_Ctr$nt_boot_betweenness, Root_As$nt_boot_betweenness)
ks.test(Root_Ctr$nt_boot_betweenness, Root_Bx$nt_boot_betweenness)
ks.test(Root_Ctr$nt_boot_betweenness, Root_Tb$nt_boot_betweenness)
ks.test(Root_As$nt_boot_betweenness, Root_Bx$nt_boot_betweenness)
ks.test(Root_As$nt_boot_betweenness, Root_Tb$nt_boot_betweenness)
ks.test(Root_Bx$nt_boot_betweenness, Root_Tb$nt_boot_betweenness)
ks.test(Root_Ctr$nt_boot_closeness, Root_As$nt_boot_closeness)
ks.test(Root_Ctr$nt_boot_closeness, Root_Bx$nt_boot_closeness)
ks.test(Root_Ctr$nt_boot_closeness, Root_Tb$nt_boot_closeness)
ks.test(Root_As$nt_boot_closeness, Root_Bx$nt_boot_closeness)
ks.test(Root_As$nt_boot_closeness, Root_Tb$nt_boot_closeness)
ks.test(Root_Bx$nt_boot_closeness, Root_Tb$nt_boot_closeness)
ks.test(Root_Ctr$nt_boot_transitivity, Root_As$nt_boot_transitivity)
ks.test(Root_Ctr$nt_boot_transitivity, Root_Bx$nt_boot_transitivity)
ks.test(Root_Ctr$nt_boot_transitivity, Root_Tb$nt_boot_transitivity)
ks.test(Root_As$nt_boot_transitivity, Root_Bx$nt_boot_transitivity)
ks.test(Root_As$nt_boot_transitivity, Root_Tb$nt_boot_transitivity)
ks.test(Root_Bx$nt_boot_transitivity, Root_Tb$nt_boot_transitivity)

##########Hub Scores
Root_Ctr_net.hs <- hub_score(Root.Ctr.Root.ig)$vector
Root_Ctr_net.hs.sort <- sort(Root_Ctr_net.hs, decreasing = TRUE)

Root_As_net.hs <- hub_score(Root.As.Root.ig)$vector
Root_As_net.hs.sort <- sort(Root_As_net.hs, decreasing = TRUE)

Root_Bx_net.hs <- hub_score(Root.Bx.Root.ig)$vector
Root_Bx_net.hs.sort <- sort(Root_Bx_net.hs, decreasing = TRUE)

Root_Tb_net.hs <- hub_score(Root.Tb.Root.ig)$vector
Root_Tb_net.hs.sort <- sort(Root_Tb_net.hs, decreasing = TRUE)

##########power law testing Root-Ctr
Root.ig= readRDS("CMRE5_Rootb_Ctrc_Root_ig.RDS")

Root.net <- asNetwork(Root.ig)
Root.ctr.spiec.deg <- degree(Root.net)

Root.ctr.spiec.deg_PL = displ$new(Root.ctr.spiec.deg)# power law distribution
Root.ctr.est = estimate_xmin(Root.ctr.spiec.deg_PL)
Root.ctr.spiec.deg_PL$xmin <- Root.ctr.est$xmin
Root.ctr.spiec.deg_PL$pars <- Root.ctr.est$pars

set.seed(101)
Root.ctr.bs <- bootstrap_p(Root.ctr.spiec.deg_PL, threads=4, no_of_sims = 1000)
Root.ctr.spiec.deg_N <- dislnorm$new(Root.ctr.spiec.deg)# Normal distribution
Root.ctr.spiec.deg_N$xmin <- Root.ctr.est$xmin
Root.ctr.spiec.deg_N$pars <- estimate_pars(Root.ctr.spiec.deg_N)
Root.ctr.bs <- bootstrap_p(Root.ctr.spiec.deg_N, threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.ctr.spiec.deg_PL, Root.ctr.spiec.deg_N)
Root.ctr.spiec.deg_E <- disexp$new(Root.ctr.spiec.deg)#  Exponential distribution
Root.ctr.spiec.deg_E$xmin <- Root.ctr.est$xmin
Root.ctr.spiec.deg_E$pars <- estimate_pars(Root.ctr.spiec.deg_E)

Root.ctr.bs <- bootstrap_p(Root.ctr.spiec.deg_E, threads=4, no_of_sims = 1000)
Root.comp <- compare_distributions(Root.ctr.spiec.deg_PL, Root.ctr.spiec.deg_E)
Root.ctr.spiec.deg_P <- dispois$new(Root.ctr.spiec.deg)# Poisson distribution
Root.ctr.spiec.deg_P$xmin <- Root.ctr.est$xmin
Root.ctr.spiec.deg_P$pars <- estimate_pars(Root.ctr.spiec.deg_P)
Root.ctr.bs <- bootstrap_p(Root.ctr.spiec.deg_P, threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.ctr.spiec.deg_PL, Root.ctr.spiec.deg_P)
Root.compNE <- compare_distributions(Root.ctr.spiec.deg_N, Root.ctr.spiec.deg_E)
Root.compEP <- compare_distributions(Root.ctr.spiec.deg_E, Root.ctr.spiec.deg_P)
Root.compNP <- compare_distributions(Root.ctr.spiec.deg_N, Root.ctr.spiec.deg_P)

##########power law testing Root-As
Root.ig= readRDS("CMRE5_Rootb_Asc_Root_ig.RDS")

Root.net <- asNetwork(Root.ig)
Root.As.spiec.deg <- degree(Root.net)

Root.As.spiec.deg_PL = displ$new(Root.As.spiec.deg)# power law distribution
Root.As.est = estimate_xmin(Root.As.spiec.deg_PL)
Root.As.spiec.deg_PL$setXmin(Root.As.est)
Root.As.spiec.deg_PL$xmin <- Root.As.est$xmin
Root.As.spiec.deg_PL$pars <- Root.As.est$pars

set.seed(101)
bs <- bootstrap_p(Root.As.spiec.deg_PL, threads=4, no_of_sims = 1000)

Root.As.spiec.deg_N <- dislnorm$new(Root.As.spiec.deg)# Normal distribution
Root.As.spiec.deg_N$xmin <- Root.As.est$xmin
Root.As.spiec.deg_N$pars <- estimate_pars(Root.As.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Root.As.spiec.deg_N, threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.As.spiec.deg_PL, Root.As.spiec.deg_N)
Root.As.spiec.deg_E <- disexp$new(Root.As.spiec.deg)#  Exponential distribution
Root.As.spiec.deg_E$xmin <- Root.As.est$xmin
Root.As.spiec.deg_E$pars <- estimate_pars(Root.As.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Root.As.spiec.deg_E, threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.As.spiec.deg_PL, Root.As.spiec.deg_E)
Root.As.spiec.deg_P <- dispois$new(Root.As.spiec.deg)# Poisson distribution
Root.As.spiec.deg_P$xmin <- Root.As.est$xmin
Root.As.spiec.deg_P$pars <- estimate_pars(Root.As.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Root.As.spiec.deg_P, threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.As.spiec.deg_PL, Root.As.spiec.deg_P)
Root.compNE <- compare_distributions(Root.As.spiec.deg_N, Root.As.spiec.deg_E)
Root.compEP <- compare_distributions(Root.As.spiec.deg_E, Root.As.spiec.deg_P)
Root.compNP <- compare_distributions(Root.As.spiec.deg_N, Root.As.spiec.deg_P)

##########power law testing Root-Bx
Root.ig= readRDS("CMRE5_Rootb_Bxc_Root_ig.RDS")
Root.net <- asNetwork(Root.ig)
Root.Bx.spiec.deg <- degree(Root.net)

Root.Bx.spiec.deg_PL = displ$new(Root.Bx.spiec.deg)# power law distribution
Root.Bx.est = estimate_xmin(Root.Bx.spiec.deg_PL)
Root.Bx.spiec.deg_PL$xmin <- Root.Bx.est$xmin
Root.Bx.spiec.deg_PL$pars <- Root.Bx.est$pars

set.seed(101)
bs <- bootstrap_p(Root.Bx.spiec.deg_PL,threads=4, no_of_sims = 1000)
Root.Bx.spiec.deg_N <- dislnorm$new(Root.Bx.spiec.deg)# Normal distribution
Root.Bx.spiec.deg_N$xmin <- Root.Bx.est$xmin
Root.Bx.spiec.deg_N$pars <- estimate_pars(Root.Bx.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Root.Bx.spiec.deg_N,threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.Bx.spiec.deg_PL, Root.Bx.spiec.deg_N)
Root.Bx.spiec.deg_E <- disexp$new(Root.Bx.spiec.deg)#  Exponential distribution
Root.Bx.spiec.deg_E$xmin <- Root.Bx.est$xmin
Root.Bx.spiec.deg_E$pars <- estimate_pars(Root.Bx.spiec.deg_E)

set.seed(101)
bs <- bootstrap_p(Root.Bx.spiec.deg_E,threads=4, no_of_sims = 1000)
Root.comp <- compare_distributions(Root.Bx.spiec.deg_PL, Root.Bx.spiec.deg_E)
Root.Bx.spiec.deg_P <- dispois$new(Root.Bx.spiec.deg)# Poisson distribution
Root.Bx.spiec.deg_P$xmin <- Root.Bx.est$xmin
Root.Bx.spiec.deg_P$pars <- estimate_pars(Root.Bx.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Root.Bx.spiec.deg_P,threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.Bx.spiec.deg_PL, Root.Bx.spiec.deg_P)
Root.compNE <- compare_distributions(Root.Bx.spiec.deg_N, Root.Bx.spiec.deg_E)
Root.compEP <- compare_distributions(Root.Bx.spiec.deg_E, Root.Bx.spiec.deg_P)
Root.compNP <- compare_distributions(Root.Bx.spiec.deg_N, Root.Bx.spiec.deg_P)

##########power law testing Root-Tb
Root.ig= readRDS("CMRE5_Rootb_Tbc_Root_ig.RDS")

Root.net <- asNetwork(Root.ig)
Root.Tb.spiec.deg <- degree(Root.net)

Root.Tb.spiec.deg_PL = displ$new(Root.Tb.spiec.deg)# power law distribution
est = estimate_xmin(Root.Tb.spiec.deg_PL)
Root.Tb.spiec.deg_PL$xmin <- est$xmin
Root.Tb.spiec.deg_PL$pars <- est$pars

set.seed(101)
bs <- bootstrap_p(Root.Tb.spiec.deg_PL,threads=4, no_of_sims = 1000)
Root.Tb.spiec.deg_N <- dislnorm$new(Root.Tb.spiec.deg)# Normal distribution
Root.Tb.spiec.deg_N$xmin <- est$xmin
Root.Tb.spiec.deg_N$pars <- estimate_pars(Root.Tb.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Root.Tb.spiec.deg_N,threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.Tb.spiec.deg_PL, Root.Tb.spiec.deg_N)
Root.Tb.spiec.deg_E <- disexp$new(Root.Tb.spiec.deg)#  Exponential distribution
Root.Tb.spiec.deg_E$xmin <- est$xmin
Root.Tb.spiec.deg_E$pars <- estimate_pars(Root.Tb.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Root.Tb.spiec.deg_E,threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.Tb.spiec.deg_PL, Root.Tb.spiec.deg_E)
Root.Tb.spiec.deg_P <- dispois$new(Root.Tb.spiec.deg)# Poisson distribution
Root.Tb.spiec.deg_P$xmin <- est$xmin
Root.Tb.spiec.deg_P$pars <- estimate_pars(Root.Tb.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Root.Tb.spiec.deg_P,threads=4, no_of_sims = 1000)

Root.comp <- compare_distributions(Root.Tb.spiec.deg_PL, Root.Tb.spiec.deg_P)
Root.compNE <- compare_distributions(Root.Tb.spiec.deg_N, Root.Tb.spiec.deg_E)
Root.compEP <- compare_distributions(Root.Tb.spiec.deg_E, Root.Tb.spiec.deg_P)
Root.compNP <- compare_distributions(Root.Tb.spiec.deg_N, Root.Tb.spiec.deg_P)

########################################################################################################
########################################################################################################
########################################################################################################

######------------Water------------############

CMRE5_Water <- subset_samples(CMRE5, Sampletype == "Water")
CMRE5_Watera = prune_taxa(taxa_sums(CMRE5_Water) > 0, CMRE5_Water) 
CMRE5_Waterb <- prune_taxa(taxa_sums(CMRE5_Watera) > 100, CMRE5_Watera)
##################################################
######------------Water-Ctr------------############
CMRE5_Waterb_Ctr <- subset_samples(CMRE5_Waterb, Treatment == "Ctr" & Timepoint %in%c("1", "7"))
CMRE5_Waterb_Ctra = prune_taxa(taxa_sums(CMRE5_Waterb_Ctr) > 0, CMRE5_Waterb_Ctr) #
##########Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(CMRE5_Waterb_Ctra),
               MARGIN = ifelse(taxa_are_rows(CMRE5_Waterb_Ctra), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(CMRE5_Waterb_Ctra),
                    tax_table(CMRE5_Waterb_Ctra))
prevalenceThreshold = 0.15 * nsamples(CMRE5_Waterb_Ctra) #
Waterb_Ctr_keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
CMRE5_Waterb_Ctrb = prune_taxa(Waterb_Ctr_keepTaxa, CMRE5_Waterb_Ctra)

##########Add taxonomic classification to OTU ID
CMRE5_Waterb_Ctrc <- format_to_besthit1(CMRE5_Waterb_Ctrb)# run this function instead (see end of this script page)
CMRE5_Waterb_Ctrc_otu <- t(otu_table(CMRE5_Waterb_Ctrc)@.Data) #extract the otu table from phyloseq object
CMRE5_Waterb_Ctrc_tax <- as.data.frame(tax_table(CMRE5_Waterb_Ctrc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Waterb_Ctrc_otu_net <- spiec.easi(CMRE5_Waterb_Ctrc_otu, method='mb', lambda.min.ratio=1e-3, nlambda=50, pulsar.params=list(rep.num=99)) # increaseing permutations can decrease the number of edges (removes unstable edges)

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Ctrc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
Water.ig= readRDS("CMRE5_Waterb_Ctrc_Water_ig.RDS")

nt_all_deg <- igraph::degree(Water.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Water.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Water.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Water.ig, "local", vids = V(Water.ig))
names(nt_all_transitivity)<- V(Water.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Water-As------------#############

CMRE5_Waterb_As <- subset_samples(CMRE5_Waterb, Treatment == "As" & Timepoint %in%c("1", "7"))
CMRE5_Waterb_Asa = prune_taxa(taxa_sums(CMRE5_Waterb_As) > 0, CMRE5_Waterb_As) #remove these OTUs from class object

CMRE5_Waterb_Asb = prune_taxa(Waterb_Ctr_keepTaxa, CMRE5_Waterb_Asa)
CMRE5_Waterb_Asc <- format_to_besthit1(CMRE5_Waterb_Asb)# run this function instead (see end of this script page)

CMRE5_Waterb_Asc_otu <- t(otu_table(CMRE5_Waterb_Asc)@.Data) #extract the otu table from phyloseq object
CMRE5_Waterb_Asc_tax <- as.data.frame(tax_table(CMRE5_Waterb_Asc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Waterb_Asc_otu_net <- spiec.easi(CMRE5_Waterb_Asc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Asc_otu_net)))

positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters

Water.ig= readRDS("CMRE5_Waterb_Asc_Water_ig.RDS")

nt_all_deg <- igraph::degree(Water.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Water.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Water.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Water.ig, "local", vids = V(Water.ig))
names(nt_all_transitivity)<- V(Water.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0
## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Water-Bx------------#############

CMRE5_Waterb_Bx <- subset_samples(CMRE5_Waterb, Treatment == "Bx" & Timepoint %in%c("1", "7"))
CMRE5_Waterb_Bxa = prune_taxa(taxa_sums(CMRE5_Waterb_Bx) > 0, CMRE5_Waterb_Bx) #
CMRE5_Waterb_Bxb = prune_taxa(Waterb_Ctr_keepTaxa, CMRE5_Waterb_Bxa)

CMRE5_Waterb_Bxc <- format_to_besthit1(CMRE5_Waterb_Bxb)# 
CMRE5_Waterb_Bxc_otu <- t(otu_table(CMRE5_Waterb_Bxc)@.Data) #
CMRE5_Waterb_Bxc_tax <- as.data.frame(tax_table(CMRE5_Waterb_Bxc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Waterb_Bxc_otu_net <- spiec.easi(CMRE5_Waterb_Bxc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) # 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Bxc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Water.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Water.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Water.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Water.ig, "local", vids = V(Water.ig))
names(nt_all_transitivity)<- V(Water.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Water-Tb------------#############

CMRE5_Waterb_Tb <- subset_samples(CMRE5_Waterb, Treatment == "Tb" & Timepoint %in%c("1", "7"))
CMRE5_Waterb_Tba = prune_taxa(taxa_sums(CMRE5_Waterb_Tb) > 0, CMRE5_Waterb_Tb) #

CMRE5_Waterb_Tbb = prune_taxa(Waterb_Ctr_keepTaxa, CMRE5_Waterb_Tba)
CMRE5_Waterb_Tbc <- format_to_besthit1(CMRE5_Waterb_Tbb)# 

CMRE5_Waterb_Tbc_otu <- t(otu_table(CMRE5_Waterb_Tbc)@.Data) #
CMRE5_Waterb_Tbc_tax <- as.data.frame(tax_table(CMRE5_Waterb_Tbc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Waterb_Tbc_otu_net <- spiec.easi(CMRE5_Waterb_Tbc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Waterb_Tbc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Water.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Water.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Water.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Water.ig, "local", vids = V(Water.ig))
names(nt_all_transitivity)<- V(Water.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##########Perform Kolmogorov-Smirnov test

Water_Ctr <- read.csv("Network_properties_Water_Ctr.csv")
Water_As <- read.csv("Network_properties_Water_As.csv")
Water_Bx <- read.csv("Network_properties_Water_Bx.csv")
Water_Tb <- read.csv("Network_properties_Water_Tb.csv")

ks.test(Water_Ctr$nt_boot_degree, Water_As$nt_boot_degree)
ks.test(Water_Ctr$nt_boot_degree, Water_Bx$nt_boot_degree)
ks.test(Water_Ctr$nt_boot_degree, Water_Tb$nt_boot_degree)
ks.test(Water_As$nt_boot_degree, Water_Bx$nt_boot_degree)
ks.test(Water_As$nt_boot_degree, Water_Tb$nt_boot_degree)
ks.test(Water_Bx$nt_boot_degree, Water_Tb$nt_boot_degree)
ks.test(Water_Ctr$nt_boot_betweenness, Water_As$nt_boot_betweenness)
ks.test(Water_Ctr$nt_boot_betweenness, Water_Bx$nt_boot_betweenness)
ks.test(Water_Ctr$nt_boot_betweenness, Water_Tb$nt_boot_betweenness)
ks.test(Water_As$nt_boot_betweenness, Water_Bx$nt_boot_betweenness)
ks.test(Water_As$nt_boot_betweenness, Water_Tb$nt_boot_betweenness)
ks.test(Water_Bx$nt_boot_betweenness, Water_Tb$nt_boot_betweenness)
ks.test(Water_Ctr$nt_boot_closeness, Water_As$nt_boot_closeness)
ks.test(Water_Ctr$nt_boot_closeness, Water_Bx$nt_boot_closeness)
ks.test(Water_Ctr$nt_boot_closeness, Water_Tb$nt_boot_closeness)
ks.test(Water_As$nt_boot_closeness, Water_Bx$nt_boot_closeness)
ks.test(Water_As$nt_boot_closeness, Water_Tb$nt_boot_closeness)
ks.test(Water_Bx$nt_boot_closeness, Water_Tb$nt_boot_closeness)
ks.test(Water_Ctr$nt_boot_transitivity, Water_As$nt_boot_transitivity)
ks.test(Water_Ctr$nt_boot_transitivity, Water_Bx$nt_boot_transitivity)
ks.test(Water_Ctr$nt_boot_transitivity, Water_Tb$nt_boot_transitivity)
ks.test(Water_As$nt_boot_transitivity, Water_Bx$nt_boot_transitivity)
ks.test(Water_As$nt_boot_transitivity, Water_Tb$nt_boot_transitivity)
ks.test(Water_Bx$nt_boot_transitivity, Water_Tb$nt_boot_transitivity)

##########Hub Scores
Water_Ctr_net.hs <- hub_score(Water.Ctr.Water.ig)$vector
Water_Ctr_net.hs.sort <- sort(Water_Ctr_net.hs, decreasing = TRUE)

Water_As_net.hs <- hub_score(Water.As.Water.ig)$vector
Water_As_net.hs.sort <- sort(Water_As_net.hs, decreasing = TRUE)

Water_Bx_net.hs <- hub_score(Water.Bx.Water.ig)$vector
Water_Bx_net.hs.sort <- sort(Water_Bx_net.hs, decreasing = TRUE)

Water_Tb_net.hs <- hub_score(Water.Tb.Water.ig)$vector
Water_Tb_net.hs.sort <- sort(Water_Tb_net.hs, decreasing = TRUE)

##########power law testing Water-Ctr
Water.ig= readRDS("CMRE5_Waterb_Ctrc_Water_ig.RDS")

Water.net <- asNetwork(Water.ig)
Water.ctr.spiec.deg <- degree(Water.net)

Water.ctr.spiec.deg_PL = displ$new(Water.ctr.spiec.deg)# power law distribution
Water.ctr.est = estimate_xmin(Water.ctr.spiec.deg_PL)
Water.ctr.spiec.deg_PL$xmin <- Water.ctr.est$xmin
Water.ctr.spiec.deg_PL$pars <- Water.ctr.est$pars

set.seed(101)
Water.ctr.bs <- bootstrap_p(Water.ctr.spiec.deg_PL, threads=4, no_of_sims = 1000)
Water.ctr.spiec.deg_N <- dislnorm$new(Water.ctr.spiec.deg)# Normal distribution
Water.ctr.spiec.deg_N$xmin <- Water.ctr.est$xmin
Water.ctr.spiec.deg_N$pars <- estimate_pars(Water.ctr.spiec.deg_N)
Water.ctr.bs <- bootstrap_p(Water.ctr.spiec.deg_N, threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.ctr.spiec.deg_PL, Water.ctr.spiec.deg_N)
Water.ctr.spiec.deg_E <- disexp$new(Water.ctr.spiec.deg)#  Exponential distribution
Water.ctr.spiec.deg_E$xmin <- Water.ctr.est$xmin
Water.ctr.spiec.deg_E$pars <- estimate_pars(Water.ctr.spiec.deg_E)

Water.ctr.bs <- bootstrap_p(Water.ctr.spiec.deg_E, threads=4, no_of_sims = 1000)
Water.comp <- compare_distributions(Water.ctr.spiec.deg_PL, Water.ctr.spiec.deg_E)
Water.ctr.spiec.deg_P <- dispois$new(Water.ctr.spiec.deg)# Poisson distribution
Water.ctr.spiec.deg_P$xmin <- Water.ctr.est$xmin
Water.ctr.spiec.deg_P$pars <- estimate_pars(Water.ctr.spiec.deg_P)
Water.ctr.bs <- bootstrap_p(Water.ctr.spiec.deg_P, threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.ctr.spiec.deg_PL, Water.ctr.spiec.deg_P)
Water.compNE <- compare_distributions(Water.ctr.spiec.deg_N, Water.ctr.spiec.deg_E)
Water.compEP <- compare_distributions(Water.ctr.spiec.deg_E, Water.ctr.spiec.deg_P)
Water.compNP <- compare_distributions(Water.ctr.spiec.deg_N, Water.ctr.spiec.deg_P)

##########power law testing Water-As
Water.ig= readRDS("CMRE5_Waterb_Asc_Water_ig.RDS")

Water.net <- asNetwork(Water.ig)
Water.As.spiec.deg <- degree(Water.net)

Water.As.spiec.deg_PL = displ$new(Water.As.spiec.deg)# power law distribution
Water.As.est = estimate_xmin(Water.As.spiec.deg_PL)
Water.As.spiec.deg_PL$setXmin(Water.As.est)
Water.As.spiec.deg_PL$xmin <- Water.As.est$xmin
Water.As.spiec.deg_PL$pars <- Water.As.est$pars

set.seed(101)
bs <- bootstrap_p(Water.As.spiec.deg_PL, threads=4, no_of_sims = 1000)

Water.As.spiec.deg_N <- dislnorm$new(Water.As.spiec.deg)# Normal distribution
Water.As.spiec.deg_N$xmin <- Water.As.est$xmin
Water.As.spiec.deg_N$pars <- estimate_pars(Water.As.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Water.As.spiec.deg_N, threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.As.spiec.deg_PL, Water.As.spiec.deg_N)
Water.As.spiec.deg_E <- disexp$new(Water.As.spiec.deg)#  Exponential distribution
Water.As.spiec.deg_E$xmin <- Water.As.est$xmin
Water.As.spiec.deg_E$pars <- estimate_pars(Water.As.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Water.As.spiec.deg_E, threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.As.spiec.deg_PL, Water.As.spiec.deg_E)
Water.As.spiec.deg_P <- dispois$new(Water.As.spiec.deg)# Poisson distribution
Water.As.spiec.deg_P$xmin <- Water.As.est$xmin
Water.As.spiec.deg_P$pars <- estimate_pars(Water.As.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Water.As.spiec.deg_P, threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.As.spiec.deg_PL, Water.As.spiec.deg_P)
Water.compNE <- compare_distributions(Water.As.spiec.deg_N, Water.As.spiec.deg_E)
Water.compEP <- compare_distributions(Water.As.spiec.deg_E, Water.As.spiec.deg_P)
Water.compNP <- compare_distributions(Water.As.spiec.deg_N, Water.As.spiec.deg_P)

##########power law testing Water-Bx
Water.ig= readRDS("CMRE5_Waterb_Bxc_Water_ig.RDS")
Water.net <- asNetwork(Water.ig)
Water.Bx.spiec.deg <- degree(Water.net)

Water.Bx.spiec.deg_PL = displ$new(Water.Bx.spiec.deg)# power law distribution
Water.Bx.est = estimate_xmin(Water.Bx.spiec.deg_PL)
Water.Bx.spiec.deg_PL$xmin <- Water.Bx.est$xmin
Water.Bx.spiec.deg_PL$pars <- Water.Bx.est$pars

set.seed(101)
bs <- bootstrap_p(Water.Bx.spiec.deg_PL,threads=4, no_of_sims = 1000)
Water.Bx.spiec.deg_N <- dislnorm$new(Water.Bx.spiec.deg)# Normal distribution
Water.Bx.spiec.deg_N$xmin <- Water.Bx.est$xmin
Water.Bx.spiec.deg_N$pars <- estimate_pars(Water.Bx.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Water.Bx.spiec.deg_N,threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.Bx.spiec.deg_PL, Water.Bx.spiec.deg_N)
Water.Bx.spiec.deg_E <- disexp$new(Water.Bx.spiec.deg)#  Exponential distribution
Water.Bx.spiec.deg_E$xmin <- Water.Bx.est$xmin
Water.Bx.spiec.deg_E$pars <- estimate_pars(Water.Bx.spiec.deg_E)

set.seed(101)
bs <- bootstrap_p(Water.Bx.spiec.deg_E,threads=4, no_of_sims = 1000)
Water.comp <- compare_distributions(Water.Bx.spiec.deg_PL, Water.Bx.spiec.deg_E)
Water.Bx.spiec.deg_P <- dispois$new(Water.Bx.spiec.deg)# Poisson distribution
Water.Bx.spiec.deg_P$xmin <- Water.Bx.est$xmin
Water.Bx.spiec.deg_P$pars <- estimate_pars(Water.Bx.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Water.Bx.spiec.deg_P,threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.Bx.spiec.deg_PL, Water.Bx.spiec.deg_P)
Water.compNE <- compare_distributions(Water.Bx.spiec.deg_N, Water.Bx.spiec.deg_E)
Water.compEP <- compare_distributions(Water.Bx.spiec.deg_E, Water.Bx.spiec.deg_P)
Water.compNP <- compare_distributions(Water.Bx.spiec.deg_N, Water.Bx.spiec.deg_P)

##########power law testing Water-Tb
Water.ig= readRDS("CMRE5_Waterb_Tbc_Water_ig.RDS")

Water.net <- asNetwork(Water.ig)
Water.Tb.spiec.deg <- degree(Water.net)

Water.Tb.spiec.deg_PL = displ$new(Water.Tb.spiec.deg)# power law distribution
est = estimate_xmin(Water.Tb.spiec.deg_PL)
Water.Tb.spiec.deg_PL$xmin <- est$xmin
Water.Tb.spiec.deg_PL$pars <- est$pars

set.seed(101)
bs <- bootstrap_p(Water.Tb.spiec.deg_PL,threads=4, no_of_sims = 1000)
Water.Tb.spiec.deg_N <- dislnorm$new(Water.Tb.spiec.deg)# Normal distribution
Water.Tb.spiec.deg_N$xmin <- est$xmin
Water.Tb.spiec.deg_N$pars <- estimate_pars(Water.Tb.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Water.Tb.spiec.deg_N,threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.Tb.spiec.deg_PL, Water.Tb.spiec.deg_N)
Water.Tb.spiec.deg_E <- disexp$new(Water.Tb.spiec.deg)#  Exponential distribution
Water.Tb.spiec.deg_E$xmin <- est$xmin
Water.Tb.spiec.deg_E$pars <- estimate_pars(Water.Tb.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Water.Tb.spiec.deg_E,threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.Tb.spiec.deg_PL, Water.Tb.spiec.deg_E)
Water.Tb.spiec.deg_P <- dispois$new(Water.Tb.spiec.deg)# Poisson distribution
Water.Tb.spiec.deg_P$xmin <- est$xmin
Water.Tb.spiec.deg_P$pars <- estimate_pars(Water.Tb.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Water.Tb.spiec.deg_P,threads=4, no_of_sims = 1000)

Water.comp <- compare_distributions(Water.Tb.spiec.deg_PL, Water.Tb.spiec.deg_P)
Water.compNE <- compare_distributions(Water.Tb.spiec.deg_N, Water.Tb.spiec.deg_E)
Water.compEP <- compare_distributions(Water.Tb.spiec.deg_E, Water.Tb.spiec.deg_P)
Water.compNP <- compare_distributions(Water.Tb.spiec.deg_N, Water.Tb.spiec.deg_P)


########################################################################################################
########################################################################################################
########################################################################################################

######------------Sediment------------############

CMRE5_Sediment <- subset_samples(CMRE5, Sampletype == "sedi")
CMRE5_Sedimenta = prune_taxa(taxa_sums(CMRE5_Sediment) > 0, CMRE5_Sediment) 
CMRE5_Sedimentb <- prune_taxa(taxa_sums(CMRE5_Sedimenta) > 100, CMRE5_Sedimenta)
##################################################
######------------Sediment-Ctr------------############
CMRE5_Sedimentb_Ctr <- subset_samples(CMRE5_Sedimentb, Treatment == "Ctr" & Timepoint %in%c("1", "7"))
CMRE5_Sedimentb_Ctra = prune_taxa(taxa_sums(CMRE5_Sedimentb_Ctr) > 0, CMRE5_Sedimentb_Ctr) #
##########Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(CMRE5_Sedimentb_Ctra),
               MARGIN = ifelse(taxa_are_rows(CMRE5_Sedimentb_Ctra), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(CMRE5_Sedimentb_Ctra),
                    tax_table(CMRE5_Sedimentb_Ctra))
prevalenceThreshold = 0.15 * nsamples(CMRE5_Sedimentb_Ctra) #
Sedimentb_Ctr_keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
CMRE5_Sedimentb_Ctrb = prune_taxa(Sedimentb_Ctr_keepTaxa, CMRE5_Sedimentb_Ctra)

##########Add taxonomic classification to OTU ID
CMRE5_Sedimentb_Ctrc <- format_to_besthit1(CMRE5_Sedimentb_Ctrb)# run this function instead (see end of this script page)
CMRE5_Sedimentb_Ctrc_otu <- t(otu_table(CMRE5_Sedimentb_Ctrc)@.Data) #extract the otu table from phyloseq object
CMRE5_Sedimentb_Ctrc_tax <- as.data.frame(tax_table(CMRE5_Sedimentb_Ctrc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Sedimentb_Ctrc_otu_net <- spiec.easi(CMRE5_Sedimentb_Ctrc_otu, method='mb', lambda.min.ratio=1e-3, nlambda=50, pulsar.params=list(rep.num=99)) # increaseing permutations can decrease the number of edges (removes unstable edges)

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Ctrc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
Sediment.ig= readRDS("CMRE5_Sedimentb_Ctrc_Sediment_ig.RDS")

nt_all_deg <- igraph::degree(Sediment.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Sediment.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Sediment.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Sediment.ig, "local", vids = V(Sediment.ig))
names(nt_all_transitivity)<- V(Sediment.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Sediment-As------------#############

CMRE5_Sedimentb_As <- subset_samples(CMRE5_Sedimentb, Treatment == "As" & Timepoint %in%c("1", "7"))
CMRE5_Sedimentb_Asa = prune_taxa(taxa_sums(CMRE5_Sedimentb_As) > 0, CMRE5_Sedimentb_As) #remove these OTUs from class object

CMRE5_Sedimentb_Asb = prune_taxa(Sedimentb_Ctr_keepTaxa, CMRE5_Sedimentb_Asa)
CMRE5_Sedimentb_Asc <- format_to_besthit1(CMRE5_Sedimentb_Asb)# run this function instead (see end of this script page)

CMRE5_Sedimentb_Asc_otu <- t(otu_table(CMRE5_Sedimentb_Asc)@.Data) #extract the otu table from phyloseq object
CMRE5_Sedimentb_Asc_tax <- as.data.frame(tax_table(CMRE5_Sedimentb_Asc)@.Data)#extract the taxonomy information

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Sedimentb_Asc_otu_net <- spiec.easi(CMRE5_Sedimentb_Asc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Asc_otu_net)))

positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters

Sediment.ig= readRDS("CMRE5_Sedimentb_Asc_Sediment_ig.RDS")

nt_all_deg <- igraph::degree(Sediment.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Sediment.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Sediment.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Sediment.ig, "local", vids = V(Sediment.ig))
names(nt_all_transitivity)<- V(Sediment.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0
## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Sediment-Bx------------#############

CMRE5_Sedimentb_Bx <- subset_samples(CMRE5_Sedimentb, Treatment == "Bx" & Timepoint %in%c("1", "7"))
CMRE5_Sedimentb_Bxa = prune_taxa(taxa_sums(CMRE5_Sedimentb_Bx) > 0, CMRE5_Sedimentb_Bx) #
CMRE5_Sedimentb_Bxb = prune_taxa(Sedimentb_Ctr_keepTaxa, CMRE5_Sedimentb_Bxa)

CMRE5_Sedimentb_Bxc <- format_to_besthit1(CMRE5_Sedimentb_Bxb)# 
CMRE5_Sedimentb_Bxc_otu <- t(otu_table(CMRE5_Sedimentb_Bxc)@.Data) #
CMRE5_Sedimentb_Bxc_tax <- as.data.frame(tax_table(CMRE5_Sedimentb_Bxc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Sedimentb_Bxc_otu_net <- spiec.easi(CMRE5_Sedimentb_Bxc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) # 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Bxc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Sediment.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Sediment.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Sediment.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Sediment.ig, "local", vids = V(Sediment.ig))
names(nt_all_transitivity)<- V(Sediment.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##################################################
######------------Sediment-Tb------------#############

CMRE5_Sedimentb_Tb <- subset_samples(CMRE5_Sedimentb, Treatment == "Tb" & Timepoint %in%c("1", "7"))
CMRE5_Sedimentb_Tba = prune_taxa(taxa_sums(CMRE5_Sedimentb_Tb) > 0, CMRE5_Sedimentb_Tb) #

CMRE5_Sedimentb_Tbb = prune_taxa(Sedimentb_Ctr_keepTaxa, CMRE5_Sedimentb_Tba)
CMRE5_Sedimentb_Tbc <- format_to_besthit1(CMRE5_Sedimentb_Tbb)# 

CMRE5_Sedimentb_Tbc_otu <- t(otu_table(CMRE5_Sedimentb_Tbc)@.Data) #
CMRE5_Sedimentb_Tbc_tax <- as.data.frame(tax_table(CMRE5_Sedimentb_Tbc)@.Data)#

##########SPIEC-EASI network reconstruction
set.seed(101)
CMRE5_Sedimentb_Tbc_otu_net <- spiec.easi(CMRE5_Sedimentb_Tbc_otu, method='mb', lambda.min.ratio=1e-1, nlambda=50, pulsar.params=list(rep.num=99)) 

##########Network properties
betaMat=as.matrix(symBeta(getOptBeta(CMRE5_Sedimentb_Tbc_otu_net)))
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 

##########Prepare data for plotting  
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

##########Calculate network parameters
nt_all_deg <- igraph::degree(Sediment.ig, mode="all")
nt_all_betweenness <- igraph::betweenness(Sediment.ig, normalized = TRUE)
nt_all_closeness <- igraph::closeness(Sediment.ig, normalized = TRUE)#not working
nt_all_transitivity <- igraph::transitivity(Sediment.ig, "local", vids = V(Sediment.ig))
names(nt_all_transitivity)<- V(Sediment.ig)$name
nt_all_transitivity[is.na(nt_all_transitivity)] <- 0

## Bootstrapping degree
set.seed(8046)
nt_boot_degree = replicate(10000, sample(nt_all_deg, 1, replace=TRUE))

## Bootstrapping betweenness
set.seed(8046)
nt_boot_betweenness = replicate(10000, sample(nt_all_betweenness, 1, replace=TRUE))

## Bootstrapping closeness
set.seed(8046)
nt_boot_closeness = replicate(10000, sample(nt_all_closeness, 1, replace=TRUE))

## Bootstrapping transitivity
set.seed(8046)
nt_boot_transitivity = replicate(10000, sample(nt_all_transitivity, 1, replace=TRUE))
notill_node_characterstics<-cbind(nt_boot_degree,nt_boot_betweenness,nt_boot_closeness, nt_boot_transitivity)

##########Perform Kolmogorov-Smirnov test

Sediment_Ctr <- read.csv("Network_properties_Sediment_Ctr.csv")
Sediment_As <- read.csv("Network_properties_Sediment_As.csv")
Sediment_Bx <- read.csv("Network_properties_Sediment_Bx.csv")
Sediment_Tb <- read.csv("Network_properties_Sediment_Tb.csv")

ks.test(Sediment_Ctr$nt_boot_degree, Sediment_As$nt_boot_degree)
ks.test(Sediment_Ctr$nt_boot_degree, Sediment_Bx$nt_boot_degree)
ks.test(Sediment_Ctr$nt_boot_degree, Sediment_Tb$nt_boot_degree)
ks.test(Sediment_As$nt_boot_degree, Sediment_Bx$nt_boot_degree)
ks.test(Sediment_As$nt_boot_degree, Sediment_Tb$nt_boot_degree)
ks.test(Sediment_Bx$nt_boot_degree, Sediment_Tb$nt_boot_degree)
ks.test(Sediment_Ctr$nt_boot_betweenness, Sediment_As$nt_boot_betweenness)
ks.test(Sediment_Ctr$nt_boot_betweenness, Sediment_Bx$nt_boot_betweenness)
ks.test(Sediment_Ctr$nt_boot_betweenness, Sediment_Tb$nt_boot_betweenness)
ks.test(Sediment_As$nt_boot_betweenness, Sediment_Bx$nt_boot_betweenness)
ks.test(Sediment_As$nt_boot_betweenness, Sediment_Tb$nt_boot_betweenness)
ks.test(Sediment_Bx$nt_boot_betweenness, Sediment_Tb$nt_boot_betweenness)
ks.test(Sediment_Ctr$nt_boot_closeness, Sediment_As$nt_boot_closeness)
ks.test(Sediment_Ctr$nt_boot_closeness, Sediment_Bx$nt_boot_closeness)
ks.test(Sediment_Ctr$nt_boot_closeness, Sediment_Tb$nt_boot_closeness)
ks.test(Sediment_As$nt_boot_closeness, Sediment_Bx$nt_boot_closeness)
ks.test(Sediment_As$nt_boot_closeness, Sediment_Tb$nt_boot_closeness)
ks.test(Sediment_Bx$nt_boot_closeness, Sediment_Tb$nt_boot_closeness)
ks.test(Sediment_Ctr$nt_boot_transitivity, Sediment_As$nt_boot_transitivity)
ks.test(Sediment_Ctr$nt_boot_transitivity, Sediment_Bx$nt_boot_transitivity)
ks.test(Sediment_Ctr$nt_boot_transitivity, Sediment_Tb$nt_boot_transitivity)
ks.test(Sediment_As$nt_boot_transitivity, Sediment_Bx$nt_boot_transitivity)
ks.test(Sediment_As$nt_boot_transitivity, Sediment_Tb$nt_boot_transitivity)
ks.test(Sediment_Bx$nt_boot_transitivity, Sediment_Tb$nt_boot_transitivity)

##########Hub Scores
Sediment_Ctr_net.hs <- hub_score(Sediment.Ctr.Sediment.ig)$vector
Sediment_Ctr_net.hs.sort <- sort(Sediment_Ctr_net.hs, decreasing = TRUE)

Sediment_As_net.hs <- hub_score(Sediment.As.Sediment.ig)$vector
Sediment_As_net.hs.sort <- sort(Sediment_As_net.hs, decreasing = TRUE)

Sediment_Bx_net.hs <- hub_score(Sediment.Bx.Sediment.ig)$vector
Sediment_Bx_net.hs.sort <- sort(Sediment_Bx_net.hs, decreasing = TRUE)

Sediment_Tb_net.hs <- hub_score(Sediment.Tb.Sediment.ig)$vector
Sediment_Tb_net.hs.sort <- sort(Sediment_Tb_net.hs, decreasing = TRUE)

##########power law testing Sediment-Ctr
Sediment.ig= readRDS("CMRE5_Sedimentb_Ctrc_Sediment_ig.RDS")

Sediment.net <- asNetwork(Sediment.ig)
Sediment.ctr.spiec.deg <- degree(Sediment.net)

Sediment.ctr.spiec.deg_PL = displ$new(Sediment.ctr.spiec.deg)# power law distribution
Sediment.ctr.est = estimate_xmin(Sediment.ctr.spiec.deg_PL)
Sediment.ctr.spiec.deg_PL$xmin <- Sediment.ctr.est$xmin
Sediment.ctr.spiec.deg_PL$pars <- Sediment.ctr.est$pars

set.seed(101)
Sediment.ctr.bs <- bootstrap_p(Sediment.ctr.spiec.deg_PL, threads=4, no_of_sims = 1000)
Sediment.ctr.spiec.deg_N <- dislnorm$new(Sediment.ctr.spiec.deg)# Normal distribution
Sediment.ctr.spiec.deg_N$xmin <- Sediment.ctr.est$xmin
Sediment.ctr.spiec.deg_N$pars <- estimate_pars(Sediment.ctr.spiec.deg_N)
Sediment.ctr.bs <- bootstrap_p(Sediment.ctr.spiec.deg_N, threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.ctr.spiec.deg_PL, Sediment.ctr.spiec.deg_N)
Sediment.ctr.spiec.deg_E <- disexp$new(Sediment.ctr.spiec.deg)#  Exponential distribution
Sediment.ctr.spiec.deg_E$xmin <- Sediment.ctr.est$xmin
Sediment.ctr.spiec.deg_E$pars <- estimate_pars(Sediment.ctr.spiec.deg_E)

Sediment.ctr.bs <- bootstrap_p(Sediment.ctr.spiec.deg_E, threads=4, no_of_sims = 1000)
Sediment.comp <- compare_distributions(Sediment.ctr.spiec.deg_PL, Sediment.ctr.spiec.deg_E)
Sediment.ctr.spiec.deg_P <- dispois$new(Sediment.ctr.spiec.deg)# Poisson distribution
Sediment.ctr.spiec.deg_P$xmin <- Sediment.ctr.est$xmin
Sediment.ctr.spiec.deg_P$pars <- estimate_pars(Sediment.ctr.spiec.deg_P)
Sediment.ctr.bs <- bootstrap_p(Sediment.ctr.spiec.deg_P, threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.ctr.spiec.deg_PL, Sediment.ctr.spiec.deg_P)
Sediment.compNE <- compare_distributions(Sediment.ctr.spiec.deg_N, Sediment.ctr.spiec.deg_E)
Sediment.compEP <- compare_distributions(Sediment.ctr.spiec.deg_E, Sediment.ctr.spiec.deg_P)
Sediment.compNP <- compare_distributions(Sediment.ctr.spiec.deg_N, Sediment.ctr.spiec.deg_P)

##########power law testing Sediment-As
Sediment.ig= readRDS("CMRE5_Sedimentb_Asc_Sediment_ig.RDS")

Sediment.net <- asNetwork(Sediment.ig)
Sediment.As.spiec.deg <- degree(Sediment.net)

Sediment.As.spiec.deg_PL = displ$new(Sediment.As.spiec.deg)# power law distribution
Sediment.As.est = estimate_xmin(Sediment.As.spiec.deg_PL)
Sediment.As.spiec.deg_PL$setXmin(Sediment.As.est)
Sediment.As.spiec.deg_PL$xmin <- Sediment.As.est$xmin
Sediment.As.spiec.deg_PL$pars <- Sediment.As.est$pars

set.seed(101)
bs <- bootstrap_p(Sediment.As.spiec.deg_PL, threads=4, no_of_sims = 1000)

Sediment.As.spiec.deg_N <- dislnorm$new(Sediment.As.spiec.deg)# Normal distribution
Sediment.As.spiec.deg_N$xmin <- Sediment.As.est$xmin
Sediment.As.spiec.deg_N$pars <- estimate_pars(Sediment.As.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Sediment.As.spiec.deg_N, threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.As.spiec.deg_PL, Sediment.As.spiec.deg_N)
Sediment.As.spiec.deg_E <- disexp$new(Sediment.As.spiec.deg)#  Exponential distribution
Sediment.As.spiec.deg_E$xmin <- Sediment.As.est$xmin
Sediment.As.spiec.deg_E$pars <- estimate_pars(Sediment.As.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Sediment.As.spiec.deg_E, threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.As.spiec.deg_PL, Sediment.As.spiec.deg_E)
Sediment.As.spiec.deg_P <- dispois$new(Sediment.As.spiec.deg)# Poisson distribution
Sediment.As.spiec.deg_P$xmin <- Sediment.As.est$xmin
Sediment.As.spiec.deg_P$pars <- estimate_pars(Sediment.As.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Sediment.As.spiec.deg_P, threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.As.spiec.deg_PL, Sediment.As.spiec.deg_P)
Sediment.compNE <- compare_distributions(Sediment.As.spiec.deg_N, Sediment.As.spiec.deg_E)
Sediment.compEP <- compare_distributions(Sediment.As.spiec.deg_E, Sediment.As.spiec.deg_P)
Sediment.compNP <- compare_distributions(Sediment.As.spiec.deg_N, Sediment.As.spiec.deg_P)

##########power law testing Sediment-Bx
Sediment.ig= readRDS("CMRE5_Sedimentb_Bxc_Sediment_ig.RDS")
Sediment.net <- asNetwork(Sediment.ig)
Sediment.Bx.spiec.deg <- degree(Sediment.net)

Sediment.Bx.spiec.deg_PL = displ$new(Sediment.Bx.spiec.deg)# power law distribution
Sediment.Bx.est = estimate_xmin(Sediment.Bx.spiec.deg_PL)
Sediment.Bx.spiec.deg_PL$xmin <- Sediment.Bx.est$xmin
Sediment.Bx.spiec.deg_PL$pars <- Sediment.Bx.est$pars

set.seed(101)
bs <- bootstrap_p(Sediment.Bx.spiec.deg_PL,threads=4, no_of_sims = 1000)
Sediment.Bx.spiec.deg_N <- dislnorm$new(Sediment.Bx.spiec.deg)# Normal distribution
Sediment.Bx.spiec.deg_N$xmin <- Sediment.Bx.est$xmin
Sediment.Bx.spiec.deg_N$pars <- estimate_pars(Sediment.Bx.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Sediment.Bx.spiec.deg_N,threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.Bx.spiec.deg_PL, Sediment.Bx.spiec.deg_N)
Sediment.Bx.spiec.deg_E <- disexp$new(Sediment.Bx.spiec.deg)#  Exponential distribution
Sediment.Bx.spiec.deg_E$xmin <- Sediment.Bx.est$xmin
Sediment.Bx.spiec.deg_E$pars <- estimate_pars(Sediment.Bx.spiec.deg_E)

set.seed(101)
bs <- bootstrap_p(Sediment.Bx.spiec.deg_E,threads=4, no_of_sims = 1000)
Sediment.comp <- compare_distributions(Sediment.Bx.spiec.deg_PL, Sediment.Bx.spiec.deg_E)
Sediment.Bx.spiec.deg_P <- dispois$new(Sediment.Bx.spiec.deg)# Poisson distribution
Sediment.Bx.spiec.deg_P$xmin <- Sediment.Bx.est$xmin
Sediment.Bx.spiec.deg_P$pars <- estimate_pars(Sediment.Bx.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Sediment.Bx.spiec.deg_P,threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.Bx.spiec.deg_PL, Sediment.Bx.spiec.deg_P)
Sediment.compNE <- compare_distributions(Sediment.Bx.spiec.deg_N, Sediment.Bx.spiec.deg_E)
Sediment.compEP <- compare_distributions(Sediment.Bx.spiec.deg_E, Sediment.Bx.spiec.deg_P)
Sediment.compNP <- compare_distributions(Sediment.Bx.spiec.deg_N, Sediment.Bx.spiec.deg_P)

##########power law testing Sediment-Tb
Sediment.ig= readRDS("CMRE5_Sedimentb_Tbc_Sediment_ig.RDS")

Sediment.net <- asNetwork(Sediment.ig)
Sediment.Tb.spiec.deg <- degree(Sediment.net)

Sediment.Tb.spiec.deg_PL = displ$new(Sediment.Tb.spiec.deg)# power law distribution
est = estimate_xmin(Sediment.Tb.spiec.deg_PL)
Sediment.Tb.spiec.deg_PL$xmin <- est$xmin
Sediment.Tb.spiec.deg_PL$pars <- est$pars

set.seed(101)
bs <- bootstrap_p(Sediment.Tb.spiec.deg_PL,threads=4, no_of_sims = 1000)
Sediment.Tb.spiec.deg_N <- dislnorm$new(Sediment.Tb.spiec.deg)# Normal distribution
Sediment.Tb.spiec.deg_N$xmin <- est$xmin
Sediment.Tb.spiec.deg_N$pars <- estimate_pars(Sediment.Tb.spiec.deg_N)
set.seed(101)
bs <- bootstrap_p(Sediment.Tb.spiec.deg_N,threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.Tb.spiec.deg_PL, Sediment.Tb.spiec.deg_N)
Sediment.Tb.spiec.deg_E <- disexp$new(Sediment.Tb.spiec.deg)#  Exponential distribution
Sediment.Tb.spiec.deg_E$xmin <- est$xmin
Sediment.Tb.spiec.deg_E$pars <- estimate_pars(Sediment.Tb.spiec.deg_E)
set.seed(101)
bs <- bootstrap_p(Sediment.Tb.spiec.deg_E,threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.Tb.spiec.deg_PL, Sediment.Tb.spiec.deg_E)
Sediment.Tb.spiec.deg_P <- dispois$new(Sediment.Tb.spiec.deg)# Poisson distribution
Sediment.Tb.spiec.deg_P$xmin <- est$xmin
Sediment.Tb.spiec.deg_P$pars <- estimate_pars(Sediment.Tb.spiec.deg_P)
set.seed(101)
bs <- bootstrap_p(Sediment.Tb.spiec.deg_P,threads=4, no_of_sims = 1000)

Sediment.comp <- compare_distributions(Sediment.Tb.spiec.deg_PL, Sediment.Tb.spiec.deg_P)
Sediment.compNE <- compare_distributions(Sediment.Tb.spiec.deg_N, Sediment.Tb.spiec.deg_E)
Sediment.compEP <- compare_distributions(Sediment.Tb.spiec.deg_E, Sediment.Tb.spiec.deg_P)
Sediment.compNP <- compare_distributions(Sediment.Tb.spiec.deg_N, Sediment.Tb.spiec.deg_P)


###################################################################
######------------Combined network plot------------################

combine_plot1 <- ggarrange(W_Ctr, Se_Ctr, S_Ctr, R_Ctr, M_Ctr, W_As, Se_As, S_As, R_As,  M_As,  W_Bx, Se_bx, S_Bx, R_Bx, M_Bx, W_Tb, Se_Tb, S_Tb, R_Tb, M_Tb, ncol = 5, nrow = 4)

###################################################################
######------------Degree distribution plot------------############# 

############## Mouse Degree distributions ######################

Mouse.Ctr.feces.ig= readRDS("Mouse.Ctr.feces.ig.RDS")
Mouse.As.feces.ig= readRDS("Mouse.As.feces.ig.RDS")
Mouse.Bx.feces.ig= readRDS("Mouse.Bx.feces.ig.RDS")
Mouse.Tb.feces.ig= readRDS("Mouse.Tb.feces.ig.RDS")

Mouse.Ctr_all_deg <- as.data.frame(igraph::degree(Mouse.Ctr.feces.ig, mode="all")) # one can also extract other values of network by using such
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

Root.Ctr_all_deg <- as.data.frame(igraph::degree(Root.Ctr.Root.ig, mode="all")) # one can also extract other values of network by using such
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

Soil.Ctr_all_deg <- as.data.frame(igraph::degree(Soil.Ctr.Soil.ig, mode="all")) # one can also extract other values of network by using such
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

Water.Ctr_all_deg <- as.data.frame(igraph::degree(Water.Ctr.Water.ig, mode="all")) # one can also extract other values of network by using such
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

#############################################################################################

############## Sediment Degree distributions ######################

Sediment.Ctr.Sediment.ig= readRDS("Sediment.Ctr.Sediment.ig.RDS")
Sediment.As.Sediment.ig= readRDS("Sediment.As.Sediment.ig.RDS")
Sediment.Bx.Sediment.ig= readRDS("Sediment.Bx.Sediment.ig.RDS")
Sediment.Tb.Sediment.ig= readRDS("Sediment.Tb.Sediment.ig.RDS")

Sediment.Ctr_all_deg <- as.data.frame(igraph::degree(Sediment.Ctr.Sediment.ig, mode="all")) # one can also extract other values of network by using such
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

#############################################################################################

######combined degree distribution plot ####################

combine_degree_plot1 <- ggarrange(W_Ctr, Se_Ctr, S_Ctr, R_Ctr, M_Ctr, W_As, Se_As, S_As, R_As,  M_As,  W_Bx, Se_bx, S_Bx, R_Bx, M_Bx, W_Tb, Se_Tb, S_Tb, R_Tb, M_Tb, Water_p1, Sediment_p1, Soil_p1, Root_p1, Mouse_p1, ncol = 5, nrow = 5)

################################################################################################################################################################
