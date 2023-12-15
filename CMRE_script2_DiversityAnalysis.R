
#Copyright (c) 2023 University of Bern, Bern, Switzerland. 
#Author: Wasimuddin (wasim.bt@gmail.com; https://www.researchgate.net/profile/Wasim-Uddin)


# library(dplyr); packageVersion("dplyr")
# ## [1] '0.8.5'
# library(tidyr); packageVersion("tidyr")
# ## [1] '0.8.1'
# library(ggplot2); packageVersion("ggplot2") 
# ## [1] '3.1.1'
# library(phyloseq); packageVersion("phyloseq") 
# ## [1] ‘1.26.1’
# library("data.table"); packageVersion("data.table")
# ## [1] ‘1.11.6’
# library("vegan"); packageVersion("vegan")
# ## [1] ‘2.5.5’
# library(lme4); packageVersion("lme4")
# ## [1] ‘1.1.30’
# library(car); packageVersion("car")
# ## [1] ‘3.1.0’
# library(devtools); packageVersion("devtools")
# ## [1] ‘2.4.3’
# library(pairwiseAdonis)

####################################--------rarefaction at 8100----------####################
set.seed(301)
CMRE5_r8100 = rarefy_even_depth(CMRE5, sample.size = 8100, rngseed = TRUE) #rarefaction
CMRE5_r8100_t  = transform_sample_counts(CMRE5_r8100, function(x) x / sum(x) ) # transforming sample count (converting into relative abundance)

##################################-------Alpha diversity  --------------###################

##################################-------Calculate diversity indices

CMRE5_r8100_richness <- estimate_richness(CMRE5_r8100, measures = c("Observed", "Shannon"))

###########-Here metadata table is merged with alpha diversity table in order to easy downstream analysis. 
CMRE5_r8100_richness$sample_names <- rownames(CMRE5_r8100_richness)
CMRE5_r8100_Metadata_table = data.frame(as(sample_data(CMRE5_r8100), "data.frame"))
CMRE5_r8100_Metadata_table$sample_names <- rownames(CMRE5_r8100_Metadata_table)
CMRE5_Metadata_with_AlphaDiv_r8100 <- CMRE5_r8100_Metadata_table %>% left_join(CMRE5_r8100_richness)
rownames(CMRE5_Metadata_with_AlphaDiv_r8100) <- CMRE5_Metadata_with_AlphaDiv_r8100$sample_names

##################################------Running Alpha diversity Linear models------------------#################
Observed <- (CMRE5_alpha_r8100$Observed)
Treatment <- (CMRE5_alpha_r8100$Treatment)
Timepoint <- (CMRE5_alpha_r8100$Timepoint)
Concentration <- (CMRE5_alpha_r8100$Concentration)
Lib <- (CMRE5_alpha_r8100$Lib)
Sampletype <- (CMRE5_alpha_r8100$Sampletype)
Shannon <- (CMRE5_alpha_r8100$Shannon)
Fisher <- (CMRE5_alpha_r8100$Fisher)

##################---------------Observed
lm_M1Observed_8100_I <- aov(Observed ~  Treatment * as.factor(Timepoint) * as.factor(Concentration) * Sampletype + Lib , na.action=na.fail)
Anova(lm_M1Observed_8100_I, type="II")

##############---------- Tukey on linear model
lm_M1Observed_8100Tukey <- TukeyHSD(lm_M1Observed_8100_I) 

##################---------------Shannon
lm_M1Shannon_8100_I <- aov(Shannon ~  Treatment * as.factor(Timepoint) * as.factor(Concentration) * Sampletype + Lib , na.action=na.fail)
Anova(lm_M1Shannon_8100_I, type="II")

##############---------- Tukey on linear model
lm_M1Shannon_8100Tukey <- TukeyHSD(lm_M1Shannon_8100_I) 

######################################################################################################
###########################----------------Beta diversity analysis------------########################
########################--Bray-curtis 
set.seed(301)
CMRE5_r8100_t_br <- distance(CMRE5_r8100_t, method = "bray") 
CMRE5_r8100_t_br_adonis1 <- adonis(CMRE5_r8100_t_br ~  Treatment * as.factor(Timepoint) * as.factor(Concentration) * Sampletype + Lib, as(sample_data(CMRE5_r8100_t), "data.frame"), permutations = 999) #

########################--Jaccard 
set.seed(301)
CMRE5_r8100_t_Ja <- distance(CMRE5_r8100_t, method = "jaccard") 
CMRE5_r8100_t_Ja_adonis1 <- adonis(CMRE5_r8100_t_Ja ~  Treatment * as.factor(Timepoint) * as.factor(Concentration) * Sampletype + Lib, as(sample_data(CMRE5_r8100_t), "data.frame"), permutations = 999) #

##########-- Pairwise adonis tests (Bray)----#########
CMRE5_r8100_t_metadata <- sample_data(CMRE5_r8100_t)
CMRE5_r8100_t_br_adonis2 <- pairwise.adonis(CMRE5_r8100_t_br,  CMRE5_r8100_t_metadata$Treatment,  CMRE5_r8100_t_metadata$Sampletype) #

Treatment <-  CMRE5_r8100_t_metadata$Treatment
Timepoint <-  CMRE5_r8100_t_metadata$Timepoint
Concentration <-  CMRE5_r8100_t_metadata$Concentration
Sampletype <-  CMRE5_r8100_t_metadata$Sampletype

# Now prepare the interaction terms as variable to test
TTCS <- interaction(Treatment,  as.factor(Timepoint) , as.factor(Concentration) , Sampletype )
TTC <- interaction(Treatment,  as.factor(Timepoint) , as.factor(Concentration) )
TTS <- interaction(Treatment,  as.factor(Timepoint) , Sampletype )
TCS <- interaction(Treatment, as.factor(Concentration) , Sampletype )
TCS1 <- interaction(as.factor(Timepoint) , as.factor(Concentration) , Sampletype )
TT <- interaction(Treatment,  as.factor(Timepoint))
TC <- interaction(Treatment, as.factor(Concentration) )
TC1 <- interaction(as.factor(Timepoint), as.factor(Concentration) )
TS <- interaction(Treatment,  Sampletype )
TS1 <- interaction(as.factor(Timepoint),  Sampletype )
CS <- interaction(as.factor(Concentration),  Sampletype )

CMRE5_r8100_t_br_adonis_TS <- pairwise.adonis(CMRE5_r8100_t_br,  TS)
CMRE5_r8100_t_br_adonis_TTCS <- pairwise.adonis(CMRE5_r8100_t_br,  TTCS)
CMRE5_r8100_t_br_adonis_TTC <- pairwise.adonis(CMRE5_r8100_t_br,  TTC)
CMRE5_r8100_t_br_adonis_TTS <- pairwise.adonis(CMRE5_r8100_t_br,  TTS)
CMRE5_r8100_t_br_adonis_TCS <- pairwise.adonis(CMRE5_r8100_t_br,  TCS)
CMRE5_r8100_t_br_adonis_TCS1 <- pairwise.adonis(CMRE5_r8100_t_br,  TCS1)
CMRE5_r8100_t_br_adonis_TT <- pairwise.adonis(CMRE5_r8100_t_br,  TT)
CMRE5_r8100_t_br_adonis_TC <- pairwise.adonis(CMRE5_r8100_t_br,  TC)
CMRE5_r8100_t_br_adonis_TC1 <- pairwise.adonis(CMRE5_r8100_t_br,  TC1)
CMRE5_r8100_t_br_adonis_TS1 <- pairwise.adonis(CMRE5_r8100_t_br,  TS1)
CMRE5_r8100_t_br_adonis_CS <- pairwise.adonis(CMRE5_r8100_t_br,  CS)

##########-- Pairwise adonis tests (Jaccard)----#########
CMRE5_r8100_t_metadata <- sample_data(CMRE5_r8100_t)
CMRE5_r8100_t_ja_adonis2 <- pairwise.adonis(CMRE5_r8100_t_Ja,  CMRE5_r8100_t_metadata$Treatment,  CMRE5_r8100_t_metadata$Sampletype) #

Treatment <-  CMRE5_r8100_t_metadata$Treatment
Timepoint <-  CMRE5_r8100_t_metadata$Timepoint
Concentration <-  CMRE5_r8100_t_metadata$Concentration
Sampletype <-  CMRE5_r8100_t_metadata$Sampletype

# Now prepare the interaction terms as variable to test
TTCS <- interaction(Treatment,  as.factor(Timepoint) , as.factor(Concentration) , Sampletype )
TTC <- interaction(Treatment,  as.factor(Timepoint) , as.factor(Concentration) )
TTS <- interaction(Treatment,  as.factor(Timepoint) , Sampletype )
TCS <- interaction(Treatment, as.factor(Concentration) , Sampletype )
TCS1 <- interaction(as.factor(Timepoint) , as.factor(Concentration) , Sampletype )
TT <- interaction(Treatment,  as.factor(Timepoint))
TC <- interaction(Treatment, as.factor(Concentration) )
TC1 <- interaction(as.factor(Timepoint), as.factor(Concentration) )
TS <- interaction(Treatment,  Sampletype )
TS1 <- interaction(as.factor(Timepoint),  Sampletype )
CS <- interaction(as.factor(Concentration),  Sampletype )

CMRE5_r8100_t_ja_adonis_TS <- pairwise.adonis(CMRE5_r8100_t_Ja,  TS)
CMRE5_r8100_t_ja_adonis_TTCS <- pairwise.adonis(CMRE5_r8100_t_ja,  TTCS)
CMRE5_r8100_t_ja_adonis_TTC <- pairwise.adonis(CMRE5_r8100_t_ja,  TTC)
CMRE5_r8100_t_ja_adonis_TTS <- pairwise.adonis(CMRE5_r8100_t_ja,  TTS)
CMRE5_r8100_t_ja_adonis_TCS <- pairwise.adonis(CMRE5_r8100_t_ja,  TCS)
CMRE5_r8100_t_ja_adonis_TCS1 <- pairwise.adonis(CMRE5_r8100_t_ja,  TCS1)
CMRE5_r8100_t_ja_adonis_TT <- pairwise.adonis(CMRE5_r8100_t_ja,  TT)
CMRE5_r8100_t_ja_adonis_TC <- pairwise.adonis(CMRE5_r8100_t_ja,  TC)
CMRE5_r8100_t_ja_adonis_TC1 <- pairwise.adonis(CMRE5_r8100_t_ja,  TC1)
CMRE5_r8100_t_ja_adonis_TS1 <- pairwise.adonis(CMRE5_r8100_t_ja,  TS1)
CMRE5_r8100_t_ja_adonis_CS <- pairwise.adonis(CMRE5_r8100_t_ja,  CS)

###########################--###########################--###########################--###########################--

#################################-----CAP analysis on each compartment separately----------#####################

############## Mouse 

CMRE5_r8100_mouse <- subset_samples(CMRE5_r8100, Sampletype%in%c("feces"))
CMRE5_r8100_mouse1 = prune_taxa(taxa_sums(CMRE5_r8100_mouse) > 0, CMRE5_r8100_mouse) 

CMRE5_r8100_mouse1_t  = transform_sample_counts(CMRE5_r8100_mouse1, function(x) x / sum(x) ) # transforming sample count
CMRE5_r8100_mouse1_t_br <- distance(CMRE5_r8100_mouse1_t, method = "bray") ##Calculate Bray-curtis 

###########----CAP analysis
CMRE5_r8100_cap <- ordinate(physeq = CMRE5_r8100_mouse1_t, 
                            method = "CAP",
                            distance = CMRE5_r8100_mouse1_t_br,
                            formula = ~ Treatment)
anova(CMRE5_r8100_cap)

###########----CAP plot
CMRE5_r8100_plot_bray_cap_mouse = plot_ordination(CMRE5_r8100_mouse1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_mouse$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_mouse = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")

###########----dispersion using betadisper (PERMDISP)# Bray-Curtis
CMRE5_r8100_t_disp <- as(sample_data(CMRE5_r8100_mouse1_t), "data.frame")
groups <- CMRE5_r8100_t_disp[["Treatment"]]
mod1 <- betadisper(CMRE5_r8100_mouse1_t_br, groups)

############## Water 
CMRE5_r8100_Water <- subset_samples(CMRE5_r8100, Sampletype%in%c("water"))
CMRE5_r8100_Water1 = prune_taxa(taxa_sums(CMRE5_r8100_Water) > 0, CMRE5_r8100_Water) 

CMRE5_r8100_Water1_t  = transform_sample_counts(CMRE5_r8100_Water1, function(x) x / sum(x) ) # transforming sample count CMRE5_r8100_Water1_t_br <- distance(CMRE5_r8100_Water1_t, method = "bray") ##Calculate Bray-curtis 

###########----CAP analysis
CMRE5_r8100_cap <- ordinate(physeq = CMRE5_r8100_Water1_t, 
                            method = "CAP",
                            distance = CMRE5_r8100_Water1_t_br,
                            formula = ~ Treatment)
anova(CMRE5_r8100_cap)

###########----CAP plot
CMRE5_r8100_plot_bray_cap_Water = plot_ordination(CMRE5_r8100_Water1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Water$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Water = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")

###########----dispersion using betadisper (PERMDISP)# Bray-Curtis
CMRE5_r8100_t_disp <- as(sample_data(CMRE5_r8100_Water1_t), "data.frame")
groups <- CMRE5_r8100_t_disp[["Treatment"]]
mod1 <- betadisper(CMRE5_r8100_Water1_t_br, groups)

############## Sediment 
CMRE5_r8100_Sediment <- subset_samples(CMRE5_r8100, Sampletype%in%c("feces"))
CMRE5_r8100_Sediment1 = prune_taxa(taxa_sums(CMRE5_r8100_Sediment) > 0, CMRE5_r8100_Sediment) 

CMRE5_r8100_Sediment1_t  = transform_sample_counts(CMRE5_r8100_Sediment1, function(x) x / sum(x) ) # transforming sample count 
CMRE5_r8100_Sediment1_t_br <- distance(CMRE5_r8100_Sediment1_t, method = "bray") ##Calculate Bray-curtis 

###########----CAP analysis
CMRE5_r8100_cap <- ordinate(physeq = CMRE5_r8100_Sediment1_t, 
                            method = "CAP",
                            distance = CMRE5_r8100_Sediment1_t_br,
                            formula = ~ Treatment)
anova(CMRE5_r8100_cap)

###########----CAP plot
CMRE5_r8100_plot_bray_cap_Sediment = plot_ordination(CMRE5_r8100_Sediment1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Sediment$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Sediment = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")

###########----dispersion using betadisper (PERMDISP)# Bray-Curtis
CMRE5_r8100_t_disp <- as(sample_data(CMRE5_r8100_Sediment1_t), "data.frame")
groups <- CMRE5_r8100_t_disp[["Treatment"]]
mod1 <- betadisper(CMRE5_r8100_Sediment1_t_br, groups)

############## Soil 
CMRE5_r8100_Soil <- subset_samples(CMRE5_r8100, Sampletype%in%c("feces"))
CMRE5_r8100_Soil1 = prune_taxa(taxa_sums(CMRE5_r8100_Soil) > 0, CMRE5_r8100_Soil) 

CMRE5_r8100_Soil1_t  = transform_sample_counts(CMRE5_r8100_Soil1, function(x) x / sum(x) ) # transforming sample count 
CMRE5_r8100_Soil1_t_br <- distance(CMRE5_r8100_Soil1_t, method = "bray") ##Calculate Bray-curtis 

###########----CAP analysis
CMRE5_r8100_cap <- ordinate(physeq = CMRE5_r8100_Soil1_t, 
                            method = "CAP",
                            distance = CMRE5_r8100_Soil1_t_br,
                            formula = ~ Treatment)
anova(CMRE5_r8100_cap)

###########----CAP plot
CMRE5_r8100_plot_bray_cap_Soil = plot_ordination(CMRE5_r8100_Soil1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Soil$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Soil = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")

###########----dispersion using betadisper (PERMDISP)# Bray-Curtis
CMRE5_r8100_t_disp <- as(sample_data(CMRE5_r8100_Soil1_t), "data.frame")
groups <- CMRE5_r8100_t_disp[["Treatment"]]
mod1 <- betadisper(CMRE5_r8100_Soil1_t_br, groups)

############## Root 
CMRE5_r8100_Root <- subset_samples(CMRE5_r8100, Sampletype%in%c("feces"))
CMRE5_r8100_Root1 = prune_taxa(taxa_sums(CMRE5_r8100_Root) > 0, CMRE5_r8100_Root) 

CMRE5_r8100_Root1_t  = transform_sample_counts(CMRE5_r8100_Root1, function(x) x / sum(x) ) # transforming sample count 
CMRE5_r8100_Root1_t_br <- distance(CMRE5_r8100_Root1_t, method = "bray") ##Calculate the Bray-curtis 

###########----CAP analysis
CMRE5_r8100_cap <- ordinate(physeq = CMRE5_r8100_Root1_t, 
                            method = "CAP",
                            distance = CMRE5_r8100_Root1_t_br,
                            formula = ~ Treatment)
anova(CMRE5_r8100_cap)

###########----CAP plot
CMRE5_r8100_plot_bray_cap_Root = plot_ordination(CMRE5_r8100_Root1_t, CMRE5_r8100_cap, color="Treatment")
orddata <- CMRE5_r8100_plot_bray_cap_Root$data
centroids <- aggregate(cbind(CAP1,CAP2)~Treatment,data=orddata,mean)
gg <- merge(orddata,centroids,by="Treatment",suffixes=c("",".centroid"))
gg$Treatment<- factor(gg$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_plot_bray_cap_Root = ggplot(gg, aes(x=CAP1, y=CAP2, color=Treatment)) + theme_set(theme_bw()) + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + geom_point(size=2.00) + geom_point(data=centroid, mapping= aes(x=CAP1, y=CAP2, size=0.50)) + theme(legend.text=element_text(size=12), legend.title=element_text(size=14))  + scale_color_manual(values=c("gray50", "#E69F00", "#009E73", "#D55E00")) + theme(plot.title = element_text(size=20, hjust = 0.5))  + scale_shape_discrete(solid=T)  + theme(plot.title = element_text(size=15, hjust = 0.5))  + geom_segment(aes(x=CAP1.centroid, y=CAP2.centroid, xend=CAP1, yend=CAP2, color=Treatment), alpha=0.4) + theme(legend.position = "none")

###########----dispersion using betadisper (PERMDISP)# Bray-Curtis
CMRE5_r8100_t_disp <- as(sample_data(CMRE5_r8100_Root1_t), "data.frame")
groups <- CMRE5_r8100_t_disp[["Treatment"]]
mod1 <- betadisper(CMRE5_r8100_Root1_t_br, groups)


###############---Now make combined figure for all beta cap compartments-------##############
A <- ggarrange(CMRE5_r8100_plot_bray_cap_Water, CMRE5_r8100_plot_bray_cap_Sediment, ncol = 3) 
B <- ggarrange(CMRE5_r8100_plot_bray_cap_Soil, CMRE5_r8100_plot_bray_cap_Root, CMRE5_r8100_plot_bray_cap_mouse, ncol = 3)
ggarrange (A, B, ncol = 1)

########################################################################################################################
