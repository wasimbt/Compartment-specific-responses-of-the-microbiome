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
# library("DESeq2"); packageVersion("DESeq2")
# [1] '1.22.2'
# library("RColorBrewer"); packageVersion("RColorBrewer")
# ##[1] ‘1.1.3’
# library("scales"); packageVersion("scales")
# ##[1] ‘1.2.0’
# library("ggpubr"); packageVersion("ggpubr")
# ##[1] ‘0.4.0’

#####----------Loading and creating the object---------###########

CMRE5_High = prune_samples(sample_sums(CMRE5)>=8100, CMRE5)
CMRE5_High_p = prune_taxa(taxa_sums(CMRE5_High) > 0, CMRE5_High) # Remove taxa with 0 counts
sample_data(CMRE5_High_p)$Treatment<- factor(sample_data(CMRE5_High_p)$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

##########################---Mouse Diff Abundance Analysis-----#################

#####Subset dataset for mouse

CMRE5_High_p_mouse <- subset_samples(CMRE5_High_p, Sampletype%in%c("feces"))
CMRE5_High_mouse = prune_taxa(taxa_sums(CMRE5_High_p_mouse) > 0, CMRE5_High_p_mouse) #Remove taxa with 0 counts

#########################-------- Mouse Ctr-AS --------##########################

CMRE5_High_mouse_CtrAs <- subset_samples(CMRE5_High_mouse, Treatment %in%c("Ctr", "As"))
CMRE5_High_mouse_CtrAs1 = prune_taxa(taxa_sums(CMRE5_High_mouse_CtrAs) > 0, CMRE5_High_mouse_CtrAs) #Remove taxa with 0 counts
###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Mouse_CtrAs1 = phyloseq_to_deseq2(CMRE5_High_mouse_CtrAs1, ~ Treatment) #
diagdds_Mouse_CtrAs2 = DESeq(diagdds_Mouse_CtrAs1, test="Wald", fitType="parametric") 
res = results(diagdds_Mouse_CtrAs2, cooksCutoff = FALSE)
CMRE5_High_mouse_CtrAs1_res = res[order(res$padj, na.last=NA), ]

##########################-------- Mouse Ctr-Bx --------##########################

CMRE5_High_mouse_CtrBx <- subset_samples(CMRE5_High_mouse, Treatment %in%c("Ctr", "Bx"))
CMRE5_High_mouse_CtrBx1 = prune_taxa(taxa_sums(CMRE5_High_mouse_CtrBx) > 0, CMRE5_High_mouse_CtrBx) #Remove taxa with 0 counts

###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Mouse_CtrBx1 = phyloseq_to_deseq2(CMRE5_High_mouse_CtrBx1, ~ Treatment) #
diagdds_Mouse_CtrBx2 = DESeq(diagdds_Mouse_CtrBx1, test="Wald", fitType="parametric") ##
res = results(diagdds_Mouse_CtrBx2, cooksCutoff = FALSE)
CMRE5_High_mouse_CtrBx1_res = res[order(res$padj, na.last=NA), ]

##########################-------- Mouse Ctr-Tb --------##########################

CMRE5_High_mouse_CtrTb <- subset_samples(CMRE5_High_mouse, Treatment %in%c("Ctr", "Tb"))
CMRE5_High_mouse_CtrTb1 = prune_taxa(taxa_sums(CMRE5_High_mouse_CtrTb) > 0, CMRE5_High_mouse_CtrTb) #Remove taxa with 0 counts

###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Mouse_CtrTb1 = phyloseq_to_deseq2(CMRE5_High_mouse_CtrTb1, ~ Treatment) #
diagdds_Mouse_CtrTb2 = DESeq(diagdds_Mouse_CtrTb1, test="Wald", fitType="parametric") ##
res = results(diagdds_Mouse_CtrTb2, cooksCutoff = FALSE)
CMRE5_High_mouse_CtrTb2_res = res[order(res$padj, na.last=NA), ]

###############################################################################################

##########################---Root Diff Abundance Analysis-----#################

#####Subset dataset for Root

CMRE5_High_p_Root <- subset_samples(CMRE5_High_p, Sampletype%in%c("roots"))
CMRE5_High_Root = prune_taxa(taxa_sums(CMRE5_High_p_Root) > 0, CMRE5_High_p_Root) #Remove taxa with 0 counts
sample_data(CMRE5_High_Root)$Treatment <- factor(sample_data(CMRE5_High_Root)$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

##########################-------- Root Ctr-AS --------##########################

CMRE5_High_Root_CtrAs <- subset_samples(CMRE5_High_Root, Treatment %in%c("Ctr", "As"))
CMRE5_High_Root_CtrAs1 = prune_taxa(taxa_sums(CMRE5_High_Root_CtrAs) > 0, CMRE5_High_Root_CtrAs) #Remove taxa with 0 counts
###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Root_CtrAs1 = phyloseq_to_deseq2(CMRE5_High_Root_CtrAs1, ~ Treatment) #
diagdds_Root_CtrAs2 = DESeq(diagdds_Root_CtrAs1, test="Wald", fitType="parametric") ##
res = results(diagdds_Root_CtrAs2, cooksCutoff = FALSE)
CMRE5_High_Root_CtrAs2_res = res[order(res$padj, na.last=NA), ]

##########################-------- Root Ctr-Bx --------##########################

CMRE5_High_Root_CtrBx <- subset_samples(CMRE5_High_Root, Treatment %in%c("Ctr", "Bx"))
CMRE5_High_Root_CtrBx1 = prune_taxa(taxa_sums(CMRE5_High_Root_CtrBx) > 0, CMRE5_High_Root_CtrBx) #Remove taxa with 0 counts

###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Root_CtrBx1 = phyloseq_to_deseq2(CMRE5_High_Root_CtrBx1, ~ Treatment) #
diagdds_Root_CtrBx2 = DESeq(diagdds_Root_CtrBx1, test="Wald", fitType="parametric") ##
res = results(diagdds_Root_CtrBx2, cooksCutoff = FALSE)
CMRE5_High_Root_CtrBx2_res = res[order(res$padj, na.last=NA), ]

##########################-------- Root Ctr-Tb --------##########################

CMRE5_High_Root_CtrTb <- subset_samples(CMRE5_High_Root, Treatment %in%c("Ctr", "Tb"))
CMRE5_High_Root_CtrTb1 = prune_taxa(taxa_sums(CMRE5_High_Root_CtrTb) > 0, CMRE5_High_Root_CtrTb) #Remove taxa with 0 counts

###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Root_CtrTb1 = phyloseq_to_deseq2(CMRE5_High_Root_CtrTb1, ~ Treatment) #
diagdds_Root_CtrTb2 = DESeq(diagdds_Root_CtrTb1, test="Wald", fitType="parametric") ##
res = results(diagdds_Root_CtrTb2, cooksCutoff = FALSE)
CMRE5_High_Root_CtrTb2_res = res[order(res$padj, na.last=NA), ]

###############################################################################################

##########################---Soil Diff Abundance Analysis-----#################

#####Subset dataset for Soil

CMRE5_High_p_Soil <- subset_samples(CMRE5_High_p, Sampletype%in%c("soil"))
CMRE5_High_Soil = prune_taxa(taxa_sums(CMRE5_High_p_Soil) > 0, CMRE5_High_p_Soil) #Remove taxa with 0 counts
sample_data(CMRE5_High_Soil)$Treatment<- factor(sample_data(CMRE5_High_Soil)$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

##########################-------- Soil Ctr-AS --------##########################

CMRE5_High_Soil_CtrAs <- subset_samples(CMRE5_High_Soil, Treatment %in%c("Ctr", "As"))
CMRE5_High_Soil_CtrAs1 = prune_taxa(taxa_sums(CMRE5_High_Soil_CtrAs) > 0, CMRE5_High_Soil_CtrAs) #Remove taxa with 0 counts

###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Soil_CtrAs1 = phyloseq_to_deseq2(CMRE5_High_Soil_CtrAs1, ~ Treatment) #
diagdds_Soil_CtrAs2 = DESeq(diagdds_Soil_CtrAs1, test="Wald", fitType="parametric") ## 
res = results(diagdds_Soil_CtrAs2, cooksCutoff = FALSE)
CMRE5_High_Soil_CtrAs2_res = res[order(res$padj, na.last=NA), ]

##########################-------- Soil Ctr-Bx --------##########################

CMRE5_High_Soil_CtrBx <- subset_samples(CMRE5_High_Soil, Treatment %in%c("Ctr", "Bx"))
CMRE5_High_Soil_CtrBx1 = prune_taxa(taxa_sums(CMRE5_High_Soil_CtrBx) > 0, CMRE5_High_Soil_CtrBx) #Remove taxa with 0 counts

###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Soil_CtrBx1 = phyloseq_to_deseq2(CMRE5_High_Soil_CtrBx1, ~ Treatment) #
diagdds_Soil_CtrBx2 = DESeq(diagdds_Soil_CtrBx1, test="Wald", fitType="parametric") ##
res = results(diagdds_Soil_CtrBx2, cooksCutoff = FALSE)
CMRE5_High_Soil_CtrBx2_res = res[order(res$padj, na.last=NA), ]

##########################-------- Soil Ctr-Tb --------##########################

CMRE5_High_Soil_CtrTb <- subset_samples(CMRE5_High_Soil, Treatment %in%c("Ctr", "Tb"))
CMRE5_High_Soil_CtrTb1 = prune_taxa(taxa_sums(CMRE5_High_Soil_CtrTb) > 0, CMRE5_High_Soil_CtrTb) #Remove taxa with 0 counts

###########################-----Perform the DeSEQ2 analysis--------####################

diagdds_Soil_CtrTb1 = phyloseq_to_deseq2(CMRE5_High_Soil_CtrTb1, ~ Treatment) #
diagdds_Soil_CtrTb2 = DESeq(diagdds_Soil_CtrTb1, test="Wald", fitType="parametric") ##
res = results(diagdds_Soil_CtrTb2, cooksCutoff = FALSE)
CMRE5_High_Soil_CtrTb2_res = res[order(res$padj, na.last=NA), ]

###############################################################################################

##########################---Sediment Diff Abundance Analysis-----#################
 
####Subset dataset for Sediment
 
 CMRE5_High_p_Sediment <- subset_samples(CMRE5_High_p, Sampletype%in%c("sedi"))
 CMRE5_High_Sediment = prune_taxa(taxa_sums(CMRE5_High_p_Sediment) > 0, CMRE5_High_p_Sediment) #Remove taxa with 0 counts
 sample_data(CMRE5_High_Sediment)$Treatment <- factor(sample_data(CMRE5_High_Sediment)$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))
 
 ##########################-------- Sediment Ctr-AS --------##########################
 
 CMRE5_High_Sediment_CtrAs <- subset_samples(CMRE5_High_Sediment, Treatment %in%c("Ctr", "As"))
 CMRE5_High_Sediment_CtrAs1 = prune_taxa(taxa_sums(CMRE5_High_Sediment_CtrAs) > 0, CMRE5_High_Sediment_CtrAs) #Remove taxa with 0 counts
 
 ###########################-----Perform the DeSEQ2 analysis--------####################
 
 diagdds_Sediment_CtrAs1 = phyloseq_to_deseq2(CMRE5_High_Sediment_CtrAs1, ~ Treatment) #
 diagdds_Sediment_CtrAs2 = DESeq(diagdds_Sediment_CtrAs1, test="Wald", fitType="parametric") 
 res = results(diagdds_Sediment_CtrAs2, cooksCutoff = FALSE)
 CMRE5_High_Sediment_CtrAs2_res = res[order(res$padj, na.last=NA), ]
 
 ##########################-------- Sediment Ctr-Bx --------##########################
 
 CMRE5_High_Sediment_CtrBx <- subset_samples(CMRE5_High_Sediment, Treatment %in%c("Ctr", "Bx"))
 CMRE5_High_Sediment_CtrBx1 = prune_taxa(taxa_sums(CMRE5_High_Sediment_CtrBx) > 0, CMRE5_High_Sediment_CtrBx) #Remove taxa with 0 counts
 
 ###########################-----Perform the DeSEQ2 analysis--------####################
 
 diagdds_Sediment_CtrBx1 = phyloseq_to_deseq2(CMRE5_High_Sediment_CtrBx1, ~ Treatment) #
 diagdds_Sediment_CtrBx2 = DESeq(diagdds_Sediment_CtrBx1, test="Wald", fitType="parametric") 
 res = results(diagdds_Sediment_CtrBx2, cooksCutoff = FALSE)
 CMRE5_High_Sediment_CtrBx2_res = res[order(res$padj, na.last=NA), ]
 
 ##########################-------- Sediment Ctr-Tb --------##########################
 
 CMRE5_High_Sediment_CtrTb <- subset_samples(CMRE5_High_Sediment, Treatment %in%c("Ctr", "Tb"))
 CMRE5_High_Sediment_CtrTb1 = prune_taxa(taxa_sums(CMRE5_High_Sediment_CtrTb) > 0, CMRE5_High_Sediment_CtrTb) #Remove taxa with 0 counts
 
 ###########################-----Perform the DeSEQ2 analysis--------####################
 
 diagdds_Sediment_CtrTb1 = phyloseq_to_deseq2(CMRE5_High_Sediment_CtrTb1, ~ Treatment) #
 diagdds_Sediment_CtrTb2 = DESeq(diagdds_Sediment_CtrTb1, test="Wald", fitType="parametric") 
 res = results(diagdds_Sediment_CtrTb2, cooksCutoff = FALSE)
 CMRE5_High_Sediment_CtrTb2_res = res[order(res$padj, na.last=NA), ]

 ###############################################################################################
 
 ##########################---Water Diff Abundance Analysis-----#################
 
 #####Subset dataset for Water
 
 CMRE5_High_p_Water <- subset_samples(CMRE5_High_p, Sampletype%in%c("water"))
 CMRE5_High_Water = prune_taxa(taxa_sums(CMRE5_High_p_Water) > 0, CMRE5_High_p_Water) #Remove taxa with 0 counts
 sample_data(CMRE5_High_Water)$Treatment <- factor(sample_data(CMRE5_High_Water)$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))
 
 ##########################-------- Water Ctr-AS --------##########################
 
 CMRE5_High_Water_CtrAs <- subset_samples(CMRE5_High_Water, Treatment %in%c("Ctr", "As"))
 CMRE5_High_Water_CtrAs1 = prune_taxa(taxa_sums(CMRE5_High_Water_CtrAs) > 0, CMRE5_High_Water_CtrAs) #Remove taxa with 0 counts
 
 ###########################-----Perform the DeSEQ2 analysis--------####################
 
 diagdds_Water_CtrAs1 = phyloseq_to_deseq2(CMRE5_High_Water_CtrAs1, ~ Treatment) #
 diagdds_Water_CtrAs2 = DESeq(diagdds_Water_CtrAs1, test="Wald", fitType="parametric") ## 
 res = results(diagdds_Water_CtrAs2, cooksCutoff = FALSE)
 CMRE5_High_Water_CtrAs2_res = res[order(res$padj, na.last=NA), ]
 
 ##########################-------- Water Ctr-Bx --------##########################
 
 CMRE5_High_Water_CtrBx <- subset_samples(CMRE5_High_Water, Treatment %in%c("Ctr", "Bx"))
 CMRE5_High_Water_CtrBx1 = prune_taxa(taxa_sums(CMRE5_High_Water_CtrBx) > 0, CMRE5_High_Water_CtrBx) #Remove taxa with 0 counts
 
 ###########################-----Perform the DeSEQ2 analysis--------####################
 
 diagdds_Water_CtrBx1 = phyloseq_to_deseq2(CMRE5_High_Water_CtrBx1, ~ Treatment) #
 diagdds_Water_CtrBx2 = DESeq(diagdds_Water_CtrBx1, test="Wald", fitType="parametric") ##
 res = results(diagdds_Water_CtrBx2, cooksCutoff = FALSE)
 CMRE5_High_Water_CtrBx2_res = res[order(res$padj, na.last=NA), ]
 
##########################-------- Water Ctr-Tb --------##########################
 
 CMRE5_High_Water_CtrTb <- subset_samples(CMRE5_High_Water, Treatment %in%c("Ctr", "Tb"))
 CMRE5_High_Water_CtrTb1 = prune_taxa(taxa_sums(CMRE5_High_Water_CtrTb) > 0, CMRE5_High_Water_CtrTb) #Remove taxa with 0 counts
 
 ###########################-----Perform the DeSEQ2 analysis--------####################
 
 diagdds_Water_CtrTb1 = phyloseq_to_deseq2(CMRE5_High_Water_CtrTb1, ~ Treatment) #
 diagdds_Water_CtrTb2 = DESeq(diagdds_Water_CtrTb1, test="Wald", fitType="parametric") ##
 res = results(diagdds_Water_CtrTb2, cooksCutoff = FALSE)
 CMRE5_High_Water_CtrTb2_res = res[order(res$padj, na.last=NA), ]
 
###############################################################################################
 
#####################################---- Make MA plot for all-----################

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

##########################################################################################################################

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

##########################################################################################################################

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

##########################################################################################################################

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

##########################################################################################################################

#####################----Combined plot for all above----##################

figure <- ggarrange(New1_ASV_relative_abundance_Water_CtrAs_MAplot3, 
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




