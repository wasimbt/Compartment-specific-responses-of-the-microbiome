##########################################------Supplementary Figure S2---------###############################################

###################--Chemical validation experiment plot
---
title: "CMRE chemical data Figure v2"
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

```{r,include=FALSE}
# install.packages("sciplot")
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("pander")
# install.packages("knitr")
# install.packages("readr")
# install.packages("readxl")
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("multcomp")
# install.packages("ggpubr")
# install.packages("ggpubr")
# install.packages("ggpmisc")
# install.packages("phyloseq")
# install.packages("gplots") # to convert colour names to #hex identities
```

```{r libraries, echo=F, message=F, warning=F}
library(pander)
library(readr)
library(readxl)
library(ggpubr)
library(cowplot)
library(gplots) # to convert colour names to #hex identities
```

# colors

```{r}
col_assignments <- data.frame(color=c("grey50", "#E69F00", "#009E73", "#D55E00")) 
rownames(col_assignments) <- c("Ctr", "As", "Bx", "Tb")
# barplot(rep(1,4), col=col_assignments$color, names=col_assignments$treatment, border=NA)
plot(rep(1,4), col=col_assignments$color, cex=10, pch=19, ylim=c(.75,1.25), xlim=c(0,5), axes=F, xlab="",ylab="")
text(1:4, 1, labels=rownames(col_assignments), cex=2)
```



# treatment solution and lake water data

```{r,include=FALSE, echo=F, message=F, warning=F}
# load data
data <- as.data.frame(readxl::read_excel("input/treatment_solutions_lake_water_chem_results.xlsx", sheet="allThree") )

data$Experiment <- as.factor(data$Experiment)
data$Experiment <- factor(data$Experiment, levels(data$Experiment) [c(2,3,1)] )
levels(data$Experiment) <- c("Water/sediment", "Soil/plant", "Animal")
data$Treatment <- as.factor(data$Treatment)
# data$Treatment <- factor(data$Treatment, levels(data$Treatment) [c(1,4,3,2)] )
data$Values <- as.numeric(data$`Final Conc. (μg/L)`)


# data1 <- droplevels(data[data$Analysis=="As",])
# 
# ggplot(data1, aes(x=Treatment, y=Values)) + #ylim(-1,2000) +
#   scale_y_continuous(trans='pseudo_log', breaks=c(1,10,100,1000, 10000, 100000)) +
#   theme_bw() +
#   facet_wrap(~Experiment, nrow = 1) +  
#   geom_boxplot(show.legend=FALSE, outlier.shape=NA) +
#   scale_fill_manual(values=c("grey90", "grey90", "grey90")) +
#   geom_point(aes(colour=Treatment, size=2), position=position_jitterdodge(), alpha=0.75, show.legend=T) +
#   scale_size(guide=F) +
#   scale_color_manual(values=rep(col_assignments["As",],4), guide=F) +
#   labs(x=" ", y="concentration [µg/L]", title="As in treatment solutions")



chemical_plots1 <- function(CHEMICAL){
    DATA <- droplevels( data[data$Analysis==paste(CHEMICAL), ] )
    ggplot(DATA, aes(x=Treatment, y=Values) ) +
    scale_y_continuous(trans='pseudo_log', breaks=c(0, 10, 30, 100, 300, 1000, 3000)) +
    theme_bw() +
    facet_wrap(~Experiment, nrow = 1) +  
#    geom_boxplot(show.legend=FALSE, outlier.shape=NA) +
    scale_fill_manual(values=c("grey90", "grey90", "grey90")) +
    geom_point(aes(colour=Treatment, size=2), position=position_jitterdodge(), alpha=0.75, show.legend=F) +
    scale_size(guide=F) + 
    scale_color_manual(values=rep(col_assignments[CHEMICAL,],4), guide=F) +
    labs(x="Treatment solutions", y="concentration [µg/L]", title=paste(CHEMICAL) )
}
```

```{r, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=5}
chemical_plots1("As")
```

```{r, eval=T}
# pdf("Treatment_solutions_As.pdf", width=18/cm(1), height=9/cm(1), pointsize=2, fonts="Helvetica")
# chemical_plots1("As")
# dev.off()
```

```{r, eval=T}
chemical_plots1("Bx")

# pdf("Treatment_solutions_Bx.pdf", width=18/cm(1), height=9/cm(1), pointsize=2, fonts="Helvetica")
# chemical_plots1("Bx")
# dev.off()
```

```{r, eval=T}
chemical_plots1("Tb")

# pdf("Treatment_solutions_Tb.pdf", width=18/cm(1), height=9/cm(1), pointsize=2, fonts="Helvetica")
# chemical_plots1("Tb")
# dev.off()
```



# combined Figure 1

```{r}

chemical_plots_combined <- plot_grid(
  chemical_plots1("As"),
  chemical_plots1("Bx"),
  chemical_plots1("Tb"),
  align="hv", nrow=3, ncol=1, rel_heights=c(10, 10, 10) )

print(chemical_plots_combined)

pdf("CMRE_chemical_data_Figure_combined.pdf", width=20/cm(1), height=30/cm(1), pointsize=10, fonts="Helvetica")
print(chemical_plots_combined)
noshow <- dev.off()

```

##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure S4---------###############################################

###################--NMDS plot

levels(sample_data(CMRE5_r8100)$Sampletype) <- c("Animal", "Plant", "Sediment", "Soil","Water")
sample_data(CMRE5_r8100)$Sampletype <- factor(sample_data(CMRE5_r8100)$Sampletype, 
                                             levels = c("Water", "Sediment", "Soil", "Plant", "Animal")) # ordering samples
set.seed(100)
CMRE5_r8100_ordi <- phyloseq::ordinate(CMRE5_r8100, method="NMDS", distance="bray")
p0 <- phyloseq::plot_ordination(CMRE5_r8100, CMRE5_r8100_ordi, color="Sampletype") +
      theme_bw() +
      geom_point(size=2) + #, alpha=0.75) +
      scale_color_manual(labels = c("Water", "Sediment", "Soil", "Plant", "Animal"), 
                         values=c("steelblue3", "slategray4",  "salmon4", "darkolivegreen4", "tan3"), name="Food chain component") 
pdf("Fig_S4_betaDiversity_NMDS_v20240111.pdf", width=15/cm(1), height=10/cm(1), pointsize=2, fonts="Helvetica")
print(p0)
dev.off()
##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure S5---------#############################################

###################--Shannon final box plot--for treatement and compartment 

CMRE5_alpha_r8100$Sampletype<- factor(CMRE5_alpha_r8100$Sampletype, levels = c("Water", "Sediment", "Soil", "Root", "Mouse"))
CMRE5_alpha_r8100$Treatment<- factor(CMRE5_alpha_r8100$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

######-- Fig. S5A --####

FigS5A_data <- CMRE5_alpha_r8100 %>% 
  select(Treatment,Sampletype,Shannon,
         Concentration,Lib)

FigS5A_data_CTR <- FigS5A_data %>% filter(Treatment=="Ctr")

FigS5A_data_Obs_means <- FigS5A_data_CTR %>% 
  group_by(Sampletype) %>% 
  summarise(MeanShannon=mean(Shannon),
            QuantDShannon=quantile(Shannon,probs = 0.95))


FigS5A_data_anova <- aov(Shannon ~ Treatment * Sampletype + Lib, data = FigS5A_data)

FigS5A_tukey <- TukeyHSD(FigS5A_data_anova)# Tukey's test
#The use of letters to indicate significant differences in pairwise comparisons is called compact letter display, and can simplify the visualisation and discussion of significant differences among means. We are going to use the multcompLetters4 function from the multcompView package. The arguments are the object from an aov function and the object from the TukeyHSD function.
# compact letter display
FigS5A_cld <-as.data.frame.list(multcompView::multcompLetters4(FigS5A_data_anova, FigS5A_tukey)[[4]])

FigS5A_cld_Ctr <- FigS5A_cld %>% 
  filter(grepl("Ctr", rownames(FigS5A_cld))) #select only the controls
rownames(FigS5A_cld_Ctr) <- gsub("Ctr:","",rownames(FigS5A_cld_Ctr) )

FigS5A_cld_Tbl <- data.frame(Sampletype =rownames(FigS5A_cld_Ctr),
                            L=FigS5A_cld_Ctr$Letters)

FigS5A_data_Obs_means_letters <- merge(FigS5A_data_Obs_means,FigS5A_cld_Tbl)

plot_FigS5A <- ggbarplot(FigS5A_data %>% filter(Treatment=="Ctr"), 
                        x = "Sampletype", 
                        y = "Shannon", 
                        fill = "grey50",
                        width=0.2,
                        add = c("mean_se", "jitter"),
                        legend = "none",
                        ylim=c(0,8)
)      +
  theme_bw() +
  geom_text(data = FigS5A_data_Obs_means_letters, 
            aes(x = Sampletype, y = QuantDShannon, label = L), 
            size = 5,  vjust=-0.8, hjust =0.5)    +
  ylab("Shannon")+  xlab("") + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=14, face = "bold"), axis.text.y = element_text(size=12, face = "bold"))  + theme(plot.margin = unit(c(1,1,-1.5,0), "lines"))

plot_FigS5A

###########----Fig. S5B --############################################
#removing the mean for each group
FigS5B_data_noC <- CMRE5_alpha_r8100 %>% filter(Treatment!="Ctr") %>% select(Sampletype,Shannon,Treatment)
#removing the means
FigS5B_data_noC_diff <- merge(FigS5B_data_noC,FigS5A_data_Obs_means[,1:2])
FigS5B_data_noC_diff$Diff <- FigS5B_data_noC_diff$Shannon-FigS5B_data_noC_diff$MeanShannon

#with letters
FigS5B_cld <- FigS5A_cld %>% filter(!grepl("Ctr", rownames(FigS5A_cld))) #select all, but the controls

#Compact letter display to indicate significant differences
FigS5B_cld_Tbl <- data.frame(Sampletype =rownames(FigS5B_cld),L=FigS5B_cld$Letters)
FigS5B_cld_Tbl1 <- separate(FigS5B_cld_Tbl,col="Sampletype",sep = ":",into=c("Treatment","Sampletype"))
FigS5B_data_Obs_means <- aggregate(Shannon ~  Sampletype*Treatment, CMRE5_alpha_r8100, mean)
FigS5B_data_Obs_sd <- aggregate(Shannon ~  Sampletype*Treatment, CMRE5_alpha_r8100, sd)
colnames(FigS5B_data_Obs_sd)[3] <- "sd"
#merge
FigS5B_data_mergeX.SD <- merge(FigS5B_data_Obs_means,FigS5B_data_Obs_sd,
                              by=c("Sampletype","Treatment"))
FigS5B_data_mergeX.SD.L <- merge(FigS5B_data_mergeX.SD,
                                FigS5B_cld_Tbl1,
                                by=c("Sampletype","Treatment"))

#merging with control group mean 
FigS5B_data_mergeX.SD.L.C <- merge(FigS5B_data_mergeX.SD.L,FigS5A_data_Obs_means[,1:2],by="Sampletype")

rename(FigS5B_data_mergeX.SD.L.C,MeanControl= MeanShannon )
FigS5B_data_mergeX.SD.L.C <- FigS5B_data_mergeX.SD.L.C %>% mutate(Diff=Shannon-MeanShannon)
#removing the control
FigS5B_data_mergeX.SD.L.C_noC <- FigS5B_data_mergeX.SD.L.C %>% filter(Treatment!="Ctr") 

### rect  
FigS5B_data_mergeX.SD.L.C_noC_Sort <- FigS5B_data_mergeX.SD.L.C_noC %>% arrange(Sampletype,Treatment) %>% 
  mutate(MeanDiff=MeanShannon+Diff)

FigS5B_data_mergeX.SD.L.C_noC_Sort$N <- 1:nrow(FigS5B_data_mergeX.SD.L.C_noC_Sort)
FigS5B_data_mergeX.SD.L.C_noC_Sort

# AS"#E69F00", 
# BX"#009E73",  
# Tb"#D55E00" 
plot_FigS5B_prep <- ggplot(FigS5B_data_mergeX.SD.L.C_noC_Sort,
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

# plot_FigS5B_prep

NewY=FigS5B_data_mergeX.SD.L.C_noC_Sort$MeanDiff+ #adjusting the high of the labels
  c(0,0,-0.05,
    0,0,0,
    0.1,0.10,0.1,
    -0.1,-0.1,-0.1,
    -0.1,-0.1,-0.1
  )
plot_FigS5B <- plot_FigS5B_prep +
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
  theme_bw() + theme(panel.grid.minor.y = element_line(linewidth = 0.25, linetype = 1), legend.title = element_text(size=18), legend.text = element_text(size=16)) + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=14, face = "bold"), axis.text.y = element_text(size=12, face = "bold")) + theme(plot.margin = unit(c(1,1,-1.5,0), "lines"))


plot_FigS5B

# Combination of plots
require(patchwork)
plot_FigS5A + plot_FigS5B

ggplot2::ggsave("Fig.S5_Suppl_Shannon.tiff", 
                width = 11.80, height = 6.30, dpi=300)


##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure S6---------###############################################

###################--Beta Diversity-Box plots 

CMRE5_r8100_t  = readRDS("CMRE5_r8100_t.RDS")

CMRE5_r8100_t_feces <- subset_samples(CMRE5_r8100_t, Sampletype%in%c("feces"))
CMRE5_r8100_t_roots <- subset_samples(CMRE5_r8100_t, Sampletype%in%c("roots"))
CMRE5_r8100_t_soil <- subset_samples(CMRE5_r8100_t, Sampletype%in%c("soil"))
CMRE5_r8100_t_water <- subset_samples(CMRE5_r8100_t, Sampletype%in%c("water"))
CMRE5_r8100_t_sedi <- subset_samples(CMRE5_r8100_t, Sampletype%in%c("sedi"))

###############################--------- Bray-Curtis------------------####################

#######--For mouse -----#######

CMRE5_r8100_t_feces_As <- subset_samples(CMRE5_r8100_t_feces, Treatment%in%c("As"))
CMRE5_r8100_t_feces_Bx <- subset_samples(CMRE5_r8100_t_feces, Treatment%in%c("Bx"))
CMRE5_r8100_t_feces_Tb <- subset_samples(CMRE5_r8100_t_feces, Treatment%in%c("Tb"))
CMRE5_r8100_t_feces_Ctr <- subset_samples(CMRE5_r8100_t_feces, Treatment%in%c("Ctr"))

CMRE5_r8100_t_feces_As_br <- distance(CMRE5_r8100_t_feces_As, method = "bray") #bray-curtis 
CMRE5_r8100_t_feces_Bx_br <- distance(CMRE5_r8100_t_feces_Bx, method = "bray")
CMRE5_r8100_t_feces_Tb_br <- distance(CMRE5_r8100_t_feces_Tb, method = "bray")
CMRE5_r8100_t_feces_Ctr_br <- distance(CMRE5_r8100_t_feces_Ctr, method = "bray")

CMRE5_r8100_t_feces_As_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_feces_As_br))) ###converting distance matrix to list 
CMRE5_r8100_feces_As_br_pairlist <- CMRE5_r8100_t_feces_As_br_pairlist[as.numeric(CMRE5_r8100_t_feces_As_br_pairlist$col) > as.numeric(CMRE5_r8100_t_feces_As_br_pairlist$row),]
CMRE5_r8100_feces_As_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE5_r8100_feces_As_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE5_r8100_feces_As_br_pairlist$Treatment <- rep("As",nrow(CMRE5_r8100_feces_As_br_pairlist)) #make new column Sampletype and fill all Mouse

##

CMRE5_r8100_t_feces_Bx_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_feces_Bx_br))) ###converting distance matrix to list 
CMRE5_r8100_feces_Bx_br_pairlist <- CMRE5_r8100_t_feces_Bx_br_pairlist[as.numeric(CMRE5_r8100_t_feces_Bx_br_pairlist$col) > as.numeric(CMRE5_r8100_t_feces_Bx_br_pairlist$row),]
CMRE5_r8100_feces_Bx_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE5_r8100_feces_Bx_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE5_r8100_feces_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE5_r8100_feces_Bx_br_pairlist)) #make new column Sampletype and fill all Mouse

##

CMRE5_r8100_t_feces_Tb_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_feces_Tb_br))) ###converting distance matrix to list 
CMRE5_r8100_feces_Tb_br_pairlist <- CMRE5_r8100_t_feces_Tb_br_pairlist[as.numeric(CMRE5_r8100_t_feces_Tb_br_pairlist$col) > as.numeric(CMRE5_r8100_t_feces_Tb_br_pairlist$row),]
CMRE5_r8100_feces_Tb_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE5_r8100_feces_Tb_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE5_r8100_feces_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE5_r8100_feces_Tb_br_pairlist)) #make new column Sampletype and fill all Mouse

##

CMRE5_r8100_t_feces_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_feces_Ctr_br))) ###converting distance matrix to list 
CMRE5_r8100_feces_Ctr_br_pairlist <- CMRE5_r8100_t_feces_Ctr_br_pairlist[as.numeric(CMRE5_r8100_t_feces_Ctr_br_pairlist$col) > as.numeric(CMRE5_r8100_t_feces_Ctr_br_pairlist$row),]
CMRE5_r8100_feces_Ctr_br_pairlist$Sampletype <- rep("Mouse",nrow(CMRE5_r8100_feces_Ctr_br_pairlist)) #make new column Sampletype and fill all Mouse
CMRE5_r8100_feces_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE5_r8100_feces_Ctr_br_pairlist)) #make new column Sampletype and fill all Mouse

##################now merge all treatment csv files

CMRE5_r8100_feces_pairlist <- rbind(CMRE5_r8100_feces_Ctr_br_pairlist, CMRE5_r8100_feces_As_br_pairlist, CMRE5_r8100_feces_Bx_br_pairlist, CMRE5_r8100_feces_Tb_br_pairlist)

######################################

#######--For roots -----####### 

CMRE5_r8100_t_roots_As <- subset_samples(CMRE5_r8100_t_roots, Treatment%in%c("As"))
CMRE5_r8100_t_roots_Bx <- subset_samples(CMRE5_r8100_t_roots, Treatment%in%c("Bx"))
CMRE5_r8100_t_roots_Tb <- subset_samples(CMRE5_r8100_t_roots, Treatment%in%c("Tb"))
CMRE5_r8100_t_roots_Ctr <- subset_samples(CMRE5_r8100_t_roots, Treatment%in%c("Ctr"))

CMRE5_r8100_t_roots_As_br <- distance(CMRE5_r8100_t_roots_As, method = "bray") #bray-curtis 
CMRE5_r8100_t_roots_Bx_br <- distance(CMRE5_r8100_t_roots_Bx, method = "bray")
CMRE5_r8100_t_roots_Tb_br <- distance(CMRE5_r8100_t_roots_Tb, method = "bray")
CMRE5_r8100_t_roots_Ctr_br <- distance(CMRE5_r8100_t_roots_Ctr, method = "bray")

CMRE5_r8100_t_roots_As_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_roots_As_br))) ###converting distance matrix to list 
CMRE5_r8100_roots_As_br_pairlist <- CMRE5_r8100_t_roots_As_br_pairlist[as.numeric(CMRE5_r8100_t_roots_As_br_pairlist$col) > as.numeric(CMRE5_r8100_t_roots_As_br_pairlist$row),]
CMRE5_r8100_roots_As_br_pairlist$Sampletype <- rep("Root",nrow(CMRE5_r8100_roots_As_br_pairlist)) #make new column Sampletype and fill all Root
CMRE5_r8100_roots_As_br_pairlist$Treatment <- rep("As",nrow(CMRE5_r8100_roots_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE5_r8100_t_roots_Bx_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_roots_Bx_br))) ###converting distance matrix to list 
CMRE5_r8100_roots_Bx_br_pairlist <- CMRE5_r8100_t_roots_Bx_br_pairlist[as.numeric(CMRE5_r8100_t_roots_Bx_br_pairlist$col) > as.numeric(CMRE5_r8100_t_roots_Bx_br_pairlist$row),]
CMRE5_r8100_roots_Bx_br_pairlist$Sampletype <- rep("Root",nrow(CMRE5_r8100_roots_Bx_br_pairlist)) #make new column Sampletype and fill all Root
CMRE5_r8100_roots_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE5_r8100_roots_Bx_br_pairlist)) #make new column Sampletype and fill all Bx

##

CMRE5_r8100_t_roots_Tb_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_roots_Tb_br))) ###converting distance matrix to list 
CMRE5_r8100_roots_Tb_br_pairlist <- CMRE5_r8100_t_roots_Tb_br_pairlist[as.numeric(CMRE5_r8100_t_roots_Tb_br_pairlist$col) > as.numeric(CMRE5_r8100_t_roots_Tb_br_pairlist$row),]
CMRE5_r8100_roots_Tb_br_pairlist$Sampletype <- rep("Root",nrow(CMRE5_r8100_roots_Tb_br_pairlist)) #make new column Sampletype and fill all Root
CMRE5_r8100_roots_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE5_r8100_roots_Tb_br_pairlist)) #make new column Sampletype and fill all Tb

##

CMRE5_r8100_t_roots_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_roots_Ctr_br))) ###converting distance matrix to list 
CMRE5_r8100_roots_Ctr_br_pairlist <- CMRE5_r8100_t_roots_Ctr_br_pairlist[as.numeric(CMRE5_r8100_t_roots_Ctr_br_pairlist$col) > as.numeric(CMRE5_r8100_t_roots_Ctr_br_pairlist$row),]
CMRE5_r8100_roots_Ctr_br_pairlist$Sampletype <- rep("Root",nrow(CMRE5_r8100_roots_Ctr_br_pairlist)) #make new column Sampletype and fill all Root
CMRE5_r8100_roots_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE5_r8100_roots_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE5_r8100_roots_pairlist <- rbind(CMRE5_r8100_roots_Ctr_br_pairlist, CMRE5_r8100_roots_As_br_pairlist, CMRE5_r8100_roots_Bx_br_pairlist, CMRE5_r8100_roots_Tb_br_pairlist)

######################################

#######--For soil -----#######

CMRE5_r8100_t_soil_As <- subset_samples(CMRE5_r8100_t_soil, Treatment%in%c("As"))
CMRE5_r8100_t_soil_Bx <- subset_samples(CMRE5_r8100_t_soil, Treatment%in%c("Bx"))
CMRE5_r8100_t_soil_Tb <- subset_samples(CMRE5_r8100_t_soil, Treatment%in%c("Tb"))
CMRE5_r8100_t_soil_Ctr <- subset_samples(CMRE5_r8100_t_soil, Treatment%in%c("Ctr"))

CMRE5_r8100_t_soil_As_br <- distance(CMRE5_r8100_t_soil_As, method = "bray") #bray-curtis 
CMRE5_r8100_t_soil_Bx_br <- distance(CMRE5_r8100_t_soil_Bx, method = "bray")
CMRE5_r8100_t_soil_Tb_br <- distance(CMRE5_r8100_t_soil_Tb, method = "bray")
CMRE5_r8100_t_soil_Ctr_br <- distance(CMRE5_r8100_t_soil_Ctr, method = "bray")

CMRE5_r8100_t_soil_As_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_soil_As_br))) ###converting distance matrix to list 
CMRE5_r8100_soil_As_br_pairlist <- CMRE5_r8100_t_soil_As_br_pairlist[as.numeric(CMRE5_r8100_t_soil_As_br_pairlist$col) > as.numeric(CMRE5_r8100_t_soil_As_br_pairlist$row),]
CMRE5_r8100_soil_As_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE5_r8100_soil_As_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE5_r8100_soil_As_br_pairlist$Treatment <- rep("As",nrow(CMRE5_r8100_soil_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE5_r8100_t_soil_Bx_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_soil_Bx_br))) ###converting distance matrix to list 
CMRE5_r8100_soil_Bx_br_pairlist <- CMRE5_r8100_t_soil_Bx_br_pairlist[as.numeric(CMRE5_r8100_t_soil_Bx_br_pairlist$col) > as.numeric(CMRE5_r8100_t_soil_Bx_br_pairlist$row),]
CMRE5_r8100_soil_Bx_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE5_r8100_soil_Bx_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE5_r8100_soil_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE5_r8100_soil_Bx_br_pairlist)) #make new column Sampletype and fill all Bx
write.csv(CMRE5_r8100_soil_Bx_br_pairlist, file = "CMRE5_r8100_soil_Bx_br_pairlist.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save as CSV file
CMRE5_r8100_soil_Bx_br_pairlist = read.csv(file = "CMRE5_r8100_soil_Bx_br_pairlist.CSV", header = TRUE)

##

CMRE5_r8100_t_soil_Tb_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_soil_Tb_br))) ###converting distance matrix to list 
CMRE5_r8100_soil_Tb_br_pairlist <- CMRE5_r8100_t_soil_Tb_br_pairlist[as.numeric(CMRE5_r8100_t_soil_Tb_br_pairlist$col) > as.numeric(CMRE5_r8100_t_soil_Tb_br_pairlist$row),]
CMRE5_r8100_soil_Tb_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE5_r8100_soil_Tb_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE5_r8100_soil_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE5_r8100_soil_Tb_br_pairlist)) #make new column Sampletype and fill all Tb
write.csv(CMRE5_r8100_soil_Tb_br_pairlist, file = "CMRE5_r8100_soil_Tb_br_pairlist.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save as CSV file
CMRE5_r8100_soil_Tb_br_pairlist = read.csv(file = "CMRE5_r8100_soil_Tb_br_pairlist.CSV", header = TRUE)

##

CMRE5_r8100_t_soil_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_soil_Ctr_br))) ###converting distance matrix to list 
CMRE5_r8100_soil_Ctr_br_pairlist <- CMRE5_r8100_t_soil_Ctr_br_pairlist[as.numeric(CMRE5_r8100_t_soil_Ctr_br_pairlist$col) > as.numeric(CMRE5_r8100_t_soil_Ctr_br_pairlist$row),]
CMRE5_r8100_soil_Ctr_br_pairlist$Sampletype <- rep("Soil",nrow(CMRE5_r8100_soil_Ctr_br_pairlist)) #make new column Sampletype and fill all Soil
CMRE5_r8100_soil_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE5_r8100_soil_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE5_r8100_soil_pairlist <- rbind(CMRE5_r8100_soil_Ctr_br_pairlist, CMRE5_r8100_soil_As_br_pairlist, CMRE5_r8100_soil_Bx_br_pairlist, CMRE5_r8100_soil_Tb_br_pairlist)

######################################

#######--For water -----#######

CMRE5_r8100_t_water_As <- subset_samples(CMRE5_r8100_t_water, Treatment%in%c("As"))
CMRE5_r8100_t_water_Bx <- subset_samples(CMRE5_r8100_t_water, Treatment%in%c("Bx"))
CMRE5_r8100_t_water_Tb <- subset_samples(CMRE5_r8100_t_water, Treatment%in%c("Tb"))
CMRE5_r8100_t_water_Ctr <- subset_samples(CMRE5_r8100_t_water, Treatment%in%c("Ctr"))

CMRE5_r8100_t_water_As_br <- distance(CMRE5_r8100_t_water_As, method = "bray") #bray-curtis 
CMRE5_r8100_t_water_Bx_br <- distance(CMRE5_r8100_t_water_Bx, method = "bray")
CMRE5_r8100_t_water_Tb_br <- distance(CMRE5_r8100_t_water_Tb, method = "bray")
CMRE5_r8100_t_water_Ctr_br <- distance(CMRE5_r8100_t_water_Ctr, method = "bray")

CMRE5_r8100_t_water_As_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_water_As_br))) ###converting distance matrix to list 
CMRE5_r8100_water_As_br_pairlist <- CMRE5_r8100_t_water_As_br_pairlist[as.numeric(CMRE5_r8100_t_water_As_br_pairlist$col) > as.numeric(CMRE5_r8100_t_water_As_br_pairlist$row),]
CMRE5_r8100_water_As_br_pairlist$Sampletype <- rep("Water",nrow(CMRE5_r8100_water_As_br_pairlist)) #make new column Sampletype and fill all Water
CMRE5_r8100_water_As_br_pairlist$Treatment <- rep("As",nrow(CMRE5_r8100_water_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE5_r8100_t_water_Bx_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_water_Bx_br))) ###converting distance matrix to list 
CMRE5_r8100_water_Bx_br_pairlist <- CMRE5_r8100_t_water_Bx_br_pairlist[as.numeric(CMRE5_r8100_t_water_Bx_br_pairlist$col) > as.numeric(CMRE5_r8100_t_water_Bx_br_pairlist$row),]
CMRE5_r8100_water_Bx_br_pairlist$Sampletype <- rep("Water",nrow(CMRE5_r8100_water_Bx_br_pairlist)) #make new column Sampletype and fill all Water
CMRE5_r8100_water_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE5_r8100_water_Bx_br_pairlist)) #make new column Sampletype and fill all Bx

##

CMRE5_r8100_t_water_Tb_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_water_Tb_br))) ###converting distance matrix to list 
CMRE5_r8100_water_Tb_br_pairlist <- CMRE5_r8100_t_water_Tb_br_pairlist[as.numeric(CMRE5_r8100_t_water_Tb_br_pairlist$col) > as.numeric(CMRE5_r8100_t_water_Tb_br_pairlist$row),]
CMRE5_r8100_water_Tb_br_pairlist$Sampletype <- rep("Water",nrow(CMRE5_r8100_water_Tb_br_pairlist)) #make new column Sampletype and fill all Water
CMRE5_r8100_water_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE5_r8100_water_Tb_br_pairlist)) #make new column Sampletype and fill all Tb

##

CMRE5_r8100_t_water_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_water_Ctr_br))) ###converting distance matrix to list 
CMRE5_r8100_water_Ctr_br_pairlist <- CMRE5_r8100_t_water_Ctr_br_pairlist[as.numeric(CMRE5_r8100_t_water_Ctr_br_pairlist$col) > as.numeric(CMRE5_r8100_t_water_Ctr_br_pairlist$row),]
CMRE5_r8100_water_Ctr_br_pairlist$Sampletype <- rep("Water",nrow(CMRE5_r8100_water_Ctr_br_pairlist)) #make new column Sampletype and fill all Water
CMRE5_r8100_water_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE5_r8100_water_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE5_r8100_water_pairlist <- rbind(CMRE5_r8100_water_Ctr_br_pairlist, CMRE5_r8100_water_As_br_pairlist, CMRE5_r8100_water_Bx_br_pairlist, CMRE5_r8100_water_Tb_br_pairlist)

######################################

#######--For sediment -----#######

CMRE5_r8100_t_sedi_As <- subset_samples(CMRE5_r8100_t_sedi, Treatment%in%c("As"))
CMRE5_r8100_t_sedi_Bx <- subset_samples(CMRE5_r8100_t_sedi, Treatment%in%c("Bx"))
CMRE5_r8100_t_sedi_Tb <- subset_samples(CMRE5_r8100_t_sedi, Treatment%in%c("Tb"))
CMRE5_r8100_t_sedi_Ctr <- subset_samples(CMRE5_r8100_t_sedi, Treatment%in%c("Ctr"))

CMRE5_r8100_t_sedi_As_br <- distance(CMRE5_r8100_t_sedi_As, method = "bray") #bray-curtis 
CMRE5_r8100_t_sedi_Bx_br <- distance(CMRE5_r8100_t_sedi_Bx, method = "bray")
CMRE5_r8100_t_sedi_Tb_br <- distance(CMRE5_r8100_t_sedi_Tb, method = "bray")
CMRE5_r8100_t_sedi_Ctr_br <- distance(CMRE5_r8100_t_sedi_Ctr, method = "bray")

CMRE5_r8100_t_sedi_As_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_sedi_As_br))) ###converting distance matrix to list 
CMRE5_r8100_sedi_As_br_pairlist <- CMRE5_r8100_t_sedi_As_br_pairlist[as.numeric(CMRE5_r8100_t_sedi_As_br_pairlist$col) > as.numeric(CMRE5_r8100_t_sedi_As_br_pairlist$row),]
CMRE5_r8100_sedi_As_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE5_r8100_sedi_As_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE5_r8100_sedi_As_br_pairlist$Treatment <- rep("As",nrow(CMRE5_r8100_sedi_As_br_pairlist)) #make new column Sampletype and fill all As

##

CMRE5_r8100_t_sedi_Bx_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_sedi_Bx_br))) ###converting distance matrix to list 
CMRE5_r8100_sedi_Bx_br_pairlist <- CMRE5_r8100_t_sedi_Bx_br_pairlist[as.numeric(CMRE5_r8100_t_sedi_Bx_br_pairlist$col) > as.numeric(CMRE5_r8100_t_sedi_Bx_br_pairlist$row),]
CMRE5_r8100_sedi_Bx_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE5_r8100_sedi_Bx_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE5_r8100_sedi_Bx_br_pairlist$Treatment <- rep("Bx",nrow(CMRE5_r8100_sedi_Bx_br_pairlist)) #make new column Sampletype and fill all Bx

##

CMRE5_r8100_t_sedi_Tb_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_sedi_Tb_br))) ###converting distance matrix to list 
CMRE5_r8100_sedi_Tb_br_pairlist <- CMRE5_r8100_t_sedi_Tb_br_pairlist[as.numeric(CMRE5_r8100_t_sedi_Tb_br_pairlist$col) > as.numeric(CMRE5_r8100_t_sedi_Tb_br_pairlist$row),]
CMRE5_r8100_sedi_Tb_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE5_r8100_sedi_Tb_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE5_r8100_sedi_Tb_br_pairlist$Treatment <- rep("Tb",nrow(CMRE5_r8100_sedi_Tb_br_pairlist)) #make new column Sampletype and fill all Tb

##

CMRE5_r8100_t_sedi_Ctr_br_pairlist <- dist2list(as.dist(t(CMRE5_r8100_t_sedi_Ctr_br))) ###converting distance matrix to list 
CMRE5_r8100_sedi_Ctr_br_pairlist <- CMRE5_r8100_t_sedi_Ctr_br_pairlist[as.numeric(CMRE5_r8100_t_sedi_Ctr_br_pairlist$col) > as.numeric(CMRE5_r8100_t_sedi_Ctr_br_pairlist$row),]
CMRE5_r8100_sedi_Ctr_br_pairlist$Sampletype <- rep("Sediment",nrow(CMRE5_r8100_sedi_Ctr_br_pairlist)) #make new column Sampletype and fill all Sediment
CMRE5_r8100_sedi_Ctr_br_pairlist$Treatment <- rep("Ctr",nrow(CMRE5_r8100_sedi_Ctr_br_pairlist)) #make new column Sampletype and fill all Ctr

##################now merge all treatment csv files

CMRE5_r8100_sedi_pairlist <- rbind(CMRE5_r8100_sedi_Ctr_br_pairlist, CMRE5_r8100_sedi_As_br_pairlist, CMRE5_r8100_sedi_Bx_br_pairlist, CMRE5_r8100_sedi_Tb_br_pairlist)

######################################

##################-------now merge all Sampletype csv files-------###############

CMRE5_r8100_pairlist <- rbind(CMRE5_r8100_feces_pairlist, CMRE5_r8100_roots_pairlist, CMRE5_r8100_soil_pairlist, CMRE5_r8100_water_pairlist, CMRE5_r8100_sedi_pairlist)

######################################---Plot Bray------######################################

CMRE5_r8100_pairlist$Sampletype<- factor(CMRE5_r8100_pairlist$Sampletype, levels = c("Water",  "Sediment", "Soil", "Root",  "Mouse"))
CMRE5_r8100_pairlist$Treatment<- factor(CMRE5_r8100_pairlist$Treatment, levels = c("Ctr", "As", "Bx", "Tb"))

CMRE5_r8100_pairlist_bray_boxplot <- ggplot(CMRE5_r8100_pairlist, aes(Sampletype, value, fill = CMRE5_r8100_pairlist$Treatment)) + geom_boxplot() + scale_fill_manual(values=c("gray50","#E69F00", "#009E73", "#D55E00")) + labs(x= " ", y = "Bray-Curtis") + theme_set(theme_bw()) + theme(legend.position="none", axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x  = element_text(colour="black", vjust=0.5, size=16), axis.text.y  = element_text(colour="black", vjust=0.5, size=16)) + scale_y_continuous(breaks = seq(0.1, 0.9, by=0.15), limits=c(0.1,0.9)) 
CMRE5_r8100_pairlist_bray_boxplot1 <- CMRE5_r8100_pairlist_bray_boxplot + scale_x_discrete(labels=c("Ctr As Bx Tb \nWater", "Ctr As Bx Tb \nSediment", "Ctr As Bx Tb \nSoil", "Ctr As Bx Tb  \nPlant", "Ctr As Bx Tb \nAnimal"))
CMRE5_r8100_pairlist_bray_boxplot2 <- CMRE5_r8100_pairlist_bray_boxplot1 + stat_summary(fun.y=mean, geom="point", colour="black", shape=18, size=2, position=position_dodge(width=0.75))  

##############################################################################################

##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure S7---------###############################################

#### water 

#water_CtrAs
 
Water_CtrAs = read.csv(file = "Water_CtrAs_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created. Made new column of Tax and ASV-ID
Water_CtrAs1 <- Water_CtrAs %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Water_CtrAs1$Treatment <- rep("As",nrow(Water_CtrAs1)) #make new column Treatment and fill all As
write.csv(Water_CtrAs1, file = "Water_CtrAs_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
#Water_CtrBx ---no taxa present
#Water_CtrTb
 
Water_CtrTb = read.csv(file = "Water_CtrTb_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created 
Water_CtrTb1 <- Water_CtrTb %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Water_CtrTb1$Treatment <- rep("Tb",nrow(Water_CtrTb1)) #make new column Treatment and fill all Tb
write.csv(Water_CtrTb1, file = "Water_CtrTb_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
 #merge all  
 Water_CtrAsBxTb_deseq2_sigOTU05_3 <- rbind(Water_CtrAs1, Water_CtrTb1)

#### prepare plot#####
  
 Water_CtrAsBxTb = read.csv(file = "Water_CtrAsBxTb_deseq2_sigOTU05_3.csv", header = TRUE, sep = ',')
 Water_CtrAsBxTb1 <- transform(Water_CtrAsBxTb, ASV_ID1=paste(Treatment, ASV_ID, sep="_"))
 Water_CtrAsBxTb1$Treatment <- factor(Water_CtrAsBxTb1$Treatment, levels = c("As", "Bx", "Tb"))
 x = tapply(Water_CtrAsBxTb1$log2FoldChange, Water_CtrAsBxTb1$ASV_ID1, function(x) max(x))
 x = sort(x, TRUE)
 Water_CtrAsBxTb1$ASV_ID1 = factor(as.character(Water_CtrAsBxTb1$ASV_ID1), levels = names(x))
 ##Plot 1
 p1 <-   ggplot(Water_CtrAsBxTb1, aes(x = log2FoldChange, y = ASV_ID1)) + 
   geom_vline(xintercept = 0.0, color = "gray", linewidth = 1.5) + geom_point(data=Water_CtrAsBxTb1, aes(colour = log2FoldChange >0), size = 3, alpha=0.7) + theme(axis.title.y=element_blank()) + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "log2FoldChange", y = "ASV-ID") + theme_bw() + theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, size=14), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=8)) + labs(color='log2FC') + scale_x_continuous(limits=c(-25, 5), breaks=c(-25, -20, -15, -10, -5, 0, 5)) + theme(legend.position="none")
 ##Plot 2
 p2 <-p1 + scale_colour_manual(name = 'log2FoldChange > 0', values = setNames(c("#1465AC",  "#B31B21"),c(T, F))) + theme(axis.text.x = element_text( hjust = 0, vjust=0.5)) 
 ##Plot 3
 p3 <- p2 + scale_y_discrete(labels = function(labels) sapply(labels, function(label) {
   parts <- strsplit(as.character(label), "_")[[1]]
   if (length(parts) > 1) {
     paste(parts[-1], collapse = "_")  # Re-join parts if there's more than one underscore
   } else {
     label  # Use the original label if there's no underscore
   }
 }))
 ##Final plot
 p4 <-p3+ facet_grid(Treatment~.,scales = "free", space = "free", switch = "y")+ theme() +
   theme(strip.placement = "outside")  + theme(strip.text.y = element_text(angle = 0)) + theme(strip.text.y = element_text(angle = 180, size=11, colour="black", face = "bold")) + theme(strip.background = element_rect(fill="lightblue", colour="black"), plot.title = element_text(size=16, face="bold",hjust=0.25, vjust = 0.10)) + ggtitle("Water")

#### Sediment 

#Sediment_CtrAs
 
Sediment_CtrAs = read.csv(file = "Sediment_CtrAs_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created. Made new column of Tax and ASV-ID
Sediment_CtrAs1 <- Sediment_CtrAs %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Sediment_CtrAs1$Treatment <- rep("As",nrow(Sediment_CtrAs1)) #make new column Treatment and fill all As
write.csv(Sediment_CtrAs1, file = "Sediment_CtrAs_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
#Sediment_CtrBx 

Sediment_CtrBx = read.csv(file = "Sediment_CtrBx_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created. Made new column of Tax and ASV-ID
Sediment_CtrBx1 <- Sediment_CtrBx %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Sediment_CtrBx1$Treatment <- rep("Bx",nrow(Sediment_CtrBx1)) #make new column Treatment and fill all As
write.csv(Sediment_CtrBx1, file = "Sediment_CtrBx_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
                             
#Sediment_CtrTb
 
Sediment_CtrTb = read.csv(file = "Sediment_CtrTb_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created 
Sediment_CtrTb1 <- Sediment_CtrTb %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Sediment_CtrTb1$Treatment <- rep("Tb",nrow(Sediment_CtrTb1)) #make new column Treatment and fill all Tb
write.csv(Sediment_CtrTb1, file = "Sediment_CtrTb_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
 #merge all  
 Sediment_CtrAsBxTb_deseq2_sigOTU05_3 <- rbind(Sediment_CtrAs1, Sediment_CtrBx1, Sediment_CtrTb1)

#### prepare plot#####
  
 Sediment_CtrAsBxTb = read.csv(file = "Sediment_CtrAsBxTb_deseq2_sigOTU05_3.csv", header = TRUE, sep = ',')
 Sediment_CtrAsBxTb1 <- transform(Sediment_CtrAsBxTb, ASV_ID1=paste(Treatment, ASV_ID, sep="_"))
 Sediment_CtrAsBxTb1$Treatment <- factor(Sediment_CtrAsBxTb1$Treatment, levels = c("As", "Bx", "Tb"))
 x = tapply(Sediment_CtrAsBxTb1$log2FoldChange, Sediment_CtrAsBxTb1$ASV_ID1, function(x) max(x))
 x = sort(x, TRUE)
 Sediment_CtrAsBxTb1$ASV_ID1 = factor(as.character(Sediment_CtrAsBxTb1$ASV_ID1), levels = names(x))
 ##Plot 1
 p1 <-   ggplot(Sediment_CtrAsBxTb1, aes(x = log2FoldChange, y = ASV_ID1)) + 
   geom_vline(xintercept = 0.0, color = "gray", linewidth = 1.5) + geom_point(data=Sediment_CtrAsBxTb1, aes(colour = log2FoldChange >0), size = 3, alpha=0.7) + theme(axis.title.y=element_blank()) + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "log2FoldChange", y = "ASV-ID") + theme_bw() + theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, size=14), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=8)) + labs(color='log2FC') + scale_x_continuous(limits=c(-20, 6), breaks=c(-20, -15, -10, -5,0, 5)) + theme(legend.position="none")
 ##Plot 2
 p2 <-p1 + scale_colour_manual(name = 'log2FoldChange > 0', values = setNames(c("#1465AC",  "#B31B21"),c(T, F))) + theme(axis.text.x = element_text( hjust = 0, vjust=0.5)) 
 ##Plot 3
 p3 <- p2 + scale_y_discrete(labels = function(labels) sapply(labels, function(label) {
   parts <- strsplit(as.character(label), "_")[[1]]
   if (length(parts) > 1) {
     paste(parts[-1], collapse = "_")  # Re-join parts if there's more than one underscore
   } else {
     label  # Use the original label if there's no underscore
   }
 }))
 ##Final plot
 p4 <-p3+ facet_grid(Treatment~.,scales = "free", space = "free", switch = "y")+ theme() +
   theme(strip.placement = "outside")  + theme(strip.text.y = element_text(angle = 0)) + theme(strip.text.y = element_text(angle = 180, size=11, colour="black", face = "bold")) + theme(strip.background = element_rect(fill="lightblue", colour="black"), plot.title = element_text(size=16, face="bold",hjust=0.25, vjust = 0.10)) + ggtitle("Sediment")


#### Soil 

#Soil_CtrAs
 
Soil_CtrAs = read.csv(file = "Soil_CtrAs_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created. Made new column of Tax and ASV-ID
Soil_CtrAs1 <- Soil_CtrAs %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Soil_CtrAs1$Treatment <- rep("As",nrow(Soil_CtrAs1)) #make new column Treatment and fill all As
write.csv(Soil_CtrAs1, file = "Soil_CtrAs_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
#Soil_CtrBx 

Soil_CtrBx = read.csv(file = "Soil_CtrBx_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created. Made new column of Tax and ASV-ID
Soil_CtrBx1 <- Soil_CtrBx %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Soil_CtrBx1$Treatment <- rep("Bx",nrow(Soil_CtrBx1)) #make new column Treatment and fill all As
write.csv(Soil_CtrBx1, file = "Soil_CtrBx_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
                             
#Soil_CtrTb
 
Soil_CtrTb = read.csv(file = "Soil_CtrTb_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created 
Soil_CtrTb1 <- Soil_CtrTb %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Soil_CtrTb1$Treatment <- rep("Tb",nrow(Soil_CtrTb1)) #make new column Treatment and fill all Tb
write.csv(Soil_CtrTb1, file = "Soil_CtrTb_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
 #merge all  
 Soil_CtrAsBxTb_deseq2_sigOTU05_3 <- rbind(Soil_CtrAs1, Soil_CtrBx1, Soil_CtrTb1)

#### prepare plot#####
  
 Soil_CtrAsBxTb = read.csv(file = "Soil_CtrAsBxTb_deseq2_sigOTU05_3.csv", header = TRUE, sep = ',')
 Soil_CtrAsBxTb1 <- transform(Soil_CtrAsBxTb, ASV_ID1=paste(Treatment, ASV_ID, sep="_"))
 Soil_CtrAsBxTb1$Treatment <- factor(Soil_CtrAsBxTb1$Treatment, levels = c("As", "Bx", "Tb"))
 x = tapply(Soil_CtrAsBxTb1$log2FoldChange, Soil_CtrAsBxTb1$ASV_ID1, function(x) max(x))
 x = sort(x, TRUE)
 Soil_CtrAsBxTb1$ASV_ID1 = factor(as.character(Soil_CtrAsBxTb1$ASV_ID1), levels = names(x))
 ##Plot 1
 p1 <-   ggplot(Soil_CtrAsBxTb1, aes(x = log2FoldChange, y = ASV_ID1)) + 
   geom_vline(xintercept = 0.0, color = "gray", linewidth = 1.5) + geom_point(data=Soil_CtrAsBxTb1, aes(colour = log2FoldChange >0), size = 3, alpha=0.7) + theme(axis.title.y=element_blank()) + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "log2FoldChange", y = "ASV-ID") + theme_bw() + theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, size=14), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=8)) + labs(color='log2FC') + scale_x_continuous(limits=c(-2, 4), breaks=c(-2, 1, 0, 1, 2, 3, 4)) + theme(legend.position="none")
 ##Plot 2
 p2 <-p1 + scale_colour_manual(name = 'log2FoldChange > 0', values = setNames(c("#1465AC",  "#B31B21"),c(T, F))) + theme(axis.text.x = element_text( hjust = 0, vjust=0.5)) 
 ##Plot 3
 p3 <- p2 + scale_y_discrete(labels = function(labels) sapply(labels, function(label) {
   parts <- strsplit(as.character(label), "_")[[1]]
   if (length(parts) > 1) {
     paste(parts[-1], collapse = "_")  # Re-join parts if there's more than one underscore
   } else {
     label  # Use the original label if there's no underscore
   }
 }))
 ##Final plot
 p4 <-p3+ facet_grid(Treatment~.,scales = "free", space = "free", switch = "y")+ theme() +
   theme(strip.placement = "outside")  + theme(strip.text.y = element_text(angle = 0)) + theme(strip.text.y = element_text(angle = 180, size=11, colour="black", face = "bold")) + theme(strip.background = element_rect(fill="lightblue", colour="black"), plot.title = element_text(size=16, face="bold",hjust=0.25, vjust = 0.10)) + ggtitle("Soil")

#### Animal 

#Mouse_CtrAs
 
Mouse_CtrAs = read.csv(file = "Mouse_CtrAs_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created. Made new column of Tax and ASV-ID
Mouse_CtrAs1 <- Mouse_CtrAs %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Mouse_CtrAs1$Treatment <- rep("As",nrow(Mouse_CtrAs1)) #make new column Treatment and fill all As
write.csv(Mouse_CtrAs1, file = "Mouse_CtrAs_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
#Mouse_CtrBx 

Mouse_CtrBx = read.csv(file = "Mouse_CtrBx_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created. Made new column of Tax and ASV-ID
Mouse_CtrBx1 <- Mouse_CtrBx %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Mouse_CtrBx1$Treatment <- rep("Bx",nrow(Mouse_CtrBx1)) #make new column Treatment and fill all As
write.csv(Mouse_CtrBx1, file = "Mouse_CtrBx_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
                             
#Mouse_CtrTb
 
Mouse_CtrTb = read.csv(file = "Mouse_CtrTb_deseq2_sigOTU05_2.csv", header = TRUE, sep = ',') # these files are manually created 
Mouse_CtrTb1 <- Mouse_CtrTb %>% mutate(ASV_ID = paste(Tax, ASV.ID, sep = "_"))
Mouse_CtrTb1$Treatment <- rep("Tb",nrow(Mouse_CtrTb1)) #make new column Treatment and fill all Tb
write.csv(Mouse_CtrTb1, file = "Mouse_CtrTb_deseq2_sigOTU05_3.CSV", row.names = TRUE, sep = ',', col.names = TRUE) #save As CSV file
 
 #merge all  
 Mouse_CtrAsBxTb_deseq2_sigOTU05_3 <- rbind(Mouse_CtrAs1, Mouse_CtrBx1, Mouse_CtrTb1)

#### prepare plot#####
  
 Mouse_CtrAsBxTb = read.csv(file = "Mouse_CtrAsBxTb_deseq2_sigOTU05_3.csv", header = TRUE, sep = ',')
 Mouse_CtrAsBxTb1 <- transform(Mouse_CtrAsBxTb, ASV_ID1=paste(Treatment, ASV_ID, sep="_"))
 Mouse_CtrAsBxTb1$Treatment <- factor(Mouse_CtrAsBxTb1$Treatment, levels = c("As", "Bx", "Tb"))
 x = tapply(Mouse_CtrAsBxTb1$log2FoldChange, Mouse_CtrAsBxTb1$ASV_ID1, function(x) max(x))
 x = sort(x, TRUE)
 Mouse_CtrAsBxTb1$ASV_ID1 = factor(as.character(Mouse_CtrAsBxTb1$ASV_ID1), levels = names(x))
 ##Plot 1
 p1 <-   ggplot(Mouse_CtrAsBxTb1, aes(x = log2FoldChange, y = ASV_ID1)) + 
   geom_vline(xintercept = 0.0, color = "gray", linewidth = 1.5) + geom_point(data=Mouse_CtrAsBxTb1, aes(colour = log2FoldChange >0), size = 3, alpha=0.7) + theme(axis.title.y=element_blank()) + theme(axis.text.x = element_text(hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + labs(x= "log2FoldChange", y = "ASV-ID") + theme_bw() + theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x  = element_text(colour="black", face = "bold", vjust=0.5, size=14), axis.text.y  = element_text(colour="black", face = "bold", vjust=0.5, size=8)) + labs(color='log2FC') + scale_x_continuous(limits=c(-2, 4), breaks=c(-2, 1, 0, 1, 2, 3, 4)) + theme(legend.position="none")
 ##Plot 2
 p2 <-p1 + scale_colour_manual(name = 'log2FoldChange > 0', values = setNames(c("#1465AC",  "#B31B21"),c(T, F))) + theme(axis.text.x = element_text( hjust = 0, vjust=0.5)) 
 ##Plot 3
 p3 <- p2 + scale_y_discrete(labels = function(labels) sapply(labels, function(label) {
   parts <- strsplit(as.character(label), "_")[[1]]
   if (length(parts) > 1) {
     paste(parts[-1], collapse = "_")  # Re-join parts if there's more than one underscore
   } else {
     label  # Use the original label if there's no underscore
   }
 }))
 ##Final plot
 p4 <-p3+ facet_grid(Treatment~.,scales = "free", space = "free", switch = "y")+ theme() +
   theme(strip.placement = "outside")  + theme(strip.text.y = element_text(angle = 0)) + theme(strip.text.y = element_text(angle = 180, size=11, colour="black", face = "bold")) + theme(strip.background = element_rect(fill="lightblue", colour="black"), plot.title = element_text(size=16, face="bold",hjust=0.25, vjust = 0.10)) + ggtitle("Animal")
                              
##############################################################################################

##############################################################################################################################

##############################################################################################################################

##########################################------Supplementary Figure S9---------###############################################

###################--Network degree distribution

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
CMRE5_Sedimentb_Ctrc_Sediment_ig= readRDS("CMRE5_Sedimentb_Ctrc_Sediment_ig.RDS")
CMRE5_Sedimentb_Asc_Sediment_ig= readRDS("CMRE5_Sedimentb_Asc_Sediment_ig.RDS")
CMRE5_Sedimentb_Bxc_Sediment_ig= readRDS("CMRE5_Sedimentb_Bxc_Sediment_ig.RDS")
CMRE5_Sedimentb_Tbc_Sediment_ig= readRDS("CMRE5_Sedimentb_Tbc_Sediment_ig.RDS")

dd.Sediment.Ctr.Sediment.ig <- degree.distribution(CMRE5_Sedimentb_Ctrc_Sediment_ig)
dd.Sediment.As.Sediment.ig <- degree.distribution(CMRE5_Sedimentb_Asc_Sediment_ig)
dd.Sediment.Bx.Sediment.ig <- degree.distribution(CMRE5_Sedimentb_Bxc_Sediment_ig)
dd.Sediment.Tb.Sediment.ig <- degree.distribution(CMRE5_Sedimentb_Tbc_Sediment_ig)

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
CMRE5_Waterb_Ctrc_Water_ig= readRDS("CMRE5_Waterb_Ctrc_Water_ig.RDS")
CMRE5_Waterb_Asc_Water_ig= readRDS("CMRE5_Waterb_Asc_Water_ig.RDS")
CMRE5_Waterb_Bxc_Water_ig= readRDS("CMRE5_Waterb_Bxc_Water_ig.RDS")
CMRE5_Waterb_Tbc_Water_ig= readRDS("CMRE5_Waterb_Tbc_Water_ig.RDS")

dd.Water.Ctr.Water.ig <- degree.distribution(CMRE5_Waterb_Ctrc_Water_ig)
dd.Water.As.Water.ig <- degree.distribution(CMRE5_Waterb_Asc_Water_ig)
dd.Water.Bx.Water.ig <- degree.distribution(CMRE5_Waterb_Bxc_Water_ig)
dd.Water.Tb.Water.ig <- degree.distribution(CMRE5_Waterb_Tbc_Water_ig)

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


