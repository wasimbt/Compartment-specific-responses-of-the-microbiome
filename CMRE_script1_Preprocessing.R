
#Copyright (c) 2023 University of Bern, Bern, Switzerland. 
#Author: Wasimuddin (wasim.bt@gmail.com; https://www.researchgate.net/profile/Wasim-Uddin)



# library(dada2); packageVersion("dada2") # 
# ## [1] '1.10.1'
# library(ShortRead); packageVersion("ShortRead") # 
# ## [1] '1.40.0'
# library(Hmisc); packageVersion("Hmisc") # 
# ## [1] '4.1.1'
# library(ggplot2); packageVersion("ggplot2") # 
# ## [1] '3.1.1'
# library(phyloseq); packageVersion("phyloseq") # 
# ## [1] ‘1.26.1’
# library("data.table"); packageVersion("data.table")
# ## [1] ‘1.11.6’


#setwd("path/to/directory")### set the working directory

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # DADA2 pipeline

# Put filtered (after CUTADAPT) reads into separate sub-directories for analysis
 

datapath <- "path/to/data"
pathF <- "path/to/Fwd" # the directory containing demultiplexed, cutadapted forward fastq files
pathR <- "path/to/Rev" # the directory containing demultiplexed, cutadapted reverse fastq files

fastqFs <- list.files(pathF, pattern="fastq.gz") # CHANGE if different file extensions
fastqRs <- list.files(pathR, pattern="fastq.gz")

filtpathF <- file.path(datapath, "filteredF") # Filtered files go into the filtered/ subdirectory
filtpathR <- file.path(datapath, "filteredR") # Filtered files go into the filtered/ subdirectory


# Perform Quality Filtering, trimming etc

TLib_out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs), rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs), truncLen=c(200,105), maxEE=c(2, 2), truncQ=2, maxN=0, minLen=80, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)# 

####################### Remove samples with less than 1000 read also from sequence table 

TLib_out1 <- TLib_out[rowSums(TLib_out) >= 1000,]

################################################################################################################################

################################################################--------Learning error rates & correcting 


#######################Set up and verify the file names for the output:

# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

#learn forward error rate
set.seed(100)
errF <- learnErrors(filtFs, nbases=8e8, multithread=TRUE, verbose=TRUE, randomize = TRUE)

#learn reverse error rate
set.seed(100) 
errR <- learnErrors(filtRs, nbases=8e8, multithread=TRUE, verbose=TRUE, randomize = TRUE)


################################################################################################################################

################################################################--------Dereplication and merging
# perform the dereplication

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Do error correction on the dereplicated reads
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,  minOverlap=20, verbose=TRUE)

# Make sequence table
L_seqtab <- makeSequenceTable(mergers)


################################################################################################################################

############--Removing the samples with less than 1000 reads

L_seqtab1 <- L_seqtab[rowSums(L_seqtab) >= 1000,]


#######################--From this point tables from different runs should be merged and then downstream process should be continued


seqtab12345 <- mergeSequenceTables(L1_seqtab1, L2_seqtab1, L3_seqtab1, L4_seqtab1, L5_seqtab1)# merging sequence tables 


#######################################--------Chimera removing-------############

seqtab12345.nochim <- removeBimeraDenovo(seqtab12345, method="consensus", multithread=TRUE, verbose=TRUE)


###### ##### Remove too far length sequences

seqtab1234.nochim1 <- seqtab1234.nochim[,nchar(colnames(seqtab1234.nochim)) %in% seq(220,273)] 

######################################--------Assign taxonomy--------##########################

set.seed(100) # Initialize random number generator for reproducibility

CMRE1_taxa <- AssignTaxonomy(seqtab1234.nochim1, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)

CMRE2_taxa <- addSpecies(CMRE1_taxa, "silva_species_Assignment_v132.fa.gz", verbose=TRUE)

#################################################Merge the objects created with mapping file

CMRE1_object <- phyloseq(otu_table(seqtab1234.nochim1, taxa_are_rows=FALSE), sample_data(CMRE1234_metadata), tax_table(CMRE2_taxa))

######################--Change the sequence header to AsVs number for identity

n_seqs <- seq(ntaxa(CMRE1_object))

len_n_seqs <- nchar(max(n_seqs))

taxa_names(CMRE1_object) <- pAste("AsV", formatC(n_seqs, width = len_n_seqs, flag = "0"), sep = "_")

CMRE2 <- CMRE1_object

#######################################--Prepare the phylogenetic tree -#######################

CMRE2_random_tree = rtree(ntaxa(CMRE2), rooted=TRUE, tip.label=taxa_names(CMRE2)) ###First for unrooted tree

#########before runnging the below script, one hAs to run this function

pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch As outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}
#########then make the rooted tree from the unrooted tree

new.outgroup = pick_new_outgroup(CMRE2_random_tree)

CMRE2_rootedTree = ape::root(CMRE2_random_tree, outgroup=new.outgroup, resolve.root=TRUE)


#########Merge the created rooted tree

CMRE3 = merge_phyloseq(CMRE2, CMRE2_rootedTree)


#######################################------Removing Mt, chloroplast reads etc.---#################

CMRE3c <- subset_taxa(CMRE3, (Order!="Chloroplast") | is.na(Order)) ###### for chloroplast  

CMRE3cm <- subset_taxa(CMRE3c, (Family!="Mitochondria") | is.na(Family)) ###### for mitochondrial 

#############------Removing the taxa without Phylum Assignment (NA)

CMRE3cm_na <- subset_taxa(CMRE3cm, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#############Remove the taxa with less than 10 counts in overall dataset

CMRE4= prune_taxa(taxa_sums(CMRE3cm_na) > 10, CMRE3cm_na)

############################----------Also remove PCR Blanks and Extraction Blanks-----------######################################

CMRE4f # new object after removing blanks

#####---Here pre-processing is necessary as there many OTUs without a signle count

CMRE5 = prune_taxa(taxa_sums(CMRE4f) > 0, CMRE4f) #remove OTUs without a single count from class object

CMRE5


############################################---XXXXXXXXX------##############################################


############################################---XXXXXXXXX------##############################################





