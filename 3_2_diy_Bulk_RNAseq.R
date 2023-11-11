# (reference) https://diytranscriptomics.com/

# Introduction to the Step 1 script----
# Step 1 Learning Objectives:
# 1 - Step 1 serves as your gateway to R scripts and, as such, 
#     you will learn the proper 'anatomy' for any R script.
# 2 - Learn how to install packages and load libraries into your R environment
# 3 - Understand the various file types that describe RNAseq data and 
#     how to import these files (e.g. kallisto read mapping data) into R
# 4 - Learn basic tools for annotation


# set working directory (always check~!)
setwd("/Users/kwangwooko/Desktop/Coding_practice/diy/kallisto")

# load packages----
setRepositories()

#provides functions for handling hdf5 file formats 
# (kallisto outputs bootstraps in this format)
library(rhdf5)
# provides access to Hadley Wickham's collection of R packages for data science, 
# which we will use throughout the course
library(tidyverse) 
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl

#replace with your organism-specific database package
library(EnsDb.Hsapiens.v86) 

library(biomaRt) # an alternative for annotation
library(beepr) #just for fun
library(datapasta) # great for copy / paste data into the R environment


# read in your study design ----
#there are LOTS of ways to read data into R, but the readr package 
# (from tidyverse) is one of the simplest
# I noticed that read_tsv function is colliding with EnsDb.Hsapiens.v86 package
targets <- read_tsv("studydesign.txt")

# you can easily create file paths to the abundance files generated 
# by Kallisto using the 'file.path' function
# set file paths to your mapped data
path <- file.path(targets$sample, "abundance.tsv") 
# now check to make sure this path is correct by seeing if the files exist
all(file.exists(path)) 
# this should look familiar from your homework
# the functions which(), any() and all() take logical vectors as their argument. 
# The any() function will return TRUE if one or more of the elements 
# in the logical vector is TRUE. The all() function
# will return TRUE if every element in the logical vector is TRUE.

# get annotations using organism-specific package ----
# go to the bioconductor website and then look at many annotation packages

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx) # to make a data frame

# need to change first column name of Tx to 'target_id'
# "abundance.tsv" in kallisto has target_id
Tx <- dplyr::rename(Tx, target_id = tx_id) 
#transcrip ID needs to be the first column in the dataframe
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# OPTIONAL: get annotations using BiomaRt----
# The annotation method described in the code chunk above works great 
# if an organism-specific data base package exists for the organisms of interest
# however, this is only the case for human, mouse and rat....
# this optional shows one way to get annotation data for other target organisms
# in this example, we're retrieving 1:1 mappings between transcript identifiers 
# and gene symbols for the domesticated dog (Canis familiaris)

library(biomaRt) # an alternative for annotation

#default host is ensembl.org, and most current release of mammalian genomes
listMarts() 
#listMarts(host="parasite.wormbase.org") #access to parasite worm genomes
#listMarts(host="protists.ensembl.org") #access to protozoan genomes

#choose the 'mart' you want to work with
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
#take a look at all available datasets within the selected mart
available.datasets <- listDatasets(myMart)
#now grab the ensembl annotations for dog
dog.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                    dataset = "clfamiliaris_gene_ensembl")
dog.attributes <- listAttributes(dog.anno)
dog.filters <- listFilters(dog.anno)


Tx.dog <- getBM(attributes=c('ensembl_transcript_id_version',
                         'external_gene_name'),
            mart = dog.anno)

Tx.dog <- as_tibble(Tx.dog)
#we need to rename the two columns we just retreived from biomart
Tx.dog <- dplyr::rename(Tx.dog, target_id = ensembl_transcript_id_version, 
                    gene_name = external_gene_name)


# import Kallisto transcript counts into R using Tximport ----
# txOut (transcript is the result of Output?)
# ignoreTxVersion = TRUE :whether transcripts versions (isoforms) are ignored

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #FALSE: genes vs TRUE: transcripts
                     countsFromAbundance = "lengthScaledTPM", 
                     ignoreTxVersion = TRUE)  # transcripts: different versions.
beep(sound = 6)

#take a look at the type of object you just created
class(Txi_gene)
names(Txi_gene)

# the essentials ----
# this chunk contains the minimal essential code from this script. 
# provides access to Hadley Wickham's collection of R packages for data science
library(tidyverse) 
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
#replace with your organism-specific database package
library(EnsDb.Hsapiens.v86) 

targets <- read_tsv("studydesign.txt", "\t") # read in your study design
path <- file.path(targets$sample, "abundance.tsv") 
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE) 

Txi_transcript <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = TRUE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE) 



###############################################################################################################################################


# Introduction to the Step 2 script ----
# Now that you've read your transcript-level or gene-level data into R, you're ready to begin working with your data.

# Step 2 Learning Objectives:
# 1 - Filter and normalize your data
# 2 - use ggplot2 to visualize the impact of filtering and normalization on our data.
# 3 - understand why gene expression data is messy, and how to make it 'tidy'

# Notes:
# recall that your abundance data are TPM, while the counts are read counts mapping to each gene or transcript

# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure

# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myTPM_transcript <- Txi_transcript$abundance
myCounts <- Txi_gene$counts
myCounts_transcript <- Txi_transcript$counts

colSums(myTPM)
colSums(myTPM_transcript)
colSums(myCounts)
colSums(myCounts_transcript)

# capture sample labels from the study design file that you worked with and saved as 'targets' in step 1
targets
sampleLabels <- targets$sample

# Generate summary stats for your data ----
# 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
# then we use the 'rowSds', 'rowMeans' and 'rowMedians' functions from the matrixStats package
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# look at what you created
head(myTPM.stats)

# Create your first plot using ggplot2 ----
# produce a scatter plot of the transformed data
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=25, size=3)
# Experiment with point shape and size in the plot above
# Experiment with other plot types (e.g. 'geom_hex' instead of 'geom_point')
# Add a theme to your ggplot code above.  Try 'theme_bw()'
# How would these graphs change if you log2 converted the data?

# Let's expand on the plot above a bit more and take a look at each 'layer' of the ggplot code
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=16, size=2) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = FALSE) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  theme_classic() +
  theme_dark() + 
  theme_bw()

# Make a DGElist from your counts, and plot ----
myDGEList <- DGEList(myCounts)
# take a look at the DGEList object 
myDGEList
#DEGList objects are a good R data file to consider saving to you working directory
save(myDGEList, file = "myDGEList")
#Saved DGEList objects can be easily shared and loaded into an R environment
load(file = "myDGEList")

# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(myDGEList) 
colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
log2.cpm.df
# add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
# use the tidy package to 'pivot' your dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

# let's look at the impact of pivoting the data
log2.cpm.df.pivot

# not it is easy to plot this pivoted data
ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
  
# what do you think of the distribution of this data?
# Try using coord_flip() at the end of the ggplot code

# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(myDGEList$counts==0)==10)
# breaking down the line above is a little tricky.  Let's try:
# 1st - 'myDGEList$counts==0' returns a new 'logical matrix' where each observation (gene) is evaluated (TRUE/FALSE) for each variable (sample) as to whether it has zero counts
# 2nd - passing this logical matrix to 'rowsums' allows you to sum the total number of times an observation was 'TRUE' across all samples
# 3rd - adding the '==10' is a simple way of asking how many of the rowsums equaled 10. In other words, how many genes had 0 counts (TRUE) for all samples in our dataset
# 4th - passing all this to the 'table' function just provides a handy way to summarize the large logical produced in the previous step

# now set some cut-off to get rid of genes/transcripts with low counts
# again using rowSums to tally up the 'TRUE' results of a simple evaluation
# how many genes had more than 1 CPM (TRUE) in at least 3 samples

# The line below is important! This is where the filtering starts
# Be sure to adjust this cutoff for the number of samples in the smallest group of comparison.
keepers <- rowSums(cpm>1)>=5
# now use base R's simple subsetting method to filter your DGEList based on the logical produced above
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)


ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Normalize your data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
# take a look at this new DGEList object...how has it changed?

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
# pivot this NORMALIZED data, just as you did earlier
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)


ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# what if we wanted to put all three violin plots together?
# go back and assign each plot to a variable (rather than printing to the plots viewer)
# here we assigned the last 3 plots to p1, p2 and p3
# we'll use the 'plot_grid' function from the cowplot package to put these together in a figure
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
print("Step 2 complete!")

# the essentials ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

#p1: unfiltered, non-normalized

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

#p2: filtered, non-normalized

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


#p3: filtered, normalized

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)


p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

====================================================================================================


# Introduction to this script -----------
# the goal of this script is to identify differentially expressed genes (DEGs) 
# and differential transcript usage (DTU)
# you should already know which pairwise comparisons are most important to you
# whether you look for differential expression at the gene or transcript level 
# depends on how you read the Kallisto output into R using TxImport (Step 1)
# if you have no biological replicates, you will NOT be able to leverage 
# statistical tools for differential expression analysis
# instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' 
# we discussed in Step 3 and 4 to identify genes based on log fold-change

# Load packages -----
# limma: venerable package for differential gene expression using linear model

library(tidyverse)
library(limma) 
library(edgeR)
library(gt)
library(DT)
library(plotly)

# Set up your design matrix ----
group <- factor(targets$group)
design <- model.matrix(~0 + group) # ~0: no intercept
colnames(design) <- levels(group)

# NOTE: if you need a paired analysis (a.k.a.'blocking' design) or 
# have a batch effect, the following design is useful
# design <- model.matrix(~block + treatment)
# this is just an example. 'block' and 'treatment' would need to be objects 
# in your environment

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix ----
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", 
                      coef=1, number=40000, sort.by="logFC")
myTop_10Hits <- topTable(ebFit, adjust ="BH", 
                         coef=1, number=10, sort.by="logFC")
# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

myTop_10Hits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTop_10Hits.df)
# TopTable (from Limma) outputs a few different stats:
# logFC, AveExpr, and P.Value should be self-explanatory
# adj.P.Val is your adjusted P value, also known as an FDR 
# (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. 
# If B = 1.5, then log odds is e^1.5, where e is euler's constant (about 2.718).
# So, the odds of differential expression os about 4.8 to 1
# t statistic is ratio of the logFC to the standard error (where the error 
# has been moderated across all genes...because of Bayesian approach)

# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

# now plot
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", 
             size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, 
           alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, 
           alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Now make the volcano plot above interactive with plotly
ggplotly(vplot)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", 
                       p.value=0.01, lfc=2)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="up")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E) # E: expression data
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,] # !=0: not equal to 0
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display your DEGs ----
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

#write your DEGs to a file
# NOTE: this .txt file can be directly used for input into other clustering or 
# network analysis tools 
# (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)
write_tsv(diffGenes.df,"DiffGenes.txt") 

# OPTIONAL: differential transcript usage (DTU) analysis ----
library(IsoformSwitchAnalyzeR)

# The IsoformSwitchAnalyzeR package looks for certain column headers 
# in our study design
# So, the first step is to make sure our study design contains the following:
# unique sample IDs must be contained in column called 'sampleID'
# covariate(s) of interest must be in column labeled 'condition'
# remove extraneous columns
targets.mod <- targets %>%
  dplyr::rename(sampleID = sample, condition = group) %>%
  dplyr::select(sampleID, condition)

# import transcript Kallisto quant data
# using the same path variable we set way back in the step 1 script
Txi_trans <- importIsoformExpression(sampleVector = path)

# fix column headers of abundance and counts data to match sampleID 
# in target.mod
colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels)
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels)

# import data
mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  ignoreAfterPeriod=TRUE,
  
  # the files below must be from the same ensembl release 
  # (in this case release 108), and must match the reference release version 
  # that we originally mapped our reads to at the beginning of the course
  # you can find version 108 of the gtf file below here: 
  # https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/
  isoformExonAnnoation = "Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz",
  isoformNtFasta       = "Homo_sapiens.GRCh38.cdna.all.fa",
  showProgress = TRUE)


# We'll do the isoform analysis in one step, but there's a lot to unpack here, 
# so you should really read the package documentation at:
# https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html
# Note that without additional manual work here (beyond the scope of this class)
# , we'll only capture isoform annotations for 1) intron retention; 
# 2) ORF sequence similarity; and 3) nonsense mediate decay (NMD)

#NOTE: THIS NEXT BIT COULD TAKE A WHILE!
mySwitchList <- isoformSwitchAnalysisCombined(
  switchAnalyzeRlist   = mySwitchList,
  pathToOutput = 'isoform_output') # directory must already exist

# now look at the directory that you just created above
# in case you missed the summary output from the function above
extractSwitchSummary(mySwitchList)

# extract the top n isoform switching events
# these 'consequences' related to the annotations I reference above.
# change to TRUE if you want this list sorted by FDR-adusted Pval (q value)
extractTopSwitches(
  mySwitchList,
  filterForConsequences = TRUE, 
  n = 50,
  sortByQvals = FALSE) 

# visualize by making a 'switch plot'
switchPlot(
  mySwitchList,
  gene='FCGR3B',
  condition1 = 'disease',
  condition2 = 'healthy',
  localTheme = theme_bw())

# the essentials ----
library(tidyverse)
library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

vplot <- ggplot(myTopHits) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", 
             size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, 
           alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, 
           alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", 
                       p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)





