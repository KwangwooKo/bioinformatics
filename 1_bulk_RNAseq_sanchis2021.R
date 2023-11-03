## (reference) https://pubmed.ncbi.nlm.nih.gov/33907739/


library(DESeq2)
library(GEOquery)
library(canvasXpress)
library(ggplot2)
library(clinfun)
library(GGally)
library(factoextra)


if(!dir.exists("Data")) {
  dir.create("Data")
}


data <- getGEO (GEO = "GSE152075", destdir="Data")  

# Extract the phenotypic/ clinical data matrix from the series matrix
# clindata is metadata
clindata <- data[["GSE152075_series_matrix.txt.gz"]]@phenoData@data

head(clindata[c(1,2,8,40,39,42)])

# stringsAsFactors = FALSE : character columns will remain as character vectors,
# preserving the original txt format.
raw_counts_pre <- read.delim("Data/GSE152075_raw_counts_GEO.txt.gz", 
                             stringsAsFactors=FALSE, sep = " ")
  
head(raw_counts[c(1:10)])
head(raw_counts[, c(1:10)])

## Gene expression normalization

# This code is for changing data frame (raw_counts+pre) to matrix (raw_counts)
raw_counts <- as.matrix(raw_counts_pre)

# replace the rownames of clindata (sample ID) 
# with the same sample name (title) of raw_counts.
# This will help to match sample names in both matrixes.
rownames(clindata) <- clindata$title

# the outcome should be TRUE
# %in% : is in operator, check if the elements on the left side are present 
# in the elements on the right side
# all() function is used to check if a condition is TRUE for all elements
all(rownames(clindata) %in% colnames(raw_counts))
all(colnames(raw_counts) %in% rownames(clindata))

# Make sure that the grouping variables are factors.
# Change the original names (long and complicated) to shorter and easier names
# this is to change a column in the data frame

colnames(clindata)[colnames(clindata) == "sequencing_batch:ch1"] <- "batch"  
clindata$batch <- as.factor (clindata$batch)

# this is to change a column in the data frame
colnames(clindata)[colnames(clindata) == 
                     "sars-cov-2 positivity:ch1"] <- "positivity"  

# this is to change the character in the column
clindata$positivity[clindata$positivity == "pos"] <- "COVID19"  
clindata$positivity[clindata$positivity == "neg"] <- "HEALTHY"
clindata$positivity <- as.factor(clindata$positivity)

# Merge the read counts and clinical data matrixes into a DESeqDataSet object 
# using the DESeq2 package
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = clindata,
                              design = ~ positivity + batch)

# Normalization by estimation of size factor
# raw_counts is equal to dds@assays@data@listData[["counts"]]
dds <- estimateSizeFactors(dds)

# Create a new table with the normalized read counts (gene expression)
norm_counts <- counts(dds, normalized=TRUE) #norm_counts : matrix, array

#cf. There are alternative methods to normalize and extract the counts 
# such as "rlog" and "vst"

# stratify viral_load
colnames(clindata)[colnames(clindata) == "n1_ct:ch1"] <- "ct"
clindata$viral_load <- clindata$ct # this is just duplication, not overwrite
clindata$viral_load[clindata$viral_load == "N/A"] <- "Negative"
clindata$viral_load[clindata$viral_load >24 & clindata$viral_load != "Unknown" 
                    & clindata$viral_load != "Negative"] <- "LOW"
clindata$viral_load[clindata$viral_load <= 24 & 
                      clindata$viral_load >=19] <- "MEDIUM"
clindata$viral_load[clindata$viral_load < 19] <- "HIGH"
clindata$viral_load <- as.factor(clindata$viral_load)
clindata$viral_load <- factor(clindata$viral_load, levels = c("NEGATIVE", 
                                                              "LOW", 
                                                              "MEDIUM", 
                                                              "HIGH", 
                                                              "UNKNOWN"))
clindata$positivity <- factor(clindata$positivity, levels =c("HEALTHY", 
                                                             "COVID19"))


# cf. to see the only two columns in many others in data frame
clindata[, c("ct", "viral_load")]
# cf. to see if two columns are same or not
all(clindata$ct == clindata$viral_load)


# stratify age

colnames(clindata)[colnames(clindata) == "age:ch1"] <- "age_cat"
clindata$age_cat[clindata$age_cat <30] = "<30"
clindata$age_cat[clindata$age_cat >=30 & clindata$age_cat <40] <- "30s"
clindata$age_cat[clindata$age_cat >=40 & clindata$age_cat <50] <- "40s"
clindata$age_cat[clindata$age_cat >=50 & clindata$age_cat <60] <- "50s"
clindata$age_cat[clindata$age_cat >=60 & clindata$age_cat <70] <- "60s"
clindata$age_cat[clindata$age_cat >=70] <- "70+"
clindata$age_cat[clindata$age_cat == "Unknown"] <- NA


# cf. make a new column and then delete it selectively
# make a new column without affect the old column
clindata$submission_Date <- clindata$submission_date 
clindata <- clindata %>% 
  select(-submission_Data)


# ggplot
# t(): is the transpose function, which converts the extracted row 
# into a row vector. 
# ,used when you want to switch the rows and columns of a matrix or data frame

MX1 <- ggplot(NULL, aes(x=clindata$positivity,     
                        y=log2(t(norm_counts["MX1",]+1)))) +
  geom_jitter(aes(shape=clindata$positivity, 
                  color=clindata$positivity), size=3)+
  xlab(NULL) +
  ylab("MX1 expression \n log2 (norm counts +1)") +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title =element_text(size = 25),
        legend.position = 'none') +
  stat_summary(fun=mean,
               geom="point",
               shape= '_',
               size=14,
               colour= c('#b53d35', '#066e70'))

MX1

# cf. to check if the "MX1 row in the norm_counts exists or not
any(rownames(norm_counts) == "MX1")
any(colnames(norm_counts) == "POS_002")

# Wilcoxon test to assess the statistical significance of MX1 expression 
MX1stat <- wilcox.test(norm_counts["MX1", ] ~ clindata$positivity,
                       paired = FALSE)
MX1stat                       

# jonckheere.terpstra (Arif et al., 2015) trend test to evaluate gene
# expression trends across ordered strata.

p_trend_age <- jonckheere.test(x=log2(t(norm_counts["MX1", ] +1))[
  clindata$positivity == "COVID19"],
  g=factor(clindata$age_cat[clindata$positivity == "COVID19"], 
           ordered = TRUE), 
  alternative = "decreasing", 
  nperm=500)


## Correlation analysis
# calculate the Spearman coefficients and plot all pairwise correlations

pairwise_corr <- ggpairs(as.data.frame(log2(t(norm_counts+1))),
                         columns = c("MX1", "MX2", "ACE2", "TMPRSS2"),
                         upper = list(continuous = wrap('cor', 
                                                        method = "spearman", 
                                                        size = 3),
                                      combo = "box_no_facet", 
                                      discrete = "count",
                                      na ="na"),
                         ggplot2::aes(colour=clindata$positivity, 
                                      shape=clindata$positivity, alpha = 0.01))

pairwise_corr <- pairwise_corr + theme(strip.placement = "outside",
                                       text = element_text(size = 9 , 
                                                           face = "bold")) +
  ggtitle("Gene correlation") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5)) +
  ylab("log2(counts +1)") + 
  xlab("log2 (counts +1)")

pairwise_corr


MX1_MX2 <- ggplot(NULL, aes(x = log2(t(norm_counts["MX1",]+1)[which(clindata$positivity=="COVID19" &  
                                                                      clindata$ct != "Unknown")]), y = log2(t(norm_counts["MX2",]+1))[which(clindata$positivity=="COVID19"& clindata$ct != "Unknown")], color = as.integer(clindata$ct[(which(clindata$positivity=="COVID19" & clindata$ct != "Unknown"))]))) +
  geom_point(size = 4, na.rm = TRUE) +
  scale_color_gradientn(colours=c("red","white","blue"), name = "Viral load (ct)") +
  ylab("MX2 expression RNA-seq \n log2 (norm counts +1)") +
  xlab("MX1 expression RNA-seq \n log2 (norm counts +1)") +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title =element_text(size = 25))

MX1_MX2 

MX1_MX2stat <-cor.test(norm_counts["MX1",]  
                       [which(clindata$positivity=="COVID19")],
                       norm_counts["MX2",]
                       [which(clindata$positivity=="COVID19")],
                       method = "spearman")
MX1_MX2stat

## Patient segregation based on gene expression


# Calculate the principal component on the log2 transformed gene expression data
# The output is a table as shown in Figure 4A, 
# containing the weight of each gene in the variance of the samples 
# for Principal Component 1 
# (PC1; the largest component of variance in the data set) and 
# Principal Component 2 
# (PC2; the second most important component influencing the variance):
  
# res.pca <- prcomp(t(log2(norm_counts[c("gene1","gene2", "geneN"),]+1)), 
#                   scale = TRUE)   

res.pca <- prcomp(t(log2(norm_counts[c("MX1","ACE2", "BSG"),]+1)), 
                  scale = TRUE)

#replace gene1, gene2, geneN by the list of genes of your interest

res.pca


# Plot PC1 _vs_. PC2 (Figure 4B)

p<- fviz_pca_biplot(res.pca, col.ind = clindata$positivity,
                    geom = "point",
                    addEllipses = TRUE,
                    palette = c('#F8766D', '#00BFC4'),
                    title='Principal Component Analysis')

p

# #D graphs showing expression for the selected genes

canvasXpress(
  data=t(log2(norm_counts[c("MX1","ACE2","BSG"),]+1)),
  varAnnot=as.data.frame(clindata$positivity, 
                         row.names=rownames(clindata)),
  axisTickScaleFontFactor=0.6,
  axisTitleScaleFontFactor=0.6,
  ellipseBy="clindata$positivity",
  colorBy="clindata$positivity",
  colorKey=list("clindata$positivity"=list("COVID19"="#F8766D", 
                                           "HEALTHY"="#00BFC4")),
  graphType="Scatter3D",
  title="3D scatter plot",
  xAxis=list("ACE2"),
  yAxis=list("BSG"),
  zAxis=list("MX1"),
  showLoessFit = FALSE)








