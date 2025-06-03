
# Coral Orthogroup RNAseq Analysis

> [!WARNING]
> This Github page is a work in progress. All of the current code is valid, however there is more to be added!

## Research Questions:

### 1. What is the proximal molecular mechanism that triggers Stony Coral Tissue Losss Disease? 
### 2. What is the molecular explanation to Resilience and Susceptibility to Stony Coral Tissue Loss Disease?

Question 1 is basically "what is sctld?". A question that has been nearly unanswered for any coral disease. The difficulty will be really honing the analysis to separate response mechanisms from what changed from a productive, life-sustaining existence to one that would result in holobiont collapse? Why did the mechanisms that prevent this suddenly fail now? The difficulty is that it enters a period where things change before reaching a noticeable tipping point. Its both really tricky but totally exciting. It will take introducing new methods, working towards a fine-scale pathogenesis and then the burden of challenging those conclusions until/if it breaks.
Question 2 is everything after question 1. The response, what was different, why would that be relevant to disease?

So Here is what SCTLD is:

In corals with a SCTLD lesion, there is a signficant increase in bacteria that were not present in control samples. So this is what changed from a happy coral to one that is now vulnerable. However, the bacteria that increased in abundance in disease-infected coral are not taxonomically identical, but we do see a lot of *functional* redundancy. Two different taxonomic bacteria signficantly expressing molecular mechanisms that are indrectly antagonistic or directly virulent. One of these early molecular switches is the activation of mobile genetic elements. Mobile genetic elements are typically silenced because they remodel the bacterias genome to rapidly adapt and this results in the expression of potentially antagonistic processes that help the bacteria survive, even if its at the expense of surrounding life. So seeing this activated is a clear mechanisms of something that changes. We also see funcitonal redundancy of other indrect antagonistic mechanisms such as quorum sensing, evasion of host detection, antioxidant production and breakdown, membrane strengthening, and direct mechanisms of virulence sucha as iron virulence. From here we get into the rest of the pathogenesis, the response to these actions. We also see the host responding to these mechanisms, supported through genes that are signficantly correlated to this bacteria activity and those particular genes would be involved in the processes that would biologigcally be responding. So this begins to build the case to address Question 1 and I will try to find where these conclusions break and are inconsistent and if they dont then it builds support for the conclusion. 


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host")
knitr::opts_chunk$set(dev='ragg_png')
library(ggpubr)
library(DESeq2)
library(tibble)
library(tidyverse)
library(tximport)
#library(GOplot)
library(dplyr)
library(tidyr)
library(reshape2)
#library(ggalt)
library(ggplot2)
library(lifecycle)
library(stringr)
library(utils)
library(corrplot)
#library(repr) 
library(data.table)
```

### Step 1: Obtain Shared Orthogroups from Orthofinder

To obtain Orthogroups, we filter our Orthofinder results to only contain Orthogroups with at least one sequence per species. 

This can be done in EXCEL like so:
1. Open Orthogroups.tsv in Excel from Orthofinder results
2. Click *Find & Select*
3. Click *Go to Special*
4. Choose *Blanks*
5. Click OK and then all the *blank rows/cells will be highlighted*
6. Choose the *Delete under Cells* section on the Home Tab
7. Click *Delete Sheet Rows*
8. Steps 1-8 were repeated 4 times so that all blank rows were eventually removed. I think the amount of data excel needed to process required these steps to be repeated. You can ensure you have the accurate final number of orthogroups by comparing your row count to the "Statistics_Overall" file in /Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Comparative_Genomics_Statistics/ specifically the value in row "Number of orthogroups with all species present" 
9. Save as "host_shared_orthogroups.csv"

### Make tx2gene files using shared orthogroups
> This is taking the shared orthogroups file and separating it by the columns of transcipts to their respective species. For example one of these columns would be called "Acer_Host_reference_proteome_AllORF_SingleBestOnly".
```{r Host Orthogroups}

Orthogroups <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/host_shared_orthogroups_SingleBestORF.csv")

# Acer Orthogroups
Acer_orthogroups <- Orthogroups[,c(1,2)]
Acer_orthogroups <- Acer_orthogroups %>% separate_rows(Acer_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Acer_orthogroups <- Acer_orthogroups[,c(2,1)]
colnames(Acer_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Acer_host_Orthogroup_Transcripts <- Acer_orthogroups[c(1)]
names(Acer_host_Orthogroup_Transcripts)
write.table(Acer_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Acer_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Acer_orthogroups$Transcript <- gsub("\\..*","",Acer_orthogroups$Transcript)
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Acerhost_orthogroup_tx2gene.csv",row.names = FALSE)

# Mcav Orthogroups
Mcav_orthogroups <- Orthogroups[,c(1,3)]
Mcav_orthogroups <- Mcav_orthogroups %>% separate_rows(Mcav_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Mcav_orthogroups <- Mcav_orthogroups[,c(2,1)]
colnames(Mcav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Mcav_host_Orthogroup_Transcripts <- Mcav_orthogroups[c(1)]
names(Mcav_host_Orthogroup_Transcripts)
write.table(Mcav_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Mcav_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Mcav_orthogroups$Transcript <- gsub("\\..*","",Mcav_orthogroups$Transcript)
Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Mcavhost_orthogroup_tx2gene.csv",row.names = FALSE)

# Ofav Orthogroups
Ofav_orthogroups <- Orthogroups[,c(1,4)]
Ofav_orthogroups <- Ofav_orthogroups %>% separate_rows(Ofav_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Ofav_orthogroups <- Ofav_orthogroups[,c(2,1)]
colnames(Ofav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Ofav_host_Orthogroup_Transcripts <- Ofav_orthogroups[c(1)]
names(Ofav_host_Orthogroup_Transcripts)
write.table(Ofav_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Ofav_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Ofav_orthogroups$Transcript <- gsub("\\..*","",Ofav_orthogroups$Transcript)
Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Ofavhost_orthogroup_tx2gene.csv",row.names = FALSE)

# Past Orthogroups
Past_orthogroups <- Orthogroups[,c(1,5)]
Past_orthogroups <- Past_orthogroups %>% separate_rows(Past_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Past_orthogroups <- Past_orthogroups[,c(2,1)]
colnames(Past_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Past_host_Orthogroup_Transcripts <- Past_orthogroups[c(1)]
names(Past_host_Orthogroup_Transcripts)
write.table(Past_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Past_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Past_orthogroups$Transcript <- gsub("\\..*","",Past_orthogroups$Transcript)
Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Pasthost_orthogroup_tx2gene.csv",row.names = FALSE)

```


### Count Matrix Generation

Read in the raw counts from Salmon for all coral samples and assign Orthogroups
In terminal
```{r}
# Acer
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Acer

# Mcav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Mcav

# Ofav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Ofav

# Past
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Past

```

### Combine Orthogroups with Salmon counts.
> [!NOTE]
> Make sure Acer_samples.csv is the path to your salmon quants.
>
> Here is what Acer_samples.csv looks like:

![image](https://github.com/user-attachments/assets/7b2eb069-01ef-45fd-80a2-69784a6fd35f)

> in those paths, there will be a file called "quant.sf" for each sample that these paths will be extracting.
>

 

### A. cervicornis
```{r,results='hide',tidy=TRUE}
setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/host/Acer/")
Acer_samples <- read.csv("~/Desktop/Microbial\ Metatranscriptomics/salmon/host/Acer/Acer_samples.csv", header = TRUE)
Acer_files <- file.path(Acer_samples$sample, "quant.sf")
names(Acer_files) <- paste0("sample", 1:13)
all(file.exists(Acer_files))
# Do not continue if output shows [FALSE]
Acer_tx2gene <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Acer_orthogroup_tx2gene_SingleBestORF.csv")
Acer_txi <- tximport(Acer_files, type = 'salmon',tx2gene=Acer_tx2gene)
dim(Acer_txi$counts)

# 6,955 Orthogroups - SingleBestORF

Acer_counts <- Acer_txi$counts
Acer_counts <- as.data.frame(Acer_counts)

# Extract sample names from file paths
sample_names <- Acer_samples$sample
sample_names_cleaned <- basename(sample_names)

# Use the cleaned sample names to rename the columns of Acer_counts
colnames(Acer_counts) <- c(sample_names_cleaned)

Acer_counts <- tibble::rownames_to_column(Acer_counts,"Orthogroup")
Acer_counts[2:14] <-  round(Acer_counts[2:14], digits=0)
Acer_counts <- as.data.frame(Acer_counts)
```

### M. cavernosa
```{r,results='hide',tidy=TRUE}
setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/host/Mcav/")
Mcav_samples <- read.csv("~/Desktop/Microbial Metatranscriptomics/salmon/host/Mcav/Mcav_samples.csv")
Mcav_files <- file.path(Mcav_samples$sample, "quant.sf")
names(Mcav_files) <- paste0("sample", 1:18)
all(file.exists(Mcav_files))
# Do not continue if output shows [FALSE]
Mcav_tx2gene <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Mcav_orthogroup_tx2gene_SingleBestORF.csv")
Mcav_txi <- tximport(Mcav_files, type = 'salmon',tx2gene=Mcav_tx2gene)
dim(Mcav_txi$counts)

# 6,955 Orthogroups - SingleBestORF

Mcav_counts <- Mcav_txi$counts
Mcav_counts <- as.data.frame(Mcav_counts)

# Extract sample names from file paths
sample_names <- Mcav_samples$sample
sample_names_cleaned <- basename(sample_names)

# Use the cleaned sample names to rename the columns of Mcav_counts
colnames(Mcav_counts) <- c(sample_names_cleaned)

Mcav_counts <- tibble::rownames_to_column(Mcav_counts,"Orthogroup")
Mcav_counts[2:19] <-  round(Mcav_counts[2:19], digits=0)
Mcav_counts <- as.data.frame(Mcav_counts)

```

### O. faveolata
```{r,results='hide',tidy=TRUE}
setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/host/Ofav/")
Ofav_samples <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/host/Ofav/Ofav_samples.csv")
Ofav_files <- file.path(Ofav_samples$sample, "quant.sf")
names(Ofav_files) <- paste0("sample", 1:29)
all(file.exists(Ofav_files))
# Do not continue if output shows [FALSE]
Ofav_tx2gene <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Ofav_orthogroup_tx2gene_SingleBestORF.csv")
Ofav_txi <- tximport(Ofav_files, type = 'salmon',tx2gene=Ofav_tx2gene)

dim(Ofav_txi$counts)

# 6,955 Orthogroups - SingleBestORF

Ofav_counts <- Ofav_txi$counts
Ofav_counts <- as.data.frame(Ofav_counts)

# Extract sample names from file paths
sample_names <- Ofav_samples$sample
sample_names_cleaned <- basename(sample_names)

# Use the cleaned sample names to rename the columns of Ofav_counts
colnames(Ofav_counts) <- c(sample_names_cleaned)

Ofav_counts <- tibble::rownames_to_column(Ofav_counts,"Orthogroup")
Ofav_counts[2:30] <-  round(Ofav_counts[2:30], digits=0)
Ofav_counts <- as.data.frame(Ofav_counts)

```

### P. astreoides
```{r,results='hide',tidy=TRUE}
setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/host/Past/")
Past_samples <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/host/Past/Past_samples.csv")
Past_files <- file.path(Past_samples$sample, "quant.sf")
names(Past_files) <- paste0("sample", 1:17)
all(file.exists(Past_files))
# Do not continue if output shows [FALSE]
Past_tx2gene <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Past_orthogroup_tx2gene_SingleBestORF.csv")
Past_txi <- tximport(Past_files, type = 'salmon',tx2gene=Past_tx2gene)
dim(Past_txi$counts)

# 6,955 Orthogroups - SingleBestORF

Past_counts <- Past_txi$counts
Past_counts <- as.data.frame(Past_counts)

# Extract sample names from file paths
sample_names <- Past_samples$sample
sample_names_cleaned <- basename(sample_names)

# Use the cleaned sample names to rename the columns of Past_counts
colnames(Past_counts) <- c(sample_names_cleaned)

Past_counts <- tibble::rownames_to_column(Past_counts,"Orthogroup")
Past_counts[2:18] <-  round(Past_counts[2:18], digits=0)
Past_counts <- as.data.frame(Past_counts)
```


## Analysis of Shared Orthogroups
### Merge the count data
```{r,results='hide',tidy=TRUE}
Acer_Mcav <- merge(Acer_counts,Mcav_counts, by ="Orthogroup")
Acer_Mcav_Past <- merge(Acer_Mcav,Past_counts,by="Orthogroup")
Acer_Mcav_Past_Ofav <- merge(Acer_Mcav_Past,Ofav_counts,by="Orthogroup")
#CladeA_CladeB_CladeC_CladeD <- merge(CladeA_CladeB_CladeC,CladeD_counts,by="Orthogroup")
#all_species <- merge(CladeA_CladeB_CladeC_CladeD,pstr_counts,by="Orthogroup")
all_species <- Acer_Mcav_Past_Ofav

# 6,955 Orthogroups SingleBestORF
```

### Formatting Sample Names to make datasets agreeable
```{r,results='hide',tidy=TRUE}
countData <- all_species
countData <- tibble::column_to_rownames(countData,"Orthogroup")
colData <- read.csv("/Users/nicholas.macknight/Desktop/Autumn16s/Autumn16s_05082023/230613_M02476_0578_000000000-L3NTN/20230615_222347/metadata.csv")

# So countData and colData have some slight name formatting differences that we need to make consistent.

# Remove "_host_quant" from the end of the column names
colnames(countData) <- gsub("_host_quant$", "", colnames(countData))

# Step 1: Remove the first two underscores in colData[[1]]
#colData[[1]] <- gsub("^(.*?)_(.*?)_(.*?)_", "\\1", colData[[1]])
colData[[1]] <- gsub("_", "", colData[[1]])  # Remove everything after the underscore
colData[[1]] <- gsub("Control", "Control_", colData[[1]])  # Remove everything after the underscore
colData[[1]] <- gsub("Disease", "Disease_", colData[[1]])  # Remove everything after the underscore

colData[[1]]

# Step 2: Remove everything after the first underscore and replace dashes with underscores in colnames(countData)
colnames(countData) <- gsub("_(.*)", "", colnames(countData))  # Remove everything after the underscore
colnames(countData) <- gsub("-", "_", colnames(countData))     # Replace dashes with underscores
colnames(countData) <- gsub("17_7", "17-7", colnames(countData))  # Remove everything after the underscore

# Step 3: Subset colData to only include rows where the first column matches the colnames in countData
colData_subset <- colData[colData[[1]] %in% colnames(countData), ]

#countData_subset <- colnames(countData) %in% colData_subset[colData_subset[[1]], ]

# View the subsetted colData
colData_subset
colData_subset[[1]]

# Remove the column named "PastPA3Disease_5_2" from countData
countData <- countData[, !colnames(countData) %in% "PastPA3Disease_5_2"]
```


# Deseq2
Differential Abundace

### ~ No Design
> I use this no design for EVE input because I do not want to be controlling for coral.species or treatment affects as this would likely skew the eve model. I simply want rlog normalized data for eve input.
```{r}
# Create a DESeqDataSet with a placeholder design
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData_subset, 
                              design = ~ 1) # No design formula

# Filter out rows with low counts
dds <- dds[rowMeans(counts(dds)) > 10, ]

# Perform rlog normalization
rld <- rlog(dds, blind = TRUE)

# Extract the rlog-transformed data
rlog_data <- assay(rld)

# Transpose rlog_data
rlog_data_t <- as.data.frame(t(rlog_data))

# Add a column for sample names to match with colData_subset
rlog_data_t$SampleID <- rownames(rlog_data_t)

# Merge the datasets based on the SampleID column
rlog_metadata_nodesign <- merge(colData_subset, rlog_data_t, by = "SampleID")

# Save the merged dataset if needed
write.csv(rlog_metadata_nodesign, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/rlog_metadata_nodesign_singlebestorf.csv", row.names = FALSE)

# Display a preview
head(rlog_metadata_nodesign)

```

### ~Coral.Species+Treatment
> DESeq2 will control for differences between coral species before assessing the effect of Treatment by specifying Coral.Species first in the design ~Coral.Species+Treatment.
```{r,results='hide',tidy=TRUE}
# Set the design
dds_treatment <- DESeqDataSetFromMatrix(countData = countData, 
                                        colData = colData_subset, 
                                        design = ~Coral.Species+Treatment)

# Filter out columns with an average number of reads less than 10
dds_treatment <- dds_treatment[ rowMeans(counts(dds_treatment)) > 10, ] 

# Run DESeq2
test <- DESeq(dds_treatment)
rld <- rlog(dds_treatment, blind=FALSE) 
rrld <- assay(rld)

write.csv(rrld, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/treatment_rlog_orthogroups_singlebestorf.csv")  

# Extract the results for the comparison of interest
# For example, comparing "TreatmentA" to "TreatmentB":
results_treatment <- results(test, contrast = c("Treatment", "Control", "Disease"))

# Order results by padj
results_treatment <- results_treatment[order(results_treatment$padj), ]

# Exclude rows with NA in padj before filtering for significance
significant_results <- results_treatment[!is.na(results_treatment$padj) & results_treatment$padj < 0.1, ]

# Save the significant results to a CSV file
 write.csv(as.data.frame(significant_results), file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/significant_treatment_comparison_singlebestorf.csv")

# View summary of significant results
summary(results_treatment)

```

out of 713 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 25, 3.5%
LFC < 0 (down)     : 84, 12%
outliers [1]       : 0, 0%
low counts [2]     : 83, 12%
(mean count < 55)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

### ~Treatment+Coral.Species
```{r,results='hide',tidy=TRUE}
colData_subset$Coral.Species <- as.factor(colData_subset$Coral.Species)
colData_subset$Treatment <- as.factor(colData_subset$Treatment)


# Set the design
dds_Species <- DESeqDataSetFromMatrix(countData = countData, 
                                        colData = colData_subset, 
                                        design = ~Treatment + Coral.Species)

# Filter out rows with an average number of reads less than 10
dds_Species <- dds_Species[ rowMeans(counts(dds_Species)) > 10, ] 

# Run DESeq2
test_Species <- DESeq(dds_Species)

# Perform rlog normalization
rld_Species <- rlog(test_Species, blind = FALSE)

# Save rlog-transformed data for downstream analysis
rlog_data <- assay(rld_Species)
#write.csv(rlog_data, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/rlog_normalized_counts_BySpecies_BestORF.csv")

# Perform pairwise comparisons for Coral.Species
res_Acer_vs_Mcav <- results(test_Species, contrast = c("Coral.Species", "Acer", "Mcav"))
res_Acer_vs_Ofav <- results(test_Species, contrast = c("Coral.Species", "Acer", "Ofav"))
res_Acer_vs_Past <- results(test_Species, contrast = c("Coral.Species", "Acer", "Past"))
res_Mcav_vs_Ofav <- results(test_Species, contrast = c("Coral.Species", "Mcav", "Ofav"))
res_Mcav_vs_Past <- results(test_Species, contrast = c("Coral.Species", "Mcav", "Past"))
res_Ofav_vs_Past <- results(test_Species, contrast = c("Coral.Species", "Ofav", "Past"))

summary(res_Acer_vs_Mcav)
summary(res_Acer_vs_Ofav)
summary(res_Acer_vs_Past)
summary(res_Mcav_vs_Ofav)
summary(res_Mcav_vs_Past)
summary(res_Ofav_vs_Past)

```
summary(res_Acer_vs_Mcav)
out of 6944 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1953, 28%
LFC < 0 (down)     : 2665, 38%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

summary(res_Acer_vs_Ofav)
out of 6944 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1786, 26%
LFC < 0 (down)     : 2838, 41%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

summary(res_Acer_vs_Past)
out of 6944 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 193, 2.8%
LFC < 0 (down)     : 2161, 31%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

summary(res_Mcav_vs_Ofav)
out of 6944 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1864, 27%
LFC < 0 (down)     : 2265, 33%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

summary(res_Mcav_vs_Past)
out of 6944 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 2152, 31%
LFC < 0 (down)     : 2602, 37%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

summary(res_Ofav_vs_Past)
out of 6944 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 2234, 32%
LFC < 0 (down)     : 2556, 37%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


### Running DESeq2 on Each Species
```{r}
# Load necessary libraries
library(DESeq2)
library(VennDiagram)
library(ggplot2)
library(tibble)

# List of species datasets
species_datasets <- list("Acer" = Acer_counts, 
                         "Mcav" = Mcav_counts,
                         "Ofav" = Ofav_counts,
                         "Past" = Past_counts)

# Prepare colData for all species (assuming same metadata structure for each)
colData <- read.csv("/Users/nicholas.macknight/Desktop/Autumn16s/Autumn16s_05082023/230613_M02476_0578_000000000-L3NTN/20230615_222347/metadata.csv")

# Remove the first two underscores in colData, keep both lines, they are necessary.
colData[[1]] <- sub("_", "", colData[[1]])
colData[[1]] <- sub("_", "", colData[[1]])

# Create an empty list to store significant orthogroups for each species
significant_orthogroups <- list()

# Iterate through each species, run DESeq2, and collect significant orthogroups
for(species in names(species_datasets)) {
  
  # Prepare count data for the species
  countData <- species_datasets[[species]]
  countData <- tibble::column_to_rownames(countData, "Orthogroup")
  
  # Clean column names (if necessary)
  colnames(countData) <- gsub("_host_quant$", "", colnames(countData))
  colnames(countData) <- gsub("_(.*)", "", colnames(countData))
  colnames(countData) <- gsub("-", "_", colnames(countData))
  colnames(countData) <- gsub("17_7", "17-7", colnames(countData))  # Remove everything after the underscore

  # Subset colData for the current species
  colData_subset <- colData[colData[[1]] %in% colnames(countData), ]

  # Check if Coral.Genotype has more than one unique value for the current species
  if(length(unique(colData_subset$Coral.Genotype)) > 1) {
    # Include Coral.Genotype in the design if it has more than one level
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData_subset,
                                  design = ~Coral.Genotype + Treatment)
  } else {
    # Use Treatment only if Coral.Genotype has only one level
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData_subset,
                                  design = ~Treatment)
  }
  
  # Filter rows with low counts
  dds <- dds[rowMeans(counts(dds)) > 10, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results for the comparison of interest (e.g., Control vs Disease)
  res <- results(dds, contrast = c("Treatment", "Control", "Disease"))
  
  # Filter for significant orthogroups (adjust based on your threshold)
  sig_orthogroups <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 0.1)]
  
  # Store the significant orthogroups for this species
  significant_orthogroups[[species]] <- sig_orthogroups
}

# Print significant orthogroups for each species
print(significant_orthogroups)


```

$Acer
  [1] "OG0000002" "OG0000046" "OG0000095" "OG0000115" "OG0000142" "OG0000183" "OG0000195" "OG0000223" "OG0000248" "OG0000314" "OG0000337" "OG0000378" "OG0000408" "OG0000413" "OG0000432" "OG0000435" "OG0000454" "OG0000457" "OG0000515" "OG0000535" "OG0000560" "OG0000571"
 [23] "OG0000593" "OG0000596" "OG0000609" "OG0000648" "OG0000665" "OG0000706" "OG0000716" "OG0000781" "OG0000786" "OG0000792" "OG0000808" "OG0000818" "OG0000830" "OG0000833" "OG0000904" "OG0000905" "OG0000907" "OG0000912" "OG0000913" "OG0000926" "OG0000938" "OG0000950"
 [45] "OG0000960" "OG0000962" "OG0000972" "OG0000975" "OG0000981" "OG0000988" "OG0001054" "OG0001055" "OG0001065" "OG0001069" "OG0001081" "OG0001088" "OG0001110" "OG0001118" "OG0001133" "OG0001140" "OG0001166" "OG0001170" "OG0001181" "OG0001188" "OG0001192" "OG0001235"
 [67] "OG0001241" "OG0001250" "OG0001272" "OG0001303" "OG0001347" "OG0001373" "OG0001397" "OG0001404" "OG0001412" "OG0001413" "OG0001496" "OG0001519" "OG0001583" "OG0001596" "OG0001599" "OG0001636" "OG0001646" "OG0001648" "OG0001664" "OG0001666" "OG0001668" "OG0001717"
 [89] "OG0001727" "OG0001742" "OG0001749" "OG0001753" "OG0001852" "OG0001853" "OG0001903" "OG0001907" "OG0001940" "OG0001954" "OG0001979" "OG0001995" "OG0001996" "OG0002062" "OG0002081" "OG0002095" "OG0002099" "OG0002104" "OG0002105" "OG0002138" "OG0002160" "OG0002165"
[111] "OG0002167" "OG0002192" "OG0002193" "OG0002194" "OG0002200" "OG0002343" "OG0002374" "OG0002378" "OG0002423" "OG0002430" "OG0002458" "OG0002459" "OG0002461" "OG0002495" "OG0002504" "OG0002539" "OG0002540" "OG0002546" "OG0002560" "OG0002565" "OG0002567" "OG0002570"
[133] "OG0002571" "OG0002598" "OG0002602" "OG0002614" "OG0002671" "OG0002686" "OG0002688" "OG0002694" "OG0002697" "OG0002717" "OG0002740" "OG0002747" "OG0002764" "OG0002767" "OG0002771" "OG0002775" "OG0002984" "OG0003017" "OG0003023" "OG0003035" "OG0003036" "OG0003073"
[155] "OG0003092" "OG0003096" "OG0003101" "OG0003102" "OG0003115" "OG0003126" "OG0003130" "OG0003141" "OG0003145" "OG0003174" "OG0003181" "OG0003210" "OG0003212" "OG0003264" "OG0003302" "OG0003311" "OG0003335" "OG0003364" "OG0003382" "OG0003394" "OG0003406" "OG0003417"
[177] "OG0003425" "OG0003459" "OG0003463" "OG0003496" "OG0003520" "OG0003541" "OG0003544" "OG0003552" "OG0003568" "OG0003573" "OG0003581" "OG0003590" "OG0003591" "OG0003622" "OG0003626" "OG0003636" "OG0003639" "OG0003640" "OG0003691" "OG0003711" "OG0003715" "OG0003727"
[199] "OG0004005" "OG0004015" "OG0004029" "OG0004048" "OG0004061" "OG0004088" "OG0004095" "OG0004102" "OG0004128" "OG0004141" "OG0004198" "OG0004205" "OG0004213" "OG0004228" "OG0004239" "OG0004241" "OG0004249" "OG0004256" "OG0004257" "OG0004271" "OG0004282" "OG0004311"
[221] "OG0004331" "OG0004380" "OG0004390" "OG0004423" "OG0004448" "OG0004456" "OG0004461" "OG0004528" "OG0004546" "OG0004560" "OG0004563" "OG0004585" "OG0004591" "OG0004650" "OG0004658" "OG0004677" "OG0004692" "OG0004697" "OG0004717" "OG0004768" "OG0004786" "OG0004788"
[243] "OG0004793" "OG0004795" "OG0004837" "OG0004842" "OG0004857" "OG0004876" "OG0004901" "OG0004916" "OG0004920" "OG0004933" "OG0004981" "OG0004983" "OG0005003" "OG0005025" "OG0005053" "OG0005077" "OG0005100" "OG0005148" "OG0005158" "OG0005176" "OG0005194" "OG0005202"
[265] "OG0005217" "OG0005218" "OG0005669" "OG0005680" "OG0005684" "OG0005721" "OG0005750" "OG0005758" "OG0005761" "OG0005814" "OG0005839" "OG0005847" "OG0005861" "OG0005880" "OG0005881" "OG0005908" "OG0005911" "OG0005955" "OG0005990" "OG0006014" "OG0006028" "OG0006040"
[287] "OG0006064" "OG0006096" "OG0006156" "OG0006196" "OG0006203" "OG0006217" "OG0006252" "OG0006262" "OG0006284" "OG0006298" "OG0006332" "OG0006339" "OG0006370" "OG0006412" "OG0006435" "OG0006515" "OG0006522" "OG0006541" "OG0006600" "OG0006611" "OG0006614" "OG0006619"
[309] "OG0006661" "OG0006683" "OG0006702" "OG0006707" "OG0006727" "OG0006756" "OG0006766" "OG0006789" "OG0006833" "OG0006849" "OG0006871" "OG0006877" "OG0006896" "OG0006897" "OG0006960" "OG0006964" "OG0006974" "OG0006979" "OG0007065" "OG0007078" "OG0007085" "OG0007101"
[331] "OG0007122" "OG0007124" "OG0007137" "OG0007159" "OG0007167" "OG0007199" "OG0007282" "OG0007302" "OG0007310" "OG0007359" "OG0007362" "OG0007365" "OG0007401" "OG0007407" "OG0007411" "OG0007415" "OG0007423" "OG0007444" "OG0007471" "OG0007474" "OG0007478" "OG0007480"
[353] "OG0007492" "OG0007506" "OG0007507" "OG0007514" "OG0007530" "OG0007542" "OG0007556" "OG0007571" "OG0007605" "OG0007632" "OG0007634" "OG0007647" "OG0007658" "OG0007660" "OG0007668" "OG0007671" "OG0007675" "OG0007676" "OG0008340" "OG0008395" "OG0008462" "OG0008557"
[375] "OG0008597" "OG0008624" "OG0008630" "OG0008633" "OG0008692" "OG0008694" "OG0008752" "OG0008771" "OG0008809" "OG0008832" "OG0008839" "OG0008851" "OG0008857" "OG0008880" "OG0008887" "OG0008923" "OG0009001" "OG0009048" "OG0009058" "OG0009082" "OG0009121" "OG0009134"
[397] "OG0009142" "OG0009169" "OG0009196" "OG0009206" "OG0009209" "OG0009210" "OG0009221" "OG0009324" "OG0009335" "OG0009378" "OG0009381" "OG0009407" "OG0009422" "OG0009423" "OG0009424" "OG0009434" "OG0009454" "OG0009478" "OG0009500" "OG0009528" "OG0009531" "OG0009532"
[419] "OG0009533" "OG0009535" "OG0009559" "OG0009575" "OG0009645" "OG0009661" "OG0009697" "OG0009698" "OG0009745" "OG0009795" "OG0009804" "OG0009808" "OG0009894" "OG0009927" "OG0009941" "OG0009942" "OG0010014" "OG0010023" "OG0010052" "OG0010063" "OG0010091" "OG0010128"
[441] "OG0010206" "OG0010264" "OG0010266" "OG0010267" "OG0010268" "OG0010288" "OG0010290" "OG0010343" "OG0010356" "OG0010371" "OG0010386" "OG0010393" "OG0010399" "OG0010404" "OG0010425" "OG0010443" "OG0010465" "OG0010469" "OG0010479" "OG0010496" "OG0010515" "OG0010535"
[463] "OG0010546" "OG0010563" "OG0010645" "OG0010668" "OG0010721" "OG0010740" "OG0010747" "OG0010753" "OG0010758" "OG0010835" "OG0010874" "OG0010890" "OG0010894" "OG0010904" "OG0010927" "OG0010956" "OG0010969" "OG0010973" "OG0010981" "OG0011007" "OG0011017" "OG0011052"
[485] "OG0011057" "OG0011071" "OG0011085" "OG0011106" "OG0011118" "OG0011129" "OG0011157" "OG0011165" "OG0011169" "OG0011187" "OG0011218" "OG0011227" "OG0011247" "OG0011259" "OG0011278" "OG0011304" "OG0011331" "OG0011345" "OG0011385" "OG0011396" "OG0011410" "OG0011417"
[507] "OG0011427" "OG0011443" "OG0011452" "OG0011482" "OG0011492" "OG0011496" "OG0011503" "OG0011531" "OG0011539" "OG0011574" "OG0011642" "OG0011663" "OG0011676" "OG0011683" "OG0011689" "OG0011698" "OG0011719" "OG0011743" "OG0011767" "OG0011770" "OG0011774" "OG0011801"

$Mcav
 [1] "OG0000716" "OG0001160" "OG0001384" "OG0001713" "OG0003111" "OG0004481" "OG0004641" "OG0005774" "OG0007642" "OG0011530"

$Ofav
 [1] "OG0000142" "OG0000277" "OG0000389" "OG0000594" "OG0000600" "OG0000639" "OG0000782" "OG0000805" "OG0001540" "OG0001939" "OG0001970" "OG0002504" "OG0003177" "OG0005113" "OG0006222" "OG0007474" "OG0007579" "OG0007609" "OG0007614" "OG0007670" "OG0008479" "OG0008913"
[23] "OG0009113" "OG0009169" "OG0010803" "OG0010950" "OG0011375"

$Past
 [1] "OG0000593" "OG0000688" "OG0000711" "OG0000845" "OG0001358" "OG0001729" "OG0002002" "OG0002547" "OG0003243" "OG0003260" "OG0003729" "OG0004127" "OG0005036" "OG0005697" "OG0005954" "OG0006297" "OG0006915" "OG0007291" "OG0007510" "OG0007580" "OG0008388" "OG0008773"
[23] "OG0008910" "OG0009455" "OG0009691" "OG0009702" "OG0010250" "OG0010402" "OG0010452" "OG0010741" "OG0011133" "OG0011707"

### Venn Diagram
```{r}
# Install the ggvenn package if needed
#install.packages("ggvenn")

# Load ggvenn
library(ggvenn)

# Create the Venn diagram using ggvenn
ggvenn_plot <- ggvenn(
  significant_orthogroups,
  fill_color = c("red", "blue", "green","yellow"),
  text_size = 5
)

# Print the plot in R Markdown
print(ggvenn_plot)
```

![image](https://github.com/user-attachments/assets/cab38579-8b41-4498-87b5-4828141bf633)

# EVE

```{r EVE_phyloTree,fig.align='center',fig.height=1.5,fig.width=1.5,dpi=100}
library(ape)

# Trialing out rlog no design data in eve
NodesignTbl <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/rlog_metadata_nodesign.csv")
NodesignTbl.t <- t(NodesignTbl)
colnames(NodesignTbl.t) <- NodesignTbl.t[1,]
NodesignTbl.t <- NodesignTbl.t[-(1:10), ]
colnames(NodesignTbl.t)
# Convert row names to a column
# Convert the dataset to numeric
NodesignTbl.t <- as.data.frame(NodesignTbl.t)

NodesignTbl.t <- cbind(Orthogroup = rownames(NodesignTbl.t), NodesignTbl.t)
# Convert all other columns to numeric
NodesignTbl.t[, -1] <- lapply(NodesignTbl.t[, -1], function(x) as.numeric(as.character(x)))
NodesignTbl.t <- as.matrix(NodesignTbl.t[,-1])
dim(NodesignTbl.t)
# 9376 x 76 ORF
# 6944 x 75 ORF

exprTbl <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/treatment_rlog_orthogroups_LongestORF.csv")
colnames(exprTbl)[colnames(exprTbl)=="X"] <- "Orthogroup"

# The table needs to be converted to a matrix:
exprMat <- as.matrix(exprTbl[,-1])
rownames(exprMat) <- exprTbl$Orthogroup
exprMat
dim(exprMat)  
#Should return 1766 x 52 - KB Orthologs
# Should return 7060 x 76 - NM OrthoGroups
# Should return 5485 x 76 - NM LongestORF Orthogroups

# Species phylogeny can be obtained from OrthoFinder output in DIR=/OrthoFinder/Results/Species_Tree/SpeciesTree_rooted.txt

# I find it is best to run OrthoFinder with an outgroup species to get a better species tree. I included a Cassiopea proteome (generated the same way as the other proteomes) as my outgroup to obtain the following rooted species tree:
# (Cassiopea_proteome:0.405291,(past_reference_proteome:0.207024,((mcav_reference_proteome:0.0938366,oann_reference_proteome:0.0723831)0.262775:0.0844971,(pstr_reference_proteome:0.0954775,cnat_reference_proteome:0.106944)0.514097:0.105377)0.45022:0.114876)1:0.405291); 
# I then remove the out group from the rooted tree and abbreviate the species names
#speciesTree <-read.tree(text="((Acer:0.384136,Past:0.335512)0.691517:0.126116,(Mcav:0.15129,Ofav:0.220088)0.691517:0.126116);
#;")

speciesTree <-read.tree(text="((Mcav:0.19085,Ofav:0.23287)0.637884:0.139991,(Acer:0.426372,Past:0.385854)0.637884:0.139991);")
speciesTree <-read.tree(text="(Acer:0.201967,(Past:0.348625,(Mcav:0.171507,Ofav:0.232107)0.669734:0.262447)1:0.201967);") # SingleBestORF


# Check validity of species tree
plot(speciesTree,
     type = "phylogram",
     use.edge.length = TRUE,
     node.pos = NULL,
     show.tip.label = TRUE,
     root.edge = TRUE)
```

![image](https://github.com/user-attachments/assets/89c4cdef-2251-4f44-88fe-ea65e563c450)


EVE needs to know which species each column in the expression matrix belongs to. This is done by creating a vector of species names corresponding to the tip labels in the species tree

```{r,results='hide',tidy=TRUE}
speciesTree$tip.label
colnames(NodesignTbl.t)
colSpecies <- substr(colnames(NodesignTbl.t), 1, 4)
colSpecies

library(evemodel)
# Run EVE
EVE.NoDesign_SingleBestORF <- betaSharedTest(tree = speciesTree, gene.data = NodesignTbl.t, colSpecies = colSpecies)
```

> [!IMPORTANT]
> The sharedBeta value ( which is 0.2839223) is an integral value in the EVE results. We will log transform it  for interpretation use with the rest of the eve results.
> The sharedBeta is the cutoff between classifying genes as lineage specific (low beta) or highly variable (high beta), but the p-value is required to classify it as significantly either of those categories.
```{r,results='hide',tidy=TRUE}
# P-value can then be calculated using:
pval = pchisq(EVE.NoDesign_SingleBestORF$LRT,df = 1,lower.tail = F)

# The shared beta:
EVE.NoDesign_SingleBestORF$sharedBeta
# [1] 0.2839223 EVE No Design - SingleBestORF

log(0.2839223)
# [1] -1.259055

#Combine LRT, beta, theta, sigma2, alpha:
head(cbind(EVE.NoDesign_SingleBestORF$LRT,EVE.NoDesign_SingleBestORF$indivBetaRes$par,pval))
colnames(NodesignTbl.t)
dim(NodesignTbl.t)
# [1] 6944   75

a <- cbind(EVE.NoDesign_SingleBestORF$LRT,EVE.NoDesign_SingleBestORF$indivBetaRes$par,pval,NodesignTbl.t)
rownames(a) <- rownames(NodesignTbl.t)

# Give Column Name to LRT column
colnames(a)
colnames(a)[colnames(a)==""] <- "LRT"
head(a)
EVE_results <- as.data.frame(a)
EVE_results <- tibble::rownames_to_column(EVE_results, "Orthogroup")

EVE_results$type <- ifelse(EVE_results$beta<0.2839223,"Lineage Specific","Highly Variable")
EVE_results$significant <- ifelse(EVE_results$pval<=0.05,"Significant","Not Significant")
EVE_results$category <- ifelse(EVE_results$significant == "Significant",EVE_results$type,"NS")
colnames(EVE_results)
EVE_results <- EVE_results %>% relocate(type, .after = pval)
EVE_results <- EVE_results %>% relocate(significant, .after = type)
EVE_results <- EVE_results %>% relocate(category, .after = significant)

write.csv(EVE_results, file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/EVE_results_Host_Orthogroups_SingleBestOnly.csv",row.names = FALSE)

#Visualize LRT v Beta by volcano plot for gene data
Shared_EVE <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/EVE_results_Host_Orthogroups_SingleBestOnly.csv")
Shared_EVE_sig <- Shared_EVE %>% filter(pval <= 0.05)
write.csv(Shared_EVE_sig, file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/Signficant_EVE_results_Host_Orthogroups_SingleBestOnly.csv")
# Make a basic volcano plot
Shared_EVE_LS <- Shared_EVE_sig %>% filter(log(beta)< (-1.259055))
dim(Shared_EVE_LS)  #277 #76 Longest ORF # 332 ORF # 332 SingleBestORF
Shared_EVE_HV <- Shared_EVE_sig %>% filter(log(beta)> (-1.259055))
dim(Shared_EVE_HV)  #651 #186 Longest ORF # 588 ORF # # 588 SingleBestORF
# 928 total #262 total Longest ORF # 920 total ORF
332+588 #920 total
```

```{r EVE_volcano,fig.height=2,fig.width=2,dpi=300}
p <- ggplot(data = Shared_EVE,
       aes(x=log(beta),y=-log(pval),color=category))+
  geom_point()+
  scale_color_manual(values=c("#f37225","#00b9e3","black"))+
  geom_vline(xintercept = -1.259055, col="black", linetype="dashed")+
  annotate("text", x = 0.9, y = 9.5, label = "Shared beta", size = 3)+
  geom_hline(yintercept = -log(0.05), col="black",linetype = "dashed")+
  annotate("text", x = -4, y = 2.5, label = "p = 0.05", size = 3)+
  xlim(-5,5)+
  theme_light()
p + labs(color = "EVE Gene Category",title = "EVE Genes Volcano Plot",subtitle = "332 Lineage-Specific and 588 Highly Variable Genes") +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
```

![image](https://github.com/user-attachments/assets/767cf84a-c9d0-40ce-bcdb-c6cd9d73a98c)

## Highly Variable Genes

```{r,results='hide',tidy=TRUE}
sig_EVE <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/Signficant_EVE_results_Host_Orthogroups_SingleBestOnly.csv")
colnames(sig_EVE)
sig_EVE <- sig_EVE[,c(2,7,8,11)]   #Orthogroup,beta, pval, category
dim(sig_EVE)   ##262 Sig EVE genes ORF # 928 significant EVE genes 
names(sig_EVE)
rlogs <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/rlog_metadata_nodesign.csv")
rlogs <- as.data.frame(NodesignTbl.t)
rlogs$Orthogroup <- rownames(rlogs)

sig_EVE_rlogs <- left_join(sig_EVE,rlogs,by="Orthogroup")    #920
dim(sig_EVE_rlogs)
names(sig_EVE_rlogs)
write.csv(sig_EVE_rlogs,file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/Signficant_EVE_results_Host_Orthogroups_rlogs_SingleBestOnly.csv")

plasticGenes <- subset(sig_EVE_rlogs, category == "Highly Variable")
names(plasticGenes)
plasticGenes <- plasticGenes[,-c(2:4)]

dim(plasticGenes)
# [1] 651  77
# [1] 186  77
# [1] 716  76
names(plasticGenes)
t_plasticGenes <- setNames(data.frame(t(plasticGenes[,-1])), plasticGenes[,1])
row.names(t_plasticGenes)
t_plasticGenes <- tibble::rownames_to_column(t_plasticGenes, "SampleID")

label <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/MetaData.csv")
colnames(label)
colnames(t_plasticGenes)
#label <- label[-c(53),c(3,4,5,7,8)] # remove pstr_d8; keep name,treatment,status,LGR,dominant_clade
HV_metadata <- merge(label,t_plasticGenes,by="SampleID")
HV_metadata
dim(HV_metadata)
# [1]  70 668
# [1]  70 203 ORF
# [1]  70 473 ORF pval 0.1
# [1]  70 605
# [1]  69 726

write.csv(HV_metadata,file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/highlyVariable_rlogs_orthogroup_SingleBestOnly.csv")
```
# HV PCA
```{r}
# Load required libraries
library(ggplot2)
library(dplyr)

# Assume HV_metadata is already loaded and cleaned
# Select only numeric columns for PCA (remove metadata columns)
numeric_data <- HV_metadata[,11:length(HV_metadata)]

# Perform PCA
pca_results <- prcomp(numeric_data, scale. = TRUE)

# Calculate the percentage of variation explained by each principal component
percent_variance <- round(100 * (pca_results$sdev^2 / sum(pca_results$sdev^2)), 1)

# Create a data frame with the PCA results and Coral.Species
pca_data <- as.data.frame(pca_results$x)  # Extract the PCA coordinates
pca_data$Coral.Species <- HV_metadata$Coral.Species  # Add Coral.Species for coloring
pca_data$SampleID <- HV_metadata$SampleID  # Optionally add SampleID for labeling

# Plot the PCA results using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Coral.Species)) +
  geom_point(size = 3, alpha = 0.8) +  # Points representing samples
  stat_ellipse(level = 0.95, linetype = 2) +  # Add ellipses around species groups
  labs(
    title = "PCA of HV_metadata", 
    x = paste0("PC1: ", percent_variance[1], "% variance"), 
    y = paste0("PC2: ", percent_variance[2], "% variance")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_color_brewer(palette = "Set2")  # Use a nice color palette for species

```
> The purpose in producing a PCA of Highly variable genes is to observe that samples are not separating by their coral species. Overall Samples are overlapping by species, demonstrating the expression of these highly variable genes is not being driven by species as expected.

![image](https://github.com/user-attachments/assets/6cfd5d27-4de6-4d87-9a20-5519c2a74188)


## Lineage-Specific Genes

```{r,results='hide',tidy=TRUE}

LSGenes <- subset(sig_EVE_rlogs, sig_EVE$category == "Lineage Specific")
names(LSGenes)
LSGenes <- LSGenes[,-c(2:4)]
dim(LSGenes)
# [1] 277  77
names(LSGenes)
t_LSGenes <- setNames(data.frame(t(LSGenes[,-1])), LSGenes[,1])
row.names(t_LSGenes)
t_LSGenes <- tibble::rownames_to_column(t_LSGenes, "SampleID")

label <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/MetaData.csv")
colnames(label)
colnames(t_LSGenes)
#label <- label[-c(53),c(3,4,5,7,8)] # remove pstr_d8; keep name,treatment,status,LGR,dominant_clade
LS_metadata <- merge(label,t_LSGenes,by="SampleID")
LS_metadata
dim(LS_metadata)
# [1]  70 34
write.csv(LS_metadata,file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/LineageSpecific_rlogs_orthogroup_singlebestORF.csv")
```
### LS PCA
```{r}
# Load required libraries
library(ggplot2)
library(dplyr)

# Assume HV_metadata is already loaded and cleaned
# Select only numeric columns for PCA (remove metadata columns)
numeric_data <- LS_metadata[,11:length(LS_metadata)]

# Perform PCA
pca_results <- prcomp(numeric_data, scale. = TRUE)

# Calculate the percentage of variation explained by each principal component
percent_variance <- round(100 * (pca_results$sdev^2 / sum(pca_results$sdev^2)), 1)

# Create a data frame with the PCA results and Coral.Species
pca_data <- as.data.frame(pca_results$x)  # Extract the PCA coordinates
pca_data$Coral.Species <- LS_metadata$Coral.Species  # Add Coral.Species for coloring
pca_data$SampleID <- LS_metadata$SampleID  # Optionally add SampleID for labeling

# Plot the PCA results using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Coral.Species)) +
  geom_point(size = 3, alpha = 0.8) +  # Points representing samples
  stat_ellipse(level = 0.95, linetype = 2) +  # Add ellipses around species groups
  labs(
    title = "PCA of LS_metadata", 
    x = paste0("PC1: ", percent_variance[1], "% variance"), 
    y = paste0("PC2: ", percent_variance[2], "% variance")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_color_brewer(palette = "Set2")  # Use a nice color palette for species

```
> Because these genes were classified as Lineage-specific, we would expect them to not separate by their coral species. As we can see the samples do separate by their species.
![image](https://github.com/user-attachments/assets/b2778473-ec3b-4321-a643-c08ae06560cf)


## Annotating Orthogroups by Representative Annotations

### Read in Annotated Transcripts
```{r}
# In terminal: 
scp nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/*_host_orthologs_annotated_e-5.txt /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/R/Coral\ Host/

# update working directory to where the annotated orthologs are stored:
setwd("~/Desktop/Microbial Metatranscriptomics/R/Coral Host")

# Load required libraries
library(dplyr)
library(readr)

# List of input file prefixes
file_prefixes <- c("Acer", "Mcav", "Ofav", "Past")

# Function to process a single file
process_file <- function(prefix) {
  # Construct file paths
  input_file <- paste0(prefix, "_host_orthologs_annotated_e-5.txt")
  output_file <- paste0(prefix, "_host_orthologs_annotated_e-5_formatted.txt")
  
  # Check if the input file exists
  if (!file.exists(input_file)) {
    message(paste("File not found:", input_file))
    return(NULL)
  }
  
  # Load the data
  data <- read.table(input_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
  
  # Format the data
  formatted_data <- data %>%
    mutate(Entry = sub(".*\\|(.*)\\|.*", "\\1", V1)) %>%
    select(Entry, Transcript = V2, Evalue = V3)
  
  # Write the output file
  write.table(formatted_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # Return the formatted data for verification
  return(formatted_data)
}

# Process each file and store results in a list
formatted_results <- lapply(file_prefixes, process_file)

# Optionally, print the first few rows of each processed file
names(formatted_results) <- file_prefixes
lapply(formatted_results, function(x) {
  if (!is.null(x)) head(x) else NULL
})

Acer_annot_transcripts <- formatted_results$Acer
Acer_annot_transcripts$Transcript <- gsub("\\..*","",Acer_annot_transcripts$Transcript) # Remove ORF Identifier

Mcav_annot_transcripts <- formatted_results$Mcav
Mcav_annot_transcripts$Transcript <- gsub("\\..*","",Mcav_annot_transcripts$Transcript)

Ofav_annot_transcripts <- formatted_results$Ofav
Ofav_annot_transcripts$Transcript <- gsub("\\..*","",Ofav_annot_transcripts$Transcript)

Past_annot_transcripts <- formatted_results$Past
Past_annot_transcripts$Transcript <- gsub("\\..*","",Past_annot_transcripts$Transcript)

```

### Annotate Orthogroups
```{r}
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Acerhost_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Mcavhost_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Ofavhost_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Pasthost_orthogroup_tx2gene_annot.csv",row.names = FALSE)


```

### Annotating Orthogroups based on the most Common Entry

```{r}
# Load required libraries
library(dplyr)

# Define file paths for the annotated files
files <- c("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Acerhost_orthogroup_tx2gene_annot.csv",
           "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Mcavhost_orthogroup_tx2gene_annot.csv",
           "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Ofavhost_orthogroup_tx2gene_annot.csv",
           "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Pasthost_orthogroup_tx2gene_annot.csv")

# Read all files and combine them into a single data frame with file origin
all_orthogroups <- lapply(files, function(file) {
  data <- read.csv(file)
  data$SourceFile <- basename(file) # Add the file name as a new column
  return(data)
}) %>% bind_rows()

# Identify the most common Entry for each Orthogroup
representative_entries <- all_orthogroups %>%
  group_by(Orthogroup) %>%
  count(Entry, sort = TRUE) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>% # Get the most common Entry per Orthogroup
  ungroup() %>%
  select(Orthogroup, RepresentativeEntry = Entry)

# Merge representative entries back to the original data for reference
all_orthogroups_with_representative <- all_orthogroups %>%
  left_join(representative_entries, by = "Orthogroup")

# Save results
write.csv(all_orthogroups_with_representative,
          file = "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Orthogroups_with_Representative_Entry.csv",
          row.names = FALSE)

# Save representative entries only
write.csv(representative_entries,
          file = "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Representative_Entries.csv",
          row.names = FALSE)

# Print the first few rows of the representative entries
head(representative_entries)

```

### Apply Representative Orthogroups Entry Annotations to all_species

```{r}

all_species_annot <- merge(representative_entries,all_species, by ="Orthogroup")

uniprot <- read_csv("~/Desktop/Microbial Metatranscriptomics/References/uniprot.csv")

all_species_annot <- merge(uniprot,all_species_annot, by.x="Entry", by.y="RepresentativeEntry")
```

## EVE Gene Box and Whisker Plots
These are for loop chunks that will create box and whisker plots of all EVE genes and then save them into folders in your computer.

### HV Box and Whisker Plots
```{r}
library(dplyr)
library(ggplot2)
library(broom)

# Custom fill colors
custom_colors <- c("Control" = "cornflowerblue", "Healthy" = "darkseagreen", "Disease" = "salmon")

# Step 1: Extract Orthogroup expression data
orthogroup_data <- HV_metadata[, 11:ncol(HV_metadata)]

# Step 2: ANOVA + Tukey to find significant orthogroups
signif_results <- list()
anova_summary_list <- list()

for (og in colnames(orthogroup_data)) {
  temp_data <- HV_metadata %>%
    select(Coral.Species, Outcome, all_of(og)) %>%
    rename(Expression = all_of(og))
  
  anova_result <- aov(Expression ~ Outcome, data = temp_data)
  p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
  
  if (p_value < 0.1) {
    tukey_result <- TukeyHSD(anova_result)
    signif_results[[og]] <- list(
      orthogroup = og,
      p_value = p_value,
      tukey = tidy(tukey_result)
    )
    
    # Extract annotation
    og_annot <- all_species_annot[all_species_annot$Orthogroup == og, ]
    entry <- if ("Entry" %in% colnames(og_annot)) as.character(og_annot$Entry[1]) else NA
    protein_name <- if ("Protein names" %in% colnames(og_annot)) as.character(og_annot$`Protein names`[1]) else NA

    # Add to summary list
    anova_summary_list[[og]] <- data.frame(
      Orthogroup = og,
      P_value = p_value,
      Entry = entry,
      Protein_Name = protein_name,
      stringsAsFactors = FALSE
    )
  }
}

# Step 3: Flatten results and sort
sig_orthogroups <- do.call(rbind, lapply(signif_results, function(x) data.frame(Orthogroup = x$orthogroup, x$tukey)))
anova_summary_df <- do.call(rbind, anova_summary_list)
anova_summary_df <- anova_summary_df[order(anova_summary_df$P_value), ]

# Step 4: Create output directory if it doesn't exist
output_dir <- "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/HV_BoxandWhisker"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save ANOVA summary
write.csv(anova_summary_df, file = file.path(output_dir, "anova_summary_HV.csv"), row.names = FALSE)

# Step 5: Plot with annotation
for (og in unique(sig_orthogroups$Orthogroup)) {
  
  # Get annotation again for plotting
  og_annot <- all_species_annot[all_species_annot$Orthogroup == og, ]
  entry <- if ("Entry" %in% colnames(og_annot)) as.character(og_annot$Entry[1]) else "Entry unavailable"
  protein_name <- if ("Protein names" %in% colnames(og_annot)) as.character(og_annot$`Protein names`[1]) else "Protein name unavailable"
  
  subtitle_text <- paste(entry, "-", protein_name)

  # Filter metadata
  temp_data <- HV_metadata %>%
    select(Coral.Species, Outcome, all_of(og)) %>%
    rename(Expression = all_of(og))

  # Recalculate ANOVA p-value for caption
  anova_p <- summary(aov(Expression ~ Outcome, data = temp_data))[[1]]$`Pr(>F)`[1]
  caption_text <- paste("ANOVA p-value:", signif(anova_p, 3))

  # Set ordering
  temp_data$Coral.Species <- factor(temp_data$Coral.Species, levels = c("Mcav", "Past", "Ofav", "Acer"))
  temp_data$Outcome <- factor(temp_data$Outcome, levels = c("Control", "Healthy", "Disease"))

  # Species plot
  species_plot <- ggplot(temp_data, aes(x = Coral.Species, y = Expression, fill = Outcome)) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_hline(yintercept = mean(temp_data$Expression, na.rm = TRUE),
               linetype = "dashed", color = "red", size = 0.8) +
    labs(
      title = paste("Expression of", og, "by Coral Species and Outcome"),
      subtitle = subtitle_text,
      caption = caption_text,
      x = "Coral Species",
      y = "rlog Expression",
      fill = "Outcome"
    ) +
    theme_minimal() +
    scale_fill_manual(values = custom_colors)

  print(species_plot)
  ggsave(filename = file.path(output_dir, paste0(og, "_species.pdf")),
         plot = species_plot, device = "pdf", width = 8, height = 6)

  # Outcome plot (boxplot by outcome only)
  outcome_plot <- ggplot(temp_data, aes(x = Outcome, y = Expression, fill = Outcome)) +
    geom_boxplot(width = 0.6) +
    geom_hline(yintercept = mean(temp_data$Expression, na.rm = TRUE),
               linetype = "dashed", color = "red", size = 0.8) +
    labs(
      title = paste("Expression of", og, "by Outcome"),
      subtitle = subtitle_text,
      caption = caption_text,
      x = "Outcome",
      y = "rlog Expression",
      fill = "Outcome"
    ) +
    theme_minimal() +
    scale_fill_manual(values = custom_colors)

  print(outcome_plot)
  ggsave(filename = file.path(output_dir, paste0(og, "_outcome.pdf")),
         plot = outcome_plot, device = "pdf", width = 8, height = 6)
}

```

### LS Box and Whisker Plots

```{r}
library(dplyr)
library(ggplot2)
library(broom)

# Custom fill colors
custom_colors <- c("Control" = "cornflowerblue", "Healthy" = "darkseagreen", "Disease" = "salmon")

# Step 1: Extract Orthogroup expression data
orthogroup_data <- LS_metadata[, 11:ncol(LS_metadata)]

# Step 2: ANOVA + Tukey to find significant orthogroups
signif_results <- list()

for (og in colnames(orthogroup_data)) {
  temp_data <- LS_metadata %>%
    select(Coral.Species, Outcome, all_of(og)) %>%
    rename(Expression = all_of(og))

  anova_result <- aov(Expression ~ Coral.Species, data = temp_data)
  p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]

  if (p_value < 0.1) {
    tukey_result <- TukeyHSD(anova_result)
    signif_results[[og]] <- list(
      orthogroup = og,
      p_value = p_value,
      tukey = tidy(tukey_result)
    )
  }
}

# Step 3: Flatten result and get only significant orthogroups
sig_orthogroups <- do.call(rbind, lapply(signif_results, function(x) data.frame(Orthogroup = x$orthogroup, x$tukey)))

# Step 4: Plot with annotation
for (og in unique(sig_orthogroups$Orthogroup)) {

  # Extract annotation
  og_annot <- all_species_annot[all_species_annot$Orthogroup == og, ]
  entry <- if ("Entry" %in% colnames(og_annot)) as.character(og_annot$Entry[1]) else "Entry unavailable"
  protein_name <- if ("Protein names" %in% colnames(og_annot)) as.character(og_annot$`Protein names`[1]) else "Protein name unavailable"

  subtitle_text <- paste(entry, "-", protein_name)

  # Filter metadata for this orthogroup
  temp_data <- LS_metadata %>%
    select(Coral.Species, Outcome, all_of(og)) %>%
    rename(Expression = all_of(og))

  # Recalculate ANOVA p-value for caption
  anova_p <- summary(aov(Expression ~ Coral.Species, data = temp_data))[[1]]$`Pr(>F)`[1]
  caption_text <- paste("ANOVA p-value:", signif(anova_p, 3))

  # Set custom factor levels for ordering
  temp_data$Coral.Species <- factor(temp_data$Coral.Species, levels = c("Mcav", "Past", "Ofav", "Acer"))
  temp_data$Coral.Species <- factor(temp_data$Coral.Species, levels = c("Mcav", "Past", "Ofav","Acer"))

  # Plot 1: by Coral Species (ordered + grouped by Outcome)
  species_plot <- ggplot(temp_data, aes(x = Coral.Species, y = Expression, fill = Outcome)) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_hline(yintercept = mean(temp_data$Expression, na.rm = TRUE),
               linetype = "dashed", color = "red", size = 0.8) +
    labs(
      title = paste("Expression of", og, "by Coral Species and Outcome"),
      subtitle = subtitle_text,
      caption = caption_text,
      x = "Coral Species",
      y = "rlog Expression",
      fill = "Outcome"
    ) +
    theme_minimal() +
    scale_fill_manual(values = custom_colors)

  print(species_plot)

  # Plot 2: Outcome boxplot (not grouped by Coral.Species)
  temp_data$Coral.Species <- factor(temp_data$Coral.Species, levels = c("Mcav", "Past", "Ofav","Acer"))

  outcome_plot <- ggplot(temp_data, aes(x = Coral.Species, y = Expression, fill = Outcome)) +
    geom_boxplot(width = 0.6) +
    geom_hline(yintercept = mean(temp_data$Expression, na.rm = TRUE),
               linetype = "dashed", color = "red", size = 0.8) +
    labs(
      title = paste("Expression of", og, "by Coral.Species"),
      subtitle = subtitle_text,
      caption = caption_text,
      x = "Coral.Species",
      y = "rlog Expression",
      fill = "Coral.Species"
    ) +
    theme_minimal() +
    scale_fill_manual(values = custom_colors)

  print(outcome_plot)
}


```
## HV DESEQ2
### ~Coral.Species+Treatment
```{r,results='hide',tidy=TRUE}

HV_metadata

# Extract count data
HV.countData <- HV_metadata[,c(1,11:length(HV_metadata))]

# Transpose the matrix
HV.countData.t <- t(HV.countData)

# Make the first row the column names
colnames(HV.countData.t) <- HV.countData.t[1, ]

# Remove the first row (now column names)
HV.countData.t <- as.data.frame(HV.countData.t[-1, ], stringsAsFactors = FALSE)

# Extract only the samples (rows) in HV.countData.t from the raw integer count matrix
matched_counts <- countData[rownames(HV.countData.t), ]

# (Optional) Check if row names matched correctly
all(rownames(matched_counts) == rownames(HV.countData.t))  # should return TRUE

# Extract Metadata
HV.colData <- HV_metadata[,1:10]

# Samples in matched_counts not present in HV.colData
setdiff(colnames(matched_counts), HV.colData$SampleID)

# Keep only columns (samples) in matched_counts that are in HV.colData$SampleID
matched_counts_filtered <- matched_counts[, colnames(matched_counts) %in% HV.colData$SampleID]


# Set the design
# **Pull orthogroups names from raw data as input for HV.countData instead
dds_treatment.HV <- DESeqDataSetFromMatrix(
  countData = matched_counts_filtered,
  colData = HV.colData,
  design = ~Coral.Species + Treatment
)


# Filter out columns with an average number of reads less than 10
dds_treatment.HV <- dds_treatment.HV[ rowMeans(counts(dds_treatment.HV)) > 10, ] 

# Run DESeq2
test.HV <- DESeq(dds_treatment.HV)
rld.HV <- rlog(dds_treatment.HV, blind=FALSE) 
rrld.HV <- assay(rld.HV)
# 9376 Orthogroups after filtering out low counts
write.csv(rrld.HV, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/dds_treatment.HV_rlog_orthogroups_ORF.csv")  
#write.csv(rrld, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/Species_rlog_orthogroups_ORF.csv")  

# Extract the results for the comparison of interest
# For example, comparing "TreatmentA" to "TreatmentB":
results_treatment.HV <- results(test.HV, contrast = c("Treatment", "Control", "Disease"))

# Order results by padj
results_treatment.HV <- results_treatment.HV[order(results_treatment.HV$padj), ]

# Exclude rows with NA in padj before filtering for significance
significant_results.HV <- results_treatment.HV[!is.na(results_treatment.HV$padj) & results_treatment.HV$padj < 0.1, ]

# View summary of significant results
summary(results_treatment.HV)

# Save the significant results to a CSV file
 write.csv(as.data.frame(significant_results.HV), file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/significant_dds_treatment.HV_comparison.csv")


```
out of 713 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 25, 3.5%
LFC < 0 (down)     : 84, 12%
outliers [1]       : 0, 0%
low counts [2]     : 83, 12%
(mean count < 55)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


## HV Heatmap 
```{r heatmap-preview, fig.height=16, fig.width=8} 

library(grid) 
p <- pheatmap( rld_grouped_t, scale = "row", cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete", show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 12, fontsize_row = 6, cellwidth = 20, cellheight = 8, main = "Significant HV Orthogroups (Preview)", silent = TRUE ) 
grid.newpage() 
grid.draw(p$gtable)
library(pheatmap)
library(dplyr)

# Step 1: Extract rlog matrix and subset to significant HV orthogroups
rld_mat <- assay(rld.HV)
sig_genes <- rownames(significant_results.HV)
rld_sig <- rld_mat[sig_genes, ]

# Step 2: Transpose to make samples rows and orthogroups columns
rld_sig_t <- as.data.frame(t(rld_sig))
rld_sig_t$SampleID <- rownames(rld_sig_t)

# Step 3: Merge with metadata to get Outcome
rld_sig_annotated <- merge(HV_metadata[, c("SampleID", "Outcome")], rld_sig_t, by = "SampleID")

# Step 4: Aggregate by Outcome (mean expression per orthogroup)
rld_grouped <- aggregate(. ~ Outcome, data = rld_sig_annotated[,-1], FUN = mean)

# Step 5: Set Outcome as rownames and transpose so that Outcome becomes columns
rownames(rld_grouped) <- rld_grouped$Outcome
rld_grouped <- rld_grouped[ , -1]  # remove Outcome column
rld_grouped_t <- t(as.matrix(rld_grouped))

# Step 6: Reorder columns manually (Control, Healthy, Disease)
desired_order <- c("Control", "Healthy", "Disease")
rld_grouped_t <- rld_grouped_t[, desired_order]

# Step 7: Replace rownames (Orthogroups) with truncated Protein names + OG
orthogroups <- rownames(rld_grouped_t)
protein_names <- sapply(orthogroups, function(og) {
  match_row <- all_species_annot[all_species_annot$Orthogroup == og, ]
  if ("Protein names" %in% colnames(match_row) && nrow(match_row) > 0) {
    name <- substr(as.character(match_row$`Protein names`[1]), 1, 50)
    paste0(name, " - ", og)
  } else {
    og
  }
})
rownames(rld_grouped_t) <- protein_names

# Step 8: Export expression metrics to CSV
heatmap_df <- as.data.frame(rld_grouped_t)
heatmap_df$Orthogroup <- orthogroups
heatmap_df$Protein_name <- sapply(orthogroups, function(og) {
  row <- all_species_annot[all_species_annot$Orthogroup == og, ]
  if (nrow(row) > 0 && "Protein names" %in% colnames(row)) {
    substr(as.character(row$`Protein names`[1]), 1, 50)
  } else {
    NA
  }
})
heatmap_df$Entry <- sapply(orthogroups, function(og) {
  row <- all_species_annot[all_species_annot$Orthogroup == og, ]
  if (nrow(row) > 0 && "Entry" %in% colnames(row)) {
    as.character(row$Entry[1])
  } else {
    NA
  }
})

# Reorder columns for export
heatmap_export <- heatmap_df[, c("Orthogroup", "Protein_name", "Entry", desired_order)]
write.csv(heatmap_export, "HV heatmap DESEq sig.csv", row.names = FALSE)

# Step 9: Save heatmap to a PDF file
pdf("HV_heatmap_DESeq_sig.pdf", width = 10, height = 12)
p <- pheatmap(
  rld_grouped_t,
  scale = "row",
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 12,
  fontsize_row = 6,
  cellwidth = 25,
  cellheight = 10,
  main = "Significant HV Orthogroups"
)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()


# Step 10: Render again inline with slightly reduced height for Rmd preview
pheatmap(
  rld_grouped_t,
  scale = "row",
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 12,
  fontsize_row = 6,
  cellwidth = 25,
  cellheight = 8,
  main = "Significant HV Orthogroups"
)


```
![image](https://github.com/user-attachments/assets/75494925-8480-4225-943c-9ec310e2d3d9)


## LS DESEQ2
### ~Treatment+Coral.Species
```{r,results='hide',tidy=TRUE}

LS_metadata

# Extract count data
LS.countData <- LS_metadata[,c(1,11:length(LS_metadata))]

# Transpose the matrix
LS.countData.t <- t(LS.countData)

# Make the first row the column names
colnames(LS.countData.t) <- LS.countData.t[1, ]

# Remove the first row (now column names)
LS.countData.t <- as.data.frame(LS.countData.t[-1, ], stringsAsFactors = FALSE)

# Extract only the samples (rows) in LS.countData.t from the raw integer count matrix
matched_counts <- countData[rownames(LS.countData.t), ]

# (Optional) Check if row names matched correctly
all(rownames(matched_counts) == rownames(LS.countData.t))  # should return TRUE

# Extract Metadata
LS.colData <- LS_metadata[,1:10]

# Samples in matched_counts not present in LS.colData
setdiff(colnames(matched_counts), LS.colData$SampleID)

# Keep only columns (samples) in matched_counts that are in LS.colData$SampleID
matched_counts_filtered <- matched_counts[, colnames(matched_counts) %in% LS.colData$SampleID]


# Set the design
# **Pull orthogroups names from raw data as input for LS.countData instead
dds_Coral.Species.LS <- DESeqDataSetFromMatrix(
  countData = matched_counts_filtered,
  colData = LS.colData,
  design = ~Coral.Species + Coral.Species
)


# Filter out columns with an average number of reads less than 10
dds_Coral.Species.LS <- dds_Coral.Species.LS[ rowMeans(counts(dds_Coral.Species.LS)) > 10, ] 

# Run DESeq2
test.LS <- DESeq(dds_Coral.Species.LS)
rld.LS <- rlog(dds_Coral.Species.LS, blind=FALSE) 
rrld.LS <- assay(rld.LS)
# 9376 Orthogroups after filtering out low counts
write.csv(rrld.LS, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/dds_Coral.Species.LS_rlog_orthogroups_ORF.csv")  
#write.csv(rrld, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/Species_rlog_orthogroups_ORF.csv")  

# Extract the results for the comparison of interest
# Perform pairwise comparisons for Coral.Species
res_LS_Acer_vs_Mcav <- results(test.LS, contrast = c("Coral.Species", "Acer", "Mcav"))
res_LS_Acer_vs_Ofav <- results(test.LS, contrast = c("Coral.Species", "Acer", "Ofav"))
res_LS_Acer_vs_Past <- results(test.LS, contrast = c("Coral.Species", "Acer", "Past"))
res_LS_Mcav_vs_Ofav <- results(test.LS, contrast = c("Coral.Species", "Mcav", "Ofav"))
res_LS_Mcav_vs_Past <- results(test.LS, contrast = c("Coral.Species", "Mcav", "Past"))
res_LS_Ofav_vs_Past <- results(test.LS, contrast = c("Coral.Species", "Ofav", "Past"))

summary(res_LS_Acer_vs_Mcav)
summary(res_LS_Acer_vs_Ofav)
summary(res_LS_Acer_vs_Past)
summary(res_LS_Mcav_vs_Ofav)
summary(res_LS_Mcav_vs_Past)
summary(res_LS_Ofav_vs_Past)

# Order results by padj
res_LS_Acer_vs_Mcav <- res_LS_Acer_vs_Mcav[order(res_LS_Acer_vs_Mcav$padj), ]
res_LS_Acer_vs_Ofav <- res_LS_Acer_vs_Ofav[order(res_LS_Acer_vs_Ofav$padj), ]
res_LS_Acer_vs_Past <- res_LS_Acer_vs_Past[order(res_LS_Acer_vs_Past$padj), ]
res_LS_Mcav_vs_Ofav <- res_LS_Mcav_vs_Ofav[order(res_LS_Mcav_vs_Ofav$padj), ]
res_LS_Mcav_vs_Past <- res_LS_Mcav_vs_Past[order(res_LS_Mcav_vs_Past$padj), ]
res_LS_Ofav_vs_Past <- res_LS_Ofav_vs_Past[order(res_LS_Ofav_vs_Past$padj), ]

# Exclude rows with NA in padj before filtering for significance
significant_res_LS_Acer_vs_Mcav <- res_LS_Acer_vs_Mcav[!is.na(res_LS_Acer_vs_Mcav$padj) & res_LS_Acer_vs_Mcav$padj < 0.1, ]
significant_res_LS_Acer_vs_Ofav <- res_LS_Acer_vs_Ofav[!is.na(res_LS_Acer_vs_Ofav$padj) & res_LS_Acer_vs_Ofav$padj < 0.1, ]
significant_res_LS_Acer_vs_Past <- res_LS_Acer_vs_Past[!is.na(res_LS_Acer_vs_Past$padj) & res_LS_Acer_vs_Past$padj < 0.1, ]
significant_res_LS_Mcav_vs_Ofav <- res_LS_Mcav_vs_Ofav[!is.na(res_LS_Mcav_vs_Ofav$padj) & res_LS_Mcav_vs_Ofav$padj < 0.1, ]
significant_res_LS_Mcav_vs_Past <- res_LS_Mcav_vs_Past[!is.na(res_LS_Mcav_vs_Past$padj) & res_LS_Mcav_vs_Past$padj < 0.1, ]
significant_res_LS_Ofav_vs_Past <- res_LS_Ofav_vs_Past[!is.na(res_LS_Ofav_vs_Past$padj) & res_LS_Ofav_vs_Past$padj < 0.1, ]


```


### LS Heatmap
```{r, fig.height=16, fig.width=8}
library(pheatmap)

# Step 1: Extract rlog matrix and subset to significant LS orthogroups
# Combine all significant orthogroups from all pairwise species comparisons
sig_orthos <- unique(c(
  rownames(significant_res_LS_Acer_vs_Mcav),
  rownames(significant_res_LS_Acer_vs_Ofav),
  rownames(significant_res_LS_Acer_vs_Past),
  rownames(significant_res_LS_Mcav_vs_Ofav),
  rownames(significant_res_LS_Mcav_vs_Past),
  rownames(significant_res_LS_Ofav_vs_Past)
))

# Subset the rlog-transformed matrix
rld_mat <- assay(rld.LS)
rld_sig <- rld_mat[sig_orthos, ]

# Step 2: Transpose to make samples rows and orthogroups columns
rld_sig_t <- as.data.frame(t(rld_sig))
rld_sig_t$SampleID <- rownames(rld_sig_t)

# Step 3: Merge with metadata to get Coral.Species
rld_sig_annotated <- merge(LS_metadata[, c("SampleID", "Coral.Species")], rld_sig_t, by = "SampleID")

# Step 4: Aggregate by Coral.Species (mean expression per orthogroup)
rld_grouped <- aggregate(. ~ Coral.Species, data = rld_sig_annotated[,-1], FUN = mean)
rownames(rld_grouped) <- rld_grouped$Coral.Species
rld_grouped <- rld_grouped[, -1]  # remove "Coral.Species" column
rld_grouped_t <- t(as.matrix(rld_grouped))  # genes = rows, Coral.Speciess = columns

# Step 5: Reorder columns
desired_order <- c("Mcav", "Past", "Ofav","Acer")
rld_grouped_t <- rld_grouped_t[, desired_order]

# Step 6: Replace rownames with Protein Name - Orthogroup (using sig_orthos for matching)
orthogroups <- rownames(rld_grouped_t)

# Ensure all_species_annot is in character format
all_species_annot$Orthogroup <- as.character(all_species_annot$Orthogroup)

# Make a lookup table: Orthogroup → Protein Name
protein_lookup <- sapply(orthogroups, function(og) {
  match_row <- all_species_annot[all_species_annot$Orthogroup == og, ]
  if (nrow(match_row) > 0 && "Protein names" %in% colnames(match_row)) {
    protein_name <- as.character(match_row$`Protein names`[1])
    protein_name <- ifelse(is.na(protein_name) || protein_name == "", "Unknown protein", protein_name)
    paste0(substr(protein_name, 1, 50), " - ", og)
  } else {
    paste0("Unknown protein - ", og)
  }
})

# Update rownames
rownames(rld_grouped_t) <- protein_lookup


# Step 7: Plot heatmap
pheatmap(
  rld_grouped_t,
  scale = "row",
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 10,
  fontsize_row = 6,
  cellwidth = 15,
  cellheight = 6,
  main = "DESeq Significant LS Orthogroups"
)


# ----- Export the heatmap to PDF -----
# Export heatmap
pdf("LS_ProteinName_Orthogroup_Heatmap.pdf", width = 10, height = 12)
pheatmap(
  rld_grouped_t,
  scale = "row",
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 10,
  fontsize_row = 6,
  cellwidth = 15,
  main = "Significant LS Orthogroups\nProtein Name - Orthogroup"
)
dev.off()

# ----- Create and Export Annotation + Expression Data Table -----
# Export annotated CSV
# Step 1: Convert to dataframe and preserve custom rownames
rld_grouped_df <- as.data.frame(rld_grouped_t)
rld_grouped_df$Row_Label <- rownames(rld_grouped_t)

# Step 2: Extract Orthogroup from rownames (after the last " - ")
rld_grouped_df$Orthogroup <- sub(".* - (OG[0-9]+)$", "\\1", rld_grouped_df$Row_Label)

# Step 3: Join with Entry info from all_species_annot
entry_lookup <- all_species_annot[, c("Orthogroup", "Entry")]
entry_lookup$Orthogroup <- as.character(entry_lookup$Orthogroup)

# Step 4: Merge in Entry info
rld_grouped_with_entry <- merge(rld_grouped_df, entry_lookup, by = "Orthogroup", all.x = TRUE)

# Step 5 (optional): Reorder columns
rld_grouped_with_entry <- rld_grouped_with_entry[, c("Row_Label", "Orthogroup", "Entry", setdiff(colnames(rld_grouped_with_entry), c("Row_Label", "Orthogroup", "Entry")))]

# Preview result
head(rld_grouped_with_entry)

write.csv(rld_grouped_with_entry, "LS_ProteinName_Orthogroup_Heatmap_DataResults.csv", row.names = FALSE)

```
[LS_ProteinName_Orthogroup_Heatmap.pdf](https://github.com/user-attachments/files/20524719/LS_ProteinName_Orthogroup_Heatmap.pdf)



# WGCNA
```{r, fig.width=14, fig.height=6}
#install.packages('BiocManager')
#library(BiocManager)
#BiocManager::install('WGCNA')
#BiocManager::install('flashClust')

library(WGCNA)
library(flashClust)
library(curl)

head(HV_metadata)

# Retaining only Expression data
expression.data <- HV_metadata[,11:length(HV_metadata)] #removing variables not holding expression data
row.names(expression.data) <- HV_metadata[,1]
names(expression.data)


# Identifying Good Genes
gsg <-goodSamplesGenes(expression.data)
summary(gsg)

#If the allOK object returns true, which it does in this case, there are no outliers present. If it doesn’t return true, we need to filter out the outliers manually using the following code --> See "https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html" for WGCNA tutorial.
gsg$allOK


# Another way to identify outliers is to visualize by hierarchical tree
sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

# You can remove the outlier using a cutree function.
#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 50, col = "red"); # 50 is an arbitrart cutoff, however it is visually guided.

# 43 and 57 appear to be outliers and can be removed by setting the cut height to 50 here:
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 50, minSize = 10) #returns numeric vector
#Remove outlier
expression.data <- expression.data[cut.sampleTree==1, ]


#Pick the Soft Threshold Power
spt <- pickSoftThreshold(expression.data) 
spt

# Plot the r2 values as a function of the soft thresholds
#We should be maximizing the r2 value and minimizing mean connectivity.
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

# Plot the Mean Connectivity (This was not working in rmakrdonw I needed to copy and paiste into the console to work. Odd)
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

#You can determine the soft power threshold should be set to 5 as it is the spt that retains the highest mean connectivity while reaching an r2 value above 0.80.
#NOTE: the higher the value, the stronger the connection strength will be of highly correlated gene expression profiles and the more devalued low correlations will be.

# Calling the Adjacency Function
#Now that you have the soft threshold power determined you can call on the adjacency() function of the WGCNA package.

#REMINDER: This function calculates the similarity measurement and transforms the similarity by the adjacency function and generates a weighted network adjacency matrix.

softPower <- 5
adjacency <- adjacency(expression.data, power = softPower)

# Module Construction
TOM <- TOMsimilarity(adjacency)

# To convert this matrix into a dissimilarity matrix you can subtract the TOM object from 1.
TOM.dissimilarity <- 1-TOM

#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
labels = FALSE, hang = 0.04)

# Identify Modules - Set ClusterSize 30 is default.
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)

table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module. 

#Plot the module assignment under the gene dendrogram for visualization.
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")


# Module Eigengene Identification
#A ME (Module Eigengene) is the standardized gene expression profile for a given module.(Summary statistic of Expression)
#To identify the Module Eigengene you can call on the expression data into the moduleEigengenes() function.
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

# Module Merging (optional)

#To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations and merge the modules that have similar expression profiles.

#REMINDER: An eigengene is the gene whose expression is representative of the the majority of genes expressed within a module.

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

#Plot
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75


merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)

# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# The similar modules are now merged! Let’s compare them with the original modules by creating another dendrogram
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
c("Original Module", "Merged Module"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors for original and merged modules")


# External Trait Matching
#Once you have constructed the network (the adjacency matrix) and divided it into meaningful modules, you can begin to relate the network to external traits.

#Start by reading in Trait Data Table
#allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA.csv") # I took the file MetaData.csv and converted some of the categorical phenotype data into continuous data and added relative risk and disease prevalence. We can also add bacteria abundances, bacteria gene expression, and algal gene expression in time.

# Combine module eigengenes with traits
#allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA_ME.csv")
#allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA_ME_Bac16s.csv")
allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA_ME_Bac16s_Net_5.19.csv")


# Match Sample IDs and retain only relevant rows
Samples <- rownames(expression.data)
traitRows <- match(Samples, allTraits$SampleID)
datTraits <- allTraits[traitRows,]
rownames(datTraits) <- allTraits[traitRows, 1]

# Ensure mergedMEs has rownames as SampleIDs (should match expression.data)
mergedMEs$SampleID <- rownames(mergedMEs)

# Clean join: no need to add SampleID column if it's already there
merged_data <- mergedMEs %>%
  left_join(datTraits, by = "SampleID")

# Remove NAs caused by missing eigengenes
merged_data_clean <- na.omit(merged_data)

# Separate back into matrices
mergedMEs_clean <- merged_data_clean %>% select(starts_with("ME")) %>% as.data.frame()
datTraits_clean <- merged_data_clean %>% select(-SampleID, -starts_with("ME")) %>% as.data.frame()

# Run correlations
nSamples <- nrow(mergedMEs_clean)
module.trait.correlation <- cor(mergedMEs_clean, datTraits_clean, use = "p")
module.trait.Pvalue <- corPvalueStudent(module.trait.correlation, nSamples)

# Prepare text matrix for heatmap
textMatrix <- paste(signif(module.trait.correlation, 2), "\n(",
                    signif(module.trait.Pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module.trait.correlation)

# Plot heatmap

labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits_clean),
               yLabels = names(mergedMEs_clean),
               ySymbols = names(mergedMEs_clean),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Coral Host Highly Variable Module-trait relationships"))

```

![image](https://github.com/user-attachments/assets/4751d05c-6690-4ecc-a563-fe219c9243ae)


# Bact_ME_CellvibrionalesFiltered_brown - red Module
```{r}

# Target Gene Identification
'''
You can use the gene significance along with the genes intramodular connectivity to identify potential target genes associated with a particular trait of interest. For this analysis weight will be the clinical trait.

Connectivity - how connected a speficic node is in the network (how many nodes have high correlation with that node). High connectivity indicates a hub gene (central to many nodes). Whole Network connectivity - a measure for how well the node is connected throughout the entire system Intramodular connectivity - a measure for how well the node is connected within its assigned module. Also an indicator for how well that node belongs to its module. This is also known as module membership.

The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. This quantifies the similarity of all genes on the array to every module.
'''
# Define variable Bact_ME_CellvibrionalesFiltered_brown containing the Bact_ME_CellvibrionalesFiltered_brown column of datTrait
Bact_ME_CellvibrionalesFiltered_brown = as.data.frame(datTraits$Bact_ME_CellvibrionalesFiltered_brown)
names(Bact_ME_CellvibrionalesFiltered_brown) = "Bact_ME_CellvibrionalesFiltered_brown"

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, Bact_ME_CellvibrionalesFiltered_brown, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Bact_ME_CellvibrionalesFiltered_brown), sep="")
names(GSPvalue) = paste("p.GS.", names(Bact_ME_CellvibrionalesFiltered_brown), sep="")
head(GSPvalue)

# Using the gene significance you can identify genes that have a high significance for Bact_ME_CellvibrionalesFiltered_brown Using the module membership measures you can identify genes with high module membership in interesting modules.

# As an example, you can look at the red module as it has the highest significant association with weight (.59).
#Plot a scatter plot of gene significance vs. module membership in the red module.

par(mar=c(1,1,1,1))
module = "red"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
abs(geneTraitSignificance[moduleGenes,1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for SCTLD Bact_ME_CellvibrionalesFiltered_brown",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

'''
The red gene significance and module membership have a positive correlation of .73 with a very significant p-value. This indicates that the genes that are highly significantly associated with the trait (high gene significance) are also the genes that are the most connected within their module (high module membership). Therefore genes in the red module could be potential target genes when looking at Disease Bact_ME_CellvibrionalesFiltered_brown.
'''
```
### Get a list of genes correlated to our trait of interest
```{r}
#Summary output of network analysis results

names(expression.data)
names(expression.data)[mergedColors=="red"]

all_species_annot
dim(all_species_annot)
names(all_species_annot)
probes = names(expression.data)[mergedColors=="red"]
probes2annot = match(probes, all_species_annot$Orthogroup)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0. It Doesnt. lol. Good luck.

# Create the starting data frame
geneInfo0 = data.frame(
  Orthogroup = probes,
  geneSymbol = all_species_annot$`Protein names`[probes2annot],
  LocusLinkID = all_species_annot$Entry[probes2annot],
  moduleColor = (ModuleColors[moduleGenes]),
  GS.Bact_ME_CellvibrionalesFiltered_brown = geneTraitSignificance[moduleGenes, 1],
  p.GS.Bact_ME_CellvibrionalesFiltered_brown = GSPvalue[moduleGenes, 1]
)

# Order modules by their significance with Bact_ME_CellvibrionalesFiltered_brown
modOrder = order(-abs(cor(MEs, Bact_ME_CellvibrionalesFiltered_brown, use = "p")))

# Loop through modules and append membership and p-values
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  
  # Append module membership and its p-value
  geneInfo0 = data.frame(
    geneInfo0,
    MM = geneModuleMembership[moduleGenes, modOrder[mod]],
    pMM = MMPvalue[moduleGenes, modOrder[mod]]
  )
  
  # Rename new columns with module names
  names(geneInfo0) = c(
    oldNames,
    paste("MM.", modNames[modOrder[mod]], sep = ""),
    paste("p.MM.", modNames[modOrder[mod]], sep = "")
  )
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Bact_ME_CellvibrionalesFiltered_brown));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo_redModule_HV_Bact_ME_CellvibrionalesFiltered_brown.csv")

geneInfo$geneSymbol
```

### Functional Redundancy

> These next three code chunks will extract the genes that are signficantly correlated to traits of interest (our bacteria from 16s ancombc2 that sig increased in disease treated coral). We will identify the genes that are the most frequently correlated to these bacteria, identify the bacteria they are correlated to and then finally include the GS values which tells us if that correlation is positive or negative between a specific gene and specific bacteria.
>
> ```{r, fig.height=16, fig.width=12}
# Load libraries
library(WGCNA)
library(dplyr)
library(ggplot2)

# ---- INPUTS ----
traitList <- c("JTB23_f87fd1fab019a118b3a6505944739d22_JTB23",
               "Verrucomicrobiales_99810b14f7cf45a9a4b2ddca790e5519_Verrucomicrobiales",
               "Rhodobacterales_48a6a0e363dfb3823afdb092a0bb834d_Rhodobacterales",
               "Cellvibrionales_bb2b2e03048bf07414f3dfb6a85ca835_Cellvibrionales",
               "Vibrionales_d0c7ebb307c0dae2345e684d9904b117_Vibrionales",
               "Peptostreptococcales.Tissierellales_1219316521169aa7d8fa7e220ffeb483_Peptostreptococcales.Tissierellales",
               "Desulfovibrionales_e9590b36781b69df7b2f7958c8003724_Desulfovibrionales",
               "Campylobacterales_07420998b15f4c1a118e96daa30c2f18_Campylobacterales",
               "Bacteroidales_f5f370a81cdaa64e8fd99723045ebe33_Bacteroidales",
               "Francisellales_f6878fcb75e7a023b6012c9b023c6e35_Francisellales"
)

# ---- CALCULATE MODULE MEMBERSHIP ----
modNames <- substring(names(mergedMEs), 3) # Remove 'ME' prefix

geneModuleMembership <- as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# ---- CALCULATE GENE SIGNIFICANCE TO BACTERIA ----
bacteria_matrix <- datTraits[, traitList]
geneTraitSignificance <- as.data.frame(cor(expression.data, bacteria_matrix, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitList, sep = "")
names(GSPvalue) <- paste("p.GS.", traitList, sep = "")

# ---- BUILD GENE INFO TABLE BY MODULE ----
all_geneInfo <- list()

# Define full column structure to standardize
trait_cols <- unlist(lapply(traitList, function(trait) c(paste0("GS.", trait), paste0("p.GS.", trait))))
mm_cols <- unlist(lapply(modNames, function(mod) c(paste0("MM.", mod), paste0("p.MM.", mod))))
full_colnames <- c("Orthogroup", "geneSymbol", "LocusLinkID", "moduleColor", trait_cols, mm_cols)

# Define empty list to hold gene info tables
all_geneInfo <- list()

for (mod in modNames) {
  cat("Processing module:", mod, "\n")
  
  moduleGenes <- mergedColors == mod
  probes <- names(expression.data)[moduleGenes]
  probes2annot <- match(probes, all_species_annot$Orthogroup)

  geneInfo0 <- data.frame(
    Orthogroup = probes,
    geneSymbol = all_species_annot$`Protein names`[probes2annot],
    LocusLinkID = all_species_annot$Entry[probes2annot],
    moduleColor = mod,
    stringsAsFactors = FALSE
  )

  # Add gene significance and p-value columns
  for (i in seq_along(traitList)) {
    trait <- traitList[i]
    gs_col <- paste0("GS.", trait)
    pgs_col <- paste0("p.GS.", trait)

    geneInfo0[[gs_col]] <- geneTraitSignificance[moduleGenes, i]
    geneInfo0[[pgs_col]] <- GSPvalue[moduleGenes, i]
  }

  # Add module membership only for *this* module
  mm_col <- paste0("MM.", mod)
  pmm_col <- paste0("p.MM.", mod)
  geneInfo0[[mm_col]] <- geneModuleMembership[moduleGenes, mm_col]
  geneInfo0[[pmm_col]] <- MMPvalue[moduleGenes, pmm_col]

  # Store this table
  all_geneInfo[[mod]] <- geneInfo0
}

# ---- FINAL GENE INFO TABLE ----
final_geneInfo <- do.call(rbind, all_geneInfo)

# Optional: Save to CSV
# write.csv(final_geneInfo, "geneInfo_correlated_to_bacteria.csv", row.names = FALSE)

# ---- FIND SIGNIFICANT CORRELATIONS TO BACTERIA (p < 0.1) ----
library(tidyr)
library(dplyr)

sig_genes_long <- final_geneInfo %>%
  select(Orthogroup, geneSymbol, starts_with("GS."), starts_with("p.GS.")) %>%
  pivot_longer(
    cols = starts_with("p.GS."),
    names_to = "Trait",
    values_to = "p_value"
  ) %>%
  mutate(Trait = gsub("p.GS.", "", Trait)) %>%
  filter(p_value < 0.1)

# ---- FREQUENCY TABLE OF SIGNIFICANT GENE APPEARANCES ----
freq_table <- sig_genes_long %>%
  group_by(geneSymbol) %>%
  summarise(Frequency = n()) %>%
  arrange(desc(Frequency)) %>%
  filter(!is.na(geneSymbol))

# ---- PLOT TOP FUNCTIONALLY REDUNDANT GENES ----
library(ggplot2)
ggplot(freq_table[1:100, ], aes(x = reorder(geneSymbol, Frequency), y = Frequency)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Most Frequently Correlated Genes Across Bacteria",
       x = "Protein Name", y = "Count of Significant Correlations")

freq_table[1:100, ]$geneSymbol
```

```{r}
# ---- FIND SIGNIFICANT CORRELATIONS TO BACTERIA (p < 0.1) ----
library(tidyr)
library(dplyr)

# Pivot to long format, filter for significant correlations (p < 0.1)
sig_genes_long <- final_geneInfo %>%
  select(Orthogroup, geneSymbol, starts_with("GS."), starts_with("p.GS.")) %>%
  pivot_longer(cols = starts_with("p.GS."),
               names_to = "Trait",
               values_to = "p_value") %>%
  mutate(Trait = gsub("p.GS.", "", Trait)) %>%
  filter(p_value < 0.1)

# Attach GO terms and Entry ID
sig_genes_long <- sig_genes_long %>%
  left_join(all_species_annot %>%
              select(Orthogroup, Entry, `Gene ontology IDs`),
            by = "Orthogroup")

# ---- CREATE TABLE OF GENES AND THEIR CORRELATED BACTERIA ----
# This will show each significant gene–bacteria correlation
sig_gene_bacteria_table <- sig_genes_long %>%
  select(geneSymbol, Trait, Entry, `Gene ontology IDs`) %>%
  distinct() %>%
  arrange(geneSymbol)

# Optional: Save table for GOMWU
# write.table(sig_gene_bacteria_table, "sig_gene_bacteria_GO.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# ---- PLOT GENE-BY-BACTERIA COUNTS (instead of just frequency) ----
# This gives number of *bacteria* traits each gene was significantly correlated with
freq_table <- sig_gene_bacteria_table %>%
  group_by(geneSymbol) %>%
  summarise(Significant_Bacteria = paste(unique(Trait), collapse = ", "),
            Count = n()) %>%
  arrange(desc(Count)) %>%
  filter(!is.na(geneSymbol))

# ---- PLOT ----
ggplot(freq_table[1:100, ], aes(x = reorder(geneSymbol, Count), y = Count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Genes Correlated with Multiple Bacteria",
       x = "Protein Name", y = "Number of Significant Bacterial Correlations")

# Output table
head(freq_table, 20)


write.csv(freq_table, file="Host_HV_freqtabletoPathogens.csv")

```
The first few rows of what freq_table looks like:
![image](https://github.com/user-attachments/assets/32a2899d-6515-4c89-a2a3-283a7a83f049)


### Include GS values
```{r}
# ---- LOAD LIBRARIES ----
library(dplyr)
library(tidyr)

# ---- STEP 1: Filter for Significant Correlations ----
gs_sig_long <- final_geneInfo %>%
  select(geneSymbol, starts_with("GS."), starts_with("p.GS.")) %>%
  pivot_longer(cols = starts_with("p.GS."),
               names_to = "Trait",
               values_to = "p_value") %>%
  mutate(Trait = gsub("p.GS.", "", Trait)) %>%
  filter(p_value < 0.1)

# Attach GS values corresponding to each Trait
gs_sig_long <- gs_sig_long %>%
  rowwise() %>%
  mutate(GS_value = final_geneInfo[[paste0("GS.", Trait)]][match(geneSymbol, final_geneInfo$geneSymbol)]) %>%
  ungroup()

# Remove potential duplicates (if any)
gs_sig_long_clean <- gs_sig_long %>%
  distinct(geneSymbol, Trait, GS_value)

# ---- STEP 2: Pivot to Wide Format ----
sig_gs_wide <- gs_sig_long_clean %>%
  pivot_wider(names_from = Trait, values_from = GS_value)

# ---- STEP 3: Add Frequency Column ----
gene_freq <- gs_sig_long_clean %>%
  group_by(geneSymbol) %>%
  summarise(Frequency = n(), .groups = "drop")

# ---- STEP 4: Merge and Rank ----
sig_gs_wide_ranked <- gene_freq %>%
  left_join(sig_gs_wide, by = "geneSymbol") %>%
  arrange(desc(Frequency))

# View top rows
sig_gs_wide_ranked


# ---- Optional: Save ----
#write.csv(sig_gs_wide_ranked, "Host_HV_freqtabletoPathogens_GS.csv")

```
The first few rows of what sig_gs_wide_ranked looks like:
![image](https://github.com/user-attachments/assets/52386ff7-5270-402a-a056-da3bf8caa604)


## WGCNA - LS_metadata
```{r, fig.width=14, fig.height=6}
#install.packages('BiocManager')
#library(BiocManager)
#BiocManager::install('WGCNA')
#BiocManager::install('flashClust')

library(WGCNA)
library(flashClust)
library(curl)

head(LS_metadata)

# Retaining only Expression data
expression.data <- LS_metadata[,12:length(LS_metadata)] #removing variables not holding expression data
row.names(expression.data) <- LS_metadata[,1]
names(expression.data)


# Identifying Good Genes
gsg <-goodSamplesGenes(expression.data)
summary(gsg)

#If the allOK object returns true, which it does in this case, there are no outliers present. If it doesn’t return true, we need to filter out the outliers manually using the following code --> See "https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html" for WGCNA tutorial.
gsg$allOK


# Another way to identify outliers is to visualize by hierarchical tree
sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

# You can remove the outlier using a cutree function.
#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 50, col = "red"); # 50 is an arbitrart cutoff, however it is visually guided.

# 43 and 57 appear to be outliers and can be removed by setting the cut height to 50 here:
#cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 50, minSize = 10) #returns numeric vector
#Remove outlier
#expression.data <- expression.data[cut.sampleTree==1, ]


#Pick the Soft Threshold Power
spt <- pickSoftThreshold(expression.data) 
spt

# Plot the r2 values as a function of the soft thresholds
#We should be maximizing the r2 value and minimizing mean connectivity.
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

# Plot the Mean Connectivity (This was not working in rmakrdonw I needed to copy and paiste into the console to work. Odd)
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

#You can determine the soft power threshold should be set to 5 as it is the spt that retains the highest mean connectivity while reaching an r2 value above 0.80.
#NOTE: the higher the value, the stronger the connection strength will be of highly correlated gene expression profiles and the more devalued low correlations will be.

# Calling the Adjacency Function
#Now that you have the soft threshold power determined you can call on the adjacency() function of the WGCNA package.

#REMINDER: This function calculates the similarity measurement and transforms the similarity by the adjacency function and generates a weighted network adjacency matrix.

softPower <- 15
adjacency <- adjacency(expression.data, power = softPower)

# Module Construction
TOM <- TOMsimilarity(adjacency)

# To convert this matrix into a dissimilarity matrix you can subtract the TOM object from 1.
TOM.dissimilarity <- 1-TOM

#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
labels = FALSE, hang = 0.04)

# Identify Modules - Set ClusterSize 30 is default.
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 10)

table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module. 

#Plot the module assignment under the gene dendrogram for visualization.
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")


# Module Eigengene Identification
#A ME (Module Eigengene) is the standardized gene expression profile for a given module.(Summary statistic of Expression)
#To identify the Module Eigengene you can call on the expression data into the moduleEigengenes() function.
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

# Module Merging (optional)

#To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations and merge the modules that have similar expression profiles.

#REMINDER: An eigengene is the gene whose expression is representative of the the majority of genes expressed within a module.

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

#Plot
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75


merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)

# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# The similar modules are now merged! Let’s compare them with the original modules by creating another dendrogram
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
c("Original Module", "Merged Module"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors for original and merged modules")


# External Trait Matching
#Once you have constructed the network (the adjacency matrix) and divided it into meaningful modules, you can begin to relate the network to external traits.

#Start by reading in Trait Data Table
#allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA.csv") # I took the file MetaData.csv and converted some of the categorical phenotype data into continuous data and added relative risk and disease prevalence. We can also add bacteria abundances, bacteria gene expression, and algal gene expression in time.

# Combine module eigengenes with traits
#allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA_ME.csv")
# Load metadata and remove first two columns (assumed to be unnamed row numbers or irrelevant columns)
#allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA_ME_Bac16s_Net.csv")
allTraits <- read.csv("~/Desktop/Microbial Metatranscriptomics/R/MetaData_WGCNA_ME_Bac16s_Net_5.19.csv")


# Keep SampleID and all trait columns
#allTraits <- allTraits[, -c(1)]  # You can remove this line if SampleID is not duplicated or has been moved up

# Extract SampleIDs from expression data to ensure matching
Samples <- rownames(expression.data)

# Make sure SampleID column exists in allTraits (sometimes it's rownames instead)
if (!"SampleID" %in% colnames(allTraits)) {
  allTraits <- allTraits %>% rownames_to_column("SampleID")
}

# Keep only the traits for matched samples (in expression data)
datTraits <- allTraits %>% filter(SampleID %in% Samples)
rownames(datTraits) <- datTraits$SampleID

# Add SampleID column to mergedMEs if it's not already there
if (!"SampleID" %in% colnames(mergedMEs)) {
  mergedMEs <- mergedMEs %>% rownames_to_column("SampleID")
}

# Perform clean merge
merged_data <- mergedMEs %>%
  left_join(datTraits, by = "SampleID")


# Remove NAs caused by missing eigengenes
merged_data_clean <- merged_data
#merged_data_clean <- na.omit(merged_data)

# Separate back into matrices
mergedMEs_clean <- merged_data_clean %>% select(starts_with("ME")) %>% as.data.frame()
datTraits_clean <- merged_data_clean %>% select(-SampleID, -starts_with("ME")) %>% as.data.frame()

# Run correlations
nSamples <- nrow(mergedMEs_clean)
module.trait.correlation <- cor(mergedMEs_clean, datTraits_clean, use = "p")
module.trait.Pvalue <- corPvalueStudent(module.trait.correlation, nSamples)

# Prepare text matrix for heatmap
textMatrix <- paste(signif(module.trait.correlation, 2), "\n(",
                    signif(module.trait.Pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module.trait.correlation)

# Plot heatmap
#par(mar = c(15, 12, 0.1, 5))
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits_clean),
               yLabels = names(mergedMEs_clean),
               ySymbols = names(mergedMEs_clean),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Lineage-Specific Host Module-trait relationships"))


pdf("Host_LS_WGCNAHeatmap.pdf", width = 18, height = 8)
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits_clean),
               yLabels = names(mergedMEs_clean),
               ySymbols = names(mergedMEs_clean),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Lineage-Specific Host Module-trait relationships"))
dev.off()

```

![image](https://github.com/user-attachments/assets/e3cccabb-7d1c-4010-b436-3f69a0b4dc47)


