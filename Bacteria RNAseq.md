<div align="center">
  <img width="75" style="margin: 0 10px;" src="https://github.com/user-attachments/assets/b17ae489-cb10-46f0-992e-d436cc4ac872" alt="image" />
  <img width="150" style="margin: 0 10px;" src="https://github.com/user-attachments/assets/c1ffc07b-a4b2-4769-95d4-f6cc1b1e2d74" alt="image" />
  <img width="75" style="margin: 0 10px;" src="https://github.com/user-attachments/assets/aba7fade-7367-4282-ab30-92096ea1afa8" alt="image" />
  
</div>

# Bacteria Orthogroup RNAseq Analysis

> [!WARNING]
> This Github page is a work in progress. All of the current code is valid, however there is more to be added!

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria")
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

## Obtain Orthologs from Orthofinder

To obtain orthologs, we filter our Orthofinder results to only contain
Orthogroups with at least one sequence per species.

This can be done in excel like so: 
1. Open Orthogroups.tsv in Excel
2. Click *Find & Select*
3. Click *Go to Special*
4. Choose *Blanks*
5. Click OK and then all the *blank rows/cells will be highlighted*
6. Choose the *Delete under Cells* section on the Home Tab
7. Click *Delete Sheet Rows*
8. Save as "host_shared_orthogroups.csv"
9. Steps 1-8 were repeated 4 times so that all blank rows were eventually removed. I think the amount of data excel needed to process required these steps to be
repeated. You can ensure you have the accurate final number of orthogroups by comparing your row count to the "Statistics_Overall" file
in /Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/Comparative_Genomics_Statistics/
specifically the value in row "Number of orthogroups with all species present"

### Make tx2gene files using shared orthogroups

For this analysis we will use "Orthogroup_Sequences.csv" from the
"Command_Line_Code.rmd" file to create new tx2gene files for each
species

```{r Host Orthogroups}

Orthogroups <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/Bacteria_shared_orthogroups_SingleBestORF.csv")

# Acer Orthogroups
Acer_orthogroups <- Orthogroups[,c(1,2)]
Acer_orthogroups <- Acer_orthogroups %>% separate_rows(Acer_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Acer_orthogroups <- Acer_orthogroups[,c(2,1)]
colnames(Acer_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Acer_bac_Orthogroup_Transcripts <- Acer_orthogroups[c(1)]
names(Acer_bac_Orthogroup_Transcripts)
write.table(Acer_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Acer_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Acer_orthogroups$Transcript <- gsub("\\..*","",Acer_orthogroups$Transcript)
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/AcerBac_orthogroup_tx2gene.csv",row.names = FALSE)

# Mcav Orthogroups
Mcav_orthogroups <- Orthogroups[,c(1,3)]
Mcav_orthogroups <- Mcav_orthogroups %>% separate_rows(Mcav_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Mcav_orthogroups <- Mcav_orthogroups[,c(2,1)]
colnames(Mcav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Mcav_bac_Orthogroup_Transcripts <- Mcav_orthogroups[c(1)]
names(Mcav_bac_Orthogroup_Transcripts)
write.table(Mcav_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Mcav_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Mcav_orthogroups$Transcript <- gsub("\\..*","",Mcav_orthogroups$Transcript)
Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/McavBac_orthogroup_tx2gene.csv",row.names = FALSE)

# Ofav Orthogroups
Ofav_orthogroups <- Orthogroups[,c(1,4)]
Ofav_orthogroups <- Ofav_orthogroups %>% separate_rows(Ofav_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Ofav_orthogroups <- Ofav_orthogroups[,c(2,1)]
colnames(Ofav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Ofav_bac_Orthogroup_Transcripts <- Ofav_orthogroups[c(1)]
names(Ofav_bac_Orthogroup_Transcripts)
write.table(Ofav_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Ofav_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Ofav_orthogroups$Transcript <- gsub("\\..*","",Ofav_orthogroups$Transcript)
Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/OfavBac_orthogroup_tx2gene.csv",row.names = FALSE)

# Past Orthogroups
Past_orthogroups <- Orthogroups[,c(1,5)]
Past_orthogroups <- Past_orthogroups %>% separate_rows(Past_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Past_orthogroups <- Past_orthogroups[,c(2,1)]
colnames(Past_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Past_bac_Orthogroup_Transcripts <- Past_orthogroups[c(1)]
names(Past_bac_Orthogroup_Transcripts)
write.table(Past_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Past_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Past_orthogroups$Transcript <- gsub("\\..*","",Past_orthogroups$Transcript)
Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/PastBac_orthogroup_tx2gene.csv",row.names = FALSE)

```

### Read in Annotated Transcripts
```{r}
# In terminal: 
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/*_Bac_orthologs_annotated_e-5.txt /Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/
  
# Load required libraries
library(dplyr)
library(readr)

# List of input file prefixes
file_prefixes <- c("Acer", "Mcav", "Ofav", "Past")

# Function to process a single file
process_file <- function(prefix) {
  # Construct file paths
  input_file <- paste0(prefix, "_Bac_orthologs_annotated_e-5.txt")
  output_file <- paste0(prefix, "_Bac_orthologs_annotated_e-5_formatted.txt")
  
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

## Annotate Orthogroups
```{r}
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/AcerBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/McavBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/OfavBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/PastBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)

```

### Import Salmon Counts
Read in the raw counts from Salmon for all coral samples and assign Orthogroups
In terminal:
```{r}
# Acer
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Acer_bac_ORF

# Mcav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Mcav_bac_ORF

# Ofav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Ofav_bac_ORF

# Past
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Past_bac_ORF
```

### A. cervicornis Salmon Quants
```{r,results='hide',tidy=TRUE}

setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/bacteria/Acer/")
Acer_samples <- read.csv("~/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Acer/Acer_samples.csv", header = TRUE)
Acer_files <- file.path(Acer_samples$sample, "quant.sf")
names(Acer_files) <- paste0("sample", 1:13)
all(file.exists(Acer_files))
# Do not continue if output shows [FALSE]
Acer_tx2gene <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/bacteria/SingleBestORF/AcerBac_orthogroup_tx2gene_annot.csv")
Acer_txi <- tximport(Acer_files, type = 'salmon',tx2gene=Acer_tx2gene)
dim(Acer_txi$counts)

# 1313 Orthogroups Single Best ORF

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

### M. cavernosa Salmon Quants
```{r,results='hide',tidy=TRUE}
setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/bacteria/Mcav/")
Mcav_samples <- read.csv("~/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Mcav/Mcav_samples.csv", header = TRUE)
Mcav_files <- file.path(Mcav_samples$sample, "quant.sf")
names(Mcav_files) <- paste0("sample", 1:18)
all(file.exists(Mcav_files))
# Do not continue if output shows [FALSE]
Mcav_tx2gene <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/bacteria/SingleBestORF/McavBac_orthogroup_tx2gene_annot.csv")
Mcav_txi <- tximport(Mcav_files, type = 'salmon',tx2gene=Mcav_tx2gene)
dim(Mcav_txi$counts)

# 1313 Orthogroups Single Best ORF

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

### O. faveolata Salmon Quants
```{r,results='hide',tidy=TRUE}

setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/bacteria/Ofav/")
Ofav_samples <- read.csv("~/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Ofav/Ofav_samples.csv", header = TRUE)
Ofav_files <- file.path(Ofav_samples$sample, "quant.sf")
names(Ofav_files) <- paste0("sample", 1:29)
all(file.exists(Ofav_files))
# Do not continue if output shows [FALSE]
Ofav_tx2gene <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/OfavBac_orthogroup_tx2gene_annot.csv")
Ofav_txi <- tximport(Ofav_files, type = 'salmon',tx2gene=Ofav_tx2gene)
dim(Ofav_txi$counts)

# 1313 Orthogroups Single Best ORF

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

### P. astreoides Salmon Quants
```{r,results='hide',tidy=TRUE}

setwd("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/salmon/bacteria/Past/")
Past_samples <- read.csv("~/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Past/Past_samples.csv", header = TRUE)
Past_files <- file.path(Past_samples$sample, "quant.sf")
names(Past_files) <- paste0("sample", 1:17)
all(file.exists(Past_files))
# Do not continue if output shows [FALSE]
Past_tx2gene <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/bacteria/SingleBestORF/PastBac_orthogroup_tx2gene_annot.csv")
Past_txi <- tximport(Past_files, type = 'salmon',tx2gene=Past_tx2gene)
dim(Past_txi$counts)

# 1313 Orthogroups Single Best ORF

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
Past_counts
```

Merge the count data
```{r,results='hide',tidy=TRUE}
Acer_Mcav <- merge(Acer_counts,Mcav_counts, by ="Orthogroup")
Acer_Mcav_Past <- merge(Acer_Mcav,Past_counts,by="Orthogroup")
Acer_Mcav_Past_Ofav <- merge(Acer_Mcav_Past,Ofav_counts,by="Orthogroup")
#CladeA_CladeB_CladeC_CladeD <- merge(CladeA_CladeB_CladeC,CladeD_counts,by="Orthogroup")
#all_species <- merge(CladeA_CladeB_CladeC_CladeD,pstr_counts,by="Orthogroup")
all_species <- Acer_Mcav_Past_Ofav

# 1,313 Orthogroups Single Best ORF
```

## Annotating Orthogroups based on the most Common Entry
```{r}
# Load required libraries
library(dplyr)

# Define file paths for the annotated files
files <- c("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/AcerBac_orthogroup_tx2gene_annot.csv",
           "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/McavBac_orthogroup_tx2gene_annot.csv",
           "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/OfavBac_orthogroup_tx2gene_annot.csv",
           "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/PastBac_orthogroup_tx2gene_annot.csv")

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
          file = "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/Orthogroups_with_Representative_Entry.csv",
          row.names = FALSE)

# Save representative entries only
write.csv(representative_entries,
          file = "~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/Representative_Entries.csv",
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

# DESeq2

Convert counts to integers. Run DESeq2 on Orthogroups to get Rlog
normalized expression.

```{r,results='hide',tidy=TRUE, fig.height=6, fig.width=4}

countData <- all_species_annot
countData <- tibble::column_to_rownames(countData,"Orthogroup")

# Remove some annotation column names to run DESeq2
colnames(countData)
countData <- countData[,c(9:length(countData))]


colData <- read.csv("/Users/nicholas.macknight/Desktop/Autumn16s/Autumn16s_05082023/230613_M02476_0578_000000000-L3NTN/20230615_222347/metadata.csv")

# So countData and colData have some slight name formatting differences that we need to make consistent.

# Remove "_host_quant" from the end of the column names
colnames(countData) <- gsub("_Bacteria_quant$", "", colnames(countData))

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


# View the subsetted colData
colData[[1]]
colnames(countData)
colData_subset
colData_subset[[1]]

# Remove the column named "PastPA3Disease_5_2" from countData
countData <- countData[, !colnames(countData) %in% "PastPA3Disease_5_2"]


colData_subset[[1]]
colData
```

### ~Coral.Species + Treatment

```{r}
# Set the design
dds_treatment <- DESeqDataSetFromMatrix(countData = countData, 
                                        colData = colData_subset, 
                                        design = ~Coral.Species + Treatment )

# Filter out columns with an average number of reads less than 10
dds_treatment <- dds_treatment[ rowMeans(counts(dds_treatment)) > 10, ] 

# Run DESeq2
test <- DESeq(dds_treatment)
rld <- rlog(dds_treatment, blind=FALSE) 
rrld <- assay(rld)
# 5,125 after filtering out low counts
# 371 after filtering out low counts
write.csv(rrld, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/treatment_rlog_SingleBestORF.csv")  


resCvD <- results(test,contrast=c("Treatment", "Control", "Disease"))
head(resCvD)
resorderdCvD <-resCvD[order(resCvD$padj),] 
head(resorderdCvD, 12)
summary(resorderdCvD)

# Annotate and save each results table
annotate_and_save(resorderdCvD, 
                  "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_CvD_BestORF_Annotated.csv", 
                  all_species_annot)

write.csv(resorderdCvD, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_CvD_BestORF.csv")

#subset <- read.csv("Resistant_subset.csv") #EDA counts
DEG_EDA <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_CvD_BestORF.csv") #DEG output
colnames(DEG_EDA)[colnames(DEG_EDA)=="X"] <- "Orthogroups"
countData$Orthogroups <- rownames(countData)

DEG_EDA <- merge(DEG_EDA,countData,by="Orthogroups")
```
<pre>
> summary(resorderdCvD)
out of 739 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 9, 1.2%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 26, 3.5%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
</pre>


### ~Coral.Species + Outcome
```{r}
# Set the design
dds_Outcome <- DESeqDataSetFromMatrix(countData = countData, 
                                        colData = colData_subset, 
                                        design = ~Coral.Species + Outcome )

# Filter out columns with an average number of reads less than 10
dds_Outcome <- dds_Outcome[ rowMeans(counts(dds_Outcome)) > 10, ] 

# Run DESeq2
test <- DESeq(dds_Outcome)

resO.CvD <- results(test,contrast=c("Outcome", "Control", "Disease"))
resO.CvH <- results(test,contrast=c("Outcome", "Control", "Healthy"))
resO.HvD <- results(test,contrast=c("Outcome", "Healthy", "Disease"))

head(resO.CvD)
resorderdO.CvD <-resO.CvD[order(resO.CvD$padj),] 
head(resorderdO.CvD, 12)
summary(resorderdO.CvD)
write.csv(resorderdO.CvD, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.CvD_BestORF.csv")

head(resO.CvH)
resorderdO.CvH <-resO.CvH[order(resO.CvH$padj),] 
head(resorderdO.CvH, 12)
summary(resorderdO.CvH)
write.csv(resorderdO.CvH, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.CvH_LORF.csv")

head(resO.HvD)
resorderdO.HvD <-resO.HvD[order(resO.HvD$padj),] 
head(resorderdO.HvD, 12)
summary(resorderdO.HvD)
write.csv(resorderdO.HvD, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.HvD_LORF.csv")

#subset <- read.csv("Resistant_subset.csv") #EDA counts
DEG_HvD <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.CvH_LORF.csv") #DEG output
colnames(DEG_HvD)[colnames(DEG_HvD)=="X"] <- "Orthogroups"
countData$Orthogroups <- rownames(countData)

DEG_HvD <- merge(DEG_HvD,countData,by="Orthogroups")











# Load required libraries
library(dplyr)
library(tibble)

# Set the design
dds_Outcome <- DESeqDataSetFromMatrix(countData = countData, 
                                      colData = colData_subset, 
                                      design = ~Coral.Species + Outcome)

# Filter out rows with an average number of reads less than 10
dds_Outcome <- dds_Outcome[rowMeans(counts(dds_Outcome)) > 10, ]

# Run DESeq2
test <- DESeq(dds_Outcome)

# Generate results for Outcome contrasts
resO.CvD <- results(test, contrast = c("Outcome", "Control", "Disease"))
resO.CvH <- results(test, contrast = c("Outcome", "Control", "Healthy"))
resO.HvD <- results(test, contrast = c("Outcome", "Healthy", "Disease"))

# Order results by adjusted p-value
resorderdO.CvD <- resO.CvD[order(resO.CvD$padj), ]
resorderdO.CvH <- resO.CvH[order(resO.CvH$padj), ]
resorderdO.HvD <- resO.HvD[order(resO.HvD$padj), ]

# Peak at the results
head(resorderdO.CvD, 12)
summary(resorderdO.CvD)

head(resorderdO.CvH, 12)
summary(resorderdO.CvH)

head(resorderdO.HvD, 12)
summary(resorderdO.HvD)

# Load annotation file
annotation_file <- all_species_annot  # Update with actual path
all_species_annot <- read.csv(annotation_file)

# Ensure the annotation file has columns: "Orthogroup" and "Protein.names"
colnames(all_species_annot) <- c("Orthogroup", "Protein.names")

# Function to annotate and save results
annotate_and_save <- function(res_ordered, output_file, annotation_data) {
  annotated <- as.data.frame(res_ordered) %>%
    rownames_to_column("Orthogroup") %>%
    left_join(annotation_data, by = "Orthogroup") %>%
    column_to_rownames("Orthogroup")
  
  # Write annotated results to file
  write.csv(annotated, file = output_file, row.names = TRUE)
}

# Annotate and save each results table
annotate_and_save(resorderdO.CvD, 
                  "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.CvD_BestORF_Annotated.csv", 
                  all_species_annot)

annotate_and_save(resorderdO.CvH, 
                  "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.CvH_LORF_Annotated.csv", 
                  all_species_annot)

annotate_and_save(resorderdO.HvD, 
                  "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.HvD_LORF_Annotated.csv", 
                  all_species_annot)

# Process DEG_HvD data
DEG_HvD <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_O.HvD_LORF_Annotated.csv")

# Rename column for consistency
colnames(DEG_HvD)[colnames(DEG_HvD) == "X"] <- "Orthogroups"

# Add Orthogroups to countData
countData$Orthogroups <- rownames(countData)

# Merge DEG_HvD with countData by Orthogroups
DEG_HvD <- merge(DEG_HvD, countData, by = "Orthogroups")

# Save the final merged file
write.csv(DEG_HvD, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_HvD_Annotated.csv", row.names = FALSE)

```

<pre>
> summary(resorderdO.CvD)
out of 739 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 12, 1.6%
LFC < 0 (down)     : 4, 0.54%
outliers [1]       : 23, 3.1%
low counts [2]     : 171, 23%
(mean count < 8)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(resorderdO.CvH)
out of 739 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 30, 4.1%
LFC < 0 (down)     : 13, 1.8%
outliers [1]       : 23, 3.1%
low counts [2]     : 227, 31%
(mean count < 10)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(resorderdO.HvD)
out of 739 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1, 0.14%
LFC < 0 (down)     : 2, 0.27%
outliers [1]       : 23, 3.1%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

</pre>

### ~Treatment + Coral.Species

```{r}
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
write.csv(rlog_data, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/rlog_normalized_counts_BySpecies_BestORF.csv")

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

# Save results for each comparison
write.csv(res_Acer_vs_Mcav[order(res_Acer_vs_Mcav$padj), ], 
          file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_Acer_vs_Mcav.csv")

write.csv(res_Acer_vs_Ofav[order(res_Acer_vs_Ofav$padj), ], 
          file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_Acer_vs_Ofav.csv")

write.csv(res_Acer_vs_Past[order(res_Acer_vs_Past$padj), ], 
          file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_Acer_vs_Past.csv")

write.csv(res_Mcav_vs_Ofav[order(res_Mcav_vs_Ofav$padj), ], 
          file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_Mcav_vs_Ofav.csv")

write.csv(res_Mcav_vs_Past[order(res_Mcav_vs_Past$padj), ], 
          file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_Mcav_vs_Past.csv")

write.csv(res_Ofav_vs_Past[order(res_Ofav_vs_Past$padj), ], 
          file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/DEG_Ofav_vs_Past.csv")

# Merge rlog-transformed data with DEG results if needed
rlog_data_df <- as.data.frame(rlog_data)
rlog_data_df$Orthogroups <- rownames(rlog_data_df)

res_Acer_vs_Mcav_df <- as.data.frame(res_Acer_vs_Mcav)
res_Acer_vs_Mcav_df$Orthogroups <- rownames(res_Acer_vs_Mcav_df)

merged_Acer_vs_Mcav <- merge(res_Acer_vs_Mcav_df, rlog_data_df, by = "Orthogroups")




res_CvD_test <- results(test_Species, contrast = c("Treatment", "Control", "Disease"))
summmary(res_CvD_test)

write.csv(res_CvD_test, file="res_CvD_test.csv")
#subset <- read.csv("Resistant_subset.csv") #EDA counts
DEG_CvD_test <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/res_CvD_test.csv") #DEG output
colnames(DEG_CvD_test)[colnames(DEG_CvD_test)=="X"] <- "Orthogroups"
countData$Orthogroups <- rownames(countData)

DEG_CvD_test <- merge(DEG_CvD_test,countData,by="Orthogroups")


```
<pre>
> summary(res_Acer_vs_Mcav)
out of 704 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 136, 19%
LFC < 0 (down)     : 333, 47%
outliers [1]       : 24, 3.4%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(res_Acer_vs_Ofav)
out of 704 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 146, 21%
LFC < 0 (down)     : 270, 38%
outliers [1]       : 24, 3.4%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(res_Acer_vs_Past)
out of 704 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 23, 3.3%
LFC < 0 (down)     : 44, 6.2%
outliers [1]       : 24, 3.4%
low counts [2]     : 28, 4%
(mean count < 5)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(res_Mcav_vs_Ofav)
out of 704 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 265, 38%
LFC < 0 (down)     : 131, 19%
outliers [1]       : 24, 3.4%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(res_Mcav_vs_Past)
out of 704 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 317, 45%
LFC < 0 (down)     : 112, 16%
outliers [1]       : 24, 3.4%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(res_Ofav_vs_Past)
out of 739 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 217, 29%
LFC < 0 (down)     : 155, 21%
outliers [1]       : 26, 3.5%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
</pre>


### ~ No Design

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
write.csv(rlog_metadata_nodesign, file = "~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/rlog_metadata_nodesign.csv", row.names = FALSE)

# Display a preview
head(rlog_metadata_nodesign)



# Run DESeq2
test_NoDesign <- DESeq(dds)

# Perform rlog normalization
rld_NoDesign <- rlog(test_NoDesign, blind = FALSE)

# Save rlog-transformed data for downstream analysis
rld_NoDesign_data <- assay(rld_NoDesign)
write.csv(rld_NoDesign_data, file = "/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/rld_NoDesign_data.csv")

res_CvD_NoDesign <- results(test_NoDesign, contrast = c("Treatment", "Control", "Disease"))

summary(res_CvD_NoDesign)


```

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
colData_subset$Treatment <- factor(colData_subset$Treatment)
colData_subset$Coral.Genotype <- factor(colData_subset$Coral.Genotype)

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
  colnames(countData) <- gsub("_Bacteria_quant$", "", colnames(countData))
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
  sig_orthogroups <- rownames(res)[which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.1)]
  
  # Store the significant orthogroups for this species
  significant_orthogroups[[species]] <- sig_orthogroups
}

# Print significant orthogroups for each species
print(significant_orthogroups)


```
$Acer
[1] "OG0000020" "OG0000353" "OG0000509"

$Mcav
character(0)

$Ofav
 [1] "OG0000018" "OG0000070" "OG0000101" "OG0000123" "OG0000534" "OG0000570" "OG0000694" "OG0000769" "OG0000811" "OG0000818" "OG0001380"

$Past
character(0)


### Visualize Overlaps Using a Venn Diagram
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

![image](https://github.com/user-attachments/assets/31c7f75e-b7f3-4212-8396-afe72670bb40)


# EVE

```{r EVE_phyloTree,fig.align='center',fig.height=1.5,fig.width=1.5,dpi=100}
library(ape)


exprTbl <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/treatment_rlog_SingleBestORF.csv")
colnames(exprTbl)[colnames(exprTbl)=="X"] <- "Orthogroup"
# The table needs to be converted to a matrix:
exprMat <- as.matrix(exprTbl[,-1])
rownames(exprMat) <- exprTbl$Orthogroup
exprMat
dim(exprMat)  
# 739 x 76 - Best ORF

# Species phylogeny can be obtained from OrthoFinder output in DIR=/OrthoFinder/Results/Species_Tree/SpeciesTree_rooted.txt

speciesTree <-read.tree(text="((Mcav:0.14721,Acer:0.102176)0.430312:0.0317775,(Ofav:0.164205,Past:0.0930845)0.430312:0.0317775);") # Single Best ORF

# Check validity of species tree
plot(speciesTree,
     type = "phylogram",
     use.edge.length = TRUE,
     node.pos = NULL,
     show.tip.label = TRUE,
     root.edge = TRUE)


#ggsave(file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Bacteria/spciesTree.pdf", tree)

```

<div align="center">
  <img src="https://github.com/user-attachments/assets/e27dc222-4efb-4758-bd75-c29a99e583fc" alt="image" />
</div>


EVE needs to know which species each column in the expression matrix
belongs to. This is done by creating a vector of species names
corresponding to the tip labels in the species tree

```{r,results='hide',tidy=TRUE}
speciesTree[[1]]$tip.label
hostTree[[1]]$tip.label
colnames(exprMat)
colSpecies <- substr(colnames(exprMat), 1, 4)
colSpecies

library(evemodel)
# Run EVE
EVE.BestORF <- betaSharedTest(tree = hostTree[[1]], gene.data = exprMat, colSpecies = colSpecies)

```

```{r,results='hide',tidy=TRUE}
# P-value can then be calculated using:
pval = pchisq(EVE.BestORF$LRT,df = 1,lower.tail = F)

# The shared beta:
EVE.BestORF$sharedBeta
# [1] 0.4121188 # BEst ORF

log(0.4121188)
# [1] -0.8864436 # Best ORF

#Combine LRT, beta, theta, sigma2, alpha:
head(cbind(EVE.BestORF$LRT,EVE.BestORF$indivBetaEVE$par,pval))
colnames(exprMat)
dim(exprMat)
# [1]  371   76
a <- cbind(EVE.BestORF$LRT,EVE.BestORF$indivBetaRes$par,pval,exprMat)
#rownames(a) <- NodesignTbl.t$Orthogroup

# Give Column Name to LRT column
colnames(a)
colnames(a)[colnames(a)==""] <- "LRT"
head(a)
EVE_results <- as.data.frame(a)
EVE_results <- tibble::rownames_to_column(EVE_results, "Orthogroup")

EVE_results$type <- ifelse(EVE_results$beta<0.4121188,"Lineage Specific","Highly Variable")
EVE_results$significant <- ifelse(EVE_results$pval<=0.05,"Significant","Not Significant")
EVE_results$category <- ifelse(EVE_results$significant == "Significant",EVE_results$type,"NS")
colnames(EVE_results)
EVE_results <- EVE_results %>% relocate(type, .after = pval)
EVE_results <- EVE_results %>% relocate(significant, .after = type)
EVE_results <- EVE_results %>% relocate(category, .after = significant)

write.csv(EVE_results, file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/EVE_results_Bacteria_Orthogroups_BestORF_HostTree.csv",row.names = FALSE)

#Visualize LRT v Beta by volcano plot for gene data
Shared_EVE <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/EVE_results_Bacteria_Orthogroups_BestORF_HostTree.csv")
Shared_EVE_sig <- Shared_EVE %>% filter(pval <= 0.05)
write.csv(Shared_EVE_sig, file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Signficant_EVE_results_Bacteria_Orthogroups_BestORF.csv")
# Make a basic volcano plot
Shared_EVE_LS <- Shared_EVE_sig %>% filter(log(beta)< (-0.8864436))
dim(Shared_EVE_LS)  #37
Shared_EVE_HV <- Shared_EVE_sig %>% filter(log(beta)> (-0.8864436))
dim(Shared_EVE_HV)  #112
  #  total
```


```{r EVE_volcano,fig.height=3,fig.width=3,dpi=300}
p <- ggplot(data = Shared_EVE,
       aes(x=log(beta),y=-log(pval),color=category))+
  geom_point()+
  scale_color_manual(values=c("#f37225","#00b9e3","black"))+
  geom_vline(xintercept = -0.8864436, col="black", linetype="dashed")+
  annotate("text", x = 0.9, y = 9.5, label = "Shared beta", size = 3)+
  geom_hline(yintercept = -log(0.05), col="black",linetype = "dashed")+
  annotate("text", x = -4, y = 2.5, label = "p = 0.05", size = 3)+
  xlim(-5,5)+
  theme_light()
p + labs(color = "EVE Gene Category",title = "EVE Genes Volcano Plot",subtitle = "37 Lineage-Specific and 112 Highly Variable Genes") +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
```

<div align="center">
  <img src="https://github.com/user-attachments/assets/f872eaaf-e62e-4f9c-909a-d1d8c9f595bc" alt="image" />
</div>


### Highly Variable Genes

```{r,results='hide',tidy=TRUE}
sig_EVE <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Signficant_EVE_results_Bacteria_Orthogroups_BestORF.csv")
colnames(sig_EVE)
sig_EVE <- sig_EVE[,c(2,7,8,11)]   #Orthogroup,beta, pval, category
dim(sig_EVE)   #928 significant EVE genes
names(sig_EVE)
rlogs <- read.csv("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/treatment_rlog_SingleBestORF.csv")
names(rlogs)
colnames(rlogs)[1] <- "Orthogroup"
sig_EVE_rlogs <- left_join(sig_EVE,rlogs,by="Orthogroup")    #39
dim(sig_EVE_rlogs)
names(sig_EVE_rlogs)
write.csv(sig_EVE_rlogs,file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/Signficant_EVE_results_Bacteria_Orthogroups_BestORF_rlogs.csv")

plasticGenes <- subset(sig_EVE_rlogs, sig_EVE$category == "Highly Variable")
names(plasticGenes)
plasticGenes <- plasticGenes[,-c(2:4)]
dim(plasticGenes)
# [1] 30 77
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
# [1] 70 47
write.csv(HV_metadata,file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/highlyVariable_rlogs_orthogroup_SingleBestORF.csv")
```


### HV PCA

```{r}
# Load required libraries
library(ggplot2)
library(dplyr)

# Assume HV_metadata is already loaded and cleaned
# Select only numeric columns for PCA (remove metadata columns)
numeric_data <- HV_metadata[,18:length(HV_metadata)]

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
  scale_color_brewer(palette = "Set2") +  # Use a nice color palette for species
  xlim(-10, 10) +  # Set x-axis limits
  ylim(-15, 12)      # Set y-axis limits # Use a nice color palette for species

```

<div align="center">
  <img src="https://github.com/user-attachments/assets/b4b0cc20-94ac-485e-907b-82b79677ab74" alt="image" />
</div>


### Lineage-Specific Genes

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
# [1]  70 294
# [1] 70 26
write.csv(LS_metadata,file="/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Coral Host/LineageSpecific_rlogs_orthogroup_LongestORF.csv")
```

### LS PCA

```{r}
# Load required libraries
library(ggplot2)
library(dplyr)

# Assume HV_metadata is already loaded and cleaned
# Select only numeric columns for PCA (remove metadata columns)
numeric_data <- LS_metadata[,18:length(LS_metadata)]

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

<div align="center">
  <img src="https://github.com/user-attachments/assets/5f5b73ea-1192-442e-b623-b6f263952288" alt="image" />
</div>


### EVE - Combined Plot

```{r,fig.height=6,fig.width=3,dpi=300}
# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggdendro)
library(dendextend)

# Create the hierarchical tree plot
numeric_data <- LS_metadata[,18:length(LS_metadata)]
numeric_data_scaled <- scale(numeric_data)
distance_matrix <- dist(numeric_data_scaled, method = "euclidean")
hclust_results <- hclust(distance_matrix, method = "average")
dend <- as.dendrogram(hclust_results)

# Color branches by Coral.Species
species_colors <- as.factor(LS_metadata$Coral.Species)
dend <- color_branches(dend, k = length(unique(species_colors)), groupLabels = FALSE)

# Label the dendrogram with Species_Outcome
labels(dend) <- LS_metadata$Species_Outcome

# Convert to a plot-ready object using ggdendro
ggdend <- ggdendro::dendro_data(dend)

# Plot the dendrogram using ggplot2
tree_plot <- ggplot(ggdend$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = ggdend$labels, aes(x = x, y = y, label = label, color = as.factor(species_colors)), hjust = 1, size = 3) +
  labs(title = "Dendrogram - Average Linkage", y = "Height") +
  geom_text(data = ggdend$labels, aes(x = x, y = y, label = label, color = as.factor(species_colors)),
            hjust = 1, size = 3, angle = 30) +  # This line adds angled labels
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Set2")

# Display the plot
print(tree_plot)


# LS PCA plot
pca_results <- prcomp(numeric_data, scale. = TRUE)
percent_variance <- round(100 * (pca_results$sdev^2 / sum(pca_results$sdev^2)), 1)
pca_data <- as.data.frame(pca_results$x)
pca_data$Coral.Species <- LS_metadata$Coral.Species
LS_PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Coral.Species)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(
    title = "Lineage-Specific Genes", 
    x = paste0("PC1: ", percent_variance[1], "% variance"), 
    y = paste0("PC2: ", percent_variance[2], "% variance")
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

# Volcano plot
volcano_plot <- ggplot(data = Shared_EVE, aes(x=log(beta), y=-log(pval), color=category)) +
  geom_point() +
  scale_color_manual(values=c("#f37225","#00b9e3","black")) +
  geom_vline(xintercept = -0.8864436, col="black", linetype="dashed") +
  annotate("text", x = 0.9, y = 9.5, label = "Shared beta", size = 3) +
  geom_hline(yintercept = -log(0.05), col="black", linetype = "dashed") +
  annotate("text", x = -4, y = 2.5, label = "p = 0.05", size = 3) +
  xlim(-5, 5) +
  theme_light() +
  labs(color = "EVE Gene Category", title = "EVE Genes Volcano Plot", subtitle = "37 Lineage-Specific and 112 Highly Variable Genes") +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9)) +
  theme(legend.position = "none")  # This line removes the legend


# HV PCA plot
numeric_data_HV <- HV_metadata[,18:length(HV_metadata)]
pca_results_HV <- prcomp(numeric_data_HV, scale. = TRUE)
percent_variance_HV <- round(100 * (pca_results_HV$sdev^2 / sum(pca_results_HV$sdev^2)), 1)
pca_data_HV <- as.data.frame(pca_results_HV$x)
pca_data_HV$Coral.Species <- HV_metadata$Coral.Species
HV_PCA_plot <- ggplot(pca_data_HV, aes(x = PC1, y = PC2, color = Coral.Species)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(
    title = "Highly Variable Genes", 
    x = paste0("PC1: ", percent_variance_HV[1], "% variance"), 
    y = paste0("PC2: ", percent_variance_HV[2], "% variance")
  ) +
  theme_minimal() +
  theme(legend.position = "none")+
  scale_color_brewer(palette = "Set2")

# Arrange the plots using gridExtra
combined_plot <- grid.arrange(
  tree_plot, 
  arrangeGrob(LS_PCA_plot, volcano_plot, HV_PCA_plot, ncol = 3), 
  heights = c(1, 2),
  ncol = 1
)

combined_plot <- grid.arrange(
 # tree_plot, 
  arrangeGrob(volcano_plot, LS_PCA_plot, HV_PCA_plot, nrow = 3), 
  #heights = c(1, 2),
  ncol = 1
)


# Save the combined figure (optional)
#ggsave("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/EVE_Summary_BacteriaOrthogroups_Column_BestORF.png", combined_plot, width = 3, height = 8)

#ggsave("/Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/EVE_Summary_BacteriaOrthogroups_Column_BestORF.pdf", combined_plot, width = 3, height = 8)

```

<div align="center">
  <img src="https://github.com/user-attachments/assets/50259de4-e694-4f33-bc22-ee7ef1c5a315" alt="image" />
</div>


