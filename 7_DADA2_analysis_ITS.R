## DADA2 ITS pipeline via Bioconductor Pipeline 
## for the manuscript Timing matters: Land use determines temporal dynamics in structure and function of fungal communities in streams
## written by R. Salis, Department of Biology and Environmental Science, Linnaeus University, Sweden


## Install DADA2 ##
install.packages("Rcpp")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("dada2")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")


## Read files ##
path = "/raw" #path to raw fastq files
list.files(path)

# extract the file names as a vector
f.names = as.vector(list.files(path, pattern = "_1.fastq.gz",full.names = F))
r.names = as.vector(list.files(path, pattern = "_2.fastq.gz",full.names = F))

fnFs = sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Identify primers
FWD <- "GTGARTCATCGAATCTTTG"  ## fITS7 forward primer sequence 19 bp
REV <- "TCCTCCGCTTATTGATATGC"  ## ITS4 reverse primer sequence 20 bp

#verify the presence and orientation of the primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# “pre-filter” the sequences just to remove those with ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#Remove Primers - with cutadapt
cutadapt <- "/miniconda3/bin/cutadapt" # the cutadapt path
system2(cutadapt, args = "--version") 

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, 
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], 
                             fnFs.filtN[i], fnRs.filtN[i])) 
}

# check presence of primers in cutadapt-ed samples 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#read in the names of the cutadapt-ed FASTQ files and get the matched lists of forward and reverse fastq files.
# Forward and reverse fastq filenames have the format:
#path.cut = "/Volumes/mana_seagate/WP2_ITS/raw/cutadapt"
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))

# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

## Plot quality scores from .fastq files ##
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

## Filter and trim (quality control) ##
#Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#Filtering paraments: maxN=0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE=2 (this means a maximum number of “expected errors” is allowed in a read)
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
saveRDS(out, "out.rds")

# Learn forward error rates
errF <- learnErrors(filtFs, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, multithread=TRUE)
#visualise
plotErrors(errF, nominalQ = TRUE)

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Generate sequence table ##
seqtab = makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, "seqtab.rds")
table(nchar(getSequences(seqtab)))

## Remove chimeras ##
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

#Inspect distribution of sequence lengths:
table(nchar(getSequences(seqtab.nochim)))

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
saveRDS(track, "track.rds")

#################################################################
## Assign taxonomy
unite.ref <- "sh_general_release_dynamic_s_02.02.2019.fasta"

taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
saveRDS(taxa, "taxa.rds")

#view the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#################################################################
## Data analysis

## ITS analysis via Phyloseq and Ampvis

##Construct Phyloseq Object 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("phyloseq")
BiocManager::install("S4Vectors")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(S4Vectors); packageVersion("S4Vectors")
library(DESeq2); packageVersion("DESeq2")



theme_set(theme_bw())

taxa <- readRDS("taxa.rds")
seqtab.nochim <- readRDS("seqtab.nochim.rds")
DFG.ITS.2 = phyloseq(tax_table(taxa), otu_table(seqtab.nochim, taxa_are_rows = FALSE))
DFG.ITS.2
rownames(otu_table(DFG.ITS.2))
rownames(tax_table(DFG.ITS.2))
## Upload metadata ##
sd.DFG.ITS1 = read.csv("8_sample_data.csv")
head(sd.DFG.ITS1)
###Rename samples for metadata
sd.DFG.ITS <- sd.DFG.ITS1[, -1]
row.names(sd.DFG.ITS)<- sd.DFG.ITS1$x #sample data row names must align with dada2 rowname outputs
sd.DFG.ITS = as.data.frame(sd.DFG.ITS)
head(sd.DFG.ITS)
DFG.ITS.ex = phyloseq(tax_table(taxa), otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(sd.DFG.ITS))
DFG.ITS.ex
###Rename sequence variants
a.vec = as.vector(1:2046)  #total ASVs
a.nam = cbind("ASV", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(DFG.ITS.ex) = asv.names$asv.names
taxa_names(DFG.ITS.ex)
DFG.ITS.ex

#subset PCR negative control
DFG.ITS.ex.neg = subset_samples(DFG.ITS.ex, Extraction_ID == "PCR_NEG")
DFG.ITS.ex.neg 
DFG.ITS.ex.neg <- filter_taxa(DFG.ITS.ex.neg, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
DFG.ITS.ex.neg
tax_table(DFG.ITS.ex.neg)
#remove PCR negative control
DFG.ITS.ex.noneg = subset_samples(DFG.ITS.ex, Extraction_ID != "PCR_NEG")
DFG.ITS.ex.noneg 
###taxonomic filtering
rank_names(DFG.ITS.ex.noneg)
#table(tax_table(DFG.ITS.ex.noneg)[, "Kingdom"], exclude = NULL)
#DFG.ITS.ex.noneg
#just filter fungi
#DFG.ITS.ex.nonegF <- subset_taxa(DFG.ITS.ex.noneg, Kingdom %in% "k__Fungi")
#DFG.ITS.ex.nonegF
#DFG.ITS.ex.noneg_other <- subset_taxa(DFG.ITS.ex.noneg, !Kingdom %in% "k__Fungi")
#DFG.ITS.ex.noneg_other
#table(tax_table(DFG.ITS.ex.nonegF)[, "Kingdom"], exclude = NULL)
#table(tax_table(DFG.ITS.ex.noneg_other)[, "Kingdom"], exclude = NULL)
#create a table of read counts for each Phylum present in the dataset
table(tax_table(DFG.ITS.ex.noneg)[, "Phylum"], exclude = NULL)
#remove taxa for which Phylum is NA
DFG.ITS.ex.noneg0 <- subset_taxa(DFG.ITS.ex.noneg, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table(tax_table(DFG.ITS.ex.noneg0)[, "Phylum"], exclude = NULL)
DFG.ITS.ex.noneg0 <- subset_taxa(DFG.ITS.ex.noneg0, (Phylum!="p__unidentified") | is.na(Phylum))
DFG.ITS.ex.noneg0 
#remove mock sample 
ITS.nomock  = subset_samples(DFG.ITS.ex.noneg0, Extraction_ID != "Mock")
ITS.nomock <- filter_taxa(ITS.nomock, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
ITS.nomock
# Transform data to proportions 
ITS.prop <- transform_sample_counts(ITS.nomock, function(otu) otu/sum(otu))
saveRDS(ITS.prop, "ITS.prop.rds")

ITS.prop.asv <- data.frame(otu_table(ITS.prop))
ITS.sd <- data.frame(sample_data(ITS.prop))
ITS.tax <- data.frame(tax_table(ITS.prop))
library(dplyr)
ITSdata <- ITS.sd %>% select(ID, Stream, treatment, type, timepoint, Miseq.run)
head(ITSdata)
ITSdataprop <- merge(ITSdata, ITS.prop.asv, by=0)

write.csv(ITSdataprop, file = "5_ITS_DataProp.csv")
write.csv(ITS.tax, file = "9_ITS_Taxonomy.csv")

## annotate Aquatic Hyphomycetes

# List of genera to check
AH_genera_list <- c("g__Alatospora", "g__Anguillospora", "g__Articulospora", "g__Campylospora",
                 "g__Clavatospora", "g__Cudoniella", "g__Culicidospora", "g__Dendrospora",
                 "g__Dimorphospora", "g__Filosporella", "g__Flagellospora", "g__Gyoerffyella",
                 "g__Hydrocina", "g__Lemonniera", "g__Lunulospora", "g__Mycoarthris",
                 "g__Mycocentrospora", "g__Mycofalcella", "g__Tetrachaetum", "g__Tetracladium",
                 "g__Tricladium", "g__Triscelophorus", "g__Tumularia", "g__Varicosporium",
                 "g__Hymenoscyphus s__tetracladius", "g__Hymenoscyphus s__varicosporoides")

# Add the AQH column
ITS.tax_AH <- ITS.tax %>%
  mutate(AQH = ifelse(Genus %in% AH_genera_list | 
                        (Genus == "g__Hymenoscyphus" & Species %in% c("s__tetracladius", "s__varicosporoides")), "Y", ""))
write.csv(ITS.tax_AH, file = "6_ITS_Taxonomy_AH.csv")

