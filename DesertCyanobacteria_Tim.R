library(plyr)
library(data.table)
library(tidyverse)
library(tidyr)
library(stringr)
library(readxl)
library(readr)
library(ggplot2)
library(markdown)
library(purrr)
library(Hmisc)
library(Formula)
library(latticeExtra)
library(htmlTable)
library(ggeasy)
library(forcats)
library(directlabels)
library(DECIPHER)
library(dplyr)
library(plyr)
library(gplots)
library(ggplot2)
library(gtools)
library(phangorn)
library(dada2); packageVersion("dada2") 
library(ShortRead); packageVersion("ShortRead")
library(plotly)
library(phyloseq)

setwd("/Users/timozcelik/Desktop/Sequencing")
path <- "/Users/timozcelik/Desktop/Sequencing"
project.fp <- path
list.files(path)

# sort fwd and rvs file names, and then extract sample names
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

save.image(file = "desert_cyanobacteria.RData")

preprocess.fp <- file.path(project.fp, "01_preprocess")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter")
table.fp <- file.path(project.fp, "03_tabletax")

fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

pre_out <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0)
pre_out

save.image(file = "desert_cyanobacteria.RData")

cutadapt <- "/Users/timozcelik/mambaforge/envs/cutadaptenv/bin/cutadapt"
system2(cutadapt, args = "--version")

FWD <- "GTGCCAGCMGCCGCGGTAA"
REV <- "GGACTACHVGGGTWTCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer) # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) # Convert back to character vector
}


FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

### Save the primer orientations to pass to cutadapt in fwd and reverse compliment orientations
FWD2 <- FWD.orients[["RevComp"]]
REV2 <- REV.orients[["RevComp"]]

### Save the primer orientations to pass to cutadapt
FWD.orients.2 <- allOrients(FWD2)
REV.orients.2 <- allOrients(REV2)

### sanity check
FWD.orients
REV.orients
FWD.orients.2
REV.orients.2

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN))

rbind(FWD.ForwardReads = sapply(FWD.orients.2, primerHits, fn = fnFs.filtN),
      FWD.ReverseReads = sapply(FWD.orients.2, primerHits, fn = fnRs.filtN),
      REV.ForwardReads = sapply(REV.orients.2, primerHits, fn = fnFs.filtN),
      REV.ReverseReads = sapply(REV.orients.2, primerHits, fn = fnRs.filtN))

### Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

### Save the reverse complements of the primers to variables
FWD.RC.2 <- dada2:::rc(FWD2) #FWD primer
REV.RC.2 <- dada2:::rc(REV2) #Rev primer

save.image(file = "desert_cyanobacteria.RData")

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
R1.flags.2 <- paste("-g", FWD2, "-a", REV.RC.2)
R2.flags.2 <- paste("-G", REV2, "-A", FWD.RC.2)

if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)

fnFs.cut <- file.path(trimmed.fp, basename(fnFs))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

# running cutadapt
for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags,R1.flags.2, R2.flags.2, "-n", 5, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             "--minimum-length", 150, #reads below 100 base pairs are removed
                             "-e", 0, # reduce maximum allowed error rate to 0 instead of default 0.1%
                             "--report", "minimal", # one line summary
                             "-j", 0,	
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut))

rbind(FWD.ForwardReads = sapply(FWD.orients.2, primerHits, fn = fnFs.cut),
      FWD.ReverseReads = sapply(FWD.orients.2, primerHits, fn = fnRs.cut),
      REV.ForwardReads = sapply(REV.orients.2, primerHits, fn = fnFs.cut),
      REV.ReverseReads = sapply(REV.orients.2, primerHits, fn = fnRs.cut))

#### all primers should be zero

save.image(file = "desert_cyanobacteria.RData")

dir.create(filter.fp)

subF.fp <- file.path(filter.fp, "preprocessed_F")
subR.fp <- file.path(filter.fp, "preprocessed_R")
dir.create(subF.fp)
dir.create(subR.fp)

subF.fp
subR.fp

fnFs.Q <- file.path(subF.fp, basename(fnFs))
fnRs.Q <- file.path(subR.fp, basename(fnRs))

#list files, sanity check
file.copy(from = fnFs.cut, to = fnFs.Q)
file.copy(from = fnRs.cut, to = fnRs.Q)

fastqFs <- sort(list.files(subF.fp, pattern="_1.fastq.gz",full.names = TRUE))
fastqRs <- sort(list.files(subR.fp, pattern="_2.fastq.gz",full.names = TRUE))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

fastqFs
fastqRs

rand_samples <- sample(size = 4, 1:length(fastqFs)) # grab samples to plot 

### check chosen samples

rand_samples

fwd.qual.plot <- plotQualityProfile(fastqFs[rand_samples]) + labs(x = "Sequence Position")

rev.qual.plot <- plotQualityProfile(fastqRs[rand_samples]) + labs(x = "Sequence Position")

### save plots

pdf(file="/Users/timozcelik/Desktop/Sequencing/Rout/fwd.qual.plot.pdf", width = 15, height = 8.5)
fwd.qual.plot
dev.off()

pdf(file="/Users/timozcelik/Desktop/Sequencing/Rout/rev.qual.plot.pdf", width = 15, height = 8.5) 
rev.qual.plot 
dev.off()

saveRDS(fwd.qual.plot, "/Users/timozcelik/Desktop/Sequencing/Rout/fwd.qual.plot.rds")
saveRDS(rev.qual.plot, "/Users/timozcelik/Desktop/Sequencing/Rout/rev.qual.plot.rds")

save.image(file = "desert_cyanobacteria.RData")

# Filter and Trim

filtpathF <- file.path(filter.fp, "filtered_F")
filtpathR <- file.path(filter.fp, "filtered_R")
dir.create(filtpathF)
dir.create(filtpathR)

filtpathF.file <- file.path(filtpathF, basename(fnFs))
filtpathR.file <- file.path(filtpathR, basename(fnRs))

### TRIMMING

filt_out <- filterAndTrim(fastqFs, filtpathF.file, fastqRs, filtpathR.file,
                          maxEE=c(2,2), truncQ=2, maxN=0, truncLen=c(200,200),rm.phix=TRUE,
                          compress=FALSE, verbose=TRUE, multithread=24,matchIDs=TRUE)

filt_out

filt_out_summary <- filt_out %>% data.frame() %>% 
  mutate(Samples = rownames(.), percent_kept = 100*(reads.out/reads.in)) %>% 
  select(Samples, everything())

filt_out_summary # summary of samples in filt_out by percentage

filtFs <- sort(list.files(filtpathF, pattern="_1.fastq.gz",full.names = TRUE))
filtRs <- sort(list.files(filtpathR, pattern="_2.fastq.gz",full.names = TRUE))

filtFs
filtRs

# Plot the quality of the filtered fastq files
rand_samples_filt <- sample(size = 4, 1:length(filtFs))
fwd.fil.plot <- plotQualityProfile(filtFs[rand_samples_filt]) + labs(x = "Sequence Position")
rev.fil.plot <- plotQualityProfile(filtRs[rand_samples_filt]) + labs(x = "Sequence Position")

pdf(file="/Users/timozcelik/Desktop/Sequencing/Rout/fwd.fil.plot.pdf", width = 15, height = 8.5)
fwd.fil.plot
dev.off()

pdf(file="/Users/timozcelik/Desktop/Sequencing/Rout/rev.fil.plot.pdf", width = 15, height = 8.5) 
rev.fil.plot
dev.off()

saveRDS(fwd.fil.plot, "/Users/timozcelik/Desktop/Sequencing/Rout/fwd.fil.plot.rds")
saveRDS(rev.fil.plot, "/Users/timozcelik/Desktop/Sequencing/Rout/rev.fil.plot.rds")

save.image(file = "desert_cyanobacteria.RData")

### Learn Error Rates
# check files listing and order
filtFs
filtRs

sample.namesF <- basename(filtFs) # doesn't drop fastq.gz 
sample.namesF

sample.namesF <- gsub("_1.fastq.gz", "", sample.namesF) 
sample.namesF

sample.namesR <- basename(filtRs) #doesn't drop fastq.gz 
sample.namesR

sample.namesR <- gsub("_2.fastq.gz", "", sample.namesR)
sample.namesR

if(!identical(sample.namesF, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR

names(filtFs)
names(filtRs)

loessErrfun_mod3 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } 
    } 
  }
  
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  return(err)
}

errF_3 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

errR_3 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

# Plot errors

pdf(file="/Users/timozcelik/Desktop/Sequencing/Rout/errF_3.pdf", width = 15, height = 8.5)
plotErrors(errF_3, nominalQ=TRUE)
dev.off()

pdf(file="/Users/timozcelik/Desktop/Sequencing/Rout/errR_3.pdf", width = 15, height = 8.5)
plotErrors(errR_3, nominalQ=TRUE)
dev.off()

saveRDS(errF_3, "/Users/timozcelik/Desktop/Sequencing/Rout/errF_3.rds")
saveRDS(errR_3, "/Users/timozcelik/Desktop/Sequencing/Rout/errR_3.rds")

save.image(file = "desert_cyanobacteria.RData")

# Dereplicate reads

derepF <- derepFastq(filtFs, verbose = TRUE)
derepR <- derepFastq(filtRs, verbose = TRUE)

derepF # outputs the summary of dereplication for each file 
derepR # outputs the summary of dereplication for each file

names(derepF) <- sample.namesF 
names(derepR) <- sample.namesR

save.image(file = "desert_cyanobacteria.RData")

# Infer sequence variants, Merge paired-end reads and create sequence table

dadaF.3 <- dada(derepF, err = errF_3, multithread = TRUE, pool = "pseudo")
dadaR.3 <- dada(derepR, err = errR_3, multithread = TRUE, pool = "pseudo")

dadaF.3[[1]]
dadaR.3[[1]]

saveRDS(dadaF.3, "/Users/timozcelik/Desktop/Sequencing/Rout/dadaF.3.rds")
saveRDS(dadaR.3, "/Users/timozcelik/Desktop/Sequencing/Rout/dadaR.3.rds")

mergers.3 <- mergePairs(dadaF.3, derepF, dadaR.3, derepR, verbose = TRUE)

head(mergers.3[[1]])

saveRDS(mergers.3, "/Users/timozcelik/Desktop/Sequencing/Rout/mergers.3.rds")

# onstruct an amplicon sequence variant table (ASV) table

seqtab.3 <- makeSequenceTable(mergers.3)

dim.seqtab.3 <- dim(seqtab.3)
dim.seqtab.3

saveRDS(dim.seqtab.3, "/Users/timozcelik/Desktop/Sequencing/Rout/dim.seqtab.3.rds")

## Inspect distribution of sequence(ASVs) lengths
seqlen.3 <- table(nchar(getSequences(seqtab.3)))
seqlen.3

saveRDS(seqlen.3, "/Users/timozcelik/Desktop/Sequencing/Rout/seqlen.3.rds")

### Histogram of read distribution
hist.3 <- hist(nchar(getSequences(seqtab.3)), main="Distribution of sequence lengths using Option 3 alter loess function (weights only) and enforce monotonicity")

saveRDS(hist.3, "/Users/timozcelik/Desktop/Sequencing/Rout/hist.3.rds")
saveRDS(seqtab.3, "/Users/timozcelik/Desktop/Sequencing/Rout/seqtab.3.rds")

save.image(file = "desert_cyanobacteria.RData")

# Remove Chimera and Assign Taxonomy

SeqTable.all.3  <- readRDS("/Users/timozcelik/Desktop/Sequencing/Rout/seqtab.3.rds")

### Remove chimeras

seqtab.nochim.3 <- removeBimeraDenovo(SeqTable.all.3, method="pooled", multithread= TRUE, verbose = TRUE)

### Save results to disk

saveRDS(seqtab.nochim.3, "/Users/timozcelik/Desktop/Sequencing/Rout/seqtab.nochim.3.rds")

### Check chimera table & Print percentage of our seqences that were not chimeric.

## read abundances of non-chimeric variants
total.reads.3 <- sum(seqtab.3)
total.nonchimeric.reads.3 <- sum(seqtab.nochim.3)

saveRDS(total.reads.3, "/Users/timozcelik/Desktop/Sequencing/Rout/total.reads.3.rds")
saveRDS(total.nonchimeric.reads.3, "/Users/timozcelik/Desktop/Sequencing/Rout/total.nonchimeric.reads.3.rds")

print( paste("Total reads 3: ", total.reads.3))
print( paste("Total non-chimeric reads 3: ", total.nonchimeric.reads.3))

nonchimera.3 <- 100*sum(seqtab.nochim.3)/sum(seqtab.3)

saveRDS(nonchimera.3, "/Users/timozcelik/Desktop/Sequencing/Rout/nonchimera.3.rds")

print( paste("% Abundance of non-chimeric reads 3: ", nonchimera.3))


## gives total number of non-chimeric ASVs
dim.seqtab.nochime.3 <- dim(seqtab.nochim.3)
dim.seqtab.nochime.3

saveRDS(dim.seqtab.nochime.3, "/Users/timozcelik/Desktop/Sequencing/Rout/dim.seqtab.nochime.3.rds")

print( paste("Number of samples & non-chimeric varients 3: ", dim.seqtab.nochime.3 ))

### TRACK READS THROUGH THE PIPELINE - DADA2 Tutorial: https://benjjneb.github.io/dada2/tutorial_1_8.html

getN <- function(x) sum(getUniques(x))

### filt_out returns 2 columnns: reads_in and reads_our

### track reads using  Traditional DADA2 error model
track.3 <- cbind(pre_out, filt_out, sapply(dadaF.3, getN), sapply(dadaR.3, getN), sapply(mergers.3, getN), rowSums(seqtab.nochim.3))
colnames(track.3) <- c("pre.DADA2.input", "pre.DADA2.filtered", "DADA2.input", "DADA2.filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.3) <- sample.names
track.3

### Save results to disk

saveRDS(track.3, "/Users/timozcelik/Desktop/Sequencing/Rout/track.3.rds")

### Assign Taxonomy using naive Bayesian classifier method 

###  After removing chimeras, use a taxonomy database to train a classifer-algorithm to assign names to our sequence variants

### Download DADA2 formated database from https://benjjneb.github.io/dada2/training.html: 
#   NOTE: As of Silva version 138, the official DADA2-formatted reference fastas are optimized for classification of Bacteria and Archaea, 
#   and are not suitable for classifying Eukaryotes.

### Assign taxonomy with Silva db v138 using silva_nr99_v138.1_train_set.fa : with taxonomy till genus level (tutorial) and
#   add species assignment with silva_species_assignment_v128.fa.gz

### Assign taxonomy at genus and species level for reads obtained from Traditional DADA2 error model
taxa.nb.3 <- assignTaxonomy(seqtab.nochim.3, "/Users/timozcelik/Desktop/Sequencing/tax/silva_nr99_v138.1_train_set.fa.gz", tryRC = TRUE,  multithread=TRUE, outputBootstraps = TRUE, verbose = TRUE)
taxa.nb.3.1 <- addSpecies(taxa.nb.3[[1]], "/Users/timozcelik/Desktop/Sequencing/tax/silva_species_assignment_v138.1.fa.gz")
head(taxa.nb.3.1)

saveRDS(taxa.nb.3.1, "/Users/timozcelik/Desktop/Sequencing/Rout/taxa.nb.3.1.rds")

saveRDS(taxa.nb.3[[2]], "/Users/timozcelik/Desktop/Sequencing/Rout/taxa.nb.3.2.rds")

###Genus and Species assignment only 

gs.nb.3 <- assignSpecies(seqtab.nochim.3, "/Users/timozcelik/Desktop/Sequencing/tax/silva_species_assignment_v138.1.fa.gz", tryRC = TRUE, allowMultiple=TRUE)
head(gs.nb.3)

saveRDS(gs.nb.3, "/Users/timozcelik/Desktop/Sequencing/Rout/gs.nb.3.rds")


################################################################################################################################################

### Check Taxonomy

################################################################################################################################################

### Assign taxonomy at genus and species level for reads obtained from Traditional DADA2 error model
taxa.nb.3.print <- taxa.nb.3.1
rownames(taxa.nb.3.print) <- NULL
head(taxa.nb.3.print)

saveRDS(taxa.nb.3.print, "/Users/timozcelik/Desktop/Sequencing/Rout/taxa.nb.3.print.rds")


### save all created objects in the R environemnt. Loading this file will load all created objects so far
## Change the name of the file accordingly
save.image(file = "desert_cyanobacteria.RData")

###############################################################################################################################

### All DADA2 steps are now completed #####

### Move onto Phyloseq to create a phyloseq object for further data analyses and visualization ####

### To create a phyloseq object you will need:
#   OTU table : already created with DADA2
#   Tax table : already created with DADA2
#   Metadata : sample data table with information about samples such as Sample_ids, experiment type. 
#   The metadata must be created such that each row contains infomration unique to each sample
#   Samples are ordered in rows and every column 


seqtab.nochim.3 <- readRDS("/Users/timozcelik/Desktop/Sequencing/Rout/seqtab.nochim.3.rds")

taxa.nb.3.1 <- readRDS("/Users/timozcelik/Desktop/Sequencing/Rout/taxa.nb.3.1.rds")



# Constructing phylogenetic tree

sequences<-getSequences(seqtab.nochim.3)
names(sequences)<-sequences

alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

# metadata <- read_csv("metadata.xls")
metadata <- as.data.frame(read.csv("metadata.csv", sep = ";", row.names = 1))

ps.3 <- phyloseq(otu_table(seqtab.nochim.3, taxa_are_rows=FALSE), 
                 sample_data(metadata), 
                 tax_table(taxa.nb.3.1),
                 phy_tree(fitGTR$tree))

### Use short names for ASVs instead of using entire DNA sequence. http://benjjneb.github.io/dada2/tutorial.html
#   Store the DNA sequences of our ASVs in the refseq slot of the phyloseq object, and then rename our taxa to a short string. That way, the short new taxa names will appear in tables and plots, and we can still recover the DNA sequences corresponding to each ASV as needed with refseq(ps).
dna.3 <- Biostrings::DNAStringSet(taxa_names(ps.3))
names(dna.3) <- taxa_names(ps.3)
ps.3 <- merge_phyloseq(ps.3, dna.3)
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))
ps.3

### Graphs time!

# remove NA
ps.3.clean <- phyloseq::tax_glom(ps.3, "Genus", NArm = TRUE)

# Absolute abundance bar plots by genus, too messy
plot_bar(ps.3.clean, x="Culture_ID", fill="Genus")

# Isolate cyanobacteria
ps.3.cyanobacteria = subset_taxa(ps.3.clean, Phylum == "Cyanobacteria")

# Cyanobacteria abundance bar plots by genus
plot_bar(ps.3.cyanobacteria, x="Culture_ID", fill="Genus")

# Facet plots by sample
plot_bar(ps.3.cyanobacteria, "Genus", fill="Genus", facet_grid=~Culture_ID)

install.packages("wesanderson")
library(wesanderson)
library(RColorBrewer)
display.brewer.all()

# Calculate relative abundance by phylum, and plot
ps.3.glom.phylum <- tax_glom(ps.3.clean, "Phylum")
ps.3.glom.phylum1 <- transform_sample_counts(ps.3.glom.phylum, function(x) x / sum(x))
ps.3.glom.phylum2 <- merge_samples(ps.3.glom.phylum1, "Culture_ID")
ps.3.glom.phylum3 <- transform_sample_counts(ps.3.glom.phylum2, function(x) x / sum(x))
RelativeAbundance_Phylum <- plot_bar(ps.3.glom.phylum3, fill="Phylum")
RelativeAbundance_Phylum + ggtitle("Relative Abundance, Phyla") + xlab("Culture ID") 


# Relative abundance by genus, very messy
ps.3.glom <- tax_glom(ps.3.clean, "Genus")
ps.3.glom1 <- transform_sample_counts(ps.3.glom, function(x) x / sum(x))
ps.3.glom2 <- merge_samples(ps.3.glom1, "Culture_ID")
ps.3.glom3 <- transform_sample_counts(ps.3.glom2, function(x) x / sum(x))
RelativeAbundance_Genus <- plot_bar(ps.3.glom3, fill="Genus")
RelativeAbundance_Genus + ggtitle("Relative Abundance, Genera") + xlab("Culture ID")

# Relative abundance of cyanobacterial genera
ps.3.cyanobacteria.glom <- tax_glom(ps.3.cyanobacteria, "Genus")
ps.3.cyanobacteria.glom1 <- transform_sample_counts(ps.3.cyanobacteria.glom, function(x) x / sum(x))
ps.3.cyanobacteria.glom2 <- merge_samples(ps.3.cyanobacteria.glom1, "Culture_ID")
ps.3.cyanobacteria.glom3 <- transform_sample_counts(ps.3.cyanobacteria.glom2, function(x) x / sum(x))
RelativeAbundance_Cyanobacteria <- plot_bar(ps.3.cyanobacteria.glom3, fill="Genus"
RelativeAbundance_Cyanobacteria + ggtitle("Relative Abundance, Cyanobacteria") + xlab("Culture ID") 

# ASV table
write.table(ps.3 %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              arrange(OTU) %>% rename(ASV = OTU) %>% 
              select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>%
              spread(Sample, Abundance), 
              file = "ps.relative_abundance.all.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

ASV_abundance <- as.data.frame(ps.3 %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% arrange(OTU) %>% rename(ASV = OTU) %>% select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>% spread(Sample, Abundance))

# Phylotree phyla
plot_tree(ps.3.clean, size="abundance", color="Culture_ID", label.tips="Genus", base.spacing=0.1, ladderize="TRUE")
plot_tree(ps.3.clean, size="abundance", color="Culture_ID", label.tips="Genus", base.spacing=0.05, min.abundance=10000000, ladderize="TRUE")

# Phylotree cyanobacterial genera
Phylotree_Cyanobacteria <- plot_tree(ps.3.cyanobacteria, size="abundance", color="Culture_ID", label.tips="Genus", base.spacing=0.06, plot.margin=0.4, ladderize="TRUE")
Phylotree_Cyanobacteria + ggtitle("Phylogenetic Tree, Cultured Cyanobacteria Strains")

plot_tree(ps.3.cyanobacteria.glom3, size="abundance", color="Culture_ID", label.tips="Genus", base.spacing=0.1, ladderize="TRUE") ## Doesn't work, no colour
