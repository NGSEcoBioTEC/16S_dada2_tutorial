# https://benjjneb.github.io/dada2/tutorial.html
# https://github.com/benjjneb/dada2/issues/311
# https://stacks.stanford.edu/file/druid:mh194vj6733/MouseFeces_Final.html
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# Load library
library(dada2)

set.seed(57)

# Which directory
getwd()
# [1] "/Volumes/Seagate Backup Plus Drive/SOI 2018/16S"

# Define the path
path <- "/Volumes/Seagate Backup Plus Drive/SOI 2018/16S/Biocant_files/raw_fastq" #these files are as raw as we can get
list.files(path)

# read the names in the fastaq files
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

# Plot quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# gray-scale - heat map of the frequency of each quality score at each base position
# green line - median quality score at each position
# orange lines - quartiles of the quality score distribution
# red line - scaled proportion of reads that extend to at least that position

# Filter and trim 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)


# standard filtering parameters: maxN=0, maxEE=2, truncQ=2
# maxEE can be relaxed to allow more expected errors - see https://academic.oup.com/bioinformatics/article/31/21/3476/194979
# truncQ - leave at its default of truncQ=2 and let the maxEE parameter be the primary quality filter. We've found expected errors to generally be the best quality filtering method -https://github.com/benjjneb/dada2/issues/232

head(out) # number of reads passing seems ok, it did not eliminate most of the reads

# Look at the filtered resullts
pathFilt <- "/Volumes/Seagate Backup Plus Drive/SOI 2018/16S/Biocant_files/raw_fastq/filtered" 
list.files(pathFilt)

fnFsFilt <- sort(list.files(pathFilt, pattern="_F_filt.fastq.gz", full.names = TRUE))
fnRsFilt <- sort(list.files(pathFilt, pattern="_R_filt.fastq.gz", full.names = TRUE))

plotQualityProfile(fnFsFilt[1:2])
plotQualityProfile(fnRsFilt[11:12])

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# points - the observed error rates for each consensus quality score
# black line - the estimated error rates
# red line - the error rates expected under the nominal definition of the Q-score

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# inspecte the returned dada-class object
dadaFs[[1]] # - dada returns the result for each individual sample, in this case we see sample 1 
dadaRs[[1]] 

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

head(mergers[[1]])

# Construct sequence table - We can now construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
## [1] 38 1115

# The sequence table is a matrix with rows corresponding to (and named by) the samples
# and columns corresponding to (and named by) the sequence variants

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## !! the lenght is not the expected given the truncLen applied

# Track reads through the pipeline - look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab")
rownames(track) <- sample.names
track

# outside of filtering there should be no step in which a majority of reads are lost
# have at least thousands of sequences per-sample, to have good resolution down to 1% frequency - https://github.com/benjjneb/dada2/issues/232

saveRDS(seqtab, "seqtab")

# Remove chimeras - Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 3101 bimeras out of 12686 input sequences.

dim(seqtab.nochim)
## [1]  38 9585
sum(seqtab.nochim)/sum(seqtab)
## [1] 0.9391222 - the frequency of non-chimeric sequences (~94%)

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim),reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","perc_reads")
rownames(track) <- sample.names
track

# keep a good proportion of my reads (i.e. >>10%) all the way through the pipeline
# sample 5418  !! is 9.6% retained seqs, different from everything else
# use mergePairs(..., just Concatenate=TRUE) - just for this sample?

# Create an ASV table
write.csv(seqtab.nochim, "asv_soi16S.csv")

# Create fasta table
library(seqinr)

name.seq <- paste0("ASV", seq(ncol(seqtab.nochim)))
uniqueSeqs <- as.list(colnames(seqtab.nochim))
write.fasta(uniqueSeqs, name.seq, "/Users/MBaptista/Desktop/SOI16/uniqueSeqs_SOI16S.fasta")

######################### Assign taxonomy #############################

# The paper introducing the IDTAXA algorithm reports classification performance that is better than the long-time standard set by the naive Bayesian classifier.
# Here we include a code block that allows you to use IdTaxa as a drop-in replacement for assignTaxonomy

#http://www2.decipher.codes/Classification/TrainingSets/

library(DECIPHER)

load("/Users/MBaptista/Desktop/SOI16/RDP_v16-mod_March2018.RData")
#load("/Users/MBaptista/Desktop/SOI16/GTDB_r86-mod_September2018.RData") #- choose accordingly

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

# Create a taxa table
write.csv(taxid, "taxa_soi16S_GTDB.csv")

######## Handoff to phyloseq

library(phyloseq)

ps_soi16_rdp <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxid))

save(ps_soi16_rdp, file = "ps_soi16_rdp.rds")

ps_soi16_gtdb <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                     tax_table(taxid))

save(ps_soi16_gtdb, file = "ps_soi16_gtdb.rds")

######################### Tree #############################

library(phangorn)

### Align sequences
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)

compact.align <- unique(phangAlign)
length(compact.align) # 9585
length(phangAlign) # 9585

#Use only NNI rearrangements first
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))

#This is only 5 times instead of (a maximum) of 100 rearrangements
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic", ratchet.par = list(iter = 5L, maxit = 5L, prop = 1/3), multicore = TRUE)

# and in the end optimize all the other parameters again (they should not change much)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))

save(fitGTR, file = "fitGTR_soi16.rds")

###################

ps_soi16_rdp_tree <- phyloseq(otu_table(ps_soi16_rdp), 
                    tax_table(ps_soi16_rdp),
                    phy_tree(fitGTR$tree))

save(ps_soi16_rdp_tree, file = "ps_soi16_rdp_tree.rds")

ps_soi16_gtdb_tree <- phyloseq(otu_table(ps_soi16_gtdb), 
                              tax_table(ps_soi16_gtdb),
                              phy_tree(fitGTR$tree))

save(ps_soi16_gtdb_tree, file = "ps_soi16_gtdb_tree.rds")
