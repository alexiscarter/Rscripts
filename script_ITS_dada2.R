# DADA2 for ITS Illumina sequences from SBL project
# Alexis Carteron
# 23/11/2017

# Tutorial
# Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), 581-583.
# http://benjjneb.github.io/dada2/tutorial.html

# This workflow assumes that the data you are starting with meets certain criteria:
# -Non-biological nucleotides have been removed (primers/adapters/barcodes…) (can be removed in filterandtrim function)
# -Samples are demultiplexed (split into individual per-sample fastqs)
# -If paired-end sequencing, the forward and reverse fastqs contain reads in matched order

# Run on Intel® Core™ i7-6900K CPU @ 3.20GHz × 16 with memory of 62.8 GiB 

#### Load packages and data ####
library(dada2); packageVersion("dada2")
# ‘1.5.8’

#~-~-~-~-~-~-~-~-~~-~-~-~-~-~-#
###### ITS PAIRED-END #########
#~-~-~-~-~-~-~-~-~-~-~-~-~~-~-#

# rename file in the terminal
# for f in *; do mv "$f" "${f:12}"; done

path <- "/home/udem/Documents/dada2/ITS"
list.files(path)

#### Filter and Trim ####

# Read in the names of the fastq files, and perform some string manipulation to get lists of the forward and reverse fastq files in matched order:
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names 

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#### Examine quality profiles of forward and reverse reads ####

## looking at data, starting by visualizing the quality profiles by sample and then samples combined
plotQualityProfile(fnFs[1:81], n = 5e+05) # forward reads
plotQualityProfile(fnRs[1:81], n = 5e+05) # reverse reads

plotQualityProfile(fnFs, n = 5e+06, aggregate = TRUE) 
plotQualityProfile(fnRs, n = 5e+06, aggregate = TRUE) 
# In the figures: the mean is in green, the median the solid orange line and the quartiles are the dotted orange lines.

# Phred score of 10 => 90 % precision
# Phred score of 20 => 99 % precision
# Phred score of 30 => 99.9 % precision

#### Perform filtering and trimming ####

filt_path <- file.path(path, "filtered_pairedend") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

## Standard filtering parameters
# maxN=0 (DADA2 requires no Ns), 
# truncQ=2, 
# rm.phix=TRUE
# maxEE=2, sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
# /!\ For common ITS amplicon strategies, it is undesirable to truncate reads to a fixed length due to the large amount of length variation at that locus.

## ITS primers infos
# ITS3_KYO2	GATGAAGAACGYAGYRAA = 18bp (forward)
# ITS4	TCCTCCGCTTATTGATATGC = 20bp (reverse)


# Filter the forward and reverse reads:
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     minLen = 50,
                     truncLen=c(290,290),
                     maxN=0, 
                     maxEE=c(2,2), 
                     truncQ=6, # no difference up to 6 (probably mean it is already filtered). When = 7 only 0.2707456 % kept
                     rm.phix=TRUE,
                     trimLeft=c(18,20),
                     compress=TRUE, 
                     multithread=TRUE)

# Advice: In addition trim the first 10 nucleotides of each read based on empirical observations 
# Across many Illumina datasets that these base positions are particularly likely to contain pathological errors. (Bioconductor workflow)

# Check pourcentage of discarded reads
pourc <- cbind(out[,2]/out[,1])
plot(out)
plot(pourc)
pourc_disc <- cbind(out, pourc)
pourc_disc 

# overall mean of discarded reads
mean(out[,2])/mean(out[,1])
# 0.4361927

#### Learn the Error Rates ####
# It learns the error model from the data, by alternating estimation of the error rates and inference of sample
errF <- learnErrors(filtFs, nreads = 1e+07, # using all the reads. 
                    multithread=TRUE)
#Convergence after  6  rounds.
#Total reads used:  2590495
#total reads doesn't change if re-reun
#If nreads = 1e+06, convergence is after 5 rounds and Total reads used is 1013240. PlotErros output very very similar.

errR <- learnErrors(filtRs, nreads = 1e+07, 
                    multithread=TRUE)
#Convergence after  5  rounds.
#Total reads used:  2590495

# Visualize the estimated error rates (sanity check)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# The error rates for each possible transition (eg. A->C, A->G, …) are shown. 
# Points are the observed error rates for each consensus quality score. 
# The black line shows the estimated error rates after convergence. 
# The red line shows the error rates expected under the nominal definition of the Q-value (for Illumina technology?). 
# /!\ If the plotted error model does not look like a good fit, try increasing nreads() (1M by default)

#### Dereplication ####

#combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”
#reduces computation time by eliminating redundant comparisons
derepFs <- derepFastq(filtFs, 
                      n = 1e+06) # no change obseved when n = 1e+07

derepRs <- derepFastq(filtRs, 
                      n = 1e+06)
# If dataset exceeds available RAM, it is preferable to process samples one-by-one in a streaming fashion
# The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads.

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### Sample Inference ####
# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, 
               err=errF, 
               BAND_SIZE = 16, # Default is 16, based on 16S amplicon data that has lower indel rates than ITS. Only cost: run-time is longer
               # no differences observed when = 32
               KDIST_CUTOFF = 0.42, # Default is 0.42, based on 16S amplicon data
               # no differences observed when = .6
               multithread=TRUE)

dadaFs.pool <- dada(derepFs, err=errF, 
                    BAND_SIZE = 16, # Default is 16, based on 16S amplicon data that has lower indel rates than ITS. Only cost: run-time is longer
                    KDIST_CUTOFF = 0.42, # Default is 0.42, based on 16S amplicon data 
                    multithread=TRUE, 
                    pool=TRUE)

#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, (sapply(dadaFs, getN) - sapply(dadaFs.pool, getN)))

dadaRs <- dada(derepRs, 
               err=errR, 
               multithread=TRUE)
dadaFs[[1]]
# 128 sample sequences were inferred from 20805 input unique sequences.
# for other infos check help("dada-class")
dadaRs[[1]]
# 124 sample sequences were inferred from 19404 input unique sequences.

#### Merging ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      minOverlap = 20, 
                      maxMismatch = 0, 
                      returnRejects = FALSE, 
                      propagateCol = character(0),
                      justConcatenate = FALSE, 
                      trimOverhang = FALSE)
# That way paired reads that did not exactly overlap were removed

# if returnRejects = TRUE
# pairs that that were rejected based on mismatches in the overlap region are retained BUT with accept = FALSE in table

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
tail(mergers[[1]])

#### Construct sequence table ####

# SV version of the OTU table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 81 4161
# 81 3845 the 17/10/2017

# Inspect distribution of sequence lengths
#seqtab.length.var <- table(nchar(getSequences(seqtab)))
# huge variance in length, from 272 to 512 (as expected ?)
# write.table(length.var, file="dada2/saved_table/seqtab.length.var.ITS.paired.txt", row.names=TRUE, col.names=TRUE)
pdf("dada2/saved_table/seqtab.length.var.plot.pdf")
seqtab.length.var.plot <- plot(table(nchar(getSequences(seqtab))))
dev.off()

#### Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method="pooled", # The samples in the sequence table are all pooled together for bimera identification. Default = consensus
                                    multithread=TRUE, verbose=TRUE) 
# Identified 623 bimeras out of 3845 input sequences
# method="pooled" is longer but identified 623 bimeras vs 502 for for method="consensus" (!!!)

dim(seqtab.nochim)
# 81 3222
sum(seqtab.nochim)/sum(seqtab) #  percentage of the total sequence reads
# 0.9692966

# save seq table
# write.table(seqtab.nochim, file="dada2/saved_table/seqtab.nochim.ITS.txt", row.names=TRUE, col.names=TRUE)
# save(seqtab.nochim, file="dada2/saved_table/seqtab.nochim.ITS.rdata")

#### Assign taxonomy ####

taxa.paired <- assignTaxonomy(seqtab.nochim, "reference_database/sh_general_release_dynamic_10.10.2017.fasta", 
                              minBoot = 50, #Default 50. The minimum bootstrap confidence for assigning a taxonomic level.
                              multithread=TRUE, verbose=TRUE)
unname(taxa.paired)
unique(unname(taxa.paired[,7])) # 458 species identified (vs 461 found on a previous run) ###(Minus 1 for NA)###
unique(unname(taxa.paired[,6])) # 353 genus identified
unique(unname(taxa.paired[,5])) # 207 family identified
unique(unname(taxa.paired[,4])) # 93 order identified
unique(unname(taxa.paired[,3])) # 40 class identified
unique(unname(taxa.paired[,2])) # 15 Phylum identified

# save it
# write.table(taxa.paired, file = "dada2/saved_table/assigntaxaITS.paired.txt", row.names=TRUE, col.names=TRUE)

# NB1: dada2 does not throw away singleton reads. 
# However, it does not infer biological sequence variants that are only supported by a single read - singletons are assumed too difficult to differentiate from errors. 
# Hence no singletons in the output table of amplicon sequence variants.

# NB2: The ASVs with no species assignment do not match the same species in over 50% of the bootstrap replicate kmer-based assignments 
# (see Wang et al. 2007 for more info on the naive Bayesian classifier method).

#~-~-~-~-~-~-~-~-~~-~-~-~-~-~-#
####### ITS FWD ONLY ##########
#~-~-~-~-~-~-~-~-~~-~-~-~-~-~-#

#### Construct sequence table ####
# ESV table for forward only version of the OTU table
seqtabFs <- makeSequenceTable(dadaFs)
dim(seqtabFs)
# 81 3956

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtabFs)))
# all 272 long

#### Remove chimeras ####
seqtab.nochimFs <- removeBimeraDenovo(seqtabFs, 
                                      method="pooled", 
                                      multithread=TRUE, 
                                      verbose=TRUE)
# Identified 705 bimeras out of 3956 input sequences.
dim(seqtab.nochimFs)
# 81 3251
sum(seqtab.nochimFs)/sum(seqtabFs) #  percentage of the total sequence reads
# 0.9649712

# save seq table
# write.table(seqtab.nochimFs, file="dada2/saved_table/seqtab.nochimFs.ITS.txt", row.names=TRUE, col.names=TRUE)
# save(seqtab.nochimFs, file="dada2/saved_table/seqtab.nochimFs.ITS.rdata")

#### Assign taxonomy ####
taxa.Fs <- assignTaxonomy(seqtab.nochimFs, "reference_database/sh_general_release_dynamic_10.10.2017.fasta", multithread=TRUE, verbose=TRUE)
unname(taxa.Fs)
unique(unname(taxa.Fs[,7])) # 460 species found
length(unique(unname(taxa.Fs[,6]))) # 364 genus identified
length(unique(unname(taxa.Fs[,5]))) # 209 family identified
length(unique(unname(taxa.Fs[,4]))) # 93 order identified
length(unique(unname(taxa.Fs[,3]))) # 41 class identified
length(unique(unname(taxa.Fs[,2]))) # 15 Phylum identified
# save it
#write.table(taxa.Fs, file = "dada2/saved_table/assigntaxaITS.Fs.txt", row.names=TRUE, col.names=TRUE)

#~-~-~-~-~-~-~-~-~~-~-~-~-~-~-#
####### ITS REV ONLY ##########
#~-~-~-~-~-~-~-~-~~-~-~-~-~-~-#

#### Construct sequence table ####
# ESV table for forward only version of the OTU table
seqtabRs <- makeSequenceTable(dadaRs)
dim(seqtabRs)
# 81 4267

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtabRs)))
# all 270 long

#### Remove chimeras ####
seqtab.nochimRs <- removeBimeraDenovo(seqtabRs, 
                                      method="pooled", 
                                      multithread=TRUE, 
                                      verbose=TRUE)
# Identified 748 bimeras out of 4267 input sequences.
dim(seqtab.nochimRs)
# 81 3519
sum(seqtab.nochimRs)/sum(seqtabRs) #  percentage of the total sequence reads
# 0.9660918

# save seq table
# write.table(seqtab.nochimRs, file="dada2/saved_table/seqtab.nochimRs.ITS.txt", row.names=TRUE, col.names=TRUE)
# save(seqtab.nochimRs, file="dada2/saved_table/seqtab.nochimRs.ITS.rdata")

#### Assign taxonomy ####

taxa.Rs <- assignTaxonomy(seqtab.nochimRs, "reference_database/sh_general_release_dynamic_10.10.2017.fasta", multithread=TRUE, verbose=TRUE)
unname(taxa.Rs)
unique(unname(taxa.Rs[,7])) # 1 species identified 9???)
length(unique(unname(taxa.Rs[,6]))) # 1 genus identified
length(unique(unname(taxa.Rs[,5]))) # 1 family identified
length(unique(unname(taxa.Rs[,4]))) # 1 order identified
length(unique(unname(taxa.Rs[,3]))) # 2 class identified
length(unique(unname(taxa.Rs[,2]))) # 2 Phylum identified

# save it
write.table(taxa.Rs, file = "dada2/saved_table/assigntaxaITS.Rs.txt", row.names=TRUE, col.names=TRUE)

#### Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim), 
               rowSums(seqtab.nochimFs),
               rowSums(seqtab.nochimRs))

colnames(track) <- c("input", 
                     "filtered", 
                     "denoised F", 
                     "denoised R", 
                     "merged", 
                     "nonchim paired",
                     "nonchim F",
                     "nonchim R")
rownames(track) <- sample.names
head(track)
# save it
# write.table(track, file = "dada2/saved_table/track_ITS_17_10_2017.txt", row.names=TRUE, col.names=TRUE)

#~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-#
#### Comparing Forward only, Reverse only and Paired ####
#~-~-~-~-~-~-~-~-~~-~-~-~-~-~-~-~-~-~-~-~-~-~-~~-~-~-~-~-~-#

# Comparing ESV between Fs and paired but sequences are not the same length...
ESV.Fs <- colnames(seqtab.nochimFs)
ESV.paired <- colnames(seqtab.nochim)
length(Reduce(intersect, list(ESV.Fs,ESV.paired))) 

# from taxa file find how many shared identified taxa between paired reads and forward only and the name of the not shared
# At the species level
sp.Fs <- unique(unname(taxa.Fs[,7])) # 460 species found
sp.paired <- unique(unname(taxa.paired[,7])) # 458 species found
length(Reduce(intersect, list(sp.Fs,sp.paired))) # 415 shared species
setdiff(sp.Fs,sp.paired)

# At the genus level
g.Fs <- unique(unname(taxa.Fs[,6])) # 364 genus identified
g.paired <- unique(unname(taxa.paired[,6])) # 353 genus identified
length(Reduce(intersect, list(g.Fs,g.paired))) # 333 shared genus
setdiff(g.Fs,g.paired)

# At the family level
f.Fs <- unique(unname(taxa.Fs[,5])) # 209 family identified
f.paired <- unique(unname(taxa.paired[,5])) # 207 family identified
nrow(Reduce(intersect, list(f.Fs,f.paired))) # 204 shared genus
setdiff(f.Fs,f.paired)

# Venn diagram
library(eulerr)
v <- euler(c(Forward=460, Paired=458, "Forward&Paired"=415))
plot(v, fill_alpha = 0.4, counts = TRUE, 
     fill = c("darkgreen", "darkred"), border = "transparent",
     auto.key = list(space = "right", columns = 1),
     main = "Shared identified species")

