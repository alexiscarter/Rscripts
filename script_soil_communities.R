# Annalysis of fungal community and environment at the SBL permanent plot
# Alexis Carteron 
# Created 23/11/2017

# Package
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(ade4)
library(vegan)
library(cluster)
library(dplyr)
library(data.table)
library(clue)
library(gclus)
library(GGally)
library(labdsv)
library(GGally)
library(RVAideMemoire)
library(ggmap)
library(ggsn)
source('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:coldiss?do=export_code&codeblock=1')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Data import and processing ####

# SBL data
label <- read.csv("data/label_soil_sample.csv", sep = ";")
load("data/soil3.rda")

# Dada2 data
taxa_table <- read.table("dada2/saved_table/assigntaxaITS.Fs.txt")
asv_table <- read.table("dada2/saved_table/seqtab.nochimFs.ITS.txt")

# Discard soil samples
# S40 because M_AS_FS_03_Ah06 replicate of M_AS_FS_03_Ah
# all the Ah horizons: S69, S70, S57, S32, S38
remove_tmp <- c("S40", "S69", "S70", "S57", "S32", "S38")
asv_table <- asv_table[!rownames(asv_table) %in% remove_tmp, ]

# Make a data.frame holding the sample data
soil3$myco <- gsub("mixed", "Mix", soil3$myco)
soil3$sampleID <- gsub("1/2", "Mix", soil3$sampleID)
soil <- soil3 %>%
  left_join(label, by = 'sampleID')
rownames(soil) <- soil$sample

# Standardize environmental data
keep <- c("totalC", "totalN", "totalP", "inorgP", "orgP")
soil[,env] <- decostand(soil[,keep], method = "standardize")

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(asv_table, taxa_are_rows=FALSE), 
               sample_data(soil), 
               tax_table(as.matrix(taxa_table)),
               phy_tree(fitGTR$tree))
ps

# Reorder factors
sample_data(ps)$horiz = factor(sample_data(ps)$horiz, levels = c('L', 'F', 'H', 'Ae', 'B'))
sample_data(ps)$myco = factor(sample_data(ps)$myco, levels = c('AM', 'Mix', 'ECM'))

# save ps object
# save(ps, file = "phyloseq/saved/ps.Fs.ITS.rda")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Filtering and transormation ####

# delete ESV with zero abundance
ps = filter_taxa(ps, function(x) sum(x) > 0, TRUE)

# Subset of the dataset with only fungi
ps.fungi = subset_taxa(ps, Kingdom=="k__Fungi")

# relative abundance (fractional abundance)
ps.fungi.rf  = transform_sample_counts(ps.fungi, function(x) x/sum(x))

# ESV with a mean greater than 10e-6 are kept
ps.fungi.rf = filter_taxa(ps.fungi.rf, function(x) mean(x) > 1e-6, TRUE) # why is there ESV with zero abundance ???

ntaxa(ps.fungi.rf)/ntaxa(ps.fungi)*100 # 99.01861 % of ASVs kept

# Hellinger transformation
ps.fungi.hel = transform_sample_counts(ps.fungi.rf, function(x) sqrt(x)) # why is there ESV with zero abundance ???

ntaxa(ps.fungi.hel) # 2926
# Agglomerate taxa
# at the species level
ps.fungi.hel = tax_glom(ps.fungi.hel, "Species", NArm = FALSE)
ntaxa(ps.fungi.species.hel) # 775

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Map of the locations of the sites ####

# Geographic coordinates x and y from the plots data frame
load("data/plots.rda")
plots <- plots[1:15,]

plots$plot <- gsub("1/2", "Mix", plots$plot)
plots$plot <- gsub("FS_", "", plots$plot)
plots$plot <- gsub("M_", "", plots$plot)

#correct longitudes
plots$long <- plots$long * -1

# get the map
boite <- make_bbox(lon=-73.995, lat=45.989)
SBL <- get_map(location=boite, maptype="satellite", source="google", zoom=14)
sbl <- ggmap(SBL)

# add points
sbl2 <- sbl + 
  geom_point(data = plots, aes(long, lat, shape = myco), size = 2, color = 'lightgrey') +
  xlab("Longitude (°)") +
  ylab("Latitude (°)") +
  geom_label_repel(data = plots,
                   aes(long, lat, fill = block, label = plot),
                   fontface = 'bold', color = 'lightgrey',
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = 'lightgrey') +
  scale_shape_discrete((name="Mycorrhizal\nType"), labels=c("AM", "EcM", "Mixed")) +
  scale_fill_discrete(name="Block") +
  scalebar(location="topright", y.min=45.940, y.max=45.973, 
           x.min=-73.98, x.max=-73.97, dist=.5, dd2km= TRUE, model='WGS84',
           st.dist=.04)
north2(sbl2, x=.22, y=.7, symbol=4, scale = .09) # north arrow



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
# Make a table of the chemical characteristics with the mean and standard deviaiton per block

# Make a data.frame holding the sample data
soil3$myco <- gsub("mixed", "Mix", soil3$myco)
soil3$sampleID <- gsub("1/2", "Mix", soil3$sampleID)
soil3$plot <- gsub("1/2", "Mix", soil3$plot)
soil <- soil3 %>%
  left_join(label, by = 'sampleID')
rownames(soil) <- soil$sample

soil$horiz = factor(soil$horiz, levels = c('L', 'F', 'H', 'Ae', 'B'))

# Summarize chemical variables by myco*horiz in a table
soil <- soil[,c("sp", "horizon", "totalC", "totalN", "totalP", "inorgP", "orgP")]

env.chem.mean <- soil %>%
  group_by(sp, horizon) %>%
  summarise_all(funs(mean, sd))

env.mean.sd <- env.chem.mean[,c("sp", "horizon", "totalC_mean", "totalC_sd", "totalN_mean", "totalN_sd", "totalP_mean",  "totalP_sd", "inorgP_mean", "inorgP_sd", "orgP_mean","orgP_sd")]
soil$horiz = factor(soil$horiz, levels = c('L', 'F', 'H', 'Ae', 'B'))
# write.csv(env.mean.sd, file="exploration_R/saved/env.mean.sd.csv")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Exploration ####

# Sample summary
# mean, max and min of sample read counts and total sum
min(sample_sums(ps))
# 5262
mean(sample_sums(ps))
# 31034.68
max(sample_sums(ps))
# 55562
sum(sample_sums(ps))
# 2296566

# getting various info on the phyloseq object
nsamples(ps)
sample_names(ps)
rank_names(ps)
sample_variables(ps)
otu_table(ps)[1:2, 1:2]

# Sequencing depth per samples
plot_bar(ps, x = "ID") +
  ggtitle("Distribution sequencing depth per sample")

plot_bar(ps.fungi.hel, x = "ID") +
  ggtitle("Distribution sequencing depth per sample after Hellinger transformation and grouped at the species level")

# Reads per kingdom 
plot_bar(ps, fill = "Kingdom", x = "ID")

# Check dataset with fungi only
round(sum(sample_sums(ps.fungi)) / sum(sample_sums(ps)) * 100 ,1) # 91.7% of reads kept
ntaxa(ps) # 3251
ntaxa(ps.fungi) # 3075
ntaxa(ps.fungi.hel) # 775

# reads per phylum
plot_bar(ps.fungi, fill = "Phylum") + facet_wrap(~horiz+myco, scales="free_x", nrow = 5, ncol = 3)

# Relative abundance of known taxa at the species level per sample
plot_bar(ps.fungi.hel, fill = "Species", x = "ID") +
  theme(legend.position="none")

# plot ESV-abundance curve
# this converts taxa counts in each sample to a percentage
#phyloTemp = transform_sample_counts(ps.fungi, function(x) 1e+02 * x/sum(x))
clusterData = psmelt(ps.fungi.hel)
clusterData = filter(clusterData,Abundance > 0)
# this is where the mean is calculated and the taxa to display is chosen
clusterAgg = aggregate(Abundance ~ OTU + Phylum,data=clusterData,sum)
# filtering and picking the number to display
clusterAgg = clusterAgg[order(-clusterAgg$Abundance),][1:100,]
ggplot(clusterAgg,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point(aes(color=Phylum),size=3) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10()

# Values for most abundant species
clusterAggSpecies = aggregate(Abundance ~ OTU + Genus + Species,data=clusterData,sum) # sum or mean
clusterAggSpecies[order(-clusterAggSpecies$Abundance),][1:50,]

# Plot observed richness (for exploration only)
plot_richness(ps.fungi, measures="Observed") + 
  ylab("Observed richness") +
  facet_wrap(~horiz+myco, scales="free_x", nrow = 5, ncol = 3)# S8, S19 very low richness

plot_richness(ps.fungi.hel, measures="Shannon") + 
  ylab("Shannon index") +
  facet_wrap(~horiz+myco, scales="free_x", nrow = 5, ncol = 3) # S8, S19 very low richness

# Plot with most abundance species per samples
t = 20
top <- names(sort(taxa_sums(ps.fungi.hel), decreasing=TRUE))[1:t]
ps.top <- prune_taxa(top, ps.fungi.hel)
plot_bar(ps.top, fill = "Species") + facet_wrap(~horiz+myco, scales="free_x", nrow = 5, ncol = 3)
# ggsave("phyloseq/saved/SpeciesTop20.Fs.ITS.fungi.pdf", width = 15, height = 10)

# visualize networks
# among samples
# Hellinger data and Bray-Curtis distance
ps.nw.sample <- make_network(ps.fungi.hel, type = "samples", max.dist = .7, keep.isolates = T, distance = "bray")
plot_network(ps.nw.sample, ps.fungi.hel, color = "horiz", label = "sample", shape = "myco")

# Bray-Curtis distance
ps.nw.sample <- make_network(ps.fungi.rf, type = "samples", max.dist = .9, keep.isolates = T, distance = "bray")
plot_network(ps.nw.sample, ps.fungi.rf, color = "horiz", label = "sample", shape = "myco")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Clustering ####

#ESV data
# ESV with a mean greater than 0 are kept
ps.fungi.f = filter_taxa(ps.fungi, function(x) mean(x) > 0, TRUE) # why is there ESV with zero abundance ???
ntaxa(ps.fungi.f)/ntaxa(ps.fungi)*100 # 96.09756 % of ESVs kept

# Agglomerate taxa at the species level
ps.fungi.species.f = tax_glom(ps.fungi.f, "Species", NArm = FALSE)
ntaxa(ps.fungi.species.f) # 780

# Change ESV DNA sequence for just SV
taxa_names(ps.fungi.species.f) <- paste0("SV", seq(ntaxa(ps.fungi.species.f)))

# Extract taxonomy table with SV label
SV_taxo <- as.data.frame(tax_table(ps.fungi.species.f))

# Extract abundance matrix from the phyloseq object
seqtab1 = as(otu_table(ps.fungi.species.f), "matrix")

# Remove samples
remove_tmp <- c("S40", "S69", "S70", "S57", "S32", "S38")
seqtab1 <- seqtab1[!rownames(seqtab1) %in% remove_tmp, ]

seqtabdf = as.data.frame(seqtab1)

seqtab2 <- setDT(seqtabdf, keep.rownames = T)[]
colnames(seqtab2)[1] <- 'sample'

seqtab3 <- label %>%
  right_join(seqtab2, by = 'sample')
rownames(seqtab3) <- seqtab3$ID

# drop unsued levels in factors
seqtab3 <- droplevels(seqtab3)
seqtab <- seqtab3[c(-1:-9)]

### Environmental data
sample.data <- sample_data(ps.fungi.species.f)
rownames(sample.data) <- sample.data$ID
env.chem <- sample.data[,c("totalC", "totalN", "totalP", "inorgP", "orgP")]

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
############## Q-mode ####
# For a critical view on the use of dissimilarity and distance, see Warton et al. (2012).
# Q mode analysis refers to an analysis type which focuses on relationships between objects.
# Dissimilarity coefficients indicate the degree to which objects do not resemble one another, they reach their maxima when objects share no similar variable values. 

# Untrasnformed data
OTU_unt <- dist(seqtab)

# Bray-Curtis dissimilarity : Measures for equally weighted, raw abundance data VS Hellinger distance (or X2 metric and distance) : Measures for differentially weighted, raw abundance data
# transform the data with Hellinger distance (metric, asymetrical quantitative coefficients)
# is an asymmetric distance. No weights are applied, the square roots of conditional probabilities are used as variance-stabilising data transformations. Variables with few non-zero counts are given lower weights.
OTU_hel <- decostand(seqtab, method = "hellinger")
OTU_dh <- dist(OTU_hel)

# Bray-Curtis dissimilarity
# is an asymmetrical measure often used for raw count data. Treats differences between high and low variable values equally.
OTU_bc <- vegdist(seqtab, method="bray")	

# Bray-Curtis dissimilarity matrix on log-transformed abundances
OTU_logbc <- vegdist(log1p(seqtab), method="bray")

# Bray-Curtis dissimilarity on Hellinger 
OTU_helbc <- vegdist(OTU_hel, method="bray")	

# "Anderson Log" dissimilarity matrix on log-transformed abundances
OTU_logA <- decostand(seqtab,  method = "log")
OTU_dlogA <- dist(OTU_logA)

# Check association matrices
coldiss(OTU_unt, byrank=TRUE, diag=FALSE)
coldiss(OTU_unt, byrank=FALSE, diag=FALSE)

coldiss(OTU_dh, byrank=TRUE, diag=FALSE)
coldiss(OTU_dh, byrank=FALSE, diag=FALSE)

coldiss(OTU_bc, byrank=TRUE, diag=FALSE)
coldiss(OTU_bc, byrank=FALSE, diag=FALSE)

coldiss(OTU_logbc, byrank=TRUE, diag=FALSE)
coldiss(OTU_logbc, byrank=FALSE, diag=FALSE)

coldiss(OTU_bc, byrank=TRUE, diag=FALSE)
coldiss(OTU_bc, byrank=FALSE, diag=FALSE)

coldiss(OTU_dlogA, byrank=TRUE, diag=FALSE)
coldiss(OTU_dlogA, byrank=FALSE, diag=FALSE)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Hierarchical agglomerative clustering of the ESV data ####

# Compare different transformation
OTU_ward_unt <- hclust(OTU_unt, method="ward.D2")
plot(OTU_ward_unt, hang =-1)

OTU_ward_dh <- hclust(OTU_dh, method="ward.D2")
plot(OTU_ward_dh, hang =-1)

OTU_ward_bc <- hclust(OTU_bc, method="ward.D2")
plot(OTU_ward_bc, hang =-1)

OTU_ward_logbc <- hclust(OTU_logbc, method="ward.D2")
plot(OTU_ward_logbc, hang =-1)

OTU_ward_helbc <- hclust(OTU_helbc, method="ward.D2")
plot(OTU_ward_helbc, hang =-1)

OTU_ward_dlogA <- hclust(OTU_dlogA, method="ward.D2")
plot(OTU_ward_dlogA, hang =-1)

# Compute UPGMA agglomerative clustering. Unweighted arithmetic average clustering or Unweighted Pair-Group Method using Arithmetic averages
# Use of unweighted methods because we have a systematic sampling designs so it gives equal weights to the original distances and allow us to extrapolate results to a "larger population" (Legendre et Legendre 2012, chapter 8.5)
# On Hellinger- transformed species abundances.
OTU_UPGMA <- hclust(OTU_dh, method="average")
#OTU_UPGMAbc <- hclust(OTU_logbc, method="average")

# Compute Ward's minimum variance clustering
OTU_ward <- hclust(OTU_dh, method="ward.D2")
#OTU_wardbc <- hclust(OTU_logbc, method="ward.D2")

# Compute clustering
OTU_single <- hclust(OTU_dh, method="single")

# Compute clustering
OTU_complete <- hclust(OTU_dh, method="complete")

# Compute clustering
OTU_UPGMC <- hclust(OTU_dh, method="centroid")

# Why 2 methods ? to be compared
# Plot the dendrograms
plot(OTU_UPGMA, hang =-1)
plot(OTU_ward, hang =-1)

# plot(OTU_UPGMAbc, hang=-1)
# rect.hclust(OTU_UPGMAbc, k=10)
# plot(OTU_wardbc, hang=-1)
# rect.hclust(OTU_wardbc, k=10)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Algorithm distortion measures ####

# Cophenetic matrices: Used to cpompare different classifications for the same objects (Legendre et Legendre 2012, chapter 8)
OTU_UPGMA_coph <- cophenetic(OTU_UPGMA)
OTU_ward_coph <- cophenetic(OTU_ward)
OTU_UPGMC_coph <- cophenetic(OTU_UPGMC)
OTU_single_coph <- cophenetic(OTU_single)
OTU_complete_coph <- cophenetic(OTU_complete)

#2norm, better than cophenetic correlation coefficients for mathematical reason (Mérigot et al. 201?)
norm2_UPGMA <- cl_dissimilarity(OTU_dh, OTU_UPGMA_coph, method = 'spectral')
norm2_ward <- cl_dissimilarity(OTU_dh, OTU_ward_coph, method = 'spectral')
norm2_UPGMC <- cl_dissimilarity(OTU_dh, OTU_UPGMC_coph, method = 'spectral')
norm2_single <- cl_dissimilarity(OTU_dh, OTU_single_coph, method = 'spectral')
norm2_complete <- cl_dissimilarity(OTU_dh, OTU_complete_coph, method = 'spectral')

# The algorithm that least distort the data is UPGMA with the lowest 2norm value of:
min(c(norm2_UPGMA, norm2_ward, norm2_UPGMC, norm2_single, norm2_complete))

#Cophenetic correlation coefficients (Legendre et Legendre 2012, chapter 8)
#shapiro.test(OTU_UPGMA_coph)
cor(OTU_dh, OTU_UPGMA_coph, method="spearman")
#shapiro.test(OTU_ward_coph)
cor(OTU_dh, OTU_ward_coph, method="spearman")

# However if one is interested at minimizing the sum of within-group sums of squares, Ward's method is a good candidate (Legendre et Legendre 2012, chapter 8.5). It has statistical basis. 
# Miniminzing overall variance better to look at k-means
# It is worth noting that low values of the D distance matrix (closely related objects) are less distorted (see Shepard-like diagrams below) 
#UPGMA is sequential

# Shepard-like diagrams
# représentativité des distances entre les objets  after transformation
plot(OTU_dh, OTU_UPGMA_coph, xlab="Hellinger distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,2), ylim=c(0,3),
     main=c("UPGMA", paste("Cophenetic correlation =", round(cor(OTU_dh, OTU_UPGMA_coph, method="spearman"),2)), paste("2norm =", round((norm2_UPGMA),2))))
abline(0,1)

plot(OTU_dh, OTU_ward_coph, xlab="Hellinger distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,2), ylim=c(0,3),
     main=c("Ward", paste("Cophenetic correlation =", round(cor(OTU_dh, OTU_ward_coph, method="spearman"),2)), paste("2norm =", round((norm2_ward),2))))
abline(0,1)

# ward algorithm distort the U matrix more than UPGMA as shown by the cophenetic correlation coefficients and the 2-norm

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Fusion level ####
plot(OTU_UPGMA$height, nrow(seqtab):2, type="S", 
     main="Fusion levels - Hellinger - UPGMA", 
     ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(OTU_UPGMA$height, nrow(seqtab):2, nrow(seqtab):2, col="red", cex=0.8)

plot(OTU_ward$height, nrow(seqtab):2, type="S", 
     main="Fusion levels - Hellinger - Ward", 
     ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(OTU_ward$height, nrow(seqtab):2, nrow(seqtab):2, col="red", cex=0.8)

#### Cut the trees to obtain k groups and compare the group contents using contingency tables ####
k=3 # common number of groups, small jump in the plot of fusion levels for UPGMA and ward
OTU_UPGMA_g <- cutree(OTU_UPGMA, k)
OTU_ward_g <- cutree(OTU_ward, k)

# Compare classifications by constructing contingency tables, UPGMA vs Ward
table(OTU_UPGMA_g, OTU_ward_g)
# different classif!

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Optimal number of clusters according to silhouette widths (Rousseeuw quality index) ####

#UPGMA
# Plot average silhouette widths (using UPGMA clustering) for all partitions except for the trivial partition in a single group (k=1)
asw <- numeric(nrow(seqtab)) # empty vector in which the asw values will be written
for (k in 2:(nrow(seqtab)-1)) {
  sil <- silhouette(cutree(OTU_UPGMA, k=k), OTU_dh)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)
plot(1:nrow(seqtab), asw, type="h", 
     main="Silhouette-optimal number of clusters, UPGMA, from k = 2 to n-1", 
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
sort(asw); order(-asw)

# Optimum is 3 but 10 and 2 and 10 have good value as well
# Silhouette-optimal number of clusters k = 3 
# with an average silhouette width of 0.1053669 

#WARD
# Plot average silhouette widths (using Ward clustering) for all partitions except for the trivial partition in a single group (k=1)
asw <- numeric(nrow(seqtab)) # empty vector in which the asw values will be written
for (k in 2:(nrow(seqtab)-1)) {
  sil <- silhouette(cutree(OTU_ward, k=k), OTU_dh)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)
plot(1:nrow(seqtab), asw, type="h", 
     main="Silhouette-optimal number of clusters, Ward, from k = 2 to n-1", 
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
sort(asw); order(-asw)

# Optimum is 2 but 15 has good value as well
# Silhouette-optimal number of clusters k = 2 
# with an average silhouette width of 0.1295747 

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Optimal number of clusters according to Mantel statistic (Pearson) ####

# Function to compute a binary distance matrix from groups
grpdist <- function(X)
{
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}

#UPGMA
# Run based on the UPGMA clustering
kt <- data.frame(k=1:nrow(seqtab), r=0)
for (i in 2:(nrow(seqtab)-1)) {
  gr <- cutree(OTU_UPGMA, i)
  distgr <- grpdist(gr)
  mt <- cor(OTU_dh, distgr, method="pearson")
  kt[i,2] <- mt
}
# kt
(k.best <- which.max(kt$r))

plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - UPGMA", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
sort(kt$r); order(-kt$r)
# Optimal number of clusters for UPGMA is 11 but 10 and 12 values around those are good as well in term of r
# Mantel-optimal number of clusters k = 11 
# with a matrix linear correlation of 0.5963291 

#WARD
# Run based on the Ward clustering
kt <- data.frame(k=1:nrow(seqtab), r=0)
for (i in 2:(nrow(seqtab)-1)) {
  gr <- cutree(OTU_ward, i)
  distgr <- grpdist(gr)
  mt <- cor(OTU_dh, distgr, method="pearson")
  kt[i,2] <- mt
}
#kt
(k.best <- which.max(kt$r))

plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Ward", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
sort(kt$r); order(-kt$r)
# Optimal number of clusters for Ward is 5
# Mantel-optimal number of clusters k = 5 
# with a matrix linear correlation of 0.5476189

# If we compare both methods (silhouette widths or Mantel test) 
# optimum number of cluster is different for UPGMA. 10 groups seem to be the best compromise
# optimum number of cluster is different for Ward. 3 or 10 groups seem to be the best compromise

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Silhouette plot of the final partition ####
# Choose the number of clusters
k1 <- 3
k2 <- 10

# Silhouette plot for Ward with k=5
cutg <- cutree(OTU_ward, k=k1)
sil <- silhouette(cutg, OTU_dh)
rownames(sil) <- row.names(seqtab)
plot(sil, main="Silhouette plot - Ward - k=3", 
     cex.names=0.8, col=2:(k1+1), nmax=100)

# Silhouette plot for Ward with k=10
cutg <- cutree(OTU_ward, k=k2)
sil <- silhouette(cutg, OTU_dh)
rownames(sil) <- row.names(seqtab)
plot(sil, main="Silhouette plot - Ward - k=10", 
     cex.names=0.8, col=2:(k2+1), nmax=100)

# Silhouette plot for UPGMA with k=10
cutg <- cutree(OTU_UPGMA, k=k2)
sil <- silhouette(cutg, OTU_dh)
rownames(sil) <- row.names(seqtab)
plot(sil, main="Silhouette plot - UPGMA - k=10", 
     cex.names=0.8, col=2:(k2+1), nmax=100)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Final dendrogram with the selected groups ####
OTU_UPGMA_final <- reorder(OTU_UPGMA, OTU_dh)
OTU_ward_final <- reorder(OTU_ward, OTU_dh)

# Plot reordered dendrogram for Ward
plot(OTU_ward_final, hang=-1, xlab="Groups", sub="", 
     ylab="Height", main="Hellinger - Ward (reordered)")
rect.hclust(OTU_ward_final, k=5)

# Plot reordered dendrogram for Ward
plot(OTU_ward_final, hang=-1, xlab="Groups", sub="", 
     ylab="Height", main="Hellinger - Ward (reordered)")
rect.hclust(OTU_ward_final, k=k1)

# Plot reordered dendrogram for Ward
plot(OTU_ward_final, hang=-1, xlab="Groups", sub="", 
     ylab="Height", main="Hellinger - Ward (reordered)")
rect.hclust(OTU_ward_final, k=k2)

# Plot reordered dendrogram for UPGMA
plot(OTU_UPGMA_final, hang=-1, xlab="Groups", sub="", 
     ylab="Height", main="Hellinger - UPGMA (reordered)")
rect.hclust(OTU_UPGMA_final, k=11)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Heat map ####
# Heat map of the distance matrix ordered with the dendrogram
dend_ward <- as.dendrogram(OTU_ward_final)
heatmap(as.matrix(OTU_dh), Rowv=dend_ward, symm=TRUE, margin=c(4,4), main = "Heatmap - sites - Ward")

dend_UPGMA <- as.dendrogram(OTU_UPGMA_final)
heatmap(as.matrix(OTU_dh), Rowv=dend_UPGMA, symm=TRUE, margin=c(4,4), main = "Heatmap - sites - UPGMA")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### k-means partitioning of the pre-transformed species data ####
# = Non-hierarchical methods
# With 3 or 10 groups
OTU_kmeans3 <- kmeans(OTU_hel, centers=3, nstart=100)
OTU_kmeans10 <- kmeans(OTU_hel, centers=10, nstart=100)

# Comparison with the k-group classification derived from Ward clustering:
table(OTU_kmeans3$cluster, OTU_ward_g) # almost the same grouping

# k-means partitioning, 2 to 15 groups
OTU.KM.cascade <- cascadeKM(OTU_hel, inf.gr=2, sup.gr=15, iter=1000, # iter = 100000, need high number of iteration for reproducibility (but time-consuming)
                            criterion="ssi")
OTU.KM.cascade$results
plot(OTU.KM.cascade, sortg=TRUE, main = 'cascadeKM')

# Reorder the sites according to the k-means result
# seqtab[order(OTU_kmeans10$cluster),]
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Species indicator values (Dufrene and Legendre) ####

#If the site classification vector is obtained independently of species data, the significance of statistical tests carried out on the indicator species will be meaningful. For example, one could classify the sites using environmental data before indicator species analysis. An example is found in Borcard et al. 2011
#library(indicspecies)
#iva = multipatt(seqtab, sample.data$mycohoriz, max.order = 1, control = how(nperm=999)) # Dufrêne and Legendre (1997) calculates the IndVal index but with combinations of site groups, as explained in De Cáceres et al. (2010)
#summary(indicval)

# Indicator species for this typology of the sites
iva <- indval(seqtab, seqtab3$myco_horiz)
#iva <- indval(seqtab, env_ward_5) 
#iva <- indval(seqtab, env_kmeans7) 
#iva <- indval(seqtab, sample.data$horizon)


# Table of the significant indicator species
gr <- iva$maxcls[iva$pval <= 0.05]
iv <- iva$indcls[iva$pval <= 0.05]
pv <- iva$pval[iva$pval <= 0.05]
fr <- apply(seqtab > 0, 2, sum)[iva$pval <= 0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg

fidg$SV <- rownames(fidg)
SV_taxo$SV <- rownames(SV_taxo)

groups <- data.frame(group = as.integer(rep(1:15)),
                     myco_horiz=levels(seqtab3$myco_horiz))

indic.species <- fidg %>%
  left_join(SV_taxo, by = 'SV') %>%
  left_join(groups, by = 'group')
indic.species

indic.max.species <- indic.species %>% 
  group_by(group) %>% 
  top_n(1, indval)
indic.max.species

# Export the result to a CSV file (to be opened in a spreadsheet)
# write.csv(indic.species, "exploration_R/saved/cluster/IndVal_myco_horiz.csv")
# write.csv(indic.max.species, "exploration_R/saved/cluster/IndValmax_myco_horiz.csv")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
############## R-mode ####

# Center and scale = standardize variables (z-scores)
env.z <- decostand(env.chem, "standardize")

# Euclidean distance matrix of the standardized soil chemsitry data
env.de <- dist(scale(env.z))
coldiss(env.de, nc=6, diag=FALSE) # we can see two big cluster and 

# Pearson r linear correlation among environmental variables
env.spearman <- cor(env.z, method = "spearman")
round(env.spearman, 2)

ggpairs(env.z, upper = list(continuous = wrap('cor', method = "spearman")))
# we decided not to keep totalP as it higlhy correlated to orgP (.978) to avoid variance-inflation (Borcard et al. 2011)
# see VIF (later in the script)
env.z<- env.z[,c("totalC", "totalN", "inorgP", "orgP")]

env_UPGMA <- hclust(env.de, method="average")
env_ward <- hclust(env.de, method="ward.D2")

#### Optimal number of clusters according to Mantel statistic (Pearson) ####
grpdist <- function(X)
{
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}

# Run based on the Ward clustering
kt <- data.frame(k=1:nrow(env.z), r=0)
for (i in 2:(nrow(env.z)-1)) {
  gr <- cutree(env_ward, i)
  distgr <- grpdist(gr)
  mt <- cor(env.de, distgr, method="pearson")
  kt[i,2] <- mt
}
#kt
(k.best <- which.max(kt$r))

plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - ward", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
# Mantel-optimal number of clusters k = 2 
# with a matrix linear correlation of 0.7989963 

env_ward_best <- cutree(env_ward, k.best)
env_ward_5 <- cutree(env_ward, 5)

# Run based on the UPGMA clustering
kt <- data.frame(k=1:nrow(env.z), r=0)
for (i in 2:(nrow(env.z)-1)) {
  gr <- cutree(env_UPGMA, i)
  distgr <- grpdist(gr)
  mt <- cor(env.de, distgr, method="pearson")
  kt[i,2] <- mt
}
#kt
(k.best <- which.max(kt$r))

plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - UPGMA", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
# Mantel-optimal number of clusters k = 3 
# with a matrix linear correlation of 0.8139792 

plot(env_UPGMA, hang =-1)
rect.hclust(env_UPGMA, k=2)

plot(env_ward, hang =-1)
rect.hclust(env_ward, k=2)
# Grouping that separate among organic and mineral horizons (except for a few H, Ah)
rect.hclust(env_ward, k=5)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### k-means partitioning  ####
# Non-hierarchical methods

# k-means partitioning, 2 to 15 groups
env.KM.cascade <- cascadeKM(env.z, inf.gr=2, sup.gr=20, iter=100000, # need high number of iteration for reproducibility (but time-consuming)
                            criterion="ssi")
env.KM.cascade$results # 13 groups
plot(env.KM.cascade, sortg=TRUE, main = 'cascadeKM')

env_kmeans7 <- kmeans(env.z, centers=7, nstart=100)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
############ variance paritioning ####
# matrices Y
env.chem <- sample.data[,c("totalC", "totalN", "inorgP", "orgP")]
env.z <- decostand(env.chem, "standardize")

#spa <- sample.data[,c("elev", "slope")]
#spa.z <- decostand(spa, "standardize")

horiz <- as(sample.data[,"horiz"], "data.frame")
myco <- as(sample.data[,"sp"], "data.frame")
block <- as(sample.data[,"block.x"], "data.frame")

# matrix X
OTU_hel <- decostand(seqtab, method = "hellinger")
OTU_hel <- vegdist(OTU_hel, method = "bray")

# Converting Categorical Columns into Multiple Binary Columns and trimming latent variables
env.horiz <- as.data.frame(model.matrix(~horiz -1, data=horiz))
env.myco <- as.data.frame(model.matrix(~sp -1, data=myco))
env.block <- as.data.frame(model.matrix(~block.x -1, data=block))

# deleting redundant variables
env.horiz <- env.horiz[,1:4]
env.myco <- env.myco[,1:2]
env.block <- env.block[,1:4]

# Appying mutliple partial RDA to partitionate the variance
var.part.all <- varpart(OTU_hel, env.z, env.block, env.horiz, env.myco)

plot(var.part.all, bg = 1:4, digits = 2, cutoff = 0.0001)

# Should be same results than var.part
rda.all <- rda(OTU_hel, cbind(env.z, env.block, env.horiz, env.myco)) # Inertia is variance, proportion is r.squared (in vegan)
RsquareAdj(rda.all) # same than var.part
anova(rda.all)

# Use function ‘rda’ to test significance of fractions of interest. Semi-Partial correlation
prda.X1 <- rda(OTU_hel, env.z, cbind(env.block, env.horiz, env.myco))
RsquareAdj(prda.X1)
anova(prda.X1)

prda.X2 <- rda(OTU_hel, env.block , cbind(env.z, env.horiz, env.myco))
RsquareAdj(prda.X2)
anova(prda.X2)

prda.X3 <- rda(OTU_hel, env.horiz , cbind(env.z, env.block, env.myco))
RsquareAdj(prda.X3)
anova(prda.X3)

prda.X4 <- rda(OTU_hel, env.myco , cbind(env.z, env.block, env.horiz))
RsquareAdj(prda.X4)
anova(prda.X4)

# Partial correlation
rda.X1 <- rda(OTU_hel, env.z)
RsquareAdj(rda.X1)
anova(rda.X1)

rda.X2 <- rda(OTU_hel, env.block)
RsquareAdj(rda.X2)
anova(rda.X2)

rda.X3 <- rda(OTU_hel, env.horiz)
RsquareAdj(rda.X3)
anova(rda.X3)

rda.X4 <- rda(OTU_hel, env.myco)
RsquareAdj(rda.X4)
anova(rda.X4)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Ordination ####

# rename mycorrhizal types in plots
myco_names <- c(`AM` = "AM", `Mix` = "Mixed", `ECM` = "EcM")

### PCoA on Hellinger data with Bray dissimilarity
ord.pcoa.hel <- ordinate(ps.fungi.hel, method="PCoA", distance = "bray")
#trace = 22.2362
evals <- ord.pcoa.hel$values$Relative_eig

# Selection of axes with broken stick model
n <- length(evals)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n)
{
  bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
}
bsm$p <- 100*bsm$p/n

# Plot eigenvalues and % of variance for each axis
barplot(evals, main="PCA Eigenvalues", col="bisque")
abline(h=mean(evals), col="red")		# average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(100*evals/sum(evals),bsm$p[n:1])), beside=TRUE, 
        main="% variance", col=c("bisque",2))
legend("topright", c("% eigenvalue", "Broken stick model"), 
       pch=15, col=c("bisque",2), bty="n")
# we should keep the first 6 PC. Might be due to the high number of PC
# keep PCo1 and PCo2 for simplicity purpose

plot_ordination(ps.fungi.hel, ord.pcoa.hel, type="biplot", color = "horiz") + 
  scale_shape_manual(values=c(4,15,16,17)) +
  geom_point(size=3) +
  coord_fixed(sqrt(evals[2] / evals[1]))

plot_ordination(ps.fungi.hel, ord.pcoa.hel, type="samples", color = "horiz") + 
  geom_point(size=2) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  facet_wrap(~myco) +
  stat_ellipse(type = "t", level = 0.9) + # 90% confidenc einterval, assumes a multivariate t-distribution
  scale_color_discrete(name="Horizon") +
  facet_wrap(~myco, labeller = as_labeller(myco_names))+
  xlab("PCo 1 (19.7% of variation)") +
  ylab("PCo 2 (8.6% of variation)")

# we have to take into account that the second eigenvalue is always smaller than the first, 
# sometimes considerably so, thus we normalize the axis norm ratios to the relevant eigenvalue ratios (coord_fixed).

### PCoA on presence-absence (Sorensen)
ord.pcoa.bin <- ordinate(ps.fungi.hel, method="PCoA", distance = "bray", binary = TRUE)
# trace = 20.83582
evals <- ord.pcoa.hel$values$Relative_eig
plot_ordination(ps.fungi.hel, ord.pcoa.bin, type="samples", color = "horiz") + 
  geom_point(size=2) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  facet_wrap(~myco) +
  stat_ellipse(type = "t", level = 0.9) + # 90% confidenc einterval
  scale_color_discrete(name="Horizon") +
  facet_wrap(~myco, labeller = as_labeller(myco_names)) +
  xlab("PCo 1 (21.8% of variation)") +
  ylab("PCo 2 (9.3% of variation)")

### NMDS
ord.nmds.hel <- ordinate(ps.fungi.hel, method="NMDS", try = 100, distance = "bray")
# stress 0.1938509
plot_ordination(ps.fungi.hel, ord.nmds.hel, type="biplot", color = "horiz", shape = "myco") + 
  scale_shape_manual(values=c(4,15,16,17)) +
  geom_point(size=4)

### PCA
ord.pca.hel <- ordinate(ps.fungi.hel, method="RDA")
# inertia = 1.442
plot_ordination(ord.pca.hel, type="samples", color = "horiz", shape = "myco") + 
  geom_point(size=4)

### CCA on Hellinger transformation
ord.cca.hel <- ordinate(ps.fungi.hel, formula = (ps.fungi.hel ~ totalC + totalN + inorgP + orgP), method="CCA")

# Plot eigenvalues CCA and CA
plot_ordination(ps.fungi.hel, ord.cca.hel, type="scree")
# get eigenvalues
evals <- ord.cca.hel$CCA$eig

# variance inflation factors 
vif.cca(ord.cca.hel)
#  totalC   totalN   inorgP     orgP 
#6.978471 7.221011 1.425696 4.356632 

# Add the environmental variables as arrows
arrowmat = vegan::scores(ord.cca.hel, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(arrowmat)
rownames(arrowdf) <- c("Total C", "Total N", "Pi", "Po")
labels = c("Total C", "Total N", "Pi", "Po")
arrowdf <- data.frame(labels = c("Total C", "Total N", "Pi", "Po"), arrowmat)

plot_ordination(ps.fungi.hel, ord.cca.hel, type="samples",  color = "horiz") + 
  geom_point(size = 2) +
  stat_ellipse(type = "t", level = 0.9) + # 90% confidence interval
  coord_fixed(sqrt(evals[2] / evals[1])) +
  #  scale_color_manual(values=c("khaki4","sienna4","black", "grey60", "brown3", "firebrick1")) +
  geom_segment(aes(xend = CCA1, yend = CCA2, x = 0, y = 0, shape = NULL, color = NULL, label = labels), 
               size = 0.6, 
               data = arrowdf, color = "black", arrow = arrow(length = unit(0.025, "npc"))) + 
  geom_text(aes(x = CCA1 - 0.4, y = CCA2*1.45, shape = NULL, color = NULL, label = labels), 
            size = 4, 
            data = arrowdf) +
  scale_color_discrete(name="Horizon") +
  facet_wrap(~myco, labeller = as_labeller(myco_names)) +
  xlab("CCA Axis 1 (5% of variation)") +
  ylab("CCA Axis 2 (2.6% of variation)")

RsquareAdj(ord.cca.hel) # 0.05421339
anova.cca(ord.cca.hel) # pr 0.001 

plot_ordination(ps.fungi.hel, ord.cca.hel, type="samples",  color = "myco") + 
  geom_point(size = 3) +
  theme(legend.position="none") +
  coord_fixed() +
  stat_ellipse()

### partial CCA
ord.pcca.hel <- ordinate(ps.fungi.hel, formula = (ps.fungi.hel ~ totalC + totalN + inorgP + orgP + Condition(horiz)), method="CCA")

# Plot eigenvalues
plot_ordination(ps.fungi.hel, ord.pcca.hel, type="scree")

plot_ordination(ps.fungi.hel, ord.pcca.hel, type="samples",  color = "myco") + 
  geom_point(size = 3) +
  theme(legend.position="none") +
  coord_fixed() +
  stat_ellipse()

RsquareAdj(ord.pcca.hel) # 0.006750063
anova.cca(ord.pcca.hel) # pr 0.077

### RDA
ord.rda.hel <- ordinate(ps.fungi.hel, formula = (ps.fungi.hel ~ totalC + totalN + inorgP + orgP), method="RDA")

# Plot eigenvalues
plot_ordination(ps.fungi.hel, ord.rda.hel, type="scree")
# get eigenvalues
evals <- ord.rda.hel$CCA$eig

# Add the environmental variables as arrows
arrowmat = vegan::scores(ord.rda.hel, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

plot_ordination(ps.fungi.hel, ord.rda.hel, type="biplot",  color = "horiz", shape = "myco") + 
  geom_point(size = 2) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  #  scale_color_manual(values=c("khaki4","sienna4","black", "grey60", "brown3", "firebrick1")) +
  scale_shape_manual(values=c(1,15,16,17)) +
  geom_segment(aes(xend = RDA1, yend = RDA2, x = 0, y = 0, shape = NULL, color = NULL, label = labels), 
               size = 0.8, 
               data = arrowdf, color = "black", arrow = arrow(length = unit(0.025, "npc"))) + 
  geom_text(aes(x = RDA1 - 0.1, y = RDA2, shape = NULL, color = NULL, label = labels), 
            size = 6, 
            data = arrowdf)

RsquareAdj(ord.rda.hel) # 0.1572035
anova.cca(ord.rda.hel) # pr 0.001

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#
#### Permanova ####
metadata <- as(sample_data(ps.fungi.hel), "data.frame")

# horizon and myco effect with interaction
# method = euclidean because data already Hellinger transformed
adonis(distance(ps.fungi.hel, method="euclidean") ~ horiz * myco, strata = metadata$block.y,
       data = metadata,
       permutations = 9999)

# Sorensen presence-absence
adonis(distance(ps.fungi.hel, method = "bray", binary =TRUE) ~ horiz * myco, strata = metadata$block.y, 
       data = metadata,
       permutations = 9999)

# Homogeneity of dispersion test for Hellinger data
betadisp <- betadisper(distance(ps.fungi.hel, method="euclidean"), metadata$myco)
permutest(betadisp, permutations = 9999) # p-value 0.2031
betadisp <- betadisper(distance(ps.fungi.hel, method="euclidean"), metadata$horiz)
permutest(betadisp, permutations = 9999) # p-value 0.9929

# Homogeneity of dispersion test for Sorensen data
betadisp <- betadisper(distance(ps.fungi.hel, method = "bray", binary =TRUE), metadata$myco)
permutest(betadisp, permutations = 9999) # p-value 0.1869
betadisp <- betadisper(distance(ps.fungi.hel, method = "bray", binary =TRUE), metadata$horiz)
permutest(betadisp, permutations = 9999) # p-value 0.009

## pairwise comparison for Sorensen data
pairwise.perm.manova(distance(ps.fungi.hel, method = "bray", binary =TRUE), metadata$myco, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method

# pairwise.perm.manova(distance(ps.fungi.hel, method = "bray", binary =TRUE), metadata$horiz, nperm=99999,
#                      p.method = "fdr")

## pairwise comparaison for Hellinger data for horizon within forest
# subset by forest
ps.AM <- subset_samples(ps.fungi.hel, myco == "AM")
meta.AM <- as(sample_data(ps.AM), "data.frame")

ps.ECM <- subset_samples(ps.fungi.hel, myco == "ECM")
meta.ECM <- as(sample_data(ps.ECM), "data.frame")

ps.Mix <- subset_samples(ps.fungi.hel, myco == "Mix")
meta.Mix <- as(sample_data(ps.Mix), "data.frame")

pairwise.perm.manova(distance(ps.AM, method = "euclidean"), meta.AM$horiz, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.AM, method = "euclidean"), meta.AM$horiz)
permutest(betadisp) # p-value 0.712

pairwise.perm.manova(distance(ps.ECM, method = "euclidean"), meta.ECM$horiz, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.ECM, method = "euclidean"), meta.ECM$horiz)
permutest(betadisp) # p-value 0.907

pairwise.perm.manova(distance(ps.Mix, method = "euclidean"), meta.Mix$horiz, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.Mix, method = "euclidean"), meta.Mix$horiz)
permutest(betadisp) # p-value 0.829

# pairwise comparaison for Hellinger data for forest among horizon
ps.L <- subset_samples(ps.fungi.hel, horiz == "L")
meta.L <- as(sample_data(ps.L), "data.frame")
ps.F <- subset_samples(ps.fungi.hel, horiz == "F")
meta.F <- as(sample_data(ps.F), "data.frame")
ps.H <- subset_samples(ps.fungi.hel, horiz == "H")
meta.H <- as(sample_data(ps.H), "data.frame")
ps.Ae <- subset_samples(ps.fungi.hel, horiz == "Ae")
meta.Ae <- as(sample_data(ps.Ae), "data.frame")
ps.B <- subset_samples(ps.fungi.hel, horiz == "B")
meta.B <- as(sample_data(ps.B), "data.frame")

pairwise.perm.manova(distance(ps.L, method = "euclidean"), meta.L$myco, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.L, method = "euclidean"), meta.L$myco)
permutest(betadisp) # p-value 0.61

pairwise.perm.manova(distance(ps.F, method = "euclidean"), meta.F$myco, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.F, method = "euclidean"), meta.F$myco)
permutest(betadisp) # p-value 0.891

pairwise.perm.manova(distance(ps.H, method = "euclidean"), meta.H$myco, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.F, method = "euclidean"), meta.F$myco)
permutest(betadisp) # p-value 0.89

pairwise.perm.manova(distance(ps.Ae, method = "euclidean"), meta.Ae$myco, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.Ae, method = "euclidean"), meta.Ae$myco)
permutest(betadisp) # p-value 0.624

pairwise.perm.manova(distance(ps.B, method = "euclidean"), meta.B$myco, nperm=99999,
                     p.method = "fdr") # Benjamini-Hochberg adjustment method
betadisp <- betadisper(distance(ps.B, method = "euclidean"), meta.B$myco)
permutest(betadisp) # p-value 0.696


# horizon and myco effect with interaction by GENUS
ps.fungi.genus.hel = tax_glom(ps.fungi.hel, "Genus", NArm = FALSE)
adonis(distance(ps.fungi.genus.hel, method="euclidean") ~ horiz * myco, strata = metadata$block.y, 
       data = metadata,
       permutations = 9999)

# horizon and myco effect with interaction by Family
ps.fungi.family.hel = tax_glom(ps.fungi.hel, "Family", NArm = FALSE)
adonis(distance(ps.fungi.family.hel, method="euclidean") ~ horiz * myco, strata = metadata$block.y, 
       data = metadata,
       permutations = 9999)

# horizon and myco effect with interaction by Family
ps.fungi.order.hel = tax_glom(ps.fungi.hel, "Order", NArm = FALSE)
adonis(distance(ps.fungi.order.hel, method="euclidean") ~ horiz * myco, strata = metadata$block.y, 
       data = metadata,
       permutations = 9999)
# interaction horiz:myco no more significant


