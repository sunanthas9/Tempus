####################################################################################################################
####################################################################################################################
#### Tempus bioinformatics technical challenge: Variant annotation
#### Author: Sunantha Sethuraman
#### Date: 09/07/2017
#### Use: Annotates the given vcf file with information from online resources
####################################################################################################################
####################################################################################################################

# Set working directory

setwd("/ufrc/renne/sunantha.s/research/VCF")

# Load bioconductor and update default packages
source("https://bioconductor.org/biocLite.R")
biocLite()

# Load necessary packages
library(VariantAnnotation)
library(cgdv17)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(PolyPhen.Hsapiens.dbSNP131)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(VariantFiltering)
library(MafDb.ExAC.r1.0.hs37d5)
library(MafDb.1Kgenomes.phase1.hs37d5)
library(PolyPhen.Hsapiens.dbSNP131)
library(SIFT.Hsapiens.dbSNP137)
library(phastCons100way.UCSC.hg19)
library(GenomicRanges)

##################################### Preprocessing #########################################################

# The input file requires preprocessing since it does not contain the ID information. ID is required later to 
# obtain ExAc information. The dbSNP ids (aka'rs' ids) are added to the input file.

# Obtain basic information from the file
temp <- readLines("Challenge_data_edited.vcf")
start <- grep("CHROM", temp)
Information <- read.csv("Challenge_data_edited.vcf", skip=start-1, sep="\t")
Information <- Information[,c(1,2,4,5)]

# Obtain the position information from the file
vcf <- readVcf("Challenge_data_edited.vcf") 
locs <- rowRanges(vcf)
seqlevels(locs) <- seqlevels(locs)[1:25]
newStyle= mapSeqlevels(seqlevels(locs), style = "NCBI")
locs <- renameSeqlevels(locs, newStyle)

# Fetch the dbSNP ids, called 'rs ids' using the SNPlocs.Hsapiens.dbSNP144.GRCh37 package
rsids <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, locs, type="start",drop.rs.prefix=FALSE) 
rsid_df <- data.frame(matrix(0, nrow=length(rsids),ncol=2))
colnames(rsid_df) <- c("POS","RSIDS")
rsid_df$RSIDS <- values(rsids)$RefSNP_id
rsid_df$POS <- start(ranges(rsids))
RowNameDf <- merge(rsid_df, Information, by="POS", all.y=T)
RowNameDf <- RowNameDf[order(RowNameDf$X.CHROM),]
RowNameDf <- RowNameDf[!duplicated(RowNameDf$POS,fromLast=T),] #Retain only the first rs id if two or more are found for the same position
Unknown <- is.na(RowNameDf$RSIDS) # For SNPs without reference, use the VariantFiltering norm of "Chr:POS_REF" format for naming
RowNameDf$RSIDS[Unknown] <- paste0(RowNameDf$X.CHROM[Unknown],":",RowNameDf$POS[Unknown],"_",RowNameDf$REF[Unknown])

# Add rs ids to file and write vcf
names(rowRanges(vcf)) <- RowNameDf$RSIDS
writeVcf(vcf,"Preprocessed.vcf", row.names=F)

##################################### Annotation #########################################################

# Perform variant annotation using the VariantFiltering package available on bioconductor
bgzip("Preprocessed.vcf", overwrite=T)
vfpar <- VariantFilteringParam("Preprocessed.vcf.bgz", otherAnnotations = c("MafDb.ExAC.r1.0.hs37d5")) # Use ExAC annotation package in addition to default packages
uind <- unrelatedIndividuals(vfpar)

# Write out a summary of annotated variants and create an output file of annotated variants
write.csv(summary(uind), "Summary.csv", row.names = F)
reportVariants(uind, type="tsv", "Variants.tsv")

##################################### Clean-up ############################################################

# Clean-up the output to present information in an user-friendly way
data <- read.csv("Variants.tsv", sep="\t")
Anno <- merge(Information, data, by="POS")
Anno$PER_REF <- 100*Anno$RO/Anno$DP
Anno$PER_ALT <- 100*Anno$AO/Anno$DP
Anno <- Anno[,c("CHR","POS","REF","ALT","DP","RO","AO","PER_REF", "PER_ALT","SAMPLEID","TYPE","TYPE.1","AFExAC","LOCATION","GENE","TXNAME",
          "CONSEQUENCE")]

colnames(Anno) <- c("Chromosome", "Position","ReferenceAllele", "VariantAllele","SequencingDepth", "ReferenceReads",
                    "VariantReads", "PercentageReferenceReads", "PercentageVariantReads", "Sample", "Type", "Class", "AlleleFrequency_ExAC", "Location", "GeneName",
                    "TranscriptID_UCSC", "ProteinConsequence" )

write.csv(Anno, "AnnotatedTable.csv", row.names=F)

###################################### END OF SCRIPT #########################################################