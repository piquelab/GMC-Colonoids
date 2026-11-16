##Make QC Plots for GMCcol_deep
##Script original from Sarah Dubasi
##Edited by Shreya Nirmalan on 5/16/2023
library(ggplot2)


### Plot1: RawReads vs Unique Reads
########### 
RawReads <-read.table("/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/bams/total_reads.txt",header=F, stringsAsFactors=F)

AllReads<-read.table("/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/bams/all_counts.txt",header=F, stringsAsFactors=F)

colnames(RawReads) <- c("Filename", "RawReads")

colnames(AllReads) <- c("Filename", "Sorted", "Quality", "Clean")

##Bringing treatment information so plots can be organized by treatment
cov.file <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/covariates/GMCcol_deep_covariates.txt"
cv <- read.table(cov.file,  sep="\t", header=T, quote='"', comment="")
cv$SampleID <- gsub("_","-", cv$SampleID)
ratio_file <- read.table("/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/counts/counts2/GMCcol_deep_alignment_HTseq.txt", stringsAsFactors=F, header=T, comment="")
cv <- cbind(cv, fraction_mapped = ratio_file$fraction_mapped)
colnames(cv)[1] <- "Filename"


temp <- merge(RawReads,AllReads, "Filename")
Merged <- merge(temp,cv, "Filename")
Merged$RawReads <- Merged$RawReads * 2
Merged$NonDup <- Merged$Clean / Merged$RawReads #Proportion of reads non-duplicates

pdf("./GMCcol_deep_QC_Unique_vs_RawReads_All2_By_Line.pdf")
p <- ggplot(Merged, aes(log10(RawReads), NonDup, color= Colonoid_line))
PBP_Unique_vs_RawReads_All_Plate <- p + geom_point(aes(shape= Treatment)) + ylim(0, 1.0) + labs(title = "Unique Read Fraction vs Raw Reads", x = "log10(Raw Read Count)", y = "Unique Read Fraction")
print(PBP_Unique_vs_RawReads_All_Plate)
dev.off()

##################################################################################################################
##### Plot2: 

counts <- AllReads

#colnames(counts) <- c("Filename", "Sorted", "Quality", "Clean")
counts$Sorted.adj <- counts$Sorted - counts$Quality
counts$Quality.adj <- counts$Quality - counts$Clean

total_reads <- RawReads
colnames(total_reads) <- c("Filename", "Total.Reads")

library(reshape2)
countsMelted <- melt(counts, id.vars = "Filename", measure.vars = c( "Sorted.adj", "Quality.adj", "Clean"))

#countsMelted <- melt(counts, id.vars = "Filename", measure.vars = c( "Quality.adj", "Clean"))

##library(ggplot2)

levels(countsMelted$variable)

countsMelted$variable <- factor(countsMelted$variable, levels = c("Clean", "Quality.adj", "Sorted.adj"))

#countsMelted$variable <- factor(countsMelted$variable, levels = c("Clean", "Quality.adj"))

countsMelted<-merge(countsMelted, total_reads, by="Filename")

counts_Melted<-rbind(countsMelted[countsMelted$variable == "Clean",], countsMelted[countsMelted$variable == "Quality.adj",])
countsMelted$Total.Reads<-countsMelted$Total.Reads *2 #Because Hisat reports fragments instead of reads

pdf("./GMCcol_deep_QC_Melted.pdf", height = 20, width = 20)
p <- ggplot(data = countsMelted, aes(x = Filename, y = value, fill = factor(variable))) + geom_bar(stat="identity", position= position_stack (reverse=TRUE)) +
  geom_point(data=countsMelted, aes(y = Total.Reads, x=Filename))+ theme(text = element_text(size=36))+
  coord_flip()

p
dev.off()


######### Plot 3: Quality versus Total Reads ####################################################################
tempTable <- merge(RawReads, AllReads, "Filename")
tempTable <- merge(tempTable, cv, "Filename")
tempTable$RawReads <- tempTable$RawReads * 2
tempTable$QvTR <- tempTable$Quality / tempTable$RawReads #proportion of quality reads to total reads
pdf("./GMCcol_deep_QC_Quality_versus_RawReads_By_Line.pdf")
p <- ggplot(tempTable, aes(log10(RawReads), QvTR, color=Colonoid_line))
QC_Quality_v_RawReads <- p + geom_point(aes(shape= Treatment)) + ylim(0, 1.0) + labs(title = "Quality Read Fraction vs Raw Reads", x = "log10(Raw Read Count)", y = "Quality Read Fraction")
print(QC_Quality_v_RawReads)
dev.off()

######## Plot 4: Clean versus Quality Reads #####################################################################

tempTableCleanvQuality <- merge(RawReads, AllReads, "Filename")
tempTableCleanvQuality <- merge(tempTableCleanvQuality, cv, "Filename")
tempTableCleanvQuality$CleanVQuality <- tempTableCleanvQuality$Clean / tempTableCleanvQuality$Quality
pdf("./GMCcol_deep_QC_Clean_Versus_QualityReads_By_Line.pdf")
p <- ggplot(tempTableCleanvQuality, aes(log10(Quality), CleanVQuality, color=Colonoid_line))
QC_Clean_v_Quality <- p + geom_point(aes(shape = Treatment)) + ylim(0, 1.0) + labs(title = "Clean Versus Quality Reads", x = "log10(Quality Read Count)", y = "Clean Versus Quality Read Fraction")
print(QC_Clean_v_Quality)
dev.off()



########### Plot 5: Fraction mapped versus Raw Reads ###############################################################

tempTableRatio <- merge(RawReads, cv, "Filename")
tempTableRatio$RawReads <- tempTableRatio$RawReads * 2
pdf("GMCcol_deep_fraction_mapped_versus_Raw_reads_By_Line.pdf")
p <- ggplot(tempTableRatio, aes(log10(RawReads), fraction_mapped, color=Colonoid_line)) +
      geom_point(aes(shape=Treatment)) +
      ylim(0, 1.0) +
      labs(title = "Fraction Mapped Versus Raw Reads", x = "log10(Raw Read Count)", y = "Fraction Mapped to Exons")
print(p)
dev.off()


##make covariates file with ratio and read info
library(tidyverse)
temp <- merge(RawReads,AllReads, "Filename")
Merged <- merge(temp,cv, "Filename")
Merged$TotalReads <- Merged$RawReads * 2
Merged$NonDuplicated_reads <- Merged$Clean / Merged$TotalReads
Merged <- rename(Merged, Ratio = fraction_mapped)
Merged <- rename(Merged, SampleID = Filename)
write.table(Merged, file="/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/covariates/GMCcol_deep_covariates.txt", quote=F, sep="\t")