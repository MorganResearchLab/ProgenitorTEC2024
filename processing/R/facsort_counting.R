library(ggplot2)
library(ggsci)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)

# compute proportions of TEC in different sort fractions
# Read in cell numbers from 10 replicates in each fraction
# Use these to back-calculate the expected numbers of 
# Ly6d+/- intertypical TEC

cell.counts <- read.csv2("/Volumes/Single_Cell_Genetics-CF10180-29/AgeingExperiment/Perinatal_manuscript/data/FAC_counts.csv",
                         header=TRUE, sep=",")
# convert to wide format
cell.df <- cell.counts %>% 
  pivot_wider(names_from=c("TECSort", "Fraction"), values_from=Ncells) %>%
  group_by(Day, Replicate) %>% 
  mutate("Total.Cells"=sum(c(cTEC_ZsGp, cTEC_ZsGn, mTEC_ZsGp, mTEC_ZsGn)), .keep="all") %>%
  mutate("p.cTEC_ZsGpos"=cTEC_ZsGp/Total.Cells, "p.mTEC_ZsGpos"=mTEC_ZsGp/Total.Cells,
         "p.cTEC_ZsGneg"=cTEC_ZsGn/Total.Cells, "p.mTEC_ZsGneg"=mTEC_ZsGn/Total.Cells) %>%
  group_by(Day) %>% 
  summarise("Mean.cTEC_ZsGpos"=mean(p.cTEC_ZsGpos), "Mean.mTEC_ZsGpos"=mean(p.mTEC_ZsGpos),
            "Mean.cTEC_ZsGneg"=mean(p.cTEC_ZsGneg), "Mean.mTEC_ZsGneg"=mean(p.mTEC_ZsGneg)) %>%
  pivot_longer(cols=!Day, values_to="Freq", names_to="Sample")
cell.df$Fraction <- gsub(cell.df$Sample, pattern="(Mean\\.)(\\S+)_(ZsG\\S+)", replacement="\\3")
cell.df$TECSort <- gsub(cell.df$Sample, pattern="(Mean\\.)(\\S+)_(ZsG\\S+)", replacement="\\2")
  
ly6d.freqs <- read.table("/Volumes/Single_Cell_Genetics-CF10180-29/AgeingExperiment/Perinatal_manuscript/data/Ly6d_freqs_TECSort.csv",
                        sep=",", header=TRUE)
ly6d.freqs$Fraction <- ifelse(ly6d.freqs$Fraction == "ZsGp", "ZsGpos", "ZsGneg")
ly6d.freqs$Day <- gsub(ly6d.freqs$Day, pattern="d", replacement="D")
ly6d.df <- ly6d.freqs %>% group_by(Day, Fraction, TECtype) %>% summarise("Mean"=mean(Freq))
ly6d.tecsort.df <- ly6d.freqs %>% group_by(Day, Fraction, TECtype, TECSort) %>% summarise("Mean"=mean(Freq))

ly6d.tecsort.df[ly6d.tecsort.df$Day %in% c("D2") & ly6d.tecsort.df$TECSort %in% c("mTEC"),] %>%
  group_by(TECtype, Fraction) %>% summarise("Mean"=sum(Mean))



# merge.df <- merge(ly6d.freqs, cell.df, by=c("Day", "Fraction"))



cell.sum <- cell.counts %>% group_by(Day, Replicate) %>% summarise("Total.Cells"=sum(Ncells))
frac.sum <- cell.counts %>% mutate(.by=c(Day, Replicate), .keep="all", )