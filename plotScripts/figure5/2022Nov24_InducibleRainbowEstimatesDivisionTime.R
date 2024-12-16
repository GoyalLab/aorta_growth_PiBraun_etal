####script for PiBraun et al. paper analysis of inducible rainbow clones
####Started by Yogesh Goyal on Nov 27, 2022 (Last modified: Nov 27, 2022)
####input data provided by Jonas Braun

library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)

##defining the directories
dataDirectory <- "data/inducible rainbow/"
dataDirectory <- "data/inducible rainbow/"
plotDirectory <- 'plots/'
plotDirectory <- "plots/"

##loading the data
allInducibleClonesRed = as_tibble(read.table(paste0(dataDirectory, 'allClusterSizesRed.csv'), stringsAsFactors=F, header = T, sep = ",")) %>% filter(age == "5-10")
allInducibleClonesOrange = as_tibble(read.table(paste0(dataDirectory, 'clusterSizesOrangeRevised.csv'), stringsAsFactors=F, header = T, sep = ","))
ageListOrange = allInducibleClonesOrange %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "10-60", "30-60", "5-10","5-30", "5-60" ))
allInducibleClonesOrange = inner_join(allInducibleClonesOrange, ageListOrange, by = "age") %>% select(-age) %>% rename(age = ageReal) %>% filter(age == "5-10") 

duration = 10-5
allInducibleClonesRed1 = allInducibleClonesRed %>% gather(key = "cloneSize", value = "frequency", 3:15)
allInducibleClonesRed1$cloneSize = as.integer(gsub('X', '', allInducibleClonesRed1$cloneSize))
allInducibleClonesRed2 = allInducibleClonesRed1 %>% group_by(age, aorta) %>% mutate(cumulativeCounts = cloneSize*frequency) %>% filter(cloneSize >3) %>% mutate(tag = "red")

allInducibleClonesOrange1 = allInducibleClonesOrange %>% gather(key = "cloneSize", value = "frequency", 2:11)
allInducibleClonesOrange1$cloneSize = as.integer(gsub('X', '', allInducibleClonesOrange1$cloneSize))
allInducibleClonesOrange2 = allInducibleClonesOrange1 %>% group_by(age, aorta) %>% mutate(cumulativeCounts = cloneSize*frequency) %>% filter(cloneSize > 3) %>% mutate(tag = "orange")

allInducibleCloneSummary = bind_rows(allInducibleClonesRed2, allInducibleClonesOrange2) %>% mutate(cellDivisionTime_hours_Cal1 = duration/log(cloneSize,2), cellDivisionTime_hours_Cal2 = duration/(cloneSize-1))
allInducibleCloneSummary$tag = factor(allInducibleCloneSummary$tag, levels = c("red", "orange"))
allInducibleSummary.summary = allInducibleCloneSummary %>% group_by(aorta, tag) %>% summarise(meanCellDivision_Cal1 = sum(24*frequency*cellDivisionTime_hours_Cal1)/sum(frequency), meanCellDivision_Cal2 = sum(24*frequency*cellDivisionTime_hours_Cal2)/sum(frequency))
allInducibleSummary.summary = allInducibleSummary.summary %>% gather(key = "divisionType", value = "value", 3:4)

allInducibleSummary.summary.summary = allInducibleSummary.summary %>% group_by(tag,divisionType) %>% summarise(sdValue = sd(value), value = mean(value))

plot1 = ggplot(allInducibleSummary.summary, aes(x = divisionType, y = value)) +
  geom_point(aes(color = tag), shape = 16, size = 1.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  geom_pointrange(aes(ymin = value-sdValue, ymax = value+sdValue, color = tag), size = 1, position = position_dodge(width = 0.5), data = allInducibleSummary.summary.summary) +
  theme_classic((base_size = 18))+
  ylim(0,60) +
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot1, file = paste0(plotDirectory, 'cellCycleTime.svg'), width = 2, height = 3)

