####script for PiBraun et al. paper analysis of inducible rainbow clones
####Started by Yogesh Goyal on Nov 24, 2022 (Last modified: Nov 27, 2022)
####input data provided by Jonas Braun

library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)

##defining the directories
dataDirectory <- "data/inducible rainbow/"
#dataDirectory <- "data/inducible rainbow/"
plotDirectory <- 'plots/'

##loading the data
allInducibleClonesRed = as_tibble(read.table(paste0(dataDirectory, 'allClusterSizesRed.csv'), stringsAsFactors=F, header = T, sep = ","))
allInducibleClonesOrange = as_tibble(read.table(paste0(dataDirectory, 'clusterSizesOrangeRevised.csv'), stringsAsFactors=F, header = T, sep = ","))

ageListOrange = allInducibleClonesOrange %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "10-60", "30-60", "5-10","5-30", "5-60" ))
allInducibleClonesOrange = inner_join(allInducibleClonesOrange, ageListOrange, by = "age") %>% select(-age) %>% rename(age = ageReal) %>% filter(age != "5-60") 


allInducibleClonesOrange = allInducibleClonesOrange %>% rowwise() %>% mutate(totalTags = sum(c_across(X1:X10)), fractionDividing = sum(c_across(X2:X10))/sum(c_across(X1:X10)), tag = "orange")
allInducibleClonesRed = allInducibleClonesRed %>% rowwise() %>% mutate(totalTags = sum(c_across(X1:X13)), fractionDividing = sum(c_across(X2:X13))/sum(c_across(X1:X13)), tag = "red")

allInducibleSummary = bind_rows(allInducibleClonesRed %>% select(age, totalTags, fractionDividing, tag) %>% filter(age %in% c("10-21", "10-30", "10-60", "30-60", "5-10","5-30")), allInducibleClonesOrange %>% select(age, totalTags, fractionDividing, tag))
allInducibleSummary$age =  factor(allInducibleSummary$age, levels=c("5-10", "5-30", "10-21", "10-30", "10-60", "30-60"))
allInducibleSummary$tag =  factor(allInducibleSummary$tag, levels=c("red", "orange"))

allInducibleSummary.summary = allInducibleSummary %>% group_by(age,tag) %>% summarise(sdFractionDividing = sd(fractionDividing), fractionDividing = mean(fractionDividing), sdTotalTags = sd(totalTags), totalTags = mean(totalTags))
allInducibleSummary.summary$tag =  factor(allInducibleSummary.summary$tag, levels=c("red", "orange"))

plot1 = ggplot(allInducibleSummary, aes(x = age, y = totalTags)) +
  geom_point(aes(color = tag), shape = 16, size = 1.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  geom_pointrange(aes(ymin = totalTags-sdTotalTags, ymax = totalTags+sdTotalTags, color = tag), size = 1, position = position_dodge(width = 0.5), data = allInducibleSummary.summary) +
  theme_classic((base_size = 18))+
  ylim(0,900) +
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot1, file = paste0(plotDirectory, 'inducibleTotalTags.svg'), width = 6, height = 4)

plot2 = ggplot(allInducibleSummary, aes(x = age, y = fractionDividing)) +
  geom_point(aes(color = tag), shape = 16, size = 1.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  geom_pointrange(aes(ymin = fractionDividing-sdFractionDividing, ymax = fractionDividing+sdFractionDividing, color = tag), size = 1, position = position_dodge(width = 0.5), data = allInducibleSummary.summary) +
  theme_classic((base_size = 18))+
  ylim(0,0.6) +
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot2, file = paste0(plotDirectory, 'inducibleFractionDividing.svg'), width = 6, height = 4)

plot2a = ggplot(allInducibleSummary, aes(x = age, y = fractionDividing)) +
  geom_point(aes(color = tag), shape = 16, size = 1.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  geom_pointrange(aes(ymin = fractionDividing-sdFractionDividing, ymax = fractionDividing+sdFractionDividing, color = tag), size = 1, position = position_dodge(width = 0.5), data = allInducibleSummary.summary) +
  geom_path(aes(x = age, y = fractionDividing, group=tag, color = tag), size = 1, position = position_dodge(width = 0.5), data = allInducibleSummary.summary) +
  theme_classic((base_size = 18))+
  ylim(0,0.6) +
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot2a, file = paste0(plotDirectory, 'inducibleFractionDividingConnected.svg'), width = 6, height = 4)


#####
allInducibleClonesRed1 = allInducibleClonesRed %>% gather(key = "cloneSize", value = "frequency", 3:15)
allInducibleClonesRed1$cloneSize = as.integer(gsub('X', '', allInducibleClonesRed1$cloneSize))
allInducibleClonesRed2 = allInducibleClonesRed1 %>% group_by(age, aorta) %>% mutate(cumulativeCounts = cloneSize*frequency)

allInducibleClonesOrange1 = allInducibleClonesOrange %>% gather(key = "cloneSize", value = "frequency", 2:11)
allInducibleClonesOrange1$cloneSize = as.integer(gsub('X', '', allInducibleClonesOrange1$cloneSize))
allInducibleClonesOrange2 = allInducibleClonesOrange1 %>% group_by(age, aorta) %>% mutate(cumulativeCounts = cloneSize*frequency)

cloneDataSummaryRed = allInducibleClonesRed2 %>% group_by(age,aorta) %>% summarise(meanCloneSize = sum(cumulativeCounts)/sum(frequency)) %>% mutate(tag = "red")
cloneDataSummaryOrange = allInducibleClonesOrange2 %>% group_by(age,aorta) %>% summarise(meanCloneSize = sum(cumulativeCounts)/sum(frequency)) %>% mutate(tag = "orange")

allInducibleCloneSummary = bind_rows(cloneDataSummaryRed %>% select(age, meanCloneSize, tag) %>% filter(age %in% c("10-21", "10-30", "10-60", "30-60", "5-10","5-30")), cloneDataSummaryOrange %>% select(age, meanCloneSize,tag))
allInducibleCloneSummary$age = factor(allInducibleCloneSummary$age, levels = c("5-10", "5-30", "10-21", "10-30", "10-60", "30-60"))
allInducibleCloneSummary$tag =  factor(allInducibleCloneSummary$tag, levels=c("red", "orange"))

allInducibleSummary.summary = allInducibleCloneSummary %>% group_by(age,tag) %>% summarise(sdMeanCloneSize = sd(meanCloneSize), meanCloneSize = mean(meanCloneSize))

plot3 = ggplot(allInducibleCloneSummary, aes(x = age, y = meanCloneSize)) +
  geom_point(aes(color = tag), shape = 16, size = 1.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  geom_path(aes(x = age, y = meanCloneSize, group=tag, color = tag), size = 1, position = position_dodge(width = 0.5), data = allInducibleSummary.summary) +
  geom_pointrange(aes(ymin = meanCloneSize-sdMeanCloneSize, ymax = meanCloneSize+sdMeanCloneSize, color = tag), size = 1, position = position_dodge(width = 0.5), data = allInducibleSummary.summary) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot3, file = paste0(plotDirectory, 'inducibleMeanCloneSize.svg'), width = 6, height = 4)

####Getting Distributions

allInducibleClonesRedSelectAges = allInducibleClonesRed1 %>% filter(age %in% c("5-30", "30-60")) %>% group_by(age,cloneSize) %>% summarise(cumulativeFrequency = sum(frequency)) %>% mutate(tag = "red", normalizedFrequency = 100*cumulativeFrequency/sum(cumulativeFrequency))
allInducibleClonesOrangeSelectAges = allInducibleClonesOrange1 %>% filter(age %in% c("5-30", "30-60")) %>% group_by(age,cloneSize) %>% summarise(cumulativeFrequency = sum(frequency)) %>% mutate(tag = "orange", normalizedFrequency = 100*cumulativeFrequency/sum(cumulativeFrequency))

allInducibleClonesSelectAgesAll = bind_rows(allInducibleClonesRedSelectAges, allInducibleClonesOrangeSelectAges)
allInducibleClonesSelectAgesAll$tag = factor(allInducibleClonesSelectAgesAll$tag, levels = c("red","orange"))
allInducibleClonesSelectAgesAllSummary = allInducibleClonesSelectAgesAll %>% group_by(age, tag) %>% summarise(totalClones = sum(cumulativeFrequency))

plot4 = ggplot(allInducibleClonesSelectAgesAll, aes(x = cloneSize, y = normalizedFrequency)) +
  geom_bar(aes(fill = tag), stat="identity") +
  facet_wrap(~tag + age) +
  xlim(0,14) +
  ylim(0,100) +
  theme_classic((base_size = 18)) +
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot4, file = paste0(plotDirectory, 'distributionsAges.svg'), width = 6, height = 4.5)

###########################################################################################################
####################################### Analysis of the full dataset ######################################
###########################################################################################################



