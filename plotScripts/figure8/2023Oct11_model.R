library(ggplot2)
library(tidyverse)
library(svglite)
library(tibble)
library(dplyr)
library(formattable)
library(reshape2) 

dataDirectory <- "data/model/from Sayantan/"
dataDirectory <- "data/model/from Sayantan/"
plotDirectory <- "plots/figure7/"


changePruning <- tibble(read.table(file = paste0(dataDirectory, "pruning_vs_nopruning.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
changePruning2 <- tibble(read.table(file = paste0(dataDirectory, "pruningTall.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

changePruning = changePruning %>% mutate(ratio = noPruning/pruning)

changePruning2$type = factor(changePruning2$type, levels=c("pruning", "noPruning"))

plot1 = ggplot(data=changePruning2, aes(x=type, y=value, group=parameterValue)) + 
  geom_point() + geom_line() +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  ylab("Normalized Cell Numbers")

ggsave(plot2, file = paste0(plotDirectory, 'pruning_p19_EdU.svg'), width = 5, height = 4)



#########################
####re-running simulations

newTime <- tibble(read.table(file = paste0(dataDirectory, "T1_data_tall_top15.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

newTime2 <- tibble(read.table(file = paste0(dataDirectory, "T1_data.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

newTime$type = factor(newTime$type, levels=c("pruning", "noPruning"))
newTime2$type = factor(newTime2$type, levels=c("pruning", "noPruning"))

newTime.summary = newTime %>% group_by(type) %>% dplyr::summarise(sdValue = sd(value),
                                                                  meanValue = median(value))

newTime2.summary = newTime2 %>% group_by(type) %>% dplyr::summarise(sdValue = sd(value),
                                                                  meanValue = median(value))

ggplot(data=newTime, aes(x=type, y=value)) + 
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 16) +
  geom_pointrange(aes(y = meanValue, ymin = meanValue - sdValue, ymax = meanValue + sdValue), data = newTime.summary)+
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  ylab("Cell Cycle Time (days)")

ggplot(data=newTime2, aes(x=type, y=value)) + 
  #geom_jitter(width = 0.1, alpha = 0.3) +
  geom_boxplot(width = 0.4) + 
  theme_classic(base_size = 16) +
  #geom_pointrange(aes(y = meanValue, ymin = meanValue - sdValue, ymax = meanValue + sdValue), data = newTime2.summary, color = '#E69F00')+
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  ylab("Cell Cycle Time (days)")
