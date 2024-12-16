####script for PiBraun et al. paper analysis of SpindleAngles and amounts
####Started by Yogesh Goyal on Nov 24, 2022 (Last modified: Nov 24, 2022)
####input data provided by Jonas Braun

library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)
library(ggdist)

##defining the directories
dataDirectory <- "data/spindleAngles/"
plotDirectory <- 'plots/'

##loading the data
allAngles = as_tibble(read.table(paste0(dataDirectory, 'spindleAnglesFinal_all.tsv'), stringsAsFactors=F, header = T))
allAngles = allAngles %>% mutate(normalizedAngle = 180-angle) %>% filter(age != "19")
allAngles$age = as.character(allAngles$age)
allAngles$age = factor(allAngles$age, levels = c("5", "7","10", "19", "20", "30", "60"))

allAngles.summary = allAngles %>% group_by(age) %>% summarise(sd = sd(normalizedAngle, na.rm = TRUE), normalizedAngle = mean(normalizedAngle, na.rm = TRUE))

plot1 = ggplot(data = allAngles, aes(x = age, y = normalizedAngle, color = replicate)) +
  geom_jitter(width = 0.15, shape =16) +
  geom_hline(yintercept=0, size=1) +
  geom_pointrange(aes(ymin = normalizedAngle-sd, ymax = normalizedAngle+sd), size = 1, color = "black", data = allAngles.summary) +
  theme_classic(base_size = 18)+
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-90, 90), breaks= c(-90, -60, -30, 0, 30, 60, 90))
ggsave(plot1, file = paste0(plotDirectory, 'spindleAnglesAll_Cartesian.svg'), width = 5, height = 3.5)

allAngles7 = allAngles %>% filter(age %in% c("7"))
allAngles7[which(allAngles7$normalizedAngle<0),4] = allAngles7[which(allAngles7$normalizedAngle<0),4] + 360

allAngles20 = allAngles %>% filter(age %in% c("20"))
allAngles20[which(allAngles20$normalizedAngle<0),4] = allAngles20[which(allAngles20$normalizedAngle<0),4] + 360

plot2 = ggplot(allAngles7, aes(x = normalizedAngle)) +
  geom_histogram(binwidth = 5, boundary = 0, fill = "hotpink3") +
  scale_x_continuous(breaks = seq(0, 360, 30)) +
  coord_polar(theta= "x", start = 0) +
  xlab(NULL)+ylab(NULL) +
  theme_bw()
ggsave(plot2, file = paste0(plotDirectory, 'spindleAnglesP7_Radial.svg'), width = 3, height = 3)


plot3 = ggplot(allAngles20, aes(x = normalizedAngle)) +
  geom_histogram(binwidth = 5, boundary = 0, fill = "hotpink3") +
  scale_x_continuous(breaks = seq(0, 360, 30)) +
  coord_polar(theta= "x", start = 0) +
  xlab(NULL)+ylab(NULL) +
  theme_bw()
ggsave(plot3, file = paste0(plotDirectory, 'spindleAnglesP20_Radial.svg'), width = 3, height = 3)

###########################################################################################################
########################################## Spindle Angles plots ###########################################
###########################################################################################################

data2Directory <- "data/spindleAngles/"
data2Directory <- "data/spindleAngles/"

allSpindle <- tibble(read.table(file = paste0(data2Directory, "actmitotSpindleWOP19.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

allSpindle <- allSpindle %>% mutate(percentActive = fractionActive*100)

allSpindle.summary = allSpindle %>% group_by(age) %>% summarise(sd = sd(percentActive, na.rm = TRUE), mean = mean(percentActive, na.rm = TRUE))

plot4 <- ggplot(data = allSpindle.summary, aes(x = age)) +  
  geom_line(aes(y = mean), size = 0.5, color = "black") + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  geom_point(data = allSpindle, mapping = aes(x=age, y=percentActive, color = aorta),size=2, shape = 16) +
  theme_bw() + ylab(expression(active~mitotic~nuclei~(per~mille))) + xlab("age (days)") +  
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) +
  scale_x_continuous(limits = c(0,60.5), breaks = c(0,5,7,10, 20, 30, 60)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot4, file = paste0(plotDirectory, 'spindleFraction.svg'), width = 6, height = 3)



allSpindle1.summary = allSpindle %>% group_by(age) %>% summarise(sd = sd(active, na.rm = TRUE), mean = mean(active, na.rm = TRUE))

plot5 <- ggplot(data = allSpindle1.summary, aes(x = age)) +  
  geom_line(aes(y = mean), size = 0.5, color = "black") + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  geom_point(data = allSpindle, mapping = aes(x=age, y=active, color = aorta),size=2, shape = 16) +
  theme_bw() + ylab(expression(active~mitotic~nuclei~(per~mille))) + xlab("age (days)") +  
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) +
  scale_x_continuous(limits = c(0,60.5), breaks = c(0,5,7,10, 20, 30, 60)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot5, file = paste0(plotDirectory, 'spindleTotal.svg'), width = 6, height = 3)


##### Comparisons for p7 vs p20
allSpindlep5p7p20 = allSpindle %>% filter(age %in% c("5", "7"))

allSpindlep5p7 = allSpindle %>% filter(age %in% c("5", "7"))
allSpindlep7 = allSpindle %>% filter(age == "7")

allSpindlep20 = allSpindle %>% filter(age== "20")

wilcox.test(allSpindlep5p7$active, allSpindlep20$active, alternative = "two.sided")
wilcox.test(allSpindlep5p7$percentActive, allSpindlep20$percentActive, alternative = "two.sided")

plot5 <- ggplot(data = allSpindle1.summary, aes(x = age)) +  
  geom_line(aes(y = mean), size = 0.5, color = "black") + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  geom_point(data = allSpindle, mapping = aes(x=age, y=active, color = aorta),size=2, shape = 16) +
  theme_bw() + ylab(expression(active~mitotic~nuclei~(per~mille))) + xlab("age (days)") +  
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) +
  scale_x_continuous(limits = c(0,60.5), breaks = c(0,5,7,10, 20, 30, 60)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

