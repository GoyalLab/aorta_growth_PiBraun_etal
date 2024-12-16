####script for PiBraun et al. paper analysis of kNN between inducible rainbow clones
####Started by Yogesh Goyal on Nov 24, 2022 (Last modified: Feb 22, 2023)
####input data provided by Jonas Braun

library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)

##defining the directories
dataDirectory <- "data/inducible rainbow/spatialPattern/kNN/"
#dataDirectory <- " data/inducible rainbow/spatialPattern/kNN/"
plotDirectory <- 'plots/Supplementary/'

##loading the data
allInducible = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNormClusterskNN.csv'), stringsAsFactors=F, header = T, sep = ","))

ageList = allInducible %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "10-60","21-30","30-60","5-10","5-30", "5-60" ))

allInducibleAgeCorrected = inner_join(allInducible, ageList, by = "age") %>% select(-age) %>% rename(age = ageReal, distance = kNNType)

allInducibleAgeCorrected$age =  factor(allInducibleAgeCorrected$age, levels=c("5-10", "5-30", "5-60", "10-21", "10-30", "10-60", "21-30", "30-60"))
allInducibleAgeCorrected$aortaNumber =  factor(allInducibleAgeCorrected$aortaNumber, levels=c("aorta07", "aorta06", "aorta05", "aorta04", "aorta03", "aorta02", "aorta01"))

allInducibleAgeCorrected.summary = allInducibleAgeCorrected %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))

allInducible1NN.summary = allInducibleAgeCorrected %>% filter(distance == "1NN") %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))

plot1 = ggplot(allInducibleAgeCorrected %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_jitter(width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.5, color = "gray58") +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducible1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot1, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1ALL.svg'), width = 6, height = 4)

plot1a = ggplot(allInducibleAgeCorrected %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_jitter(aes(color = aortaNumber), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducible1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot1a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1ALLColored.svg'), width = 6, height = 4)
ggsave(plot1a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1ALLColoredLegend.svg'), width = 6, height = 4)


plot2 = ggplot(allInducibleAgeCorrected, aes(x = age, y = value)) +
  geom_jitter(color = "gray58", width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducibleAgeCorrected.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot2, file = paste0(plotDirectory, 'inducibleRandomBaseline_knnAll.svg'), width = 6, height = 4)


plot2a = ggplot(allInducibleAgeCorrected, aes(x = age, y = value)) +
  geom_jitter(aes(color = distance), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducibleAgeCorrected.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot2a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knnAllColored.svg'), width = 6, height = 4)
ggsave(plot2a, file = paste0(plotDirectory, 'inducibleRandomBaseline_knnAllColoredLegend.svg'), width = 6, height = 4)


###################
### Plotting Red and Orange channels on the same plot
allInducibleRed = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNormClusterskNN.csv'), stringsAsFactors=F, header = T, sep = ","))
allInducibleOrange = as_tibble(read.table(paste0(dataDirectory, 'markedNormalized_kNN/', 'markedNormClustersOrangekNN.csv'), stringsAsFactors=F, header = T, sep = ","))

ageListRed = allInducibleRed %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "10-60","21-30","30-60","5-10","5-30", "5-60" ))
ageListOrange = allInducibleOrange %>% select(age) %>% unique() %>% mutate(ageReal = c("10-21", "10-30", "5-10","5-30" ))

allInducibleAgeCorrectedRed = inner_join(allInducibleRed, ageListRed, by = "age") %>% select(-age) %>% rename(age = ageReal, distance = kNNType) %>% mutate(color = "red") %>% filter(age %in% c("10-21", "10-30", "5-10","5-30" ))
allInducibleAgeCorrectedOrange = inner_join(allInducibleOrange, ageListOrange, by = "age") %>% select(-age) %>% rename(age = ageReal, distance = kNNType) %>% mutate(color = "orange")

allInducibleAgeCorrectedOrangeRed = bind_rows(allInducibleAgeCorrectedRed,allInducibleAgeCorrectedOrange)
allInducibleAgeCorrectedOrangeRed$age =  factor(allInducibleAgeCorrectedOrangeRed$age, levels=c("5-10", "5-30", "10-21", "10-30"))

allInducibleAgeCorrectedOrangeRed.summary = allInducibleAgeCorrectedOrangeRed %>% group_by(age, color) %>% summarise(sd = sd(value), value = mean(value))
allInducibleAgeCorrectedOrangeRed1NN.summary = allInducibleAgeCorrectedOrangeRed %>% filter(distance == "1NN") %>% group_by(age,color) %>% summarise(sd = sd(value), value = mean(value))

plot3 = ggplot(allInducibleAgeCorrectedOrangeRed %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_jitter(aes(color = color), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allInducibleAgeCorrectedOrangeRed1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot3, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1RedOrange.svg'), width = 3, height = 4)


plot4 = ggplot(allInducibleAgeCorrectedOrangeRed %>% filter(distance == "1NN"), aes(x = age, y = value)) +
  geom_point(aes(color = color), shape = 16, size = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd, color = color), size = 1, position = position_dodge(width = 0.5), data = allInducibleAgeCorrectedOrangeRed1NN.summary) +
  ylim(-2.5,2.5) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot4, file = paste0(plotDirectory, 'inducibleRandomBaseline_knn1RedOrangeSideBySide.svg'), width = 5, height = 4)


###########################################################################################################
####################################### checking for division bias ########################################
###########################################################################################################

dataDirectory <- "spatialPattern/kNN/"
plotDirectory <- 'plots/Supplementary/'

dataDirectory <- "data/inducible rainbow/spatialPattern/kNN/"
plotDirectory <- 'plots/'


allDividingInducible = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividingClusterskNN_v3.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingInducible = allDividingInducible %>% mutate(type = "realData")

allDividingPos1 = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividedPositiveControlMidkNN1Run.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingPos2 = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividedPositiveControlOutkNN1Run.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingPos1= allDividingPos1 %>% mutate(type = "PositiveControl1")
allDividingPos2 = allDividingPos2 %>% mutate(type = "PositiveControl2")

allDividingNeg = as_tibble(read.table(paste0(dataDirectory, 'dividingClusterskNN/', 'dividedNegativeControlkNN1Run.csv'), stringsAsFactors=F, header = T, sep = ","))
allDividingNeg = allDividingNeg %>% mutate(type = "NegativeControl")


allDividingInducibleWithControls = bind_rows(allDividingInducible,allDividingPos1,allDividingPos2,allDividingNeg)

allDividingInducible.summary = allDividingInducible %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))
allInducible1NN.summary = allDividingInducible %>% filter(kNNType == "1NN") %>% group_by(age) %>% summarise(sd = sd(value), value = mean(value))

allDividingInducibleWithControls.summary = allDividingInducibleWithControls %>% group_by(type) %>% summarise(sd = sd(value), value = mean(value))

#type                sd    value
#<fct>            <dbl>    <dbl>
#  1 NegativeControl  0.381  0.00273
#2 PositiveControl1 0.205 -0.446  
#3 PositiveControl2 0.290 -0.322  
# 4 realData         0.409  0.00791

allDividingInducibleWithControls$type = factor(allDividingInducibleWithControls$type, levels=c("realData", "NegativeControl","PositiveControl1", "PositiveControl2"))
allDividingInducibleWithControls.summary$type = factor(allDividingInducibleWithControls.summary$type, levels=c("realData", "NegativeControl","PositiveControl1", "PositiveControl2"))

allDividingInducible$age = factor(allDividingInducible$age, levels = c( "0-5", "0-10", "5-10","5-30", "5-60","10-21", "10-30", "10-60","21-30","30-60"))
allDividingInducible.summary$age = factor(allDividingInducible.summary$age, levels = c( "0-5", "0-10", "5-10","5-30", "5-60","10-21", "10-30", "10-60","21-30","30-60"))


plot5 = ggplot(allDividingInducible, aes(x = age, y = value)) +
  geom_jitter(aes(color = aorta), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allDividingInducible.summary) +
  ylim(-4,4) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

ggsave(plot5, file = paste0(plotDirectory, 'inducibleDividingRedAll.svg'), width = 8, height = 5)

set.seed(2059)

allDividingInducibleWithControls1 = sample_n(allDividingInducibleWithControls,25000) ## subsampling some rows otherwise the plots are pretty large size
allDividingInducibleWithControls.summary$type = factor(allDividingInducibleWithControls.summary$type, levels=c("realData", "NegativeControl","PositiveControl1", "PositiveControl2"))

real = 

wilcox.test(allDividingInducible$value, allDividingPos1$value ,
            alternative = "greater", paired = FALSE, conf.int = TRUE, conf.level = 0.95)

wilcox.test(allDividingInducible$value, allDividingPos2$value ,
            alternative = "greater", paired = FALSE, conf.int = TRUE, conf.level = 0.95)

wilcox.test(allDividingInducible$value, allDividingNeg$value ,
            alternative = "greater", paired = FALSE, conf.int = TRUE, conf.level = 0.95)


plot6 = ggplot(allDividingInducibleWithControls1, aes(x = type, y = value)) +
  geom_jitter(aes(color = type), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allDividingInducibleWithControls.summary) +
  ylim(-4,4) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot6, file = paste0(plotDirectory, 'DividingWControls.svg'), width = 5, height = 4)

plot7 = ggplot(allDividingInducibleWithControls, aes(x = type, y = value)) +
  geom_jitter(aes(color = age), width = 0.2, height = 0.1, shape = 16, size = 0.5, alpha = 0.4) +
  geom_hline(yintercept=0, color = "hotpink3", size=1) +
  geom_pointrange(aes(ymin = value-sd, ymax = value+sd), size = 1, color = "black", data = allDividingInducibleWithControls.summary) +
  ylim(-4,4) +
  theme_classic((base_size = 18))+
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot7, file = paste0(plotDirectory, 'DividingWControls_includingAgeColors.svg'), width = 5, height = 4)
