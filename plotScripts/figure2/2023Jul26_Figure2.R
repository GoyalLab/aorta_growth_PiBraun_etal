library(ggplot2)
library(tidyverse)
library(svglite)
library(tibble)
library(dplyr)
library(formattable)
library(reshape2) 

dataDirectory <- "data/figure2/"
dataDirectory <- "data/figure2/"

plotDirectory <- "plots/figure2/"
plotDirectory <- "plots/figure2/"

########################################################################################
###################################### Cell numbers #####################################
########################################################################################


cellNumbers <- tibble(read.table(file = paste0(dataDirectory, "nucNumberWOP50.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

cellNumbers.summary <- cellNumbers %>% group_by(age) %>% dplyr::summarise(sd = sd(numberNuclei),
                                                                          mean = mean(numberNuclei)) 

plot1 <- ggplot(data = cellNumbers.summary, aes(x = age, y= mean)) +  
  geom_jitter(data = cellNumbers, aes(x=age, y=numberNuclei, color = aorta),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  geom_line() + 
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot1, file = paste0(plotDirectory, 'totalCellNumbers.svg'), width = 5, height = 3)


########################################################################################
###################################### major and minor axis numbers #####################################
########################################################################################

shapeNumbers <- tibble(read.table(file = paste0(dataDirectory, "allCellShapesUpdate_MicroCTAspect_NEWFINAL.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

shapeNumbers.summary <- shapeNumbers %>% group_by(age, aorta) %>% dplyr::summarise(sdMinor = sd(minor_axis_length_Micro),
                                                                            sdMajor = sd(major_axis_length_MicroCT),
                                                                            sdAspect = sd(aspectRatio),
                                                                            sdArea = sd(areaAdjusted_MicroCT),
                                                                            meanMinor = mean(minor_axis_length_Micro),
                                                                            meanMajor = mean(major_axis_length_MicroCT),
                                                                            meanAspect = mean(aspectRatio),
                                                                            meanArea = mean(areaAdjusted_MicroCT))

shapeNumbers.summary.summary <- shapeNumbers.summary %>% group_by(age) %>% dplyr::summarise(sdMinor = sd(meanMinor),
                                                                          sdMajor = sd(meanMajor),
                                                                          sdAspect = sd(meanAspect),
                                                                          sdArea = sd(meanArea),
                                                                          meanMinor = mean(meanMinor),
                                                                          meanMajor = mean(meanMajor),
                                                                          meanAspect = mean(meanAspect),
                                                                          meanArea = mean(meanArea))
                                                                         

plot2 <- ggplot(data = shapeNumbers.summary.summary, aes(x = age, y= meanMajor)) +  
  geom_jitter(data = shapeNumbers.summary, aes(x=age, y=meanMajor),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = meanMajor, ymin = meanMajor - sdMajor, ymax = meanMajor + sdMajor), alpha = .2) +
  geom_line() + 
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot2, file = paste0(plotDirectory, 'meanMajorAxis_v1.svg'), width = 5, height = 3)

plot3 <- ggplot(data = shapeNumbers.summary.summary, aes(x = age, y= meanMinor)) +  
  geom_jitter(data = shapeNumbers.summary, aes(x=age, y=meanMinor),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = meanMinor, ymin = meanMinor - sdMinor, ymax = meanMinor + sdMinor), alpha = .2) +
  geom_line() + 
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot3, file = paste0(plotDirectory, 'meanMinorAxis_v1.svg'), width = 5, height = 3)

plot4 <- ggplot(data = shapeNumbers.summary.summary, aes(x = age, y= meanAspect)) +  
  geom_jitter(data = shapeNumbers.summary, aes(x=age, y=meanAspect),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = meanAspect, ymin = meanAspect - sdAspect, ymax = meanAspect + sdAspect), alpha = .2) +
  geom_line() + 
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot4, file = paste0(plotDirectory, 'meanAspect_v1.svg'), width = 5, height = 3)

plot5 <- ggplot(data = shapeNumbers.summary.summary, aes(x = age, y= meanArea)) +  
  geom_jitter(data = shapeNumbers.summary, aes(x=age, y=meanArea),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = meanArea, ymin = meanArea - sdArea, ymax = meanArea + sdArea), alpha = .2) +
  geom_line() + 
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot5, file = paste0(plotDirectory, 'meanArea_v1.svg'), width = 5, height = 3)

#################################################################################################################
###################################### comparisons with global measurements #####################################
#################################################################################################################
data1Directory <- "data/general_microCT/"
data2Directory <- "data/circumference/"

microCT <- tibble(read.table(file = paste0(data1Directory, "totalLengthAortas.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

microCT.summary <- microCT  %>% group_by(age) %>% dplyr::summarise(length = mean(length_1_8),
                                                                   sd = sqrt(var(length_1_8))/sqrt(length(length_1_8)))

circumferenceIF <- tibble(read.table(file = paste0(data2Directory, "230703_circumference.csv"), header = FALSE, stringsAsFactors=F, sep = ",", fill = TRUE))

circumferenceIF = circumferenceIF %>% rename(age = V1,
                                             aortaID = V2,
                                             circumference = V3)

circumferenceIF$age <- gsub("P", "", circumferenceIF$age)
circumferenceIF$age <- as.integer(circumferenceIF$age)

circumferenceIF.summary <- circumferenceIF  %>% group_by(age) %>% dplyr::summarise(sdCircumference = sd(circumference),
                                                                                   circumference = mean(circumference))

normalizedLength = microCT.summary %>% mutate(normalizedLength = (length - min(length))/(max(length) - min(length)))
normalizedCircumference = circumferenceIF.summary %>% mutate(normalizedCircumference = (circumference - min(circumference))/(max(circumference) - min(circumference)))

normalizedMajor = shapeNumbers.summary.summary %>% mutate(normalizedMajor = (meanMajor - min(meanMajor))/(max(meanMajor) - min(meanMajor)))
normalizedMinor = shapeNumbers.summary.summary %>% mutate(normalizedMinor = (meanMinor - min(meanMinor))/(max(meanMinor) - min(meanMinor)))

plot6 <- ggplot() +  
  geom_line(data = normalizedLength, aes(x = age, y= normalizedLength)) +
  geom_line(data = normalizedMajor, aes(x = age, y= normalizedMajor)) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot6, file = paste0(plotDirectory, 'NormalizedMajorLength.svg'), width = 5, height = 3)

plot7 <- ggplot() +  
  geom_line(data = normalizedCircumference, aes(x = age, y= normalizedCircumference), color = "blue") +
  geom_line(data = normalizedMinor, aes(x = age, y= normalizedMinor), color = "blue") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot7, file = paste0(plotDirectory, 'NormalizedMinorCircumference.svg'), width = 5, height = 3)

#################################################################################################################
###################################### comparisons with global measurements V1 #####################################
#################################################################################################################

dataDirectory <- "data/figure2/"

circumferenceAll <- tibble(read.table(file = paste0(dataDirectory, "nucAlongCircumferenceMean.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
lengthAll<- tibble(read.table(file = paste0(dataDirectory, "nucAlongLengthMean.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

circumferenceAll = circumferenceAll %>% mutate(error = 100*(circumference_estimation/circumference_real - 1))
circumferenceAll.summary = circumferenceAll %>% group_by(age) %>% dplyr::summarise(meanCircumference_real = mean(circumference_real),
                                                                                   meanCircumference_estimation = mean(circumference_estimation),
                                                                                   meanError = mean(error),
                                                                                   sdCircumference_real = sd(circumference_real),
                                                                                   sdCircumference_estimation = sd(circumference_estimation),
                                                                                   sdError = sd(error))

lengthAll = lengthAll %>% mutate(error = 100*(length_estimation/length_real - 1))
lengthAll.summary = lengthAll %>% group_by(age) %>% dplyr::summarise(meanLength_real = mean(length_real),
                                                                                   meanLength_estimation = mean(length_estimation),
                                                                                   meanError = mean(error),
                                                                                   sdLength_real = sd(length_real),
                                                                                   sdLength_estimation = sd(length_estimation),
                                                                                   sdError = sd(error))

plot8 <- ggplot(data = circumferenceAll.summary) + 
  geom_line(aes(x = age, y = meanCircumference_real), color = "blue") +
  geom_line(aes(x = age, y = meanCircumference_estimation), color = "gray") +
  geom_pointrange(aes(x = age, y =meanCircumference_real, ymin = meanCircumference_real-sdCircumference_real, ymax = meanCircumference_real+sdCircumference_real), color = "blue") +
  geom_pointrange(aes(x = age, y =meanCircumference_estimation, ymin = meanCircumference_estimation-sdCircumference_estimation, ymax = meanCircumference_estimation+sdCircumference_estimation), color = "gray") +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot8, file = paste0(plotDirectory, 'CircumferenceComparison_V1.svg'), width = 5, height = 3)

plot9 <- ggplot(data = lengthAll.summary) + 
  geom_line(aes(x = age, y = meanLength_real), color = "blue") +
  geom_line(aes(x = age, y = meanLength_estimation), color = "gray") +
  geom_pointrange(aes(x = age, y =meanLength_real, ymin = meanLength_real-sdLength_real, ymax = meanLength_real+sdLength_real), color = "blue") +
  geom_pointrange(aes(x = age, y =meanLength_estimation, ymin = meanLength_estimation-sdLength_estimation, ymax = meanLength_estimation+sdLength_estimation), color = "gray") +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot9, file = paste0(plotDirectory, 'LengthComparison_V1.svg'), width = 5, height = 3)

plot10 <- ggplot(data = lengthAll.summary) + 
  geom_line(aes(x = age, y = meanError), color = "blue") +
  geom_line(aes(x = age, y = meanError), color = "gray", data = circumferenceAll.summary) +
  geom_pointrange(aes(x = age, y =meanError, ymin = meanError-sdError, ymax = meanError+sdError), color = "blue") +
  geom_pointrange(aes(x = age, y =meanError, ymin = meanError-sdError, ymax = meanError+sdError), color = "gray", data = circumferenceAll.summary) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") +
  ylim (-100,100)

ggsave(plot10, file = paste0(plotDirectory, 'ErrorComparison_V1.svg'), width = 5, height = 3)

nucCircumferenceAll <- tibble(read.table(file = paste0(dataDirectory, "nucAlongCircumferenceMean.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
nucLengthAll <- tibble(read.table(file = paste0(dataDirectory, "nucAlongLengthMean.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

nucCircumferenceAll.summary = nucCircumferenceAll %>% group_by(age) %>% dplyr::summarise(meanNucsMean = mean(nucsMean),
                                                                                   meanAxisLength = mean(axisLength),
                                                                                   sdNucsMean = sd(nucsMean),
                                                                                   sdAxisLength = sd(axisLength))

nucLengthAll.summary = nucLengthAll %>% group_by(age) %>% dplyr::summarise(meanNucsMean = mean(nucsMean),
                                                                                         meanAxisLength = mean(axisLength),
                                                                                         sdNucsMean = sd(nucsMean),
                                                                                         sdAxisLength = sd(axisLength))
plot11 <- ggplot(data = nucLengthAll.summary) + 
  geom_line(aes(x = age, y = meanNucsMean), color = "blue") +
  geom_line(aes(x = age, y = meanNucsMean), color = "gray", data = nucCircumferenceAll.summary) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") 

ggsave(plot11, file = paste0(plotDirectory, 'meanCellNumbers_lengthCircumference.svg'), width = 5, height = 3)

plot12 <- ggplot(data = nucLengthAll.summary) + 
  geom_line(aes(x = age, y = meanAxisLength), color = "blue") +
  geom_line(aes(x = age, y = meanAxisLength), color = "gray", data = nucCircumferenceAll.summary) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") 

ggsave(plot12, file = paste0(plotDirectory, 'meanAxisLength_lengthCircumference.svg'), width = 5, height = 3)


########################################################################################
###################################### Cell density #####################################
##########################################################################################

cellDensityLength <- tibble(read.table(file = paste0(dataDirectory, "densitiesLengthLocationMicroCT_FINAL.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

cellDensityLength.summary <- cellDensityLength %>% group_by(age, aorta, location) %>% dplyr::summarise(mean = mean(densities_MicroCT))


cellDensityLength.summary.summary <- cellDensityLength.summary %>% group_by(location,age) %>% dplyr::summarise(sd = sd(mean),
                                                                          mean = mean(mean)) 

cellDensityLength.summary.summary$ageCategory <- NA

cellDensityLength.summary.summary$ageCategory[cellDensityLength.summary.summary$age <8] = "early"
cellDensityLength.summary.summary$ageCategory[cellDensityLength.summary.summary$age >8] = "mid"
cellDensityLength.summary.summary$ageCategory[cellDensityLength.summary.summary$age >20] = "late"

cellDensityLength.summary.summary.summary = cellDensityLength.summary.summary %>% group_by(ageCategory, location) %>% dplyr::summarise(sd = sd(mean),
                                                                                                                                       mean = mean(mean))

plot8 <- ggplot(data = cellDensityLength.summary.summary.summary, aes(x = location, y= mean, fill=as.factor(ageCategory))) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  geom_line() +
  theme_classic() +
  ylim(0,0.005) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot8, file = paste0(plotDirectory, 'DensityLength_v1.svg'), width = 5, height = 3)

###########
cellDensityCircum <- tibble(read.table(file = paste0(dataDirectory, "densitiesCircumferenceLocationMicroCT_FINAL.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

cellDensityCircum.summary <- cellDensityCircum %>% group_by(age, aorta, location) %>% dplyr::summarise(mean = mean(densities_MicroCT))


cellDensityCircum.summary.summary <- cellDensityCircum.summary %>% group_by(location,age) %>% dplyr::summarise(sd = sd(mean),
                                                                                                               mean = mean(mean)) 

cellDensityCircum.summary.summary$ageCategory <- NA

cellDensityCircum.summary.summary$ageCategory[cellDensityCircum.summary.summary$age <8] = "early"
cellDensityCircum.summary.summary$ageCategory[cellDensityCircum.summary.summary$age >8] = "mid"
cellDensityCircum.summary.summary$ageCategory[cellDensityCircum.summary.summary$age >20] = "late"

cellDensityCircum.summary.summary.summary = cellDensityCircum.summary.summary %>% group_by(ageCategory, location) %>% dplyr::summarise(sd = sd(mean),
                                                                                                                                       mean = mean(mean))

plot9 <- ggplot(data = cellDensityCircum.summary.summary.summary, aes(x = location, y= mean, fill=as.factor(ageCategory))) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  geom_line() +
  theme_classic() +
  ylim(0,0.005) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none")

ggsave(plot9, file = paste0(plotDirectory, 'DensityCircumference_v1.svg'), width = 5, height = 3)


########################################################################################
###################################### Nuclei Statistics #####################################
##########################################################################################

dataDirectory <- "data/cellShape/"

nucleiAll <- tibble(read.table(file = paste0(dataDirectory, "statsNuclei_AortaLevel.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

nucleiAllSummary = nucleiAll %>% group_by(age) %>% dplyr::summarise(meanArea = mean(area),
                                                                    sdArea = sd(area),
                                                                    meanMajorAxisLength = mean(major_axis_length),
                                                                    sdMajorAxisLength = sd(major_axis_length),
                                                                    meanMinorAxisLength = mean(minor_axis_length),
                                                                    sdMinorAxisLength = sd(minor_axis_length))

nucleiAllSummary$age = as.numeric(gsub("P", "", nucleiAllSummary$age))


plot11 <- ggplot() + 
  geom_line(aes(x = age, y = meanArea), data = nucleiAllSummary) +
  geom_ribbon(aes(x = age, y = meanArea, ymin = meanArea - sdArea, ymax = meanArea + sdArea), alpha = .2, data = nucleiAllSummary) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") 
ggsave(plot11, file = paste0(plotDirectory, 'nucleiArea.svg'), width = 5, height = 3)

plot12 <- ggplot() + 
  geom_line(aes(x = age, y = meanMajorAxisLength), data = nucleiAllSummary) +
  geom_ribbon(aes(x = age, y = meanMajorAxisLength, ymin = meanMajorAxisLength - sdMajorAxisLength, ymax = meanMajorAxisLength + sdMajorAxisLength), alpha = .2, data = nucleiAllSummary) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") 
ggsave(plot12, file = paste0(plotDirectory, 'nucleiMajorAxis.svg'), width = 5, height = 3)

plot13 <- ggplot() + 
  geom_line(aes(x = age, y = meanMinorAxisLength), data = nucleiAllSummary) +
  geom_ribbon(aes(x = age, y = meanMinorAxisLength, ymin = meanMinorAxisLength - sdMinorAxisLength, ymax = meanMinorAxisLength + sdMinorAxisLength), alpha = .2, data = nucleiAllSummary) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="none") 
ggsave(plot13, file = paste0(plotDirectory, 'nucleiMinorAxis.svg'), width = 5, height = 3)
