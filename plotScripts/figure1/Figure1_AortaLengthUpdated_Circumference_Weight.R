library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2) 

dataDirectory <- "data/general_microCT/"
data2Directory <- "data/circumference/"
data3Directory <- "data/weightBackbone/"
data4Directory <- "data/general_microCT/"

dataDirectory <- "data/general_microCT/"
data2Directory <- "data/circumference/"
data3Directory <- "data/weightBackbone/"
data4Directory <- "data/general_microCT/"

plotDirectory <- "plots/"
plot2Directory <- "plots/figure1/"

plotDirectory <- "plots/"
plot2Directory <- "plots/figure1/"

########################################################################################
######################################Length plots #####################################
########################################################################################

microCT <- tibble(read.table(file = paste0(dataDirectory, "totalLengthAortas.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

microCT.summary <- microCT  %>% group_by(age) %>% dplyr::summarise(length = mean(length_1_8),
                                                                     sd = sqrt(var(length_1_8))/sqrt(length(length_1_8)))

micro3_7 <- filter(microCT, age %in% c(3,5,7))
micro7_15 <- filter(microCT, age %in% c(7,10,14))
micro7_21 <- filter(microCT, age %in% c(7,10,14, 21))
micro7_30 <- filter(microCT, age %in% c(7,10,14,21,30))

#get slope of linear model
slope3_7 = coef(lm(micro3_7$length_1_8~ micro3_7$age))[2]
slope7_30 = coef(lm(micro7_30$length_1_8~ micro7_30$age))[2]
slope7_15 = coef(lm(micro7_15$length_1_8~ micro7_15$age))[2]
slope7_21 = coef(lm(micro7_21$length_1_8~ micro7_21$age))[2]

plot1 <- ggplot(data = microCT.summary, aes(x = age, y= length)) +  
  geom_jitter(data = microCT, aes(x=age, y=length_1_8, color = sex),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = length, ymin = length - sd, ymax = length + sd), alpha = .2) +
  geom_line() + 
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot1, file = paste0(plotDirectory, 'microCT_length1-8.svg'), width = 5, height = 3)


plot2 <- ggplot(data = microCT.summary, aes(x = age, y= length)) +  
  geom_jitter(data = microCT, aes(x=age, y=length_1_8),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  stat_smooth(data = micro3_7 ,aes(x=age, y=length_1_8), method = "lm", se=FALSE) +
  stat_smooth(data = micro7_15 ,aes(x=age, y=length_1_8), method = "lm", se=FALSE) +
  theme_classic(base_size = 14) +
  xlim(0,15) +
  ylim(5,8.5) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot2, file = paste0(plotDirectory, 'microCT_length1-8_inset_3_7_7_15_slope_30.svg'), width = 4, height = 3)

plot3 <- ggplot(data = microCT.summary, aes(x = age, y= length)) +  
  geom_jitter(data = microCT, aes(x=age, y=length_1_8),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  stat_smooth(data = micro3_7 ,aes(x=age, y=length_1_8), method = "lm", se=FALSE) +
  stat_smooth(data = micro7_21 ,aes(x=age, y=length_1_8), method = "lm", se=FALSE) +
  theme_classic(base_size = 14) +
  xlim(0,22) +
  ylim(5,10.5) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot3, file = paste0(plotDirectory, 'microCT_length1-8_inset_3_7_7_21_slope_20.svg'), width = 4, height = 3)

####################################Individual microCT


individualMicroCT = tibble(read.table(file = paste0(dataDirectory, "intercostalLengths1_8.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))


individualMicroCT.summary = individualMicroCT %>% group_by(age,intercostal)%>% dplyr::summarise(sd = sqrt(var(length))/sqrt(length(length)),
                                                                                                length = mean(length))

AverageIndividualMicroCT.summary = individualMicroCT %>% group_by(age) %>% dplyr::summarise(sd = sqrt(var(length))/sqrt(length(length)),
                                                                                                 length = mean(length))

plot1 <- ggplot(data = individualMicroCT.summary, aes(x = age, y= length)) +  
  geom_line(aes(color = intercostal)) +
  geom_ribbon(aes(y = length, ymin = length - sd, ymax = length + sd, color = intercostal), alpha = .2) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plotX, file = paste0(plotDirectory, 'microCT_length1-8_Individual.svg'), width = 5, height = 3)



############################################ beyond 1-8 intercostal

microCTlength_lsc_ct.summary <- microCT  %>% group_by(age) %>% dplyr::summarise(length = mean(length_lsc_ct),
                                                                                sd = sqrt(var(length_lsc_ct))/sqrt(length(length_lsc_ct)))
plotX <- ggplot(data = microCTlength_lsc_ct.summary, aes(x = age, y= length)) +  
  geom_jitter(data = microCT, aes(x=age, y=length_lsc_ct, color = sex),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = length, ymin = length - sd, ymax = length + sd), alpha = .2) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plotX, file = paste0(plotDirectory, 'microCT_lengthSubclavianCeliacTrunk.svg'), width = 5, height = 3)

########################################################################################
###################################Circumference plots #################################
########################################################################################

circumferenceIF <- tibble(read.table(file = paste0(data2Directory, "230703_circumference.csv"), header = FALSE, stringsAsFactors=F, sep = ",", fill = TRUE))

circumferenceIF = circumferenceIF %>% rename(age = V1,
                                             aortaID = V2,
                                             circumference = V3)

circumferenceIF$age <- gsub("P", "", circumferenceIF$age)
circumferenceIF$age <- as.integer(circumferenceIF$age)

circumferenceIF.summary <- circumferenceIF  %>% group_by(age) %>% dplyr::summarise(sdCircumference = sd(circumference),
                                                                                   circumference = mean(circumference))
circum3_7 <- circumferenceIF %>% filter(age %in% c(3,5,7)) 
circum7_15 <- circumferenceIF %>% filter(age %in% c(7,10,14))
circum3_30 <- circumferenceIF %>% filter(age %in% c(3,5,7,10,14,21,30))

#get slope of linear model
slope3_7 = coef(lm(circum3_7$circumference~ circum3_7$age))[2]
slope7_15 = coef(lm(circum7_15$circumference~ circum7_15$age))[2]
slope3_30 = coef(lm(circum3_30$circumference~ circum3_30$age))[2]


plot1 <- ggplot(data = circumferenceIF.summary, aes(x = age, y= circumference)) +  
  geom_jitter(data = circumferenceIF, aes(x=age, y=circumference),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  geom_ribbon(aes(y = circumference, ymin = circumference - sdCircumference, ymax = circumference + sdCircumference), alpha = .2) +
  geom_line() + 
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  ylim(500,2000)

ggsave(plot1, file = paste0(plotDirectory, 'circumferenceIF.svg'), width = 5, height = 3)

plot2 <- ggplot(data = circumferenceIF.summary, aes(x = age, y= circumference)) +  
  geom_jitter(data = circumferenceIF, aes(x=age, y=circumference),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  stat_smooth(data = circum3_7 ,aes(x=age, y=circumference), method = "lm", se=FALSE) +
  stat_smooth(data = circum7_15 ,aes(x=age, y=circumference), method = "lm", se=FALSE) +
  theme_classic(base_size = 14) +
  xlim(0,15) +
  ylim(500,1500) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot2, file = paste0(plotDirectory, 'circumferenceIF_inset_3_7_7_15_slope_.svg'), width = 4, height = 3)

########################################################################################
######################################Weight plot #####################################
########################################################################################

animalWeight <- tibble(read.table(file = paste0(data3Directory, "animalWeights.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE)) 
animalWeight = animalWeight %>% mutate(testColumn = "mouseWeight")

animalWeightMelt = melt(animalWeight, id.vars  = 25) %>% select(-testColumn) %>% rename(age = variable,
                                                                                        weight = value) 
animalWeightMelt <- na.omit(animalWeightMelt)
animalWeightMelt$age =  gsub("X", "", animalWeightMelt$age)
animalWeightMelt$age <- as.integer(animalWeightMelt$age)

animalWeightMelt.summary = animalWeightMelt %>% group_by(age) %>% dplyr::summarise(sdWeight = sd(weight),
                                                                                   weight = mean(weight))

weight2_8 <- animalWeightMelt %>% filter(age %in% c(2:8)) 
weight8_16 <- animalWeightMelt %>% filter(age %in% c(8:16))

#get slope of linear model
slope2_8 = coef(lm(weight2_8$weight~ weight2_8$age))[2]
slope8_16 = coef(lm(weight8_16$weight~ weight8_16$age))[2]


plot1 <- ggplot(data = animalWeightMelt.summary, aes(x = age, y= weight)) +  
  geom_jitter(data = animalWeightMelt, aes(x=age, y=weight),size=2, shape = 16, width = 0.3, height = 0.3)+ 
  geom_ribbon(aes(y = weight, ymin = weight - sdWeight, ymax = weight + sdWeight), alpha = .2) +
  geom_line() +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot1, file = paste0(plotDirectory, 'animalWeight.svg'), width = 4, height = 3)

plot2 <- ggplot(data = animalWeightMelt.summary, aes(x = age, y= weight)) +  
  geom_jitter(data = animalWeightMelt, aes(x=age, y=weight),size=2, shape = 16, width = 0.3, height = 0.3)+ 
  stat_smooth(data = weight2_8 ,aes(x=age, y=weight), method = "lm", se=FALSE) +
  stat_smooth(data = weight8_16 ,aes(x=age, y=weight), method = "lm", se=FALSE) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  xlim(1,17) +
  ylim(1,8)

ggsave(plot2, file = paste0(plotDirectory, 'animalWeightMelt_inset_2_8_8_16_slope_.svg'), width = 4, height = 3)

########################################################################################################
##########################################Vertebrate Length ############################################
########################################################################################################

vertebrateLength <- tibble(read.table(file = paste0(data4Directory, "vertibrateLength.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE)) 

vertebrateLength.summary = vertebrateLength %>% group_by(age) %>% dplyr::summarise(sdVertebrateLength = sd(vertibrateLength_mm),
                                                                                   meanVertebrateLength = mean(vertibrateLength_mm))

vert3_7 <- filter(vertebrateLength, age %in% c(3,5,7))
vert7_15 <- filter(vertebrateLength, age %in% c(7,10,14))

slope3_7 = coef(lm(vert3_7$vertibrateLength_mm~ vert3_7$age))[2]
slope7_15 = coef(lm(vert7_15$vertibrateLength_mm~ vert7_15$age))[2]

plot1 <- ggplot(data = vertebrateLength.summary, aes(x = age, y= meanVertebrateLength)) +  
  geom_jitter(data = vertebrateLength, aes(x=age, y=vertibrateLength_mm),size=2, shape = 16, width = 0.3, height = 0.3)+ 
  geom_ribbon(aes(y = meanVertebrateLength, ymin = meanVertebrateLength - sdVertebrateLength, ymax = meanVertebrateLength + sdVertebrateLength), alpha = .2) +
  geom_line() +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
ggsave(plot1, file = paste0(plot2Directory, 'vertebrateLength.svg'), width = 4, height = 3)

plot2 <- ggplot() +
  geom_jitter(data = vertebrateLength, aes(x=age, y=vertibrateLength_mm),size=2, shape = 16, width = 0.3, height = 0.0)+ 
  stat_smooth(data = vert3_7 ,aes(x=age, y=vertibrateLength_mm), method = "lm", se=FALSE) +
  stat_smooth(data = vert7_15 ,aes(x=age, y=vertibrateLength_mm), method = "lm", se=FALSE) +
  theme_classic(base_size = 14) +
  xlim(0,15) +
  ylim(9.5,17.5) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot2, file = paste0(plot2Directory, 'vertebrateLength_inset_3_7_7_15.svg'), width = 4, height = 3)

########################################################################################################
######################################Plotting All on the same plot #####################################
########################################################################################################

normalizedLength = microCT.summary %>% mutate(normalizedLength = length/min(length))
normalizedCircumference = circumferenceIF.summary %>% mutate(normalizedCircumference = circumference/min(circumference))
normalizedWeight = animalWeightMelt.summary %>% filter(age != 0) %>% mutate(normalizedWeight = weight/min(weight))
normalizedVertebrate = vertebrateLength.summary %>% mutate(normalizedVertebrate = meanVertebrateLength/min(meanVertebrateLength))

plot1 <- ggplot() +  
  geom_line(data = normalizedWeight, aes(x = age, y= normalizedWeight), color = "hotpink3") +
  geom_line(data = normalizedCircumference, aes(x = age, y= normalizedCircumference), color = "blue") +
  geom_line(data = normalizedLength, aes(x = age, y= normalizedLength)) +
  geom_line(data = normalizedVertebrate, aes(x = age, y= normalizedVertebrate), color = "red") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())


plot1a <- ggplot() +  
  geom_line(data = normalizedCircumference, aes(x = age, y= normalizedCircumference), color = "blue") +
  geom_line(data = normalizedLength, aes(x = age, y= normalizedLength)) +
  geom_line(data = normalizedVertebrate, aes(x = age, y= normalizedVertebrate), color = "red") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot1, file = paste0(plotDirectory, 'fig1AllCombined_v1.svg'), width = 4, height = 3)
ggsave(plot1, file = paste0(plot2Directory, 'fig1AllCombined_v1.svg'), width = 4, height = 3)
ggsave(plot1a, file = paste0(plot2Directory, 'fig1AllCombinedNoWeight_v1.svg'), width = 4, height = 3)
####pinching both sides

normalizedLength = microCT.summary %>% mutate(normalizedLength = (length - min(length))/(max(length) - min(length)))
normalizedCircumference = circumferenceIF.summary %>% mutate(normalizedCircumference = (circumference - min(circumference))/(max(circumference) - min(circumference)))
normalizedWeight = animalWeightMelt.summary %>% filter(age != 0) %>% mutate(normalizedWeight = (weight - min(weight))/(max(weight) - min(weight)))
normalizedVertebrate = vertebrateLength.summary %>% mutate(normalizedVertebrate = (meanVertebrateLength- min(meanVertebrateLength))/(max(meanVertebrateLength) - min(meanVertebrateLength)))

plot2 <- ggplot() +  
  geom_line(data = normalizedWeight, aes(x = age, y= normalizedWeight), color = "hotpink3") +
  geom_line(data = normalizedCircumference, aes(x = age, y= normalizedCircumference), color = "blue") +
  geom_line(data = normalizedLength, aes(x = age, y= normalizedLength)) +
  geom_line(data = normalizedVertebrate, aes(x = age, y= normalizedVertebrate), color = "red") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot2, file = paste0(plotDirectory, 'fig1AllCombined_pinchTopBottom_v1.svg'), width = 4, height = 3)
ggsave(plot2, file = paste0(plot2Directory, 'fig1AllCombined_pinchTopBottom_v1.svg'), width = 4, height = 3)

