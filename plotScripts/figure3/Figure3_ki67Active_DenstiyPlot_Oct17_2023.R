
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library("ggsci")

home1Directory <- "data/ki67/"
plotDirectory <- "plots/figure3/"



df <- tibble(read.table(file = paste0(home1Directory, 'densityki67Circumference.csv'), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

dfMean <- df %>% group_by(age, positionFinal) %>% dplyr::summarise(mean = mean(averageki67Active)) 
dfSD <- df %>% group_by(age, positionFinal) %>% dplyr::summarise(sd = sd(averageki67Active)) 

df_join <- dfMean %>% 
  left_join(dfSD)

  
cbPalett <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 

plot1 <- ggplot(data = df_join, aes(x=positionFinal, y=mean, fill=as.factor(age), color = age)) +
  geom_line(data = df_join, mapping = aes(x=positionFinal, y=mean)) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  theme_bw() + ylab(expression(density)) + xlab(expression(position~circumference)) +
  scale_colour_manual(values=cbPalett, name = "ages") + 
  scale_fill_manual(values=cbPalett, name = "ages") +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
ggsave(plot1, file = paste0(plotDirectory, 'ki67_alongCircumference.svg'), width = 6, height = 3)

####

df <- tibble(read.table(file = paste0(home1Directory, 'densityki67Length.csv'), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

dfMean <- df %>% group_by(age, positionFinal) %>% dplyr::summarise(mean = mean(averageki67Active)) 
dfSD <- df %>% group_by(age, positionFinal) %>% dplyr::summarise(sd = sd(averageki67Active)) 

df_join <- dfMean %>% 
  left_join(dfSD)


cbPalett <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 

plot2 <- ggplot(data = df_join, aes(x=positionFinal, y=mean, fill=as.factor(age), color = age)) +
  geom_line(data = df_join, mapping = aes(x=positionFinal, y=mean)) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = .2) +
  theme_bw() + ylab(expression(density)) + xlab(expression(position~length)) +
  scale_colour_manual(values=cbPalett, name = "ages") + 
  scale_fill_manual(values=cbPalett, name = "ages") +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot2, file = paste0(plotDirectory, 'ki67_alongLength.svg'), width = 6, height = 3)

