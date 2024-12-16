
library(ggplot2)
library(dplyr)
library(tidyr)

dataDirectory <- "results/inducedRainbow/"
plotDirectory <- "plots/inducedRainbow/gini/"

dataDirectory <- "data/inducible rainbow/"
plotDirectory <- "plots/inducibleRainbowGini/"

dataDirectory <- "data/inducible rainbow/"
plotDirectory <- "plots/inducibleRainbowGini/"


df <- tibble(read.table(file = paste0(dataDirectory, "clusterSizesStatistics.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

df_depr <- df %>%
  mutate(group = case_when(
    endsWith(age, "10") ~ "Wave1",
    startsWith(age, "P10") ~ "Wave2",
    age == "P21-30" ~ "short",
    age == "P30-60" ~ "short",
    age == "P0-5" ~ "Wave1",
    age == "P5-30" ~ "Wave1_2",
    age == "P5-60" ~ "Wave1_2",
  ))

dfWave1 <- filter(df_depr, age %in%  c("P0-5", "P0-10", "P5-10"))
dfWave2 <- filter(df_depr, age %in%  c("P10-21", "P10-30", "P10-60"))
dfWave1_2 <- filter(df_depr, age %in%  c("P5-30", "P5-60"))

cbPalett <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", 171)


plot <- ggplot(data = df, aes(x=divided, y=gini, colour = age)) +
  geom_point(data = df, mapping = aes(x=divided, y=gini, colour = age)) + 
  theme_bw() + ylab(expression(gini)) + xlab(expression(divided/marked)) +
  ylim(0,0.4) + xlim(0, 0.6) + 
  scale_colour_manual(values=cbPalett, name = "ages") + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'gini.svg'), width = 10, height = 5)

plot <- ggplot(data = df_depr, aes(x=divided, y=gini, colour = group)) +
  geom_point(data = df_depr, mapping = aes(x=divided, y=gini, colour = group)) + 
  geom_smooth(data = dfWave1, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#E69F00", se=FALSE) + 
  geom_smooth(data = dfWave2, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#009E73", se=FALSE) + 
  geom_smooth(data = dfWave1_2, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#56B4E9", se=FALSE)+ 
  theme_bw() + ylab(expression(gini)) + xlab(expression(divided/marked)) +
  scale_colour_manual(values=cbPalett, name = "group") + 
  ylim(0,0.4) + xlim(0, 0.6) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'giniGroupedWithGrayAndBlue.svg'), width = 10, height = 5)

plot <- ggplot(data = df_depr, aes(x=divided, y=gini)) +
  geom_point(data = dfWave1, mapping = aes(x=divided, y=gini), colour = "#E69F00") + 
  geom_point(data = dfWave2, mapping = aes(x=divided, y=gini), colour = "#009E73") +
  geom_smooth(data = dfWave1, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#E69F00", se=FALSE) + 
  geom_smooth(data = dfWave2, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#009E73", se=FALSE) + 
  theme_bw() + ylab(expression(gini)) + xlab(expression(divided/marked)) +
  ylim(0,0.4) + xlim(0, 0.6) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'giniGroupedWave1And2Comparison.svg'), width = 10, height = 5)

plot <- ggplot(data = df_depr, aes(x=divided, y=gini)) +
  geom_point(data = dfWave1, mapping = aes(x=divided, y=gini), colour = "#E69F00") + 
  geom_point(data = dfWave2, mapping = aes(x=divided, y=gini), colour = "#009E73") +
  geom_smooth(data = dfWave1, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#E69F00", se=FALSE) + 
  geom_smooth(data = dfWave2, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#009E73", se=FALSE) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'giniGroupedWave1And2ComparisonZoomed.svg'), width = 10, height = 5)

plot <- ggplot(data = df_depr, aes(x=divided, y=gini)) +
  geom_point(data = dfWave1, mapping = aes(x=divided, y=gini), colour = "#E69F00", shape = 16) + 
  geom_point(data = dfWave2, mapping = aes(x=divided, y=gini), colour = "#009E73", shape = 16) +
  geom_point(data = dfWave1_2, mapping = aes(x=divided, y=gini), colour = "#56B4E9", shape = 16) +
  geom_smooth(data = dfWave1, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#E69F00", se=FALSE) + 
  geom_smooth(data = dfWave2, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#009E73", se=FALSE) +
  geom_smooth(data = dfWave1_2, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#56B4E9", se=FALSE)+ 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'giniGroupedWave1And2BothComparisonZoomed.svg'), width = 10, height = 5)

dfWave2 <- filter(df_depr, age %in%  c("P10-21", "P10-30", "P10-60", "P21-30", "P30-60"))

plot <- ggplot(data = df_depr, aes(x=divided, y=gini)) +
  geom_point(data = dfWave1, mapping = aes(x=divided, y=gini), colour = "#E69F00") + 
  geom_point(data = dfWave2, mapping = aes(x=divided, y=gini), colour = "#009E73") +
  geom_smooth(data = dfWave1, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#E69F00", se=FALSE) + 
  geom_smooth(data = dfWave2, mapping = aes(x=divided, y=gini), method=lm, size = 1, color = "#009E73", se=FALSE) + 
  theme_bw() + ylab(expression(gini)) + xlab(expression(divided/marked)) +
  ylim(0,0.4) + xlim(0, 0.6) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'giniGroupedWave1And2ComparisonExtend.svg'), width = 10, height = 5)





