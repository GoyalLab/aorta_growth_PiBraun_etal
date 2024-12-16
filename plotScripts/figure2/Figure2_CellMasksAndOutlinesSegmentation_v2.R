library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library("ggsci")

#plot segmentation masks
dataDirectory <- "data/cellSegmentationVisualization/"

dataDirectory <- "data/cellSegmentationVisualization/"

plotDirectory <- "images/cellSize/Segmented/"
plotDirectory <- "plots/figure2/"

plotDirectory <- "plots/figure2/"

fullMasks <- tibble(read.table(file = paste0(dataDirectory, "segP60Masks.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

#df <- tibble(read.table(file = path, header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
data_tablePlot <- subset(fullMasks, value != 0)

unique_values <- unique(data_tablePlot$value)

#mapping from unique values to unqiue colors (here: 51)
mapping <- data.table(
  from_values = unique_values,
  to_values = sample(unique_values,length(unique_values) , replace = FALSE)[1:51] 
)

dataPlot <- merge(data_tablePlot, mapping, by.x = "value", by.y = "from_values", all.x = TRUE)

plot1 <- ggplot(data = dataPlot , aes(x=x, y=y_new, color=as.factor(to_values)))  + geom_point(size = 0.05)+
  scale_color_igv() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),  
        axis.text.x = element_blank(),   
        axis.text.y = element_blank(),    
        panel.border=element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none") 

ggsave(plot1, file = paste0(plotDirectory, 'P60Fill.png'), width = 6, height = 6)

#plot outlines with single cell in same color

fullMasksSC = tibble(read.table(file = paste0(dataDirectory, "segP60Masks_SCImages.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
sc <- subset(fullMasksSC, valueSC != 0)
numberSC <- unique(sc$value)[1]
df_singleCell <- dataPlot %>% filter(value == numberSC)

fullMasks <- tibble(read.table(file = paste0(dataDirectory, "segP60Outlines.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
data_tablePlot <- subset(fullMasks, value != 0)

data_tablePlot$to_values <- 0
df_singleCell$to_values <- 1
df_combined <- rbind(data_tablePlot, df_singleCell)

plot <- ggplot(data = df_combined, aes(x=x, y=y_new))  + geom_point(size = 0.4, color= "hotpink3")+
  theme_classic() +
  theme(axis.title.x = element_blank(),   
        axis.title.y = element_blank(),  
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border=element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none") 

ggsave(plot, file = paste0(plotDirectory, 'P60Outlines_fillingSame.png'), width = 6, height = 6)

#plot outlines with single cell in black

plot <- ggplot(data = df_combined, aes(x=x, y=y_new, color = as.factor(to_values)))  + geom_point(size = 0.4)+
  scale_color_manual(values = c('0' = 'hotpink3', '1' = 'black')) + 
  theme_classic() +
  theme(axis.title.x = element_blank(),   
        axis.title.y = element_blank(),  
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border=element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none") 

ggsave(plot, file = paste0(plotDirectory, 'P60Outlines_fillingBlack.png'), width = 6, height = 6)

#plot outlines no filling

fullMasks <- tibble(read.table(file = paste0(dataDirectory, "segP60Outlines.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
data_tablePlot <- subset(fullMasks, value != 0)

plot <- ggplot(data = data_tablePlot, aes(x=x, y=y_new))  + geom_point(size = 0.4, color= "hotpink3")+
  theme_classic() +
  theme(axis.title.x = element_blank(),   
        axis.title.y = element_blank(),  
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border=element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none") 

ggsave(plot, file = paste0(plotDirectory, 'P60Outlines.png'), width = 6, height = 6)



###
fullMasksSC = tibble(read.table(file = paste0(dataDirectory, "segP10Masks_SCImages.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
data_tablePlot <- subset(fullMasksSC, value != 0)

unique_values <- unique(data_tablePlot$value)

#mapping from unique values to unqiue colors (here: 51)
mapping <- data.table(
  from_values = unique_values,
  to_values = sample(unique_values,length(unique_values) , replace = FALSE)[1:51] 
)

dataPlot <- merge(data_tablePlot, mapping, by.x = "value", by.y = "from_values", all.x = TRUE)

dataPlotSC <- subset(dataPlot, valueSC != 0)

plot2 <- ggplot(data = dataPlotSC , aes(x=x, y=y_new, color=as.factor(to_values)))  + geom_point(size = 0.05)+
  scale_color_igv() +
  xlim(0, max(dataPlot$x)) +
  ylim(0, max(dataPlot$y_new)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none")
ggsave(plot2, file = paste0(plotDirectory, "P10FillNEWSC.png"), width = 6, height = 6)


########################################################################################################
######################################## Figure 2A #####################################################
########################################################################################################

dataDirectory <- "schematic/"
plotDirectory <- "plots/figure2/"

fullMasks <- tibble(read.table(file = paste0(dataDirectory, "masksShape.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
data_tablePlot <- subset(fullMasks, value != 0)
unique_values <- unique(data_tablePlot$value)

mapping <- data.table(
  from_values = unique_values,
  to_values = sample(unique_values,length(unique_values) , replace = FALSE)[1:50] 
)

dataPlot <- merge(data_tablePlot, mapping, by.x = "value", by.y = "from_values", all.x = TRUE)


plot1 <- ggplot(data = dataPlot , aes(x=x, y=y_new, color=as.factor(to_values)))  + geom_point(size=0.0000001)+
  scale_color_igv() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),  
        axis.text.x = element_blank(),   
        axis.text.y = element_blank(),    
        panel.border=element_blank(), 
        axis.line = element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none") 

ggsave(plot1, file = paste0(plotDirectory, "cellsMasks.png"), width = 20.5, height = 5)


fullMasks <- tibble(read.table(file = paste0(dataDirectory, "masksNuclei.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
data_tablePlot <- subset(fullMasks, value != 0)
unique_values <- unique(data_tablePlot$value)

mapping <- data.table(
  from_values = unique_values,
  to_values = sample(unique_values,length(unique_values) , replace = FALSE)[1:50] 
)

dataPlot <- merge(data_tablePlot, mapping, by.x = "value", by.y = "from_values", all.x = TRUE)


plot2 <- ggplot(data = dataPlot , aes(x=x, y=y_new, color=as.factor(to_values)))  + geom_point(size=0.0000000001)+
  scale_color_igv() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),  
        axis.text.x = element_blank(),   
        axis.text.y = element_blank(),    
        panel.border=element_blank(), 
        axis.line = element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none") 

ggsave(plot2, file = paste0(plotDirectory, "NucleiMasks_v1.png"), width = 20.5, height = 5)
