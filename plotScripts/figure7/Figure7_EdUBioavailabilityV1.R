
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library("ggsci")

home1Directory <- "data/"
plotDirectory <- "plots/figure5/"


## EduData for bioavailability
EdU_Ki67 <- tibble(read.table(file = paste0(home1Directory, "231005_EdU bioavailability quantification.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

EdU_Ki67 <- EdU_Ki67 %>% mutate(EdUpercent = 100*EdUpercent, Ki67percent = 100*Ki67percent) %>% select(-Ki67, -Replicate)

EdU_Ki67_tall <- gather(EdU_Ki67, "type", "percentPositive", Ki67percent, EdUpercent)

EdU_Ki67_tall <- EdU_Ki67_tall %>% group_by(Time) 

EdU_Ki67_tall$Time = factor(EdU_Ki67_tall$Time, c("fresh", "20", "60", "180"))
EdU_Ki67_tall$type = factor(EdU_Ki67_tall$type, c("Ki67percent", "EdUpercent"))


EdU_Ki67_tall.summary = EdU_Ki67_tall %>% group_by(Time, type) %>% summarise(mean = mean(percentPositive))


plot <- ggplot(data=EdU_Ki67_tall, aes(x= Time, y=percentPositive, fill=type)) +
  geom_point(pch = 21, position = position_jitterdodge()) +
  geom_boxplot(outlier.size = 0) +
  theme_classic()

ggsave(plot, file = paste0(plotDirectory, "SuppFigure_EdU_Bioavailability.svg"), width = 6, height = 4)
