library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library("ggsci")
library(forcats)

dataDirectory <- "data/inducible rainbow/"

plotDirectory = "plots/"

df <- tibble(read.table(file = "data/inducible rainbow/allClusterSizesRed.csv", header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

df_help <- df %>%
  select(-average, -maximum)

long_df <- melt(data = df_help, id.vars = c("age", "aorta"), 
                variable.name = "cluster_size", value.name = "count")
long_df$cluster_size <- as.integer(sub("X", "", long_df$cluster_size))

ages_split <- strsplit(long_df$age, "-")
long_df$start <- as.integer(sapply(ages_split, function(x) as.integer(sub("P", "", x[1]))))
long_df$end <- as.integer(sapply(ages_split, function(x) as.integer(x[2])))

long_df <- long_df %>%
  group_by(age, aorta) %>%
  mutate(fractionAorta = count / sum(count)) %>%
  ungroup()


df_aorta <- long_df %>%
  filter(!is.na(cluster_size)) %>% 
  group_by(age, aorta, cluster_size) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

df_aorta <- df_aorta %>%
  mutate(clusterSizeAggregated = case_when(
    cluster_size == 1 ~ "1",
    cluster_size == 2 ~ "2",
    cluster_size == 3 ~ "3",
    cluster_size == 4 ~ "4",
    cluster_size >= 5 ~ "5+"
  ))

summary_df_aorta <- df_aorta %>%
  group_by(age, aorta, clusterSizeAggregated) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

summary_df_aorta$clusterSizeAggregated <- factor(summary_df_aorta$clusterSizeAggregated, levels = c("1", "2", "3", "4", "5+"))
summary_df_aorta <- summary_df_aorta %>%
  group_by(age, aorta) %>%
  mutate(fractionCluster = count / sum(count)) %>%
  ungroup()

summary_df_aorta$clusterSizeAggregated <- factor(summary_df_aorta$clusterSizeAggregated)
summary_df_aorta$clusterSizeAggregated <- fct_rev(summary_df_aorta$clusterSizeAggregated)


df_age <- long_df %>%
  filter(!is.na(cluster_size)) %>% 
  select(-aorta) %>%
  group_by(age, cluster_size) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

df_age <- df_age %>%
  mutate(clusterSizeAggregated = case_when(
    cluster_size == 1 ~ "1",
    cluster_size == 2 ~ "2",
    cluster_size == 3 ~ "3",
    cluster_size == 4 ~ "4",
    cluster_size >= 5 ~ "5+"
  ))

summary_df <- df_age %>%
  group_by(age, clusterSizeAggregated) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

summary_df$clusterSizeAggregated <- factor(summary_df$clusterSizeAggregated, levels = c("1", "2", "3", "4", "5+"))
  

summary_df <- summary_df %>%
  group_by(age) %>%
  mutate(fractionCluster = count / sum(count)) %>%
  ungroup()

summary_df$age <- factor(summary_df$age, levels = c("P30-60","P21-30",
                                                    "P10-60", "P10-30", "P10-21", 
                                                    "P5-60","P5-30", "P5-10", 
                                                    "P0-10", "P0-5"))
summary_df$clusterSizeAggregated <- factor(summary_df$clusterSizeAggregated)
summary_df$clusterSizeAggregated <- fct_rev(summary_df$clusterSizeAggregated)

plot <- ggplot(summary_df, aes(x = age, y = fractionCluster, fill = clusterSizeAggregated)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_igv() + 
  labs(x = "Age", y = "Composition of cluster size", fill = "cluster Size") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22), legend.title = element_text(size = 22), legend.text = element_text(size =18)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  coord_flip()


ggsave(plot, file = paste0(plotDirectory, 'indRainbow_clusterSize_proliferatingCellsOnly.svg'), width = 12, height = 8)


###cluster sizes aorta level
unique_ages <- unique(summary_df_aorta$age)

for (current_age in unique_ages) {
  # Filter the data for the current age
  age_df <- summary_df_aorta %>% filter(age == current_age)
  
  plot <- ggplot(age_df, aes(x = aorta, y = fractionCluster, fill = clusterSizeAggregated)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_igv() + 
    labs(x = "Aorta", y = "Composition of cluster size", fill = "cluster Size") +
    theme_minimal() +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text=element_text(size=18),axis.title=element_text(size=22), legend.title = element_text(size = 22), legend.text = element_text(size =18)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    coord_flip()
  ggsave(
    filename = paste0(plotDirectory, "clusterSize_aorta_", current_age, ".svg"),
    plot = plot,
    width = 12,
    height = 8
  )
}



####Max and Average plot

df_stats <- df[, c("age", "aorta", "average", "maximum")]


library(dplyr)

# Remove rows where count == 0
long_df_filt <- long_df %>% filter(count != 0)

df_stats <- long_df_filt  %>%
  group_by(age, aorta) %>%
  summarise(
    average = sum(cluster_size * count) / sum(count),  # Weighted average
    start = first(start),  # Assuming 'start' is consistent within groups
    end = first(end),      # Assuming 'end' is consistent within groups
    maximum = max(cluster_size)  # For the fill scale
  ) %>%
  ungroup()



df_statsAge <- df_stats %>%
  select(-aorta) %>%
  group_by(age) %>%
  summarise(average = mean(average), maximum = max(maximum))

ages_split <- strsplit(df_statsAge$age, "-")
df_statsAge$start <- as.integer(sapply(ages_split, function(x) as.integer(sub("P", "", x[1]))))
df_statsAge$end <- as.integer(sapply(ages_split, function(x) as.integer(x[2])))

df_statsAge$minAv = df_statsAge$average - 0.015
df_statsAge$maxAv = df_statsAge$average + 0.015

plot <- ggplot() +
  geom_rect(data = df_statsAge, aes(xmin = start, xmax = end, ymin = minAv, ymax = maxAv, fill = maximum), alpha = 0.5) +
  scale_x_continuous(limits = c(0,60.5), breaks = c(0,5,10, 21, 30, 60)) +
  ylim(2, 3) +  xlab("age (days)") + ylab("average cluster size") + 
  scale_fill_material("purple", limits = c(min(df_statsAge$maximum), max(df_statsAge$maximum))) +
  theme(axis.text=element_text(size=23),axis.title=element_text(size=28)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, 'indRainbow_clusterSize_stats_proliferatingCellsOnly.svg'), width = 8, height = 5)





####Log2

df_ageLog2 <- df_age
df_ageLog2$cluster_size <- ceiling(log2(df_ageLog2$cluster_size))

df_ageLog2 <- df_ageLog2  %>%
  filter(!is.na(cluster_size)) %>% 
  group_by(age, cluster_size) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(age) %>%
  mutate(fractionCluster = count / sum(count)) %>%
  ungroup()



  





