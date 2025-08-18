setwd('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/Telomeres')

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(gridExtra)
freq_file='genome_telomere_unit.tsv'
freq_file='telomere.tsv'
df <- read_tsv('genome_telomere_unit.tsv', col_types = cols())
df <- read_tsv(freq_file, col_types = cols())

# Define columns for G and C strands
g_columns <- c('1_T', '2_T', '3_A', '4_G', '5_G', '6_G')
c_columns <- c('-1_C', '-2_C', '-3_C', '-4_T', '-5_A', '-6_A')

# Prepare long format data for G strand
data_g <- df %>%
  filter(Total_G > 0) %>%
  select(Sample, Celltype, all_of(g_columns), Total_G) %>%
  pivot_longer(cols = all_of(g_columns), names_to = "Location", values_to = "Count") %>%
  mutate(
    Frequency = Count / Total_G,
    Strand = "G",
    Location = sub("_.*", "", Location)
  ) %>%
  select(Sample, Celltype, Location, Frequency, Strand)

# Prepare long format data for C strand
data_c <- df %>%
  filter(Total_C > 0) %>%
  select(Sample, Celltype, all_of(c_columns), Total_C) %>%
  pivot_longer(cols = all_of(c_columns), names_to = "Location", values_to = "Count") %>%
  mutate(
    Frequency = Count / Total_C,
    Strand = "C",
    # Convert location index to numeric + 7 as in Python code
    Location = as.character(7 + as.integer( sub("_.*", "", Location)))
  ) %>%
  select(Sample, Celltype, Location, Frequency, Strand)

# Combine data
data <- bind_rows(data_g, data_c)

# Preserve the order of Celltype as in data
data$Celltype <- factor(data$Celltype, levels = unique(data$Celltype))

# Define color palette
base_colors <- c("A" = "#D55E00", "C" = "#0072B2" , "G" = "#F0E442", "T" = "#009E73")

# Palette for G strand locations
palette_g <- sapply(g_columns, function(b) {
  base <- sub(".*_", "", b)
  base_colors[base]
})
names(palette_g) <- sub("_.*", "", g_columns)

# Palette for C strand locations
palette_c <- sapply(c_columns, function(b) {
  base <- sub(".*_", "", b)
  base_colors[base]
})
names(palette_c) <- as.character(7 + as.integer(sub("_.*", "", c_columns)))

summary_data <- data %>%
  group_by(Celltype, Location, Strand) %>%
  summarise(
    mean_freq = mean(Frequency, na.rm = TRUE),
    se_freq = sd(Frequency, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )


# Plot for G strand
p_g <- ggplot(summary_data %>% filter(Strand == "G"), 
              aes(x = Celltype, y = mean_freq, fill = Location)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_freq - se_freq, ymax = mean_freq + se_freq),
                position = position_dodge(width = 0.9), width = 0.5) +
  geom_jitter(data = data %>% filter(Strand == "G"),
              aes(x = Celltype, y = Frequency, color = Location),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9),
              size = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = palette_g) +
  scale_color_manual(values = rep("black", length(palette_g))) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2), expand = c(0, 0))+
  theme_classic(base_size = 16) +
  labs(title = " ", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 20),
    legend.position = "none"
  )

# Plot for C strand
p_c <- ggplot(summary_data %>% filter(Strand == "C"), 
              aes(x = Celltype, y = mean_freq, fill = Location)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_freq - se_freq, ymax = mean_freq + se_freq),
                position = position_dodge(width = 0.9), width = 0.5) +
  geom_jitter(data = data %>% filter(Strand == "C"),
              aes(x = Celltype, y = Frequency, color = Location),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9),
              size = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = palette_c) +
  scale_color_manual(values = rep("black", length(palette_c))) +
  scale_y_reverse(limits = c(0.9, 0), breaks = seq(0, 0.8, 0.2), expand = c(0, 0))+
  scale_x_discrete(position = "top")+
  theme_classic(base_size = 16) +
  labs(title = " ", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 20),
    legend.position = "none"
  )

# Arrange plots vertically
combined_plot <- grid.arrange(p_g, p_c, ncol = 1, heights = c(1, 1))

# Save plot if output filename provided

ggsave(paste(strsplit(freq_file, '\\.')[[1]][1], ".png", sep=""), combined_plot, dpi = 1000, width = 9, height = 6)

