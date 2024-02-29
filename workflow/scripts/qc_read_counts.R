library(tidyverse)

data <- read_csv(snakemake@input[[1]])

samples <- data |>
  filter(str_ends(Annotation, "reads")) |>
  mutate(
    direction = word(Annotation, 1, sep = " "),
    type = word(Annotation, 2, -1)
    )

p1 <- ggplot(samples, aes(x = Count, fill = Annotation), alpha=0.5)+
  geom_histogram() +
  theme_classic() +
  xlab("Number of reads") +
  ylab("Number of samples") +
  facet_wrap(~type, scale = "free")+
  theme(legend.position = "bottom") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))

ggsave(snakemake@output[[1]], plot = p1, device = "pdf", width = 11, height = 8.5, units = "in", dpi = 300)

custom_order <- c("raw reads", "trimmed reads", "SSU reads", "rRNA assembly", "not rRNA reads", "mRNA assembly", "Filtered mRNA assembly")
agg <- samples |>
  group_by(type) |>
  summarize(Count = sum(Count)) |>
  bind_rows(
    data |>
      filter(str_ends(Annotation, "assembly")) |>
      rename(type = Annotation)
  ) |>
  mutate(type = factor(type, levels = custom_order, labels = paste0(type, " (", format(Count, scientific = TRUE), ")")))

p2 <- agg |>
 ggplot(aes(x = Count, y = type, fill = type))+
  geom_col()+
  xlab("Total number of sequences") +
  ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), trans = "sqrt") + 
  theme_classic() +
  theme(legend.position = "none")

ggsave(snakemake@output[[2]], plot = p2, device = "pdf", width = 11, height = 8.5, units = "in", dpi = 300)
