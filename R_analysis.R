```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# Load libraries:

library(tidyverse);packageVersion("tidyverse")
library(ggplot2);packageVersion("ggplot2")
library(phyloseq);packageVersion("phyloseq") 
library(RColorBrewer);packageVersion("RColorBrewer") 
library(ComplexHeatmap);packageVersion("ComplexHeatmap") 
library(ggpubr);packageVersion("ggpubr") 
library(readxl);packageVersion("readxl") 
library(ADImpute);packageVersion("ADImpute") 
library(patchwork);packageVersion("patchwork")
library(cowplot);packageVersion("cowplot") 
library(rstatix);packageVersion("rstatix")
library(svglite);packageVersion("svglite")
library(circlize);packageVersion("circlize")
library(reshape2);packageVersion("reshape2")
library(hilldiv);packageVersion("hilldiv")
```

```{r Load virome data}
#metadata
metadata_cf <- read.csv("~/Desktop/01_CrappyFish_manuscript/data/virome_data/metadata.csv", sep = ",")  #%>% filter(!Sample_name %in% c("S02", "H04")) #failed samples removed - only 38 used downstream
rownames(metadata_cf) <- metadata_cf$Sample_name %>% as.matrix()

# virus df
virus_tax_cf <-  read_excel("/Users/zts270/Desktop/01_CrappyFish_manuscript/data/move_to_dori/vOTU_files/tax_vOTU_overview.xlsx", sheet = "Sheet1") 
virus_abun_cf <- read.csv("/Users/zts270/Desktop/01_CrappyFish_manuscript/data/move_to_dori/rpkm_0.02_abundance_table.csv", sep = ";") 
virus_length_cf <- read.csv("/Users/zts270/Desktop/01_CrappyFish_manuscript/data/move_to_dori/vOTU_files/vOTU_length.tsv", sep = "\t")

```


```{r Load bacterial 16S rRNA data}
# bacteria df
bacteria_16S <- read.csv("~/Desktop/01_CrappyFish_manuscript/data/16S_data/Bozzi_OTUs.csv", sep = ";") %>% as.data.frame() 
bac.tax <- read.csv("~/Desktop/01_CrappyFish_manuscript/data/16S_data/Bozzi_bacteria_taxonomy.csv")
bac.meta <- read.csv("~/Desktop/01_CrappyFish_manuscript/data/16S_data/Bozzi_bacteria_metadata.csv", sep = ";") 

virome.bac <- read.csv("/Users/zts270/Desktop/01_CrappyFish_manuscript/data/16S_data/rpkm_0.02_abundance_table_bacteria.csv", sep = ";") # virome data but only in samples overlapping with samples from which 16S has been generated

bacteria_MAGs <- read.csv("/Users/zts270/Desktop/01_CrappyFish_manuscript/MAGs/CF_mags_no2/MAG_cov_aliMyco.csv", sep = ",")  

otus_in_otu_table <- colnames(bacteria_16S)[-1]

bac.tax <- bac.tax[bac.tax$Column1 %in% otus_in_otu_table, ]
```

```{r average weighted mean length by abundance}
average_abun <- virus_abun_cf
colnames(average_abun)[colnames(average_abun) == "vOTU"] <- "seq_name"

average_abun$total_abundance <- rowSums(average_abun[, -1])

average_length <- mean(virus_length_cf$length)

merged_df <- merge(virus_length_cf, average_abun, by = "seq_name")

weighted_avg_length <- weighted.mean(merged_df$length, merged_df$total_abundance)
weighted_avg_length

```

### Relative abudance plot of 16S rRNA data
```{r}
Myco <- bac.tax %>%
  select(Column1,Family) %>%
  filter(Family == "Mycoplasmataceae") %>%
  .$Column1
Samples <- bacteria_16S$SampleID


vib <- bac.tax %>%
 select(Column1,Genus) %>%
 filter(Genus == "Vibrio") %>%
 .$Column1

Ali <- bac.tax %>%
 select(Column1,Genus) %>%
 filter(Genus == "Aliivibrio") %>%
 .$Column1

Myco.df <- bacteria_16S %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  hilldiv::tss() %>%
  t() %>%
  as_tibble() %>%
  select(all_of(Myco[1:8])) %>%
  mutate(SampleID = Samples, .before = 1) %>%
  reshape2::melt() %>%
  group_by(SampleID) %>%
  summarise(Mycoplasma = sum(value))

vib.df <- bacteria_16S %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  hilldiv::tss() %>%
  t() %>%
  as_tibble() %>%
  select(all_of(vib[1:1])) %>%
  mutate(SampleID = Samples, .before = 1) %>%
  reshape2::melt() %>%
  group_by(SampleID) %>%
  summarise(Vibrio = sum(value))


Ali.df <- bacteria_16S %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  hilldiv::tss() %>%
  t() %>%
  as_tibble() %>%
  select(all_of(Ali[1:5])) %>%
  mutate(SampleID = Samples, .before = 1) %>%
  reshape2::melt() %>%
  group_by(SampleID) %>%
  summarise(Aliivibrio = sum(value))

all_OTUs <- bacteria_16S %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  hilldiv::tss() %>%
  t() %>%
  as_tibble() %>%
  mutate(SampleID = Samples, .before = 1) %>%
  reshape2::melt() %>%
  group_by(SampleID) %>%
  summarise(Total = sum(value))

combined_df <- all_OTUs %>%
  left_join(Myco.df, by = "SampleID") %>%
  left_join(Ali.df, by = "SampleID") %>%
  left_join(vib.df, by = "SampleID") %>%
  mutate(Other = Total - Mycoplasma - Aliivibrio - Vibrio) %>%
  select(-Total) %>%
  arrange(desc(Ali.df$Aliivibrio))

combined_df <- combined_df %>%
  gather(key = "Genus", value = "Value", -SampleID)  

combined_df$SampleID <- factor(combined_df$SampleID, levels = unique(combined_df$SampleID))


library(colorspace)

plot_16S <-ggplot(combined_df, aes(x = SampleID, y = Value, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.3), 
        legend.text = element_text(size = 10, angle = ),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  labs(y = "Relative Abundance", title = "16S rRNA microbiota profile", x = "Samples") +
  scale_fill_manual(values = c(
    scales::alpha("goldenrod1", 0.9), 
     scales::alpha("salmon2", 0.9),
     scales::alpha("tomato4", 0.8), 
     scales::alpha("goldenrod3", 0.8))) + 
  geom_tile() +
  theme(axis.line = element_blank(),
        panel.grid = element_blank())

plot_16S

# svglite("plot_16S.svg", width = 18, height = 12)
# plot(plot_16S)
# dev.off()
```

```{r MAG vs 16S relative abundance profile}
Myco <- bac.tax %>%
  select(Column1, Family) %>%
  filter(Family == "Mycoplasmataceae") %>%
  pull(Column1)

vib <- bac.tax %>%
  select(Column1, Genus) %>%
  filter(Genus == "Vibrio") %>%
  pull(Column1)

Ali <- bac.tax %>%
  select(Column1, Genus) %>%
  filter(Genus == "Aliivibrio") %>%
  pull(Column1)

Samples_16S <- bacteria_16S$SampleID
Samples_MAGs <- bacteria_MAGs$SampleID

# extract genus function
extract_genus_df <- function(bacteria_df, genus_list, genus_name, sample_ids) {
  bacteria_df %>%
    column_to_rownames("SampleID") %>%
    t() %>%
    hilldiv::tss() %>%
    t() %>%
    as_tibble() %>%
    select(all_of(genus_list)) %>%
    mutate(SampleID = sample_ids, .before = 1) %>%
    melt() %>%
    group_by(SampleID) %>%
    summarise("{genus_name}" := sum(value))
}

# 16S setup
Myco.df.16S <- extract_genus_df(bacteria_16S, Myco[1:8], "Mycoplasma", Samples_16S)
Ali.df.16S  <- extract_genus_df(bacteria_16S, Ali[1:5], "Aliivibrio", Samples_16S)
Vib.df.16S  <- extract_genus_df(bacteria_16S, vib[1:1], "Vibrio", Samples_16S)

all_OTUs.16S <- bacteria_16S %>%
  column_to_rownames("SampleID") %>%
  t() %>%
  hilldiv::tss() %>%
  t() %>%
  as_tibble() %>%
  mutate(SampleID = Samples_16S, .before = 1) %>%
  melt() %>%
  group_by(SampleID) %>%
  summarise(Total = sum(value))

combined_df_16S <- all_OTUs.16S %>%
  left_join(Myco.df.16S, by = "SampleID") %>%
  left_join(Ali.df.16S, by = "SampleID") %>%
  left_join(Vib.df.16S, by = "SampleID") %>%
  mutate(Other = Total - Mycoplasma - Aliivibrio - Vibrio) %>%
  select(-Total) %>%
  gather(key = "Genus", value = "Value", -SampleID)

# MAG setup
Myco.df.MAGs <- extract_genus_df(bacteria_MAGs, Myco[1:1], "Mycoplasma", Samples_MAGs)
Ali.df.MAGs  <- extract_genus_df(bacteria_MAGs, Ali[1:1], "Aliivibrio", Samples_MAGs)

all_OTUs.MAGs <- bacteria_MAGs %>%
  column_to_rownames("SampleID") %>%
  t() %>%
  hilldiv::tss() %>%
  t() %>%
  as_tibble() %>%
  mutate(SampleID = Samples_MAGs, .before = 1) %>%
  melt() %>%
  group_by(SampleID) %>%
  summarise(Total = sum(value))

combined_df_MAGs <- all_OTUs.MAGs %>%
  left_join(Myco.df.MAGs, by = "SampleID") %>%
  left_join(Ali.df.MAGs, by = "SampleID") %>%
  mutate(Other = Total - Mycoplasma - Aliivibrio) %>%
  select(-Total) %>%
  gather(key = "Genus", value = "Value", -SampleID)

ordered_from_16S <- combined_df_16S %>%
  filter(Genus == "Aliivibrio") %>%
  arrange(desc(Value)) %>%
  pull(SampleID) %>%
  unique()

mag_only_samples <- setdiff(combined_df_MAGs$SampleID, combined_df_16S$SampleID)

sample_order <- c(ordered_from_16S, mag_only_samples)

all_samples <- union(combined_df_16S$SampleID, combined_df_MAGs$SampleID)
all_genera_16S <- unique(combined_df_16S$Genus)
all_genera_MAGs <- unique(combined_df_MAGs$Genus)

full_grid_16S <- expand.grid(SampleID = all_samples, Genus = all_genera_16S)
full_grid_MAGs <- expand.grid(SampleID = all_samples, Genus = all_genera_MAGs)

combined_df_16S <- full_grid_16S %>%
  left_join(combined_df_16S, by = c("SampleID", "Genus")) %>%
  mutate(Value = replace_na(Value, 0))

combined_df_MAGs <- full_grid_MAGs %>%
  left_join(combined_df_MAGs, by = c("SampleID", "Genus")) %>%
  mutate(Value = replace_na(Value, 0))

combined_df_16S$SampleID <- factor(combined_df_16S$SampleID, levels = sample_order)
combined_df_MAGs$SampleID <- factor(combined_df_MAGs$SampleID, levels = sample_order)

plot_16S <- ggplot(combined_df_16S, aes(x = SampleID, y = Value, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1)) +
  labs(title = "16S rRNA microbiota profile", y = "Relative Abundance", x = NULL) +
  scale_fill_manual(values = c("goldenrod1", "salmon2", "tomato4", "goldenrod3"))

plot_MAGs <- ggplot(combined_df_MAGs, aes(x = SampleID, y = Value, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1)) +
  labs(title = "MAGs microbiota profile", y = "Relative Abundance", x = "Samples") +
  scale_fill_manual(values = c("darkgoldenrod1", "salmon2", "tomato4"))

mag16plot <- plot_grid(plot_16S, plot_MAGs, ncol = 1, align = "v")
mag16plot

# svglite("Mmag16plot.svg", width = 18, height = 12)
# plot(mag16plot)
# dev.off()


```

### Phyloseq object generation:
```{r phyloseq object setup}
virus_abun_cf2 <- virus_abun_cf %>%
  filter(vOTU %in% virus_length_cf$seq_name) %>%
 column_to_rownames(var = "vOTU") 

virus_abun_cf[is.na(virus_abun_cf)] <- 0
virus_abun_cf <- virus_abun_cf %>% as_tibble()

physeq <- phyloseq(otu_table((virus_abun_cf), taxa_are_rows = TRUE),
                             tax_table(as.matrix(virus_tax_cf)),
                             sample_data(metadata_cf))

subset_alii <- subset_taxa(physeq, Host_family == "Vibrionaceae")

```


```{r}
barplot_all <- subset_samples(physeq, Category %in% c("Mycoplasma sp.", "Aliivibrio sp.", "Unknown microbiota"))
otu_table(barplot_all) <- otu_table(barplot_all) + 1
sample_order <- sample_names(barplot_all)[order(sample_data(barplot_all)$Category)]

plot_1 <- barplot_all
taxa_names(plot_1) <- tax_table(plot_1)[, "vOTU"]

p1 <- plot_heatmap(plot_1,
                   method = "NMDS",
                   distance = "bray",
                   sample.label = "Sample_name",
                   taxa.label = "vOTU")

ord_sample_order <- unique(p1$data$Sample)

p1 <- plot_heatmap(plot_1,
                   method = NULL,
                   sample.order = ord_sample_order,
                   sample.label = "Sample_name",
                   taxa.label = "vOTU") +
  scale_fill_gradientn(colors = c("white", "dodgerblue2", "lightblue")) +
  labs(y = "vOTU") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size=8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10)) +
  ggforce::facet_row(~Category, scales = "free_x", space = "free")

ord_sample_order <- as.character(ord_sample_order)

plot_2 <- subset_samples(subset_alii, Category %in% c("Mycoplasma sp.", "Aliivibrio sp.", "Unknown microbiota"))


plot_2 <- prune_samples(ord_sample_order, plot_2)

p2 <- plot_heatmap(plot_2,
                   method = NULL,
                   sample.order = ord_sample_order,
                   sample.label = "Sample_name",
                   taxa.label = "vOTU") +
  scale_fill_gradientn(colors = c("white", "goldenrod1", "orange")) +
  labs(y = "vOTU") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size=8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10)) +
  ggforce::facet_row(~Category, scales = "free_x", space = "free")

x_title <- ggdraw() + draw_label("Samples", size = 12, hjust = 0.5)

p1 <- cowplot::plot_grid(p1, x_title, ncol = 1, rel_heights = c(1, 0.05))
p2 <- cowplot::plot_grid(p2, x_title, ncol = 1, rel_heights = c(1, 0.05))

combined <- cowplot::plot_grid(p1, p2, ncol = 1, labels = c("A", "B"), label_size = 10, rel_heights = c(26,11))
plot(combined)


# svglite("heatmaps_vOTUs_figure3.svg", width = 13, height = 7)
# plot(combined)
# dev.off()

```

### Alpha diversity analysis between the three microbiota groups (Aliivibrio, Mycoplasma or no 16S)
```{r}
physeq.mycoAlii.nn <- subset_samples(physeq, Category %in% c("Mycoplasma sp.","Aliivibrio sp.", "Unknown microbiota"))

richness <- estimate_richness(physeq, measures = c("Observed","Shannon", "Simpson"))
md.richness <- metadata_cf[match(rownames(richness),metadata_cf$Sample_name),]
richness <- cbind(md.richness, richness)

y_labels <- c("vOTU richness", "Shannon Index") 

plot_list <- list()
for (i in c("Observed", "Shannon")) { 
    divtestdata <- data.frame(vOTU = richness[, i], Category = richness$Category)

    shapiro <- shapiro.test(divtestdata$vOTU)
    shapiro$p.value
    
    if (shapiro$p.value > 0.05) {
        test <- c("TukeyHSD")
        stat.test <- divtestdata %>%
            group_by("Category") %>%
            rstatix::tukey_hsd(vOTU ~ Category) %>%
            adjust_pvalue(method = "fdr") %>%
            add_significance("p.adj")
    } else {
        test <- c("NP Dunn")
        stat.test <- divtestdata %>%
            group_by("Category") %>%
            dunn_test(vOTU ~ Category) %>%
            adjust_pvalue(method = "fdr") %>%
            add_significance("p.adj")
    }

    stat.test$p.adj <- round(stat.test$p.adj, 4) 

    stat.test <- stat.test %>%
        add_x_position(x = "Category", dodge = 0.8) %>%
        add_y_position()

    level_colors <- list(Bacterial_profile = c(Aliivibrio = "goldenrod1", Mycoplasma = "salmon2", none = "tomato4"))

    divtestdata$Myco90 <- factor(divtestdata$Category, levels = c("Aliivibrio sp.", "Mycoplasma sp.", "Unknown microbiota"))

    plot <- ggboxplot(
        divtestdata, x = "Category", y = "vOTU", 
        color = "black",
        fill = "Category",
        outlier.shape = 8,
        size = 0.5,
        title = ""
    ) +
    scale_color_manual(values = c("goldenrod1", "salmon2", "tomato4")) + 
    stat_pvalue_manual(
        stat.test, label = paste(test, ": p.adj={p.adj}", sep = ""), tip.length = 0.05, label.size = 4,
        step.increase = 0.04
    ) +
    geom_boxplot(aes(col = Myco90)) +
    geom_point(aes(col = Myco90), alpha = 0.5) +
    theme_minimal() +
    xlab("") +
    ylab(y_labels[which(c("Observed", "Shannon") == i)]) +  
    theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90, hjust = 1)) +
      theme_bw()

    plot_list[[i]] <- plot}

plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "none")

legend_b <- get_legend(
    plot_list[[2]] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "none"))

plot_list[[2]] <- plot_list[[2]] + theme(legend.position = "none") 

combined_plot <- cowplot::plot_grid(
    plot_list[[1]], plot_list[[2]], 
    nrow = 1)

diversity_plot <- cowplot::plot_grid(combined_plot, 
    ncol = 1, 
    rel_heights = c(0.1, 1))

plot(diversity_plot)

# svglite("alpha_div_plot3_new.svg", width = 8, height = 7)
# plot(diversity_plot)
# dev.off()
```


### MAG module completion
```{r, print=F}
module_completion <- read.table(file='/Users/zts270/Desktop/01_CrappyFish_manuscript/MAGs/CF_mags_no2/metabolic_independence/metabolism-module_pathwise_completeness-MATRIX.txt', header = TRUE, sep = "\t")

module_info <- read.csv("~/Desktop/01_CrappyFish_manuscript/MAGs/CF_mags_no2/metabolic_independence/modules_info.txt", sep = "\t")

external_genomes <- read.table(file='/Users/zts270/Desktop/01_CrappyFish_manuscript/MAGs/CF_mags_no2/metabolic_independence/external-genomes.txt', header = TRUE, sep = "\t")

external_genomes <- external_genomes %>%
  dplyr::mutate(Group = c(rep("Mycoplasma", 21), rep("Aliivibrio", 11)))

df <- melt(module_completion)
colnames(df) <- c('module', 'genome', 'completion')

df$group <- external_genomes$Group[match(df$genome, external_genomes$name)]
df$individual <- external_genomes$individual[match(df$genome, external_genomes$name)]

p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")}}

groups_order <- c("Mycoplasma", "Aliivibrio")
group_colors <- c("goldenrod1", "salmon2") 

p_value <- wilcox.test(df[df$group == "Mycoplasma", ]$completion, df[df$group == "Aliivibrio", ]$completion, exact = FALSE)$p.value

plot_all_genomes <- ggplot(data=df, aes(x=group, y=completion, group=group)) +
    geom_violin(aes(fill=group), color="black", fill="grey", alpha=0.3) +
    geom_jitter(aes(shape=group, color=group), width=0.2, height=0.01, size=2.5, alpha=0.1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_shape_manual(values = c(16, 17),
  labels = c("Mycoplasma" = expression(italic("Mycoplasma")), "Aliivibrio" = expression(italic("Aliivibrio")))) +  
  annotate("text", x=1.5, y=1.03, label=paste(p_to_stars(p_value)), size=3.5) +
  geom_segment(aes(x=1, xend=2, y=1.025, yend=1.025), size=0.5) +  
  geom_segment(aes(x=1, xend=1, y=1.025, yend=1.02), size=0.5) +   
  geom_segment(aes(x=2, xend=2, y=1.025, yend=1.02), size=0.5) +   
  theme_classic() +
  theme(legend.position="none",
        axis.title.y = element_text(size=14, vjust=0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14)) +
  ylab("Estimated KEGG module completion") +
  xlab("Pangenomes") +
  scale_x_discrete(limits = groups_order, labels = c("Mycoplasma" = expression(italic("Mycoplasma")), 
                                                     "Aliivibrio" = expression(italic("Aliivibrio")))) +  
  scale_color_manual(values = group_colors,
                     labels = c("Mycoplasma" = expression(italic("Mycoplasma")), 
                                "Aliivibrio" = expression(italic("Aliivibrio")))) 


df_filt <- df %>% filter(genome %in% c("CF_MAG_00002", "CF_MAG_00001_contigs"))
p_value_filt <- wilcox.test(df_filt[df_filt$group == "Mycoplasma", ]$completion, df_filt[df_filt$group == "Aliivibrio", ]$completion, exact = FALSE)$p.value

plot_CF_genomes <- df %>% 
  filter(genome %in% c("CF_MAG_00002", "CF_MAG_00001_contigs")) %>%
  ggplot(aes(x=group, y=completion, group=group)) +
    geom_violin(aes(fill=group), color="black", fill="grey", alpha=0.3) +
    geom_jitter(aes(shape=group, color=group), width=0.2, height=0.01, size=2.5, alpha=0.3) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_shape_manual(values = c(16, 17)) +  # Circle for Mycoplasma, triangle for Aliivibrio
  annotate("text", x=1.5, y=1.03, label=paste(p_to_stars(p_value_filt)), size=3.5) +
  geom_segment(aes(x=1, xend=2, y=1.025, yend=1.025), size=0.5) +   
  geom_segment(aes(x=1, xend=1, y=1.025, yend=1.02), size=0.5) +  
  geom_segment(aes(x=2, xend=2, y=1.025, yend=1.02), size=0.5) +   
  theme_classic() +
  theme(legend.position="right",
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=14, vjust=0.5),
        axis.title.x = element_text(size = 14),
        plot.margin = unit(c(2,2,2,2), "cm")) +
  ylab(NULL) +
  xlab("Study Genomes") +
  labs(fill = "Genus", shape = "Genus", color = "Genus") +
   guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
  scale_x_discrete( limits = groups_order, labels = c("Mycoplasma" = "SSM", "Aliivibrio" = "SSA")) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
   guides(shape = guide_legend(override.aes = list(alpha = 1, size = 5)), 
         color = guide_legend(override.aes = list(alpha = 1))) +
  scale_fill_manual(values = group_colors,
                  labels = c("Mycoplasma" = expression(italic("Mycoplasma")), 
                             "Aliivibrio" = expression(italic("Aliivibrio")))) +
scale_color_manual(values = group_colors,
                   labels = c("Mycoplasma" = expression(italic("Mycoplasma")), 
                              "Aliivibrio" = expression(italic("Aliivibrio")))) +
scale_shape_manual(values = c(16, 17),
                   labels = c("Mycoplasma" = expression(italic("Mycoplasma")), 
                              "Aliivibrio" = expression(italic("Aliivibrio"))))


df %>% filter(completion == 1) %>% group_by(group) %>% summarise(n = n())

MoI <- df %>% filter(genome %in% c("CF_MAG_00002", "CF_MAG_00001_contigs")) %>%
  filter(completion == 1) %>% 
  group_by(module, genome) %>% 
  summarise(completion) %>% 
  summarise(n = n()) %>% 
  pull(module)

MoI <- df %>% filter(genome %in% c("CF_MAG_00002", "CF_MAG_00001_contigs")) %>%
  filter(module %in% MoI) %>% arrange(module) %>% 
  left_join(module_info, by = "module")

MoI_plot <- MoI %>% 
  ggplot(aes(x = group, subcategory, fill = completion)) + 
  geom_tile(color = "white", size = 0.02) +
  scale_fill_gradient(low = "white", high = "darkslategrey") +
  scale_x_discrete( limits = groups_order, labels = c("Mycoplasma" = "SSM", "Aliivibrio" = "SSA")) +
  labs(y = "",
       x = "",
       fill = "Completion") +
  theme_minimal() +
  theme(axis.line = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(margin = margin(b = 20), size=14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


metabolic_plot <- ((plot_all_genomes + plot_CF_genomes) | MoI_plot + plot_layout(guides = "collect", widths = c(3, 1)))

metabolic_plot


# svglite("fig5_metabolic_plot3.svg", width = 17.5, height = 9)
# plot(metabolic_plot)
# dev.off()

```

### Plot vOTU annotations
```{r}
data <- read.table("/Users/zts270/Desktop/01_CrappyFish_manuscript/data/move_to_dori/pharokka_cds_functions.tsv", 
                   header=TRUE, sep="\t") 

heatmap_data <- dcast(data, Description ~ vOTU, value.var = "Count", fill = 0) 

rownames(heatmap_data) <- heatmap_data$Description
heatmap_data <- heatmap_data[, -1]

heatmap_data[heatmap_data > 0] <- 1 

function_descriptions <- rownames(heatmap_data)

function_categories <- sapply(function_descriptions, function(x) {
  if (grepl("Head|Tail|Connector", x, ignore.case = TRUE)) {
    return("Structural")
  } else if (grepl("DNA|RNA|Nucleotide|Transcription", x, ignore.case = TRUE)) {
    return("Phage Metabolism")
  } else if (grepl("Lysis", x, ignore.case = TRUE)) {
    return("Lytic")
  } else if (grepl("Integration and Excision", x, ignore.case = TRUE)) {
    return("Temperate")
  } else if (grepl("tRNA|CRISPR|tmRNA|Moron|Virulence|AMR", x, ignore.case = TRUE)) {
    return("Accessory Genes")
  } else if (grepl("Other|CDS|Unknown", x, ignore.case = TRUE)) {
    return("Unknown/CDS") 
  } 
})

row_annotation <- data.frame(Category = factor(function_categories))
rownames(row_annotation) <- function_descriptions

row_annotation <- as.data.frame(row_annotation)
rownames(row_annotation) <- rownames(heatmap_data)

sorted_indices <- order(row_annotation$Category)
heatmap_data <- heatmap_data[sorted_indices, ]
row_annotation <- row_annotation[sorted_indices, , drop = FALSE]

cds_counts <- aggregate(data$Count, by = list(vOTU = data$vOTU), sum)
colnames(cds_counts) <- c("vOTU", "CDS_count")

column_annotation <- data.frame(CDS_count = cds_counts$CDS_count)
rownames(column_annotation) <- cds_counts$vOTU

column_annotation <- column_annotation[match(colnames(heatmap_data), rownames(column_annotation)), , drop = FALSE]

category_colors <- c(
  "Structural" = "deepskyblue3",
  "Phage Metabolism" = "gold2",
  "Lytic" = "orangered4",
  "Temperate" = "aquamarine3",
  "Accessory Genes" = "tan",
  "Unknown/CDS" = "darkorange3"
)

row_anno <- rowAnnotation(
  "Annotation Category" = factor(row_annotation$Category), 
  col = list("Annotation Category" = category_colors),
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 8, fontface = "plain"),  
  annotation_legend_param = list(
    "Annotation Category" = list(  
      title_gp = gpar(fontface = "plain", fontsize = 10)  
    )
  )
)

color_palette <- c("0" = "ghostwhite", "1" = "darkslategrey")


anno_heatmap <- Heatmap(
  as.matrix(heatmap_data),
  name = "Presence/Absence",
  col = color_palette,
  top_annotation = HeatmapAnnotation("CDS count" = anno_barplot(cds_counts$CDS_count, 
                                                               gp = gpar(fill = "lightblue4", col = NA), 
                                                               axis = TRUE),
                                     annotation_name_gp = gpar(fontsize = 10)),
  right_annotation = row_anno,  
  cluster_rows = FALSE,          
  cluster_columns = FALSE,      
  show_row_names = TRUE,         
  show_column_names = TRUE,      
  row_names_gp = gpar(fontsize = 10),  
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title_gp = gpar(fontface = "plain", fontsize = 10)
  )
)

anno_heatmap

# svglite("pharokka_heatmap_fig4.svg", width = 9, height = 5)
# plot(anno_heatmap)
# dev.off()

```

### Cycle threshold correlation analysis with Aliivibrio vOTUs.
```{r}
CTs <- read.csv("/Users/zts270/Bioinformatics/R_data_analyses/PhD_viruses_and_metagenomes_R/CrappyFish/ct_plot_input_file.csv") %>% as.data.frame()
md <- read.csv("/Users/zts270/Bioinformatics/R_data_analyses/PhD_viruses_and_metagenomes_R/CrappyFish/Metadata_filtered.csv")
CSW_df <- read.csv("/Users/zts270/Bioinformatics/R_data_analyses/PhD_viruses_and_metagenomes_R/CrappyFish/CSW_Bozzi_Samples.csv", sep=",", na.strings = "") %>% as.data.frame()

CTs <- CTs %>%
  full_join(md, by = "Sample") %>%
  na.omit()

CSW_df <- CSW_df %>% 
  rename(Sample = Bozzi_et_al) %>%
  na.omit() %>%
  full_join(CTs, by = "Sample")

physeq_ct <- phyloseq(otu_table((virus_abun_cf), taxa_are_rows = TRUE),
                             tax_table(as.matrix(virus_tax_cf)),
                             sample_data(metadata_cf))

richness <- estimate_richness(physeq_ct, measures = "Observed") 
richness$Sample_names <- rownames(richness)  

merged_data <- merge(richness, CSW_df, by = "Sample_names") 

aliivibrio_data <- merged_data %>% filter(Myco90 == "Aliivibrio")

ct_vs_vOTU <- ggplot(aliivibrio_data, aes(x = Observed, y = Ct, color = Myco90)) +  
  geom_point(color = "goldenrod1") +
  geom_smooth(method = "lm", color = "goldenrod1") +  
  labs(x = "vOTU Richness", y = "Ct Value", title = "Bacterial biomass and vOTU richness") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., "~italic(P)*`=`~", ..p.value..)), 
           label.x = Inf, label.y = Inf, 
           hjust = 1.15, vjust = 1.8, 
           method = "spearman", color = "black")

ct_vs_vOTU

# svglite("ct_vs_vOTU_figure3.svg", width = 5, height = 4.5)
# plot(ct_vs_vOTU)
# dev.off()

```
