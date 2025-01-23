---
title: "crappyfish_2.0_clean"
author: "Caroline Winther-Have"
date: "2024-01-17"
output: html_document
---

# Load libraries:
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(tidyverse);packageVersion("tidyverse")
library(ggplot2);packageVersion("ggplot2")
library(ape);packageVersion("ape")
library(phyloseq);packageVersion("phyloseq") 
library(survival);packageVersion("survival")
library(Hmisc);packageVersion("Hmisc") 
library(RColorBrewer);packageVersion("RColorBrewer") 
library(ComplexHeatmap);packageVersion("ComplexHeatmap") 
library(ggpubr);packageVersion("ggpubr") 
library(ggrepel);packageVersion("ggrepel") 
library(readxl);packageVersion("readxl") 
library(ADImpute);packageVersion("ADImpute") 
library(DESeq2);packageVersion("DESeq2") 
library(corrplot);packageVersion("corrplot")
library(corrr);packageVersion("corrr") 
library(patchwork);packageVersion("patchwork")
library(cowplot);packageVersion("cowplot") 
library(Matrix);packageVersion("Matrix") 
library(psych);packageVersion("psych")
library(forcats);packageVersion("forcats") 
library(rstatix);packageVersion("rstatix")
library(vegan);packageVersion("vegan") 
library(ggsignif);packageVersion("ggsignif")
library(svglite);packageVersion("svglite")
library(circlize);packageVersion("circlize")
library(reshape2);packageVersion("reshape2")
library(MASS);packageVersion("MASS")
library(ggridges);packageVersion("ggridges")
library(gridExtra);packageVersion("gridExtra")
library(reshape)packageVersion("reshape")
```

### Load virome data
```{r}
#metadata
metadata_cf <- read.csv("/metadata_cf.csv", sep = ";")  %>% filter(!Sample_name %in% c("S02", "H04")) #failed samples removed - only 38 used downstream
rownames(metadata_cf) <- metadata_cf$Sample_name %>% as.matrix()

# virus df
virus_tax_cf <-  read_excel("/tax_vOTU_overview.xlsx", sheet = "Sheet1") 
virus_abun_cf <- read.csv("/rpkm_0.02_abundance_table.csv", sep = ";") 
virus_length_cf <- read.csv("/vOTU_files/vOTU_length.tsv", sep = "\t")

# Calculate average length weighted by abundance
average_abun <- virus_abun_cf
colnames(average_abun)[colnames(average_abun) == "vOTU"] <- "seq_name"

average_abun$total_abundance <- rowSums(average_abun[, -1])

average_length <- mean(virus_length_cf$length)

# Merge the data frames
merged_df <- merge(virus_length_cf, average_abun, by = "seq_name")

# Calculate weighted mean
weighted_avg_length <- weighted.mean(merged_df$length, merged_df$total_abundance)

coverage <- rowSums(average_abun > 0) / ncol(average_abun)
coverage_df <- data.frame(vOTU = rownames(average_abun), Coverage = coverage)
```

### Load bacterial 16S rRNA data
```{r}
bacteria <- read.csv("/Bozzi_OTUs.csv", sep = ";") %>% as.data.frame() 
  
virome.bac <- read.csv("/rpkm_0.02_abundance_table_bacteria.csv", sep = ";") # virome data but only in samples overlapping with samples from which 16S has been generated
bac.tax <- read.csv("/Bozzi_bacteria_taxonomy.csv")
bac.meta <- read.csv("/Bozzi_bacteria_metadata.csv", sep = ";") 

otus_in_otu_table <- colnames(bacteria)[-1]

# Filter taxonomy table
bac.tax <- bac.tax[bac.tax$Column1 %in% otus_in_otu_table, ]
```

### Relative abudance plot of 16S rRNA data
```{r}
Myco <- bac.tax %>%
  select(Column1,Family) %>%
  filter(Family == "Mycoplasmataceae") %>%
  .$Column1
Samples <- bacteria$SampleID

Ali <- bac.tax %>%
 select(Column1,Genus) %>%
 filter(Genus == "Aliivibrio") %>%
 .$Column1

Myco.df <- bacteria %>%
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

Ali.df <- bacteria %>%
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

all_OTUs <- bacteria %>%
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
  mutate(Other = Total - Mycoplasma - Aliivibrio) %>%
  select(-Total) %>%
  arrange(desc(Ali.df$Aliivibrio))

#Combine all data into one dataframe
combined_df <- combined_df %>%
  gather(key = "Genus", value = "Value", -SampleID)  

combined_df$SampleID <- factor(combined_df$SampleID, levels = unique(combined_df$SampleID))

# plot
plot_16S <-ggplot(combined_df, aes(x = SampleID, y = Value, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.3), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  labs(y = "Relative Abundance", title = "16S rRNA microbiota profile", x = "Samples") +
  scale_fill_manual(values = c("goldenrod1", "salmon2", "tomato4")) + 
  geom_tile() +
  theme(axis.line = element_blank(),
        panel.grid = element_blank())

plot_16S
```

### Phyloseq object generation:
```{r}
# Subset vOTU virome table based on vOTU IDs in size_virome and set vOTU as row names
virus_abun_cf <- virus_abun_cf %>%
  filter(vOTU %in% virus_length_cf$seq_name) %>%
 column_to_rownames(var = "vOTU") 

# Normalize TPM for virome data using contig length data
colnames(virus_length_cf) <- c("hgnc_symbol","transcript_length")
vContig.tmp <- NormalizeTPM(virus_abun_cf, tr_length = virus_length_cf, scale = 1e+06) 

# Extract vOTU IDs 
virome_vcontig <- rownames(vContig.tmp)
virus_abun_c_uhmf <- rownames(virus_abun_cf)
virus_abun_cf[is.na(virus_abun_cf)] <- 0
virus_abun_cf <- virus_abun_cf %>% as_tibble()

vContig.tmp[is.na(vContig.tmp)] <- 0

vContig.tmp <- vContig.tmp %>% as_tibble() 

# physeq object
physeq.tmp <- phyloseq(otu_table(vContig.tmp, taxa_are_rows = TRUE),
               tax_table(as.matrix(virus_tax_cf)),
               sample_data(metadata_cf))

#Non-normalised phyloseq object
physeq.nn <- phyloseq(otu_table((virus_abun_cf), taxa_are_rows = TRUE),
                             tax_table(as.matrix(virus_tax_cf)),
                             sample_data(metadata_cf))

subset_alii <- subset_taxa(physeq.nn, Host_family == "Vibrionaceae")

# Set up data for heatmap plots
gpt <- subset_samples(physeq.tmp, Myco90 %in% c("Mycoplasma","Aliivibrio sp.", "unknown"))

# Transformation to the OTU table 
otu_table(gpt) <- otu_table(gpt) + 1

# Subset the samples based on the microbiota column
gpt <- subset_samples(physeq.nn, Category %in% c("Mycoplasma sp.", "Aliivibrio sp.", "Sick fish, no 16S"))

otu_table(gpt) <- otu_table(gpt) + 1

# Reorder samples based on the 'Sample_name' column
sample_data <- data.frame(sample_data(gpt))
sample_order <- sample_data$Sample_name[order(sample_data$Category)]

taxa_names(gpt) <- tax_table(gpt)[, "vOTU"] 

# Generate heatmap plot 
p_heatmap <- plot_heatmap(gpt, method = "NMDS", distance = "bray", sample.label = "Sample_name", taxa.label = "vOTU")

# Reorder the samples 
p_heatmap$data$Sample <- factor(p_heatmap$data$Sample, levels = sample_order)

p_heatmap <- plot_heatmap(gpt, method = "NMDS", distance = "bray", sample.label = "Sample_name", taxa.label = "vOTU") +
  scale_fill_gradientn(colors = c("white", "dodgerblue2", "lightblue")) +  
  labs(y = "vOTU", x = "Samples") +
  theme(axis.text.y = element_text(size = 7, colour = "black"), 
        axis.text.x = element_text(angle = 45),  
        axis.ticks.x = element_blank(),  
        panel.spacing = unit(0, "lines"),
        axis.title.x = element_blank()) +
  ggforce::facet_row(~Category, scales = "free_x", space = "free") 
  
x_title_samples <- ggdraw() + 
  draw_label("Samples", size = 12, hjust = 0.5)

p_heatmap <- cowplot::plot_grid(
  p_heatmap, x_title_samples, 
  ncol = 1, 
  rel_heights = c(1, 0.05) )

p_heatmap

# vOTUs that had a Vibrionaceae host predicted by iPHoP
gpt <- subset_samples(subset_alii, Category %in% c("Mycoplasma sp.", "Aliivibrio sp.", "Sick fish, no 16S"))

# Transformation to the OTU table 
otu_table(gpt) <- otu_table(gpt) + 1

# Reorder the samples based on the 'Sample_name' column
sample_data <- data.frame(sample_data(gpt))
sample_order <- sample_data$Sample_name[order(sample_data$Category)]

taxa_names(gpt) <- tax_table(gpt)[, "vOTU"]  

p_heatmap <- plot_heatmap(gpt, method = "NMDS", distance = "bray", sample.label = "Sample_name", taxa.label = "vOTU")

# Reorder the samples 
p_heatmap$data$Sample <- factor(p_heatmap$data$Sample, levels = sample_order)

p_heatmap_ali <- plot_heatmap(gpt, method = "NMDS", distance = "bray", sample.label = "Sample_name", taxa.label = "vOTU") +
  scale_fill_gradientn(colors = c("white", "goldenrod1", "orange")) +  
  labs(y = "vOTU", x = "Samples") +
  theme(axis.text.y = element_text(size = 7, colour = "black"), 
        axis.text.x = element_text(angle = 45), 
        axis.ticks.x = element_blank(), 
        panel.spacing = unit(0, "lines"),
        axis.title.x = element_blank()) +
  ggforce::facet_row(~Category, scales = "free_x", space = "free") 
  
x_title_samples <- ggdraw() + 
  draw_label("Samples", size = 12, hjust = 0.5)

p_heatmap_ali <- cowplot::plot_grid(
  p_heatmap_ali, x_title_samples, 
  ncol = 1, 
  rel_heights = c(1, 0.05))

p_heatmap_ali
```

### Alpha diversity analysis between the three microbiota groups (Aliivibrio, Mycoplasma or no 16S)
```{r}
# Subset physeq object
physeq.mycoAlii.nn <- subset_samples(physeq.nn, Category %in% c("Mycoplasma sp.","Aliivibrio sp.", "Sick fish, no 16S"))

#Alpha diversity using Richness, Shannon and Simpson indices
richness <- estimate_richness(physeq.nn, measures = c("Observed","Shannon", "Simpson"))
md.richness <- metadata_cf[match(rownames(richness),metadata_cf$Sample_name),]
richness <- cbind(md.richness, richness)

y_labels <- c("vOTU richness", "Shannon index") 

plot_titles <- c("Observed", "Shannon")

# Loop stats and plots
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

    divtestdata$Myco90 <- factor(divtestdata$Category, levels = c("Aliivibrio sp.", "Mycoplasma sp.", "Sick fish, no 16S"))

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
    ggtitle(plot_titles[which(c("Observed", "Shannon") == i)]) +  
    xlab("") +
    ylab(y_labels[which(c("Observed", "Shannon") == i)]) +  
    theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90, hjust = 1))

    plot_list[[i]] <- plot}

# Adjust legend 
plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "none")

# Extract the legend
legend_b <- get_legend(
    plot_list[[2]] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "none"))

plot_list[[2]] <- plot_list[[2]] + theme(legend.position = "none") 

# Combine individual plots into single row
combined_plot <- cowplot::plot_grid(
    plot_list[[1]], plot_list[[2]], 
    nrow = 1)

# Combine the plots 
diversity_plot <- cowplot::plot_grid(combined_plot, 
    ncol = 1, 
    rel_heights = c(0.1, 1))

diversity_plot
```


### MAG module completion
```{r, print=F}
# Load data 
module_completion <- read.table(file='/metabolism-module_pathwise_completeness-MATRIX.txt', header = TRUE, sep = "\t")

# Load metabolic completion values across genomes
module_info <- read.csv("/modules_info.txt", sep = "\t")

# Load external genomes (pangenomes)
external_genomes <- read.table(file='/external-genomes.txt', header = TRUE, sep = "\t")

external_genomes <- external_genomes %>%
  dplyr::mutate(Group = c(rep("Mycoplasma", 21), rep("Aliivibrio", 11)))

# df generation
df <- melt(module_completion)
colnames(df) <- c('module', 'genome', 'completion')

# Associate genomes with groups 
df$group <- external_genomes$Group[match(df$genome, external_genomes$name)]
df$individual <- external_genomes$individual[match(df$genome, external_genomes$name)]

p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")}}

groups_order <- c("Mycoplasma", "Aliivibrio")
group_colors <- c("goldenrod1", "salmon2") 

# Wilcox test
p_value <- wilcox.test(df[df$group == "Mycoplasma", ]$completion, df[df$group == "Aliivibrio", ]$completion, exact = FALSE)$p.value

# Plot
plot_all_genomes <- ggplot(data=df, aes(x=group, y=completion, group=group)) +
    geom_violin(aes(fill=group), color="black", fill="grey", alpha=0.3) +
    geom_jitter(aes(shape=group, color=group), width=0.2, height=0.01, size=2.5, alpha=0.1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_shape_manual(values = c(16, 17)) +  
  annotate("text", x=1.5, y=1.03, label=paste(p_to_stars(p_value)), size=3.5) +
  geom_segment(aes(x=1, xend=2, y=1.025, yend=1.025), size=0.5) +  
  geom_segment(aes(x=1, xend=1, y=1.025, yend=1.02), size=0.5) +   
  geom_segment(aes(x=2, xend=2, y=1.025, yend=1.02), size=0.5) +   
  theme_classic() +
  theme(legend.position="none",
        axis.title.y = element_text(size=13, angle=90, vjust=0.5),
        axis.ticks.y = element_blank()) +
  ylab("Estimated KEGG module completion") +
  xlab("Pangenomes") +
  scale_x_discrete(limits = groups_order) +
  scale_color_manual(values = group_colors)


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
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=13, angle=90, vjust=0.5),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2,2,2,2), "cm")) +
    ylab(NULL) +
  xlab("Study Genomes") +
  labs(fill = "Genus", shape = "Genus", color = "Genus", size=13) +
   guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
  scale_x_discrete(limits = groups_order) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
   guides(shape = guide_legend(override.aes = list(alpha = 1, size = 5)), 
         color = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()

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
  ggplot(aes(x = group, y = subcategory, fill = completion)) + 
  geom_tile(color = "white", size = 0.02) +
  scale_fill_gradient(low = "white", high = "darkslategrey") +
  labs(y = "",
       x = "",
       fill = "Completion") +
  theme_minimal() +
  theme(axis.line = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12))


(plot_all_genomes + plot_CF_genomes) | MoI_plot + plot_layout(guides = "collect")
```

### Plot vOTU annotations
```{r}
# Data
data <- read.table("/Users/zts270/Desktop/01_CrappyFish_manuscript/data/move_to_dori/pharokka_cds_functions.tsv", 
                   header=TRUE, sep="\t")

# Reshape the data for heatmap
heatmap_data <- dcast(data, Description ~ vOTU, value.var = "Count", fill = 0)

# Set row names to be the descriptions
rownames(heatmap_data) <- heatmap_data$Description
heatmap_data <- heatmap_data[, -1]

heatmap_data[heatmap_data > 0] <- 1  # Convert to binary 

function_descriptions <- rownames(heatmap_data)

# Assign functions to broad categories
function_categories <- sapply(function_descriptions, function(x) {
  if (grepl("head|tail|connector", x, ignore.case = TRUE)) {
    return("Structural")
  } else if (grepl("DNA|RNA|nucleotide|transcription", x, ignore.case = TRUE)) {
    return("Phage Metabolism")
  } else if (grepl("lysis", x, ignore.case = TRUE)) {
    return("Lysis")
  } else if (grepl("integration and excision", x, ignore.case = TRUE)) {
    return("Temperate")
  } else if (grepl("tRNA|CRISPR|tmRNA|moron|VFDB|CARD", x, ignore.case = TRUE)) {
    return("Accessory Genes")
  } else if (grepl("other|CDS|unknown", x, ignore.case = TRUE)) {
    return("Unknown/CDS") 
  } 
})

# Create a df for row annotation
row_annotation <- data.frame(Category = factor(function_categories))
rownames(row_annotation) <- function_descriptions

# Aligned
row_annotation <- as.data.frame(row_annotation)
rownames(row_annotation) <- rownames(heatmap_data)

# Sort heatmap data by functional categories 
sorted_indices <- order(row_annotation$Category)
heatmap_data <- heatmap_data[sorted_indices, ]
row_annotation <- row_annotation[sorted_indices, , drop = FALSE]

# Extract the CDS counts
cds_counts <- aggregate(data$Count, by = list(vOTU = data$vOTU), sum)
colnames(cds_counts) <- c("vOTU", "CDS_count")

# Column annotation for CDS count
column_annotation <- data.frame(CDS_count = cds_counts$CDS_count)
rownames(column_annotation) <- cds_counts$vOTU

# Align 
column_annotation <- column_annotation[match(colnames(heatmap_data), rownames(column_annotation)), , drop = FALSE]

# Category colours
category_colors <- c(
  "Structural" = "deepskyblue3",
  "Phage Metabolism" = "gold2",
  "Lysis" = "orangered4",
  "Temperate" = "aquamarine3",
  "Accessory Genes" = "tan",
  "Unknown/CDS" = "darkorange3"
)

# Row annotation 
row_anno <- rowAnnotation(
  Category = factor(row_annotation$Category), 
  col = list(Category = category_colors) )

color_palette <- colorRamp2(c(0, 1), c("white", "darkslategrey"))

# Plot heatmap
anno_heatmap <- Heatmap(
  as.matrix(heatmap_data),
  name = "Presence/Absence",
  col = color_palette,
  top_annotation = HeatmapAnnotation(CDS_count = anno_barplot(cds_counts$CDS_count, 
                                                               gp = gpar(fill = "lightgray"), 
                                                               axis = TRUE)),
  right_annotation = row_anno,  
  cluster_rows = FALSE,          
  cluster_columns = FALSE,      
  show_row_names = TRUE,         
  show_column_names = TRUE,      
  row_names_gp = gpar(fontsize = 10),  
  column_names_gp = gpar(fontsize = 10) 
)

anno_heatmap
```
