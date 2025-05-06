
```{r load libraries}
library(tidyverse)
library(phyloseq)
library(readxl)
library(ggpubr)
library(rstatix)

```


```{r load data}
# metadata and vOTU table
metadata_cf <- read.csv("/metadata_cf.csv", sep = ";")  %>% filter(!Sample_name %in% c("S02", "H04"))
rownames(metadata_cf) <- metadata_cf$Sample_name %>% as.matrix()

virus_tax_cf <-  read_excel("/tax_vOTU_overview.xlsx", sheet = "Sheet1") 
virus_length_cf <- read.csv("/vOTU_length.tsv", sep = "\t")
virus_abun_cf <- read.csv("/rpkm_0.02_abundance_table.csv", sep = ";") %>%
  filter(vOTU %in% virus_length_cf$seq_name) %>%
 column_to_rownames(var = "vOTU") %>% as_tibble()


# CT data
CTs <- read.csv("/ct_plot_input_file.csv") %>% as.data.frame()
md <- read.csv("/Metadata_filtered.csv")

CTs <- CTs %>%
  full_join(md, by = "Sample") %>%
  na.omit()

CSW_df <- read.csv("/CSW_Bozzi_Samples.csv", sep=",", na.strings = "") %>% as.data.frame()

CSW_df <- CSW_df %>% 
  rename(Sample = Bozzi_et_al) %>%
  na.omit() %>%
  full_join(CTs, by = "Sample")
```

```{r data wrangling}
# phyloseq object
physeq <- phyloseq(otu_table((virus_abun_cf), taxa_are_rows = TRUE),
                             tax_table(as.matrix(virus_tax_cf)),
                             sample_data(metadata_cf))

# Calculate richness per sample (total number of vOTUs)
richness <- estimate_richness(physeq, measures = "Observed")  
richness$Sample_names <- rownames(richness)  
merged_data <- merge(richness, CSW_df, by = "Sample_names")  %>% filter(!OTU %in% "contaminants")
```

```{r plot regression}
# Plot with Spearman correlation
ggplot(merged_data, aes(x = Observed, y = Ct, color = Myco90)) +  
  geom_point() +
  geom_smooth(method = "lm") +  # Add regression line
  labs(x = "vOTU Richness", y = "Ct Value", title = "Richness vs. Ct with Regression Lines") +
  theme_minimal() +
  facet_wrap(~Myco90, scales = "free") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = Inf, label.y = Inf, 
           hjust = 1.2, vjust = 2,  # Align to upper right corner
           method = "spearman", color = "black")  # Use Spearman correlation

aliivibrio_data <- merged_data %>% filter(Myco90 == "Aliivibrio")

ct_vs_vOTU <- ggplot(aliivibrio_data, aes(x = Observed, y = Ct, color = Myco90)) +  
  geom_point(color = "goldenrod1") +
  geom_smooth(method = "lm", color = "goldenrod1") +  # Add regression line
  labs(x = "vOTU Richness", y = "Ct Value", title = "Bacterial biomass and vOTU richness") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = Inf, label.y = Inf, 
           hjust = 1.15, vjust = 1.8, 
           method = "spearman", color = "black")

plot(ct_vs_vOTU)
```


```{r box plot sample difference}
stat.test <- CSW_df %>%
  filter(!is.na(Wolters_et_al),
         OTU == "OTU2") %>%
  dplyr::mutate(Group = ifelse(Proportion >0.9, "SSM", "SSA")) %>%
  wilcox_test(Ct~Group)

boxplotA <- CSW_df %>%
  filter(!is.na(Wolters_et_al),
         OTU == "OTU2") %>%
  dplyr::mutate(Group = ifelse(Proportion >0.9, "SSM", "SSA")) %>%
  ggboxplot(x = "Group", y = "Ct",
            fill = "Group") +
  scale_y_reverse() +
  scale_fill_manual(values = c("SSM" = "salmon2", "SSA" = "goldenrod1")) +
  stat_pvalue_manual(stat.test, label = "Wilcoxon Rank Sum Test: {p}", y.position = -10) +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 60, hjust =1)) +
  theme(axis.text.x.bottom = element_text()) +
  theme(legend.position = "right") +
  guides(fill=guide_legend(title="Microbiota profile")) +
  xlab("Microbiota profile")


stat.test <- CSW_df %>%
  filter(!is.na(Wolters_et_al),
         OTU == "OTU2") %>%
  dplyr::mutate(Health_state = ifelse(Health_state == "H", "Healthy", "Sick")) %>%
  wilcox_test(Ct~Health_state)

boxplotB <- CSW_df %>%
  filter(!is.na(Wolters_et_al),
         OTU == "OTU2") %>%
  dplyr::mutate(Health_state = factor(ifelse(Health_state == "H", "Healthy", "Sick"), levels = c("Sick", "Healthy"))) %>%
  ggboxplot(x = "Health_state", y = "Ct",
            fill = "Health_state") +
  scale_y_reverse() +
  scale_fill_manual(values = c("Healthy" = "dodgerblue3", "Sick" = "lightblue3")) +
  stat_pvalue_manual(stat.test, label = "Wilcoxon Rank Sum Test: {p}", y.position = -10) +
  theme(axis.text.x.bottom = element_text()) +
  theme(legend.position = "right") + 
  guides(fill=guide_legend(title="Fish phenotype")) +
  xlab("Fish phenotype")

boxesAB <- cowplot::plot_grid(boxplotA,boxplotB, nrow = 1)
boxesAB
```
