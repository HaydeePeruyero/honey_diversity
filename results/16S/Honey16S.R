library(dada2); packageVersion("dada2")
library(ggplot2)
library(phyloseq)
library(vegan)
library(Biostrings)
library(dplyr)
library(tidyr)

path <- "/d/Users/hayde/Downloads/OneDrive_1_1-7-2025/samples/"  # Cambia a tu ruta
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# 3. Filtrar y recortar
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)
#4. Inferencia con DADA2
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#5. Crear tabla de secuencias
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#6. Asignar taxonomía (SILVA o Greengenes)
# Descarga primero el archivo (ej: SILVA v138)
# https://zenodo.org/record/4587955/files/silva_nr99_v138_train_set.fa.gz
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/fastq/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

#7. Verificar secuencias nuevas
# Busca si aparecen secuencias sin clasificar a nivel de género o especie
table(is.na(taxa[, "Genus"]))

# ¿Y si hay especies nuevas?    
uniquesToFasta(seqtab.nochim, fout="~/Documents/fastq/asv_sequences.fasta", ids=colnames(seqtab.nochim))
               
table(taxa[, "Genus"] == "Acetilactobacillus")

## phyloseq
# Asegúrate de que la tax_table tenga nombres correctos
rank_names(ps)  # Debería incluir "Phylum"



############ Class

#8. Crear objeto phyloseq (opcional para visualización)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),   tax_table(taxa))

# Agregamos el nivel de phylum a las OTUs
ps_class <- tax_glom(ps, taxrank = "Class")

# Transformamos a proporciones relativas por muestra
ps_class_rel <- transform_sample_counts(ps_class, function(x) x / sum(x))

# Extraer tabla como dataframe
df <- psmelt(ps_class_rel)  # convierte phyloseq a tidy dataframe

top10 <- df %>%
  group_by(Class) %>%
  summarize(Abundance = sum(Abundance)) %>%
  top_n(10, Abundance) %>%
  pull(Class)

# Filtrar top10 y renombrar el resto como "Other"
df$Class <- ifelse(df$Class %in% top10, as.character(df$Class), "Other")

ggplot(df, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Class composition per sample")
############ Family
#8. Crear objeto phyloseq (opcional para visualización)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),   tax_table(taxa))

# Agregamos el nivel de phylum a las OTUs
ps_family <- tax_glom(ps, taxrank = "Family")

# Transformamos a proporciones relativas por muestra
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))

# Extraer tabla como dataframe
df <- psmelt(ps_family_rel)  # convierte phyloseq a tidy dataframe
top10 <- df %>%
  group_by(Family) %>%
  summarize(Abundance = sum(Abundance)) %>%
  top_n(10, Abundance) %>%
  pull(Family)

# Filtrar top10 y renombrar el resto como "Other"
df$Family <- ifelse(df$Family %in% top10, as.character(df$Family), "Other")

ggplot(df, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Family composition per sample")


############ Genus
#8. Crear objeto phyloseq (opcional para visualización)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),   tax_table(taxa))

# Agregamos el nivel de phylum a las OTUs
ps_genus <- tax_glom(ps, taxrank = "Genus")

# Transformamos a proporciones relativas por muestra
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Extraer tabla como dataframe
df <- psmelt(ps_genus_rel)  # convierte phyloseq a tidy dataframe
top10 <- df %>%
  group_by(Genus) %>%
  summarize(Abundance = sum(Abundance)) %>%
  top_n(10, Abundance) %>%
  pull(Genus)

# Filtrar top10 y renombrar el resto como "Other"
df$Genus <- ifelse(df$Genus %in% top10, as.character(df$Genus), "Other")

ggplot(df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Genus composition per sample")


############

# Extraer tabla de OTUs (ASVs en filas)
otu <- as.data.frame(otu_table(ps))

if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}

otu <- as.data.frame(otu)
otu$ASV <- rownames(otu)

tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)

# Combinar usando la columna ASV
otu_taxa <- merge(otu, tax, by = "ASV")
nrow(otu_taxa)  # ¿sigue dando 0?

head(otu_taxa)
write.csv(otu_taxa, "otu_table_with_taxonomy.csv", row.names = FALSE)

####### ahora con green genes 2024
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/fastq/gg2_2024_09_toSpecies_trainset.fa.gz", multithread=TRUE)

#7. Verificar secuencias nuevas
# Busca si aparecen secuencias sin clasificar a nivel de género o especie
table(is.na(taxa[, "Genus"]))

# ¿Y si hay especies nuevas?    
uniquesToFasta(seqtab.nochim, fout="~/Documents/fastq/asv_sequences.fasta", ids=colnames(seqtab.nochim))

table(taxa[, "Genus"] == "Acetilactobacillus")

############ Family
#8. Crear objeto phyloseq (opcional para visualización)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),   tax_table(taxa))

# Agregamos el nivel de phylum a las OTUs
ps_family <- tax_glom(ps, taxrank = "Family")

# Transformamos a proporciones relativas por muestra
ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))

# Extraer tabla como dataframe
df <- psmelt(ps_family_rel)  # convierte phyloseq a tidy dataframe
top10 <- df %>%
  group_by(Family) %>%
  summarize(Abundance = sum(Abundance)) %>%
  top_n(10, Abundance) %>%
  pull(Family)

# Filtrar top10 y renombrar el resto como "Other"
df$Family <- ifelse(df$Family %in% top10, as.character(df$Family), "Other")

ggplot(df, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Family composition per sample")

############

# Extraer tabla de OTUs (ASVs en filas)
otu <- as.data.frame(otu_table(ps))

if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}

otu <- as.data.frame(otu)
otu$ASV <- rownames(otu)

tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)

# Combinar usando la columna ASV
otu_taxa <- merge(otu, tax, by = "ASV")
nrow(otu_taxa)  # ¿sigue dando 0?

head(otu_taxa)
write.csv(otu_taxa, "otu_table_with_taxonomy.csv", row.names = FALSE)
########## Generar fasta
# Instalar Biostrings si no está instalado
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

library(Biostrings)

# Crear los IDs: ASV1, ASV2, ...
otu_taxa$ID <- paste0("ASV", seq_len(nrow(otu_taxa)))

# Crear objeto DNAStringSet para escribir FASTA
asv_fasta <- DNAStringSet(otu_taxa$ASV)
names(asv_fasta) <- otu_taxa$ID

# Guardar archivo FASTA
writeXStringSet(asv_fasta, filepath = "~/Documents/fastq/ASVs_from_otu_taxa.fasta")

#######After blast 
blast <- read.table("D:/Users/hayde/Documents/R_sites/honey_diversity/data/blast_results.tsv", header = FALSE, sep = "\t")

blast <- as.data.frame(blast)

colnames(blast) <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)


# Asegurar tipos
blast$evalue <- as.numeric(blast$evalue)
blast$bitscore <- as.numeric(blast$bitscore)

best_hits <- blast %>%
  group_by(qseqid) %>%
  arrange(evalue, desc(bitscore)) %>%
  slice(1) %>%         # solo el mejor hit
  ungroup()

best_hits <- blast %>%
  group_by(qseqid) %>%
  arrange(desc(pident), desc(bitscore)) %>%
  slice(1) %>%         # solo el mejor hit
  ungroup()



library(readr)
otu_taxa_hits <- merge(otu_taxa, best_hits, by.x = "ID", by.y = "qseqid", all.x = TRUE)

  
  
otu_taxa_hits <- read.csv("D:/Users/hayde/Documents/R_sites/honey_diversity/data/otu_taxa_hits_16s.csv", stringsAsFactors = FALSE)

#colnames(otu_taxa_hits)[3:8] <- c("10.S_Melli", "12.S_Melli", "14.S_Scapto", 
                                  "15.S_Melli", "16.S_Melli", "4.S_Scapto")
otu_taxa_hits

otu_taxa_hits$genome <- sub("\\|.*", "", otu_taxa_hits$sseqid)
otu_taxa_hits <- otu_taxa_hits %>%
  mutate(genome = case_when(
    pident < 95 ~ "other, less than 95%",
    TRUE ~ genome  # si no, conserva el valor original
  ))


# Asegúrate de que los nombres de columnas sean strings
colnames(otu_taxa_hits) <- as.character(colnames(otu_taxa_hits))

# Seleccionar columnas de muestra (solo números)
sample_cols <- colnames(otu_taxa_hits)[grepl("^X\\d+$", colnames(otu_taxa_hits))]

sample_cols <- colnames(otu_taxa_hits)[3:8]

#Convertir a relativas por muestra (por columna)
otu_taxa_hits_rel <- otu_taxa_hits %>%
  mutate(across(all_of(sample_cols), ~ .x / sum(.x, na.rm = TRUE)))

# Pivotear a formato largo
otu_long <- otu_taxa_hits_rel %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  )

ggplot(otu_long, aes(x = Sample, y = Abundance, fill = genome)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=c("#4F81BD", "gray"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Genome Best Hits per Sample")+
  theme( legend.text = element_text(size = 16, face = "italic"), 
         legend.key.size = unit(0.4, 'cm'), 
         axis.title = element_text(size = 16), 
         axis.text.x = element_text(size = 16, angle = 90),
         axis.text.y = element_text(size = 16), 
         plot.title = element_text(size = 20), 
         axis.title.x = element_text(size = 16, face = "bold"),
         axis.title.y = element_text(size = 16, face = "bold"), 
         legend.title = element_text(size = 16),
         legend.position = "bottom")

table(otu_taxa_hits$genome)

#### for each qseqid, diference:
blast_subset <- blast %>%
  filter(sseqid %in% c("Melli_45_016-contigs|k141_1186275",
                       "Melli_41_001-contigs|k141_1172984"))

blast_best <- blast_subset %>%
  group_by(qseqid, sseqid) %>%
  slice_max(pident, n = 1) %>%
  ungroup()

blast_wide <- blast_best %>%
  select(qseqid, sseqid, pident) %>%
  pivot_wider(names_from = sseqid, values_from = pident)

blast_wide_filtered <- blast_wide %>%
  filter(
    `Melli_45_016-contigs|k141_1186275` > 95,
    `Melli_41_001-contigs|k141_1172984` > 95
  )


blast_wide_filtered <- blast_wide_filtered %>%
  mutate(pident_diff = `Melli_45_016-contigs|k141_1186275` - `Melli_41_001-contigs|k141_1172984`)

####### mismatches
blast_subset <- blast %>%
  filter(sseqid %in% c("Melli_45_016-contigs|k141_1186275",
                       "Melli_41_001-contigs|k141_1172984"))

blast_best <- blast_subset %>%
  group_by(qseqid, sseqid) %>%
  slice_max(mismatch, n = 1, with_ties = FALSE) %>%
  ungroup()

blast_wide <- blast_best %>%
  select(qseqid, sseqid, mismatch) %>%
  pivot_wider(names_from = sseqid, values_from = mismatch)

blast_wide_filtered <- blast_wide %>%
  filter(
    `Melli_45_016-contigs|k141_1186275` < 6,
    `Melli_41_001-contigs|k141_1172984` < 6
  )


blast_wide_filtered <- blast_wide_filtered %>%
  mutate(mismatch_diff = `Melli_45_016-contigs|k141_1186275` - `Melli_41_001-contigs|k141_1172984`)

blast_wide_filtered
###########
target_genomes <- c("Melli_45_016-contigs|k141_1186275", 
                    "Melli_41_001-contigs|k141_1172984")

blast_filtered <- blast %>%
  filter(sseqid %in% target_genomes, pident > 95)

summary_stats <- blast_filtered %>%
  group_by(sseqid) %>%
  summarise(
    mean_pident = mean(pident, na.rm = TRUE),
    sd_pident = sd(pident, na.rm = TRUE),
    
    mean_mismatch = mean(mismatch, na.rm = TRUE),
    sd_mismatch = sd(mismatch, na.rm = TRUE),
    
    mean_length = mean(length, na.rm = TRUE),
    sd_length = sd(length, na.rm = TRUE),
    
    total_hits = n()
  )
print(summary_stats)
