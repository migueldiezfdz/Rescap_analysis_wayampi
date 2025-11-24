# ============================================================================
# LOAD REQUIRED PACKAGES
# ============================================================================

# Core data manipulation package suite (includes dplyr, tidyr, readr, etc.)
library("tidyverse")

# Ecological analysis package - used for rarefaction (rrarefy function)
library("vegan")

# Parallel processing packages for computationally intensive rarefaction
library('foreach')      # Enables parallel loops with %dopar%
library('doParallel')   # Backend for parallel processing with foreach


# ============================================================================
# DATA IMPORT
# ============================================================================

# Load main database table containing template matching results
# This table has sequence alignment statistics from a bioinformatics pipeline
Tabla_inicial <- read_tsv("/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/all_db_22_fulldb.tsv", 
                          col_names = c("template","score","expected","template_length",
                                        "template_identity","template_coverage",
                                        "query_identity","query_coverage","depth",
                                        "q_value","p_value","sample"),
                          comment = "#") %>%
  na.omit() %>%  # Remove rows with any missing values
  separate(sample, c("dataset", "kk"), sep = "\\.") %>%  # Split sample name by "."
  select(-kk) %>%  # Remove temporary column
  filter(p_value < 0.05) %>%  # Keep only statistically significant matches
  filter(dataset != "BichatK12" & dataset != "BichatK58")  # Exclude specific datasets

# Load supplementary data tables
Tabla_equivalencias <- read_csv("/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/equivalencias_guayana.csv")
tabla_bowtie <- read_csv('/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/bowtie_1_lines.txt', col_names = F)
tabla_lines <- read_csv('/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/bowtie_names.txt', col_names = F)

# Sample ID mapping table
tabla_label <- read_tsv('/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/name_rosetta.tsv')

# Metadata tables (original and binary versions for antibiotic data)
tabla_metadata <- read_csv('/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/metadata_original_sin_titulo_antibioticos.csv')
tabla_metadata_bin <- read_csv('/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/metadata_original_sin_titulo_antibioticos_cerounos.csv')

# Database cross-references
ID_EMBL <- read_csv('/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/ID_EMBL_todo.csv')

# BacMet2 database: predicted biocide and metal resistance genes
Bac_met_pred <- read_tsv("/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/BacMet2_PRE.155512.mapping.txt")

# ResFinder database: antibiotic resistance gene phenotypes
Resfinder_phenotype <- read_tsv("/Users/miguel/MEGA/Documentos_dropbox/Datos_investigacion/guayana/tablas/phenotypes.txt")


# ============================================================================
# DATA CLEANING AND STANDARDIZATION
# ============================================================================

# Standardize dataset naming conventions
Tabla_inicial$dataset <- gsub("MTC", "HE", Tabla_inicial$dataset)  # Replace MTC prefix with HE
Tabla_inicial$dataset <- gsub("HE0", "HE", Tabla_inicial$dataset)  # Remove leading zero

# Merge metadata with sample labels
tabla_metadata <- rename(tabla_metadata, dataset = id)
tabla_metadata <- left_join(tabla_metadata, tabla_label, by = 'dataset') %>% 
  select(-dataset)

tabla_metadata_bin <- rename(tabla_metadata_bin, dataset = id)
tabla_metadata_bin <- left_join(tabla_metadata_bin, tabla_label, by = 'dataset') %>% 
  select(-dataset)

# Fix typo in ResFinder phenotype class
Resfinder_phenotype <- rename(Resfinder_phenotype, template = `Gene_accession no.`)
Resfinder_phenotype$Class <- gsub("Tettracycline", "Tetracycline", Resfinder_phenotype$Class)


# ============================================================================
# PROCESS BACMET2 DATABASE (BIOCIDE/METAL RESISTANCE)
# ============================================================================

# Clean GenBank IDs by removing version numbers
Bac_met_pred$GenBank_ID <- gsub("\\..*", "", Bac_met_pred$GenBank_ID)
Bac_met_pred <- Bac_met_pred %>% rename(ID = GenBank_ID)

# Separate rows where multiple compounds are listed together
Bac_met_pred <- Bac_met_pred %>% 
  separate_rows(Compound, sep = ", ") %>% 
  distinct()

# Special handling for SDS compound: add class information
Bac_met_pred <- Bac_met_pred %>% 
  filter(Compound == 'SDS') %>% 
  mutate(Compound = "Sodium Dodecyl Sulfate (SDS) [Class_compound: Organo-sulfate]") %>% 
  full_join(Bac_met_pred %>% filter(Compound != "SDS"))

# Clean compound names and extract class information
Bac_met_pred$Compound <- str_remove(Bac_met_pred$Compound, "and ")
Bac_met_pred$Compound <- str_remove(Bac_met_pred$Compound, "\\[class")
Bac_met_pred$Compound <- str_remove(Bac_met_pred$Compound, "\\]")
Bac_met_pred <- Bac_met_pred %>% 
  separate(Compound, c("Compound", "Class_compound"), sep = ' : ')

# Fix typos in compound classes
Bac_met_pred$Class_compound <- gsub("Diamindine", "Diamidine", Bac_met_pred$Class_compound)
Bac_met_pred$Class_compound <- gsub("Phenolic Compounds", "Phenolic compounds", Bac_met_pred$Class_compound)
Bac_met_pred$Class_compound <- gsub("Phenanthridine ", "Phenanthridine", Bac_met_pred$Class_compound)

# Assign "Metal" class to entries without a class
Bac_met_pred <- Bac_met_pred %>% 
  filter(is.na(Class_compound)) %>% 
  mutate(Class_compound = "Metal") %>% 
  full_join(Bac_met_pred %>% filter(!is.na(Class_compound)))

# Merge EMBL IDs and prepare template gene identifiers
Bac_met_pred <- full_join(ID_EMBL, Bac_met_pred, by = "ID")
Bac_met_pred <- rename(Bac_met_pred, template = EMBL) %>% 
  select(-GI_number)
Bac_met_pred$template_gene <- Bac_met_pred$template
Bac_met_pred <- Bac_met_pred %>% 
  unite(template_gene, template, Gene_name, remove = F)


# ============================================================================
# ANNOTATE MAIN TABLE WITH METADATA AND DATABASE CLASSIFICATIONS
# ============================================================================

# Add sample labels to main table
Tabla_inicial <- left_join(Tabla_inicial, tabla_label, by = 'dataset')

# Classify samples by origin (Europe vs. Guyana)
Tabla_inicial$origin <- "europa"
Tabla_inicial$origin[grepl("^HE", Tabla_inicial$dataset)] <- "guayana"

# Merge with equivalence table and create combined stage/country identifier
Tabla_inicial <- left_join(Tabla_inicial, Tabla_equivalencias) %>% 
  unite(together, c(Stage, Country), sep = "_", remove = FALSE)

# Remove treatment group samples
Tabla_inicial <- Tabla_inicial %>% filter(Treatment != 'Treat')

# Classify genes by database type based on template ID prefix
Tabla_inicial$database <- "ABR"  # Default: Antibiotic Resistance
Tabla_inicial$database[grepl("^ENA", Tabla_inicial$template)] <- "BAC"  # Biocide/Metal
Tabla_inicial$database[grepl("^BAC", Tabla_inicial$template)] <- "BAC"
Tabla_inicial$database[grepl("^MOB", Tabla_inicial$template)] <- "REL"  # Mobile genetic elements


# ============================================================================
# FILTER AND ANNOTATE BY DATABASE TYPE
# ============================================================================

# ANTIBIOTIC RESISTANCE GENES (ABR): Apply 90% identity threshold
Tabla_inicial_ABR_filt <- Tabla_inicial %>% 
  filter(database == "ABR") %>% 
  filter(query_identity >= 90.00) %>% 
  filter(template_identity >= 90.00)
Tabla_inicial_ABR_filt$template_gene <- Tabla_inicial_ABR_filt$template
Tabla_inicial_ABR_filt <- left_join(Tabla_inicial_ABR_filt, Resfinder_phenotype, by = 'template')

# BIOCIDE/METAL RESISTANCE GENES (BAC): Apply 95% identity threshold (more stringent)
Tabla_inicial_BAC_pred_filt <- Tabla_inicial %>% 
  filter(database == "BAC") %>% 
  filter(query_identity >= 90.00) %>% 
  filter(template_identity >= 90.00)

# Parse template IDs (format: prefix|ID|suffix)
Tabla_inicial_BAC_pred_filt <- rename(Tabla_inicial_BAC_pred_filt, template_full = template)
Tabla_inicial_BAC_pred_filt <- Tabla_inicial_BAC_pred_filt %>% 
  separate(template_full, c("kk", "template", "template.1"), sep = "\\|", remove = F) %>% 
  distinct() %>% 
  select(-kk, -template.1)

# Merge with BacMet annotations
Tabla_inicial_BAC_pred_filt <- inner_join(Tabla_inicial_BAC_pred_filt, Bac_met_pred, by = 'template')

# Classify as Biocide or Metal resistance
Tabla_inicial_BAC_pred_filt$BAC <- "Biocide"
Tabla_inicial_BAC_pred_filt$BAC[grepl("^Metal", Tabla_inicial_BAC_pred_filt$Class_compound)] <- "Metal"
Tabla_inicial_BAC_pred_filt <- Tabla_inicial_BAC_pred_filt %>% 
  unite("database_BAC", database, BAC, sep = "_", remove = FALSE) %>% 
  select(-BAC)

# MOBILE GENETIC ELEMENTS (REL): Apply 90% identity threshold
Tabla_inicial_REL <- Tabla_inicial %>% 
  filter(database == "REL") %>% 
  filter(query_identity >= 90.00) %>% 
  filter(template_identity >= 90.00)
# Extract gene name (remove everything after first underscore)
Tabla_inicial_REL$template_gene <- sub("_.*", " ", Tabla_inicial_REL$template)

# Save unfiltered table for reference
Tabla_inicial_prev <- Tabla_inicial

# Combine all filtered tables
Tabla_inicial <- bind_rows(Tabla_inicial_ABR_filt, Tabla_inicial_BAC_pred_filt, Tabla_inicial_REL) %>% 
  distinct()

# Create detailed information table
Tabla_inicial_info <- Tabla_inicial %>% distinct()


# ============================================================================
# AGGREGATE DATA BY GENE AND SAMPLE
# ============================================================================

# For each gene-sample combination, keep maximum depth and minimum p-value
Tabla_inicial <- Tabla_inicial %>% 
  group_by(template_gene, label, origin, database) %>% 
  summarise(depth = max(depth), 
            p_value = min(p_value), 
            q_value = max(q_value)) %>% 
  ungroup()

# Rename template_gene back to template for consistency
Tabla_inicial <- rename(Tabla_inicial, template = template_gene)

# Create unique identifier combining gene and sample
Tabla_inicial <- Tabla_inicial %>% 
  unite(unique_id, c(template, label), sep = "_", remove = FALSE)

# Extract metal/biocide classification for later use
Tabla_metal_biocide <- Tabla_inicial_info %>% 
  filter(database == "BAC") %>% 
  select(template_gene, database_BAC) %>% 
  distinct() %>% 
  rename(template = template_gene)


# ============================================================================
# Plot : Plot Rarefaction curves
# ============================================================================

otu_matrix <- Tabla_inicial %>%
  select(template, label, depth) %>%
  pivot_wider(names_from = label, values_from = depth, values_fill = 0) %>%
  column_to_rownames("template")
otu_matrix_int <- round(as.matrix(otu_matrix))
str(otu_matrix_int)
rarecurve(otu_matrix_int, step=20, sample=min(rowSums(otu_matrix)), col="blue", cex=0.6)



# ============================================================================
# RAREFACTION: NORMALIZE READ COUNTS ACROSS SAMPLES
# ============================================================================

# Convert data to wide format matrix (samples × genes)
Tabla_rare <- Tabla_inicial %>% 
  select(template, depth, label) %>%
  pivot_wider(names_from = template, values_from = depth, values_fill = 0) %>%  
  column_to_rownames("label")

# Ensure all values are integers
Tabla_rare[] <- lapply(Tabla_rare, as.integer)

# Determine rarefaction depth (minimum library size across samples)
(raremax <- min(rowSums(Tabla_rare)))

# Set up parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)  # Use all but one core
registerDoParallel(cl)

# Perform rarefaction 1000 times and average (reduces stochastic variation)
# This randomly subsamples each sample to the same depth (raremax)
raro <- foreach(i = 1:1000, .combine = '+') %dopar% {
  rarox <- vegan::rrarefy(Tabla_rare, raremax) 
  rarox
}

# Stop parallel cluster
stopCluster(cl)

# Calculate average across 1000 rarefactions
raro <- raro / 1000


# ============================================================================
# FINALIZE RAREFIED DATA
# ============================================================================

# Convert rarefied matrix back to long format
tabla_rarefaction <- as.data.frame(raro) %>% 
  rownames_to_column(var = "label") %>% 
  pivot_longer(cols = 2:ncol(.), names_to = 'template', values_to = 'depth')

# Remove genes with zero depth after rarefaction
tabla_rarefaction <- tabla_rarefaction %>% filter(depth > 0)

# Create unique IDs
tabla_rarefaction <- tabla_rarefaction %>% 
  unite(unique_id, c(template, label), sep = "_", remove = FALSE)

# Merge rarefied data with original metadata
Tabla_inicial_minimal <- Tabla_inicial %>% select(unique_id, origin, database)
Tabla_inicial <- left_join(tabla_rarefaction, Tabla_inicial_minimal, by = 'unique_id')


# ============================================================================
# RECONSTRUCT DETAILED ANNOTATION TABLE
# ============================================================================

# Extract gene annotation information from original detailed table
Tabla_inicial_info_old <- Tabla_inicial_info %>% 
  select(template, template_gene:NCBI_annotation) %>% 
  distinct()

# Annotate ABR genes
Tabla_inicial_ABR_info <- Tabla_inicial %>% 
  filter(database == "ABR") %>% 
  left_join(., Tabla_inicial_info_old, by = 'template')

# Annotate BAC genes (requires renaming for proper join)
Tabla_inicial_BAC_pred_filt_info <- Tabla_inicial %>% 
  filter(database == "BAC")
Tabla_inicial_BAC_pred_filt_info <- rename(Tabla_inicial_BAC_pred_filt_info, template_gene = template)
Tabla_inicial_BAC_pred_filt_info <- Tabla_inicial_BAC_pred_filt_info %>% 
  left_join(., Tabla_inicial_info_old, by = 'template_gene')

# Annotate REL genes
Tabla_inicial_REL_info <- Tabla_inicial %>% 
  filter(database == "REL")

# Combine all annotated tables into final comprehensive dataset
Tabla_inicial_info <- bind_rows(Tabla_inicial_ABR_info, 
                                Tabla_inicial_BAC_pred_filt_info, 
                                Tabla_inicial_REL_info) %>% 
  distinct()
