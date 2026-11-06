# Copy this file to config/paths.R and edit paths as needed. Alternatively,
# set AML_PROJECT_ROOT or AML_PATH_CONFIG in the shell environment.

PROJECT_ROOT <- normalizePath(
  Sys.getenv("AML_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)
DATA_ROOT <- file.path(PROJECT_ROOT, "data")
RESULTS_ROOT <- file.path(PROJECT_ROOT, "results")

# Bulk transcriptomics (GSE114868)
GSE114868_DATA_DIR <- file.path(DATA_ROOT, "GSE114868")
BULK_DEG_RESULTS_DIR <- file.path(RESULTS_ROOT, "01_bulk_deg")

# Script 02 writes the full WGCNA analysis here. Script 03 and the auxiliary
# soft-power script read the module-trait and soft-power subdirectories.
BULK_WGCNA_RESULTS_DIR <- file.path(RESULTS_ROOT, "02_bulk_wgcna")

# Single-cell transcriptomics (GSE130756)
SCRNA_AML_TXT_DIR <- file.path(DATA_ROOT, "GSE130756", "AML_txt")
SCRNA_NORMAL_TXT_DIR <- file.path(DATA_ROOT, "GSE130756", "normal_txt")
SCRNA_QC_DIR <- file.path(RESULTS_ROOT, "04_scRNA_preprocessing")
SCRNA_ANNOTATION_DIR <- file.path(RESULTS_ROOT, "05_scRNA_annotation")
SCRNA_ANNOTATION_REFERENCE_FILE <- file.path(
  DATA_ROOT,
  "annotation",
  "ACT_Annotation results_top 1.txt"
)
SCRNA_ANNOTATED_OBJECT_FILE <- file.path(
  SCRNA_ANNOTATION_DIR,
  "2.resolution_0.5_annotation_1.qs"
)
SCRNA_HDWGCNA_INPUT_FILE <- file.path(
  SCRNA_ANNOTATION_DIR,
  "2.resolution_0.5_annotation.qs"
)
SCRNA_CELL_PROPORTION_DIR <- file.path(RESULTS_ROOT, "06_scRNA_cell_proportions")
SCRNA_HDWGCNA_DIR <- file.path(RESULTS_ROOT, "07_scRNA_hdWGCNA")
SCRNA_CYTOTRACE_DIR <- file.path(RESULTS_ROOT, "08_scRNA_CytoTRACE2")
SCRNA_PSEUDOBULK_DIR <- file.path(RESULTS_ROOT, "09_scRNA_pseudobulk")
SCRNA_CELLCHAT_DIR <- file.path(RESULTS_ROOT, "10_scRNA_CellChat")

# Integrated core-gene analysis
CORE_GENE_RESULTS_DIR <- file.path(RESULTS_ROOT, "11_core_gene_integration")
CORE_BULK_DEG_FILE <- file.path(
  BULK_DEG_RESULTS_DIR,
  "11_limma_all_results_AML_vs_Healthy.tsv"
)
CORE_BULK_WGCNA_GENE_FILE <- file.path(
  BULK_WGCNA_RESULTS_DIR,
  "05_hub_genes",
  "01_gene_module_membership_GS_DEG.tsv"
)
CORE_SCRNA_DEG_FILE <- file.path(
  SCRNA_PSEUDOBULK_DIR,
  "Tab04_AML_malignant_vs_NormalBM_pseudobulk_DEG_all.tsv"
)
CORE_SCRNA_MODULE_FILE <- file.path(SCRNA_HDWGCNA_DIR, "Tab01_moduleGene.csv")
CORE_SCRNA_HUB_FILE <- file.path(SCRNA_HDWGCNA_DIR, "Tab02_hubGene.csv")
CORE_SCRNA_EXPRESSION_PCT_FILE <- file.path(
  DATA_ROOT,
  "derived_inputs",
  "scRNA_malignant_NormalBM_expression_pct.tsv"
)

# Perturbational drug screening
DRUG_CMAP_FILE <- file.path(DATA_ROOT, "drug_screening", "cMAPnew2.xlsx")
DRUG_L1000_FILE <- file.path(DATA_ROOT, "drug_screening", "L1000fwd.csv")
DRUG_PRIORITIZATION_DIR <- file.path(RESULTS_ROOT, "12_drug_prioritization")

# Directional transcriptional reversal
REVERSAL_CORE_SIGNATURE_FILE <- file.path(
  CORE_GENE_RESULTS_DIR,
  "02_Final_core_signature_v2.tsv"
)
REVERSAL_EXPRESSION_FILE <- file.path(DATA_ROOT, "drug_rnaseq", "expr_matrix.txt")
REVERSAL_METADATA_FILE <- file.path(DATA_ROOT, "drug_rnaseq", "sample_metadata.txt")
REVERSAL_RESULTS_DIR <- file.path(RESULTS_ROOT, "13_core_gene_reversal")
