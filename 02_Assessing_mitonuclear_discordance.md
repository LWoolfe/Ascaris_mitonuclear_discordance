# ======================================================================
# MITONUCLEAR TREE COMPARISON
# Overall discordance and discordance within host/geographic groups
# ======================================================================
#
# Purpose:
#   Compare nuclear and mitochondrial neighbour-joining (NJ) trees
#   generated from the same 141 Ascaris samples.
#
# Analyses:
#   1. Confirm identical samples are represented in both trees.
#   2. Compare overall tree topology using Robinson-Foulds (RF) distance.
#   3. Quantify the number and proportion of shared internal splits.
#   4. Compare patristic distances using Pearson and Spearman correlation.
#   5. Test correspondence between distance matrices using a Mantel test.
#   6. Determine whether host, region and country groups form discrete
#      splits in each tree.
#   7. Quantify mitonuclear discordance within host, region and country
#      groups.
#
# Notes:
#   - Formal tree comparisons are performed on unrooted trees because
#     NJ trees are fundamentally unrooted.
#   - Rooting used for visualisation/tanglegram construction is therefore
#     not included in the formal topology comparison.
#
# ======================================================================


# ======================================================================
# 1. WORKING DIRECTORY AND PACKAGES
# ======================================================================

setwd(
  "S:/Research/Lauren_Ascaris_2/mtDNA_mitoz_mitos_analyses/arnoud_temp"
)

library(ape)
library(phangorn)
library(vegan)


# ======================================================================
# 2. INPUT FILES
# ======================================================================

NUCLEAR_TREE <- "Ascaris.nuclear.NJ.n141.100k.nwk"
MTDNA_TREE   <- "mtDNA_all12genes_v1.nwk"

# Change this line only if the metadata filename changes.
METADATA <- "tanglegram_metadata_LW.csv"


# Metadata column names
SAMPLE_COL  <- "sample_id"
HOST_COL    <- "host"
REGION_COL  <- "region"
COUNTRY_COL <- "country"


# ======================================================================
# 3. READ TREES AND METADATA
# ======================================================================

nuc <- read.tree(NUCLEAR_TREE)
mt  <- read.tree(MTDNA_TREE)

meta <- read.csv(
  METADATA,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("============================================================\n")
cat("INPUT DATA\n")
cat("============================================================\n")

cat("\nMetadata columns:\n")
print(names(meta))

cat("\nMetadata rows:", nrow(meta), "\n")
cat("Nuclear tree tips:", Ntip(nuc), "\n")
cat("mtDNA tree tips:", Ntip(mt), "\n")


# ======================================================================
# 4. CHECK REQUIRED METADATA COLUMNS
# ======================================================================

required_cols <- c(
  SAMPLE_COL,
  HOST_COL,
  REGION_COL,
  COUNTRY_COL
)

missing_cols <- setdiff(
  required_cols,
  names(meta)
)

if (length(missing_cols) > 0) {
  
  stop(
    paste(
      "Required metadata columns not found:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ======================================================================
# 5. CLEAN METADATA
# ======================================================================
#
# Excel/CSV files can contain:
#
#   "Pig"
#   "Pig "
#   " Pig"
#
# or non-breaking spaces which visually look identical.
#
# These would otherwise be interpreted by R as different groups.
#
# ======================================================================


# ----------------------------------------------------------------------
# Function for cleaning metadata text
# ----------------------------------------------------------------------

clean_metadata_text <- function(x) {
  
  x <- as.character(x)
  
  # Replace non-breaking spaces that can be introduced by Excel
  x <- gsub(
    "\u00A0",
    " ",
    x,
    fixed = TRUE
  )
  
  # Remove leading and trailing whitespace
  x <- trimws(x)
  
  # Collapse repeated internal whitespace
  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )
  
  # Convert empty strings to NA
  x[x == ""] <- NA_character_
  
  return(x)
}


# ----------------------------------------------------------------------
# Clean sample IDs
# ----------------------------------------------------------------------

meta[[SAMPLE_COL]] <- clean_metadata_text(
  meta[[SAMPLE_COL]]
)


# ----------------------------------------------------------------------
# Clean categorical metadata
# ----------------------------------------------------------------------

meta[[HOST_COL]] <- clean_metadata_text(
  meta[[HOST_COL]]
)

meta[[REGION_COL]] <- clean_metadata_text(
  meta[[REGION_COL]]
)

meta[[COUNTRY_COL]] <- clean_metadata_text(
  meta[[COUNTRY_COL]]
)


# ----------------------------------------------------------------------
# Check categories AFTER cleaning
# ----------------------------------------------------------------------

cat("\n============================================================\n")
cat("CLEANED METADATA CATEGORIES\n")
cat("============================================================\n")

cat("\nHost:\n")
print(
  table(
    meta[[HOST_COL]],
    useNA = "ifany"
  )
)

cat("\nRegion:\n")
print(
  table(
    meta[[REGION_COL]],
    useNA = "ifany"
  )
)

cat("\nCountry:\n")
print(
  table(
    meta[[COUNTRY_COL]],
    useNA = "ifany"
  )
)


# ======================================================================
# 6. CHECK FOR DUPLICATE SAMPLE IDs
# ======================================================================

if (anyDuplicated(meta[[SAMPLE_COL]]) > 0) {
  
  cat("\nDuplicated metadata sample IDs:\n")
  
  print(
    meta[[SAMPLE_COL]][
      duplicated(meta[[SAMPLE_COL]])
    ]
  )
  
  stop(
    "Duplicate sample IDs detected after metadata cleaning."
  )
}


# ======================================================================
# 7. CHECK METADATA AGAINST TREE TIP LABELS
# ======================================================================

cat("\n============================================================\n")
cat("SAMPLE ID CHECKS\n")
cat("============================================================\n")


cat("\nSamples in nuclear tree but absent from metadata:\n")

print(
  setdiff(
    nuc$tip.label,
    meta[[SAMPLE_COL]]
  )
)


cat("\nSamples in mtDNA tree but absent from metadata:\n")

print(
  setdiff(
    mt$tip.label,
    meta[[SAMPLE_COL]]
  )
)


cat("\nMetadata samples absent from both trees:\n")

print(
  setdiff(
    meta[[SAMPLE_COL]],
    union(
      nuc$tip.label,
      mt$tip.label
    )
  )
)


# ======================================================================
# 8. VERIFY THAT BOTH TREES CONTAIN IDENTICAL SAMPLES
# ======================================================================

if (!setequal(
  nuc$tip.label,
  mt$tip.label
)) {
  
  cat("\nNuclear-only tree tips:\n")
  
  print(
    setdiff(
      nuc$tip.label,
      mt$tip.label
    )
  )
  
  cat("\nmtDNA-only tree tips:\n")
  
  print(
    setdiff(
      mt$tip.label,
      nuc$tip.label
    )
  )
  
  stop(
    "Nuclear and mitochondrial trees do not contain identical samples."
  )
}

cat(
  "\nPASS: Nuclear and mtDNA trees contain identical sample IDs.\n"
)


# ======================================================================
# 9. RESTRICT METADATA TO TREE SAMPLES
# ======================================================================

meta <- meta[
  meta[[SAMPLE_COL]] %in% nuc$tip.label,
]

cat(
  "\nMetadata records retained for analysis:",
  nrow(meta),
  "\n"
)

if (nrow(meta) != 141) {
  
  stop(
    paste(
      "Expected metadata for 141 samples but retained",
      nrow(meta)
    )
  )
}


# ----------------------------------------------------------------------
# Save cleaned metadata actually used for the analysis
# ----------------------------------------------------------------------

write.csv(
  meta,
  "tanglegram_metadata_n141_cleaned.csv",
  row.names = FALSE
)


# ======================================================================
# 10. UNROOT BOTH TREES
# ======================================================================
#
# NJ trees are fundamentally unrooted.
#
# Any rooting used for visualisation of the final trees/tanglegram
# is therefore ignored for formal topology comparison.
#
# ======================================================================

nuc_unroot <- unroot(nuc)
mt_unroot  <- unroot(mt)


cat("\n============================================================\n")
cat("UNROOTED TREE STRUCTURE\n")
cat("============================================================\n")

cat(
  "\nNuclear:",
  Ntip(nuc_unroot), "tips |",
  Nnode(nuc_unroot), "internal nodes |",
  "binary =", is.binary(nuc_unroot),
  "\n"
)

cat(
  "mtDNA:",
  Ntip(mt_unroot), "tips |",
  Nnode(mt_unroot), "internal nodes |",
  "binary =", is.binary(mt_unroot),
  "\n"
)


# Shared-split calculation below assumes fully bifurcating trees.

if (
  !is.binary(nuc_unroot) ||
  !is.binary(mt_unroot)
) {
  
  stop(
    paste(
      "At least one tree is not fully bifurcating.",
      "Shared-split calculations require binary trees."
    )
  )
}


# ======================================================================
# 11. OVERALL ROBINSON-FOULDS DISTANCE
# ======================================================================
#
# Robinson-Foulds distance compares the internal splits present in
# the two tree topologies.
#
# Branch lengths are ignored.
#
# Normalised RF:
#
#   0 = identical topology
#   1 = maximally different topology
#
# ======================================================================

rf_raw <- RF.dist(
  nuc_unroot,
  mt_unroot,
  normalize = FALSE,
  check.labels = TRUE,
  rooted = FALSE
)

rf_norm <- RF.dist(
  nuc_unroot,
  mt_unroot,
  normalize = TRUE,
  check.labels = TRUE,
  rooted = FALSE
)


# ----------------------------------------------------------------------
# Number of shared non-trivial internal splits
# ----------------------------------------------------------------------
#
# A fully bifurcating unrooted tree with n tips contains:
#
#                  n - 3
#
# non-trivial internal splits.
#
# For two binary trees containing the same tips:
#
# shared splits =
#
#       [2(n - 3) - RF] / 2
#
# ----------------------------------------------------------------------

n_tips <- Ntip(nuc_unroot)

n_splits <- n_tips - 3

shared_splits <- (
  (2 * n_splits) - rf_raw
) / 2

prop_shared <- (
  shared_splits / n_splits
)


cat("\n============================================================\n")
cat("OVERALL TOPOLOGICAL COMPARISON\n")
cat("============================================================\n")

cat(
  "\nRaw Robinson-Foulds distance:",
  rf_raw,
  "\n"
)

cat(
  "Normalized Robinson-Foulds distance:",
  rf_norm,
  "\n"
)

cat(
  "Non-trivial splits per tree:",
  n_splits,
  "\n"
)

cat(
  "Shared non-trivial splits:",
  shared_splits,
  "\n"
)

cat(
  "Proportion of splits shared:",
  prop_shared,
  "\n"
)


# ======================================================================
# 12. OVERALL PATRISTIC-DISTANCE MATRICES
# ======================================================================
#
# Patristic distance is the total branch length separating two tips
# within a tree.
#
# Rerooting does not alter pairwise patristic distances.
#
# ======================================================================

nuc_pat <- cophenetic.phylo(
  nuc_unroot
)

mt_pat <- cophenetic.phylo(
  mt_unroot
)


# ----------------------------------------------------------------------
# Put mtDNA matrix into exactly the same sample order as nuclear matrix
# ----------------------------------------------------------------------

mt_pat <- mt_pat[
  rownames(nuc_pat),
  colnames(nuc_pat)
]


stopifnot(
  identical(
    rownames(nuc_pat),
    rownames(mt_pat)
  ),
  identical(
    colnames(nuc_pat),
    colnames(mt_pat)
  )
)


# Convert matrices to distance objects

nuc_dist <- as.dist(
  nuc_pat
)

mt_dist <- as.dist(
  mt_pat
)


# ======================================================================
# 13. CORRELATION BETWEEN PATRISTIC DISTANCES
# ======================================================================

pearson_r <- cor(
  as.vector(nuc_dist),
  as.vector(mt_dist),
  method = "pearson"
)

spearman_rho <- cor(
  as.vector(nuc_dist),
  as.vector(mt_dist),
  method = "spearman"
)


cat("\n============================================================\n")
cat("OVERALL PATRISTIC-DISTANCE CORRELATION\n")
cat("============================================================\n")

cat(
  "\nPearson r:",
  pearson_r,
  "\n"
)

cat(
  "Spearman rho:",
  spearman_rho,
  "\n"
)


# ======================================================================
# 14. MANTEL TEST
# ======================================================================
#
# Tests whether the relationship between the nuclear and mitochondrial
# pairwise distance matrices is stronger than expected after random
# permutation of sample labels.
#
# ======================================================================

set.seed(20260825)

mantel_result <- mantel(
  nuc_dist,
  mt_dist,
  method = "pearson",
  permutations = 9999
)


cat("\n============================================================\n")
cat("MANTEL TEST\n")
cat("============================================================\n")

cat(
  "\nMantel r:",
  unname(mantel_result$statistic),
  "\n"
)

cat(
  "Mantel P:",
  mantel_result$signif,
  "\n"
)


# ======================================================================
# 15. SAVE OVERALL RESULTS
# ======================================================================

overall_results <- data.frame(
  
  comparison = c(
    "Raw Robinson-Foulds distance",
    "Normalized Robinson-Foulds distance",
    "Non-trivial splits per tree",
    "Shared non-trivial splits",
    "Proportion of splits shared",
    "Patristic Pearson correlation",
    "Patristic Spearman correlation",
    "Mantel statistic",
    "Mantel P-value"
  ),
  
  value = c(
    rf_raw,
    rf_norm,
    n_splits,
    shared_splits,
    prop_shared,
    pearson_r,
    spearman_rho,
    unname(mantel_result$statistic),
    mantel_result$signif
  )
)


write.csv(
  overall_results,
  "mitonuclear_tree_comparison_overall.csv",
  row.names = FALSE
)


cat("\nOverall results:\n")
print(overall_results)


# ======================================================================
# 16. FUNCTION: DOES A METADATA GROUP FORM A DISTINCT TREE SPLIT?
# ======================================================================
#
# This asks whether all members of a particular metadata group can
# be separated from the remainder of the tree by cutting one branch.
#
# Example:
#
#   Do all pig-derived samples form a single distinct grouping?
#
# To use is.monophyletic(), the unrooted tree is temporarily rooted
# using one sample outside the focal group.
#
# This temporary root has NO biological interpretation.
#
# ======================================================================

forms_split <- function(
    tree,
    tips
) {
  
  tips <- intersect(
    tips,
    tree$tip.label
  )
  
  outside <- setdiff(
    tree$tip.label,
    tips
  )
  
  
  # Groups with fewer than two samples cannot be meaningfully tested
  
  if (length(tips) < 2) {
    
    return(NA)
  }
  
  
  # Cannot test if the focal group contains the entire tree
  
  if (length(outside) < 1) {
    
    return(NA)
  }
  
  
  # Temporarily root using one individual outside the focal group
  
  tr <- root(
    tree,
    outgroup = outside[1],
    resolve.root = TRUE
  )
  
  
  return(
    is.monophyletic(
      tr,
      tips
    )
  )
}


# ======================================================================
# 17. FUNCTION: TEST ALL GROUPS WITHIN A METADATA VARIABLE
# ======================================================================

test_group_splits <- function(
    variable
) {
  
  groups <- sort(
    unique(
      meta[[variable]][
        !is.na(meta[[variable]])
      ]
    )
  )
  
  
  out <- lapply(
    groups,
    function(g) {
      
      
      # Samples belonging to focal group
      
      group_rows <- (
        !is.na(meta[[variable]]) &
          meta[[variable]] == g
      )
      
      tips <- meta[[SAMPLE_COL]][
        group_rows
      ]
      
      
      # Retain samples present in the nuclear tree
      
      tips <- intersect(
        tips,
        nuc_unroot$tip.label
      )
      
      
      data.frame(
        
        Variable = variable,
        
        Group = g,
        
        n = length(tips),
        
        Nuclear_split = forms_split(
          nuc_unroot,
          tips
        ),
        
        mtDNA_split = forms_split(
          mt_unroot,
          tips
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  return(
    do.call(
      rbind,
      out
    )
  )
}


# ======================================================================
# 18. TEST MAJOR HOST, REGION AND COUNTRY GROUPINGS
# ======================================================================

host_splits <- test_group_splits(
  HOST_COL
)

region_splits <- test_group_splits(
  REGION_COL
)

country_splits <- test_group_splits(
  COUNTRY_COL
)


cat("\n============================================================\n")
cat("HOST GROUPINGS\n")
cat("============================================================\n")

print(
  host_splits
)


cat("\n============================================================\n")
cat("REGIONAL GROUPINGS\n")
cat("============================================================\n")

print(
  region_splits
)


cat("\n============================================================\n")
cat("COUNTRY GROUPINGS\n")
cat("============================================================\n")

print(
  country_splits
)


# ======================================================================
# 19. FUNCTION: COMPARE NUCLEAR AND mtDNA TREES WITHIN EACH GROUP
# ======================================================================
#
# Each focal host/region/country is isolated by pruning both trees to
# contain only samples belonging to that group.
#
# The nuclear and mtDNA subtrees are then compared using:
#
#   - Raw RF distance
#   - Normalised RF distance
#   - Number of shared internal splits
#   - Proportion of internal splits shared
#   - Pearson correlation of patristic distances
#   - Spearman correlation of patristic distances
#
# At least four tips are required for an informative unrooted RF
# comparison because an unrooted tree with fewer than four tips
# contains no non-trivial internal split.
#
# ======================================================================

compare_within_groups <- function(
    variable
) {
  
  groups <- sort(
    unique(
      meta[[variable]][
        !is.na(meta[[variable]])
      ]
    )
  )
  
  
  out <- lapply(
    groups,
    function(g) {
      
      
      # ------------------------------------------------------------
      # Identify samples belonging to this group
      # ------------------------------------------------------------
      
      group_rows <- (
        !is.na(meta[[variable]]) &
          meta[[variable]] == g
      )
      
      tips <- meta[[SAMPLE_COL]][
        group_rows
      ]
      
      
      # Retain only tips present in BOTH trees
      
      tips <- intersect(
        tips,
        intersect(
          nuc_unroot$tip.label,
          mt_unroot$tip.label
        )
      )
      
      
      n <- length(tips)
      
      
      # ------------------------------------------------------------
      # Fewer than four samples:
      # no informative unrooted RF comparison
      # ------------------------------------------------------------
      
      if (n < 4) {
        
        return(
          
          data.frame(
            
            Variable = variable,
            
            Group = g,
            
            n = n,
            
            RF_raw = NA,
            
            RF_normalized = NA,
            
            Shared_splits = NA,
            
            Total_splits = NA,
            
            Proportion_shared = NA,
            
            Pearson_r = NA,
            
            Spearman_rho = NA,
            
            stringsAsFactors = FALSE
          )
        )
      }
      
      
      # ============================================================
      # Prune both trees to samples in the focal group
      # ============================================================
      
      nuc_sub <- keep.tip(
        nuc_unroot,
        tips
      )
      
      mt_sub <- keep.tip(
        mt_unroot,
        tips
      )
      
      
      # Treat both resulting trees as unrooted
      
      nuc_sub <- unroot(
        nuc_sub
      )
      
      mt_sub <- unroot(
        mt_sub
      )
      
      
      # ============================================================
      # Robinson-Foulds distance
      # ============================================================
      
      rf_raw_sub <- RF.dist(
        
        nuc_sub,
        
        mt_sub,
        
        normalize = FALSE,
        
        check.labels = TRUE,
        
        rooted = FALSE
      )
      
      
      rf_norm_sub <- RF.dist(
        
        nuc_sub,
        
        mt_sub,
        
        normalize = TRUE,
        
        check.labels = TRUE,
        
        rooted = FALSE
      )
      
      
      # ============================================================
      # Shared internal splits
      # ============================================================
      
      total_splits <- n - 3
      
      
      if (
        is.binary(nuc_sub) &&
        is.binary(mt_sub)
      ) {
        
        shared <- (
          (2 * total_splits) -
            rf_raw_sub
        ) / 2
        
        
        prop_shared <- (
          shared /
            total_splits
        )
        
      } else {
        
        shared <- NA
        
        prop_shared <- NA
      }
      
      
      # ============================================================
      # Patristic-distance correspondence
      # ============================================================
      
      nuc_pat_sub <- cophenetic.phylo(
        nuc_sub
      )
      
      mt_pat_sub <- cophenetic.phylo(
        mt_sub
      )
      
      
      # Reorder mtDNA matrix to match nuclear matrix
      
      mt_pat_sub <- mt_pat_sub[
        rownames(nuc_pat_sub),
        colnames(nuc_pat_sub)
      ]
      
      
      nuc_vec <- as.vector(
        as.dist(
          nuc_pat_sub
        )
      )
      
      mt_vec <- as.vector(
        as.dist(
          mt_pat_sub
        )
      )
      
      
      pearson_sub <- cor(
        nuc_vec,
        mt_vec,
        method = "pearson"
      )
      
      
      spearman_sub <- cor(
        nuc_vec,
        mt_vec,
        method = "spearman"
      )
      
      
      # ============================================================
      # Return results
      # ============================================================
      
      data.frame(
        
        Variable = variable,
        
        Group = g,
        
        n = n,
        
        RF_raw = rf_raw_sub,
        
        RF_normalized = rf_norm_sub,
        
        Shared_splits = shared,
        
        Total_splits = total_splits,
        
        Proportion_shared = prop_shared,
        
        Pearson_r = pearson_sub,
        
        Spearman_rho = spearman_sub,
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  return(
    do.call(
      rbind,
      out
    )
  )
}


# ======================================================================
# 20. RUN WITHIN-GROUP COMPARISONS
# ======================================================================

host_comparison <- compare_within_groups(
  HOST_COL
)

region_comparison <- compare_within_groups(
  REGION_COL
)

country_comparison <- compare_within_groups(
  COUNTRY_COL
)


cat("\n============================================================\n")
cat("DISCORDANCE WITHIN HOST GROUPS\n")
cat("============================================================\n")

print(
  host_comparison
)


cat("\n============================================================\n")
cat("DISCORDANCE WITHIN REGIONS\n")
cat("============================================================\n")

print(
  region_comparison
)


cat("\n============================================================\n")
cat("DISCORDANCE WITHIN COUNTRIES\n")
cat("============================================================\n")

print(
  country_comparison
)


# ======================================================================
# 21. SAVE GROUPING RESULTS
# ======================================================================

write.csv(
  host_splits,
  "mitonuclear_major_splits_host.csv",
  row.names = FALSE
)

write.csv(
  region_splits,
  "mitonuclear_major_splits_region.csv",
  row.names = FALSE
)

write.csv(
  country_splits,
  "mitonuclear_major_splits_country.csv",
  row.names = FALSE
)


write.csv(
  host_comparison,
  "mitonuclear_discordance_within_host.csv",
  row.names = FALSE
)

write.csv(
  region_comparison,
  "mitonuclear_discordance_within_region.csv",
  row.names = FALSE
)

write.csv(
  country_comparison,
  "mitonuclear_discordance_within_country.csv",
  row.names = FALSE
)


# ======================================================================
# 22. SAVE SESSION INFORMATION
# ======================================================================
#
# Records the versions of R and packages used for reproducibility.
#
# ======================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  "mitonuclear_tree_comparison_sessionInfo.txt"
)


# ======================================================================
# 23. FINAL SUMMARY
# ======================================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("\nOverall normalized RF:", rf_norm, "\n")
cat(
  "Shared internal splits:",
  shared_splits,
  "of",
  n_splits,
  "\n"
)

cat(
  "Proportion of internal splits shared:",
  prop_shared,
  "\n"
)

cat(
  "Patristic Pearson r:",
  pearson_r,
  "\n"
)

cat(
  "Patristic Spearman rho:",
  spearman_rho,
  "\n"
)

cat(
  "Mantel r:",
  unname(mantel_result$statistic),
  "\n"
)

cat(
  "Mantel P:",
  mantel_result$signif,
  "\n"
)

cat("\nOutput files written successfully.\n")
cat("============================================================\n")


# ============================================================
# Plot mitonuclear discordance by population group
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# Read analysis outputs
# ------------------------------------------------------------

overall <- read.csv(
  "mitonuclear_tree_comparison_overall.csv",
  stringsAsFactors = FALSE
)

host <- read.csv(
  "mitonuclear_discordance_within_host.csv",
  stringsAsFactors = FALSE
)

region <- read.csv(
  "mitonuclear_discordance_within_region.csv",
  stringsAsFactors = FALSE
)

country <- read.csv(
  "mitonuclear_discordance_within_country.csv",
  stringsAsFactors = FALSE
)


# ============================================================
# Extract overall values
# ============================================================

overall_rf <- overall$value[
  overall$comparison == "Normalized Robinson-Foulds distance"
]

overall_spearman <- overall$value[
  overall$comparison == "Patristic Spearman correlation"
]


overall_plot <- data.frame(
  Category = "Overall",
  Group = "All samples",
  n = 141,
  RF_normalized = overall_rf,
  Spearman_rho = overall_spearman
)


# ============================================================
# Select main host groups
# ============================================================

host_plot <- host %>%
  filter(Group %in% c("Human", "Pig")) %>%
  transmute(
    Category = "Host",
    Group = Group,
    n = n,
    RF_normalized = RF_normalized,
    Spearman_rho = Spearman_rho
  )


# ============================================================
# Select main regional group
# ============================================================
#
# Africa is retained in the main figure because it contains
# 123 samples. Asia and Europe contain considerably fewer
# samples and are retained in the complete results table.
#
# ============================================================

region_plot <- region %>%
  filter(Group == "Africa") %>%
  transmute(
    Category = "Region",
    Group = Group,
    n = n,
    RF_normalized = RF_normalized,
    Spearman_rho = Spearman_rho
  )


# ============================================================
# Select well-sampled country populations
# ============================================================
#
# Ethiopia and Kenya were the two country-level groups with
# sufficiently large sample numbers for meaningful comparison.
#
# ============================================================

country_plot <- country %>%
  filter(Group %in% c("Ethiopia", "Kenya")) %>%
  transmute(
    Category = "Country",
    Group = Group,
    n = n,
    RF_normalized = RF_normalized,
    Spearman_rho = Spearman_rho
  )


# ============================================================
# Combine groups
# ============================================================

plot_data <- bind_rows(
  overall_plot,
  host_plot,
  region_plot,
  country_plot
)


# Add sample size to labels

plot_data <- plot_data %>%
  mutate(
    Group_label = paste0(Group, " (n = ", n, ")")
  )


# Set plotting order

group_order <- c(
  "All samples (n = 141)",
  "Human (n = 124)",
  "Pig (n = 17)",
  "Africa (n = 123)",
  "Ethiopia (n = 56)",
  "Kenya (n = 60)"
)

plot_data$Group_label <- factor(
  plot_data$Group_label,
  levels = rev(group_order)
)


# ============================================================
# Convert to long format
# ============================================================

plot_long <- plot_data %>%
  select(
    Group_label,
    RF_normalized,
    Spearman_rho
  ) %>%
  pivot_longer(
    cols = c(
      RF_normalized,
      Spearman_rho
    ),
    names_to = "Metric",
    values_to = "Value"
  )


# Give panels clearer names

plot_long$Metric <- factor(
  plot_long$Metric,
  levels = c(
    "RF_normalized",
    "Spearman_rho"
  ),
  labels = c(
    "Normalized RF\nHigher = greater topological discordance",
    "Spearman \u03c1\nHigher = greater concordance"
  )
)


# ============================================================
# Create dot plot
# ============================================================

p <- ggplot(
  plot_long,
  aes(
    x = Value,
    y = Group_label
  )
) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", Value)
    ),
    nudge_x = 0.035,
    size = 3.5,
    hjust = 0
  ) +
  facet_wrap(
    ~ Metric,
    nrow = 1
  ) +
  scale_x_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Value",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      size = 10
    ),
    plot.margin = margin(
      10, 25, 10, 10
    )
  )

print(p)



# ======================================================================
# 24. PLOT MITONUCLEAR DISCORDANCE BY POPULATION GROUP
# ======================================================================

library(ggplot2)
library(dplyr)
library(tidyr)


# ----------------------------------------------------------------------
# Read analysis outputs
# ----------------------------------------------------------------------

overall <- read.csv(
  "mitonuclear_tree_comparison_overall.csv",
  stringsAsFactors = FALSE
)

host <- read.csv(
  "mitonuclear_discordance_within_host.csv",
  stringsAsFactors = FALSE
)

region <- read.csv(
  "mitonuclear_discordance_within_region.csv",
  stringsAsFactors = FALSE
)

country <- read.csv(
  "mitonuclear_discordance_within_country.csv",
  stringsAsFactors = FALSE
)


# ======================================================================
# 25. EXTRACT OVERALL VALUES
# ======================================================================

overall_rf <- overall$value[
  overall$comparison == "Normalized Robinson-Foulds distance"
]

overall_spearman <- overall$value[
  overall$comparison == "Patristic Spearman correlation"
]


overall_plot <- data.frame(
  Category = "Overall",
  Group = "All samples",
  n = 141,
  RF_normalized = overall_rf,
  Spearman_rho = overall_spearman
)


# ======================================================================
# 26. SELECT HOST GROUPS
# ======================================================================

host_plot <- host %>%
  filter(Group %in% c("Human", "Pig")) %>%
  transmute(
    Category = "Host",
    Group = Group,
    n = n,
    RF_normalized = RF_normalized,
    Spearman_rho = Spearman_rho
  )


# ======================================================================
# 27. SELECT REGIONAL GROUP
# ======================================================================

region_plot <- region %>%
  filter(Group == "Africa") %>%
  transmute(
    Category = "Region",
    Group = Group,
    n = n,
    RF_normalized = RF_normalized,
    Spearman_rho = Spearman_rho
  )


# ======================================================================
# 28. SELECT COUNTRY GROUPS
# ======================================================================

country_plot <- country %>%
  filter(Group %in% c("Ethiopia", "Kenya")) %>%
  transmute(
    Category = "Country",
    Group = Group,
    n = n,
    RF_normalized = RF_normalized,
    Spearman_rho = Spearman_rho
  )


# ======================================================================
# 29. COMBINE GROUPS
# ======================================================================

plot_data <- bind_rows(
  overall_plot,
  host_plot,
  region_plot,
  country_plot
)


plot_data <- plot_data %>%
  mutate(
    Group_label = paste0(
      Group,
      " (n = ",
      n,
      ")"
    )
  )


group_order <- c(
  "All samples (n = 141)",
  "Human (n = 124)",
  "Pig (n = 17)",
  "Africa (n = 123)",
  "Ethiopia (n = 56)",
  "Kenya (n = 60)"
)


plot_data$Group_label <- factor(
  plot_data$Group_label,
  levels = rev(group_order)
)


# ======================================================================
# 30. CONVERT TO LONG FORMAT
# ======================================================================

plot_long <- plot_data %>%
  select(
    Group_label,
    RF_normalized,
    Spearman_rho
  ) %>%
  pivot_longer(
    cols = c(
      RF_normalized,
      Spearman_rho
    ),
    names_to = "Metric",
    values_to = "Value"
  )


plot_long$Metric <- factor(
  plot_long$Metric,
  levels = c(
    "RF_normalized",
    "Spearman_rho"
  ),
  labels = c(
    "Normalized RF\nHigher = greater topological discordance",
    "Spearman \u03c1\nHigher = greater concordance"
  )
)


# ======================================================================
# 31. CREATE MAIN DISCORDANCE PLOT
# ======================================================================

p <- ggplot(
  plot_long,
  aes(
    x = Value,
    y = Group_label
  )
) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Value
      )
    ),
    nudge_x = 0.035,
    size = 3.5,
    hjust = 0
  ) +
  facet_wrap(
    ~ Metric,
    nrow = 1
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1.05
    ),
    breaks = seq(
      0,
      1,
      0.2
    )
  ) +
  labs(
    x = "Value",
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      size = 10
    ),
    plot.margin = margin(
      10,
      25,
      10,
      10
    )
  )


print(p)


ggsave(
  "mitonuclear_discordance_by_group.png",
  plot = p,
  width = 9,
  height = 5,
  dpi = 300
)


# ======================================================================
# 32. PREPARE METADATA FOR ETHIOPIA CLADE ANALYSIS
# ======================================================================
#
# The metadata currently contains an unnamed column.
#
# Do NOT delete this column yet because it may contain mitochondrial
# clade membership.
#
# Instead, give unnamed columns temporary valid names so that dplyr
# can operate on the dataframe.
#
# ======================================================================

bad_names <- (
  is.na(names(meta)) |
    trimws(names(meta)) == ""
)


if (any(bad_names)) {
  
  number_bad <- sum(bad_names)
  
  names(meta)[bad_names] <- paste0(
    "unnamed_",
    seq_len(number_bad)
  )
}


cat("\n============================================================\n")
cat("METADATA COLUMNS AFTER FIXING UNNAMED COLUMNS\n")
cat("============================================================\n")

print(
  names(meta)
)


# ======================================================================
# 33. ETHIOPIA: TEST MITONUCLEAR DISCORDANCE USING KNOWN
#     MITOCHONDRIAL CLADE ASSIGNMENTS
# ======================================================================
#
# The original mitochondrial analysis assigned 61 Ethiopian samples to:
#
#   Clade D = 51
#   Clade A =  8
#   Clade B =  2
#
# Rather than inferring mitochondrial groups using hierarchical
# clustering, the known clade assignments are used directly here.
#
# This allows us to:
#
#   1. Determine exactly which Ethiopian samples survived into the
#      matched nuclear-mitochondrial dataset.
#
#   2. Determine whether either of the two Clade B samples was lost.
#
#   3. Test whether the previously observed two-cluster mitochondrial
#      division represented Clade D versus Clades A/B.
#
#   4. Compare nuclear and mitochondrial distances separately for:
#
#         Within A
#         Within B
#         Within D
#         A vs B
#         A vs D
#         B vs D
#
#   5. Test the broader contrast:
#
#         Clade D versus non-D (A + B)
#
#      This is useful because the main mitochondrial subdivision within
#      Ethiopia appears to separate the divergent D lineage from the
#      smaller A/B component.
#
# ======================================================================


# ======================================================================
# 34. MAKE SURE METADATA COLUMN NAMES ARE VALID
# ======================================================================
#
# The original metadata contained an unnamed column.
#
# We do not need that column for this analysis, but dplyr will fail if
# the dataframe still contains a blank column name. Therefore any blank
# column is renamed rather than deleted.
#
# ======================================================================

bad_names <- (
  is.na(names(meta)) |
    trimws(names(meta)) == ""
)

if (any(bad_names)) {
  
  names(meta)[bad_names] <- paste0(
    "unnamed_",
    seq_len(sum(bad_names))
  )
}


cat("\n============================================================\n")
cat("METADATA COLUMNS\n")
cat("============================================================\n")

print(names(meta))


# ======================================================================
# 35. ENTER THE ORIGINAL ETHIOPIAN MITOCHONDRIAL CLADE ASSIGNMENTS
# ======================================================================
#
# These are the 61 Ethiopian samples from the original mitochondrial
# dataset.
#
# Keeping this lookup inside the script means that the current analysis
# is reproducible and does not depend on the previous arbitrary k = 2
# clustering step.
#
# ======================================================================

eth_clade_lookup <- read.table(
  text = "
sample_id clade
SRR31631925 D
SRR31675196 D
SRR31675197 D
SRR31675198 B
SRR31675199 B
SRR31675200 D
SRR31675201 D
SRR31675202 A
SRR31675203 D
SRR31675204 D
SRR31692059 D
SRR31692061 D
SRR31692062 A
SRR31692063 D
SRR31692064 D
SRR31692065 D
SRR31692066 D
SRR31692067 D
SRR31692068 D
SRR31692069 D
SRR31692070 D
SRR31692071 D
SRR31692072 D
SRR31692073 D
SRR31692074 A
SRR31692075 D
SRR31692076 D
SRR31692077 D
SRR31692078 D
SRR31757516 A
SRR31757517 D
SRR31757518 D
SRR31757519 A
SRR31757520 A
SRR31757521 D
SRR31757522 D
SRR31757523 D
SRR31757524 D
SRR31757525 D
SRR31757526 D
SRR31757527 A
SRR31757528 D
SRR31757529 D
SRR31757530 D
SRR31757531 D
SRR31757532 D
SRR31757533 D
SRR31757534 D
SRR31757535 D
SRR31757536 D
SRR31757537 D
SRR31757538 D
SRR31757539 D
SRR31757540 D
SRR31757541 D
SRR31757542 D
SRR31757543 D
SRR31757544 D
SRR31757545 D
SRR31757546 A
SRR31757547 D
",
  header = TRUE,
  stringsAsFactors = FALSE
)


# ======================================================================
# 36. VERIFY THAT THE LOOKUP REPRODUCES THE ORIGINAL RESULTS
# ======================================================================
#
# This is an important error-checking step.
#
# The original dataset should contain:
#
#   A = 8
#   B = 2
#   D = 51
#   Total = 61
#
# If this check fails, the script stops rather than silently analysing
# an incorrectly transcribed clade table.
#
# ======================================================================

original_clade_counts <- table(
  factor(
    eth_clade_lookup$clade,
    levels = c("A", "B", "D")
  )
)


cat("\n============================================================\n")
cat("ORIGINAL ETHIOPIAN mtDNA CLADE COUNTS\n")
cat("============================================================\n")

print(original_clade_counts)

cat(
  "\nTotal Ethiopian samples:",
  nrow(eth_clade_lookup),
  "\n"
)


expected_counts <- c(
  A = 8,
  B = 2,
  D = 51
)


if (
  nrow(eth_clade_lookup) != 61 ||
  !all(
    as.integer(original_clade_counts) ==
    as.integer(expected_counts)
  )
) {
  
  stop(
    paste(
      "Original Ethiopian clade assignments do not match",
      "the expected 8 A, 2 B and 51 D samples."
    )
  )
}


cat(
  "\nPASS: Original Ethiopian dataset contains",
  "8 Clade A, 2 Clade B and 51 Clade D samples.\n"
)


# ======================================================================
# 37. DETERMINE WHICH ORIGINAL ETHIOPIAN SAMPLES SURVIVED
#     INTO THE MATCHED DATASET
# ======================================================================
#
# The mitochondrial dataset originally contained 61 Ethiopian samples.
#
# The mitonuclear analysis contains fewer Ethiopian samples because only
# individuals present in both the nuclear and mitochondrial analyses are
# retained.
#
# Here we explicitly check every original Ethiopian sample against:
#
#   - the cleaned metadata
#   - the nuclear tree
#   - the mitochondrial tree
#
# This tells us exactly which five samples were lost.
#
# ======================================================================

eth_status <- eth_clade_lookup


eth_status$in_metadata <- (
  eth_status$sample_id %in%
    meta[[SAMPLE_COL]]
)


eth_status$in_nuclear_tree <- (
  eth_status$sample_id %in%
    nuc_unroot$tip.label
)


eth_status$in_mtDNA_tree <- (
  eth_status$sample_id %in%
    mt_unroot$tip.label
)


eth_status$retained_matched <- (
  eth_status$in_metadata &
    eth_status$in_nuclear_tree &
    eth_status$in_mtDNA_tree
)


# ----------------------------------------------------------------------
# Add the country recorded in the matched metadata as a cross-check
# ----------------------------------------------------------------------

country_lookup <- setNames(
  meta[[COUNTRY_COL]],
  meta[[SAMPLE_COL]]
)


eth_status$metadata_country <- unname(
  country_lookup[
    eth_status$sample_id
  ]
)


# ----------------------------------------------------------------------
# Check that retained samples are actually labelled Ethiopia
# ----------------------------------------------------------------------

wrong_country <- eth_status[
  eth_status$retained_matched &
    !is.na(eth_status$metadata_country) &
    eth_status$metadata_country != "Ethiopia",
  ,
  drop = FALSE
]


if (nrow(wrong_country) > 0) {
  
  cat(
    "\nWARNING: The following samples are not labelled Ethiopia",
    "in the matched metadata:\n"
  )
  
  print(wrong_country)
}


# ======================================================================
# 38. REPORT RETAINED AND EXCLUDED CLADE COUNTS
# ======================================================================

matched_eth <- eth_status[
  eth_status$retained_matched,
  ,
  drop = FALSE
]


excluded_eth <- eth_status[
  !eth_status$retained_matched,
  ,
  drop = FALSE
]


matched_clade_counts <- table(
  factor(
    matched_eth$clade,
    levels = c("A", "B", "D")
  )
)


excluded_clade_counts <- table(
  factor(
    excluded_eth$clade,
    levels = c("A", "B", "D")
  )
)


cat("\n============================================================\n")
cat("ETHIOPIAN SAMPLES RETAINED FOR MITONUCLEAR ANALYSIS\n")
cat("============================================================\n")

print(matched_clade_counts)

cat(
  "\nTotal retained:",
  nrow(matched_eth),
  "\n"
)


cat("\n============================================================\n")
cat("ETHIOPIAN SAMPLES LOST DURING MATCHING / QC\n")
cat("============================================================\n")

print(
  excluded_eth[
    ,
    c(
      "sample_id",
      "clade",
      "in_metadata",
      "in_nuclear_tree",
      "in_mtDNA_tree"
    )
  ]
)


cat("\nExcluded clade counts:\n")

print(excluded_clade_counts)


# ======================================================================
# 39. SPECIFICALLY CHECK THE TWO CLADE B SAMPLES
# ======================================================================
#
# Original Clade B samples:
#
#   SRR31675198
#   SRR31675199
#
# This directly answers whether Clade B remained represented after
# matching to the nuclear dataset.
#
# ======================================================================

original_B <- eth_status[
  eth_status$clade == "B",
  ,
  drop = FALSE
]


retained_B <- original_B[
  original_B$retained_matched,
  ,
  drop = FALSE
]


excluded_B <- original_B[
  !original_B$retained_matched,
  ,
  drop = FALSE
]


cat("\n============================================================\n")
cat("CLADE B CHECK\n")
cat("============================================================\n")


cat("\nOriginal Clade B samples:\n")

print(
  original_B[
    ,
    c(
      "sample_id",
      "retained_matched"
    )
  ]
)


cat(
  "\nNumber of Clade B samples retained:",
  nrow(retained_B),
  "of",
  nrow(original_B),
  "\n"
)


if (nrow(retained_B) > 0) {
  
  cat("\nClade B samples retained:\n")
  
  print(
    retained_B$sample_id
  )
}


if (nrow(excluded_B) > 0) {
  
  cat("\nClade B samples excluded:\n")
  
  print(
    excluded_B$sample_id
  )
}


# ======================================================================
# 40. CREATE MATCHED ETHIOPIAN CLADE TABLE
# ======================================================================

eth_matched <- matched_eth[
  ,
  c(
    "sample_id",
    "clade"
  ),
  drop = FALSE
]


eth_tips <- eth_matched$sample_id


if (length(eth_tips) < 4) {
  
  stop(
    "Too few matched Ethiopian samples for analysis."
  )
}


cat("\n============================================================\n")
cat("MATCHED ETHIOPIAN DATASET\n")
cat("============================================================\n")

cat(
  "\nNumber of samples:",
  length(eth_tips),
  "\n"
)

print(
  table(
    eth_matched$clade
  )
)


# ======================================================================
# 41. EXTRACT ETHIOPIAN NUCLEAR AND mtDNA DISTANCE MATRICES
# ======================================================================
#
# The same Ethiopian individuals are extracted from the nuclear and
# mitochondrial patristic-distance matrices.
#
# This ensures that every pairwise comparison refers to exactly the
# same pair of worms in both genomic compartments.
#
# ======================================================================

eth_nuc_pat <- nuc_pat[
  eth_tips,
  eth_tips,
  drop = FALSE
]


eth_mt_pat <- mt_pat[
  eth_tips,
  eth_tips,
  drop = FALSE
]


# ======================================================================
# 42. REPRODUCE THE OVERALL ETHIOPIAN CORRELATION
# ======================================================================
#
# This should reproduce approximately:
#
#   Pearson r   = 0.060
#   Spearman rho = 0.039
#
# Any small difference would indicate that the samples being used here
# differ from those used in the previous country-level comparison.
#
# ======================================================================

eth_nuc_vec <- as.vector(
  as.dist(
    eth_nuc_pat
  )
)


eth_mt_vec <- as.vector(
  as.dist(
    eth_mt_pat
  )
)


eth_pearson <- cor(
  eth_nuc_vec,
  eth_mt_vec,
  method = "pearson"
)


eth_spearman <- cor(
  eth_nuc_vec,
  eth_mt_vec,
  method = "spearman"
)


cat("\n============================================================\n")
cat("OVERALL ETHIOPIAN PATRISTIC-DISTANCE CORRELATION\n")
cat("============================================================\n")


cat(
  "\nPearson r:",
  eth_pearson,
  "\n"
)


cat(
  "Spearman rho:",
  eth_spearman,
  "\n"
)


# ======================================================================
# 43. CHECK WHAT THE PREVIOUS k = 2 mtDNA CLUSTERING REPRESENTED
# ======================================================================
#
# Previously, hierarchical clustering divided the 56 Ethiopian samples
# into groups of 47 and 9.
#
# Now that we have the real clade labels, we can directly cross-tabulate
# that two-cluster solution against A/B/D.
#
# If the previous interpretation is correct, one cluster should contain
# almost/all D samples and the other should contain the A/B samples.
#
# This clustering is ONLY used as a diagnostic check.
# It is NOT used for the biological analysis below.
#
# ======================================================================

eth_mt_hclust <- hclust(
  as.dist(
    eth_mt_pat
  ),
  method = "average"
)


eth_k2 <- cutree(
  eth_mt_hclust,
  k = 2
)


k2_check <- data.frame(
  sample_id = names(eth_k2),
  k2_cluster = as.integer(eth_k2),
  clade = eth_matched$clade[
    match(
      names(eth_k2),
      eth_matched$sample_id
    )
  ],
  stringsAsFactors = FALSE
)


k2_table <- table(
  k2_check$k2_cluster,
  k2_check$clade
)


cat("\n============================================================\n")
cat("PREVIOUS k = 2 CLUSTERING vs TRUE CLADE LABELS\n")
cat("============================================================\n")


print(k2_table)


# ======================================================================
# 44. CREATE ALL UNIQUE PAIRWISE COMPARISONS
# ======================================================================
#
# For n Ethiopian samples there are:
#
#                   n(n - 1) / 2
#
# unique sample pairs.
#
# Each pair receives:
#
#   - nuclear patristic distance
#   - mitochondrial patristic distance
#   - mitochondrial clade of sample 1
#   - mitochondrial clade of sample 2
#
# ======================================================================

pair_indices <- which(
  upper.tri(
    eth_nuc_pat
  ),
  arr.ind = TRUE
)


eth_pairs <- data.frame(
  
  sample1 = rownames(
    eth_nuc_pat
  )[
    pair_indices[, 1]
  ],
  
  sample2 = colnames(
    eth_nuc_pat
  )[
    pair_indices[, 2]
  ],
  
  nuclear_distance = eth_nuc_pat[
    pair_indices
  ],
  
  mtDNA_distance = eth_mt_pat[
    pair_indices
  ],
  
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# Add actual mitochondrial clade assignments
# ----------------------------------------------------------------------

matched_clade_lookup <- setNames(
  eth_matched$clade,
  eth_matched$sample_id
)


eth_pairs$clade1 <- unname(
  matched_clade_lookup[
    eth_pairs$sample1
  ]
)


eth_pairs$clade2 <- unname(
  matched_clade_lookup[
    eth_pairs$sample2
  ]
)


# ======================================================================
# 45. CLASSIFY PAIRS USING THE TRUE A/B/D CLADE LABELS
# ======================================================================
#
# This gives six possible pair categories:
#
#   Within A
#   Within B
#   Within D
#   A vs B
#   A vs D
#   B vs D
#
# Depending on whether one of the two B samples was excluded, there may
# be no "Within B" pair in the matched dataset.
#
# ======================================================================

same_clade <- (
  eth_pairs$clade1 ==
    eth_pairs$clade2
)


eth_pairs$pair_type_exact <- ifelse(
  
  same_clade,
  
  paste0(
    "Within ",
    eth_pairs$clade1
  ),
  
  paste0(
    pmin(
      eth_pairs$clade1,
      eth_pairs$clade2
    ),
    " vs ",
    pmax(
      eth_pairs$clade1,
      eth_pairs$clade2
    )
  )
)


eth_pairs$pair_type_exact <- factor(
  eth_pairs$pair_type_exact,
  levels = c(
    "Within A",
    "Within B",
    "Within D",
    "A vs B",
    "A vs D",
    "B vs D"
  )
)


cat("\n============================================================\n")
cat("EXACT A/B/D PAIR TYPES\n")
cat("============================================================\n")


print(
  table(
    eth_pairs$pair_type_exact,
    useNA = "ifany"
  )
)


# ======================================================================
# 46. SUMMARISE DISTANCES FOR EACH EXACT CLADE COMPARISON
# ======================================================================
#
# This is useful because it allows us to see whether Clade B behaves
# more like A or D rather than automatically combining A and B.
#
# ======================================================================

exact_distance_summary <- eth_pairs %>%
  filter(
    !is.na(pair_type_exact)
  ) %>%
  group_by(
    pair_type_exact
  ) %>%
  summarise(
    
    n_pairs = n(),
    
    mean_mtDNA = mean(
      mtDNA_distance,
      na.rm = TRUE
    ),
    
    median_mtDNA = median(
      mtDNA_distance,
      na.rm = TRUE
    ),
    
    mean_nuclear = mean(
      nuclear_distance,
      na.rm = TRUE
    ),
    
    median_nuclear = median(
      nuclear_distance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("DISTANCES BY TRUE mtDNA CLADE COMPARISON\n")
cat("============================================================\n")


print(
  exact_distance_summary
)


# ======================================================================
# 47. CREATE THE BROADER D versus A/B COMPARISON
# ======================================================================
#
# The previous two-cluster analysis suggested that the dominant
# mitochondrial division is between:
#
#       Clade D
#
# and:
#
#       non-D samples (Clades A and B)
#
# We therefore retain the full three-clade results above, but also test
# this broader division explicitly.
#
# Importantly, A and B are NOT being claimed to be the same clade.
# They are combined here only to ask whether the major D/non-D
# mitochondrial subdivision has an equivalent nuclear counterpart.
#
# ======================================================================

eth_pairs$broad1 <- ifelse(
  eth_pairs$clade1 == "D",
  "D",
  "A/B"
)


eth_pairs$broad2 <- ifelse(
  eth_pairs$clade2 == "D",
  "D",
  "A/B"
)


eth_pairs$pair_type_broad <- ifelse(
  
  eth_pairs$broad1 ==
    eth_pairs$broad2,
  
  ifelse(
    eth_pairs$broad1 == "D",
    "Within D",
    "Within A/B"
  ),
  
  "D vs A/B"
)


eth_pairs$pair_type_broad <- factor(
  eth_pairs$pair_type_broad,
  levels = c(
    "Within D",
    "Within A/B",
    "D vs A/B"
  )
)


# ======================================================================
# 48. SUMMARISE D versus A/B DISTANCES
# ======================================================================

broad_distance_summary <- eth_pairs %>%
  group_by(
    pair_type_broad
  ) %>%
  summarise(
    
    n_pairs = n(),
    
    mean_mtDNA = mean(
      mtDNA_distance,
      na.rm = TRUE
    ),
    
    median_mtDNA = median(
      mtDNA_distance,
      na.rm = TRUE
    ),
    
    mean_nuclear = mean(
      nuclear_distance,
      na.rm = TRUE
    ),
    
    median_nuclear = median(
      nuclear_distance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("D versus A/B DISTANCE SUMMARY\n")
cat("============================================================\n")


print(
  broad_distance_summary
)


# ======================================================================
# 49. CALCULATE BETWEEN / WITHIN DISTANCE RATIOS
# ======================================================================
#
# For each genomic compartment:
#
#   ratio =
#
#       mean distance between D and A/B
#       --------------------------------
#       mean distance within broad groups
#
#
# A large mitochondrial ratio combined with a nuclear ratio close to 1
# would demonstrate that D versus A/B is strongly differentiated in the
# mitochondrial genome but not equivalently differentiated in the
# nuclear genome.
#
# ======================================================================

broad_same <- (
  eth_pairs$broad1 ==
    eth_pairs$broad2
)


mean_mt_within <- mean(
  eth_pairs$mtDNA_distance[
    broad_same
  ],
  na.rm = TRUE
)


mean_mt_between <- mean(
  eth_pairs$mtDNA_distance[
    !broad_same
  ],
  na.rm = TRUE
)


mean_nuc_within <- mean(
  eth_pairs$nuclear_distance[
    broad_same
  ],
  na.rm = TRUE
)


mean_nuc_between <- mean(
  eth_pairs$nuclear_distance[
    !broad_same
  ],
  na.rm = TRUE
)


mt_between_within_ratio <- (
  mean_mt_between /
    mean_mt_within
)


nuc_between_within_ratio <- (
  mean_nuc_between /
    mean_nuc_within
)


cat("\n============================================================\n")
cat("D versus A/B: BETWEEN / WITHIN DISTANCE RATIOS\n")
cat("============================================================\n")


cat(
  "\nMean mtDNA distance within broad groups:",
  mean_mt_within,
  "\n"
)


cat(
  "Mean mtDNA distance between D and A/B:",
  mean_mt_between,
  "\n"
)


cat(
  "mtDNA between/within ratio:",
  mt_between_within_ratio,
  "\n"
)


cat(
  "\nMean nuclear distance within broad groups:",
  mean_nuc_within,
  "\n"
)


cat(
  "Mean nuclear distance between D and A/B:",
  mean_nuc_between,
  "\n"
)


cat(
  "Nuclear between/within ratio:",
  nuc_between_within_ratio,
  "\n"
)


# ======================================================================
# 50. TEST WHETHER TRUE mtDNA CLADES FORM TREE SPLITS
# ======================================================================
#
# Both trees are pruned to the matched Ethiopian samples.
#
# We then ask whether:
#
#   Clade A
#   Clade B
#   Clade D
#   non-D (A/B)
#
# form discrete splits in the nuclear and mitochondrial topologies.
#
# The key expectation for mitonuclear discordance is that Clade D or
# the broader D/non-D division forms a clear mitochondrial split without
# an equivalent nuclear split.
#
# ======================================================================

eth_nuc_tree <- keep.tip(
  nuc_unroot,
  eth_tips
)


eth_mt_tree <- keep.tip(
  mt_unroot,
  eth_tips
)


eth_nuc_tree <- unroot(
  eth_nuc_tree
)


eth_mt_tree <- unroot(
  eth_mt_tree
)


clade_A_tips <- eth_matched$sample_id[
  eth_matched$clade == "A"
]


clade_B_tips <- eth_matched$sample_id[
  eth_matched$clade == "B"
]


clade_D_tips <- eth_matched$sample_id[
  eth_matched$clade == "D"
]


non_D_tips <- eth_matched$sample_id[
  eth_matched$clade %in%
    c("A", "B")
]


split_groups <- list(
  "Clade A" = clade_A_tips,
  "Clade B" = clade_B_tips,
  "Clade D" = clade_D_tips,
  "Non-D (A/B)" = non_D_tips
)


split_results <- do.call(
  rbind,
  lapply(
    names(split_groups),
    function(group_name) {
      
      tips <- split_groups[[group_name]]
      
      data.frame(
        
        Group = group_name,
        
        n = length(tips),
        
        Nuclear_split = forms_split(
          eth_nuc_tree,
          tips
        ),
        
        mtDNA_split = forms_split(
          eth_mt_tree,
          tips
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


rownames(split_results) <- NULL


cat("\n============================================================\n")
cat("TRUE CLADE SPLIT TESTS\n")
cat("============================================================\n")


print(
  split_results
)


# ======================================================================
# 51. CREATE ONE ETHIOPIA FIGURE
# ======================================================================
#
# One figure is retained rather than creating several separate boxplots.
#
# Each point is one pair of Ethiopian samples.
#
# The shapes indicate whether that pair is:
#
#   within A
#   within B
#   within D
#   A vs B
#   A vs D
#   B vs D
#
# If D versus A/B drives discordance, A-D and B-D pairs should have
# large mitochondrial distances without an equivalent shift along the
# nuclear-distance axis.
#
# ======================================================================

p_eth <- ggplot(
  eth_pairs,
  aes(
    x = nuclear_distance,
    y = mtDNA_distance,
    shape = pair_type_exact
  )
) +
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  # geom_smooth(
  #   aes(
  #     group = 1
  #   ),
  #   method = "lm",
  #   se = FALSE,
  #   linewidth = 0.7
  # ) +
  labs(
    x = "Nuclear patristic distance",
    y = "Mitochondrial patristic distance",
    shape = "Mitochondrial clade comparison",
    title = "Mitonuclear correspondence within Ethiopia"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    legend.position = "right"
  )


print(
  p_eth
)


# ----------------------------------------------------------------------
# Overwrite the existing Ethiopia figure
# ----------------------------------------------------------------------

ggsave(
  "ethiopia_nuclear_vs_mtDNA_patristic_distances.png",
  plot = p_eth,
  width = 8,
  height = 5.5,
  dpi = 600
)


# ======================================================================
# 52. SAVE ONE SAMPLE-LEVEL FILE
# ======================================================================
#
# This file records all 61 original Ethiopian samples and makes it clear:
#
#   - original mitochondrial clade
#   - whether present in metadata
#   - whether present in each tree
#   - whether retained in the matched analysis
#
# This overwrites the previous file of the same name.
#
# ======================================================================

write.csv(
  eth_status,
  "ethiopia_samples_mitochondrial_clades.csv",
  row.names = FALSE
)


# ======================================================================
# 53. CREATE ONE COMPACT SUMMARY FILE
# ======================================================================
#
# Instead of writing multiple separate CSVs for:
#
#   clade counts
#   distance summaries
#   correlations
#   split tests
#   two-cluster checks
#
# they are combined into one long-format summary file.
#
# ======================================================================


# ----------------------------------------------------------------------
# Matched clade counts
# ----------------------------------------------------------------------

summary_counts <- data.frame(
  
  Section = "Matched clade counts",
  
  Group = names(
    matched_clade_counts
  ),
  
  Metric = "n_samples",
  
  Value = as.character(
    as.integer(
      matched_clade_counts
    )
  ),
  
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# Overall correlations
# ----------------------------------------------------------------------

summary_cor <- data.frame(
  
  Section = "Overall Ethiopia",
  
  Group = "All matched samples",
  
  Metric = c(
    "Pearson_r",
    "Spearman_rho"
  ),
  
  Value = as.character(
    c(
      eth_pearson,
      eth_spearman
    )
  ),
  
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# Exact A/B/D distance summaries
# ----------------------------------------------------------------------

exact_long <- exact_distance_summary %>%
  pivot_longer(
    cols = c(
      n_pairs,
      mean_mtDNA,
      median_mtDNA,
      mean_nuclear,
      median_nuclear
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  transmute(
    
    Section = "Exact A/B/D pair distances",
    
    Group = as.character(
      pair_type_exact
    ),
    
    Metric = Metric,
    
    Value = as.character(
      Value
    )
  )


# ----------------------------------------------------------------------
# Broad D versus A/B summaries
# ----------------------------------------------------------------------

broad_long <- broad_distance_summary %>%
  pivot_longer(
    cols = c(
      n_pairs,
      mean_mtDNA,
      median_mtDNA,
      mean_nuclear,
      median_nuclear
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  transmute(
    
    Section = "D versus A/B pair distances",
    
    Group = as.character(
      pair_type_broad
    ),
    
    Metric = Metric,
    
    Value = as.character(
      Value
    )
  )


# ----------------------------------------------------------------------
# Key ratios
# ----------------------------------------------------------------------

ratio_summary <- data.frame(
  
  Section = "D versus A/B explanation",
  
  Group = "Between versus within",
  
  Metric = c(
    "mtDNA_between_within_ratio",
    "nuclear_between_within_ratio"
  ),
  
  Value = as.character(
    c(
      mt_between_within_ratio,
      nuc_between_within_ratio
    )
  ),
  
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# Split results
# ----------------------------------------------------------------------

split_long <- split_results %>%
  pivot_longer(
    cols = c(
      Nuclear_split,
      mtDNA_split
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  transmute(
    
    Section = "Tree split tests",
    
    Group = paste0(
      Group,
      " (n=",
      n,
      ")"
    ),
    
    Metric = Metric,
    
    Value = as.character(
      Value
    )
  )


# ----------------------------------------------------------------------
# Previous k = 2 clustering diagnostic
# ----------------------------------------------------------------------

k2_long <- as.data.frame(
  k2_table,
  stringsAsFactors = FALSE
)


names(k2_long) <- c(
  "Cluster",
  "Clade",
  "Freq"
)


k2_long <- k2_long[
  k2_long$Freq > 0,
  ,
  drop = FALSE
]


k2_long <- k2_long %>%
  transmute(
    
    Section = "Previous k=2 diagnostic",
    
    Group = paste0(
      "Cluster ",
      Cluster,
      " / Clade ",
      Clade
    ),
    
    Metric = "n_samples",
    
    Value = as.character(
      Freq
    )
  )


# ----------------------------------------------------------------------
# Combine all summaries
# ----------------------------------------------------------------------

ethiopia_summary <- bind_rows(
  summary_counts,
  summary_cor,
  exact_long,
  broad_long,
  ratio_summary,
  split_long,
  k2_long
)


# ----------------------------------------------------------------------
# Overwrite the existing summary file
# ----------------------------------------------------------------------

write.csv(
  ethiopia_summary,
  "ethiopia_clade_explanation_summary.csv",
  row.names = FALSE
)


# ======================================================================
# 54. FINAL CONSOLE SUMMARY
# ======================================================================

cat("\n============================================================\n")
cat("ETHIOPIA CLADE ANALYSIS COMPLETE\n")
cat("============================================================\n")


cat("\nOriginal mitochondrial dataset:\n")

print(
  original_clade_counts
)


cat("\nMatched mitonuclear dataset:\n")

print(
  matched_clade_counts
)


cat("\nSamples lost during matching/QC:\n")

print(
  excluded_eth[
    ,
    c(
      "sample_id",
      "clade"
    )
  ]
)


cat("\nClade B status:\n")

print(
  original_B[
    ,
    c(
      "sample_id",
      "retained_matched"
    )
  ]
)


cat("\nPrevious k = 2 clustering versus true clades:\n")

print(
  k2_table
)


cat("\nOverall Ethiopia correlation:\n")

cat(
  "Pearson r =",
  eth_pearson,
  "\n"
)

cat(
  "Spearman rho =",
  eth_spearman,
  "\n"
)


cat("\nExact A/B/D distance summary:\n")

print(
  exact_distance_summary
)


cat("\nD versus A/B distance summary:\n")

print(
  broad_distance_summary
)


cat("\nBetween/within distance ratios:\n")

cat(
  "mtDNA =",
  mt_between_within_ratio,
  "\n"
)

cat(
  "Nuclear =",
  nuc_between_within_ratio,
  "\n"
)


cat("\nTree split tests:\n")

print(
  split_results
)


cat("\nFiles overwritten:\n")

cat(
  "1. ethiopia_samples_mitochondrial_clades.csv\n"
)

cat(
  "2. ethiopia_clade_explanation_summary.csv\n"
)

cat(
  "3. ethiopia_nuclear_vs_mtDNA_patristic_distances.png\n"
)


cat("\n============================================================\n")
