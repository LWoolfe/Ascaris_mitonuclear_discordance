#!/bin/bash

# ======================================================================
# Prepare matched mitochondrial and nuclear sample sets for tanglegram
# ======================================================================
#
# Purpose:
#   Identify individuals represented in both the mitochondrial dataset
#   from Chapter 4 and the QC-filtered nuclear SNP dataset from Chapter 5.
#
#   Known differences in sample naming between datasets are first
#   harmonised. Samples are then intersected to generate a single list
#   of individuals to be used for both mitochondrial and nuclear NJ trees.
#
# Input:
#   - Concatenated alignment of the 12 mitochondrial protein-coding genes
#   - Final QC-filtered 146-sample autosomal nuclear SNP VCF
#
# Output:
#   - Harmonised mitochondrial FASTA
#   - Lists of mitochondrial and nuclear samples
#   - List of samples shared between datasets
#   - 141-sample mitochondrial FASTA
#
# ======================================================================

set -euo pipefail

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

WORKDIR=/parallel_scratch/lw01165/00_mitogenome/tanglegram

MTDNA_FASTA=${WORKDIR}/ascaris_12genes_concatenated_nucleotide.fasta

NUCLEAR_VCF=/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz

cd "${WORKDIR}"

echo "============================================================"
echo "Preparing matched samples for mitonuclear comparison"
echo "Started: $(date)"
echo "============================================================"


# ======================================================================
# STEP 1: Extract mitochondrial sample IDs
# ======================================================================

# FASTA headers are assumed to begin with the sample ID.
# Any text after the first whitespace is ignored.

grep "^>" "${MTDNA_FASTA}" \
| sed 's/^>//' \
| cut -d' ' -f1 \
> mtDNA_samples.original.txt

sort mtDNA_samples.original.txt \
> mtDNA_samples.original.sorted.txt

echo
echo "Number of samples in original mitochondrial dataset:"
wc -l < mtDNA_samples.original.txt


# ======================================================================
# STEP 2: Extract nuclear sample IDs
# ======================================================================

bcftools query -l "${NUCLEAR_VCF}" \
> nuclear_samples.txt

sort nuclear_samples.txt \
> nuclear_samples.sorted.txt

echo
echo "Number of samples in QC-filtered nuclear dataset:"
wc -l < nuclear_samples.txt


# ======================================================================
# STEP 3: Identify initial differences between datasets
# ======================================================================

echo
echo "Mitochondrial samples without an exact nuclear ID:"
comm -23 \
mtDNA_samples.original.sorted.txt \
nuclear_samples.sorted.txt \
> mtDNA_only.before_rename.txt

cat mtDNA_only.before_rename.txt

echo
echo "Nuclear samples without an exact mitochondrial ID:"
comm -13 \
mtDNA_samples.original.sorted.txt \
nuclear_samples.sorted.txt \
> nuclear_only.before_rename.txt

cat nuclear_only.before_rename.txt


# ======================================================================
# STEP 4: Harmonise known sample-name differences
# ======================================================================
#
# Six individuals had "_merge" appended to their identifiers in the
# nuclear dataset but not in the mitochondrial dataset.
#
# The mitochondrial identifiers are changed to match the nuclear VCF.
# The original FASTA is retained unchanged.
#
#   mtDNA ID          Nuclear ID
#   ------------------------------------------------
#   AL17SLB           AL17SLB_merge
#   AL2SLB            AL2SLB_merge
#   BunawanS10        BunawanS10_merge
#   BunawanS9         BunawanS9_merge
#   DE_P17a           DE_P17a_merge
#   TrentoS7          TrentoS7_merge
#
# ======================================================================

MTDNA_RENAMED=${WORKDIR}/ascaris_12genes_concatenated_nucleotide.tanglegram_names.fasta

cp "${MTDNA_FASTA}" "${MTDNA_RENAMED}"

sed -i \
-e 's/^>AL17SLB\( .*\)\?$/>AL17SLB_merge\1/' \
-e 's/^>AL2SLB\( .*\)\?$/>AL2SLB_merge\1/' \
-e 's/^>BunawanS10\( .*\)\?$/>BunawanS10_merge\1/' \
-e 's/^>BunawanS9\( .*\)\?$/>BunawanS9_merge\1/' \
-e 's/^>DE_P17a\( .*\)\?$/>DE_P17a_merge\1/' \
-e 's/^>TrentoS7\( .*\)\?$/>TrentoS7_merge\1/' \
"${MTDNA_RENAMED}"


# ======================================================================
# STEP 5: Regenerate mitochondrial sample list after renaming
# ======================================================================

grep "^>" "${MTDNA_RENAMED}" \
| sed 's/^>//' \
| cut -d' ' -f1 \
> mtDNA_samples.txt

sort mtDNA_samples.txt \
> mtDNA_samples.sorted.txt


# ======================================================================
# STEP 6: Identify samples present in both datasets
# ======================================================================

comm -12 \
mtDNA_samples.sorted.txt \
nuclear_samples.sorted.txt \
> tanglegram_shared_samples.txt

N_SHARED=$(wc -l < tanglegram_shared_samples.txt)

echo
echo "Number of samples shared between datasets:"
echo "${N_SHARED}"

if [[ "${N_SHARED}" -ne 141 ]]; then
echo "ERROR: Expected 141 shared samples but found ${N_SHARED}."
exit 1
fi


# ======================================================================
# STEP 7: Record samples found in only one dataset
# ======================================================================

comm -23 \
mtDNA_samples.sorted.txt \
nuclear_samples.sorted.txt \
> mtDNA_only.final.txt

comm -13 \
mtDNA_samples.sorted.txt \
nuclear_samples.sorted.txt \
> nuclear_only.final.txt

echo
echo "Mitochondrial-only sequences:"
cat mtDNA_only.final.txt

echo
echo "Nuclear-only samples:"
cat nuclear_only.final.txt


# ======================================================================
# STEP 8: Filter mitochondrial FASTA to the 141 shared individuals
# ======================================================================

MTDNA_N141=${WORKDIR}/ascaris_12genes_concatenated_nucleotide.tanglegram.n141.fasta

seqkit grep \
-f tanglegram_shared_samples.txt \
"${MTDNA_RENAMED}" \
> "${MTDNA_N141}"


# ======================================================================
# STEP 9: Verify filtered mitochondrial dataset
# ======================================================================

N_MTDNA=$(grep -c "^>" "${MTDNA_N141}")

echo
echo "Sequences in filtered mitochondrial FASTA:"
echo "${N_MTDNA}"

if [[ "${N_MTDNA}" -ne 141 ]]; then
echo "ERROR: Expected 141 mitochondrial sequences."
exit 1
fi

grep "^>" "${MTDNA_N141}" \
| sed 's/^>//' \
| cut -d' ' -f1 \
| sort \
> mtDNA_n141.sorted.txt

if diff -q \
mtDNA_n141.sorted.txt \
tanglegram_shared_samples.txt > /dev/null
then
echo "PASS: Mitochondrial FASTA contains exactly the 141 shared samples."
else
  echo "ERROR: Mitochondrial sample identities do not match shared sample list."
exit 1
fi

echo
echo "============================================================"
echo "Sample preparation completed successfully"
echo "Completed: $(date)"
echo "============================================================"

#!/bin/bash
#SBATCH --job-name=tangle_subset_vcf
#SBATCH --partition=general
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=tangle_subset_vcf.%j.out
#SBATCH --error=tangle_subset_vcf.%j.err


# ======================================================================
# Subset nuclear VCF to samples shared with mitochondrial dataset
# ======================================================================
#
# Purpose:
#   Extract the 141 individuals represented in both the mitochondrial
#   and nuclear datasets.
#
# Input:
#   Final QC-filtered Chapter 5 autosomal SNP VCF (n = 146)
#
# Output:
#   Nuclear SNP VCF containing exactly the 141 individuals shared
#   with the mitochondrial dataset.
#
# ======================================================================

set -euo pipefail

echo "============================================================"
echo "Creating nuclear VCF for mitonuclear comparison"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"


# ----------------------------------------------------------------------
# Software
# ----------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

echo
bcftools --version | head -n 1


# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

WORKDIR=/parallel_scratch/lw01165/00_mitogenome/tanglegram

INPUT_VCF=/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz

SAMPLE_LIST=${WORKDIR}/tanglegram_shared_samples.txt

OUTPUT_VCF=${WORKDIR}/Ascaris.autosomes.SNPs.tanglegram.n141.vcf.gz


# ======================================================================
# STEP 1: Check input files
# ======================================================================

if [[ ! -f "${INPUT_VCF}" ]]; then
echo "ERROR: Input VCF not found:"
echo "${INPUT_VCF}"
exit 1
fi

if [[ ! -f "${SAMPLE_LIST}" ]]; then
echo "ERROR: Shared sample list not found:"
echo "${SAMPLE_LIST}"
exit 1
fi

echo
echo "Input nuclear VCF:"
echo "${INPUT_VCF}"

echo
echo "Number of requested samples:"
wc -l < "${SAMPLE_LIST}"


# ======================================================================
# STEP 2: Subset VCF
# ======================================================================

echo
echo "Subsetting nuclear VCF..."

bcftools view \
--samples-file "${SAMPLE_LIST}" \
--threads 2 \
--output-type z \
--output "${OUTPUT_VCF}" \
"${INPUT_VCF}"


# ======================================================================
# STEP 3: Index subsetted VCF
# ======================================================================

echo
echo "Indexing output VCF..."

bcftools index \
--tbi \
--threads 2 \
"${OUTPUT_VCF}"


# ======================================================================
# STEP 4: Verify sample number
# ======================================================================

N_INPUT=$(bcftools query -l "${INPUT_VCF}" | wc -l)
N_REQUESTED=$(wc -l < "${SAMPLE_LIST}")
N_OUTPUT=$(bcftools query -l "${OUTPUT_VCF}" | wc -l)

echo
echo "Samples in original nuclear VCF: ${N_INPUT}"
echo "Samples requested:               ${N_REQUESTED}"
echo "Samples in tanglegram VCF:       ${N_OUTPUT}"

if [[ "${N_OUTPUT}" -ne "${N_REQUESTED}" ]]; then
echo "ERROR: Output sample number does not match requested sample number."
exit 1
fi


# ======================================================================
# STEP 5: Verify exact sample identities
# ======================================================================

bcftools query -l "${OUTPUT_VCF}" \
| sort \
> "${WORKDIR}/nuclear_tanglegram_samples.sorted.txt"

sort "${SAMPLE_LIST}" \
> "${WORKDIR}/tanglegram_shared_samples.sorted.txt"

if diff -q \
"${WORKDIR}/tanglegram_shared_samples.sorted.txt" \
"${WORKDIR}/nuclear_tanglegram_samples.sorted.txt" > /dev/null
then
echo "PASS: Nuclear VCF contains exactly the requested samples."
else
  echo "ERROR: Nuclear sample identities do not match."
diff \
"${WORKDIR}/tanglegram_shared_samples.sorted.txt" \
"${WORKDIR}/nuclear_tanglegram_samples.sorted.txt"
exit 1
fi


echo
echo "Output VCF:"
echo "${OUTPUT_VCF}"

echo
echo "Index:"
echo "${OUTPUT_VCF}.tbi"

echo
echo "============================================================"
echo "Completed successfully: $(date)"
echo "============================================================"

#!/bin/bash
#SBATCH --job-name=tangle_nuclear_prep
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=tangle_nuclear_prep.%j.out
#SBATCH --error=tangle_nuclear_prep.%j.err


# ======================================================================
# Prepare nuclear SNP dataset for mitonuclear NJ comparison
# ======================================================================
#
# Purpose:
#   Re-filter the Chapter 5 nuclear SNP dataset after reducing the
#   sample set from 146 to the 141 individuals represented in the
#   mitochondrial dataset.
#
# Steps:
#   1. Remove variants that became monomorphic following sample removal.
#   2. Convert the remaining variants to PLINK format.
#   3. Perform LD pruning using the same parameters as Chapter 5:
#          --indep-pairwise 50 10 0.15
#   4. Select an exact random subset of 100,000 LD-pruned SNPs.
#
# A fixed random seed is used to make SNP selection reproducible.
#
# ======================================================================

set -euo pipefail

echo "============================================================"
echo "Preparing nuclear SNP dataset for Ascaris tanglegram"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"


# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

WORKDIR=/parallel_scratch/lw01165/00_mitogenome/tanglegram

INPUT_VCF=${WORKDIR}/Ascaris.autosomes.SNPs.tanglegram.n141.vcf.gz

POLY_VCF=${WORKDIR}/Ascaris.autosomes.SNPs.tanglegram.n141.polymorphic.vcf.gz

# File containing the autosomal chromosome identifiers used in
# the Chapter 5 analysis.
CHROM_FILE=/parallel_scratch/lw01165/chrom_order.txt

PLINK_BASE=${WORKDIR}/Ascaris.tanglegram.n141.polymorphic

PRUNE_OUT=${WORKDIR}/Ascaris.tanglegram.n141.LD015

PRUNED_BASE=${WORKDIR}/Ascaris.tanglegram.n141.LD015.pruned

FINAL_BASE=${WORKDIR}/Ascaris.tanglegram.n141.LD015.100k


# ======================================================================
# STEP 1: Check input
# ======================================================================

if [[ ! -f "${INPUT_VCF}" ]]; then
echo "ERROR: Input VCF not found:"
echo "${INPUT_VCF}"
exit 1
fi

echo
echo "Input VCF:"
echo "${INPUT_VCF}"


# ======================================================================
# STEP 2: Remove SNPs that became monomorphic
# ======================================================================
#
# Removing individuals can cause SNPs that were polymorphic in the
# original n=146 dataset to become fixed within the n=141 dataset.
#
# Requiring a minor allele count >=1 retains only sites that remain
# polymorphic in the reduced dataset.
#
# ======================================================================

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

echo
echo "bcftools version:"
bcftools --version | head -n 1

N_BEFORE=$(bcftools view -H "${INPUT_VCF}" | wc -l)

echo
echo "Variants before monomorphic-site removal: ${N_BEFORE}"

bcftools view \
-c 1:minor \
--threads 8 \
-Oz \
-o "${POLY_VCF}" \
"${INPUT_VCF}"

bcftools index \
--tbi \
--threads 8 \
"${POLY_VCF}"

N_POLY=$(bcftools view -H "${POLY_VCF}" | wc -l)

echo "Variants remaining polymorphic:  ${N_POLY}"
echo "Variants removed as monomorphic: $((N_BEFORE - N_POLY))"


# ======================================================================
# STEP 3: Convert polymorphic VCF to PLINK format
# ======================================================================

conda activate NJ_env

echo
echo "PLINK2 version:"
plink2 --version

CHR_LIST=$(paste -sd, "${CHROM_FILE}")

echo
echo "Converting VCF to PLINK binary format..."

plink2 \
--vcf "${POLY_VCF}" \
--chr "${CHR_LIST}" \
--make-bed \
--allow-extra-chr \
--set-all-var-ids @:# \
  --out "${PLINK_BASE}"

echo
echo "Samples in PLINK dataset:"
wc -l < "${PLINK_BASE}.fam"

echo
echo "Variants in PLINK dataset:"
wc -l < "${PLINK_BASE}.bim"


# ======================================================================
# STEP 4: LD pruning
# ======================================================================
#
# Parameters reproduce the approach used for the Chapter 5
# population-structure analyses:
#
#   Window size = 50 SNPs
#   Step size   = 10 SNPs
#   r2 threshold = 0.15
#
# ======================================================================

echo
echo "Performing LD pruning using: 50 10 0.15"

plink2 \
--bfile "${PLINK_BASE}" \
--allow-extra-chr \
--indep-pairwise 50 10 0.15 \
--out "${PRUNE_OUT}"

N_PRUNE_IN=$(wc -l < "${PRUNE_OUT}.prune.in")

echo
echo "Variants retained after LD pruning:"
echo "${N_PRUNE_IN}"


# ----------------------------------------------------------------------
# Extract LD-pruned SNPs
# ----------------------------------------------------------------------

plink2 \
--bfile "${PLINK_BASE}" \
--extract "${PRUNE_OUT}.prune.in" \
--allow-extra-chr \
--make-bed \
--out "${PRUNED_BASE}"

N_PRUNED=$(wc -l < "${PRUNED_BASE}.bim")

echo
echo "LD-pruned variants in extracted dataset:"
echo "${N_PRUNED}"


# ======================================================================
# STEP 5: Select exactly 100,000 SNPs
# ======================================================================
#
# An exact random subset of 100,000 LD-pruned SNPs is selected.
#
# A fixed seed is provided so that the same SNP subset can be
# reproduced in future runs.
#
# ======================================================================

if [[ "${N_PRUNED}" -lt 100000 ]]; then
echo "ERROR: Fewer than 100,000 SNPs remain after LD pruning."
echo "Cannot generate requested 100,000-SNP dataset."
exit 1
fi

echo
echo "Selecting exact random subset of 100,000 SNPs..."

plink2 \
--bfile "${PRUNED_BASE}" \
--allow-extra-chr \
--thin-count 100000 \
--seed 20260825 \
--make-bed \
--out "${FINAL_BASE}"


# ======================================================================
# STEP 6: Verify final dataset
# ======================================================================

N_FINAL=$(wc -l < "${FINAL_BASE}.bim")
N_SAMPLES=$(wc -l < "${FINAL_BASE}.fam")

echo
echo "============================================================"
echo "FINAL NUCLEAR DATASET"
echo "============================================================"
echo "Samples: ${N_SAMPLES}"
echo "SNPs:    ${N_FINAL}"
echo
echo "Dataset prefix:"
echo "${FINAL_BASE}"

if [[ "${N_SAMPLES}" -ne 141 ]]; then
echo "ERROR: Final dataset does not contain 141 samples."
exit 1
fi

if [[ "${N_FINAL}" -ne 100000 ]]; then
echo "ERROR: Final dataset does not contain exactly 100,000 SNPs."
exit 1
fi

echo
echo "PASS: Final dataset contains 141 samples and 100,000 SNPs."

echo
echo "Completed: $(date)"
echo "============================================================"

#!/bin/bash
#SBATCH --job-name=tangle_nuclear_NJ
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tangle_nuclear_NJ.%j.out
#SBATCH --error=tangle_nuclear_NJ.%j.err


# ======================================================================
# Construct nuclear neighbour-joining tree
# ======================================================================
#
# Purpose:
#   Calculate genome-wide pairwise genetic distances among the 141
#   shared individuals using the exact random subset of 100,000
#   LD-pruned nuclear SNPs.
#
#   Genetic distance is calculated as:
#
#       1 - IBS
#
#   where IBS is identity-by-state between each pair of individuals.
#
#   A neighbour-joining tree is then constructed from the resulting
#   pairwise distance matrix using the R package ape.
#
# Output:
#   - Pairwise 1-IBS distance matrix
#   - Corresponding sample-ID file
#   - Nuclear NJ tree in Newick format
#
# ======================================================================

set -euo pipefail

echo "============================================================"
echo "Nuclear NJ tree for Ascaris mitonuclear comparison"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"


# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

WORKDIR=/parallel_scratch/lw01165/00_mitogenome/tanglegram

INPUT=${WORKDIR}/Ascaris.tanglegram.n141.LD015.100k

DIST=${WORKDIR}/Ascaris.tanglegram.n141.100k.1minusIBS

OUTPUT_TREE=${WORKDIR}/Ascaris.nuclear.NJ.n141.100k.nwk

cd "${WORKDIR}"


# ======================================================================
# STEP 1: Load software
# ======================================================================

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate NJ_env

echo
echo "PLINK version:"
plink --version

echo
echo "R version:"
R --version | head -n 1


# ======================================================================
# STEP 2: Check final PLINK input dataset
# ======================================================================

if [[ ! -f "${INPUT}.bed" || \
      ! -f "${INPUT}.bim" || \
      ! -f "${INPUT}.fam" ]]; then

echo "ERROR: PLINK input dataset not found."
exit 1
fi

N_SAMPLES=$(wc -l < "${INPUT}.fam")
N_SNPS=$(wc -l < "${INPUT}.bim")

echo
echo "Input dataset:"
echo "Samples: ${N_SAMPLES}"
echo "SNPs:    ${N_SNPS}"

if [[ "${N_SAMPLES}" -ne 141 ]]; then
echo "ERROR: Expected 141 samples."
exit 1
fi

if [[ "${N_SNPS}" -ne 100000 ]]; then
echo "ERROR: Expected 100,000 SNPs."
exit 1
fi


# ======================================================================
# STEP 3: Calculate pairwise nuclear genetic distances
# ======================================================================
#
# PLINK calculates pairwise identity-by-state and reports
# genetic distance as 1 - IBS.
#
# 'square' requests a complete square distance matrix.
#
# ======================================================================

echo
echo "Calculating pairwise 1-IBS genetic distances..."

plink \
--bfile "${INPUT}" \
--distance square 1-ibs \
--allow-extra-chr \
--out "${DIST}"


# ----------------------------------------------------------------------
# Expected PLINK output:
#
#   Ascaris.tanglegram.n141.100k.1minusIBS.mdist
#       141 x 141 distance matrix
#
#   Ascaris.tanglegram.n141.100k.1minusIBS.mdist.id
#       Sample identifiers in matrix order
#
# ----------------------------------------------------------------------

if [[ ! -f "${DIST}.mdist" ]]; then
echo "ERROR: PLINK distance matrix was not created."
exit 1
fi

if [[ ! -f "${DIST}.mdist.id" ]]; then
echo "ERROR: PLINK distance ID file was not created."
exit 1
fi

echo
echo "Pairwise distance matrix created successfully."


# ======================================================================
# STEP 4: Construct neighbour-joining tree in R
# ======================================================================
#
# The square PLINK distance matrix is imported into R and converted
# to a 'dist' object.
#
# The neighbour-joining tree is constructed using ape::nj().
#
# The resulting tree is exported in Newick format for downstream
# visualisation and comparison with the mitochondrial NJ tree.
#
# ======================================================================

echo
echo "Constructing neighbour-joining tree in R..."

Rscript - <<'EOF'

# ----------------------------------------------------------------------
# Nuclear NJ tree from PLINK 1-IBS distance matrix
# ----------------------------------------------------------------------

library(ape)

workdir <- "/parallel_scratch/lw01165/00_mitogenome/tanglegram"

dist_file <- file.path(
  workdir,
  "Ascaris.tanglegram.n141.100k.1minusIBS.mdist"
)

id_file <- file.path(
  workdir,
  "Ascaris.tanglegram.n141.100k.1minusIBS.mdist.id"
)

out_tree <- file.path(
  workdir,
  "Ascaris.nuclear.NJ.n141.100k.nwk"
)


# ----------------------------------------------------------------------
# Read sample identifiers
# ----------------------------------------------------------------------
#
# PLINK .mdist.id files contain:
#
#   FID IID
#
# The individual ID (IID; column 2) is used as the tree-tip label.
#
# ----------------------------------------------------------------------

ids <- read.table(
  id_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

sample_ids <- ids[, 2]


# ----------------------------------------------------------------------
# Read square genetic-distance matrix
# ----------------------------------------------------------------------

D <- as.matrix(
  read.table(
    dist_file,
    header = FALSE,
    check.names = FALSE
  )
)


# ----------------------------------------------------------------------
# Verify matrix dimensions
# ----------------------------------------------------------------------

if (nrow(D) != length(sample_ids)) {
  stop(
    "Distance matrix dimensions do not match number of sample IDs."
  )
}

if (ncol(D) != length(sample_ids)) {
  stop(
    "Distance matrix dimensions do not match number of sample IDs."
  )
}


# ----------------------------------------------------------------------
# Assign sample names to rows and columns
# ----------------------------------------------------------------------

rownames(D) <- sample_ids
colnames(D) <- sample_ids


# ----------------------------------------------------------------------
# Confirm matrix is symmetric
# ----------------------------------------------------------------------

if (!isTRUE(all.equal(D, t(D), tolerance = 1e-10))) {
  stop("Distance matrix is not symmetric.")
}


# ----------------------------------------------------------------------
# Confirm diagonal distances are effectively zero
# ----------------------------------------------------------------------

if (any(abs(diag(D)) > 1e-10)) {
  stop("Non-zero values detected on distance-matrix diagonal.")
}


# ----------------------------------------------------------------------
# Convert matrix to R distance object
# ----------------------------------------------------------------------

Ddist <- as.dist(D)


# ----------------------------------------------------------------------
# Construct neighbour-joining tree
# ----------------------------------------------------------------------

nuclear_nj <- nj(Ddist)


# ----------------------------------------------------------------------
# Verify number of tree tips
# ----------------------------------------------------------------------

if (length(nuclear_nj$tip.label) != 141) {
  stop("NJ tree does not contain 141 tips.")
}


# ----------------------------------------------------------------------
# Export tree in Newick format
# ----------------------------------------------------------------------

write.tree(
  nuclear_nj,
  file = out_tree
)


# ----------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------

cat("\nNuclear NJ tree created successfully\n")
cat("Samples:", length(nuclear_nj$tip.label), "\n")
cat("Output:", out_tree, "\n")

EOF


# ======================================================================
# STEP 5: Check Newick output
# ======================================================================

if [[ ! -s "${OUTPUT_TREE}" ]]; then
echo "ERROR: Newick tree was not generated."
exit 1
fi

echo
echo "Newick tree:"
echo "${OUTPUT_TREE}"

echo
echo "============================================================"
echo "Nuclear NJ analysis completed successfully"
echo "Completed: $(date)"
echo "============================================================"

r