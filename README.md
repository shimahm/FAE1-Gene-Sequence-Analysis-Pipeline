# FAE1-Gene-Sequence-Analysis-Pipeline

This pipeline identifies, extracts, aligns, and constructs a phylogenetic tree for FAE1 gene homologs in *Brassica rapa* and *Brassica oleracea* using BLAST, variant-aware consensus sequence generation, multiple sequence alignment, and phylogenetic tree construction.

---

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Directory Structure](#directory-structure)
- [Pipeline Steps](#pipeline-steps)
  - [Step 0: Identify FAE1 Homologs with BLAST](#step-0-identify-fae1-homologs-with-blast)
  - [Step 1: Extract Consensus Sequences](#step-1-extract-consensus-sequences)
  - [Step 2: Combine Consensus Sequences](#step-2-combine-consensus-sequences)
  - [Step 3: Multiple Sequence Alignment and Phylogenetic Tree](#step-3-multiple-sequence-alignment-and-phylogenetic-tree)
  - [Step 4: Visualize Phylogenetic Tree](#step-4-visualize-phylogenetic-tree)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Overview

This pipeline uses genomic data from *Brassica rapa* (Chiifu v4.1) and *Brassica oleracea* (JZS v2) to:

- Identify FAE1 gene homologs using BLAST against reference genomes.
- Create consensus sequences for five samples (A4, A8, A9 for *B. rapa*; C36, C47 for *B. oleracea*), incorporating variants from VCF files.
- Combine consensus sequences into a single FASTA file.
- Perform multiple sequence alignment and construct a phylogenetic tree.
- Visualize the phylogenetic relationships using R.

The pipeline generates 12 consensus sequences (6 for *B. rapa*, 6 for *B. oleracea*) and a phylogenetic tree showing relationships among FAE1 homologs.

---

## Prerequisites

**System:** Linux (e.g., Ubuntu)

**Software:**

- Conda (for environment management)
- BLAST+ (makeblastdb, blastn)
- Samtools (`samtools faidx`)
- Bcftools (`bcftools consensus`)
- Seqtk (`seqtk seq`)
- Tabix (`tabix`)
- MAFFT (`mafft`)
- FastTree (`FastTree`)
- R with `ape` package (for tree visualization)

**Conda Environment:**
```bash
conda create -n align_env -c bioconda blast samtools bcftools seqtk htslib mafft fasttree
conda activate align_env
```

**R Setup:**
```R
install.packages("ape")
```

---

## Directory Structure

```
parent_line/
├── fae1_coordinates.bed
├── fae1_coordinates_original.bed
├── A4.bam, A8.bam, A9.bam, C36.bam, C47.bam
├── A4.vcf.gz, A8.vcf.gz, A9.vcf.gz, C36.vcf.gz, C47.vcf.gz
├── A4.vcf.gz.tbi, A8.vcf.gz.tbi, A9.vcf.gz.tbi, C36.vcf.gz.tbi, C47.vcf.gz.tbi
├── step1_extract_consensus_with_strand.sh
├── step2_combine_sequences.sh
├── step3_align_and_tree.sh
├── fae1_analysis/
│   ├── A4_FAE1_brapa_1.fasta, A4_FAE1_brapa_2.fasta, ...
│   ├── all_fae1_sequences.fasta
│   ├── all_fae1_sequences_aligned.fasta
│   ├── fae1_phylogenetic_tree.nwk
├── /home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/
│   ├── ref_a/
│   │   ├── Brapa_chiifu_v41_genome20230413.fasta
│   │   ├── Brapa_chiifu_v41_genome20230413.fasta.fai
│   ├── ref_c/
│   │   ├── Brassica_oleracea_JZS_v2.fasta
│   │   ├── Brassica_oleracea_JZS_v2.fasta.fai
```

---

## Pipeline Steps

### Step 0: Identify FAE1 Homologs with BLAST

**Objective:** Identify FAE1 gene homologs in *B. rapa* and *B. oleracea* genomes using BLAST.

#### Commands:

```bash
# Create BLAST databases
makeblastdb -in /home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/ref_a/Brapa_chiifu_v41_genome20230413.fasta -dbtype nucl -out BrapaDB
makeblastdb -in /home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/ref_c/Brassica_oleracea_JZS_v2.fasta -dbtype nucl -out BoleraceaDB

# Run BLAST for B. rapa
blastn -query AT4G34520.gene.fa -db BrapaDB -out fae1_blast_brapa.out -evalue 1e-10 -outfmt 6 -num_alignments 10

# Run BLAST for B. oleracea
blastn -query AT4G34520.gene.fa -db BoleraceaDB -out fae1_blast_boleracea.out -evalue 1e-10 -outfmt 6 -num_alignments 10
```

**Output:**
- `fae1_blast_brapa.out`: BLAST hits for *B. rapa*.
- `fae1_blast_boleracea.out`: BLAST hits for *B. oleracea*.

#### Process BLAST Results

Convert BLAST outputs to BED format with strand information:

```bash
# B. rapa BED
awk 'BEGIN{OFS="\t"} {start=($9<$10)?$9:$10; end=($9<$10)?$10:$9; strand=($9<$10)?"+":"-"; print $2, start-1, end, "FAE1_"NR, ".", strand}' fae1_blast_brapa.out > FAE1_brapa_copies.bed

# B. oleracea BED
awk 'BEGIN{OFS="\t"} {start=($9<$10)?$9:$10; end=($9<$10)?$10:$9; strand=($9<$10)?"+":"-"; print $2, start-1, end, "FAE1_"NR, ".", strand}' fae1_blast_boleracea.out > FAE1_boleracea_copies.bed
```

##### Output BED Files Example:

`FAE1_brapa_copies.bed`
```
A08	20984059	20985581	FAE1_copy1	0	-
A09	6715834	6716953	FAE1_copy2	0	+
```
`FAE1_boleracea_copies.bed`
```
C03	68196636	68198158	FAE1_copy4	0	-
C09	8280133	8281083	FAE1_copy5	0	+
C01	2453445	2453823	FAE1_copy3	0	-
```

##### Combine and Rename BED Files:

Manually create `fae1_coordinates_original.bed` and `fae1_coordinates.bed`:

`fae1_coordinates_original.bed`
```
#chrom	start	end	name	score	strand
A08	20984059	20985581	FAE1_brapa_1	.	-
A09	6715834	6716953	FAE1_brapa_2	.	+
C03	68196636	68198158	FAE1_boleracea_1	.	-
C09	8280132	8281083	FAE1_boleracea_2	.	+
C01	2453445	2453823	FAE1_boleracea_3	.	-
```

`fae1_coordinates.bed`
```
A08	20984059	20985581	FAE1_brapa_1	rapa
A09	6715834	6716953	FAE1_brapa_2	rapa
C03	68196636	68198158	FAE1_boleracea_1	oleracea
C09	8280132	8281083	FAE1_boleracea_2	oleracea
C01	2453445	2453823	FAE1_boleracea_3	oleracea
```

---

### Step 1: Extract Consensus Sequences

**Objective:** Generate consensus sequences for FAE1 homologs from BAM and VCF files, incorporating strand information.

**Script:** `step1_extract_consensus_with_strand.sh`

<details>
<summary>Click to view script</summary>

```bash
#!/bin/bash
# Input files
B_RAPA_REF="/home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/ref_a/Brapa_chiifu_v41_genome20230413.fasta"
B_OLERACEA_REF="/home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/ref_c/Brassica_oleracea_JZS_v2.fasta"
COORDINATES="fae1_coordinates.bed"
ORIGINAL_COORDINATES="fae1_coordinates_original.bed"
OUT_DIR="fae1_analysis"
mkdir -p "$OUT_DIR"
rm -f "$OUT_DIR"/*.fasta "$OUT_DIR"/*.nwk
while IFS=$'\t' read -r chrom start end gene_id sample_type; do
    if [[ "$sample_type" == "rapa" ]]; then
        REF="$B_RAPA_REF"
        SAMPLES=("A4:bam_files_a/SRR14319585_sorted.bam" "A8:bam_files_a/SRR14319582_sorted.bam" "A9:bam_files_a/SRR14319581_sorted.bam")
    elif [[ "$sample_type" == "oleracea" ]]; then
        REF="$B_OLERACEA_REF"
        SAMPLES=("C36:bam_files_c/SRR14319592_sorted.bam" "C47:bam_files_c/SRR14319586_sorted.bam")
    else
        continue
    fi
    strand=$(awk -F'[ \t]+' -v c="$chrom" -v s="$start" -v e="$end" -v g="$gene_id" '$1==c && $2==s && $3==e && $4==g {print $6}' "$ORIGINAL_COORDINATES")
    if [[ -z "$strand" ]]; then
        case "$gene_id" in
            "FAE1_brapa_1"|"FAE1_boleracea_1"|"FAE1_boleracea_3") strand="-";;
            "FAE1_brapa_2"|"FAE1_boleracea_2") strand="+";;
            *) strand="+";;
        esac
    fi
    for sample_pair in "${SAMPLES[@]}"; do
        sample=$(echo "$sample_pair" | cut -d':' -f1)
        vcf_sample=$(echo "$sample_pair" | cut -d':' -f2)
        BAM="${sample}.bam"
        VCF="${sample}.vcf.gz"
        OUT_FASTA="${OUT_DIR}/${sample}_${gene_id}.fasta"
        TEMP_FASTA="${OUT_DIR}/${sample}_${gene_id}_temp.fasta"
        if [[ ! -f "$BAM" || ! -f "$VCF" || ! -f "$REF" || ! -f "${VCF}.tbi" ]]; then
            continue
        fi
        samtools faidx "$REF" "${chrom}:${start}-${end}" | \
        bcftools consensus -H 1 -s "$vcf_sample" -f - "$VCF" | \
        sed "s/>${chrom}.*/>${sample}_${gene_id}/" > "$TEMP_FASTA"
        if [[ "$strand" == "-" ]]; then
            seqtk seq -r "$TEMP_FASTA" > "$OUT_FASTA"
            rm "$TEMP_FASTA"
        else
            mv "$TEMP_FASTA" "$OUT_FASTA"
        fi
    done
done < "$COORDINATES"
```
</details>

#### Run:

```bash
conda activate align_env
tabix -p vcf A4.vcf.gz A8.vcf.gz A9.vcf.gz C36.vcf.gz C47.vcf.gz
chmod +x step1_extract_consensus_with_strand.sh
./step1_extract_consensus_with_strand.sh
```

**Output:**  
12 FASTA files in `fae1_analysis/` (e.g., `A4_FAE1_brapa_1.fasta`, ...).

---

### Step 2: Combine Consensus Sequences

**Objective:** Combine the 12 consensus sequences into a single FASTA file.

**Script:** `step2_combine_sequences.sh`

<details>
<summary>Click to view script</summary>

```bash
#!/bin/bash
OUT_DIR="fae1_analysis"
ALL_FASTA="${OUT_DIR}/all_fae1_sequences.fasta"
FASTA_COUNT=$(ls -1 "$OUT_DIR"/*.fasta | wc -l)
if [[ "$FASTA_COUNT" -ne 12 ]]; then
    echo "Error: Expected 12 FASTA files in $OUT_DIR, found $FASTA_COUNT."; exit 1
fi
cat "$OUT_DIR"/*.fasta > "$ALL_FASTA"
if [[ ! -s "$ALL_FASTA" ]]; then
    echo "Error: Combined FASTA file is empty."; exit 1
fi
```
</details>

#### Run:

```bash
conda activate align_env
chmod +x step2_combine_sequences.sh
./step2_combine_sequences.sh
```

**Output:**  
`fae1_analysis/all_fae1_sequences.fasta` (contains all 12 consensus sequences).

---

### Step 3: Multiple Sequence Alignment and Phylogenetic Tree

**Objective:** Align the consensus sequences and construct a phylogenetic tree.

**Script:** `step3_align_and_tree.sh`

<details>
<summary>Click to view script</summary>

```bash
#!/bin/bash
OUT_DIR="fae1_analysis"
ALL_FASTA="${OUT_DIR}/all_fae1_sequences.fasta"
ALIGNED_FASTA="${OUT_DIR}/all_fae1_sequences_aligned.fasta"
TREE_FILE="${OUT_DIR}/fae1_phylogenetic_tree.nwk"
if [[ ! -s "$ALL_FASTA" ]]; then
    echo "Error: Input file $ALL_FASTA is missing or empty."; exit 1
fi
mafft --auto "$ALL_FASTA" > "$ALIGNED_FASTA"
if [[ ! -s "$ALIGNED_FASTA" ]]; then
    echo "Error: Alignment file $ALIGNED_FASTA is empty or failed."; exit 1
fi
FastTree -nt -gtr "$ALIGNED_FASTA" > "$TREE_FILE"
if [[ ! -s "$TREE_FILE" ]]; then
    echo "Error: Tree file $TREE_FILE is empty or failed."; exit 1
fi
```
</details>

#### Run:

```bash
conda activate align_env
chmod +x step3_align_and_tree.sh
./step3_align_and_tree.sh
```

**Output:**
- `fae1_analysis/all_fae1_sequences_aligned.fasta`: Aligned sequences.
- `fae1_analysis/fae1_phylogenetic_tree.nwk`: Phylogenetic tree (Newick format).

---

### Step 4: Visualize Phylogenetic Tree

**Objective:** Visualize the phylogenetic tree in R.

**Script:** `visualize_tree.R`

```R
setwd("/home/mahmoudi/Bigdata/computing/Shima/parent_line")
library(ape)
tree <- read.tree("fae1_analysis/fae1_phylogenetic_tree.nwk")
plot(tree)
```

**Run:**
```bash
Rscript visualize_tree.R
```

**Output:**  
A graphical plot of the phylogenetic tree, showing relationships among FAE1 homologs.

**Alternative:**  
Use [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) for an interactive GUI:

```bash
sudo apt-get install figtree
figtree fae1_analysis/fae1_phylogenetic_tree.nwk
```

---

## Input Files

- **Reference Genomes:**
  - `/home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/ref_a/Brapa_chiifu_v41_genome20230413.fasta`
  - `/home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/ref_c/Brassica_oleracea_JZS_v2.fasta`
- **Query:** `AT4G34520.gene.fa` (Arabidopsis FAE1 gene sequence)
- **BAM Files:** `A4.bam`, `A8.bam`, `A9.bam`, `C36.bam`, `C47.bam`
- **VCF Files:** `A4.vcf.gz`, `A8.vcf.gz`, `A9.vcf.gz`, `C36.vcf.gz`, `C47.vcf.gz`
- **BED Files:**  
  - `fae1_coordinates.bed` (sample-type annotated)
  - `fae1_coordinates_original.bed` (strand-annotated)

---

## Output Files

- **Step 0:**  
  - `fae1_blast_brapa.out`, `fae1_blast_boleracea.out`
  - `FAE1_brapa_copies.bed`, `FAE1_boleracea_copies.bed`
- **Step 1:** 12 FASTA files in `fae1_analysis/` (e.g., `A4_FAE1_brapa_1.fasta`)
- **Step 2:** `fae1_analysis/all_fae1_sequences.fasta`
- **Step 3:**  
  - `fae1_analysis/all_fae1_sequences_aligned.fasta`
  - `fae1_analysis/fae1_phylogenetic_tree.nwk`

---

## Troubleshooting

### Step 0
- **Issue:** BLAST finds no hits.  
  **Solution:** Check `AT4G34520.gene.fa` format and sequence; adjust `-evalue`; check output with `head fae1_blast_brapa.out`.

### Step 1
- **Issue:** Empty FASTA files.  
  **Solution:** Verify existence of BAM, VCF, and reference files. Check VCF sample names (`zcat A4.vcf.gz | grep '^#CHROM'`).  
  Test with:  
  ```
  samtools faidx $B_RAPA_REF A08:20984059-20985581 | bcftools consensus -H 1 -s "bam_files_a/SRR14319585_sorted.bam" -f - A4.vcf.gz
  ```
- **Issue:** Strand errors.  
  **Solution:** Check `fae1_coordinates_original.bed` format (`cat -t fae1_coordinates_original.bed`).

### Step 2
- **Issue:** Incorrect file count.  
  **Solution:** Verify 12 FASTA files: `ls -l fae1_analysis/`. Check for stray files: `ls -l fae1_analysis/*_temp.fasta`.

### Step 3
- **Issue:** Alignment or tree fails.  
  **Solution:** Verify MAFFT and FastTree installation; check alignment: `head fae1_analysis/all_fae1_sequences_aligned.fasta`.  
  Rerun with bootstrap:  
  ```
  FastTree -nt -gtr -gamma -boot 1000
  ```

---

## License

This pipeline is provided under the MIT License. See [LICENSE](LICENSE) for details.

---

