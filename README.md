## nanoASM
One has a conviction that proxy ASM might be more causal than mQTL. This repository was developed to perform such validation
![alt text](github_fig.png)

Pipeline for allele-specific methylation (ASM) analysis using nanopore long-read sequencing.

---

## Overview

This repository contains scripts for:

* **Basecalling (Dorado)**
* **Variant calling and phasing**
* **Allele-specific methylation (ASM) detection**

The workflow is designed for nanopore R10 data with modified base detection (5mCG / 5hmCG).

---

## Pipeline Structure

### 1. Basecalling

```bash
scripts/dorado_basecall_cg.sh
```

* Uses Dorado basecaller (R10.4.1 SUP model)
* Generates BAM files with CG-context modified bases
* Model: `5mCG / 5hmCG`

---

### 2. Call modified bases

Before running `bam_prep.sh`, call modified bases:

```bash
for bam in *.sorted.bam; do
  modkit call-mods -t 47 "$bam" "${bam%.sorted.bam}.called.bam"
done
```

* The resulting `*.called.bam` files are used as input for `bam_prep.sh`.

---

### 3. Variant Calling and Phasing

```bash
scripts/variant_call_phase.sh
```

* PEPPER-Margin-DeepVariant (optional / historical)
* SNP filtering via bcftools
* LongPhase phasing or population-based phasing
* Haplotagging BAM files

---

### 4. ASM Detection (Core)

```bash
scripts/asm_core.sh
```

Main steps:

* Label reads by allele using:

  ```bash
  utils/label_reads_by_variant.py
  ```
* Extract CpG methylation using `modkit`
* Build per-sample CpG `.bg` files
* Construct union CpG sites
* Filter sparse samples
* Run **metilene** for ASM detection
* Output:

  * DMR summary file
  * Optional matrix files for significant loci

---

### 5. Cohort Execution

Run ASM per cohort using Slurm:

```bash
scripts/run_asm_cancer.sh
scripts/run_asm_normal.sh
```

These scripts:

* iterate over chunked VCF inputs
* call `asm_core.sh`
* manage output directories and logs

---

## Key Dependencies

* Dorado
* samtools
* bcftools
* bedtools
* modkit
* metilene
* Python (pysam)
* R (for downstream analysis)

---

## Notes

* Scripts contain cluster-specific paths (`/share/lab_wangl/...`)
* These should be adapted for other environments
* Designed for high-depth nanopore WGS datasets

---

## Outputs

* `*.output.txt` → ASM/DMR results
* `*.matrix.txt` → per-locus methylation matrices
* log files for segmentation faults and filtering steps

---

## Author

Yijun Tian


