# ADARpipeline
### RNA Editing Variant Analysis in Myxoid Liposarcoma Cell Lines

This repository contains scripts and documentation for identifying and analyzing RNA editing events (specifically A>G and T>C changes, indicative of ADAR activity) across multiple myxoid liposarcoma (MLS) cell lines and treatment conditions.

---

## 📂 Project Overview
``` 
The workflow is designed to:
1. Identify RNA-specific variants by comparing RNA and DNA VCFs.
2. Focus on potential RNA editing sites (A>G and T>C).
3. Annotate variants using VEP https://www.ensembl.org.
4. Select variants in specific functional regions (e.g., 3'UTRs, missense).
5. Compare editing patterns between treated and control (DMSO) conditions.
6. Quantify allele frequencies (AF) and visualize overlap between conditions.
```  
---

## 🧬 Sample Information
**Cell Lines**: MLS-1765, MLS-402, MLS-Avory   
**Conditions**: DMSO (control), PANO (panobinostat), AZD and combination (HDACi and BRD4i)  
**Replicates**: 3–4 replicates per condition  
  
---

## 🛠️ Tools Used

- [BWA](https://bio-bwa.sourceforge.net/) mapping (for DNA), STAR (for RNA)
- [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) for variant calling
- [bcftools](https://samtools.github.io/bcftools/) — variant filtering, merging, and intersections   
- [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) — variant effect annotation  
- Python3 with modules `pandas`, `matplotlib`, and `venn`
  
---

## 👩‍🍳 Pre-processing pipeline  
`star.sh`  
`mutect2.sh`
```
- Reference: hg38.fa
- STAR mapping for RNAseq and BWA for DNA
- MarkDuplicates (Picard): Removes duplicate reads, outputs metrics, and creates a marked BAM file.
- SplitNCigarReads (GATK): Splits reads spanning introns to make them GATK-compatible.
- AddOrReplaceReadGroups (Picard): Adds metadata (e.g., RGID, RGSM) required for analysis.
- Base Recalibration (GATK BaseRecalibrator): Generates a recalibration table using known variant sites.
- Apply BQSR (GATK ApplyBQSR): Applies recalibration to update base quality scores.
- Variant Calling (GATK Mutect2): Detects somatic variants, generating a VCF file.
```

## 🧾 ADAR pipeline Summary

### 1. **VCF Processing (Bash Scripts)**

Located in `bashscripts/` 
Each script performs:  
- **DNA-RNA intersection** to retain RNA-only variants (`bcftools isec`)
- **Merging replicates** to create a multisample VCF (`bcftools merge`)
- **Filtering for depth and ADAR-type edits (A>G, T>C), DP>5 and in atleast 3 replicates** (`bcftools view`)
- **Optional filtering** to remove DMSO variants or common SNPs

📝 Example:
```bash
# Filter variants with depth > 5 in ALL 4 replicates and genotype present
bcftools view -i 'COUNT(GT!="." & FORMAT/DP>5)>=3' input.vcf.gz -Oz -o output.vcf.gz

# Filter for A->G and T->C variants (ADAR editing candidates)
bcftools view -i '(REF="A" & ALT="G") | (REF="T" & ALT="C")' input.gz -Oz -o output.vcf.gz
```

### 2. VEP Annotation & Region Extraction (CLI)
Annotated VCFs are parsed for:
	• 3_prime_UTR_variant
	• missense_variant
 
📝 Example commands:
```bash
grep "3_prime_UTR_variant\t" MLS-Avory_all3_r3.txt > MLS-Avory_all3_r3_3utr.txt
grep "missense_variant\t" MLS-1765_all3_r3.txt > MLS-1765_all3_r3_missense.txt
```

### 3. Allele Frequency and Visualization (Python Script)
`250703_treatmentVSdmso_ALL3.py` 

**Features:**  
``` 
• Parses .txt files for gene-position mapping  
• Compares treatment vs DMSO 3'UTR variants using Venn diagrams  
• Extracts AD and AF values from VCF  
• Computes max AF and AF ratio (treatment vs control)  
• Saves a CSV with the final results   
 ```

**Output:**  
``` 
• *_AF_Ratios_AD_3replicates_3utr.csv: Contains AF and AD data across conditions  
• Venn diagram showing overlap in RNA editing events  
```

### 📈 Output CSV Columns
``` 
• Position: Chromosomal position (e.g., chr12:112233-112233)
• GeneName: Mapped gene symbol
• AD_REF_*: Reference allele depth per sample
• AD_ALT_*: Alternate allele depth per sample
• AF_*: Allele frequency per sample
• Max_AF_T/NT: Max AF across treated/control replicates
• Max_AF_Ratio: Ratio of Max_AF_T / Max_AF_NT
```  

### 🧼 Requirements
```
• star version 2.5.2b
• picard-tools/2.9.0
• gatk version 4.6.1.0
• samtools version 2.0
• bcftools version 1.16
• htslib version 1.6
• Python3 packages:  pip install pandas matplotlib venn
	pandas==2.2.3
	matplotlib==3.10.1
	venn==0.1.3
```


### Improvments/questions:
``` 
* Is 1st replicate MLS-1765 weird?
* How to handle alternative allels, for example ALT: G,GT
* How to filter the DMSO best?
* Do you need to use QUAL filter as well?
* Use -w1 if piping output directly
* Remove -Oz -o from bcftools isec -C ... line already using -p which writes output to 0000.vcf inside a directory, so -Oz -o "$TREATMENT_SPECIFIC_VCF" can be ignored or misleading.
``` 
