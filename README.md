# genomic-variant-toolkit
Python library and parser designed to process targeted Next-Generation Sequencing (NGS) outputs from Ion Torrent / Genexus systems, handling VCFv4.1 formats and run-level quality control metrics.

## Overview
Targeted gene panels in molecular pathology generate multidimensional outputs spanning single nucleotide variants (SNVs), copy number alterations (CNVs), and RNA fusions. This repository provides robust parsing routines to extract run QC metrics, parse embedded functional annotations (FUNC JSON), and structure somatic variants for clinical evaluation.

## Features
- **Header QC Parsing:** Automatically extracts MAPD, median coverage depth, fusion read counts, and assay QC calls.
- **Multi-variant Extraction:** Distinguishes and separates SNVs/indels, CNVs (with copy number and confidence intervals), and RNA transcript alterations.
- **Robust Field Extraction:** Tolerant parsing of complex annotations and missing data attributes.
- **Automated Test Suite:** Comprehensive `pytest` coverage validating parsing accuracy across edge cases.

## Installation
```bash
git clone [https://github.com/TU_USUARIO/clinical-ngs-vcf-parser.git](https://github.com/TU_USUARIO/clinical-ngs-vcf-parser.git)
cd clinical-ngs-vcf-parser
pip install -r requirements.txt
