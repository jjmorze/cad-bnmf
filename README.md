# cad-bnmf

## Decomposing Polygenic Risk of Coronary Artery Disease into Distinct Biological Processes

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![medRxiv](https://img.shields.io/badge/medRxiv-2026.10.08.26365088-b31b1b.svg)](https://doi.org/10.64898/2026.10.08.26365088)

---

## Manuscript Information

**Title:** Decomposing polygenic risk of coronary artery disease into distinct biological processes

**Authors:** Jakub Morze, Kirk Smith, Marcin Kierczak, Göran Bergström, Hui Li, Anders Gummesson, Tove Fall, Stefan Söderberg, Tino Ebbers, J. Gustav Smith, Rikard Landberg, Paul Franks, Tuuli Lappalainen, Miriam Udler, Clemens Wittenbecher

**Corresponding Authors:**
- Clemens Wittenbecher, PhD (clemens.wittenbecher@chalmers.se)

**Affiliations:** Department of Life Sciences, Chalmers University of Technology; Massachusetts General Hospital; Uppsala University; University of Gothenburg; Umeå University; Linköping University; Lund University; KTH Royal Institute of Technology.

**Preprint:** medRxiv 2026.10.08.26365088; doi: [10.64898/2026.10.08.26365088](https://doi.org/10.64898/2026.10.08.26365088)

---

## Overview

This repository contains the analysis code for generating all results, tables, and figures presented in the manuscript. Standard polygenic risk scores (PRS) capture overall genetic liability to coronary artery disease (CAD) but do not reveal which biological processes drive risk in an individual. We applied data-driven latent factor decomposition to 526 genome-wide significant CAD-associated variants, using their associations with 78 CAD-related traits, to derive process-specific polygenic risk scores (pPRSs). These were then evaluated in UK Biobank (n=458,905), the Swedish CArdioPulmonary bioImage Study (SCAPIS; n=24,904), and the Mass General Brigham Biobank (MGBB; n=53,848).

---

## Table of Contents

1. [Key Findings](#key-findings)
2. [Analytical Workflow](#analytical-workflow)
3. [Repository Structure](#repository-structure)
4. [Data Availability](#data-availability)
5. [Contact](#contact)
6. [License](#license)

---

## Key Findings

- **Genetic decomposition:** 526 CAD variants × 78 CAD-related traits decomposed into **8 genetic components**
- **Biological processes identified:** cholesterol-enriched lipoproteins, triglyceride-rich lipoproteins, lipoprotein(a), metabolically unhealthy normal weight, obesity, endothelial dysfunction, coagulation, and cardiac remodeling
- **Beyond conventional risk factors:** pPRSs captured biological profiles not reflected by traditional CAD risk factors
- **Replication across cohorts:** pPRSs were consistently associated with clinical and subclinical CAD in UK Biobank, SCAPIS, and MGBB
- **Distinct profiles:** each pPRS showed characteristic biomarker, proteomic, and comorbidity signatures

---

## Analytical Workflow

To be extended. Main steps:

1. **Genetic decomposition** — Bayesian non-negative matrix factorization (bNMF) of CAD variant–trait associations from public GWAS summary statistics
2. **pPRS construction** — process-specific scores derived from component weights
3. **Association analyses** — pPRSs vs. clinical CAD (UK Biobank, MGBB) and subclinical atherosclerosis (SCAPIS)
4. **Phenotypic characterization** — biomarker, proteomic, and comorbidity (PheWAS) profiles of each pPRS

---

## Repository Structure

To be extended.

```
cad-bnmf/
├── cad-bnmf.Rproj
├── R/            # analysis scripts / notebooks
├── output/       # summary-level tables and figures
├── LICENSE
└── README.md
```

---

## Data Availability

Individual-level data cannot be shared publicly but are available through the respective access procedures:

- **UK Biobank:** Access via [UK Biobank Access Management System](https://www.ukbiobank.ac.uk) (Application: 570811)
- **SCAPIS:** Access via the [SCAPIS data portal](https://www.scapis.org) (researchers based in, or affiliated with an institution in, Sweden; requires Swedish Ethical Review Authority approval)
- **MGBB:** Available on request from the [Mass General Brigham Biobank](https://research.massgeneralbrigham.org/en/researchers-and-collaborators/personalized-medicine/biobank), subject to IRB approval
- **GWAS summary statistics:** Publicly available; sources cited in Supplementary Tables 1 and 2
- **Summary-level results:** Provided in the Supplementary Information and this repository

---

## Funding

- Royal Swedish Academy of Sciences (ME2024-0009)
- SciLifeLab & Wallenberg Data Driven Life Science Program (KAW 2020.0239)
- Swedish Research Council (2022-01529)
- Swedish Heart and Lung Foundation (20241115)

---

## Contact

**Jakub Morze, MD, PhD**
Department of Life Sciences
Chalmers University of Technology
Kemivägen 4, 41296 Gothenburg, Sweden
Email: jakub.morze@chalmers.se
GitHub: [@jjmorze](https://github.com/jjmorze)

---

## Citation

```bibtex
@article{morze2026pprs,
  title={Decomposing polygenic risk of coronary artery disease into distinct
         biological processes},
  author={Morze, Jakub and Smith, Kirk and Kierczak, Marcin and others},
  journal={medRxiv},
  year={2026},
  doi={10.64898/2026.10.08.26365088}
}
```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

*Code repository: https://github.com/jjmorze/cad-bnmf*
*Last updated: 2026-10-26 | Version: 1.0.0*
