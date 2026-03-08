# SIV_HIV_evolution

Data and analysis files for the manuscript:

**“A unified evolutionary framework to reconstruct the viral histories of human, simian, and prosimian immunodeficiency viruses across timescales”**

This repository contains the sequence alignments, BEAST XML files, phylogenetic tree outputs, and R code used to reconstruct the evolutionary history of HIV, SIV, and pSIV lineages using the **Prisoner of War (PoW)** framework.

## Repository contents

### Sequence alignments

- **`Fullgenome.fasta`**  
  Full-genome alignment of the complete lentivirus dataset used in the study, including HIV, SIV, and pSIVgml.

- **`hiv_clock_dataset.fasta`**  
  Time-stamped HIV-1 group M alignment used to estimate short-term substitution rates for each non-recombinant region (NRR).

---

### BEAST XML files for temporal signal and clock analyses

For each non-recombinant region (NRR), the repository contains two BEAST XML files:

- **`NRR_*_dated_iso_v1.xml`**  
  Isochronous analysis for the relevant NRR.

- **`NRR_*_dated_hetero_v1.xml`**  
  Heterochronous (tip-dated) analysis for the relevant NRR.

These XMLs were used for **Bayesian Evaluation of Temporal Signal (BETS)** via stepping-stone sampling to compare isochronous and heterochronous models and assess temporal signal in each NRR.

`*` denotes the NRR number.

---

### PoW-transformed trees

- **`MCC_filtered_PoWtransformed_NRR_*.tree`**  
  Maximum clade credibility (MCC) trees obtained after applying the PoW distance-to-time transformation to the NRRs with sufficient temporal signal.

These files correspond to the NRRs retained for downstream time-scaled phylogenetic inference.

---

### Analysis code

- **`PoW_transformation.R`**  
  R script used to:
  - read posterior tree and rate outputs,
  - apply the PoW transformation,
  - generate PoW-transformed time-scaled trees,
  - summarise node ages and divergence times used in the manuscript.

---

## Analysis overview

The study proceeds in four main steps:

1. **Dataset assembly and alignment**  
   A full lentivirus alignment was assembled, including representative HIV, SIV, and pSIV sequences.

2. **Temporal signal assessment**  
   Each NRR was tested for temporal signal using **BETS**, comparing isochronous and heterochronous BEAST models via stepping-stone sampling.

3. **Short-term rate estimation**  
   HIV-1 group M was used as the best-characterised time-stamped dataset to estimate short-term substitution rates for each NRR.

4. **PoW transformation**  
   For NRRs with sufficient temporal signal, ultrametric genetic-distance trees were transformed into **PoW time-scaled divergence trees** using the inferred short-term rate distributions.

---

## Notes

- The XML files provided here are those relevant to the analyses presented in the manuscript.
- The PoW-transformed trees are provided only for NRRs that showed evidence of temporal signal and were used in the main dating analyses.
- The repository is intended to support reproducibility of the principal results reported in the paper.
