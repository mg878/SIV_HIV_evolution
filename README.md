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

- **`NRR_*.xml`**  
  BEAST XML files used to infer the **ultrametric genetic-distance trees** for each NRR under the HKY substitution model. These trees provide the input distance trees for the downstream PoW transformation

These XML files correspond to the analyses used to generate the posterior tree files from which the PoW-transformed time-scaled divergence trees were constructed.

- **`NRR_*_dated_iso_v1.xml`**  
  Isochronous analysis for the relevant NRR.

- **`NRR_*_dated_hetero_v1.xml`**  
  Heterochronous (tip-dated) analysis for the relevant NRR.

These two sets of XML files were used for **Bayesian Evaluation of Temporal Signal (BETS)** via stepping-stone sampling to compare isochronous and heterochronous models and assess temporal signal in each NRR.

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

2. **Short-term rate estimation from HIV-1 group M**  
   A time-stamped HIV-1 group M dataset was used to estimate short-term substitution rates for each non-recombinant region (NRR), providing the short-term rate scale used in downstream analyses.

3. **Temporal signal assessment in the full HIV+SIV dataset**  
   The HIV-derived short-term rates were then used to define weakly informative priors (with standard deviations inflated 2.5× relative to the HIV posterior for each NRR) for the full HIV+SIV dataset. For each NRR, temporal signal was assessed using **Bayesian Evaluation of Temporal Signal (BETS)** by comparing isochronous and heterochronous BEAST models via stepping-stone sampling.

4. **Ultrametric distance-tree inference**  
   For each NRR, BEAST XML files (`NRR_*.xml`) were used to infer ultrametric genetic-distance trees under the HKY substitution model. These trees were used as input for the PoW transformation.

5. **PoW transformation**  
   For NRRs with sufficient temporal signal, ultrametric genetic-distance trees were transformed into **PoW time-scaled divergence trees** using the inferred short-term rate distributions.
