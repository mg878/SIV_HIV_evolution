# SIV_HIV_evolution
Ghafari et al.: Reconstructing the evolutionary history of human, simian, and prosimian immunodeficiency viruses.

## Codes

`PoW_transformation.R`: This is the code to generate PoW-transformed trees for both the pol locus and integrase gene data sets.

## Output

`pol_locus.xml`: Contains the xml file for all samples in the pol locus data set.

`integrase_gene.xml`: Contains the xml file for all the samples in the integrase gene data set.

`hiv_rate.xml`: Contains the xml file of all the time-stamped HIV samples used to infer HIV substitution rate.

`siv_rate.xml`: Contains the xml file of all the time-stamped SIV samples used to infer SIV substitution rate.

`hiv_rate.log`: Contains the posterior rate distribution of HIV generated using BEAST.

`siv_rate.log`: Contains the posterior rate distribution of SIV generated using BEAST.

`host_tree.tree`: Contains the primate host tree.

`MCC_PoWtransformed_*.tree`: Contains the maximum clade credibility tree of PoW-transformed time trees for * = pol locus (constructed using either the HIV or SIV posterior rate distributions) and integrase gene data sets.  
