############################################################
# PoW transformation of BEAST posterior trees and rates
# Manuscript:
# "A unified evolutionary framework to reconstruct the viral histories of human, simian, and prosimian immunodeficiency viruses across timescales"
# This script:
# 1. Reads posterior short-term substitution rates from BEAST logs
# 2. Reads BEAST ultrametric distance trees
# 3. Optionally adjusts the pSIVgml branch length for NRRs 2,5,6
# 4. Applies the PoW transformation to posterior trees
# 5. Writes PoW-transformed trees for downstream MCC summarisation
############################################################

############################
# Load required libraries
############################
library("stringr")
library("ape")
library("nleqslv")
library("treeio")


############################
# User-defined input/output paths
############################

# Posterior short-term substitution rate log file
# (*.log output file from BEAST v1.10.4)
filepath_rates <- "~/local/directory/NRR_1_dated_hetero_v1.log"

# Ultrametric distance trees (*.trees output file from BEAST v1.10.4)
filepath_distances <- "~/local/directory/NRR_1.trees"

# BEAST log file for the standard HKY substitution model used
# to construct the ultrametric distance trees
filepath_HKYsubstitutionModel <- "~/local/directory/NRR_1.log"

# Temporary sampled tree file
sampled_treefile <- "~/local/directory/sampled_trees.trees"

# Output path for PoW-transformed trees
output_pow_trees <- "~/local/directory/PoWtransformed_NRR_1.trees"

############################
# Settings
############################

# Number of posterior samples to draw after burn-in
sample_size <- 100

# Set seed for reproducibility of posterior subsampling
set.seed(123)

# Best-fit muMax for PoW transformation
# RNA viruses: 3.65e-2
# ssDNA viruses: 2.0e-2
# dsDNA viruses: 3.0e-3
muMax <- 3.65e-2
muMaxs <- rep(muMax, sample_size)

############################
# Read posterior logs
############################

# Read log file for inferred short-term rate
con <- readLines(filepath_rates)

# Read log file for HKY substitution model parameters
con3 <- readLines(filepath_HKYsubstitutionModel)

# Locate key columns
meanrate_column <- which(strsplit(con[5], "\t")[[1]] == "meanRate")
kappa_column <- which(strsplit(con3[3], "\t")[[1]] == "kappa")
frequencies1_column <- which(strsplit(con3[3], "\t")[[1]] == "frequencies1")

# Remove first 10% as burn-in
con <- con[5:length(con)]
con <- con[(round(length(con) / 10) + 1):length(con)]

con3 <- con3[5:length(con3)]
con3 <- con3[(round(length(con3) / 10) + 1):length(con3)]

############################
# Sample posterior states
############################

sampled_states <- sample(con, sample_size)
sampled_states2 <- sample(con3, sample_size)

############################
# Extract rates and HKY parameters
############################

rates <- c()
freq1 <- c()
freq2 <- c()
freq3 <- c()
freq4 <- c()
K <- c()

for (i in seq_along(sampled_states)) {
  rates <- c(rates, as.numeric(strsplit(sampled_states[[i]], "\t")[[1]][[meanrate_column]]))
  freq1 <- c(freq1, as.numeric(strsplit(sampled_states2[[i]], "\t")[[1]][[frequencies1_column]]))
  freq2 <- c(freq2, as.numeric(strsplit(sampled_states2[[i]], "\t")[[1]][[frequencies1_column + 1]]))
  freq3 <- c(freq3, as.numeric(strsplit(sampled_states2[[i]], "\t")[[1]][[frequencies1_column + 2]]))
  freq4 <- c(freq4, as.numeric(strsplit(sampled_states2[[i]], "\t")[[1]][[frequencies1_column + 3]]))
  K <- c(K, as.numeric(strsplit(sampled_states2[[i]], "\t")[[1]][[kappa_column]]))
}

# Quick summary of sampled rates
median(rates)
min(rates)
max(rates)

# Uncomment below if you want to add variation around muMax
# muMax_stdev <- 10.0e-3
# muMaxs <- rnorm(sample_size, muMax, muMax_stdev)

############################
# Read BEAST posterior tree file
############################

con2 <- readLines(filepath_distances)

# Locate first tree line
for (i in seq_along(con2)) {
  if (paste(strsplit(con2[[i]], "")[[1]][1:4], collapse = "") == "tree") {
    break
  }
}

# Extract posterior tree block and remove first 10% as burn-in
tree_list <- con2[i:(length(con2) - 1)]
tree_list <- tree_list[(round(length(tree_list) / 10) + 1):length(tree_list)]

# Randomly subsample trees
sampled_trees <- sample(tree_list, sample_size)

# Recreate sampled tree file
recreated_treefile <- c(con2[1:(i - 1)], sampled_trees, "End;")
write.table(recreated_treefile, sampled_treefile,
            sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Import sampled trees
importedtree <- read.nexus(sampled_treefile, tree.names = NULL, force.multi = FALSE)

############################
# Helper function
############################

# Find parent node from tip label
findParentNode <- function(tree, tipLabel) {
  if (tipLabel %in% tree$tip.label) {
    nodeId <- which(tree$tip.label == tipLabel)
    parentNode <- tree$edge[which(tree$edge[, 2] == nodeId), 1]
    return(parentNode)
  } else {
    stop("Tip label not found in the tree")
  }
}

############################################################
# pSIVgml branch adjustment
# Only run for NRRs that include pSIVgml (NRRs 2, 5, 6)
############################################################

for (el in 1:sample_size) {
  print(el)
  
  muMax <- muMaxs[[el]]
  meanRate <- rates[[el]]
  pA <- freq1[[el]]
  pC <- freq2[[el]]
  pG <- freq3[[el]]
  pT <- freq4[[el]]
  kap <- K[[el]]
  
  pY <- pT + pC
  pR <- pA + pG
  delta <- 0.1
  
  betaMax <- muMax / (2 * ((pT * pC + pA * pG) * kap + pY * pR))
  steps <- seq(0, log10(betaMax) + 9, delta)
  steps2 <- seq(1, length(steps), 1)
  
  lambda <- function(l) {
    final <- sum(exp(l * (1 / delta * steps + 1)) *
                   (2 * ((pT * pC + pA * pG) * kap + pY * pR)) *
                   10^(-9 + steps)) /
      sum(exp(l * steps2)) - meanRate
    return(final)
  }
  
  LAMBDA <- nleqslv(0, lambda)$x
  proportions <- exp(LAMBDA * (1 / delta * steps + 1)) / sum(exp(LAMBDA * steps2))
  ratespread <- 10^(-9 + steps)
  
  e2 <- function(t) exp(-ratespread * t)
  e3 <- function(t) exp(-(pR * kap + pY) * ratespread * t)
  e4 <- function(t) exp(-(pY * kap + pR) * ratespread * t)
  pTC <- function(t) pC + (pC * pR) / pY * e2(t) - pC / pY * e4(t)
  pAG <- function(t) pG + (pG * pY) / pR * e2(t) - pG / pR * e3(t)
  pTA <- function(t) pA * (1 - e2(t))
  pTG <- function(t) pG * (1 - e2(t))
  pCA <- function(t) pA * (1 - e2(t))
  pCG <- function(t) pG * (1 - e2(t))
  
  s1 <- function(t) sum(proportions * 2 * pT * pTC(t))
  s2 <- function(t) sum(proportions * 2 * pA * pAG(t))
  v  <- function(t) sum(proportions * (2 * pT * pTA(t) + 2 * pT * pTG(t) +
                                         2 * pC * pCA(t) + 2 * pC * pCG(t)))
  
  a1 <- function(t) -log(1 - (pY * s1(t)) / (2 * pT * pC) - v(t) / (2 * pY))
  a2 <- function(t) -log(1 - (pR * s2(t)) / (2 * pA * pG) - v(t) / (2 * pR))
  b  <- function(t) -log(1 - v(t) / (2 * pY * pR))
  k1 <- function(t) (a1(t) - pR * b(t)) / (pY * b(t))
  k2 <- function(t) (a2(t) - pY * b(t)) / (pR * b(t))
  
  # Conservative estimate of pSIVgml age
  tstar <- 4.2e6
  
  func_dist <- function(dt) {
    final <- 2 * pT * pC / pY * (a1(tstar) - pR * b(tstar)) +
      (2 * pA * pG / pR) * (a2(tstar) - pY * b(tstar)) +
      2 * pY * pR * b(tstar) - dt
    return(final)
  }
  
  added_branch <- nleqslv(0, func_dist)$x
  
  # Identify stem branch to SIVs
  edgeMatrix <- importedtree[[el]]$edge
  index <- which(edgeMatrix[, 1] == findParentNode(importedtree[[el]], "pSIVglm"))
  
  # Identify stem branch to pSIVglm
  pSIV_index <- which.max(importedtree[[el]]$edge.length)
  SIV_index <- setdiff(index, pSIV_index)
  
  # Extend stem branches of SIVs and pSIV
  importedtree[[el]]$edge.length[[SIV_index]]  <- importedtree[[el]]$edge.length[[SIV_index]] + added_branch / 2
  importedtree[[el]]$edge.length[[pSIV_index]] <- importedtree[[el]]$edge.length[[pSIV_index]] + added_branch / 2
}

############################################################
# Apply PoW transformation to sampled trees
############################################################

for(el in 1:sample_size){
  print(el)
  muMax = muMaxs[[el]]
  meanRate = rates[[el]]
  pA = freq1[[el]]
  pC = freq2[[el]]
  pG = freq3[[el]]
  pT = freq4[[el]]
  kap = K[[el]]
  pY = pT + pC
  pR = pA + pG
  delta = 0.1
  betaMax = muMax/(2*((pT*pC + pA*pG)*kap + pY*pR))
  steps = seq(0,log10(betaMax)+9,delta)
  steps2 = seq(1,length(steps),1)
  lambda <- function(l){
    final <- sum(exp(l*(1/delta*steps+1))*(2*((pT*pC + pA*pG)*kap + pY*pR))*10**(-9+steps))/sum(exp(l*steps2)) - meanRate
    return (final)
  }
  
  LAMBDA = nleqslv(0,lambda)$x
  proportions = exp(LAMBDA*(1/delta*steps+1))/sum(exp(LAMBDA*steps2))
  ratespread = 10**(-9+steps)
  
  e2 <- function(t) exp(-ratespread*t)
  e3 <- function(t) exp(-(pR*kap+pY)*ratespread*t)
  e4 <- function(t) exp(-(pY*kap+pR)*ratespread*t)
  pTC <- function(t) pC+(pC*pR)/pY*e2(t)-pC/pY*e4(t)
  pAG <- function(t) pG+(pG*pY)/pR*e2(t)-pG/pR*e3(t)
  pTA <- function(t) pA*(1-e2(t))
  pTG <- function(t) pG*(1-e2(t))
  pCA <- function(t) pA*(1-e2(t))
  pCG <- function(t) pG*(1-e2(t))
  s1 <- function(t) sum(proportions*2*pT*pTC(t))
  s2 <- function(t) sum(proportions*2*pA*pAG(t))
  v <- function(t) sum(proportions*(2*pT*pTA(t)+2*pT*pTG(t)+2*pC*pCA(t)+2*pC*pCG(t)))
  a1 <- function(t) -log(1-(pY*s1(t))/(2*pT*pC)-v(t)/(2*pY))
  a2 <- function(t) -log(1-(pR*s2(t))/(2*pA*pG)-v(t)/(2*pR))
  b <- function(t) -log(1-v(t)/(2*pY*pR))
  k1 <- function(t) (a1(t)-pR*b(t))/(pY*b(t))
  k2 <- function(t) (a2(t)-pY*b(t))/(pR*b(t))
  func <- function(t){
    final <- 2*pT*pC/pY*(a1(t)-pR*b(t))+(2*pA*pG/pR)*(a2(t)-pY*b(t))+2*pY*pR*b(t) - totalBranchLength
    return(final)
  }
  
  
  allBranches=mapply(c,importedtree[[el]]$edge[,1],importedtree[[el]]$edge[,2],(importedtree[[el]]$edge.length),SIMPLIFY = FALSE)
  
  Allpaths <- list(list())
  for(j in 0:importedtree[[el]]$Nnode+1){
    branches <- list()
    for(i in 1:(length(nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2))-1)){
      branches[[i]] <- c(nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2)[i+1],nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2)[i])
    }
    Allpaths[[j]] <- c(branches)
  }
  
  trackBranches <- list()
  globcount=0
  for(k in 1:length(Allpaths)){
    totalBranchLength=0
    count=0
    convertedHeight<-list()
    for(i in 1:length(Allpaths[[k]])){
      for(j in 1:length(allBranches)){
        if(sum(allBranches[[j]] %in% Allpaths[[k]][[i]],na.rm = TRUE)==2){
          matches=0
          count = count + 1
          globcount = globcount + 1
          trackBranches[[globcount]] <- c(allBranches[[j]][c(1:2)])
          for(m in 1:length(trackBranches)){
            if(sum(trackBranches[[m]] %in% Allpaths[[k]][[i]],na.rm = TRUE)==2){
              matches = matches + 1
            }
          }
          if(count==1 && matches==1){
            totalBranchLength = totalBranchLength + allBranches[[j]][[3]]
            convertedHeight[[count]] <- c(nleqslv(0,func)$x)
            importedtree[[el]]$edge.length[[j]]=convertedHeight[[count]]
          }
          if(count>1 && matches==1){
            totalBranchLength = totalBranchLength + allBranches[[j]][[3]]
            convertedHeight[[count]] <- c(nleqslv(0,func)$x)
            importedtree[[el]]$edge.length[[j]]=convertedHeight[[count]]-convertedHeight[[count-1]]
          }
          if(count>1 && matches>1){
            TRUE
          }
        }
      }
    }
  }
}

############################
# Save PoW-transformed trees
############################

ape::write.nexus(importedtree, file = output_pow_trees)


############################################################
# Example code for reading MCC tree and extracting node ages
############################################################

tips_hiv1M <- c(
  "AF077336|VI850|HIV-1-F1|BELGIUM|1993",
  "M62320|U455_U455A|HIV-1-A1|UGANDA|1985",
  "U21135|WEAU160_GHOSH|HIV-1-B|UNITED_STATES|1990",
  "U46016|ETH2220|HIV-1-C|ETHIOPIA|1986",
  "U88822|84ZR085|HIV-1-D|DEM_REP_OF_CONGO|1984"
)

tips_hiv1O <- c(
  "AY169814|99USTWLA|HIV-1-O|UNITED_STATES|1999",
  "AB485666|BCF06|HIV-1-O|CAMEROON|1994",
  "L20571|MVP5180|HIV-1-O|CAMEROON|1991",
  "AF407418|VAU|HIV-1-O|FRANCE|1992",
  "AJ302646|99SE-MP1299|HIV-1-O|SENEGAL|1999",
  "AB485668|I_2478B|HIV-1-O|UNITED_STATES|",
  "L20587|ANT70|HIV-1-O|BELGIUM|1987",
  "AY169803|96CMA102|HIV-1-O|CAMEROON|1996"
)

tips_hiv1N <- c(
  "GQ324962|U14296|HIV-1-N|CAMEROON|2006",
  "GQ324959|SJGddd|HIV-1-N|CAMEROON|2002",
  "GQ324958|U14842|HIV-1-N|CAMEROON|2006",
  "AY532635|DJO0131|HIV-1-N|CAMEROON|2002",
  "AJ006022|YBF30|HIV-1-N|CAMEROON|1995"
)

tips_hiv2A <- c(
  "AF082339|ALI|HIV-2-A|PORTUGAL|",
  "M31113|ST_JSP4_27|HIV-2-A|SENEGAL|1986",
  "M30502|BEN|HIV-2-A|GERMANY|"
)

tips_hiv2B <- c(
  "U27200|EHO|HIV-2-B|COTE_DIVOIRE|",
  "X61240|D205_ALT|HIV-2-B|GHANA|1986"
)

tips_sivcpz <- c(
  "U42720|ANT|SIV-CPZ|DEM_REP_OF_CONGO|1990",
  "KP861923|LB715|SIV-CPZ|CAMEROON|2005",
  "JX178450|LB715|SIV-CPZ|CAMEROON|2005",
  "JQ866001|BF1167|SIV-CPZ|DEM_REP_OF_CONGO|2006",
  "JQ768416|SIVcpzTAN13|SIV-CPZ|TANZANIA|2006",
  "JN091691|TAN5|SIV-CPZ|TANZANIA|2006",
  "JN091690|UG38|SIV-CPZ|TANZANIA|2009",
  "GQ217539|SIVcpzGab4|SIV-CPZ|GABON|2006",
  "EF394358|TAN3|SIV-CPZ|TANZANIA|2002",
  "EF394357|TAN2|SIV-CPZ|TANZANIA|2001",
  "EF394356|TAN1|SIV-CPZ|TANZANIA|2000",
  "AY169968|SIVcpzCAM13|SIV-CPZ|CAMEROON|2001",
  "AF447763|TAN1|SIV-CPZ|TANZANIA|2000",
  "AF382828|SIVcpzGAB2|SIV-CPZ|GABON|1988",
  "AF103818|US_Marilyn|SIV-CPZ|UNITED_STATES|1985"
)

tips_sivsmm <- c(
  "U48814|SL92E|SIV-SMM|SIERRA_LEONE|1992",
  "U48812|SL92C|SIV-SMM|SIERRA_LEONE|1992",
  "U48811|SL92b|SIV-SMM|SIERRA_LEONE|1992",
  "U17646|SL92A|SIV-SMM|SIERRA_LEONE|1992",
  "KC693496|sm487-gag|SIV-SMM|COTE_DIVOIRE|2006",
  "JX860433|SIVsmSL92B|SIV-SMM|SIERRA_LEONE|1992",
  "JX860432|SIVsmSL92A|SIV-SMM|SIERRA_LEONE|1992",
  "JX860431|SIVsmLIB1|SIV-SMM|LIBERIA|1989",
  "JX860430|SIVsmCI2|SIV-SMM|COTE_DIVOIRE|1979",
  "AY932819|TAI_37|SIV-SMM|COTE_DIVOIRE|",
  "AY932818|CI_8|SIV-SMM|COTE_DIVOIRE|",
  "AY932812|TAI_33|SIV-SMM|COTE_DIVOIRE|",
  "AY932810|TAI_29|SIV-SMM|COTE_DIVOIRE|",
  "AY932809|YAI_31|SIV-SMM|COTE_DIVOIRE|",
  "AY932806|TAI_1|SIV-SMM|COTE_DIVOIRE|",
  "AY864798|SIVsm_SL93_135|SIV-SMM|SIERRA_LEONE|",
  "AY864797|SIVsm_SL93_134|SIV-SMM|SIERRA_LEONE|",
  "AY864796|SIVsm_SL93_119|SIV-SMM|SIERRA_LEONE|",
  "AY864795|SIVsm_SL93_080|SIV-SMM|SIERRA_LEONE|",
  "AF334679|SL92B|SIV-SMM|SIERRA_LEONE|1992"
)

# Example MCC tree input
mcc_tree <- read.beast("~/local/directory/MCC_PoWtransformed_NRR_1.tree")

# Remove pSIV for some summaries
t <- drop.tip(mcc_tree, "pSIVglm")
t@data$posterior[t@data$posterior < 0.9] <- 0

write.beast(
  t,
  file = "~/local/directory/MCC_filtered_PoWtransformed_NRR_1.tree"
)

# Extract node ages
t1_hiv1M <- mcc_tree@data$height_median[match(toString(MRCA(mcc_tree, match(tips_hiv1M, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]
t2_hiv1M <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, match(tips_hiv1M, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]

t1_hiv1O <- mcc_tree@data$height_median[match(toString(MRCA(mcc_tree, match(tips_hiv1O, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]
t2_hiv1O <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, match(tips_hiv1O, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]

t1_hiv1N <- mcc_tree@data$height_median[match(toString(MRCA(mcc_tree, match(tips_hiv1N, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]
t2_hiv1N <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, match(tips_hiv1N, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]

t1_hiv2A <- mcc_tree@data$height_median[match(toString(MRCA(mcc_tree, match(tips_hiv2A, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]
t2_hiv2A <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, match(tips_hiv2A, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]

t1_hiv2B <- mcc_tree@data$height_median[match(toString(MRCA(mcc_tree, match(tips_hiv2B, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]
t2_hiv2B <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, match(tips_hiv2B, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]

t1_sivcpz <- mcc_tree@data$height_median[match(toString(MRCA(mcc_tree, match(tips_sivcpz, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]
t2_sivcpz <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, match(tips_sivcpz, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]

t1_sivsmm <- mcc_tree@data$height_median[match(toString(MRCA(mcc_tree, na.omit(match(tips_sivsmm, mcc_tree@phylo$tip.label)))), mcc_tree@data$node)]
t2_sivsmm <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, na.omit(match(tips_sivsmm, mcc_tree@phylo$tip.label)))), mcc_tree@data$node)]

t1_root <- mcc_tree@data$height_median[Ntip(mcc_tree@phylo) + Nnode(mcc_tree@phylo)]
t2_root <- mcc_tree@data$height_0.95_HPD[Ntip(mcc_tree@phylo) + Nnode(mcc_tree@phylo)]

# Without pSIV
dates <- data.frame(
  c("median", "lower_95%", "upper_95%"),
  c(t1_hiv1M, t2_hiv1M[[1]][1], t2_hiv1M[[1]][2]),
  c(t1_hiv1O, t2_hiv1O[[1]][1], t2_hiv1O[[1]][2]),
  c(t1_hiv1N, t2_hiv1N[[1]][1], t2_hiv1N[[1]][2]),
  c(t1_hiv2A, t2_hiv2A[[1]][1], t2_hiv2A[[1]][2]),
  c(t1_hiv2B, t2_hiv2B[[1]][1], t2_hiv2B[[1]][2]),
  c(t1_sivcpz, t2_sivcpz[[1]][1], t2_sivcpz[[1]][2]),
  c(t1_sivsmm, t2_sivsmm[[1]][1], t2_sivsmm[[1]][2]),
  c(t1_root, t2_root[[1]][1], t2_root[[1]][2]),
  c(0, 0, 0)
)

write.csv(dates, "~/local/directory/dates_NRR_1.csv", row.names = TRUE)

# With pSIV
tips_sivs <- setdiff(mcc_tree@phylo$tip.label, "pSIVglm")
t1_sivs <- mcc_tree@data$height[match(toString(MRCA(mcc_tree, match(tips_sivs, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]
t2_sivs <- mcc_tree@data$height_0.95_HPD[match(toString(MRCA(mcc_tree, match(tips_sivs, mcc_tree@phylo$tip.label))), mcc_tree@data$node)]

dates <- data.frame(
  c("median", "lower_95%", "upper_95%"),
  c(t1_hiv1M, t2_hiv1M[[1]][1], t2_hiv1M[[1]][2]),
  c(t1_hiv1O, t2_hiv1O[[1]][1], t2_hiv1O[[1]][2]),
  c(t1_hiv1N, t2_hiv1N[[1]][1], t2_hiv1N[[1]][2]),
  c(t1_hiv2A, t2_hiv2A[[1]][1], t2_hiv2A[[1]][2]),
  c(t1_hiv2B, t2_hiv2B[[1]][1], t2_hiv2B[[1]][2]),
  c(t1_sivcpz, t2_sivcpz[[1]][1], t2_sivcpz[[1]][2]),
  c(t1_sivsmm, t2_sivsmm[[1]][1], t2_sivsmm[[1]][2]),
  c(t1_sivs, t2_sivs[[1]][1], t2_sivs[[1]][2]),
  c(t1_root, t2_root[[1]][1], t2_root[[1]][2])
)
