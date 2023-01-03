####libraries required for the analyses
library("ape")
library("stats")
library("nleqslv")
library("scales")
library("treeio")
require(tidytree)
##########################################

####path to log files for the inferred short-term substitution rate of HIV and SIV (*.log output file from BEAST v1.10.4)

#uncomment below to use the HIV rate (for constructing a PoW-transformed HIV time tree)
#filepath_rates <- "path/to/directory/hiv_rate.log"

#uncomment below to use the SIV rate (for constructing a PoW-transformed SIV/pSIV time tree)
#filepath_rates <- "path/to/directory/siv_rate.log"

####path to ultrametric HKY distance tree (*.trees output file from BEAST v1.10.4)

#uncomment below to import the distance tree for the pol locus data set
#filepath_distances <- "path/to/directory/pol_locus.trees"

#uncomment below to import the distance tree for the integrase gene data set
#filepath_distances <- "path/to/directory/integrase_gene.trees"

#######

#######
#con reads each row in the log file, removing the first 10% of rows as burn-in
con = readLines(paste(filepath_rates, sep = ""))
con <- con[c(5:length(con))]
con <- con[c((round(length(con)/10)+1):length(con))]
#######

#sample_size sets the number of samples taken at random from the posterior rate distribution and distance trees.
sample_size=100

#######
#selecting states from the log file at random
sampled_states <- sample(con,sample_size)
#######

#######
#parameters of the substitution model

#uncomment below for pol locus if you are doing the PoW-transformation for HIV samples
#pA=0.402
#pC=0.19
#pG=0.186
#pT=0.222
#kap=8.625

#uncomment below for pol locus if you are doing the PoW-transformation for SIV/pSIV samples
#pA=0.284
#pC=0.248
#pT=0.251
#kap=4.319
#pG=0.218

#uncomment below for the integrase gene data set 
#pA=0.246
#pC=0.263
#pG=0.228
#pT=0.263
#kap=4.57
#######


#######
#reading the mean.Rate column for HIV and SIV dated samples
meanrate_column=15
rates <- c()
for(i in 1:length(sampled_states)){
  rates <- c(rates,as.numeric(strsplit(sampled_states[[i]],"\t")[[1]][[meanrate_column]]))
}
#######

#use the best fit value of muMax in RNA viruses for PoW transformation of every tree (do not change below)
muMax = 3.65*10**(-2)
muMaxs=rep(muMax,sample_size)

#uncomment below if you like to add variation around the value of muMax (not recommended using without strong evidence in favour of variation)
#muMax = 3.65*10**(-2)
#muMax_stdev=10.0*10**-3
#muMaxs=rnorm(sample_size,muMax,muMax_stdev)

#con2 reads the distance tree file
con2 = readLines(paste(filepath_distances, sep = ""))

for(i in 1:length(con2)){
  if(paste(strsplit(con2[[i]],'')[[1]][c(1:4)],collapse = '')=='tree'){
    break
  }
}

#tree_list includes a list of produced trees from BEAST, removing the first 10% as burn-in
tree_list <- con2[c(i:(length(con2)-1))] 
tree_list <- tree_list[c((round(length(tree_list)/10)+1):length(tree_list))]
sampled_trees <- sample(tree_list,sample_size)


#create a new tree file based on the subsampling in the previous step 
recreated_treefile <- c(con2[c(1:(i-1))],sampled_trees,'End;')
write.table(recreated_treefile, '/path/to/directory/sampled_trees.trees', sep = "\n", quote = FALSE, row.names = FALSE)

#import the tree file created in the previous step
filename <- "/path/to/directory/sampled_trees.trees"
importedtree <- read.nexus(filename, tree.names = NULL, force.multi = FALSE)
d <- read.beast(filename)


#only run the for loop below for the pol locus data set which includes the pSIV sample
#if running the integrase gene data set, jump to line #166
#calculate the tip distance for pSIV and re-distribute the added branch length 
for(el in 1:sample_size){
  print(el)
  muMax = muMaxs[[el]]
  meanRate = rates[[el]]
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
  
  #conservative estimate of the age of pSIVgml sample is 4.2 million years before present
  tstar<-4.2e6
  func_dist <- function(dt){
    final <- 2*pT*pC/pY*(a1(tstar)-pR*b(tstar))+(2*pA*pG/pR)*(a2(tstar)-pY*b(tstar))+2*pY*pR*b(tstar) - dt
    return(final)
  }
  added_branch<- c(nleqslv(0,func_dist)$x)
  #extend the stem branch of all SIVs
  importedtree[[el]]$edge.length[[2]]<- importedtree[[el]]$edge.length[[2]]+added_branch/2
  #extend the pSIV branch
  last_branch=1
  importedtree[[el]]$edge.length[[last_branch]]<- importedtree[[el]]$edge.length[[last_branch]]+added_branch/2
}

#transform each of the sampled distance trees using the PoW model 
for(el in 1:sample_size){
  print(el)
  muMax = muMaxs[[el]]
  meanRate = rates[[el]]
  
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

d=as.phylo(importedtree[[10]])
plot.phylo(d,show.tip.label = FALSE,use.edge.length = TRUE)

#write the PoW-transformed trees in output.trees
ape::write.nexus(importedtree, file='/path/to/directory/output.trees')

#to create the PoW-transormed tree of all SIV and HIV lineages shown in Figure 1,
#import the MCC_PoWTransformed_pol_locus_SIVrate.tree (the maximum clade credibility tree
#of PoW-transformed time tree for pol locus data set calibrated using the SIV posterior rate distribution)

poltree<-read.nexus(file="/path/to/directory/MCC_PoWTransformed_pol_SIVrate.tree")


###these are samples selected from the pol locus data set to construct the virus tree in Figure 1:

#JX860432|SIV_SMM|SIERRA_LEONE|1992
#AB485670|HIV-2_B|COTE_DIVOIRE|
#JX860430|SIV_SMM|COTE_DIVOIRE|1979
#J03654|HIV-2_A|GUINEA-BISSAU|1986
#AY965437|SIV_SMM|UNITED_STATES|1995
#AF334679|SIV_SMM|SIERRA_LEONE|1992
#KR862350|SIV_VER|SOUTH_AFRICA|2010
#LC114462|SIV_MAL|ZAMBIA|2009
#M30931|SIV_VER|GERMANY|
#U58991|SIV_TAN|UGANDA|
#M58410|SIV_GRV|ETHIOPIA|
#LM999945|SIV_TAN|CAMEROON|2001
#HQ378594|SIV_SAB|SENEGAL|2005
#KM378563|SIV_DRL|GERMANY|2011
#AF367411|SIV_MND-2|CAMEROON|1998
#AF349680|SIV_RCM|NIGERIA|
#AY169806|HIV-1_O|CAMEROON|1996
#KP004991|SIV_GOR|CAMEROON|2012
#AJ006022|HIV-1_N|CAMEROON|1995
#DQ373065|SIV_CPZ|CAMEROON|2005
#KP861923|SIV_CPZ|CAMEROON|2005
#U42720|SIV_CPZ|DEM_REP_OF_CONGO|1990
#DQ007902|HIV-1_B|CHINA|2002
#KF214240|SIV_COL|UGANDA|2010
#AF131870|SIV_SUN|GABON|1998
#AF188116|SIV_LST|DEM_REP_OF_CONGO|1988
#M27470|SIV_MND-1|GABON|
#FM165200|SIV_OLC|COTE_DIVOIRE|1997
#AM745105|SIV_WRC|COTE_DIVOIRE|1997
#KY497574|SIV_ASC|DEM_REP_OF_CONGO|2013
#AF468658|SIV_GSN|CAMEROON|1999
#AY340700|SIV_MUS-1|CAMEROON|2001
#EF070331|SIV_MUS-2|CAMEROON|2001
#KF304708|SIV_MUS-3|GABON|2011
#AY340701|SIV_MON|CAMEROON|1999
#KJ461716|SIV_ASC|UGANDA|2010
#AM182197|SIV_TAL|CAMEROON|2001
#AJ580407|SIV_DEN|DEM_REP_OF_CONGO|
#AY523865|SIV_DEB|CAMEROON|1999
#L06042|SIV_SYK|KENYA|

species<-c("AF334679","U42720","DQ373065","AJ006022","JX860432","AB485670", "JX860430", "J03654", "AY965437","KR862350",
           "LC114462", "M30931", "U58991", "M58410","LM999945", "HQ378594",
            "KM378563","AF367411","AF349680","AY169806","KP004991","KP861923",
            "DQ007902","KF214240","AF131870","AF188116","M27470","FM165200",
            "AM745105","KY497574","AF468658","AY340700","EF070331","KF304708",
            "AY340701","KJ461716","AM182197","AJ580407","AY523865","L06042")
ii<-sapply(species,grep,poltree$tip.label)

merge_ii<-c()
for (u in ii) {
  merge_ii <- c(merge_ii,u[1])
}

fullNames<-poltree$tip.label[merge_ii]

pr.tree<-drop.tip(poltree,setdiff(poltree$tip.label,fullNames))

ape::write.nexus(pr.tree, file='/path/to/directory/MCC_PoWTransformed_HIV_SIV_family.tree')
