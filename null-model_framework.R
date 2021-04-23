

### Analyze non-random patterns (A.K.A. significant species associations) from presence/absence 
### community matrix 


##########################################################################################
#                                                                                        #
#    "C-Score FRAMEWORK"                                                                 #
#                                                                                        #
# Null model analysis - C-Scores - Bayes Correction - Blois Framework                    #
#                                                                                        #
# Modified by Zach Zbinden from:                                                         #
#  D'Amen, Gotelli & Guisan: Disentangling biotic interactions, environmental filters,   #
#  and dispersal limitation as drivers of species co-occurrence                          #
##########################################################################################

### Packages required: 
install.packages("BBmisc")
install.packages("adespatial")
install.packages("usdm")
install.packages("vegan")
install.packages("ade4")
install.packages("Rmisc")



###############################
###      1.0                  #
###      Import Data          #
###############################

## read in data file
## example file contains Site Names, Geographic Coordinates, Environmental Data, and Species Presences/Absences
data <- read.csv("DATA.csv", header = T, row.names=1, sep =",")

### Example data contains data from TWO sampling periods, choose which period to analyze:
##  USE *EITHER* 1974  - OR - 2014 variables

###############################
###           1974          ###
###############################
data1974 <- data[34:66,]
species<- data1974[,c(3:26)]
env <- data1974[,c(27:48)]
coords <- data1974[,c(1:2)]
###############################


###############################
###          2014          ####
###############################
data2014 <- data[1:33,]
species <- data2014[,c(3:26)]
env <- data2014[,c(27:48)]
coords <- data2014[,c(1:2)]
###############################


###########################################
# Standardize ENV variables               #
###########################################
# change variable names from original format to be easier to understand
# not necessary in all cases
names(env)[names(env) == "dis_m3_pyr"] <- "avg.discharge"
names(env)[names(env) == "run_mm_cyr"] <- "avg.runoff"
names(env)[names(env) == "inu_pc_clt"] <- "max.innundation"
names(env)[names(env) == "ria_ha_csu"] <- "river.area.reach"
names(env)[names(env) == "ria_ha_usu"] <- "river.area.upstream"
names(env)[names(env) == "riv_tc_csu"] <- "river.vol.reach"
names(env)[names(env) == "riv_tc_usu"] <- "river.vol.upstream"
names(env)[names(env) == "gwt_cm_cav"] <- "avg.groundwater.table.depth"
names(env)[names(env) == "ele_mt_cav"] <- "avg.elevation"
names(env)[names(env) == "slp_dg_cav"] <- "avg.slope"
names(env)[names(env) == "sgr_dk_rav"] <- "stream.gradient"
names(env)[names(env) == "pre_mm_cyr"] <- "avg.annual.precipitation"
names(env)[names(env) == "glc_pc_c02"] <- "deciduous.tree.cover"
names(env)[names(env) == "glc_pc_c04"] <- "evergreen.tree.cover"
names(env)[names(env) == "glc_pc_c13"] <- "herbaceous.cover"
names(env)[names(env) == "for_pc_cse"] <- "total.forest.cover"
names(env)[names(env) == "cly_pc_cav"] <- "clay.soil"
names(env)[names(env) == "slt_pc_cav"] <- "silt.soil"
names(env)[names(env) == "snd_pc_cav"] <- "sand.soil"
names(env)[names(env) == "soc_th_cav"] <- "carbon.soil"
names(env)[names(env) == "kar_pc_cse"] <- "karst.extent"
##################################################
# standardize by mean and unit standard deviation
##################################################
library(BBmisc)
summary(env)
env <-normalize(env, method="standardize", range =c(0,1), margin=2L)
summary(env)


################################################################################################################################
#  SKIP STEP 1.1 if you wish to use Geographic Coordinates rather than distance matrix
################################################################################################################################
### 1.1: Extension for Fluvial Analysis -- distance-based Moran's eigenvector maps (dbMEM)                                     #
### OPTIONAL STEP: Read in fluvial distance matrix and decompose into eigenvectors to represent space in mutliple dimensions   #
## Required data: distance matrix representing pairwise fluvial/network distance between all samples sites                     #
## Such a matrix can be acquired by using the Network Analyst Tool in ArcGis                                                   #
## Tutorial: http://www.biology.ualberta.ca/facilities/GIS/uploads/instructions/AVRiverDistanceMatrix9.pdf                     #
## **If you wish to use LAT and LONG rather than fluvial distance, skip to step 2.0                                            #
################################################################################################################################
############################
### Reading spatial data ###
############################
### Fluvial network distance was created in ArcGIS using Network Analyst Extension
hydro_dist <- read.csv("../co-occur_workdir/fluvial_dist.csv", header = T, sep = ",", row.names = 1)

###################################
### dbMEM method with adespatial  #
###################################
hydro_dist <- as.dist(hydro_dist, diag = TRUE, upper = TRUE)
### dbMEM method with adespatial
library(adespatial)
### Spatial eigenfunction analysis of fluvial distance matrix to generate eigenvectors representing "river space"
dbmem1 <- dbmem(hydro_dist, thresh = NULL, MEM.autocor = "positive", store.listw = TRUE, silent = FALSE)
# replace coords object (default includes X&Y/LONG&LAT) with the spatial eigenvectors
coords <- dbmem1
row.names(coords) <- row.names(species)


####################################################
### check for collinearity among spatial and env  #
###################################################
#### PCA for the environmental variables
npcs <- 6 # set the number of PCs you wish to retain
pc.env<-prcomp(env, center = FALSE, scale. = FALSE) 
pc.loadings <- pc.env$rotation
write.csv(pc.loadings,"./pca_loadings.csv")
summary(pc.env) # show variance explained by each component
plot(pc.env)   
env.pcs<-pc.env$x[,1:npcs] 
library(usdm)
comb <- cbind(dbmem1, env.pcs)
# check results, may consider removing vars if VIF too high (>10)
vif(comb)


# plots



###############################################################
###############################################################
##                                                           ##
##     2.0: Create function for calculating C-scores         ##
## Pairwise CO-OCCURENCE ANALYSIS with C-score  calculation  ##
##					                                                 ##
###############################################################
###############################################################
##
## 
##        DIRECTIONS: Highlight and Run ALL CODE FOR STEP 2.0 (NO CHANGES NEEDED)
##
##
##
##  ©  C. RANDIN and M. D'Amen, Dept.of Ecology & Evolution, University of Lausanne                                       
#############################################################################################
## 
## Format required: a plots (rows) x species (columns) matrix of presences/absences
## Input matrix should have column names (species names) and row names (sampling plots)
##
## The function c.score calculates the C-score matrix to detect species association, for the whole community and for species pairs
## Randomization: column sum is fixed

## It returns the C-score index for the observed community (ObsCscoreTot), p.value (PValTot) and standardized effect size (SES.Tot). It saves also a table in the working directory where the same 
## metrics are calculated for each species pair (only the table with species pairs with significant p.values is saved in this version)

## NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli & McCabe 2002 - Ecology)
##
## Literature 
## Gotelli, N.J. and D.J. McCabe. 2002. Ecology, 83, 2091-2096.
## Gotelli, N.J. and Ulrich, W. 2010. Oecologia, 162, 463-477
## Stone, L. & Roberts, A. 1990. Oecologia, 85, 74-79


library(vegan)
library(ade4)

# ecospat.Cscore01(data.in, 100 ,outpath)

ecospat.Cscore01 <- function(data.in,npermut,outpath)
  
{
  
  # C-coef Observed matrix
  cat("Computing observed co-occurence matrix", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)  
  
  
  spec.occ <- data.matrix(data.in)
  
  ####  C score #############	
  coocc<-t(spec.occ)%*%spec.occ  # nb of checkboard units
  n.spec=dim(coocc)[1]			 
  mat1<-array(apply(spec.occ,MAR=2,sum),dim=c(n.spec,n.spec)) 
  mat2<-t(array(apply(spec.occ,MAR=2,sum),dim=c(n.spec,n.spec)))
  mat.obs.c.coef <- ((mat1 - coocc)*(mat2 - coocc))/(mat1*mat2) # observed c score 
  df.obs.c.coef <- data.frame(Col = rep(1:ncol(mat.obs.c.coef),each=ncol(mat.obs.c.coef)),Row = rep(1:nrow(mat.obs.c.coef),
                                                                                                    nrow(mat.obs.c.coef)),Sp1 = rep(colnames(mat.obs.c.coef),each=ncol(mat.obs.c.coef)),
                              Sp2 = rep(rownames(mat.obs.c.coef),nrow(mat.obs.c.coef)),Co.Occ = c(mat.obs.c.coef)) # dataframe with cscore for each species pair
  v.diago.inf <- c(rownames(df.obs.c.coef)[df.obs.c.coef[,1]>df.obs.c.coef[,2]],rownames(df.obs.c.coef)[df.obs.c.coef[,1]==df.obs.c.coef[,2]])# Remove identical combinations of species 
  df.obs.c.coef <- df.obs.c.coef[-as.numeric(v.diago.inf),]
  CscoreTot<-mean(df.obs.c.coef$Co.Occ)
  
  # Matrix to store the permuations
  mat.perm <- matrix(0,nrow(df.obs.c.coef),npermut, dimnames = list(c(paste(df.obs.c.coef[,3],df.obs.c.coef[,4])),c(1:npermut)))
  
  
  # Permutations C-coef 
  
  cat("Computing permutations", "\n",append = F)
  cat(".............", "\n",append = F)
  
  for (i in 1:npermut)
  {
    if (i == 1)
    {
      cat(npermut ," permutations to go", "\n",append = F)
      cat(".............", "\n",append = F)			
    }	
    if (i == npermut / 2)
    {
      cat(npermut / 2," permutations to go", "\n",append = F)
      cat(".............", "\n",append = F)
    }	
    
    
    spec.occ.perm1<-data.matrix(data.in)
    spec.occ.perm1 <- permatswap(spec.occ.perm1,fixedmar="both",mtype="prab",time=1) # row/column sums are preserved
    # time=1 : separate swapping sequence that always begins with the original matrix
    spec.occ.perm <- as.matrix(spec.occ.perm1[[3]][[1]] )
    
    coocc.perm <- t(spec.occ.perm)%*%spec.occ.perm 
    mat1.perm <- array(apply(spec.occ.perm,MAR=2,sum),dim=c(n.spec,n.spec))
    mat2.perm <- t(array(apply(spec.occ.perm,MAR=2,sum),dim=c(n.spec,n.spec)))
    
    mat.obs.c.coef.perm <- ((mat1.perm - coocc.perm)*(mat2.perm - coocc.perm))/(mat1.perm*mat2.perm)
    
    df.obs.c.coef.perm <- data.frame(Col = rep(1:ncol(mat.obs.c.coef.perm),each=ncol(mat.obs.c.coef.perm)),Row = rep(1:nrow(mat.obs.c.coef.perm),
                                                                                                                     nrow(mat.obs.c.coef.perm)),Sp1 = rep(colnames(mat.obs.c.coef),each=ncol(mat.obs.c.coef.perm)), 
                                     Sp2 = rep(rownames(mat.obs.c.coef),nrow(mat.obs.c.coef.perm)),Co.Occ = c(mat.obs.c.coef.perm))
    
    # Remove identical combinations of species (same Co-occ coef) and the diagonal (Co-occ coeff = 0)	
    df.obs.c.coef.perm <- df.obs.c.coef.perm[-as.numeric(v.diago.inf),]
    
    # Store result of permuation
    mat.perm[,i] <- df.obs.c.coef.perm[,5]
    
  }
  
  
  ## for the whole community
  vec.CScore.tot<-as.vector(apply(mat.perm,MAR=2,mean)) # C-score for all null communities (mean on the columns)
  SimulatedCscore<-mean(vec.CScore.tot) # mean of Simulation C-score: Simulated C-score
  sd.SimulatedCscore<-sd(vec.CScore.tot) # standard deviation of null communities
  Zscore<-(CscoreTot-SimulatedCscore)/sd.SimulatedCscore # standardized effect size
  
  randtest.less<-as.randtest(vec.CScore.tot, CscoreTot, alter="less")
  pval.less<-randtest.less$pvalue
  randtest.greater<-as.randtest(vec.CScore.tot, CscoreTot, alter="greater")
  pval.greater<-randtest.greater$pvalue
  plot(randtest.greater, xlab= "Simulated C-scores",main=paste("", sep=""))
  
  
  # Calculate P-values based on random distribution
  mat.pval <- matrix(0,nrow(mat.perm),4,dimnames = list(rownames(mat.perm),c("Obs.Co.Occ","Zscore","pval_less","pval_greater")))
  mat.pval[,1] <- df.obs.c.coef[,5]
  
  cat("Computing P-values", "\n",append = F)
  cat(".............", "\n",append = F)
  
  for (k in 1:nrow(mat.perm))
  {
    
    mat.pval[k,2]<-	(df.obs.c.coef[k,5]-mean(mat.perm[k,]))/sd(mat.perm[k,])
    
    randtest<-as.randtest(sim=mat.perm[k,], obs=df.obs.c.coef[k,5], alter="less")
    mat.pval[k,3]<-randtest$pvalue
    randtest<-as.randtest(sim=mat.perm[k,], obs=df.obs.c.coef[k,5], alter="greater")
    mat.pval[k,4]<-randtest$pvalue        
  }
  
  # Exporting Co-occ matrix
  cat("Exporting dataset", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)  
  
  hist(as.vector(mat.pval[,2]), xlab="Zscore", main = paste(""))
  abline(v=c(2,-2),col = "red")
  
  mat.pval.names<-data.frame(df.obs.c.coef[,3:4],mat.pval,df.obs.c.coef.perm[,5])
  mat.pval.names2<-data.frame(mat.pval.names[,1:3],mat.pval.names[,7],mat.pval.names[,4:6])
  names(mat.pval.names2)[3]<-"obs.C-score"
  names(mat.pval.names2)[4]<-"exp.C-score"
  write.table(mat.pval.names2,file=paste(outpath,"\\Cscores01.txt", sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  
  tab<-mat.pval.names2
  v<-c(0)
  for (i in 1:nrow(tab)){
    if (tab[i,6]<=0.05||tab[i,7]<=0.05){
      v<-c(v,i)
    }
  }
  m<-data.frame()
  for(j in 1:length(v)){
    m<-rbind(m,tab[v[j],])
  }
  
  m1<-na.omit(m)
  
  write.table(m1,file=paste(outpath,"\\Sign_Cscores01.txt",sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  write.table(mat.perm, file=paste(outpath,"\\MatrixPermutations01.txt",sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  
  l<-list(ObsCscoreTot=CscoreTot, SimCscoreTot=SimulatedCscore, PVal.less=pval.less,PVal.greater=pval.greater,Z.score=Zscore)
  
  return(l)
  
  cat("Computations finished!", "\n",append = F)
  

}
#######################################
#######################################
###### END OF FUNCTION FROM STEP 2.0  #
#######################################
#######################################








########################################################################################################################
########################################################################################################################
##                                                                                                                    ##      
##                                                                                                                    ##
##                      3.0: Calculate C-Scores & Correct for number of tests                                         ##
##                                                                                                                    ##
##  Empirical Bayes method to to identify statistically significant species pairs in a binary presence-absence matrix ##
##                                                                                                                    ##
##                                                                                                                    ##
########################################################################################################################
########################################################################################################################


##   NOTE: This step can require a long time to run and depends highly on the number of permutations chosen below. 
##     Highlight all code for step 3.0 and run, no changes are necessary


### D'Amen, Gotelli & Guisan: Disentangling biotic interactions, environmental filters, and dispersal 
### limitation as drivers of species co-occurrence

# Adapted by Manuela D'Amen from the procedures in the FORTRAN program PAIRS (Ulrich 2010)
# For a theoretical description see Gotelli & Ulrich2010 - Oecologia

## Format required: data object should be a plots (rows) x species (columns) matrix of presences/absences and
## it should have column names (species names) and row names (sampling plots)

## NOTE:
## This scritp runs after the ecospat.Cscore01.r function, using R Object names from the results table of this function
## Names should be updated when applying to co-occurrence results coming from another source  

## 3.1: 
## Apply the function from step 2 to calculate the rescaled C-score index 
## 
nperm <- 10000 # number of permutations may be adjusted here (default=9999)

## Calculate C-scores
cscore.01<-ecospat.Cscore01(species,nperm,getwd())
 


### If everything up until this point has been executed correctly, the remaining steps in section 3 
### can be highlighted an executed now. 

# extract information from the output tables, the observed C-score and the permutation table
coocc.tab.all<-read.table("Cscores01.txt",h=T)
obs.r<-data.frame(coocc.tab.all$obs.C.score)                          # Values of observed C-scores
pairs<-coocc.tab.all[,1:2]                                            # species pairs names
CooccProb<-read.table("MatrixPermutations01.txt", h=T)     # matrix with C-scores values for null communities
CooccProb.pairs<-cbind(pairs,CooccProb)                               # attribute species pairs name to these values
CooccProb.pairs.mean<-apply(CooccProb.pairs[,3:nperm],1,mean)         # calculate for each pair the mean and sd C-score across null pairs
CooccProb.pairs.sd<-apply(CooccProb.pairs[,3:nperm],1,sd)
CooccProb.stat<-data.frame(pairs,CooccProb.pairs.mean,CooccProb.pairs.sd)  

## 3.2:
## Assign each pairwise C-score to one of 22 evenly spaced bins spanning the interval from 0 to 1.
## 
library(Rmisc)
n.bin<-22

colnames(obs.r)<-"Cscore"
obs.r$bin<- cut(obs.r$Cscore, breaks=seq(0,1,by=1/n.bin),include.lowest=TRUE,labels=1:n.bin)  
obs.r<-cbind(obs.r,pairs)
## calculates how many pairs are present in each bin for oberved communities
obs.r1<-obs.r[order(obs.r[,2], obs.r[,1]),] 
obs.bin.summary<-as.data.frame(table(obs.r1$bin)) 
obs.bin.meanCS<-aggregate(obs.r1$Cscore, list(obs.r1$bin), mean)
names(obs.bin.summary)<-c("bin","n.pairs")

## 3.3:
## Calculation of the average number of species pairs from the null communities with different scores in each bin (table null.meanCS)
## This average represents the null expectation of the C-score for species pairs in each bin. 
## 

matrix.null.bins<-CooccProb 
matrix.null.bin.meanCS<-data.frame(seq(1:n.bin))
names(matrix.null.bin.meanCS)<-"Group.1"

for (i in 1:ncol(CooccProb)){
  null.bin<-cut(CooccProb[,i], breaks=seq(0,1,by=1/n.bin), include.lowest=TRUE, labels=1:n.bin)
  matrix.null.bins[,i]<-null.bin
  null.bin1<-cbind(CooccProb[,i],null.bin)   
  null.bin.meanCS<-aggregate(null.bin1[,1],list(null.bin1[,2]),mean)
  matrix.null.bin.meanCS<-merge(matrix.null.bin.meanCS, null.bin.meanCS, by="Group.1", all.x=T)}

tab.summary<-data.frame(matrix(data=NA,byrow=T,nrow=n.bin,ncol=nperm)) #table with the number of pairs in each bin for each null community

for (z in 1:nperm){                                                   
  col.summary<-as.data.frame(table(matrix.null.bins[,z]))
  tab.summary[,z]<-col.summary[,2]
  colnames(tab.summary)[z]<-z} 

matrix.null.bin.meanCS[is.na(matrix.null.bin.meanCS)] <- 0
tab.summary.meanCS<-apply(matrix.null.bin.meanCS[,2:(nperm+1)],1,mean, na.rm=T)
null.meanCS<-cbind(seq(1:n.bin),tab.summary.meanCS)
names(null.meanCS)<-c("bin", "null.meanCS")

## 3.4:
## Within each bin, calculation of the mean and the 95% confidence limit 
##

myci <- function(t) {                             # Funtion to calculate the confidence intervals (CI)
  n <- length(t) # n is the sample size
  se <- sd(t)/sqrt(n) # Find the standard error of the sample
  m <- mean(t) # Find the sample mean
  cv <- qt(0.975,df=n-1) # cv is a critical value for the t distribution. P( t > cv ) = 0.025 = P( t < -cv )
  c(m-cv*se,m+cv*se) # Return the 95% confidence interval
}

mean.pairs.bin<-(apply(tab.summary,1,mean))
CI.pairs.bin<-t(apply(tab.summary,1,myci))
sd.pairs.bin<-apply(tab.summary,1,sd)

stat.pairs.bin<-data.frame(seq(1:n.bin),mean.pairs.bin,CI.pairs.bin)
names(stat.pairs.bin)<-c("bin","MEAN","UP.CI","LOW.CI")

half.bin<-signif(((1/n.bin)/2),2)
bin.mean<-seq(0,1,by=1/n.bin)-half.bin
bin.mean<-bin.mean[2:(n.bin+1)]

#table with OBSERVED n pairs for each bin
obs.bin.summary2<-cbind(obs.bin.summary,bin.mean)
#table with MEAN and CI n pairs for each bin
stat.pairs.bin2<-cbind(stat.pairs.bin,bin.mean)

## Graph comparing how many pairs are present in each bin in average across the null communities (white cirles)
## and the number of observed pairs in the same bin (black circle)
plot(x=stat.pairs.bin2[,5],y=stat.pairs.bin2[,2], xlab="C-score",ylab="Number of pairs")
points(y=obs.bin.summary2[,2],x=obs.bin.summary2[,3], pch=16) 

## 3.5:
## Within each bin, order the species pairs by their scores and retain the species pairs with the largest scores that
## place them above the mean for the number of species pairs expected from the simulated distribution (Bayes M criterion).
## This step can be modified to apply the Bayes CL criterion using the values calculated above for the 95% confidence limit
##

tab.finalBayesM<-data.frame(matrix(ncol=ncol(obs.r1),nrow=0))
names(tab.finalBayesM)<-names(obs.r1)

for(j in 1:n.bin){
  
  obs.r2<-obs.r1[obs.r1$bin==j,]
  if((obs.bin.summary[j,2]>stat.pairs.bin[j,2])==T){   # if the observed number of species in the bin is higher 
    # than the mean number of species expected by null communities in the same bin
    exp.mean<-as.numeric(null.meanCS[j,2])             # we retain the pairs with a C-score higher than the mean C-score measured from null community in that bin
    obs.r3<-obs.r2[obs.r2$Cscore>exp.mean,]
    tab.finalBayesM<-rbind(tab.finalBayesM, obs.r3) 
  }}


## 3.6: 
## Further reduce the set of significant pairs by retaining only those that 
## are statistically significant in an individual test (simple CL criterion). 
##

BayesM_merge<-merge(tab.finalBayesM,coocc.tab.all,by.y=c("Sp1","Sp2"))              
sign.BayesM<-BayesM_merge[which(BayesM_merge$pval_less<0.05|BayesM_merge$pval_greater<0.05),]

write.table(sign.BayesM,"Sign.BayesM.txt", sep="\t")

###############################
###############################
###############################
#                             #
#  END BAYESIAN CORRECTION    #
#                             #
###############################
###############################
###############################








###########################################################
###########################################################
#########        	                      			     ########	
#########  		        4.0:  Blois framework        ########		
#########                                          ########
###########################################################
###########################################################

### D'Amen, Gotelli & Guisan: Disentangling biotic interactions, environmental filters, 
### and dispersal limitation as drivers of species co-occurrence
# Adapted by  Manuela D'Amen from : 
#
# Blois et al (2014) A framework for evaluating the influence of climate, dispersal limitation, and biotic interactions using fossil 
# pollen associations across the late Quaternary.  Ecography, 37, 1095-1108.
#
### Modified by Zach Zbinden
#### LOAD NECESSARY DATA
presabs<-species  # table with the original co-occurrence matrix
pairs<-read.table("Sign.BayesM.txt", header = TRUE, sep = "\t")   # three columns table with the significant species pairs and results of null model after Bayes correction
NullModel <- ifelse(pairs$Zscore < 0, "AGGR", "SEGR")      # column 1 ="Sp1", column 2 ="Sp2", 
pairs <- pairs[,c(1:2)]                                   # column 3 = "NullModel" containing null model result, either "SEGR" or "AGGR" for each pair
pairs <- cbind(pairs, NullModel)    

env<-as.data.frame(env)      # table for the environmental values for each site
coor<-as.data.frame(coords)    # table with coordinates for each site

#### PCA for the environmental variables
npcs <- 6# set the number of PCs you wish to retain
pc.env<-prcomp(env, center = FALSE, scale. = FALSE) 
pc.loadings <- pc.env$rotation
write.csv(pc.loadings,"./pca_loadings.csv")
summary(pc.env) # show variance explained by each component
plot(pc.env)   
comp.env<-pc.env$x[,1:npcs] 

#### preparation of the columns for store the results from the tests
pairs$coor.test<-NA
pairs$env.test<-NA


#### For each species pairs identification of the four co-occurrence classes in each site ("co00"= both absent, "co11"=both present
#### "co01" and "co10"= checkerboard distributions) and spatial configuration and environmental tests

for (i in 1:nrow(pairs)) {
  
  sp1<-as.character(pairs[i,1])
  sp2<-as.character(pairs[i,2])
  
  tab<-presabs[,c(sp1,sp2)]                
  x<-pairs[pairs$Sp1==sp1&pairs$Sp2==sp2,] 
  
  tab$id<-apply(tab,1,sum)     # identification of the four classes
  tab$id[tab$id == 2] <- "co11"
  tab$id[tab$id == 0] <- "co00"
  tab$id[tab[[1]]==1&tab[[2]]==0] <- "co10"
  tab$id[tab[[1]]==0&tab[[2]]==1] <- "co01"
  
  tab<-cbind(presabs[1],tab)  # association of the site id to the smaller table  
  
  if (x$NullModel=="SEGR") {                             # identification of the pattern of segregation or
    tab1<-tab[tab[[4]]=="co10"|tab[[4]]=="co01",]          # aggregation for the considered species pairs
  }else if (x$NullModel=="AGGR") {                       # to select the sites to test in the following lines
    tab1<-tab[tab[[4]]=="co11"|tab[[4]]=="co00",]
  }
  
  
  #### Test of spatial arrangement: MANOVA on the coordinates
  coor.tab1<-merge(coor, tab1, by="row.names",all = FALSE)
  fac1=factor(tab1$id)
  xx=as.matrix(coor.tab1[,(2:ncol(coor)+1)])
  m.coor=manova(xx~fac1)
  sum.coor=summary(m.coor)
  pairs$coor.test[i]<-round(sum.coor$stats[1,6],digits = 5) 
  
  #### Test of environmental filter: MANOVA on the PCA axes
  env.tab1<-merge(comp.env,tab1, by="row.names",all = FALSE) 
  fac1=factor(env.tab1$id)
  yy<-as.matrix(env.tab1[,2:(ncol(comp.env)+1)])
  m.env<-manova(yy~fac1)
  sum.env=summary(m.env)
  pairs$env.test[i]<-round(sum.env$stats[1,6],digits = 5)   
  
} 

pairs
### OUTOCOME: The pairs table will store p-values from the spatial test on coordinates ("coor.test" column) 

############################    END    #############################################################