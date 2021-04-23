

### Analyze non-random patterns (A.K.A. significant species associations) from presence/absence 
### community matrix 


##########################################################################################
#                                                                                        #
#    "Probabilistic FRAMEWORK"                                                           #
#                                                                                        #
# Null model analysis - Probabilistic Co-occurrence - Blois Framework                    #
#                                                                                        #
# Modified by Zach Zbinden from:                                                         #
#  D'Amen, Gotelli & Guisan: Disentangling biotic interactions, environmental filters,   #
#  and dispersal limitation as drivers of species co-occurrence                          #
##########################################################################################

### Packages required: 
install.packages("BBmisc")
install.packages("adespatial")
install.packages("usdm")
install.packages("cooccur")


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
names(env)[names(env) == "ORD_STRA"] <- "stream.order"
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
hydro_dist <- read.csv("fluvial_dist.csv", header = T, sep = ",", row.names = 1)

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



##############################################################
#                                                            #
#                                                            #
### 2.0: Calculate Probabilistic Species Co-occurrence       #
#                                                            #
#                                                            #
##############################################################

library(cooccur)

# This function runs the pairwise analysis and generates the object necessary for summary and visualization
species.t <-t(species)
CO1 <- cooccur(mat= species.t, type = "spp_site", thresh =TRUE, spp_names = TRUE, true_rand_classifier = 0.1, prob = "comb")

# Shows number of random, postive, and negative associations 
summary(CO1)
plot(CO1)


sig.pairs <- as.data.frame(print(CO1))
sig.pairs <- sig.pairs[,-c(1:7)]
NullModel <- ifelse(sig.pairs$p_gt < 0.05, "AGGR", "SEGR")      # column 1 ="Sp1", column 2 ="Sp2", 
pairs <- sig.pairs[,c(3:4)]                                   # column 3 = "NullModel" containing null model result, either "SEGR" or "AGGR" for each pair
pairs <- cbind(pairs, NullModel)  
colnames(pairs) <- c("Sp1", "Sp2", "NullModel")



###########################################################
###########################################################
#########        	                      			     ########	
#########  		        3.0:  Blois framework        ########		
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
env<-as.data.frame(env)      # table for the environmental values for each site
coor<-as.data.frame(coords)    # table with coordinates for each site

#### PCA for the environmental variables
npcs <- 6 # set the number of PCs you wish to retain
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