library(tidyverse)
library(here)
library(igraph)
library(RSiena)
library(parallel)
library(rstudioapi)
library(lattice)
library(sna)
library(ggplot2)

# Set working directory to current script
if (Sys.getenv('RSTUDIO') == 1) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory(function(){})[1])
}

# load data for participant nodes

patientdata <-read_csv(paste0(here(),'/data/rsif20210850_si_008.csv'), show_col_types = FALSE)
healthycontrol <-read_csv(paste0(here(),'/data/healthycontrol.csv'), show_col_types = FALSE)

stroke1 <-read_csv(paste0(here(),'/data/stroke1.csv'), show_col_types = FALSE)
stroke2 <-read_csv(paste0(here(),'/data/stroke2.csv'), show_col_types = FALSE)
stroke3 <- read_csv(paste0(here(),'/data/stroke3.csv'), show_col_types = FALSE)

# load edge data
week2 <- read_csv(paste0(here(),'/data/week2.csv'))
month3 <- read_csv(paste0(here(),'/data/month3.csv'))
year1 <- read_csv(paste0(here(),'/data/year1.csv'))

healthweek2 <- read_csv(paste0(here(),'/data/healthweek2.csv'))
healthmonth3<- read_csv(paste0(here(),'/data/healthmonth3.csv'))
healthyear1<- read_csv(paste0(here(),'/data/healthyear1.csv'))

#SIENA 
#Convert to dataframes
stroke_graph1 <- graph.data.frame(week2)
stroke_graph2 <- graph.data.frame(month3)
stroke_graph3 <- graph.data.frame(year1)

health_graph1<- graph.data.frame(healthweek2)
health_graph2<- graph.data.frame(healthmonth3)
health_graph3<- graph.data.frame(healthyear1)

#convert data frames to Adjacency Matrix 
matrix1 <- get.adjacency(stroke_graph1, sparse = FALSE)
matrix2 <- get.adjacency(stroke_graph2, sparse = FALSE)
matrix3 <- get.adjacency(stroke_graph3, sparse = FALSE)

health_1<- get.adjacency(health_graph1, sparse = FALSE)
health_2<- get.adjacency(health_graph2, sparse = FALSE)
health_3<- get.adjacency(health_graph3, sparse = FALSE)


#3D Matrices stacked to form array of stroke data
stroke_dv <- sienaDependent(array(c(matrix1,matrix2,matrix3), dim = c(49,49,3)))
class(stroke_dv)

health_dv <- sienaDependent(array(c(health_1,health_2,health_3), dim = c(49,49,3)))

# prepare gender as a static covariate
gender_cov <- tibble(ID = row.names(matrix2)) %>%
  left_join(patientdata %>% select(ID, Gender)) %>%
  select(-ID) %>%
  mutate(Gender = factor(Gender))
gender_cov <- coCovar(as.integer(gender_cov$Gender))

# prepare handedness as a static covariate
handedness_cov <- tibble(ID = row.names(matrix2)) %>%
  left_join(patientdata %>% select(ID,Handedness)) %>%
  select(-ID) %>%
  mutate(Handedness = factor(Handedness))
handedness_cov <- coCovar(as.integer(handedness_cov$Handedness))

# prepare lesionside as a static covariate
lesionside_cov <- tibble(ID = row.names(matrix2)) %>%
  left_join(patientdata %>% select(ID,Lesionside)) %>%
  select(-ID) %>%
  mutate(Lesionside = factor(Lesionside))
lesionside_cov <- coCovar(as.integer(lesionside_cov$Lesionside))

#Siena modelling
model_data <- sienaDataCreate(stroke_dv, health_dv, gender_cov, handedness_cov, lesionside_cov)
print01Report(model_data, modelname = 'Siena stroke')

#Effects modelling
model_effects <- getEffects(model_data)

#Check the number of default effects provided by Siena
model_effects # 
effectsDocumentation(model_effects)

#Test 1
model_effects <- includeEffects(model_effects, altX, interaction1 = "gender_cov")

#Test 2 
model_effects <- includeEffects(model_effects, simX, interaction1 = "handedness_cov")

# Test 3
model_effects <- includeEffects(model_effects, egoX, altX, egoXaltX, interaction1 = "lesionside_cov") #8

# Test 4
model_effects <- includeEffects(model_effects, egoX, altX, egoXaltX, interaction1 = "gender_cov")

#Full specification
model_effects 

#Algorithm 'stroke1'
algorithm <- sienaAlgorithmCreate(projname = 'stroke1', seed=2225)

#number of cores
num_cores = detectCores() - 1 #9 cores in my laptop
 
#Siena Simulation
model <- siena07(algorithm,data = model_data, effects = model_effects, useCluster = TRUE, nbrNodes = num_cores, returnDeps=T)
model
myRes <- model #store the results in myRes

# Time heterogeneity
tt <- sienaTimeTest(myRes)
summary(tt)

#plot ego-alter heat map
obj_n <- function(vi, vj){
  b1*(vi-v_av) + b2*(vj-v_av) + b3*(1 - abs(vi-vj)/ran_v - sim_av)
}
# Output
v_av <- 11  # average 
sim_av <- 2.44# average similarity
ran_v <-  3.44  # range of values
b1 <- 4.744 # ego
b2 <- 6.211 #  alter
b3 <- 5.74  # similarity
vv <- c(5,10,15,20)  # Define the values of v for which the table is to be given
sel_tab <- outer(vv, vv, obj_n)  # calculate the table
round(sel_tab,3)  # display the table: rows are egos, columns are alters

# Graphics
graphics.off()
levelplot(sel_tab,col.regions=colorRampPalette(c("white","red", "orange", "yellow")), 
          scales=list(col='white'), main = "Heatmap",
          xlab = "Ego's Attribute Value( Stroke)", ylab = "Alter's Attribute Value (Healthy controls)" )

filled.contour(x=vv, y=vv, z=sel_tab, zlim=c(min(sel_tab),max(sel_tab)), nlevels=50, axes=TRUE, 
               color.palette=colorRampPalette(c("white","yellow", "green", "blue")),
               main = "heatmap",
               xlab = "Ego Attribute Value (Stroke)", ylab = "Alter Attribute Value (Healthy controls)" )



#-----------------------------SIENA Analysis-----------------------------------------------#--------------------------------------------------#

# create the array using the first two waves
stroke2 <- sienaDependent(stroke_dv[,,c(1,3)])
health2 <- sienaDependent(health_dv[,,c(1,3)])

# Siena data
myDataSim <- sienaDataCreate(stroke2, health2)
myDataSim

# Default effects
rm(simEffects)
simEffects <- getEffects(myDataSim) 

# Algorithm
modelOptionsSim1 <- sienaAlgorithmCreate(projname='Stroke simulations',
                                         nsub=5, seed=78684)
# Observed Siena simulations
myResultsObs <- siena07(modelOptionsSim1, data=myDataSim,
                        effects=simEffects, returnDeps=T) 
myResultsObs 

#Get the mean values
simEffects$initialValue[simEffects$include==T] <- myResultsObs$theta
simEffects$fix[simEffects$include==T] <- TRUE
simEffects

#High influence
simEffectsHiPI <- simEffects
simEffectsHiPI <- setEffect(simEffectsHiPI, 	
                            avDegIntn, interaction1='stroke2', name='health2',
                            initialValue=12, fix=T)  
simEffectsHiPI

#Low influence
simEffectsLoPI <- simEffects
simEffectsLoPI <- setEffect(simEffectsLoPI, avDegIntn, interaction1='stroke2', name='health2',
                            initialValue=.5, fix=T)  
simEffectsLoPI

# Simulate only for the first 50 iterations
nIter <- 50
modelOptionsSim2 <- sienaAlgorithmCreate(projname='siena brainstroke',
                                         nsub=0,     # no phase 2 iterations (because parameter values are fixed)
                                         n3=nIter,   # number of iterations (i.e., number of independent runs)
                                         simOnly=T,  # skips calculation of covariance matrix
                                         seed=4)     
# Simulate for the high 
myResultsSimHiPI <- siena07(modelOptionsSim2, data=myDataSim,
                            effects=simEffectsHiPI, returnDeps=T, returnChains=T)

# Simulate for the low 
myResultsSimLoPI <- siena07(modelOptionsSim2, data=myDataSim,
                            effects=simEffectsLoPI, returnDeps=T, returnChains=T)

#   Create vectors to store means
meanSimO <- rep(0,1000)       
meanSimHiPI <- rep(0,nIter)   
meanSimLoPI <- rep(0,nIter)  

for (m in 1:1000) {
  meanSimO[m] <- mean( myResultsObs$sims[[m]][[1]]$health[[1]] )
}


#   Extract from simulations 
for (m in 1:nIter) {
  meanSimHiPI[m] <- mean( myResultsSimHiPI$sims[[m]][[1]]$health2[[1]] )
  meanSimLoPI[m] <- mean( myResultsSimLoPI$sims[[m]][[1]]$health2[[1]] )
}
meanSimHiPI
meanSimLoPI


# Mean of 50th simulation run 
mean( myResultsSimHiPI$sims[[50]][[1]]$stroke2[[1]] ) 

dim(stroke2)
colMeans(stroke2, na.rm=T)  
meanObs <- mean(stroke2[49, 49, 2])  
dim(stroke2)

# Box Plot
set.seed(9000)
toPlot <- rbind(
  data.frame(cond="Observed",mean=meanSimO),
  data.frame(cond="High",mean=meanSimHiPI),
  data.frame(cond="Low",mean=meanSimLoPI) )
boxplot(mean ~ cond, data=toPlot, main="Mean values vs Network Influence", xlab="Network Influence",
        ylab="Mean values")
abline(h=meanObs)

# plot mean  over time
simName <- myResultsSimHiPI    # look at hi PI simulations first
chainNum <- 50
datChain <- t(matrix(unlist(simName$chain[[chainNum]][[1]][[1]]), 
                     nc=length(simName$chain[[chainNum]][[1]][[1]]))) 
datBehavChain <- datChain[datChain[,2]=="1",]    # condense to behavior change opportunities
( nBehavChanges <- dim(datBehavChain)[1] )       # number of behavior change opportunities
meanSeqB <- nBehavChanges + 1    # only record behavior changes
meanSeqB[1] <- meanObs   #mean obtained
for (i in 1:nBehavChanges) {
  meanSeqB[i+1] <- meanSeqB[i] + (as.numeric(datBehavChain[i,6])/50)
}

( nChanges <- dim(datChain)[1] )
meanSeqF <- nChanges + 1         
meanSeqF[1] <- meanObs   
for (i in 1:nChanges) {
  if (datChain[i,2]=="0") { 
    meanSeqF[i+1] <- meanSeqF[i]
  }
  if (datChain[i,2]=="1") { 
    meanSeqF[i+1] <- meanSeqF[i] + (as.numeric(datChain[i,6])/50)
  }
}

simChanges <- rep(0, nIter)
for (i in 1: nIter) {simChanges[i] <- length(simName$chain[[i]][[1]][[1]]) }
maxChanges <- max(simChanges)

#Extract mean happiness values for all 50 chains, with network change opportunities 
seqs <- matrix(meanSimHiPI, nr= nIter, nc=maxChanges+1)  # fill matrix with final mean
seqs[,1] <- meanObs
for (i in 1: nIter) {
  datChain <- t(matrix(unlist(simName$chain[[i]][[1]][[1]]), nc=length(simName$chain[[i]][[1]][[1]]))) 
  for (j in 1: simChanges[i]) {
    if (datChain[j,2]=="0") { 
      seqs[i,j+1] <- seqs[i,j]
    }
    if (datChain[j,2]=="1") { 
      seqs[i,j+1] <- seqs[i,j] + (as.numeric(datChain[j,6])/50)
    }
  }
}

# plot the "loess" curve
micros <- 1:dim(seqs)[2]
lo1 <- loess(seqs[1,] ~ micros)
l1 <- predict(lo1, micros)

plot(x=1:dim(seqs)[2], y=seqs[1,], ylim=c(min(seqs),max(seqs)),type="l", 
     ylab='stroke mean values', xlab='micro step', col='white',
     main='mean stroke value change over time')	
for (i in 1:25) {
  lines(x=micros, y=predict(loess(seqs[i,] ~ micros), micros), col=colors()[i*10], lwd=3.0)
}





