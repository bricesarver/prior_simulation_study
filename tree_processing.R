##tree_processing.R
##processes the results from a folder of results for the simulation study

rm(list=ls())
##call libraries
library(ape)
library(geiger)
library(parallel)
library(laser)

##functions
rescale_posterior <- function(posterior_tree, iterator){
  posterior_tree$edge.length <- ((posterior_tree$edge.length*original_roots[iterator])/max(branching.times(posterior_tree)))
  posterior_tree
}

birthdeath_calc <- function(x){
  return(bd(branching.times(x))$r1)
}

#specify the first working directory here
#example: working <- "./100.bd.strict.subst.trees.directory"

handle <- strsplit(working, split="\\./")[[1]][2]
stem <- strsplit(handle, split="\\.subst")[[1]][1]
##set working directory for trial
setwd(working)

##initialize lists and get filename
##original_posterior <- read.nexus("../100.bd.strict.(time).trees")
original_posterior <- read.nexus(paste("../", stem, ".(time).trees", sep=""))
##sampled_trees <- read.tree("100.bd.strict.subst.trees.sampled.trees")
sampled_trees <- read.tree(paste(stem, ".subst.trees.sampled.trees", sep=""))
original_trees <- as.list(1:10)
original_roots <- 1:10

##process the simulated trees i.e. get the tree name and pull the ultrametric tree from the simulated posterior
##also, return a vector of roots for future processing
for (i in 1:length(names(sampled_trees))){
  print(paste("Simulated tree: ",i, sep=""))
  holder_name <- names(sampled_trees)[[i]]
  print(holder_name)
  original_trees[[i]] <- original_posterior[[holder_name]]
  print(original_trees)
  original_roots[i] <- max(branching.times(original_trees[[i]]))
  
}

#initialize lambda and r lists for the simulated trees
#fit birth-death and Yule models to these trees
original_lambda <- as.list(1:10)
original_r <- as.list(1:10)

for (i in 1:10){
  original_lambda[i] <- pureBirth(branching.times(original_trees[[i]]))$r1
  original_r[i] <- bd(branching.times(original_trees[[i]]))$r1
}

#some descriptive plots of the original trees
#plot(x=seq.int(from=1, to=10), y=original_lambda, col="red")
#points(x=seq.int(from=1, to=10), y=original_r, col="blue")
#plot(x=seq.int(1, 10), y=(unlist(original_lambda) - unlist(original_r)))

#begin processing the results
#initialize the lists that will store results across replicates

bd_strict_rs <- list()
length(bd_strict_rs) <- 10
bd_strict_bd <- list()
length(bd_strict_bd) <- 10
bd_strict_yule <- list()
length(bd_strict_yule) <- 10

bd_ucln_rs <- list()
length(bd_ucln_rs) <- 10
bd_ucln_bd <- list()
length(bd_ucln_bd) <- 10
bd_ucln_yule <- list()
length(bd_ucln_yule) <- 10

yule_strict_rs <- list()
length(yule_strict_rs) <- 10
yule_strict_bd <- list()
length(yule_strict_bd) <- 10
yule_strict_yule <- list()
length(yule_strict_yule) <- 10

yule_ucln_rs <- list()
length(yule_ucln_rs) <- 10
yule_ucln_bd <- list()
length(yule_ucln_bd) <- 10
yule_ucln_yule <- list()
length(yule_ucln_yule) <- 10


#Principal function calls and calculation.  These are self-contained loops per evaluation scheme.
#Kept modular on purpose for downstream modification if need be.
#mc.cores should be changed to something more appropriate for your system if needed

print("Processing Birth-Death-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_strict <- read.nexus(posterior_names[1])[1001:10001]
  bd_strict_rs[[i]] <- mclapply(bd_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_strict_rs[[i]])){
    #print(j)
    bd_strict_yule[[i]][j] <- pureBirth(branching.times(bd_strict_rs[[i]][[j]]))$r1
    
  }
  bd_strict_bd[[i]] <- mclapply(bd_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Birth-Death-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_ucln <- read.nexus(posterior_names[2])[1001:10001]
  bd_ucln_rs[[i]] <- mclapply(bd_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_ucln_rs[[i]])){
    #print(j)
    bd_ucln_yule[[i]][j] <- pureBirth(branching.times(bd_ucln_rs[[i]][[j]]))$r1
    
  }
  bd_ucln_bd[[i]] <- mclapply(bd_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_strict <- read.nexus(posterior_names[3])[1001:10001]
  yule_strict_rs[[i]] <- mclapply(yule_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for (j in 1:length(yule_strict_rs[[i]])){
    #print(j)
    yule_strict_yule[[i]][j] <- pureBirth(branching.times(yule_strict_rs[[i]][[j]]))$r1
    
  }
  yule_strict_bd[[i]] <- mclapply(yule_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_ucln <- read.nexus(posterior_names[4])[1001:10001]
  yule_ucln_rs[[i]] <- mclapply(yule_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(yule_ucln_rs[[i]])){
    #print(j)
    yule_ucln_yule[[i]][j] <- pureBirth(branching.times(yule_ucln_rs[[i]][[j]]))$r1
    
  }
  yule_ucln_bd[[i]] <- mclapply(yule_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Execution complete.")

#plotting and summary code

mean_yulestrictyule <- list()
mean_yuleuclnyule <- list()
mean_bduclnyule <- list()
mean_bdstrictyule <- list()

for (i in 1:length(yule_strict_yule)){
  mean_yulestrictyule[[i]] <- mean(yule_strict_yule[[i]])
}

for (i in 1:length(yule_ucln_yule)){
  mean_yuleuclnyule[[i]] <- mean(yule_ucln_yule[[i]])
}

for (i in 1:length(bd_ucln_yule)){
  mean_bduclnyule[[i]] <- mean(bd_ucln_yule[[i]])
}

for (i in 1:length(bd_strict_yule)){
  mean_bdstrictyule[[i]] <- mean(bd_strict_yule[[i]])
}

dat <- cbind(unlist(mean_yulestrictyule), unlist(mean_yuleuclnyule), unlist(mean_bdstrictyule), unlist(mean_bduclnyule), unlist(original_lambda))
dat <- as.data.frame(dat)
colnames(dat) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat, main="100:BD:Strict", ylab="Lambda", ylim=c(0, 2))

mean_yulestrictbd <- list()
mean_yuleuclnbd <- list()
mean_bduclnbd <- list()
mean_bdstrictbd <- list()

for (i in 1:length(yule_strict_bd)){
  mean_yulestrictbd[[i]] <- mean(as.numeric(yule_strict_bd[[i]]))
}
for (i in 1:length(yule_ucln_bd)){
  mean_yuleuclnbd[[i]] <- mean(as.numeric(yule_ucln_bd[[i]]))
}
for (i in 1:length(bd_ucln_bd)){
  mean_bduclnbd[[i]] <- mean(as.numeric(bd_ucln_bd[[i]]))
}
for (i in 1:length(bd_strict_bd)){
  mean_bdstrictbd[[i]] <- mean(as.numeric(bd_strict_bd[[i]]))
}

dat2<-cbind(unlist(mean_yulestrictbd), unlist(mean_yuleuclnbd), unlist(mean_bdstrictbd), unlist(mean_bduclnbd), unlist(original_r))
dat2 <- as.data.frame(dat2)
colnames(dat2) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat2, main="100:BD:Strict", ylab="Net Diversification Rate", ylim=c(0, 2))

all_means <- list(dat)
all_means <- append(all_means, list(dat2))
names(all_means) <- c(paste(stem, ".l", sep=""), paste(stem, ".r", sep=""))

setwd("..")

##NEXT FOLDER

#set second working directory:
#example: working <- "./100.bd.ucln.subst.trees.directory"

handle <- strsplit(working, split="\\./")[[1]][2]
stem <- strsplit(handle, split="\\.subst")[[1]][1]
##set working directory for trial
setwd(working)

##initialize lists and get filename
##original_posterior <- read.nexus("../100.bd.strict.(time).trees")
original_posterior <- read.nexus(paste("../", stem, ".(time).trees", sep=""))
##sampled_trees <- read.tree("100.bd.strict.subst.trees.sampled.trees")
sampled_trees <- read.tree(paste(stem, ".subst.trees.sampled.trees", sep=""))
original_trees <- as.list(1:10)
original_roots <- 1:10

##process the simulated trees i.e. get the tree name and pull the ultrametric tree from the simulated posterior
##also, return a vector of roots for future processing
for (i in 1:length(names(sampled_trees))){
  print(paste("Simulated tree: ",i, sep=""))
  holder_name <- names(sampled_trees)[[i]]
  print(holder_name)
  original_trees[[i]] <- original_posterior[[holder_name]]
  print(original_trees)
  original_roots[i] <- max(branching.times(original_trees[[i]]))
  
}

#initialize lambda and r lists for the simulated trees
#fit birth-death and Yule models to these trees
original_lambda <- as.list(1:10)
original_r <- as.list(1:10)

for (i in 1:10){
  original_lambda[i] <- pureBirth(branching.times(original_trees[[i]]))$r1
  original_r[i] <- bd(branching.times(original_trees[[i]]))$r1
}

#some descriptive plots of the original trees
plot(x=seq.int(from=1, to=10), y=original_lambda, col="red")
points(x=seq.int(from=1, to=10), y=original_r, col="blue")
plot(x=seq.int(1, 10), y=(unlist(original_lambda) - unlist(original_r)))

#begin processing the results
#initialize the lists that will store results across replicates

bd_strict_rs <- list()
length(bd_strict_rs) <- 10
bd_strict_bd <- list()
length(bd_strict_bd) <- 10
bd_strict_yule <- list()
length(bd_strict_yule) <- 10

bd_ucln_rs <- list()
length(bd_ucln_rs) <- 10
bd_ucln_bd <- list()
length(bd_ucln_bd) <- 10
bd_ucln_yule <- list()
length(bd_ucln_yule) <- 10

yule_strict_rs <- list()
length(yule_strict_rs) <- 10
yule_strict_bd <- list()
length(yule_strict_bd) <- 10
yule_strict_yule <- list()
length(yule_strict_yule) <- 10

yule_ucln_rs <- list()
length(yule_ucln_rs) <- 10
yule_ucln_bd <- list()
length(yule_ucln_bd) <- 10
yule_ucln_yule <- list()
length(yule_ucln_yule) <- 10


#Principal function calls and calculation.  These are self-contained loops per evaluation scheme.
#Kept modular on purpose for downstream modification if need be.

print("Processing Birth-Death-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_strict <- read.nexus(posterior_names[1])[1001:10001]
  bd_strict_rs[[i]] <- mclapply(bd_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_strict_rs[[i]])){
    #print(j)
    bd_strict_yule[[i]][j] <- pureBirth(branching.times(bd_strict_rs[[i]][[j]]))$r1
    
  }
  bd_strict_bd[[i]] <- mclapply(bd_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Birth-Death-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_ucln <- read.nexus(posterior_names[2])[1001:10001]
  bd_ucln_rs[[i]] <- mclapply(bd_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_ucln_rs[[i]])){
    #print(j)
    bd_ucln_yule[[i]][j] <- pureBirth(branching.times(bd_ucln_rs[[i]][[j]]))$r1
    
  }
  bd_ucln_bd[[i]] <- mclapply(bd_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_strict <- read.nexus(posterior_names[3])[1001:10001]
  yule_strict_rs[[i]] <- mclapply(yule_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for (j in 1:length(yule_strict_rs[[i]])){
    #print(j)
    yule_strict_yule[[i]][j] <- pureBirth(branching.times(yule_strict_rs[[i]][[j]]))$r1
    
  }
  yule_strict_bd[[i]] <- mclapply(yule_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_ucln <- read.nexus(posterior_names[4])[1001:10001]
  yule_ucln_rs[[i]] <- mclapply(yule_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(yule_ucln_rs[[i]])){
    #print(j)
    yule_ucln_yule[[i]][j] <- pureBirth(branching.times(yule_ucln_rs[[i]][[j]]))$r1
    
  }
  yule_ucln_bd[[i]] <- mclapply(yule_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Execution complete.")

#plotting and summary code

mean_yulestrictyule <- list()
mean_yuleuclnyule <- list()
mean_bduclnyule <- list()
mean_bdstrictyule <- list()

for (i in 1:length(yule_strict_yule)){
  mean_yulestrictyule[[i]] <- mean(yule_strict_yule[[i]])
}

for (i in 1:length(yule_ucln_yule)){
  mean_yuleuclnyule[[i]] <- mean(yule_ucln_yule[[i]])
}

for (i in 1:length(bd_ucln_yule)){
  mean_bduclnyule[[i]] <- mean(bd_ucln_yule[[i]])
}

for (i in 1:length(bd_strict_yule)){
  mean_bdstrictyule[[i]] <- mean(bd_strict_yule[[i]])
}

dat <- cbind(unlist(mean_yulestrictyule), unlist(mean_yuleuclnyule), unlist(mean_bdstrictyule), unlist(mean_bduclnyule), unlist(original_lambda))
dat <- as.data.frame(dat)
colnames(dat) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat, main="100:BD:UCLN", ylab="Lambda", ylim=c(0, 2))

mean_yulestrictbd <- list()
mean_yuleuclnbd <- list()
mean_bduclnbd <- list()
mean_bdstrictbd <- list()

for (i in 1:length(yule_strict_bd)){
  mean_yulestrictbd[[i]] <- mean(as.numeric(yule_strict_bd[[i]]))
}
for (i in 1:length(yule_ucln_bd)){
  mean_yuleuclnbd[[i]] <- mean(as.numeric(yule_ucln_bd[[i]]))
}
for (i in 1:length(bd_ucln_bd)){
  mean_bduclnbd[[i]] <- mean(as.numeric(bd_ucln_bd[[i]]))
}
for (i in 1:length(bd_strict_bd)){
  mean_bdstrictbd[[i]] <- mean(as.numeric(bd_strict_bd[[i]]))
}

dat2<-cbind(unlist(mean_yulestrictbd), unlist(mean_yuleuclnbd), unlist(mean_bdstrictbd), unlist(mean_bduclnbd), unlist(original_r))
dat2 <- as.data.frame(dat2)
colnames(dat2) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat2, main="100:BD:UCLN", ylab="Net Diversification Rate", ylim=c(0, 2))

all_means <- append(all_means, list(dat))
all_means <- append(all_means, list(dat2))
names(all_means)[3] <- paste(stem, ".l", sep="")
names(all_means)[4] <- paste(stem, ".r", sep="")

setwd("..")

#next directory
#example: working <- "./25.bd.strict.subst.trees.directory"
handle <- strsplit(working, split="\\./")[[1]][2]
stem <- strsplit(handle, split="\\.subst")[[1]][1]
##set working directory for trial
setwd(working)

##initialize lists and get filename
##original_posterior <- read.nexus("../100.bd.strict.(time).trees")
original_posterior <- read.nexus(paste("../", stem, ".(time).trees", sep=""))
##sampled_trees <- read.tree("100.bd.strict.subst.trees.sampled.trees")
sampled_trees <- read.tree(paste(stem, ".subst.trees.sampled.trees", sep=""))
original_trees <- as.list(1:10)
original_roots <- 1:10

##process the simulated trees i.e. get the tree name and pull the ultrametric tree from the simulated posterior
##also, return a vector of roots for future processing
for (i in 1:length(names(sampled_trees))){
  print(paste("Simulated tree: ",i, sep=""))
  holder_name <- names(sampled_trees)[[i]]
  print(holder_name)
  original_trees[[i]] <- original_posterior[[holder_name]]
  print(original_trees)
  original_roots[i] <- max(branching.times(original_trees[[i]]))
  
}

#initialize lambda and r lists for the simulated trees
#fit birth-death and Yule models to these trees
original_lambda <- as.list(1:10)
original_r <- as.list(1:10)

for (i in 1:10){
  original_lambda[i] <- pureBirth(branching.times(original_trees[[i]]))$r1
  original_r[i] <- bd(branching.times(original_trees[[i]]))$r1
}

#some descriptive plots of the original trees
plot(x=seq.int(from=1, to=10), y=original_lambda, col="red")
points(x=seq.int(from=1, to=10), y=original_r, col="blue")
plot(x=seq.int(1, 10), y=(unlist(original_lambda) - unlist(original_r)))

#begin processing the results
#initialize the lists that will store results across replicates

bd_strict_rs <- list()
length(bd_strict_rs) <- 10
bd_strict_bd <- list()
length(bd_strict_bd) <- 10
bd_strict_yule <- list()
length(bd_strict_yule) <- 10

bd_ucln_rs <- list()
length(bd_ucln_rs) <- 10
bd_ucln_bd <- list()
length(bd_ucln_bd) <- 10
bd_ucln_yule <- list()
length(bd_ucln_yule) <- 10

yule_strict_rs <- list()
length(yule_strict_rs) <- 10
yule_strict_bd <- list()
length(yule_strict_bd) <- 10
yule_strict_yule <- list()
length(yule_strict_yule) <- 10

yule_ucln_rs <- list()
length(yule_ucln_rs) <- 10
yule_ucln_bd <- list()
length(yule_ucln_bd) <- 10
yule_ucln_yule <- list()
length(yule_ucln_yule) <- 10


#Principal function calls and calculation.  These are self-contained loops per evaluation scheme.
#Kept modular on purpose for downstream modification if need be.

print("Processing Birth-Death-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_strict <- read.nexus(posterior_names[1])[1001:10001]
  bd_strict_rs[[i]] <- mclapply(bd_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_strict_rs[[i]])){
    #print(j)
    bd_strict_yule[[i]][j] <- pureBirth(branching.times(bd_strict_rs[[i]][[j]]))$r1
    
  }
  bd_strict_bd[[i]] <- mclapply(bd_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Birth-Death-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_ucln <- read.nexus(posterior_names[2])[1001:10001]
  bd_ucln_rs[[i]] <- mclapply(bd_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_ucln_rs[[i]])){
    #print(j)
    bd_ucln_yule[[i]][j] <- pureBirth(branching.times(bd_ucln_rs[[i]][[j]]))$r1
    
  }
  bd_ucln_bd[[i]] <- mclapply(bd_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_strict <- read.nexus(posterior_names[3])[1001:10001]
  yule_strict_rs[[i]] <- mclapply(yule_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for (j in 1:length(yule_strict_rs[[i]])){
    #print(j)
    yule_strict_yule[[i]][j] <- pureBirth(branching.times(yule_strict_rs[[i]][[j]]))$r1
    
  }
  yule_strict_bd[[i]] <- mclapply(yule_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_ucln <- read.nexus(posterior_names[4])[1001:10001]
  yule_ucln_rs[[i]] <- mclapply(yule_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(yule_ucln_rs[[i]])){
    #print(j)
    yule_ucln_yule[[i]][j] <- pureBirth(branching.times(yule_ucln_rs[[i]][[j]]))$r1
    
  }
  yule_ucln_bd[[i]] <- mclapply(yule_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Execution complete.")

#plotting and summary code

mean_yulestrictyule <- list()
mean_yuleuclnyule <- list()
mean_bduclnyule <- list()
mean_bdstrictyule <- list()

for (i in 1:length(yule_strict_yule)){
  mean_yulestrictyule[[i]] <- mean(yule_strict_yule[[i]])
}

for (i in 1:length(yule_ucln_yule)){
  mean_yuleuclnyule[[i]] <- mean(yule_ucln_yule[[i]])
}

for (i in 1:length(bd_ucln_yule)){
  mean_bduclnyule[[i]] <- mean(bd_ucln_yule[[i]])
}

for (i in 1:length(bd_strict_yule)){
  mean_bdstrictyule[[i]] <- mean(bd_strict_yule[[i]])
}

dat <- cbind(unlist(mean_yulestrictyule), unlist(mean_yuleuclnyule), unlist(mean_bdstrictyule), unlist(mean_bduclnyule), unlist(original_lambda))
dat <- as.data.frame(dat)
colnames(dat) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat, main="25:BD:STRICT", ylab="Lambda", ylim=c(0, 2))

mean_yulestrictbd <- list()
mean_yuleuclnbd <- list()
mean_bduclnbd <- list()
mean_bdstrictbd <- list()

for (i in 1:length(yule_strict_bd)){
  mean_yulestrictbd[[i]] <- mean(as.numeric(yule_strict_bd[[i]]))
}
for (i in 1:length(yule_ucln_bd)){
  mean_yuleuclnbd[[i]] <- mean(as.numeric(yule_ucln_bd[[i]]))
}
for (i in 1:length(bd_ucln_bd)){
  mean_bduclnbd[[i]] <- mean(as.numeric(bd_ucln_bd[[i]]))
}
for (i in 1:length(bd_strict_bd)){
  mean_bdstrictbd[[i]] <- mean(as.numeric(bd_strict_bd[[i]]))
}

dat2<-cbind(unlist(mean_yulestrictbd), unlist(mean_yuleuclnbd), unlist(mean_bdstrictbd), unlist(mean_bduclnbd), unlist(original_r))
dat2 <- as.data.frame(dat2)
colnames(dat2) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat2, main="25:BD:STRICT", ylab="Net Diversification Rate", ylim=c(0, 2))

all_means <- append(all_means, list(dat))
all_means <- append(all_means, list(dat2))
names(all_means)[5] <- c(paste(stem, ".l", sep=""))
names(all_means)[6] <- c(paste(stem, ".r", sep=""))

setwd("..")

#Final Directory:
#example: working <- "./25.bd.ucln.subst.trees.directory"
handle <- strsplit(working, split="\\./")[[1]][2]
stem <- strsplit(handle, split="\\.subst")[[1]][1]
##set working directory for trial
setwd(working)

##initialize lists and get filename
##original_posterior <- read.nexus("../100.bd.strict.(time).trees")
original_posterior <- read.nexus(paste("../", stem, ".(time).trees", sep=""))
##sampled_trees <- read.tree("100.bd.strict.subst.trees.sampled.trees")
sampled_trees <- read.tree(paste(stem, ".subst.trees.sampled.trees", sep=""))
original_trees <- as.list(1:10)
original_roots <- 1:10

##process the simulated trees i.e. get the tree name and pull the ultrametric tree from the simulated posterior
##also, return a vector of roots for future processing
for (i in 1:length(names(sampled_trees))){
  print(paste("Simulated tree: ",i, sep=""))
  holder_name <- names(sampled_trees)[[i]]
  print(holder_name)
  original_trees[[i]] <- original_posterior[[holder_name]]
  print(original_trees)
  original_roots[i] <- max(branching.times(original_trees[[i]]))
  
}

#initialize lambda and r lists for the simulated trees
#fit birth-death and Yule models to these trees
original_lambda <- as.list(1:10)
original_r <- as.list(1:10)

for (i in 1:10){
  original_lambda[i] <- pureBirth(branching.times(original_trees[[i]]))$r1
  original_r[i] <- bd(branching.times(original_trees[[i]]))$r1
}

#some descriptive plots of the original trees
plot(x=seq.int(from=1, to=10), y=original_lambda, col="red")
points(x=seq.int(from=1, to=10), y=original_r, col="blue")
plot(x=seq.int(1, 10), y=(unlist(original_lambda) - unlist(original_r)))

#begin processing the results
#initialize the lists that will store results across replicates

bd_strict_rs <- list()
length(bd_strict_rs) <- 10
bd_strict_bd <- list()
length(bd_strict_bd) <- 10
bd_strict_yule <- list()
length(bd_strict_yule) <- 10

bd_ucln_rs <- list()
length(bd_ucln_rs) <- 10
bd_ucln_bd <- list()
length(bd_ucln_bd) <- 10
bd_ucln_yule <- list()
length(bd_ucln_yule) <- 10

yule_strict_rs <- list()
length(yule_strict_rs) <- 10
yule_strict_bd <- list()
length(yule_strict_bd) <- 10
yule_strict_yule <- list()
length(yule_strict_yule) <- 10

yule_ucln_rs <- list()
length(yule_ucln_rs) <- 10
yule_ucln_bd <- list()
length(yule_ucln_bd) <- 10
yule_ucln_yule <- list()
length(yule_ucln_yule) <- 10


#Principal function calls and calculation.  These are self-contained loops per evaluation scheme.
#Kept modular on purpose for downstream modification if need be.

print("Processing Birth-Death-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_strict <- read.nexus(posterior_names[1])[1001:10001]
  bd_strict_rs[[i]] <- mclapply(bd_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_strict_rs[[i]])){
    #print(j)
    bd_strict_yule[[i]][j] <- pureBirth(branching.times(bd_strict_rs[[i]][[j]]))$r1
    
  }
  bd_strict_bd[[i]] <- mclapply(bd_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Birth-Death-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  bd_ucln <- read.nexus(posterior_names[2])[1001:10001]
  bd_ucln_rs[[i]] <- mclapply(bd_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(bd_ucln_rs[[i]])){
    #print(j)
    bd_ucln_yule[[i]][j] <- pureBirth(branching.times(bd_ucln_rs[[i]][[j]]))$r1
    
  }
  bd_ucln_bd[[i]] <- mclapply(bd_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-Strict...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_strict <- read.nexus(posterior_names[3])[1001:10001]
  yule_strict_rs[[i]] <- mclapply(yule_strict, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for (j in 1:length(yule_strict_rs[[i]])){
    #print(j)
    yule_strict_yule[[i]][j] <- pureBirth(branching.times(yule_strict_rs[[i]][[j]]))$r1
    
  }
  yule_strict_bd[[i]] <- mclapply(yule_strict_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Processing Yule-UCLN...")
for (i in 1:10){
  posterior_names <- list.files(pattern=paste("*rep_", i, "_.+.time.trees", sep=""))
  print(paste("Processing replicate ", i, sep=""))
  yule_ucln <- read.nexus(posterior_names[4])[1001:10001]
  yule_ucln_rs[[i]] <- mclapply(yule_ucln, rescale_posterior, i, mc.cores=getOption("mc.cores", 40))
  
  for(j in 1:length(yule_ucln_rs[[i]])){
    #print(j)
    yule_ucln_yule[[i]][j] <- pureBirth(branching.times(yule_ucln_rs[[i]][[j]]))$r1
    
  }
  yule_ucln_bd[[i]] <- mclapply(yule_ucln_rs[[i]], birthdeath_calc, mc.cores=getOption("mc.cores", 40))
}

print("Execution complete.")

#plotting and summary code

mean_yulestrictyule <- list()
mean_yuleuclnyule <- list()
mean_bduclnyule <- list()
mean_bdstrictyule <- list()

for (i in 1:length(yule_strict_yule)){
  mean_yulestrictyule[[i]] <- mean(yule_strict_yule[[i]])
}

for (i in 1:length(yule_ucln_yule)){
  mean_yuleuclnyule[[i]] <- mean(yule_ucln_yule[[i]])
}

for (i in 1:length(bd_ucln_yule)){
  mean_bduclnyule[[i]] <- mean(bd_ucln_yule[[i]])
}

for (i in 1:length(bd_strict_yule)){
  mean_bdstrictyule[[i]] <- mean(bd_strict_yule[[i]])
}

dat <- cbind(unlist(mean_yulestrictyule), unlist(mean_yuleuclnyule), unlist(mean_bdstrictyule), unlist(mean_bduclnyule), unlist(original_lambda))
dat <- as.data.frame(dat)
colnames(dat) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat, main="25:BD:UCLN", ylab="Lambda", ylim=c(0, 2))

mean_yulestrictbd <- list()
mean_yuleuclnbd <- list()
mean_bduclnbd <- list()
mean_bdstrictbd <- list()

for (i in 1:length(yule_strict_bd)){
  mean_yulestrictbd[[i]] <- mean(as.numeric(yule_strict_bd[[i]]))
}
for (i in 1:length(yule_ucln_bd)){
  mean_yuleuclnbd[[i]] <- mean(as.numeric(yule_ucln_bd[[i]]))
}
for (i in 1:length(bd_ucln_bd)){
  mean_bduclnbd[[i]] <- mean(as.numeric(bd_ucln_bd[[i]]))
}
for (i in 1:length(bd_strict_bd)){
  mean_bdstrictbd[[i]] <- mean(as.numeric(bd_strict_bd[[i]]))
}

dat2<-cbind(unlist(mean_yulestrictbd), unlist(mean_yuleuclnbd), unlist(mean_bdstrictbd), unlist(mean_bduclnbd), unlist(original_r))
dat2 <- as.data.frame(dat2)
colnames(dat2) <- c("Yule:Strict", "Yule:UCLN", "BD:Strict", "BD:UCLN", "Original")
boxplot(dat2, main="25:BD:UCLN", ylab="Net Diversification Rate", ylim=c(0, 2))

all_means <- append(all_means, list(dat))
all_means <- append(all_means, list(dat2))
names(all_means)[7] <- c(paste(stem, ".l", sep=""))
names(all_means)[8] <- c(paste(stem, ".r", sep=""))

setwd("..")

pdf(file="BD.results.pdf", width=12, height=8.5)
par(mfrow=c(2,4), oma=c(2, 2, 2, 2))

boxplot(all_means[[1]], ylab="Lambda", main="100:BD:Strict", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[1]]$Original))
boxplot(all_means[[2]], ylab="Net Diversification Rate", main="100:BD:Strict", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[2]]$Original))
boxplot(all_means[[3]], ylab="Lambda", main="100:BD:UCLN", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[3]]$Original))
boxplot(all_means[[4]], ylab="Net Diversification Rate", main="100:BD:UCLN", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[4]]$Original))

boxplot(all_means[[5]], ylab="Lambda", main="25:BD:Strict", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[5]]$Original))
boxplot(all_means[[6]], ylab="Net Diversification Rate", main="25:BD:Strict", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[6]]$Original))
boxplot(all_means[[7]], ylab="Lambda", main="25:BD:UCLN", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[7]]$Original))
boxplot(all_means[[8]], ylab="Net Diversification Rate", main="25:BD:UCLN", ylim=c(0,1.8), las=2)
abline(h=median(all_means[[8]]$Original))
title("Birth-Death Simulations", outer=TRUE)
dev.off()