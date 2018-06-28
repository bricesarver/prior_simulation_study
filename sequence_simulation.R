#sequence_simulation.R
#uses BEASTifier and seq-gen to simulate sequences on a set of topologies
#10 trees are sampled at random from a 'simulated posterior distribution'
#of which 50% of the samples are removed as a burn-in
#these 10 trees are written and are used as the basis for the other analyses
#5000 bp are simulated on the trees

#this is a one-shot program for each set of parameters.  filenames and paths should be
#changed to reflect what needs to be processed.  though I probably should have scripted
#this up in a way that took a vector of strings and generated files based on the values,
#this acted as a check that everything proceeded as expected

rm(list=ls())
#set working directory to a parent folder
#call it PATH; alternatively, execute this script above subdirectories

setwd(PATH)

#library() will halt on failure; require() will only return a logical if it fails
library(ape)

#principal stem of filename to process

#specify folder in which to place files
#example: filename <- "100.bd.ucln.subst.trees"

system(command=paste("mkdir", paste(filename, ".directory", sep=""),  sep=" "))
system(command=paste("cp BEASTifier config ", paste(filename, ".directory", sep=""), sep=""))

#pre-process and read in trees
trees <- read.nexus(file=filename)

trees.minusburnin <- trees[5001:10001]
pruned.trees <- sample(x=trees.minusburnin, size=10, replace=FALSE)

rescaled.trees <- pruned.trees

for (i in 1:length(rescaled.trees)){
  edges <- rescaled.trees[[i]]$edge.length*0.01
  rescaled.trees[[i]]$edge.length <- edges
  
}

setwd(paste(filename, ".directory", sep=""))
write.tree(pruned.trees, file=paste(filename, ".sampled.trees", sep=""), tree.names=TRUE)
write.tree(rescaled.trees, file=paste(filename, ".rescaled.trees", sep=""), tree.names=TRUE)

#simulation parameters defined following Weisrock, Harmon, and Larson.
numSites <- "-l 5000"
freq <- ("-f 0.1978,0.2784,0.2403,0.1835")
alpha <- ("-a 2.3592")
ratematrix <- ("-r 1.6493,2.9172,0.3969,0.9164,8.4170,1")
mods <- paste("GTR+G")

for (i in 1:length(rescaled.trees)) {
  seed <- abs(round(runif(1, max=1000000)))
  write.tree(rescaled.trees[[i]], append=FALSE, file=paste(filename,"_rep_", i, ".phy", sep=""))
  settings <- paste(" -mGTR ", numSites, " ", ratematrix, " ", freq, " ", "-on", " -z", seed, " ", alpha, " -q", sep="")
  #filename <- paste("b_", lambda, "_d_", mu, "_a_", age, "_n_", ntaxa, "_sim_", "GTR+G", "_rep_", i, ".NEX", sep="")
  print(seed)
  print(settings)
  print(filename)
  system(paste("seq-gen", settings, "<", paste(filename, "_rep_", i, ".phy", sep=""), ">",  paste(filename, "_rep_", i, ".NEX", sep=""), sep=""))
  }
allNexus <- list.files(pattern="*.NEX")
write(allNexus, file="files.txt", append=TRUE)
write(mods, file="mods.txt", append=TRUE)

#Run BEASTifier
system(command="./BEASTifier -config config")

#Make PBS scripts for submission to a distributed cluster computing system
#based on a template beast_sub.pbs, which will differ depending on your
#method of submission - just a bash loop
#system("counter=32; for i in *.xml; do echo $i; cat beast_sub.pbs | sed -e s/NUMBER/\"$counter\"/ -e s/FILENAME/\"$i\"/ > submit$counter.pbs; let counter=counter+1; echo $counter; done")