# simulate representative datasets for posterior generation
# the expectations are just calculations
# this just simulates trees under the expected values as starting points for
# creating a meaningless alignment to load into BEAUti.  BEAUti requires a NEXUS file
# with a fixed number of sequences, even though I will ultimately just sample from the
# prior to generate trees (i.e. sample without considering the data; the data is listed
# as 'NNN' for each individual in the data block of the XML file.

# set working directory
setwd("/your/path/here")

require(ape)
require(geiger)
require(TreeSim)

# exponential expectation
(log(25)-log(2))/5 # 0.5051457
2*((log(25)-log(2))/5) # 1.010291

(log(100)-(log(2)))/5 # 0.7824046
2*((log(100)-(log(2)))/5) # 1.564809

(log(500)-(log(2)))/5 # 1.104292
2*((log(500)-(log(2)))/5) # 2.208584

tree_25 <- sim.bd.taxa.age(n=25, numbsim=1, age=5, lambda=0.5051457, mu=0, frac=1.0, mrca=FALSE)
tree_100 <- sim.bd.taxa.age(n=100, numbsim=1, age=5, lambda=0.7824046, mu=0, frac=1.0, mrca=FALSE)
tree_500 <- sim.bd.taxa.age(n=500, numbsim=1, age=5, lambda=1.104, mu=0, frac=1.0, mrca=FALSE)

write.tree(phy=tree_25[[1]], file="25.tree.phy")
write.tree(phy=tree_100[[1]], file="100.tree.phy")
write.tree(phy=tree_500[[1]], file="500.tree.phy")