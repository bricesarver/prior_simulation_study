# ltt_plots.R
# make LTT plots for the simulation study

library(parallel)
library(ape)

# functions
rescale_posterior <- function(posterior_tree) {
	val <- strsplit(names(posterior_tree), split="\\.")[[1]][4]
	print(val)
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots[as.numeric(val)]/max(branching.times(posterior_tree))
	posterior_tree
}

# set working directory
# example: setwd("~/final_simulation_study/fixed_root_simulations/")
setwd("/your/path/here")
# we first need to produce MCC trees for each of the replicates in question
# specify the directories of interest sub to the main directory
directories <- c("100.bd.strict.subst.trees.directory", "100.bd.ucln.subst.trees.directory")
replicates <- seq(1, 10, 1)

setwd(directories[1])

# specify treeannPath, a path to the command line version of treeannotator
# example: treeannPath <- "~/source/BEASTv1.7.5/bin/treeannotator"
treeannPath <- "/path/to/treeannotator/"
directory <- strsplit(getwd(), split="/")[[1]][7]

annotate_em <- function(repnumber) {
	trees <- list.files(pattern=glob2rx(paste("*_rep_", repnumber, "_*time.trees", sep="")))
	stem <- strsplit(directory, split="\\.subst.")[[1]][1]
	system(paste(treeannPath, "-heights median -burnin 1000 -limit 0.5", trees[1], paste(stem, repnumber, "bd-strict.annotated.tre", sep=".")))
	system(paste(treeannPath, "-heights median -burnin 1000 -limit 0.5", trees[2], paste(stem, repnumber, "bd-ucln.annotated.tre", sep=".")))
	system(paste(treeannPath, "-heights median -burnin 1000 -limit 0.5", trees[3], paste(stem, repnumber, "yule-strict.annotated.tre", sep=".")))
	system(paste(treeannPath, "-heights median -burnin 1000 -limit 0.5", trees[4], paste(stem, repnumber, "yule-ucln.annotated.tre", sep=".")))	
}

mclapply(replicates, annotate_em, mc.cores=40)
system("mkdir annotated_trees")
system("mv *.annotated.tre annotated_trees")

setwd("..")
setwd(directories[2])
directory <- strsplit(getwd(), split="/")[[1]][7]

mclapply(replicates, annotate_em, mc.cores=40)
system("mkdir annotated_trees")
system("mv *.annotated.tre annotated_trees")

# and, the exact same as above for the Yule analyses
# example: setwd("~/final_simulation_study/fixed_root_simulations_yule/")
setwd("your/yule/path")
directories <- c("100.yule.strict.subst.trees.directory", "100.yule.ucln.subst.trees.directory")
replicates <- seq(1, 10, 1)

setwd(directories[1])

directory <- strsplit(getwd(), split="/")[[1]][7]

mclapply(replicates, annotate_em, mc.cores=40)
system("mkdir annotated_trees")
system("mv *.annotated.tre annotated_trees")

setwd("..")
setwd(directories[2])
directory <- strsplit(getwd(), split="/")[[1]][7]

mclapply(replicates, annotate_em, mc.cores=40)
system("mkdir annotated_trees")
system("mv *.annotated.tre annotated_trees")
setwd("/../..")

# specify the paths to the anntated trees (*path), plus the directories that house them (*origpath)
# example: bdstrictpath <- "~/fixed_root_simulations/100.bd.strict.subst.trees.directory/annotated_trees"
# example: bdstrictorigpath <- "~/fixed_root_simulations/"

bduclnpath <- "/your/path/here"
bduclnorigpath <- "/your/path/here"

yulestrictpath <- "/your/path/here"
yulestrictorigpath <- "/your/path/here"

yuleuclnpath <- "/your/path/here"
yuleuclnorigpath <- "/your/path/here"

original_posterior_bds <- read.nexus(paste0(bdstrictorigpath, "100.bd.strict.(time).trees"))
sampled_trees_bds <- read.tree(paste0(bdstrictorigpath, "100.bd.strict.subst.trees.directory/100.bd.strict.subst.trees.sampled.trees"))
original_trees_bds <- as.list(1:10)
original_roots_bds <- 1:10

for (i in 1:length(names(sampled_trees_bds))) {
  print(paste0("Simulated tree: ",i))
  holder_name <- names(sampled_trees_bds)[[i]]
  print(holder_name)
  original_trees_bds[[i]] <- original_posterior_bds[[holder_name]]
  print(original_trees_bds)
  original_roots_bds[i] <- max(branching.times(original_trees_bds[[i]]))
}

original_posterior_bdu <- read.nexus(paste0(bduclnorigpath, "100.bd.ucln.(time).trees"))
sampled_trees_bdu <- read.tree(paste0(bduclnorigpath, "100.bd.ucln.subst.trees.directory/100.bd.ucln.subst.trees.sampled.trees"))
original_trees_bdu <- as.list(1:10)
original_roots_bdu <- 1:10
for (i in 1:length(names(sampled_trees_bdu))) {
  print(paste0("Simulated tree: ",i))
  holder_name <- names(sampled_trees_bdu)[[i]]
  print(holder_name)
  original_trees_bdu[[i]] <- original_posterior_bdu[[holder_name]]
  print(original_trees_bdu)
  original_roots_bdu[i] <- max(branching.times(original_trees_bdu[[i]]))
}

original_posterior_ys <- read.nexus(paste0(yulestrictorigpath, "100.yule.strict.(time).trees"))
sampled_trees_ys <- read.tree(paste0(yulestrictorigpath, "100.yule.strict.subst.trees.directory/100.yule.strict.subst.trees.sampled.trees"))
original_trees_ys <- as.list(1:10)
original_roots_ys <- 1:10
for (i in 1:length(names(sampled_trees_ys))) {
  print(paste0("Simulated tree: ",i))
  holder_name <- names(sampled_trees_ys)[[i]]
  print(holder_name)
  original_trees_ys[[i]] <- original_posterior_ys[[holder_name]]
  print(original_trees_ys)
  original_roots_ys[i] <- max(branching.times(original_trees_ys[[i]])) 
}

original_posterior_yu <- read.nexus(paste0(yuleuclnorigpath, "100.yule.ucln.(time).trees"))
sampled_trees_yu <- read.tree(paste0(yuleuclnorigpath, "100.yule.ucln.subst.trees.directory/100.yule.ucln.subst.trees.sampled.trees"))
original_trees_yu <- as.list(1:10)
original_roots_yu <- 1:10
for (i in 1:length(names(sampled_trees_yu))) {
  print(paste0("Simulated tree: ",i))
  holder_name <- names(sampled_trees_yu)[[i]]
  print(holder_name)
  original_trees_yu[[i]] <- original_posterior_yu[[holder_name]]
  print(original_trees_yu)
  original_roots_yu[i] <- max(branching.times(original_trees_yu[[i]]))
}

bds_bs <- lapply(list.files(bdstrictpath, full.names=TRUE, pattern="*bd-strict.annotated.tre"), read.nexus)
names(bds_bs) <- list.files(bdstrictpath, pattern="*bd-strict.annotated.tre")
bds_bu <- lapply(list.files(bdstrictpath, full.names=TRUE, pattern="*bd-ucln.annotated.tre"), read.nexus)
names(bds_bu) <- list.files(bdstrictpath, pattern="*bd-ucln.annotated.tre")
bds_ys <- lapply(list.files(bdstrictpath, full.names=TRUE, pattern="*yule-strict.annotated.tre"), read.nexus)
names(bds_ys) <- list.files(bdstrictpath, pattern="*yule-strict.annotated.tre")
bds_yu <- lapply(list.files(bdstrictpath, full.names=TRUE, pattern="*yule-ucln.annotated.tre"), read.nexus)
names(bds_yu) <- list.files(bdstrictpath, pattern="*yule-ucln.annotated.tre")

bdu_bs <- lapply(list.files(bduclnpath, full.names=TRUE, pattern="*bd-strict.annotated.tre"), read.nexus)
names(bdu_bs) <- list.files(bduclnpath, pattern="*bd-strict.annotated.tre")
bdu_bu <- lapply(list.files(bduclnpath, full.names=TRUE, pattern="*bd-ucln.annotated.tre"), read.nexus)
names(bdu_bu) <- list.files(bduclnpath, pattern="*bd-ucln.annotated.tre")
bdu_ys <- lapply(list.files(bduclnpath, full.names=TRUE, pattern="*yule-strict.annotated.tre"), read.nexus)
names(bdu_ys) <- list.files(bduclnpath, pattern="*yule-strict.annotated.tre")
bdu_yu <- lapply(list.files(bduclnpath, full.names=TRUE, pattern="*yule-ucln.annotated.tre"), read.nexus)
names(bdu_yu) <- list.files(bduclnpath, pattern="*yule-ucln.annotated.tre")

ys_bs <- lapply(list.files(yulestrictpath, full.names=TRUE, pattern="*bd-strict.annotated.tre"), read.nexus)
names(ys_bs) <- list.files(yulestrictpath, pattern="*bd-strict.annotated.tre")
ys_bu <- lapply(list.files(yulestrictpath, full.names=TRUE, pattern="*bd-ucln.annotated.tre"), read.nexus)
names(ys_bu) <- list.files(yulestrictpath, pattern="*bd-ucln.annotated.tre")
ys_ys <- lapply(list.files(yulestrictpath, full.names=TRUE, pattern="*yule-strict.annotated.tre"), read.nexus)
names(ys_ys) <- list.files(yulestrictpath, pattern="*yule-strict.annotated.tre")
ys_yu <- lapply(list.files(yulestrictpath, full.names=TRUE, pattern="*yule-ucln.annotated.tre"), read.nexus)
names(ys_yu) <- list.files(yulestrictpath, pattern="*yule-ucln.annotated.tre")

yu_bs <- lapply(list.files(yuleuclnpath, full.names=TRUE, pattern="*bd-strict.annotated.tre"), read.nexus)
names(yu_bs) <- list.files(yuleuclnpath, pattern="*bd-strict.annotated.tre")
yu_bu <- lapply(list.files(yuleuclnpath, full.names=TRUE, pattern="*bd-ucln.annotated.tre"), read.nexus)
names(yu_bu) <- list.files(yuleuclnpath, pattern="*bd-ucln.annotated.tre")
yu_ys <- lapply(list.files(yuleuclnpath, full.names=TRUE, pattern="*yule-strict.annotated.tre"), read.nexus)
names(yu_ys) <- list.files(yuleuclnpath, pattern="*yule-strict.annotated.tre")
yu_yu <- lapply(list.files(yuleuclnpath, full.names=TRUE, pattern="*yule-ucln.annotated.tre"), read.nexus)
names(yu_yu) <- list.files(yuleuclnpath, pattern="*yule-ucln.annotated.tre")

bds_bs_rs <- list()
bds_bu_rs <- list()
bds_ys_rs <- list()
bds_yu_rs <- list()

bdu_bs_rs <- list()
bdu_bu_rs <- list()
bdu_ys_rs <- list()
bdu_yu_rs <- list()

ys_bs_rs <- list()
ys_bu_rs <- list()
ys_ys_rs <- list()
ys_yu_rs <- list()

yu_bs_rs <- list()
yu_bu_rs <- list()
yu_ys_rs <- list()
yu_yu_rs <- list()

for (i in 1:10) {
	posterior_tree <- bds_bs[[i]]
	val <- as.numeric(strsplit(names(bds_bs[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bds[val]/max(branching.times(posterior_tree))
	bds_bs_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- bds_bu[[i]]
	val <- as.numeric(strsplit(names(bds_bu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bds[val]/max(branching.times(posterior_tree))
	bds_bu_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- bds_ys[[i]]
	val <- as.numeric(strsplit(names(bds_ys[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bds[val]/max(branching.times(posterior_tree))
	bds_ys_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- bds_yu[[i]]
	val <- as.numeric(strsplit(names(bds_yu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bds[val]/max(branching.times(posterior_tree))
	bds_yu_rs[[i]] <- posterior_tree
}


for (i in 1:10) {
	posterior_tree <- bdu_bs[[i]]
	val <- as.numeric(strsplit(names(bdu_bs[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bdu[val]/max(branching.times(posterior_tree))
	bdu_bs_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- bdu_bu[[i]]
	val <- as.numeric(strsplit(names(bdu_bu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bdu[val]/max(branching.times(posterior_tree))
	bdu_bu_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- bdu_ys[[i]]
	val <- as.numeric(strsplit(names(bdu_ys[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bdu[val]/max(branching.times(posterior_tree))
	bdu_ys_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- bdu_yu[[i]]
	val <- as.numeric(strsplit(names(bdu_yu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_bdu[val]/max(branching.times(posterior_tree))
	bdu_yu_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- ys_bs[[i]]
	val <- as.numeric(strsplit(names(ys_bs[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_ys[val]/max(branching.times(posterior_tree))
	ys_bs_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- ys_bu[[i]]
	val <- as.numeric(strsplit(names(ys_bu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_ys[val]/max(branching.times(posterior_tree))
	ys_bu_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- ys_ys[[i]]
	val <- as.numeric(strsplit(names(ys_ys[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_ys[val]/max(branching.times(posterior_tree))
	ys_ys_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- ys_yu[[i]]
	val <- as.numeric(strsplit(names(ys_yu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_ys[val]/max(branching.times(posterior_tree))
	ys_yu_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- yu_bs[[i]]
	val <- as.numeric(strsplit(names(yu_bs[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_yu[val]/max(branching.times(posterior_tree))
	yu_bs_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- yu_bu[[i]]
	val <- as.numeric(strsplit(names(yu_bu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_yu[val]/max(branching.times(posterior_tree))
	yu_bu_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- yu_ys[[i]]
	val <- as.numeric(strsplit(names(yu_ys[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_yu[val]/max(branching.times(posterior_tree))
	yu_ys_rs[[i]] <- posterior_tree
}

for (i in 1:10) {
	posterior_tree <- yu_yu[[i]]
	val <- as.numeric(strsplit(names(yu_yu[i]), split="\\.")[[1]][4])
	posterior_tree$edge.length <- posterior_tree$edge.length*original_roots_yu[val]/max(branching.times(posterior_tree))
	yu_yu_rs[[i]] <- posterior_tree
}

make_ltt_gray <- function(tree) {
	bt <- branching.times(tree)
	bt <- sort(bt, decreasing=FALSE)
	names(bt) <- NULL
	dat <- cbind(bt, n=sort(seq(2, length(bt)+1, 1), decreasing=TRUE))
	points(x=dat[,1], y=log(dat[,2]), type="l", xlab="Time", ylab="log(number of lineages)", xlim=c(5, 0), col="gray", lwd=3)
}

make_ltt_black <- function(tree) {
	bt <- branching.times(tree)
	bt <- sort(bt, decreasing=FALSE)
	names(bt) <- NULL
	dat <- cbind(bt, n=sort(seq(2, length(bt)+1, 1), decreasing=TRUE))
	points(x=dat[,1], y=log(dat[,2]), type="l", xlab="Time", ylab="log(number of lineages)", xlim=c(5, 0), col="black")
}

pdf(file="ltt_plots.pdf", width=8.5, height=11)
par(mfrow=c(4,4), oma=c(3, 5, 3, 2.5), mar=c(4,2,2,0), mgp=c(3,1,0),
	xpd=NA, bty="l")
plot(xlab="", ylab="log(lineages)",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n")
axis(1, at=c(0:5), labels=NA)
lapply(original_trees_bds, make_ltt_gray)
lapply(bds_bs_rs, make_ltt_black)
plot( type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n", xlab="", ylab="")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_bds, make_ltt_gray)
lapply(bds_bu_rs, make_ltt_black)
plot( type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n", xlab="", ylab="")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_bds, make_ltt_gray)
lapply(bds_ys_rs, make_ltt_black)
plot( type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n", xlab="", ylab="")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_bds, make_ltt_gray)
lapply(bds_yu_rs, make_ltt_black)

plot(xlab="", ylab="log(lineages)",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n")
axis(1, at=c(0:5), labels=NA)
lapply(original_trees_bdu, make_ltt_gray)
lapply(bdu_bs_rs, make_ltt_black)
plot(xlab="", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_bdu, make_ltt_gray)
lapply(bdu_bu_rs, make_ltt_black)
plot(xlab="", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_bdu, make_ltt_gray)
lapply(bdu_ys_rs, make_ltt_black)
plot(xlab="", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_bdu, make_ltt_gray)
lapply(bdu_yu_rs, make_ltt_black)

plot(xlab="", ylab="log(lineages)",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_ys, make_ltt_gray)
lapply(ys_bs_rs, make_ltt_black)
plot(xlab="", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_ys, make_ltt_gray)
lapply(ys_bu_rs, make_ltt_black)
plot(xlab="", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_ys, make_ltt_gray)
lapply(ys_ys_rs, make_ltt_black)
plot(xlab="", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), xaxt="n", yaxt="n")
axis(1, at=c(0:5), labels=NA)
axis(2, at=c(0:5), labels=NA)
lapply(original_trees_ys, make_ltt_gray)
lapply(ys_yu_rs, make_ltt_black)

plot(xlab="Time", ylab="log(lineages)",  type="n", x=0:5, y=0:5, xlim=c(5,0))
lapply(original_trees_yu, make_ltt_gray)
lapply(yu_bs_rs, make_ltt_black)
plot(xlab="Time", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), yaxt="n")
lapply(original_trees_yu, make_ltt_gray)
lapply(yu_bu_rs, make_ltt_black)
plot(xlab="Time", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), yaxt="n")
lapply(original_trees_yu, make_ltt_gray)
lapply(yu_ys_rs, make_ltt_black)
plot(xlab="Time", ylab="",  type="n", x=0:5, y=0:5, xlim=c(5,0), yaxt="n")
lapply(original_trees_yu, make_ltt_gray)
lapply(yu_yu_rs, make_ltt_black)

mtext(text=c("BD:Strict", "BD:UCLN", "Yule:Strict", "Yule:UCLN"), side=2, at=c(.8875, .6375, .385, .1375) , outer=TRUE, line=2, font=2, cex=0.8)
mtext(text=c("BD:Strict", "BD:UCLN", "Yule:Strict", "Yule:UCLN"), side=3, at=c(.14, .385, .64, .89) , outer=TRUE, line=0, font=2, cex=0.8)

dev.off()