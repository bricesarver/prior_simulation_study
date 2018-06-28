# The impact of tree prior and molecular clock models on the estimate of diversificaton rates

Here, I provide code and methods to simulate trees, simulate sequence data, estimate trees 
using BEAST, and estimate diversification parameters.

## Dependencies
In addition to R, you will need:
* [BEAST v1.7.5](https://code.google.com/p/beast-mcmc/downloads/list?can=4&q=)

This approach ought to work with BEAST v2.0+, though I expect there to be differences in the XML input 
that may not be compatible with BEASTifier.

* [Seq-Gen](http://tree.bio.ed.ac.uk/software/seqgen/)
* [BEASTifier](https://github.com/josephwb/BEASTifier)

You will also need several R libraries:
* parallel
* ape
* laser
* TreeSim
* geiger

You can install these, and their dependencies, from within R:
```
install.packages("ape")
```

Alternatively, you can download the source code, compile, and install locally if you are 
working on a system where you lack sufficient privileges.


## Execution

### Step 1: Simulate Trees

BEAUti v1.7.5 requires a NEXUS file to be loaded before you can generate the BEAST XML file.  
Simulate a tree, under any conditions, with the number of tips you require.  Then, using Seq-Gen, 
simulate DNA sequence data, under any conditions, using the tree to guide the simulation.  We will 
be sampling from the prior only, so all nucleotide sequence data will be replaced with 'NNN' when 
the XML file is generated.

An example of simulating some trees is in **simulate_trees.R**.

### Step 2: Select priors

Prior distributions, and operators, can be specified within BEAUti.  If you don't want parameters 
to change from their initial value, remove the operator associated with them.

If you have an XML template, it is possible to automate this approach.  Alternatively, it can be 
done by hand.

**Caution:** BEAST will fail to execute if the starting tree is not compatible with constraints imposed 
by fixing priors.  It is advisable to incorporate some 'wiggle room' into certain parameters, just to be safe.  
I have found that fixing the root age at a strict value (by, say, using a narrow uniform prior) often halts 
execution, whereas placing a prior without hard boundary conditions (like, say, a normal prior with a mean of the value you want 
and a small standard deviation) works much better.

### Step 3: Generate a posterior distribution of trees

Let BEAST run for 10,000,000 generations (or longer!) sampling every 1000.  Be sure to name your 
output something easily identifiable.

The first few trees change drastically - I remove a 10% burnin for this reason.  Sampling from a stationary distribution 
is achieved quickly.  Keep the chronograms **AND** phylograms.

### Step 4: Select a subset of these trees and simulate nucleotide datasets.

**sequence_simulation.R** contains code that samples from the 'simulated' posterior distribution, 
rescales the tree, and simulates DNA sequences under a model using the tree.  Data is simulated using 
phylograms (i.e., chronograms with branches multiplied by their rate scalars, in this context).  We multiply branches 
by 0.01 to scale the trees to appropriate levels; alternatively, we could used prior distributions with smaller means for 
molecular clock parameters and avoid this step.

This also calls BEASTifier.  Make sure that you have appropriately generated the BEASTifier configuration file 
for your needs.

### Step 5: Estimate and process results

This step is computationally intensive.  One BEAST run is performed for each XML file.  I recommend determining 
the number of generations needed to achieve a reasonable ESS and stationarity.  For the 100 taxa datasets with the 
most complex tree prior and clock model (BD:UCLN), I found that 50,000,000 (sampling every 5000) 
was sufficient.  I performed hundreds of runs in parallel on a distributed cluster; if you do not 
have access to such a system, your analysis will be bottlenecked by the number of CPUs.  RAM should 
not be an issue, though removing R objects from memory using ```rm()``` and performing garbage collection using ```gc()``` 
can alleviate this issue once parameters are summarized.

Even with 40-60 CPUs on a single system (note that distributed systems may be bottlenecked by 
bandwidth and bus limitations), analyzing hundreds of thousands of trees is not quick.  My analyses 
took approximately 20 (wall time) minutes per simulation condition.

**tree_processing.R** contains code, using parallel processing, that rescales trees and calculates 
MLEs of speciation and net diversification for use in plots (or summary tables).

**make_ltt_plots.R** contains code to make lineage-through-time plots from maximum clade credibility 
trees summarized using TreeAnnotator.
