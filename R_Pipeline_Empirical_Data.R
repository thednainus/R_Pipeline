############################################################################################################
## Pipeline in R to claculate sytatistics described in Nascimento and Rodrigo (submitted for publicatin)  ##
## empirical data. More information will be added here when paper is accepeted for publication            ##
############################################################################################################

# Tree shape statistics
# (i) Beta-sppliting model = Beta_Index
# (ii) Ladder length = maxLadder
# (iii) "IL" branches = ILportion
# (iv) Maximum depth = maxDepth
# (v) Maximum width = maxWidth
# (vi) Maximum width over maximum depth = woverd
# (vii) Maximum difference in widths = maxDiffWidth
# (viii) Number of cherries = nCherry

# Other statistics
# (ix) Proportion of terminal branch lengths that contributed to the total tree branch length = prop
# (x) Nucleatide diversity for DNA sequence alignments

# Sets work directory to Desktop. You can change the path to any other directory
setwd("~/Desktop/sus/result/")

# Install packages below if not already installed. Use install.packages("packageNameHere").
library(apTreeshape)
library(ape)
library(phangorn)
library(phylobase) # required to install phyloTop
library(NHPoisson) # required to install phyloTop
library(igraph) # required to install phyloTop
library(phyloTop)
library(pegas)


# Read file names in a matrix. If more than one tree/DNA sequence alignment will be analyzed, make sure that the final file with the statistics values
# are in the same order. This will be important for prediction of the ERV model with RWeka (see R code about prediction with RWeka for ERV models)
  File_names_tree  <- matrix(list.files(path="~/Desktop/sus/tree/", pattern ="*.tre", recursive = FALSE, include.dirs = FALSE))
  File_names_ali <- matrix(list.files(path="~/Desktop/sus/alignment/", pattern ="*.txt", recursive = FALSE, include.dirs = FALSE))


## Function to analyse reconstructed phylogenetic trees (ML trees) in which each tree should be in a separate file.
## If trees are not rooted, you should root the trees using a midpoint root, for example.
## Please comment the line below about midpoint root if your trees are already rooted.

ML_tree_stats <- function(x) {
  
  fileName <- x
  
  # This part of the code will create a file with the name "ML_trees_statistics.txt" as below. You can change this file name.
  # You should change this part of the code in accordance to your OS.
  Path = "~/Desktop/sus/result/ML_tree_stats.txt"
  
  
  # You should paste the test of the path to "x", to read phylogenetic trees. You should change this part of the code in accordance to your OS.
  x <- paste("~/Desktop/sus/tree/", x, sep="")
  tree <- read.tree(x)
  tree <- midpoint(tree) # YOU SHOULD COMMENT THIS LINE IF YOUR TREES ARE ALREADY ROOTED
  
  # First I transfor "tree" as treeshape to calulate the Beta_Index
  Tree <- as.treeshape(tree)
  Beta_Index <- maxlik.betasplit(Tree)
  
  # Then, I transfor "tree" as "phylo4" format to calculate the other shape statistics
  tree1 <- as(tree,"phylo4")
  
  n <- nTips(tree)
  
  # Calculates the statistics as explained in the beginning. To try to account for tree size, some values were divided by the number of tips in the tree
  maxLadder <- (max(ladderNums(tree1)))/n
  ILportion <- ILnumber(tree1)/n
  maxDepth  <- (max(dists(tree1)))/n
  maxWidth <- max(widths(tree1))/n
  woverd <- max(widths(tree1))/(max(dists(tree1)))
  maxDiffWidth <- max(diff(widths(tree1)))/n
  nCherry <- (nConfig(tree1,2))/n
  
  # to calculate the proportion of terminal branch lengths that contributed to the total tree branch length 
  list_tip_branches <- which(tree$edge[,2] <= n)
  total_tip_length = 0
  
  for (i in list_tip_branches){
    total_tip_length = total_tip_length + tree$edge.length[i]
  }
  
  total_tree_length = sum(tree$edge.length)
  prop <- total_tip_length/total_tree_length
  
  Information <- data.frame(fileName, Beta_Index$max_lik, maxLadder, ILportion, maxDepth, maxWidth, woverd, maxDiffWidth, nCherry, prop)
  write.table(Information, file = Path, append = TRUE, col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
}


## Function to calculate nucleotide diversity from a DNA sequence alignment 
nuc_diversity_empirical <- function(x){
  
  fileName <- x
  
  # You should change this name is accordance to your OS. You can change the name of the file to any other meaninful name
  Path = "~/Desktop/sus/result/nuc_diversity.txt"
  
  x <- paste("~/Desktop/sus/alignment/", x, sep="") # You should change this path name is accordance to your OS
  
  Seq_ali <- read.FASTA(x) # Read DNA sequence alignment
  Nuc.Diversity <- nuc.div(Seq_ali, pairwise.deletion = TRUE) # Calculates the nucleotide diversity
  
  Information <- data.frame(fileName, Nuc.Diversity)
  write.table(Information, file = Path, append = TRUE, col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  
}


## apply to the matrix "File_names_tree" the function "ML_tree_stats"
apply(File_names_tree, 1, ML_tree_stats)

## apply to the matrix "File_names_ali" the "nuc_diversity_empirical" function
apply(File_names_ali, 1, nuc_diversity_empirical)
warnings()






 
