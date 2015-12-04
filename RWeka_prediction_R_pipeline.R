# prediction with RWeka
library(RWeka)
library(rJava)

# Add the path and file name that you wish to analyse 
tree_stats_path <- "ADD PATH HERE" # Example: "~/Desktop/ML_tree_stats.txt"
nuc_diversity_path <- "ADD PATH HERE" # Example: "~/Desktop/nuc_diversity.txt"

# After download the classifier you should load it in R using the following command for prediction of empirical data
load("classifier_1000bp.dat") # Classifier that should be used for predictions when DNA sequence alingment is approximately 1,000 bp
load("classifier_10000bp.dat") # Classifier that should be used for predictions when DNA sequence alingment is approximately 10,000 bp

tree_stats <- read.delim(tree_stats_path, header = FALSE)
nuc_diversity <- read.delim(nuc_diversity_path, header = FALSE)

# Adding names to collumns
colnames(tree_stats) <- c("file_name", "Beta_Index", "maxLadder", "ILportion", "maxDepth", "maxWidth", "woverd", "maxDiffWidth", "nCherry", "prop")
colnames(nuc_diversity) <- c("file_name", "nuc_diversity")

# Joining data to a single data frame
# MAKE SURE THAT THE FILE NAME ORDER IN "tree_stats" ARE THE SAME AS IN "nuc_diversity" (if you are analysing more than one tree/DNA sequence alignment)
allData <- cbind(tree_stats, nuc_diversity$nuc_diversity)
colnames(allData)[11] <- "nuc_diversity"

# Prediction function to be used with trees reconstructed with DNA sequence alignments of approximatelly 1,000 bp. This will give the result in probability. 
# For the example in the paper, this classifier will be used to predict data for env E and pol genes (rows 1 and 4 for the prediction results).
predict(classifier1000bp, newdata = allData, type = "probability")

# Prediction function to be used with trees reconstructed with DNA sequence alignments of approximatelly 10,000 bp. This will give the result in probability.
# For the example in the paper, this classifier will be used to predict data for all PERVs and PERV-B genomes (rows 2 and 3 for the prediction results)
predict(classifier10000bp, newdata = allData, type = "probability")




