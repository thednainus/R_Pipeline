# R-pipeline
version 1.0

R pipeline to analyse empirical data for endogenous retroviruses (ERVs) or other transposable elements (TEs). Codes were tested using R version 3.2.2 for Mac computer.

Use R code “R_Pipeline_Empirical_Data.R” to calculate all statistics for phylogenetic trees and DNA sequence alignments.

Upload file “classifier_1000bp.dat” or “classifier_10000bp.dat” to your computer, and use R code “RWeka_predication_R_pipeline.R” to predict the model of ERV/TE evolution according to simulations in Nascimento and Rodrigo (submitted for publication). 


In alignments folder you can find the alignments for PERVs used for analysis of empirical data in Nascimento and Rodrigo (submitted for publication).

-------------------------

To run “R_Pipeline_Empirical_Data.R” the following R packages a necessary:

apTreeshape version 1.4.5
ape version 3.3
phangorn version 1.99.14
phylobase version 0.8.0 (required to install phyloTop)
NHPoisson version 3.1 (required to install phyloTop)
igraph version 1.0.1 (required to install phyloTop)
phyloTop version 1.1.1
pegas version 0.8.2

To run “RWeka_predication_R_pipeline.R” the following R packages a necessary:

RWeka version 0.4.24
rJava version 0.9.7

Even though the above versions were used in the latest analyses of the data. New versions should also work.

For more information on the codes please read comments in both R codes.

-------------------------

For any questions or comments please e-mail me at thednainus@yahoo.com


-------------------------

Reproduction or any form of use of these codes should be acknowledged
