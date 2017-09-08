#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch & Simen Rod Sandve
##Script for summarizing the information from the BEAST trees
##into a single "target" tree
###################################################################



#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

# reading BEAST trees
install.packages('ips',  repos = "http://cran.rstudio.com/" )  
install.packages("tree", repos = "http://cran.rstudio.com/")
install.packages("devtools", repos = "http://cran.rstudio.com/")
install.packages("RSQLite", repos = "http://cran.rstudio.com/")
install.packages("phangorn", repos = "http://cran.rstudio.com/")
install.packages('phytools', repos = "http://cran.rstudio.com/")

library(ips)
library(tree)
library(devtools)
library(RSQLite)
library(phangorn)
library(phytools)

install_github("FabianGrammes/Ssa.RefSeq.db")
library(Ssa.RefSeq.db)

#------------------------------------------------------#
##_____ Summarize BEAST trees with TreeAnnotator _____##
#------------------------------------------------------#

#system('sbatch -a 1-633 %20 TreeAnnotator_submit.sh 10')



#-----------------------------#
##_____ Visualize trees _____##
#-----------------------------#


trees = dir(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                  'Master_internship_phylogenetics/',
                  'phylogenetics/TreeAnnotator_results/',
                  sep = ''), full.names = T)


beast.trees <- lapply(trees, read.beast)
plot(beast.tr[[1]])

#beast.trees2 <- beast.trees[1:5]
#---------------------------#
##_____ Ss4R distance _____##
#---------------------------#

# remove the trees without branchlengths
for (i in 1:length(beast.trees)){
  #print (tree)
  #i <- sum(i, 1)
  print(i)
  if (is.na(beast.trees[[i]]$rate[2])){
    beast.trees[[i]] <- NA
    }
}

all.trees <- beast.trees[!is.na(beast.trees)]


# plot(all.trees[[2]])
# nodelabels()
# tiplabels()

# find Omyk and Ssal in the tip labels and get MRCA 
Ss4R.nodes <- sapply(all.trees, function(i){
  getMRCA(i, c(substr(i$tip.label, 1, 4) == 'Ssal', 
               substr(i$tip.label, 1, 4) == 'Omyk'))
  })


# get node age of the Ss4R node
height.pos <- sapply(Ss4R.nodes, function(i){
  i }) - sapply(all.trees, function(i){
  length(i$tip.label) 
  }) 


# get all the node age estimates
# and the salmon duplicates

Ss4R.nodeage.est <- NULL
Ss4R.nodeage.min <- NULL
Ss4R.nodeage.max <- NULL
ssal.dup1 <- NULL
ssal.dup2 <- NULL

i <- 0
for (pos in height.pos){
  #cat(pos, "\n")
  i <- sum(i, 1)
  Ss4R.nodeage.est <- c(Ss4R.nodeage.est, all.trees[[i]]$height[pos])
  Ss4R.nodeage.min <- c(Ss4R.nodeage.min, all.trees[[i]]$`height_95%_HPD_MIN`[pos])
  Ss4R.nodeage.max <- c(Ss4R.nodeage.max, all.trees[[i]]$`height_95%_HPD_MAX`[pos])
  # find salmon dups in the tree
  salmon.tips <- which(substr(all.trees[[i]]$tip.label, 1, 4) == 'Ssal')
  salmon.dups <- substr(all.trees[[i]]$tip.label[salmon.tips], 6, 1000)
  ssal.dup1 <- c(ssal.dup1, salmon.dups[1])
  ssal.dup2 <- c(ssal.dup2, salmon.dups[2])
  
}


     
nodeages.Ss4R <- as.data.frame(cbind(Ss4R.nodeage.est, Ss4R.nodeage.max, Ss4R.nodeage.min,
                                     ssal.dup1, ssal.dup2))
colnames(nodeages.Ss4R) <- c('estimate', 'max', 'min', 'ssal1', 'ssal2')

write.table(nodeages.Ss4R, file = "nodeages_Ss4r.csv", append = F, 
            sep = ";", quote = F, col.names = T, row.names = F)



# tree = rtree(10)
# plot(beast.trees[[100]])
# nodelabels()
# tiplabels()
# Ancestors(beast.trees[[100]], 1:3, "all")
# Children(tree, 11)
# Descendants(tree, 11, "tips")
# Siblings(tree, 3)
#


#-------------------------------------#
##_____ duplicate chr positions _____##
#-------------------------------------#

salmon.tips <- which(substr(tree$tip.label, 1, 4) == 'Ssal')
salmon.dups <- substr(tree$tip.label[salmon.tips], 6, 1000)


#get.id('*') # get all genes/proteins/transcripts

proteinid <- get.gtf(get.id(salmon.dups, id.type = 'protein')$gene_id)

