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

# dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
# 
# install.packages('ips',  Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" ) 
# install.packages("tree", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
# install.packages("devtools", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
# install.packages("RSQLite", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
# install.packages("phangorn", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
# install.packages('phytools', Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

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
                  'phylogenetics/20170916-TreeAnnotator/',
                  sep = ''), full.names = T)


beast.trees <- lapply(trees, read.beast)
plot(beast.trees[[1]])

#beast.trees2 <- beast.trees[1:5]
#---------------------------#
##_____ Ss4R distance _____##
#---------------------------#

# remove the trees without branchlengths
beast.trees.na <- NULL
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


     
nodeages.Ss4R <- data.frame(Ss4R.nodeage.est, Ss4R.nodeage.max, Ss4R.nodeage.min,
                            ssal.dup1, ssal.dup2, stringsAsFactors = F)
colnames(nodeages.Ss4R) <- c('estimate', 'max', 'min', 'ssal1', 'ssal2')

write.table(nodeages.Ss4R, file = "20170916-nodeages_Ss4r.csv", append = F, 
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

# salmon.tips <- which(substr(all.trees[[1]]$tip.label, 1, 4) == 'Ssal')
# salmon.dups <- substr(all.trees[[1]]$tip.label[salmon.tips], 6, 1000)

# n = 1
# for (i in 1:length(beast.trees)){
#   salmon.tips <- which(substr(beast.trees[[i]]$tip.label, 1, 4) == 'Ssal')
# }
# 
# i <- 1

# find the tree tips that have salmon
salmon.tips <- lapply(all.trees, function(i) {
  which(substr(i$tip.label, 1, 4) == 'Ssal')
})

# get the names of the salmon proteins
salmon.dups <- NULL
for (i in 1:length(all.trees)){
  salmon.dups[[i]] <- substr(all.trees[[i]]$tip.label[salmon.tips[[i]]], 6, 1000)
}


# find the corresponding chromosome positions of the proteins on the genome
dup.chrom.df <- NULL
for (i in 1:length(all.trees)){
  proteinid <- get.gtf(get.id(salmon.dups[[i]], id.type = 'protein')$gene_id)
  dup.chrom.pos <- as.data.frame(c(proteinid[1,2:4], proteinid[2,2:4]), 
                                 stringsAsFactors = F)
  dup.chrom.df <- rbind(dup.chrom.df, dup.chrom.pos)
}



Ss4R.time.chrompos <- data.frame(dup.chrom.df, nodeages.Ss4R[,1:3] )
Ss4R.time.chrompos.ord <- Ss4R.time.chrompos[order(Ss4R.time.chrompos[,1], 
                                                  Ss4R.time.chrompos[,2]),]

colnames(Ss4R.time.chrompos.ord) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", 
                            "estimate", "max", "min")

write.table(Ss4R.time.chrompos.ord, file = "20170918-Ss4r_time&chropos.csv", append = F, 
            sep = ";", quote = F, col.names = T, row.names = F)



#-------------------------------------------------#
##_____ Plotting nodeages Ss4R and chro pos _____##
#-------------------------------------------------#

for (chr1 in unique(Ss4R.chrom.pos.ord$chr1)){
  cat (chr1, "\n")
  all.rows <- which(Ss4R.chrom.pos.ord$chr1 == chr1)
  chr2.df <- Ss4R.chrom.pos.ord[all.rows, 4]
  for (chr2 in unique(chr2.df)){
    cat (chr2, "\n")
    #barplot(x = Ss4R.chrom.pos.ord)
  }
  
}

#--------------------------------------#
##_____ example tree, nodelabels _____##
#--------------------------------------#


all.trees[[7]]$tip.label <- c("Rainbow trout (Oncorhynchus mykiss)", 
                              "Grayling (Thymallus Thymallus)", 
                              "Coho salmon (Oncorhynchus kisutch)",
                              "Atlantic salmon (Salmo salar)", 
                              "Atlantic salmon (Salmo salar)", 
                              "Zebrafish (Danio rerio)", 
                              "Rainbow trout (Oncorhynchus mykiss)", 
                              "Grayling (Thymallus Thymallus)", 
                              "Coho salmon (Oncorhynchus kisutch)",
                              "Northern pike (Esox Lucius)")
plot(all.trees[[7]])
node.ages <- round(all.trees[[7]]$height, digits = 1)
nodelabels(node.ages, frame = "r", bg = "white")


#--------------------------------#
##_____ average nodelabels _____##
#--------------------------------#
library(dplyr)
library(plyr)
install.packages("Stack")
library(stack)
library(ggplot2)
library(grid)
library(stringr)


Ss4R.time.chr <- read.table(file = "20170918-Ss4r_time&chropos.csv", 
                    sep = ";", stringsAsFactors = F, header = T)

Ss4R.chr <- read.table(file = "20170918-Ss4r_time_est_and_chro.csv", 
                            sep = ";", stringsAsFactors = F, header = T)
Ss4R.chr$chr1 <- as.character(str_pad(Ss4R.chr$chr1, 2, pad = "0"))



boxplot(estimate~chr1, data = Ss4R.chr, xlab = "chromosome number",
        ylab = "divergence time after Ss4R (Mya)", col = rainbow(12))

fill <- "#4271AE"
line <- "#1F3552"

boxplot.Ss4R <- ggplot(Ss4R.chr, aes(x=chr1, y=estimate)) +
  geom_boxplot(fill = fill, color = line, alpha = 0.7) +
  theme_bw() +
  scale_x_discrete(name = "Chromosome number") +
  scale_y_continuous(name = "Divergence time after Ss4R (Mya)") +
  theme(text = element_text(size = 14),
        axis.text= element_text(size = 13)) 

boxplot.Ss4R


#-----------------------------------------#
##_____ Make data for circlize plot _____##
#-----------------------------------------#

# order by chromosomeno
dup.chrom.df.ord <- dup.chrom.df[order(dup.chrom.df[,1], 
                                       dup.chrom.df[,2]),]




#get.id('*') # get all genes/proteins/transcripts
#proteinid <- get.gtf(get.id(salmon.dups, id.type = 'protein')$gene_id)
#dup.chrom.pos <- as.data.frame(c(proteinid[1,2:4], proteinid[2,2:4]))

# retrieve chromosome number
# give different color per chromosome

colors <- c('cigene_21_1', 'cigene_21_2', 'cigene_21_3', 'cigene_21_4', 'cigene_21_5', 
            'cigene_21_6', 'cigene_21_7', 'cigene_21_8', 'cigene_21_9', 'cigene_21_10',
            'cigene_21_11', 'cigene_21_12', 'cigene_21_13', 'cigene_21_14', 'cigene_21_15',
            'cigene_21_16', 'cigene_21_17', 'cigene_21_18', 'cigene_21_19', 'cigene_21_20',
            'cigene_21_21'  )

dup.chrom.df.ord$color <- NA
colorit <- 0
for (chro in unique(dup.chrom.df.ord$chr1)){
  colorit <- sum(colorit, 1)
  chro.pos <- which(dup.chrom.df.ord$chr1 == chro)
  dup.chrom.df.ord[chro.pos, 7] <- colors[colorit]
}

colorit



write.table(dup.chrom.df.ord, file = "20170918-duplicates_chrompos.txt", 
              append = F, col.names = F, row.names = F, sep = "\t", quote = FALSE)
