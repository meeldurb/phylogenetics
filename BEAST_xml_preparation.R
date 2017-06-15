#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##preparing files for BEAST analysis
###################################################################





#####___________install & load packages and functions____________####

install.packages("ape", repos="http://cran.rstudio.com/")

library("ape")
source("D:/Beast_dating_salmonids/Beast_Rscripts/orthologfinder_oct2016.R")

####____________load data______________####

# prepare homeolog table

# get homologs of Oncorhynchus mykiss (rainbow trout)
load(file.path("C:/Users/meeldurb/Google Drive",
    "Master internship phylogenetics salmonids",
    "Salmonid_genomics_resources/Orthologs_homeologs/Homeologs",
    "OmykV6_2016_best_in_homelogRegions_minpident80_mincov50.RData"))
     
# change to use env here because of same file names 
# only select first 2 columns, these contain the homologues
omyk.hom <- putative_homeologs_pepblast[,1:2] 
# add organism naming in front on homologues genes
omyk.hom = apply(omyk.hom, 2, function(i) paste('Omyk2|', i, sep=''))

# get homologs of Salmo salar (Atlantic salmon)
load(file.path("C:/Users/meeldurb/Google Drive", 
    "Master internship phylogenetics salmonids", 
    "Salmonid_genomics_resources/Orthologs_homeologs/Homeologs",
    "RefSeq_GarethLongest_2016_best_in_homelogRegions_minpident80_mincov50.RData"))
# change to use env here because of same file names 
ssal.hom <- putative_homeologs_pepblast[,1:2]
# add organism naming in front on homologues genes
ssal.hom <- apply(ssal.hom, 2, function(i) paste('Ssal|', i, sep=''))

# bind all homologs in complete dataframe
hom.table <- data.frame(rbind(omyk.hom, ssal.hom), 
                        stringsAsFactors = F)
head(hom.table)


# get clans
load(file.path("C:/Users/meeldurb/Google Drive",
"Master internship phylogenetics salmonids",
"Salmonid_genomics_resources/Orthologs_homeologs",
"Gene_trees/OMYK_Orthofinder_trees_October2016_clanfinder_clans.RData"))

clans <- Orthofinder_trees_October2016_clanfinder_clans
# add a dot behind the gene names to make it suited for unlisting
# names(clans) <- paste(names(clans), '.', sep='')
# Some labels of the genes need re-labelling
# Tthy needs to have a | behind and 
clans <- lapply(clans, function(i) {
  tip.lab <- gsub('Tthy', 'Tthy|Tthy', i$tip.label)
  i$tip.label <- tip.lab
  i
  })

## Example tree below
# select only the clans tht have at least 1 E.luc 
# and select only the clans that have 
# at least 2 of the Ssal, Omyk or Tthy species.
clans.selection = clans[sapply(clans, function(i) {
   length(grep('Eluc', i$tip.label))>0 & 
    length(grep('Ssal|Omyk|Tthy', i$tip.label))>=2
  })]
length(clans.selection)
plot(clans.selection[[1000]])

test.all = lapply(clans.selection, get.ortholog, verbouse=F)
# table(sapply(test.all, function(i) i$class))
# test.og = get.ortholog(clans.selection[[1000]])
# tree = clans.selection[[1000]]
