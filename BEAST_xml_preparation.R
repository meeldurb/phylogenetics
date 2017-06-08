#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##preparing files for BEAST analysis
###################################################################


####____________load data______________####

# prepare homeolog table

# get homologs of Oncorhynchus mykiss (rainbow trout)
load(file.path("C:/Users/meeldurb/Google Drive/",
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
    "Salmonid_genomics_resources/Orthologs_homeologs/Homeologs/",
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
# add a dot behind the gene names
names(clans) <- paste(names(clans), '.', sep='')
# Some labels of the genes need re-labelling
# Tthy needs to have a | behind and 
clans <- lapply(clans, function(i) {
  tip.lab <- gsub('Tthy', 'Tthy|Tthy', i$tip.label)
  tip.lab <- sub('\\|Gene.*', '', tip.lab)
  i$tip.label <- tip.lab
  #i 
  })
