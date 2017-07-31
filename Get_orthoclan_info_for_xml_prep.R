#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch & Simen Rod Sandve
##Script for loading and manipulating the orthoclan results
##and preparing them to be imported in the BEAUTI .xml files
###################################################################


#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

#install.packages('RCurl',  repos = "http://cran.rstudio.com/")
#install.packages("ape", repos = "http://cran.rstudio.com/")
#install.packages("phangorn", repos = "http://cran.rstudio.com/")
#install.packages("plyr", repos = "http://cran.rstudio.com/")

library(RCurl)
library(ape)
library(phangorn)
library(seqinr)
library(plyr)



source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenomics/clanfinder.R', sep = ''))
source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenomics/Phylo_functions.R', sep = ''))
source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenetics/loadRData.R', sep = ''))



#-------------------------------------------------#
##____ make clans from orthogroup trees ____##
#-------------------------------------------------#

# 1) load orthotrees

# - this is the output from protein trees made 
# from orthogroup sets of sequences (from orthofinder)
OG_trees <- loadRData(paste('C:/Users/meeldurb/Google Drive/Master internship ',
                    'phylogenetics salmonids/Salmonid_genomics_resources/',
                    'Orthologs_homeologs/orthogroups.03.06.2017/',
                    'OG_trees.30.05.2017.RData', sep = ''))

# total number of orthogroups
length(OG_trees) 
# summing up the number of specific species in the trees
table(sapply(OG_trees, function(i) sum(c('Ssal', 'Omyk') %in% substr(i$tip.label, 1, 4))))

# 2) make clans using clanfinder

OG_clanfinder = lapply(OG_trees, clanFinder, ut = c('Olat', 'Gacu', 'Drer', 'Locu', 'Mmus', 'Hsap'))

# # Play a bit with the datastructure
# names(OG_trees)
# # the names of the trees in the different clans in the OG groups are named by a integer
# names(OG_clanfinder)
# names(OG_clanfinder[1])
# names(OG_clanfinder[[1]])
# names(OG_clanfinder$OG0000000.)
# OG_clanfinder[[1]][1]

# 3) fixing names

# some of the OG clans are empty, we need to remove these
OG_clanfinder.filt <- OG_clanfinder[sapply(OG_clanfinder, length)>0]
clans.num <- as.numeric(unlist(sapply(OG_clanfinder.filt, function(i) 1:length(i))))
OG_clans <- unlist(OG_clanfinder.filt, recursive = F)
head(names(OG_clans))
# get rid of the 2 dots and make an recursive integer for the clans of the OG's 
names(OG_clans) <- paste(substr(names(OG_clans), 1, 9), paste(clans.num, '.', sep=''), sep='_')
head(names(OG_clans))
# total clans
length(OG_clans) 
table(sapply(OG_clans, function(i) sum(c('Ssal', 'Omyk') %in% substr(i$tip.label, 1, 4))))




# fix O.mykiss names in trees from proteinID to geneID
# taken from this dataframe
name.tab = read.table(paste('C:/users/meeldurb/Google Drive/', 
                            'Master internship phylogenetics salmonids/',
                            'Salmonid_genomics_resources/Orthologs_homeologs/', 
                            'orthogroups.03.06.2017/Omyk.gene2protein.table.txt', 
                            sep = ''), header = T)

OG_clans = lapply(OG_clans, function(i){
  tr = i
  tips = i$tip.label
  tips[!is.na(match(sub('Omyk\\|', '', tips), name.tab$protein))] <- paste('Omyk|', name.tab$gene_id[na.omit(match(sub('Omyk\\|', '', tips), name.tab$protein))], sep='')
  tr$tip.label <- tips
  tr
}
)

## ==> change the clan names of single clan orthogroups to clan_0
#Note - names of clans where the OG represented a single clans: .0

# is there redundancy in OG_clans??
orthogrtips = unlist(sapply(trees, function(i) i$tip.label))
table(duplicated(orthogrtips))

clantips = unlist(sapply(OG_clans, function(i) i$tip.label))
table(duplicated(clantips)) 
#NOTE: changing Omyk names generates 128 redundant sequences in clans...
#==> Gareth - figure out what Gene.ID are redundant..


#-------------------------------#
##_____ Make lookup table _____##
#-------------------------------#

tips = sapply(OG_clans, function(i) i$tip.label)
#names(tips) <- paste(names(tips), '_', sep='')
head(tips)
tips = unlist(tips)
head(tips)
lookup = data.frame(OG=substr(names(tips), 1, 9), 
                    clan = sub('\\..*', '', names(tips)),
                    gene = substr(tips, 6, 100),
                    species = substr(tips, 1, 4),
                    tip = tips, stringsAsFactors = F)
rownames(lookup) <- NULL
head(lookup, 50)

## identify clans with same Omyk tipname and remove the suckers...
clans_with_dupOmykGenes = unique(paste(lookup$clan[lookup$tip %in% 
                        lookup$tip[which(duplicated(lookup$tip))]], '.', sep=''))
length(OG_clans)-length(clans_with_dupOmykGenes)
keep.names= names(OG_clans)[!names(OG_clans) %in% clans_with_dupOmykGenes]
length(keep.names)
OG_clans_filt = OG_clans[keep.names]
length(OG_clans_filt)

##---------------------------##
#________ filter data ________#
##---------------------------##

# save(OG_clans_filt, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
#                                  'Beast_dating_salmonids/RData/',
#                                  'Ortho_clans_filtered.RData',
#                                  sep = ''))

OG_clans_filt <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                           'Beast_dating_salmonids/RData/',
                           'Ortho_clans_filtered.RData',
                           sep = ''))


# load location of fasta files
fasta.files <- paste('C:/Users/meeldurb/Google Drive/',
                     'Master internship phylogenetics salmonids/',
                     'Salmonid_genomics_resources/Orthologs_homeologs/',
                     'orthogroups.03.06.2017/Alignments/', sep = '')

# read the alignments in memory
alignments = lapply(dir(fasta.files, full.names = T), function(i){
  t = try(read.fasta(i), silent = T)
  if(class(t)=='try-error') return(NULL)
  if(class(t)!='try-error') return(t)
  }
)


# change names of the alignments. Remove .fa
names(alignments) <- sub('fa', '', dir(fasta.files))

# extract the duplicated alignments
table(duplicated(sub('fa', '', dir(fasta.files))))

# how much of the alignments are not empty
table(sapply(alignments, is.null))[1]/sum(table(sapply(alignments, is.null)))

# filter out the empty alignments
alignments_filtered = alignments[!sapply(alignments, is.null)]
length(alignments)
length(alignments_filtered)
names(alignments_filtered)

length(alignments_filtered$OG0018653.)
length(OG_clans_filt$OG0018653.$tip.label)

# check out some of the trees
par(mfrow=c(3,1))
plot(OG_trees[["OG0009052."]])
plot(OG_clans_filt$OG0009052_1.)


# not the whole tree is splitted up in clans, the mammal outgroup is removed
length(OG_trees[["OG0001162."]]$tip.label)
length(OG_clans[["OG0001162_1."]]$tip.label)
length(OG_clans[["OG0001162_2."]]$tip.label)



##---------------------------##
#_____ Find proper trees _____#
##---------------------------##

# make a dataframe of all the alignments with the names
alignments.df = ldply(alignments_filtered, function(i) {
  data.frame(tips=names(i), seqs=as.character(sapply(i, paste, collapse='')
                                              ), stringsAsFactors = F)
  })

# removing the double organism name
alignments.df$tips <- substr(alignments.df$tips, 6, 100)


# save(alignments.df, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
#                                  '/Beast_dating_salmonids/RData/',
#                                  'Clan_alignments_df_aa.RData', sep = ''))

alignments.df <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                     '/Beast_dating_salmonids/RData/',
                     'Clan_alignments_df_aa.RData', sep = ''))

head(alignments.df)
dim(alignments.df)

# find which OG clans have at least 7 tips in the alignment dataframe
idx.toanalyze = c()
for(i in 1:length(OG_clans_filt)){
  idx.toanalyze[i] <- sum(!is.na(match(OG_clans_filt[[i]]$tip.label, 
                                       alignments.df$tips)))>7
}
table(idx.toanalyze)

# only select the OG clans which have at least 7 tips that are also
# contained in the alignments dataframe
OG_clans_filt_2analyze = OG_clans_filt[idx.toanalyze]
length(OG_clans_filt_2analyze)

save(OG_clans_filt_2analyze, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
                     '/Beast_dating_salmonids/RData/',
                    'Clans_2analyze_inBeast_aa.RData', sep = ''))



##---------------------------##
#____ check for duplicates ___#
##---------------------------##

OG_clans_filt_2analyze <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                          '/Beast_dating_salmonids/RData/',
                                          'Clans_2analyze_inBeast_aa.RData', sep = ''))

# use dup tables below and check that clans contain at least one dup pair of one of the species...
#/Users/srsand/Google\ Drive/Salmonid_genomics_resources/Orthologs_homeologs/Homeologs/OmykV6_2016_best_in_homelogRegions_minpident80_mincov50_phylofiltered.RData
#/Users/srsand/Google\ Drive/Salmonid_genomics_resources/Orthologs_homeologs/Homeologs/RefSeq_GarethLongest_2016_best_in_homelogRegions_minpident80_mincov50_phylofiltered.RData 


# check if OG clans contain duplicated pair of Omyk and Ssal
Omyk.duplicates <- loadRData(paste('C:/Users/meeldurb/Google Drive/',
                            'Master internship phylogenetics salmonids/',
                            'Salmonid_genomics_resources/Orthologs_homeologs/',
                            'Homeologs/OmykV6_2016_best_in_homelogRegions_',
                            'minpident80_mincov50_phylofiltered.RData', sep = ''))[,1:2]

Ssal.duplicates <- loadRData(paste('C:/Users/meeldurb/Google Drive/',
                            'Master internship phylogenetics salmonids/',
                            'Salmonid_genomics_resources/Orthologs_homeologs/',
                            'Homeologs/RefSeq_GarethLongest_2016_best_',
                            'in_homelogRegions_minpident80_mincov50_',
                            'phylofiltered.RData', sep = ''))[,1:2]

# add organism naming in front on duplicate pairs
Omyk.dup <- as.data.frame(apply(Omyk.duplicates, 2, function(i) paste('Omyk|', i, sep='')))
Ssal.dup <- as.data.frame(apply(Ssal.duplicates, 2, function(i) paste('Ssal|', i, sep='')))



# find how many have at least 2 of Ssal or Omyk genes 
# this should resemble the duplicate check
table(sapply(OG_clans_filt_2analyze, function(i) {
  length(which(substr(i$tip.label, 1,4) == "Ssal")) >= 2 |
    length(which(substr(i$tip.label, 1,4) == "Omyk")) >= 2 }))



# check if there is at least 1 duplicate pair
# for each row in the dup.table we want to check whether
# they are contained in the tip.label
idx.dup = 0
for (idx in 1:length(OG_clans_filt_2analyze)){
  idx.dup[idx] <- sum(sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, Omyk.dup[,1]))),  
                      sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, 
                                       Omyk.dup[,2])))) >= 2 |
    sum(sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, Ssal.dup[,1]))),  
        sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, 
                         Ssal.dup[,2])))) >= 2
}

table(idx.dup)
# keep the OG groups that were TRUE
keep = which(idx.dup == 1)
length(OG_clans_filt_2analyze)
OG_clans_dupl = OG_clans_filt_2analyze[keep]
length(OG_clans_dupl)


##---------------------------------------------------##
#______ Checking if duplicate check is correct _______#
##---------------------------------------------------##

# first Omyk
idx.dup = 0
for (idx in 1:length(OG_clans_filt_2analyze)){
  idx.dup[idx] <- sum(sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, Omyk.dup[,1]))),  
        sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, 
                         Omyk.dup[,2])))) >= 2
}

table(idx.dup)
# keep the OG groups that were TRUE
keep = which(idx.dup == 1)
length(OG_clans_filt_2analyze)
OG_clans_dupl.omyk = OG_clans_filt_2analyze[keep]
length(OG_clans_dupl.omyk)

# then Ssal
idx.dup = 0
for (idx in 1:length(OG_clans_filt_2analyze)){
  idx.dup[idx] <- sum(sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, Ssal.dup[,1]))),  
                      sum(!is.na(match(OG_clans_filt_2analyze[[idx]]$tip.label, 
                                       Ssal.dup[,2])))) >= 2
}

table(idx.dup)
# keep the OG groups that were TRUE
keep = which(idx.dup == 1)
length(OG_clans_filt_2analyze)

OG_clans_dupl.ssal = OG_clans_filt_2analyze[keep]
length(OG_clans_dupl.ssal)




fulldf.check <- unique(names(c(OG_clans_dupl.omyk, OG_clans_dupl.ssal)))
length(fulldf.check)

##---------------------------##
#______ Check for pike _______#
##---------------------------##

# We also want Eluc in the trees
eluc.tips = sapply(OG_clans_dupl, function(i) sum(substr(i$tip.label, 1, 4) %in% 'Eluc'))

table(eluc.tips>0)

OG_clans_dupl_Eluc_filt = OG_clans_dupl[which(eluc.tips>0)]


##---------------------------##
#______ Write out data _______#
##---------------------------##



save(OG_clans_dupl_Eluc_filt, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                  '/Beast_dating_salmonids/RData/',
                                  'Clans_2analyze_inBeast_withduplicatesandElucfilt_aa.RData', 
                                  sep = ''))
OG_clans_dupl <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                 '/Beast_dating_salmonids/RData/',
                                 'Clans_2analyze_inBeast_withduplicatesandElucfilt_aa.RData', 
                                 sep = ''))
