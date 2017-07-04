#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch & Simen Rod Sandve
##Script for loading and manipulating the orthoclan information
###################################################################


#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

#install.packages('RCurl',  repos = "http://cran.rstudio.com/")
#install.packages("ape", repos = "http://cran.rstudio.com/")
#install.packages("phangorn", repos = "http://cran.rstudio.com/")

library(RCurl)
library(ape)
library(phangorn)
library(seqinr)



source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenomics/clanfinder.R', sep = ''))
source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenomics/Phylo_functions.R', sep = ''))

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#-------------------------------------------------#
##____ make initial clans from protein trees ____##
#-------------------------------------------------#

# 1) get clans from orthogroups

# load orthotrees - this is the output from protein trees made 
# from orthogroup sets of sequences (from orthofinder)
OG_trees <- loadRData(paste('C:/Users/meeldurb/Google Drive/Master internship ',
                    'phylogenetics salmonids/Salmonid_genomics_resources/',
                    'Orthologs_homeologs/orthogroups.03.06.2017/',
                    'OG_trees.30.05.2017.RData', sep = ''))

# total number of orthogroups
length(OG_trees) 
# summing up number of specific species in trees.
table(sapply(OG_trees, function(i) sum(c('Ssal', 'Omyk') %in% substr(i$tip.label, 1, 4))))

# 2) make clans using clanfinder

OG_clanfinder = lapply(OG_trees, clanFinder, ut = c('Olat', 'Gacu', 'Drer', 'Locu', 'Mmus', 'Hsap'))

# Play a bit with the datastructure
names(OG_trees)
# the names of the trees in the different OG groups are named by a integer
names(OG_clanfinder)
names(OG_clanfinder[1])
names(OG_clanfinder[[1]])
names(OG_clanfinder$OG0000000.)
OG_clanfinder[[1]][1]

# 3) fixing names

# some of the OG clans are empty, we need to remove these
OG_clanfinder.filt <- OG_clanfinder[sapply(OG_clanfinder, length)>0]
length(OG_clanfinder)
length(OG_clanfinder.filt)


clans.num <- as.numeric(unlist(sapply(OG_clanfinder.filt, function(i) 1:length(i))))
OG_clans <- unlist(OG_clanfinder.filt, recursive = F)
head(names(OG_clans))
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
orthogrtips = unlist(sapply(OG_trees, function(i) i$tip.label))
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
#_______ write out data ______#
##---------------------------##

#save(OG_clans_filt, file = '~/Dropbox/Work/PhDs_and_Master/ERASMUS/Melanie/Beast_dating_salmonids/RData/Ortho_clans_filtered.RData')

OG_clans_filt <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                           'Beast_dating_salmonids/RData/',
                           'Ortho_clans_filtered.RData',
                           sep = ''))


# load location of fasta files
fasta.files <- paste('C:/Users/meeldurb/Google Drive/',
                     'Master internship phylogenetics salmonids/',
                     'Salmonid_genomics_resources/Orthologs_homeologs/',
                     'orthogroups.03.06.2017/cds_pal2nal/', sep = '')

alignments = lapply(dir(fasta.files, full.names = T), function(i){
  t = try(read.fasta(i), silent = T)
  if(class(t)=='try-error') return(NULL)
  if(class(t)!='try-error') return(t)
  }
)
names(alignments) <- sub('fa', '', dir(fasta.files))
table(duplicated(sub('fa', '', dir(fasta.files))))

table(sapply(alignments, is.null))[1]/sum(table(sapply(alignments, is.null)))

alignments_filtered = alignments[!sapply(alignments, is.null)]
length(alignments_filtered)
names(alignments_filtered)

length(alignments_filtered[["OG0001162."]])
length(OG_clans[["OG0001162_1."]]$tip.label)
length(OG_clans[["OG0001162_2."]]$tip.label)

par(mfrow=c(3,1))
plot(OG_trees[["OG0001162."]])
plot(OG_clans[["OG0001162_1."]])
plot(OG_clans[["OG0001162_2."]])

length(OG_trees[["OG0001162."]]$tip.label)
length(OG_clans[["OG0001162_1."]]$tip.label)
length(OG_clans[["OG0001162_2."]]$tip.label)



##:
library(plyr)
alignments.df = ldply(alignments_filtered,function(i) data.frame(tips=names(i), seqs=as.character(sapply(i, paste, collapse='')), stringsAsFactors = F))
alignments.df$tips <- substr(alignments.df$tips, 6, 100)

save(alignments.df, file = '~/Dropbox/Work/PhDs_and_Master/ERASMUS/Melanie/Beast_dating_salmonids/RData/Clan_alignments_df.RData')
head(alignments.df)
dim(alignments.df)

idx.toanalyze = c()
for(i in 1:length(OG_clans_filt)){
  idx.toanalyze[i] <- sum(!is.na(match(OG_clans_filt[[i]]$tip.label, alignments.df$tips)))>7
}
table(idx.toanalyze)

OG_clans_filt_2analyze = OG_clans_filt[idx.toanalyze]
save(OG_clans_filt_2analyze, file = '~/Dropbox/Work/PhDs_and_Master/ERASMUS/Melanie/Beast_dating_salmonids/RData/Clans_2analyze_inBeast.RData')

# use dup tables below and check that clans contain at least one dup pair of one of the species...
#/Users/srsand/Google\ Drive/Salmonid_genomics_resources/Orthologs_homeologs/Homeologs/OmykV6_2016_best_in_homelogRegions_minpident80_mincov50_phylofiltered.RData
#/Users/srsand/Google\ Drive/Salmonid_genomics_resources/Orthologs_homeologs/Homeologs/RefSeq_GarethLongest_2016_best_in_homelogRegions_minpident80_mincov50_phylofiltered.RData 




