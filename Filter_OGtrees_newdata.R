#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch & Simen Rod Sandve
##Script for loading and manipulating the orthoclan results
##and preparing them to be imported in the BEAUTI .xml files
###################################################################


#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

install.packages('RCurl',  repos = "http://cran.rstudio.com/")
install.packages("ape", repos = "http://cran.rstudio.com/")
install.packages("phangorn", repos = "http://cran.rstudio.com/")
install.packages("seqinr", repos = "http://cran.rstudio.com/")
install.packages("plyr", repos = "http://cran.rstudio.com/")


library(RCurl)
library(ape)
library(phangorn)
library(seqinr)
library(plyr)


# source functions


# clanfinder function
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "srsand/Phylogenomics/master/clanfinder.R", 
                               sep = ""), ssl.verifypeer = FALSE)))

# phylo functions
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "srsand/Phylogenomics/master/Phylo_functions.R", 
                               sep = ""), ssl.verifypeer = FALSE)))

# loadRData in variable function
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "meeldurb/phylogenetics/master/loadRData.R", 
                               sep = ""), ssl.verifypeer = FALSE)))



#------------------------#
# ____ Load clans ______ #
#------------------------#

# loading files after running clanFinder

# file with the '.' behind clanID
OG_clans <- loadRData(paste('C:/Users/meeldurb/Google Drive/Master internship ',
                            'phylogenetics salmonids/Salmonid_genomics_resources/',
                            'Orthologs_homeologs/orthogroups_2017/',
                            'OG_trees+clans.06.08.17.RData', sep = ''))

# file without "." behind clanID
# OG_clans <- loadRData(paste('C:/Users/meeldurb/Google Drive/Master internship ',
#                                  'phylogenetics salmonids/Salmonid_genomics_resources/',
#                                  'Orthologs_homeologs/orthogroups_2017/',
#                                  'OG_clan_trees.06.08.17.RData', sep = ''))

# OG_clantable <- loadRData(paste('C:/Users/meeldurb/Google Drive/Master internship ',
#                              'phylogenetics salmonids/Salmonid_genomics_resources/',
#                              'Orthologs_homeologs/orthogroups_2017/',
#                              'OG_clan_table.06.08.17.RData', sep = ''))




# is there redundancy in OG_clans??
orthogrtips = unlist(sapply(trees, function(i) i$tip.label))
table(duplicated(orthogrtips))

clantips = unlist(sapply(OG_clans, function(i) i$tip.label))
table(duplicated(clantips)) 


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


##---------------------------##
#________ filter data ________#
##---------------------------##

# save(OG_clans_filt, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
#                                  'Beast_dating_salmonids/RData/',
#                                  'Ortho_clans_filtered.RData',
#                                  sep = ''))

# OG_clans_filt <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
#                                  'Beast_dating_salmonids/RData/',
#                                  'Ortho_clans_filtered.RData',
#                                  sep = ''))


# load location of fasta files
fasta.files <- paste('C:/Users/meeldurb/Google Drive/',
                     'Master internship phylogenetics salmonids/',
                     'Salmonid_genomics_resources/Orthologs_homeologs/',
                     'orthogroups_2017/cds_macse_nt_align/', sep = '')

# read the alignments in memory
# alignments <- NULL
# for (alignment in dir(fasta.files, full.names = T)){
#   print (alignment)
#   alignments <- c(alignments, read.fasta(alignment))
# }

# alignments = sapply(dir(fasta.files, full.names = T), function(i){
#   t = tryCatch(read.fasta(i), error = function(e) NULL)
#   })

#i <- dir(fasta.files, full.names = T)[1]

# read the alignments in memory
#fail.report <- NULL
alignments = lapply(head(dir(fasta.files, full.names = T), n = 500), function(i){
  t = read.fasta(i) })
  # if(class(t)=='try-error') {
  #   print(i)
  #   #print(substr(i, 158, 169), "has an error")
  #   #fail.report <- c(fail.report, paste(substr(i, 158, 169), "has an error", sep = ""))
  #   return(NULL)
  # }
  # if(class(t)!='try-error') {
  #   #print(substr(i, 158, 169), "was successfull")
  #   #fail.report <- c(fail.report, paste(substr(i, 158, 169), "was successfull", sep = ""))
  #   return(t)
  # }
  # })

# change names of the alignments. Remove .aln
names(alignments) <- sub('aln', '', head(dir(fasta.files), n = 500))
head(alignments)

# extract the duplicated alignments
table(duplicated(sub('aln', '', dir(fasta.files))))

# how much of the alignments are not empty
table(sapply(alignments, is.null))[1]/sum(table(sapply(alignments, is.null)))

# filter out the empty alignments
alignments_filtered = alignments[!sapply(alignments, is.null)]
length(alignments)
length(alignments_filtered)

save(alignments.df, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                 '/Beast_dating_salmonids/RData/',
                                 '20170913-Clan_aligments_nt.RData', sep = ''))

##---------------------------##
#_____ Find proper trees _____#
##---------------------------##

# make a dataframe of all the alignments with the names
alignments.df = ldply(alignments_filtered, function(i) {
  data.frame(tips=names(i), seqs=as.character(sapply(i, paste, collapse='')
  ), stringsAsFactors = F)
})


save(alignments.df, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                 '/Beast_dating_salmonids/RData/',
                                  '20170913-Clan_alignments_df_nt.RData', sep = ''))

alignments.df <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                 '/Beast_dating_salmonids/RData/',
                                 '20170913-Clan_alignments_df_nt.RData', sep = ''))

head(alignments.df)
dim(alignments.df)

# find which OG clans have at least 7 tips in the alignment dataframe
idx.toanalyze = c()
for(i in 1:length(OG_clans)){
  idx.toanalyze[i] <- sum(!is.na(match(OG_clans[[i]]$tip.label, 
                                       alignments.df$tips)))>7
}
table(idx.toanalyze)

# only select the OG clans which have at least 7 tips that are also
# contained in the alignments dataframe
OG_clans_filt = OG_clans[idx.toanalyze]
length(OG_clans_filt)

save(OG_clans_filt, file = paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                  '/Beast_dating_salmonids/RData/',
                                  '20170913-OG_Clans_forBEAST_nt.RData', sep = ''))



##---------------------------##
#____ check for duplicates ___#
##---------------------------##

OG_clans_filt <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                  '/Beast_dating_salmonids/RData/',
                                  '20170913-OG_Clans_forBEAST_nt.RData', sep = ''))

# use dup tables below and check that clans contain 
# at least one dup pair of one of the species...
#/Users/srsand/Google\ Drive/Salmonid_genomics_resources/Orthologs_homeologs/Homeologs/OmykV6_2016_best_in_homelogRegions_minpident80_mincov50_phylofiltered.RData
#/Users/srsand/Google\ Drive/Salmonid_genomics_resources/Orthologs_homeologs/Homeologs/RefSeq_GarethLongest_2016_best_in_homelogRegions_minpident80_mincov50_phylofiltered.RData 


# check if OG clans contain duplicated pair of Omyk and Ssal
Omyk.duplicates <- read.table(paste('C:/Users/meeldurb/Google Drive/',
                                    'Master internship phylogenetics salmonids/',
                                    'Salmonid_genomics_resources/Orthologs_homeologs/',
                                    'Homeologs/Omyk_RefSeq_2017_Ss4R_ohnologs_',
                                    'blastidentification_minpident80.txt', sep = ''),
                              header = T)[,1:2]

Ssal.duplicates <- read.table(paste('C:/Users/meeldurb/Google Drive/',
                                    'Master internship phylogenetics salmonids/',
                                    'Salmonid_genomics_resources/Orthologs_homeologs/',
                                    'Homeologs/Ssal_RefSeq_2017_Ss4R_ohnologs_',
                                    'blastidentification_minpident80.txt', sep = ''),
                              header = T)[,1:2]

# add organism naming in front on duplicate pairs
Omyk.dup <- as.data.frame(apply(Omyk.duplicates, 2, function(i) paste('Omyk|', i, sep='')))
Ssal.dup <- as.data.frame(apply(Ssal.duplicates, 2, function(i) paste('Ssal|', i, sep='')))


# find how many have at least 2 of Ssal or Omyk genes 
# this should resemble the duplicate check
table(sapply(OG_clans_filt, function(i) {
  length(which(substr(i$tip.label, 1,4) == "Ssal")) >= 2 |
    length(which(substr(i$tip.label, 1,4) == "Omyk")) >= 2 }))



# check if there is at least 1 duplicate pair
# for each row in the dup.table we want to check whether
# they are contained in the tip.label
idx.dup = 0
for (idx in 1:length(OG_clans_filt)){
  idx.dup[idx] <- sum(sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, Omyk.dup[,1]))),  
                      sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, 
                                       Omyk.dup[,2])))) >= 2 |
    sum(sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, Ssal.dup[,1]))),  
        sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, 
                         Ssal.dup[,2])))) >= 2
}

table(idx.dup)
# keep the OG groups that were TRUE
keep = which(idx.dup == 1)
length(OG_clans_filt)
OG_clans_dupl = OG_clans_filt[keep]
length(OG_clans_dupl)


##---------------------------------------------------##
#______ Checking if duplicate check is correct _______#
##---------------------------------------------------##

# first Omyk
idx.dup = 0
for (idx in 1:length(OG_clans_filt)){
  idx.dup[idx] <- sum(sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, Omyk.dup[,1]))),  
                      sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, 
                                       Omyk.dup[,2])))) >= 2
}

table(idx.dup)
# keep the OG groups that were TRUE
keep = which(idx.dup == 1)
length(OG_clans_filt)
OG_clans_dupl.omyk = OG_clans_filt[keep]
length(OG_clans_dupl.omyk)

# then Ssal
idx.dup = 0
for (idx in 1:length(OG_clans_filt)){
  idx.dup[idx] <- sum(sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, Ssal.dup[,1]))),  
                      sum(!is.na(match(OG_clans_filt[[idx]]$tip.label, 
                                       Ssal.dup[,2])))) >= 2
}

table(idx.dup)
# keep the OG groups that were TRUE
keep = which(idx.dup == 1)
length(OG_clans_filt)

OG_clans_dupl.ssal = OG_clans_filt[keep]
length(OG_clans_dupl.ssal)




fulldf.check <- unique(names(c(OG_clans_dupl.omyk, OG_clans_dupl.ssal)))
length(fulldf.check)
length(OG_clans_dupl)

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
                                           'Master_internship_phylogenetics/phylogenetics/RData/',
                                           '20170913-Clans_forBEAST_dupssandElucfilt_nt.RData', 
                                           sep = ''))
OG_clans_dupl <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                 'Master_internship_phylogenetics/phylogenetics/RData/',
                                 '20170918-Clans_forBEAST_dupssandElucfilt_nt.RData', 
                                 sep = ''))


print(paste("Number of OG_clans at start:", length(OG_clans), sep = ""))
print(paste("Number of alignment files at start:", length(dir(fasta.files)), sep = ""))
print(paste("OG clans number after standard filtering:", length(OG_clans_filt), sep = ""))
print(paste("OG clans number after duplicate filtering:", 
            length(OG_clans_dupl), sep = ""))
print(paste("OG clans number after duplicate and Eluc filtering:", 
            length(OG_clans_dupl_Eluc_filt), sep = ""))


