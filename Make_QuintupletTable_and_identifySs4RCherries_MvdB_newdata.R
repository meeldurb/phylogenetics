#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch & Simen Rod Sandve
##Pipeline for quality control and x-checking trees 
##with dup-predictions
###################################################################


#--------------------------------------#
##____ load libraries & functions ____##
#--------------------------------------#

library(ape)
library(phangorn)
library(plyr)
library(RCurl)


# source functions

# clanfinder function
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "srsand/Phylogenomics/master/clanfinder.R", 
                               sep = ""), ssl.verifypeer = FALSE)))

# phylo functions
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "srsand/Phylogenomics/master/Phylo_functions.R", 
                               sep = ""), ssl.verifypeer = FALSE)))

# autoroot function
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "srsand/Phylogenomics/master/auto.root_salmonid_clans.R", 
                               sep = ""), ssl.verifypeer = FALSE)))

# loadRData in variable function
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "meeldurb/phylogenetics/master/loadRData.R", 
                               sep = ""), ssl.verifypeer = FALSE)))



## NB: clanfinder_v2 er buggy! her faar vi ut masse dritt som er redundant


#---------------------#
##____ Load data ____##
#---------------------#

# 1) load OG clans

# filtered_OG_clans <- loadRData(paste('/mnt/users/melavan/RData/',
#                                  '20170913-Clans_forBEAST_dupssandElucfilt_nt.RData', 
#                                  sep = ''))


filtered_OG_clans <- loadRData(paste('C:/Users/meeldurb/Dropbox/',
                                    'Melanie/Master_internship_',
                                    'phylogenetics/phylogenetics/RData/',
                                    '20170913-Clans_forBEAST_dupssandElucfilt_nt.RData', 
                                    sep = ''))



##-------------------------------------##
## ______ Filter for duplicates ______ ##
##-------------------------------------##

# omyk.d <- read.table(paste('/mnt/users/melavan/duplicate_tables/',
#                           'Omyk_RefSeq_2017_Ss4R_ohnologs_',
#                           'blastidentification_minpident80.txt', sep = ''),
#                      header = T)[1:2]
# 
# ssal.d <-  read.table(paste('/mnt/users/melavan/duplicate_tables/',
#                             'Ssal_RefSeq_2017_Ss4R_ohnologs_',
#                             'blastidentification_minpident80.txt', sep = ''),
#                       header = T)[1:2]

omyk.d <- read.table(paste('C:/Users/meeldurb/Google Drive/',
                           'Master internship phylogenetics salmonids/',
                           'Salmonid_genomics_resources/Orthologs_homeologs/',
                           'Homeologs/Omyk_RefSeq_2017_Ss4R_ohnologs_',
                           'blastidentification_minpident80.txt', sep = ''),
                     header = T, stringsAsFactors = F)[1:2]


ssal.d <- read.table(paste('C:/Users/meeldurb/Google Drive/',
                           'Master internship phylogenetics salmonids/',
                           'Salmonid_genomics_resources/Orthologs_homeologs/',
                           'Homeologs/Ssal_RefSeq_2017_Ss4R_ohnologs_',
                           'blastidentification_minpident80.txt', sep = ''),
                     header = T, stringsAsFactors = F)[1:2]
  
  



#i <- filtered_OG_clans$OG0000019_2.

# find the number of duplicates for Ssal
dups.inclans.ssal = sapply(filtered_OG_clans, function(i) {
  tips = substr(i$tip.label, 6, 100) 
  qs = na.omit(match(tips, ssal.d$qseqid))
  ss = na.omit(match(tips, ssal.d$sseqid))
  if(length(qs)==0 | length(ss)==0 ) {
    return(NA)
  } else {
    qs.gene = ssal.d$qseqid[qs]
    ss.gene = ssal.d$sseqid[ss]
    length(c(qs.gene, ss.gene))
    }
  }
)

table(dups.inclans.ssal, useNA = 'always')
#names(which(dups.inclans.ssal >= 2))


#i <- filtered_OG_clans$OG0000019_2.

# find the number of duplicates for Omyk
dups.inclans.omyk = sapply(filtered_OG_clans, function(i) {
  tips = substr(i$tip.label, 6, 100) 
  qs = na.omit(match(tips, omyk.d$qseqid))
  ss = na.omit(match(tips, omyk.d$sseqid))
  if(length(qs)==0 | length(ss)==0 ) {
    return(NA)
  } else{
    qs.gene = omyk.d$qseqid[qs]
    ss.gene = omyk.d$sseqid[ss]
    length(c(qs.gene, ss.gene))
    }
  }
)

table(dups.inclans.omyk, useNA = 'always')
#names(which(dups.inclans.omyk >= 2))

#dups.inclans.omyk
#    2    3 <NA> 
#  950    1  977  

#dups.inclans.ssal
#     2    3 <NA> 
#  1744    1  183 


table(dups.inclans.omyk==2 & dups.inclans.ssal==2, useNA = 'always')
OG_clans_dups = filtered_OG_clans[which(dups.inclans.omyk==2 & dups.inclans.ssal==2)]

length(filtered_OG_clans)
length(OG_clans_dups)


# save(OG_clans_dups, file = paste('/mnt/users/melavan/RData/',
#                                  '20170913-Clans_forBEAST_',
#                                 'duplicates_final_nt.RData ', sep = ''))
# 
# OG_clans_dups <- loadRData(paste('/mnt/users/melavan/RData/',
#                                  '20170913-Clans_forBEAST_',
#                                  'duplicates_final_nt.RData ', sep = ''))


save(OG_clans_dups, file = paste('C:/Users/meeldurb/Dropbox/',
                                 'Melanie/Master_internship_',
                                 'phylogenetics/phylogenetics/RData/',
                                 '20170913-Clans_forBEAST_',
                                 'duplicates_final_nt.RData ', sep = ''))


OG_clans_dups <- loadRData(paste('C:/Users/meeldurb/Dropbox/',
                                 'Melanie/Master_internship_',
                                 'phylogenetics/phylogenetics/RData/',
                                 '20170913-Clans_forBEAST_',
                                 'duplicates_final_nt.RData ', sep = ''))


##------------------------------------##
## Making tables with quintuplets:    ##
## RT dups + Sal dups + pike ortholog ##
##------------------------------------##

i <- OG_clans_dups$OG0008390_1.

dup.table = lapply(OG_clans_dups, function(i) {
  # find the positions of the duplicates
  ssal.qs = na.omit(match(substr(i$tip.label, 6, 100), ssal.d$qseqid))
  ssal.ss = na.omit(match(substr(i$tip.label, 6, 100), ssal.d$sseqid))
  omyk.qs = na.omit(match(substr(i$tip.label, 6, 100), omyk.d$qseqid))
  omyk.ss = na.omit(match(substr(i$tip.label, 6, 100), omyk.d$sseqid))
  
  # take the Eluc gene
  eluc = grep('Eluc', i$tip.label, value = T)
  
  # check if duplicates are found in same position in duplicate tables
  # meaning they are found on the same line, so are true duplicates
  # and check if only one Eluc found
  if(ssal.qs==ssal.ss & omyk.qs==omyk.ss & length(eluc)==1) { 
    return(c(eluc, ssal.d$qseqid[ssal.qs], ssal.d$sseqid[ssal.qs],
             omyk.d$qseqid[omyk.ss], omyk.d$sseqid[omyk.ss]))
    }
  if(ssal.qs!=ssal.ss | omyk.qs!=omyk.ss | length(eluc)>1) {
    return('PROBLEM')
  }
})
names(dup.table) <- names(OG_clans_dups)
table(sapply(dup.table, function(i) length(grep('PROBLEM', i)))) # clans we can use

dup.table = dup.table[sapply(dup.table, function(i) length(grep('PROBLEM', i)))==0]
length(dup.table)


# make a dataframe 
dup.table = ldply(dup.table, function(i) {
  t(data.frame(i, stringsAsFactors = F))
})

save(dup.table, file = paste('C:/Users/meeldurb/Dropbox/',
                             'Melanie/Master_internship_',
                             'phylogenetics/phylogenetics/RData/',
                             '20170913_duplicate_table_aa.RData', sep = ''))


dup.table <- loadRData(paste('C:/Users/meeldurb/Dropbox/',
                             'Melanie/Master_internship_',
                             'phylogenetics/phylogenetics/RData/',
                             '20170913_duplicate_table_aa.RData', sep = ''))

##--------------------------------------------------##
##       Identify Ss4R-cherries and orthologs       ##
##--------------------------------------------------##



## for loop to check consistency between dup table and clan topology
# setting some objects

topology.results = c() 
lore.results = c()
lore.taxa.results = c()

dup_cluster_df = matrix(rep(NA), nrow = nrow(dup.table), ncol=4)

i <- 714
for(i in 1:nrow(dup.table)){
  #print(i)
  test.tree <- OG_clans_dups[[dup.table$.id[i]]]
  test.tree <- auto.root(test.tree)$rooted.clans # Alternative 1
  #test.tree <- unroot(test.tree) # Alternative 2 ==> gives 2158 good and 2985 inconsisten topologies
  #print(test.tree$tip.label)
  #plot(test.tree); nodelabels(); tiplabels()
  
  # get MRCA of 
  mrca.combo1 = getMRCA(test.tree, tip = c(paste('Ssal|', dup.table[i,3], sep=''), 
                                           paste('Omyk|', dup.table[i,5], sep='')))
  mrca.combo2 = getMRCA(test.tree, tip = c(paste('Ssal|', dup.table[i,3], sep=''), 
                                           paste('Omyk|', dup.table[i,6], sep='')))
  test.cherry1 = Descendants(test.tree, mrca.combo1, type='tips')[[1]]
  test.cherry2 = Descendants(test.tree, mrca.combo2, type='tips')[[1]]
  
  # select which combo of ssal and omyk duplicates
  # gives a mrca with smallest descedants clade size
   if(length(test.cherry1) < length(test.cherry2)) { 
    clades = list(clade1=c(3,5), clade2=c(4,6))
    } else {
      clades = list(clade1=c(3,6), clade2=c(4,5))
    }
  
  # combine the correct duplicates in a dataframe
  dup_cluster_df[i, 1] <- paste('Ssal|', dup.table[i, clades[[1]][1]], sep='')
  dup_cluster_df[i, 2] <- paste('Omyk|', dup.table[i, clades[[1]][2]], sep='')
  dup_cluster_df[i, 3] <- paste('Ssal|', dup.table[i, clades[[2]][1]], sep='')
  dup_cluster_df[i, 4] <- paste('Omyk|', dup.table[i, clades[[2]][2]], sep='')
  
  
  # testing if both cherry clades only contains 1 ssal and 1 omyk duplicate
  cherry1 = test.tree$tip.label[unlist(Descendants(test.tree, 
            getMRCA(test.tree, c(paste('Ssal|', dup.table[i, clades[[1]][1]], sep=''),
            paste('Omyk|', dup.table[i, clades[[1]][2]], sep=''))), type='tips'))]
  
  cherry2 = test.tree$tip.label[unlist(Descendants(test.tree, 
            getMRCA(test.tree, c(paste('Ssal|', dup.table[i, clades[[2]][1]], sep=''),
            paste('Omyk|', dup.table[i, clades[[2]][2]], sep=''))), type='tips'))]
  
  if(any(cherry1 %in% cherry2)){
    topology.results[i] <- 'Topology Incongruence'
  }
  if(!any(cherry1 %in% cherry2)){
    topology.results[i] <- 'Good Topology'
  }
  
  # lore.results[i] <- sum(sapply(c('Tthy', 'Ssal', 'Omyk', 'Okis'), function(i){
  #   ifelse(length(grep(i,  test.tree$tip.label))>0, 
  #          is.monophyletic(test.tree, grep(i, test.tree$tip.label)), NA)
  # } 
  # ), na.rm = T)
  # 
  # lore.taxa.names = sapply(c('Tthy', 'Ssal', 'Omyk', 'Okis'), function(i){
  #   ifelse(length(grep(i,  test.tree$tip.label))>0,  
  #          is.monophyletic(test.tree, grep(i, test.tree$tip.label)), NA)
  # } 
  # )
  # 
  # lore.taxa.results[i] <- paste(sort(names(lore.taxa.names)[lore.taxa.names]), collapse='|')
  
}


#MelB: lore.taxa.names/results are not really relevant for you at this stage 

colnames(dup_cluster_df) <- c('Ssal1', 'Omyk1', 'Ssal2', 'Omyk2')
dup_cluster_fix = cbind(dup.table[,1:2], dup_cluster_df)
dup_cluster_phylofilt = dup_cluster_fix[topology.results!='Topology Incongruence', ]

# final quintuplettable sorted based on tree topology information 
# # (ie ortholog/ohnolog relationships)
# 
# save(dup_cluster_phylofilt, file = paste('/mnt/users/melavan/RData',
#                                          '20170913_duplicate_clans_filtered_final_nt.RData', 
#                                          sep = ''))
# 
# dup_cluster_phylofilt <- loadRData(paste('/mnt/users/melavan/RData',
#                                          '20170913_duplicate_clans_filtered_final_nt.RData', 
#                                          sep = ''))
# 
# write.table(dup_cluster_phylofilt, file = paste('/mnt/users/melavan/RData',
#                                                 '20170913_duplicate_clans_filtered_final_nt.csv',
#                                                 sep = ''), 
#             append = FALSE, sep = ';', quote = FALSE, row.names = FALSE, col.names = TRUE)


save(dup_cluster_phylofilt, file = paste('C:/Users/meeldurb/Dropbox/',
                                         'Melanie/Master_internship_',
                                         'phylogenetics/phylogenetics/RData/',
                                         '20170913_duplicate_clans_filtered_final_nt.RData', 
                                         sep = ''))

dup_cluster_phylofilt <- loadRData(paste('C:/Users/meeldurb/Dropbox/',
                                         'Melanie/Master_internship_',
                                         'phylogenetics/phylogenetics/RData/',
                                         '20170913_duplicate_clans_filtered_final_nt.RData', 
                                         sep = ''))

write.table(dup_cluster_phylofilt, file = paste('C:/Users/meeldurb/Dropbox/',
                                                'Melanie/Master_internship_',
                                                'phylogenetics/phylogenetics/RData/',
                                                '20170913_duplicate_clans_filtered_final_nt.csv',
                                                sep = ''), 
            append = FALSE, sep = ';', quote = FALSE, row.names = FALSE, col.names = TRUE)

print(paste("Number of OG_clans at start: ", length(filtered_OG_clans), sep = ""))
print(paste("OG clans number after duplicates filtering: ", length(OG_clans_dups), sep = ""))
print(paste("OG clans number after linking to duplicate pairs: ", 
            length(dup_cluster_phylofilt[,1]), sep = ""))




