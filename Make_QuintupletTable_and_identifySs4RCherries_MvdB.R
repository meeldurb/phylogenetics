# pipeline for quality control and x-checking trees with dup-predictions

#--------------------------------------#
##____ load libraries & functions ____##
#--------------------------------------#

library(ape)
library(phangorn)


source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenomics/clanfinder.R', sep = ''))
source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenomics/Phylo_functions.R', sep = ''))

## NB: clanfinder_v2 er buggy! her faar vi ut masse dritt som er redundant


#-------------------------------------------------#
##____ make initial clans from protein trees ____##
#-------------------------------------------------#

# 1) load orthotrees

# - this is the output from protein trees made 
# from orthogroup sets of sequences (from orthofinder)
trees <- loadRData(paste('C:/Users/meeldurb/Google Drive/Master internship ',
                         'phylogenetics salmonids/Salmonid_genomics_resources/',
                         'Orthologs_homeologs/orthogroups.03.06.2017/',
                         'OG_trees.30.05.2017.RData', sep = ''))
# total number of orthogroups
length(trees) 
# summing up the number of specific species in the trees
table(sapply(trees, function(i) sum(c('Ssal', 'Omyk') %in% substr(i$tip.label, 1, 4))))


# 2) make clans using clanfinder
OG_clanfinder = lapply(trees, clanFinder, ut = c('Olat', 'Gacu', 'Drer', 'Locu', 'Mmus', 'Hsap'))

# 3) fixing names

# some of the OG clans are empty, we need to remove these
OG_clanfinder.filt <- OG_clanfinder[sapply(OG_clanfinder, length)>0]
clans.num <- as.numeric(unlist(sapply(OG_clanfinder.filt, function(i) 1:length(i))))
OG_clans <- unlist(OG_clanfinder.filt, recursive = F)
names(OG_clans) <- paste(substr(names(OG_clans), 1, 9), paste(clans.num, '.', sep=''), sep='_')


# fix O.mykiss names in trees from proteinID to geneID
# taken from this dataframe
name.tab = read.table(paste('C:/users/meeldurb/Google Drive/', 
                            'Master internship phylogenetics salmonids/',
                            'Salmonid_genomics_resources/Orthologs_homeologs/', 
                            'orthogroups.03.06.2017/Omyk.gene2protein.table.txt', 
                            sep = ''), header = T)

table(duplicated(name.tab$protein))
table(duplicated(name.tab$gene_id))

OG_clans = lapply(OG_clans, function(i){
  tr = i
  tips = i$tip.label
  tips[!is.na(match(sub('Omyk\\|', '', tips), name.tab$protein))] <- paste('Omyk|', name.tab$gene_id[na.omit(match(sub('Omyk\\|', '', tips), name.tab$protein))], sep='')
  tr$tip.label <- tips
  tr
}
)


#--------------------------------------------------#
##____ filter clans based on species included ____##
#--------------------------------------------------#

salmonid.tips = sapply(OG_clans, function(i) sum(unique(substr(i$tip.label, 1, 4)) %in% sal))
eluc.tips = sapply(OG_clans, function(i) sum(substr(i$tip.label, 1, 4) %in% 'Eluc'))

# Filter clans: including minimum 2 salmonid species AND Eluc
table(salmonid.tips>=2 & eluc.tips>0)
OG_clans_filt = OG_clans[which(salmonid.tips>=2 & eluc.tips>0)]






##---------------------------------------------------------------------##
## making tables with quintuplets: RT dups + Sal dups + pike orhtolog ###
##---------------------------------------------------------------------##

ssal.d <- loadRData(paste('C:/Users/meeldurb/Google Drive/',
                'Master internship phylogenetics salmonids/',
                'Salmonid_genomics_resources/Orthologs_homeologs/',
                'Homeologs/OmykV6_2016_best_in_homelogRegions_',
                'minpident80_mincov50_phylofiltered.RData', sep = ''))

omyk.d <- loadRData(paste('C:/Users/meeldurb/Google Drive/',
                          'Master internship phylogenetics salmonids/',
                          'Salmonid_genomics_resources/Orthologs_homeologs/',
                          'Homeologs/RefSeq_GarethLongest_2016_best_',
                          'in_homelogRegions_minpident80_mincov50_',
                          'phylofiltered.RData', sep = ''))




dups.inclans.ssal = sapply(OG_clans_filt, function(i) {
  tips = substr(i$tip.label, 6, 100) 
  qs = na.omit(match(tips, ssal.d$qseqid))
  ss = na.omit(match(tips, ssal.d$sseqid))
  if(length(qs)==0 | length(ss)==0 ) return(NA)
  qs.gene = ssal.d$qseqid[qs]; ss.gene = ssal.d$sseqid[ss]
  length(c(qs.gene, ss.gene))
}
)

table(dups.inclans.ssal, useNA = 'always')
#names(which(dups.inclans.ssal>2))

# omyk...
dups.inclans.omyk = sapply(OG_clans_filt, function(i) {
  tips = substr(i$tip.label, 6, 100) 
  qs = na.omit(match(tips, omyk.d$qseqid))
  ss = na.omit(match(tips, omyk.d$sseqid))
  if(length(qs)==0 | length(ss)==0 ) return(NA)
  qs.gene = omyk.d$qseqid[qs]; ss.gene = omyk.d$sseqid[ss]
  length(c(qs.gene, ss.gene))
}
)

table(dups.inclans.omyk, useNA = 'always')
#names(which(dups.inclans.omyk>2))

# dups.inclans.ssal
# 2    3    4    6    8   10 <NA> 
#   9017   14   58    3    1    1 9758 

# dups.inclans.omyk
# 2     3     4     8  <NA> 
#   6026     4    64     1 12757

# 
table(dups.inclans.omyk==2 & dups.inclans.ssal==2, useNA = 'always')
OG_clans_dups = OG_clans_filt[which(dups.inclans.omyk==2 & dups.inclans.ssal==2)]
length(dups.inclans.omyk); length(OG_clans_filt)

dup.table = lapply(OG_clans_dups, function(i) {
  ssal.qs = na.omit(match(substr(i$tip.label,  6, 100), ssal.d$qseqid))
  ssal.ss = na.omit(match(substr(i$tip.label,  6, 100), ssal.d$sseqid))
  omyk.qs = na.omit(match(substr(i$tip.label,  6, 100), omyk.d$qseqid))
  omyk.ss = na.omit(match(substr(i$tip.label,  6, 100), omyk.d$sseqid))
  
  eluc = grep('Eluc', i$tip.label, value = T)
  
  if(ssal.qs==ssal.ss & omyk.qs==omyk.ss & length(eluc)==1) { return(c(eluc, ssal.d$qseqid[ssal.qs], ssal.d$sseqid[ssal.qs],
                                                                       omyk.d$qseqid[omyk.ss], omyk.d$sseqid[omyk.ss]) )}
  if(ssal.qs!=ssal.ss | omyk.qs!=omyk.ss | length(eluc)>1) return('PROBLEM')
})
names(dup.table) <- names(OG_clans_dups)
table(sapply(dup.table, function(i) length(grep('PROBLEM', i)))) # clans we can use

dup.table = dup.table[sapply(dup.table, function(i) length(grep('PROBLEM', i)))==0]
length(dup.table)

library(plyr)
dup.table = ldply(dup.table, function(i) {
  t(data.frame(i, stringsAsFactors = F))
})


head(dup.table) ## Quintuplets table to use in code below to identify "cherries" and orthologs between ssal and omyk



##--------------------------------------------------##
##       Identify Ss4R-cherries and orthologs       ##
##--------------------------------------------------##



## for loop to check consistency between dup table and clan topology
# setting some objects
test.tree = list(); test.cherry1 = c(); test.cherry2 = c(); 
clades=list(); mrca.combo1 = c() ; mrca.combo2 = c()
cherry1 = c(); cherry2 = c()
topology.results = c(); lore.results = c(); lore.taxa.names = c(); lore.taxa.results = c()
dup_cluster_df = matrix(rep(NA), nrow = nrow(dup.table), ncol=4)

for(i in 1:nrow(dup.table)){
  #print(i)
  test.tree <- OG_clans_expr[[dup.table$.id[i]]]
  test.tree <- auto.root(test.tree)$rooted.clans # Alternative 1
  #test.tree <- unroot(test.tree) # Alternative 2 ==> gives 2158 good and 2985 inconsisten topologies
  #print(test.tree$tip.label)
  #plot(test.tree); nodelabels()
  
  mrca.combo1 = getMRCA(test.tree, tip = c(paste('Ssal|', dup.table[i,3], sep=''), paste('Omyk|', dup.table[i,5], sep='')))
  mrca.combo2 = getMRCA(test.tree, tip = c(paste('Ssal|', dup.table[i,3], sep=''), paste('Omyk|', dup.table[i,6], sep='')))
  test.cherry1 = unlist(Descendants(test.tree, mrca.combo1, type='tips'))
  test.cherry2 = unlist(Descendants(test.tree, mrca.combo2, type='tips'))
  
  # choose which combo of ssal and omyk dup gives a mrca with smallest descedants clade size
  if(length(test.cherry1)<length(test.cherry2)) { clades = list(clade1=c(3,5), clade2=c(4,6))
  } else clades = list(clade1=c(3,6), clade2=c(4,5))
  
  dup_cluster_df[i, 1] <- paste('Ssal|', dup.table[i, clades[[1]][1]], sep='')
  dup_cluster_df[i, 2] <- paste('Omyk|', dup.table[i, clades[[1]][2]], sep='')
  dup_cluster_df[i, 3] <- paste('Ssal|', dup.table[i, clades[[2]][1]], sep='')
  dup_cluster_df[i, 4] <- paste('Omyk|', dup.table[i, clades[[2]][2]], sep='')
  
  
  # test if both cherry clades only contains 1 ssal and 1 omyk dup
  cherry1 = test.tree$tip.label[unlist(Descendants(test.tree, getMRCA(test.tree, c(paste('Ssal|', dup.table[i, clades[[1]][1]], sep=''),
                                                                                   paste('Omyk|', dup.table[i, clades[[1]][2]], sep=''))), type='tips'))]
  
  cherry2 = test.tree$tip.label[unlist(Descendants(test.tree, getMRCA(test.tree, c(paste('Ssal|', dup.table[i, clades[[2]][1]], sep=''),
                                                                                   paste('Omyk|', dup.table[i, clades[[2]][2]], sep=''))), type='tips'))]
  
  if(any(cherry1 %in% cherry2))  topology.results[i] <- 'Topology Incongruence'
  if(!any(cherry1 %in% cherry2)) topology.results[i] <- 'Good Topology'
  
  lore.results[i] <- sum(sapply(c('Tthy', 'Ssal', 'Omyk', 'Okis'), function(i){
    ifelse(length(grep(i,  test.tree$tip.label))>0,  is.monophyletic(test.tree, grep(i, test.tree$tip.label)), NA)
  } 
  ), na.rm = T)
  
  lore.taxa.names = sapply(c('Tthy', 'Ssal', 'Omyk', 'Okis'), function(i){
    ifelse(length(grep(i,  test.tree$tip.label))>0,  is.monophyletic(test.tree, grep(i, test.tree$tip.label)), NA)
  } 
  )
  
  lore.taxa.results[i] <- paste(sort(names(lore.taxa.names)[lore.taxa.names]), collapse='|')
  
}


#MelB: lore.taxa.names/results are not really relevant for you at this stage 

colnames(dup_cluster_df) <- c('Ssal1', 'Omyk1', 'Ssal2', 'Omyk2')
dup_cluster_fix = cbind(dup.table[,1:2], dup_cluster_df)
dup_cluster_phylofilt = dup_cluster_fix[topology.results!='Topology Incongruence', ]
nb.dupsinclusters.filt <- nrow(dup_cluster_phylofilt)

head(nb.dupsinclusters.filt) # final quintuplettable sorted based on tree topology information (ie ortholog/ohnolog relationships)


