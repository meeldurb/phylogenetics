#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##preparing files for BEAST analysis
###################################################################





#####___________install & load packages and functions____________####

install.packages("ape", repos = "http://cran.rstudio.com/")
install.packages("phangorn", repos = "http://cran.rstudio.com/")
install.packages("plyr", repos = "http://cran.rstudio.com/")
install.packages("seqinr", repos = "http://cran.rstudio.com/")

library("ape")
library("phangorn")
library("plyr")
require("seqinr")

source(paste('C:/Users/meeldurb/Documents/WUR/Master',
             '/Internship/Beast_dating_salmonids/Beast_Rscripts/',
             'auto.root_salmonid_clans_v1.R', sep = "" ))

source(paste('C:/Users/meeldurb/Documents/WUR/Master',
             '/Internship/Beast_dating_salmonids/Beast_Rscripts/',
             'orthologfinder_oct2016.R', sep = ""))

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

####________ Example tree _______######
# select only the clans tht have at least 1 E.luc 
# and select only the clans that have 
# at least 2 of the Ssal, Omyk or Tthy species.
clans.selection = clans[sapply(clans, function(i) {
   length(grep('Eluc', i$tip.label))>0 & 
    length(grep('Ssal|Omyk|Tthy', i$tip.label))>=2
  })]
length(clans.selection)
plot(auto.root(clans.selection[[1000]])$rooted.clan)

# get orthologs of all the clans
test.all <- lapply(clans.selection, get.ortholog, verbouse=F)
table(sapply(test.all, function(i) i$class))
test.og <- get.ortholog(clans.selection[[1000]])
tree <- clans.selection[[1000]]
plot(tree)


#####________Make BEAST files________#####

# get alignment file path
alis <- file.path('C:/Users/meeldurb/Documents/WUR/Master',
                  '/Internship/Beast_dating_salmonids/',
                  'OG0000115_2_aln.fa')

# get preparation xml file
xml <- readLines(file.path('C:/Users/meeldurb/Documents/WUR/Master',
                           '/Internship/Beast_dating_salmonids/Dummy_beast_input/',
                           'test_beauti.xml'))


# get the orthologs
load(file.path('C:/Users/meeldurb/Documents/WUR/Master',
               '/Internship/Beast_dating_salmonids/',
               'orthologs_list_20161115.RData'))


count = 0

for(i in alis){
  
  #i = alis # make the xml for single alignment/OG
  
  count <- count+1 # check for loop
  temp <- xml
  # read all alingments files as fasta, but as one whole string
  d <- read.fasta(i, as.string=T)
  
  # get ortholog info ready:
  # IMPORTANT: alignments must have the name: clan_aln.fa for this to work
  # remove the _aln.fa name of the file
  clan = paste(sub('_aln.fa', '', basename(i)), '.', sep='') 

  # the orthologs.list has a "."
  #at the end of every OG (ortholog group) ID
  # we need to select this OG ID coming from
  # the alignment filename in this ortholog.list
  clan.info  = orthologs.list[[clan]]
  
  # we can plot the tree from this OG ID
  plot(clan.info$clan.tree)
  
  # get the monophyletic groups or orthogroups from the 
  # clans that were taken from the orthologs.list
  OP1 = clan.info$OP1 
  # split the names, these are separated by a "+"
  # you get a list after strsplit, you need it to be a vector
  OP1 = unlist(strsplit(names(OP1), '\\+')) 
  # OP1 should be equal to 2
  if(length(OP1) == 2) {
    # in which positions of d are these genes.
    which.OP1 <- sapply(sub('Omyk2\\||Ssal\\|', '', OP1), function(i) grep(i, names(d)))
  } else {
    which.OP1 <- NULL
  }
  # if OP2 is not nul, split by +
  if(!is.null(clan.info$OP2)) {
    OP2 = clan.info$OP2; OP2 = unlist(strsplit(names(OP2), '\\+'))
    # in which positions of d are these genes.
    which.OP2 = sapply(sub('Omyk2\\||Ssal\\|', '', OP2), function(i) grep(i, names(d)))
  }
  
  # linking the OP1 seqs to haplotype seqs
  
  #make data alignment strings line 7
  alignment.section <- c() 
  a <- c()
  b <- c()
  c <- c()
  for(n in 1:length(d)){
    a = paste('                    <sequence id="seq_', names(d)[n], '"', sep='')
    b = paste('taxon="', names(d)[n], '"', sep='')
    c = paste('totalcount="4" value="', d[[n]], '"/>', sep='')
    alignment.section[n] <- paste(a, b , c)
  }
  
  # lines to insert - wait until inserted OP1, OP2, and Eluc+ssal
  
  # make OP1 section: lines 74 + + <OP1.her>
  OP1.seqs = c()
  if(!is.null(OP1)){
    for(n in which.OP1){
      OP1.seqs = c(paste('                    <taxon id="', names(d[n]), '" spec="Taxon"/>', sep=''), OP1.seqs)
    }
  }
  
  # make OP2 section: lines 83 + + <OP2.her>
  OP2.seqs = c()
  if(!is.null(OP2)){
    for(n in which.OP2){
      OP2.seqs = c(paste('                    <taxon id="', names(d[n]), '" spec="Taxon"/>', sep=''), OP2.seqs) 
    }
  }
  
  # make Eluc+salmonids section: <esox.mono.her> line 92
  eluc = paste('                    <taxon id="', grep('Eluc', names(d), value = T), '" spec="Taxon"/>', sep='')
  others = grep('Eluc|Gacu|Olat', names(d), value = T, invert = T) # must define all non salmonids
  for(n in 1:length(others)){
    others[n] = paste('                    <taxon idref="', others[n], '"/>', sep='') # 
  }
  eluc_salmonids = c(eluc, others)
  
  # change idref to id for species not in OP1 or OP2
  #### WOW: Weird coding fuck up..
  id.previous = unlist(sapply(c(OP1, OP2), function(i) grep(i, eluc_salmonids, fixed=T)))
  id.esox = grep('Eluc', eluc_salmonids)
  swap.others = eluc_salmonids[!1:length(eluc_salmonids) %in% c(id.previous, id.esox)]
  swap.others = sub('"/>', ' spec=\"Taxon\"/>', sub('idref', 'id', others))
  eluc_salmonids[!1:length(eluc_salmonids) %in% c(id.previous, id.esox)] <- swap.others
  
  
  ### writing out
  
  # ==> insert lower lines first
  ## ==> eluc+salmonid
  xml.out = c(xml[1:6], alignment.section, xml[8:73], OP1.seqs, xml[75:82], OP2.seqs, xml[84:91], eluc_salmonids, xml[93:length(xml)])
  xml.out = sub('dummy_alignment.log', paste('BEAST_', sub('_aln.*', '', basename(i)), '.log', sep=''), xml.out)
  # NOTE: ==> remove monophyletic from OP1 and OP2 xml.out = sub(' monophyletic="true"', '', xml.out)
  # ==> FIX: output file name for trees!! 
  
  ## change out path
  out.path = file.path('C:/Users/meeldurb/Documents/WUR/Master',
                       '/Internship/Beast_dating_salmonids/Dummy_beast_input/')
  writeLines(xml.out, con = paste(out.path, paste('BEAST_', sub('_aln.*', '', basename(i)), '.xml', sep='')))
  
}
