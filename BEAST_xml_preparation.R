#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##preparing files for BEAST analysis
###################################################################



#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

install.packages("ape", repos = "http://cran.rstudio.com/")
install.packages("phangorn", repos = "http://cran.rstudio.com/")
install.packages("plyr", repos = "http://cran.rstudio.com/")
install.packages("seqinr", repos = "http://cran.rstudio.com/")

library("ape")
library("phangorn")
library("plyr")
require("seqinr")

source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/phylogenomics/',
             'auto.root_salmonid_clans.R', sep = "" ))

source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/phylogenetics/',
             'auto.root_salmonid_clans_MvdB.R', sep = "" ))

source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/phylogenetics/',
              'orthologfinder_MvdB.R', sep = ""))

source(paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'Phylogenetics/loadRData.R', sep = ''))



#-----------------------#
##_____ Load data _____##
#-----------------------#

# get clans to obtain the clan structure of the trees
OG_clans_dupl <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                 '/Beast_dating_salmonids/RData/',
                                 'Clans_2analyze_inBeast_withduplicates_aa.RData', sep = ''))


# prepare duplicated pairs table

# get duplicate pairs of Oncorhynchus mykiss (rainbow trout)
Omyk.duplicates <- loadRData(paste('C:/Users/meeldurb/Google Drive/',
                            'Master internship phylogenetics salmonids/',
                            'Salmonid_genomics_resources/Orthologs_homeologs/',
                            'Homeologs/OmykV6_2016_best_in_homelogRegions_',
                            'minpident80_mincov50_phylofiltered.RData', sep = ''))[,1:2]
  

# get duplicate pairs of Salmo salar (Atlantic salmon)
Ssal.duplicates <- loadRData(paste('C:/Users/meeldurb/Google Drive/',
                            'Master internship phylogenetics salmonids/',
                            'Salmonid_genomics_resources/Orthologs_homeologs/',
                            'Homeologs/RefSeq_GarethLongest_2016_best_',
                            'in_homelogRegions_minpident80_mincov50_',
                            'phylofiltered.RData', sep = ''))[,1:2]

# add organism naming in front on duplicate pairs
omyk.dup <- apply(Omyk.duplicates, 2, function(i) paste('Omyk2|', i, sep=''))
ssal.dup <- apply(Ssal.duplicates, 2, function(i) paste('Ssal|', i, sep=''))

# bind all homologs in complete dataframe
dup.table <- data.frame(rbind(omyk.dup, ssal.dup), 
                        stringsAsFactors = F)




#--------------------------#
##_____ Example tree _____##
#--------------------------#

# select only the OG groups with clans that have at least 1 E.luc 
# and at least 2 of the Ssal, Omyk or Tthy species.
# else no proper tree can be drawn
clans.selection <- OG_clans_dupl[sapply(OG_clans_dupl, function(i) {
  length(grep('Eluc', i$tip.label))>0 & 
    length(grep('Ssal|Omyk|Tthy', i$tip.label))>=2
  })]

length(clans.selection)
par(mfrow = c(1,1))
plot(auto.root(clans.selection[[150]])$rooted.clan)

# get orthologs of all the clans
test.all <- lapply(clans.selection, get.ortholog, verbouse=F)
table(sapply(test.all, function(i) i$class))
test.og <- get.ortholog(clans.selection[[150]])

tree <- clans.selection[[150]]
plot(tree)
nodelabels()
tiplabels()


#------------------------------#
##_____ Make BEAST files _____##
#------------------------------#


# get file path with all alignments
alignment.files <- dir(paste('C:/Users/meeldurb/Google Drive/',
                       'Master internship phylogenetics salmonids/',
                       'Salmonid_genomics_resources/Orthologs_homeologs/',
                       'orthogroups.03.06.2017/Alignments/', sep = ''), 
                       full.names = T)


# get xml file we need to import the alingment in
xml <- readLines(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                       'Master_internship_phylogenetics/', 
                       'phylogenetics/dummy_xml/secondary_constr_aa.xml', sep = ''))


# orthologs.list <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
#                                  'Beast_dating_salmonids/orthologs_list_20161115.RData',
#                                  sep = ''))

# claninformation is contained in  OG_clans_dupl

#------------------------------------------------------#
##_____ Make BEAST files, example 1 OG alingment _____##
#------------------------------------------------------#



ali = paste("C:/Users/meeldurb/Google Drive/Master internship phylogenetics ",
            "salmonids/Salmonid_genomics_resources/Orthologs_homeologs/",
            "orthogroups.03.06.2017/Alignments/OG0008397.fa", sep = "")
count <- 0

for(ali in alignment.files){
  cat(ali, '\n')
  cat(count <- count + 1, '\n')
  # checking if file is .fa file and if not empty
  if (grepl('.fa', ali)){ 
    if (!file.size(ali) == 0){
  # read all alignments of the .fa file
  ali.fasta <- read.fasta(ali)
    
  
  } else {
      print ("file is empty")
    }
  }
} 

# import ID
pat.id = ".*(OG\\d*).fa"
OG.id <- sub(pattern = pat.id, "\\1", ali)

pat.id.xml = "OG\\d*"
id <- sub(pattern = pat.id.xml, replace = OG.id, x = xml[5])

xml[5]
# change datatype depending on nt/aa

# import alignment

# import the taxids

# import monophyletic groups or priors

# change the filenames

## do it in python

convertcmd <- paste("python transform_XML.py", ali, xml, OG)
system(convertcmd)

  # get ortholog info ready:
  # IMPORTANT: alignments must have the name: clan_aln.fa for this to work
  # remove the _aln.fa name of the file
  clan = paste(sub('.fa', '', basename(ali)), '.', sep='') 

  # the orthologs.list has a "."
  #at the end of every OG (ortholog group) ID
  # we need to select this OG ID coming from
  # the alignment filename in this ortholog.list
  # I AM NOT SURE WHY WE WOULD USE THIS ORTHOLOGS LIST
  #clan.info  = orthologs.list[[clan]]
  
  clan.info = OG_clans_dupl[[clan]]
  # we can plot the tree from this OG ID
  plot(clan)
  
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
