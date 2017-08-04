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
library("RCurl")

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



#-----------------------#
##_____ Load data _____##
#-----------------------#

# get clans to obtain the clan structure of the trees
OG_clans_dupl <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                 '/Beast_dating_salmonids/RData/',
                                 'Clans_2analyze_inBeast_withduplicates_aa.RData', sep = ''))


# get file path with all corrected alignments
alignment.files <- dir(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                             'Master_internship_phylogenetics/',
                             'phylogenetics/Alignments_aa_corrected', sep = ''), 
                             full.names = T)


# xml file we need to import the alingment in
xml <- readLines(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                       'Master_internship_phylogenetics/', 
                       'phylogenetics/dummy_xml/secondary_constr_aa.xml', sep = ''))


# Table with duplicate pairs
dup_cluster_phylofilt <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                                         'Beast_dating_salmonids/RData/',
                                         '20170801_duplicate_clans_filtered_aa.RData', 
                                         sep = ''))


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




# orthologs.list <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
#                                  'Beast_dating_salmonids/orthologs_list_20161115.RData',
#                                  sep = ''))

# claninformation is contained in  OG_clans_dupl

#------------------------------------------------------#
##_____ Make BEAST files, example 1 OG alingment _____##
#------------------------------------------------------#


for (ali in alignment.files){
  convertcmd <- paste("python transform_XML.py", ali, xml, OG)
  system(convertcmd)
  
}




# loop to change the alignment files
for(clan in names(OG_clans_dupl)){
  cat(clan, '\n')
  # give the corrected alignment file a name and location it needs to be written to
  fileout <- paste(outfolder, clan, "_corr.fa", sep="")
  if (!file.exists(fileout)) {
  # get original alignment ID
  ali.id <- gsub(".*(OG\\d*)_\\d*\\.", '\\1', clan)
  # position where the alignment file is found
  ali.file.pos <- match(ali.id, gsub(".*(OG\\d*).fa", '\\1', alignment.files))
  ali.file <- alignment.files[ali.file.pos]
  # checking if file is not empty
    if (!file.size(ali.file) == 0){
      # only isolate sequences that are also in the tip.labels of the clans
      clan.select <- OG_clans_dupl[[clan]]
      clan.genes <- clan.select$tip.label
      seqs <- read.fasta(ali.file)
      # fixing names in seqs to resemble names in clans tip.labels
      # remove double species name
      names(seqs) <- lapply(names(seqs), function(i){
          gsub("\\w*_(\\w*\\|.*)", "\\1", i) })
      # fix Omyk names from proteinID to geneID
      # taken from Omyk.prot2gene
      names(seqs) <- lapply(names(seqs), function(i){
        select.row <- match(sub('Omyk\\|', '', i), Omyk.prot2gene$protein) 
        if (!is.na(select.row)) {
        i <- paste('Omyk|', Omyk.prot2gene$gene_id[select.row], sep = '')
        } else {
          i
        }
      })
      # first check length, then check if they contain the same elements
      if (length(names(seqs)) == length(clan.genes)){
        if (all.equal(sort(names(seqs)), sort(clan.genes))) {
        print ("tip.labels are the same")
        # change the alignment file for only the names
        write.fasta(seqs, names(seqs), fileout)
        }
        } else {
          print ("tip labels are not the same, making subselection of seqs")
          # select only the sequences that resemble the tip.labels
          selected.seqs <- seqs[clan.genes]
          write.fasta(selected.seqs, names(selected.seqs), fileout)
      }
      } else {
         print ("file is empty")
       }
  } else {
    print("corrected .fa file already was written")
  }
}

# change OG ID number
# import alignment
# import monophyletic groups or priors
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
  
  clan.info = OG_clans_dupl[clan]
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
