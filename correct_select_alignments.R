#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Changing OG alignments to clan alignments
##Correcting names and subselection for clans
###################################################################



#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

install.packages("Rcurl", repos = "http://cran.rstudio.com/")
install.packages("seqinr", repos = "http://cran.rstudio.com/")

require("seqinr")
library("RCurl")


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


# get file path with all alignments
alignment.files <- dir(paste('C:/Users/meeldurb/Google Drive/',
                             'Master internship phylogenetics salmonids/',
                             'Salmonid_genomics_resources/Orthologs_homeologs/',
                             'orthogroups.03.06.2017/Alignments/', sep = ''), 
                       full.names = T)

# table with Omyk gene and protein names
Omyk.prot2gene = read.table(paste('C:/users/meeldurb/Google Drive/', 
                                  'Master internship phylogenetics salmonids/',
                                  'Salmonid_genomics_resources/Orthologs_homeologs/', 
                                  'orthogroups.03.06.2017/Omyk.gene2protein.table.txt', 
                                  sep = ''), header = T)



#--------------------------------------#
##_____ Correct alignment files  _____##
#--------------------------------------#

clan <- "OG0008581_1."
clan <- "OG0008390_1."
clan <- "OG0008707_1."
# count <- 0

# creating the folder to save the data in
outfolder <- "Alignments_aa_corrected/"  
if (!file.exists(outfolder)){ 
  dir.create(outfolder)
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