#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##preparing files for BEAST analysis
###################################################################



#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

install.packages("RCurl", repos = "http://cran.rstudio.com/")


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


# get file path with all corrected alignments
alignment.files <- dir(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                             'Master_internship_phylogenetics/',
                             'phylogenetics/Alignments_aa_corrected', sep = ''), 
                             full.names = T)


# xml file we need to import the alingment in
xml <- paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/', 
             'phylogenetics/dummy_xml/secondary_constr_aa.xml', sep = '')


# Table with duplicate pairs
dup_cluster_phylofilt <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
                               'Beast_dating_salmonids/RData/',
                               '20170801_duplicate_clans_filtered_aa.RData', 
                                sep = ''))



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
  clan.id <- sub('.*(OG\\d*_\\d*.)_corr.fa', '\\1', ali)
  cat(clan.id)
  if (clan.id %in% dup_cluster_phylofilt[,1]){
  convertcmd <- paste("python transform_XML.py", xml, ali,  dup_cluster_phylofilt)
  system(convertcmd)
  } else {
    print ("clan not in duplicate table")
  }
}








