#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##preparing files for BEAST analysis
###################################################################



#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

# dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
# 
# install.packages('RCurl', Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )


library("RCurl")


# loadRData in variable function
eval(parse(text = getURL(paste("https://raw.githubusercontent.com/",
                               "meeldurb/phylogenetics/master/loadRData.R", 
                               sep = ""), ssl.verifypeer = FALSE)))



#-----------------------#
##_____ Load data _____##
#-----------------------#

# get clans to obtain the clan structure of the trees
OG_clans_dupl <- loadRData(paste('C:/Users/meeldurb/Dropbox/',
                                     'Melanie/Master_internship_',
                                     'phylogenetics/phylogenetics/RData/',
                                     '20170913-Clans_forBEAST_dupssandElucfilt_nt.RData', 
                                     sep = ''))


# get file path with all corrected alignments
alignment.files <- dir(paste('C:/Users/meeldurb/Google Drive/',
                            'Master internship phylogenetics salmonids/',
                            'Salmonid_genomics_resources/Orthologs_homeologs/',
                            'orthogroups_2017/cds_macse_nt_align', sep = ''), 
                             full.names = T)


# xml file we need to import the alingment in
xml <- paste('C:/Users/meeldurb/Dropbox/Melanie/',
             'Master_internship_phylogenetics/',
             'phylogenetics/dummy_xml/secondary_constr_nt.xml', sep = '')


# 
# # Table with duplicate pairs
# dup_cluster_phylofilt <- loadRData(paste('C:/Users/meeldurb/Dropbox/Melanie/',
#                                'Beast_dating_salmonids/RData/',
#                                '20170801_duplicate_clans_filtered_aa.RData', 
#                                 sep = ''))

# Table with duplicate pairs
dup_table <- paste('C:/Users/meeldurb/Dropbox/Melanie/',
                   'Master_internship_phylogenetics/', 
                   'phylogenetics/RData/',
                   '20170913_duplicate_clans_filtered_final_nt.csv', 
                    sep = '')

#------------------------------#
##_____ Make BEAST files _____##
#------------------------------#



#clan.id <- 'OG0000301_1.'
for (ali in alignment.files){
  clan.id <- sub('.*(OG\\d*_\\d*.+)aln', '\\1', ali)
  cat(clan.id, '\n')
  if (clan.id %in% dup.table[,1]){
  convertcmd <- paste("python 20170913-transform_XML_new.py ", xml, ' "', ali, '" ', 
                      dup_table, sep = "")
  system(convertcmd)
  } else {
    print ("clan not in duplicate table")
  }
}








