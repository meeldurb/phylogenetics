#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch & Simen Rod Sandve
##Script for installing all packages used in project
###################################################################


#------------------------------------------#
##_____ load libraries and functions _____##
#------------------------------------------#

dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

install.packages('RCurl', Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("ape", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("phangorn", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("seqinr", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("plyr", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages('ips',  Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" ) 
install.packages("tree", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("devtools", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("RSQLite", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages('phytools', Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
