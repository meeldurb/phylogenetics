#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Function for loading Rdata directly into a variable
##to avaid trivial name giving to the loaded RData files
###################################################################




loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}