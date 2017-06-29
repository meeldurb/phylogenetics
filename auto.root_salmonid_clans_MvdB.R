#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for setting the root in salmonid phylogenetic trees
###################################################################


#####___________ install & load packages and functions ____________####

install.packages("ape", repos = "http://cran.rstudio.com/")
install.packages("phangorn", repos = "http://cran.rstudio.com/")

library("ape")
library("phangorn")


######_________________ load example data _________________#####
# load tree
tr.n = "((((((((Omyk2|CIGENEomyV6.41072:0.00055,Omyk|GSONMT00002897001:0.03238)0.966:0.00292,Ssal|XP_013993288.1:0.00055):0.00159,Okis|gb|GDQG01036265.1||m.63677:0.00934)0.823:0.00152,Tthy|Tthy_00032335-RD:0.00413)0.505:0.00907,(((Omyk2|CIGENEomyV6.34600:0.00054,Omyk|GSONMT00074661001:0.00055)0.984:0.0156,Ssal|XP_014052884.1:0.0029)0.966:0.00343,Tthy|Tthy_00027454-RA:0.00951)0.969:0.00438)0.999:0.00155,Eluc|XP_010902685.1:0.01043)0.816:0.0661,Olat|ENSORLP00000000831.1:0.02215)1.000:0,Gacu|ENSGACP00000026666.1:0.01773)Root;"
tree = read.tree(text=tr.n)
plot(tree)

# load tree with mammals
tr.m = "((((((((Omyk2|CIGENEomyV6.41072:0.00055,Omyk|GSONMT00002897001:0.03238)0.966:0.00292,Ssal|XP_013993288.1:0.00055):0.00159,Okis|gb|GDQG01036265.1||m.63677:0.00934)0.823:0.00152,Tthy|Tthy_00032335-RD:0.00413)0.505:0.00907,(((Omyk2|CIGENEomyV6.34600:0.00054,Omyk|GSONMT00074661001:0.00055)0.984:0.0156,Ssal|XP_014052884.1:0.0029)0.966:0.00343,Tthy|Tthy_00027454-RA:0.00951)0.969:0.00438)0.999:0.00155,Eluc|XP_010902685.1:0.01043)0.816:0.0661,Mmus|ENSORLP00000000831.1:0.02215)1.000:0,Hsap|ENSGACP00000026666.1:0.01773)Root;"
tree = read.tree(text=tr.m)
plot(tree)


# tree = rtree(10)
# plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()
Ancestors(tree, 1:3, "all")
Children(tree, 23)
Descendants(tree, 13, "tips")
Siblings(tree, 21)
# 


input.tree = tree


######_________________ function to put the root on the trees _________________#####
auto.root <- function(input.tree, outgroup=c('Mmus', 'Hsap', 'Drer', 'Gacu', 'Olat')){
  
  if(class(input.tree)=='list'){ #the input tree that I am giving now is of class "Phylo"
    tree = input.tree[[1]]
  } else{
  tree = unroot(input.tree)
  }
  
  find.outgroup <- sapply(outgroup, function(i) grep(i, tree$tip.label))
  get.root <- unlist(find.outgroup)
  
  
  #put.root = unlist(sapply(outgroup, function(i) grep(i, tree$tip.label)))
  
  # if NO root sequences - defined by outgroup parameter - perform a midpoint rooting
  if(length(get.root) == 0){
    return(list(rooted.clans <- midpoint(tree), root.info='midpoint'), fail.root=F)
  }
  
  # fix names
  names(get.root) <- sub('[:0-10000:]', '', names(get.root))


  # check whether both mammal groups are present
  if(sum(unique(names(get.root)) %in% c('Mmus', 'Hsap')) == 2){
    cat('Minimum two mammal outgroups\n')
    mam.root = grep('Mmus|Hsap', tree$tip.label)
    
    # check whether there is a monophyletic mammal outgroup
    if(is.monophyletic(tree, mam.root, reroot = F)) { 
      cat('Monophyletic mammal outgroup\n')
      root.node = getMRCA(tree, mam.root)
      tr = try(root(tree, node = root.node, resolve.root = T), silent=F)
      if(class(tr)!='try-error') { 
        return(list(rooted.clans=tr, root.info='Monophyletic_mammal_root', fail.root=F)) 
        }
      if(class(tr)=='try-error') { 
        cat('Problems with mammal node rooting\n....rooting with human\n')
        tr = try(root(tree, outgroup=grep('Hsap', tree$tip.label)[1], resolve.root = T), silent=T)
          if(class(tr)!='try-error') return(list(rooted.clans=tr, root.info='Human_root - mammals not monophyletic', fail.root=F))
          if(class(tr)=='try-error') return(list(rooted.clans=input.tree, root.info='Human_root - mammals not monophyletic', fail.root=T))
        }
      }
    
    # or non monophyletic mammal outgroup
    if(!is.monophyletic(tree, mam.root, reroot = F)) { 
      cat('NOT Monophyletic mammal outgroup\n')
      root.node = getMRCA(tree, mam.root)
      tr = try(root(tree, node = root.node, resolve.root = T), silent = F)
      if(class(tr)!='try-error') { 
        cat('Rooting with MRCA-node of all outgroups\n')
        return(list(rooted.clans=tr, root.info='MRCA-node of all outgroups', fail.root=F))
        }
      if(class(tr)=='try-error') { 
        cat('Problems with mammal node rooting\n....rooting with a random human\n')
        tr = try(root(tree, outgroup=grep('Hsap', tree$tip.label)[1], resolve.root = T), silent=T)
        if(class(tr)!='try-error') {
          return(list(rooted.clans=tr, root.info='Random_human', fail.root=F))
        }
        if(class(tr)=='try-error') {
          return(list(rooted.clans=input.tree, root.info='Random_human', fail.root=T))
        }
      }
    }
    
  }
  # check if only one mammal is present
  if(sum(unique(names(get.root)) %in% c('Mmus', 'Hsap'))==1){
    cat('One mammal outgroup\n')
    tr = try(root(tree, outgroup= put.root[1], resolve.root = T), silent=T)
    if(class(tr)!='try-error'){
      return(list(rooted.clans=tr, root.info='Single_mammal_outgroup', fail.root=F))
    }
    if(class(tr)=='try-error') {
      return(list(rooted.clans=input.tree, root.info='Single_mammal_outgroup', fail.root=T))
    }
  }
  # checking if no mammal is found
  if(sum(unique(names(get.root)) %in% c('Mmus', 'Hsap'))==0){
    cat('No mammal outgroup\n')
    tr = try(ingroupMRCA.rooting(tree), silent = T)
    if(class(tr)!= 'try-error') {
      return(list(rooted.clans=tr, root.info='ingroup_MRCA', fail.root=F))
      }
    if(class(tr)=='try-error') {
      tr = try(root(tree, outgroup = get.root[1], resolve.root = T), silent = T)
      if(class(tr)!= 'try-error') {
        return(list(rooted.clans=tr, root.info='Most_distant_teleost', fail.root=F))
      }
      if(class(tr)=='try-error') {
        return(list(rooted.clans=input.tree, root.info='Most_distant_teleost', fail.root=T))
      }
    }
  }
}

                     

# 
# fail.root = function(tree, root.with=c('node', 'outgroup'), root=NULL){
#   if(root.with=='outgroup')
#     return(tryCatch(root(tree, root, resolve.root=TRUE), 
#                     error=function(e) try(root(unroot(tree), root, resolve.root=TRUE), 
#                                           silent=TRUE)))
#   if(root.with=='node')
#     return(tryCatch(root(tree, node=root, resolve.root=TRUE),
#                     error=function(e) try(root(unroot(tree), node=root, resolve.root=TRUE), 
#                                           silent=TRUE)))
# }

######___________ function to get the root when no mammals _____________#####
ingroupMRCA.rooting =function(tree, ingroup = c('Eluc', 'Tthy', 'Ssal', 
                                                'Omyk', 'Okis')){
  # get only the organism name
  #spec = gsub("\\|.*", "", tree$tip.label)
  spec = substr(tree$tip.label, 1, 4)
  in.tips = which(spec %in% ingroup == T)
  root.node = getMRCA(tree, in.tips)
  root(tree, node = root.node, resolve.root = T)
}


