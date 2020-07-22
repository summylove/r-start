rm(list=ls())

options(stringsAsFactors = FALSE) #禁止chr转成factor
Sys.setlocale("LC_ALL","English") #Sys.setenv(LANGUAGE = "en") 
#Sys.setlocale("LC_COLLATE", "C")
options(encoding = "utf-8")

pl=c('GEOquery','ggplot2') 
lapply(pl, function(i){
  if (!require(i,character.only = TRUE, quietly = TRUE)){
    BiocManager::install(i,character.only = TRUE)
  }  
  if (!require(i,character.only = TRUE, quietly = TRUE)){
    install.packages(i,character.only = TRUE)
  } 
  print(i,character.only = TRUE)
  require(i, character.only = TRUE)
})
(.packages())

