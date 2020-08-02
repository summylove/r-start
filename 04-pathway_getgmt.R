#https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
Homo_sapiens.gene_info=read.delim('fig5enrich/Homo_sapiens.gene_info')
Homo_sapiens.gene_info[1,]
colnames(Homo_sapiens.gene_info)
proteins=Homo_sapiens.gene_info[Homo_sapiens.gene_info$type_of_gene == "protein-coding",-c(1,4,10,16)]
table(Homo_sapiens.gene_info$type_of_gene)

cat(file='fig5enrich/Homo_sapiens.gene_info.txt',
    #paste0(Homo_sapiens.gene_info$GeneID[Homo_sapiens.gene_info$type_of_gene == "protein-coding"],collapse = ','),
    proteins$GeneID,
    sep="\n",append = T)

##block
block_length=2500
total_length=length(proteins$GeneID)
for (i in 1:(as.integer(total_length/block_length)+1)){
if(i<= as.integer(total_length/block_length)){
assign(paste0('block_',i),seq(1, total_length,  by=block_length)[i]:(seq(1, total_length, by=block_length)[i+1]-1))#i=1:7
}else{
  assign(paste0('block_',i),seq(1, total_length,  by=block_length)[i]:total_length) #i=8
}  
}

# # one file
# for(i in 1:(as.integer(total_length/block_length)+1) ){
#   cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',0,'.txt'),
#       paste0(proteins$GeneID[get(paste0('block_',i))],collapse = ' '),
#       #proteins$GeneID[get(paste0('block_',i))],
#       sep="\n",append = T)  
# }
# Homo_sapiens.gene_blocks=t(read.delim('fig5enrich/Homo_sapiens.gene_info_0.txt',sep = ' ',header = F))
# write.csv(Homo_sapiens.gene_blocks,file = 'fig5enrich/Homo_sapiens.gene_blocks.csv',na = '')

## single file
for(i in 1:(as.integer(total_length/block_length)+1) ){
  cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',i,'.txt'), 
      'gene',
  sep="\n",append = F) 
  cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',i,'.txt'),
      #paste0(proteins$GeneID[get(paste0('block_',i))],collapse = ' '),
      proteins$GeneID[get(paste0('block_',i))],
      sep="\n",append = T) 
  cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',i,'.csv'), 
      'gene',
      sep="\n",append = F) 
  cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',i,'.csv'),
      #paste0(proteins$GeneID[get(paste0('block_',i))],collapse = ''),
      proteins$GeneID[get(paste0('block_',i))],
      sep="\n",append = T) 
}

#http://metascape.org/gp/index.html#/main/step1

for (i in 1:(as.integer(total_length/block_length)+1)){
assign(paste0('GO_AllLists_block_',i), 
       read_csv(paste0("fig5enrich/",list.files('fig5enrich','^block_*')[seq(1,16,2)][i],"/Enrichment_GO/GO_AllLists.csv"))
       )
}
# 1+2
# 1+2 +3
GO_AllLists_=function(GO_AllLists_1,GO_AllLists_2){

GO_AllLists=merge(GO_AllLists_1,GO_AllLists_2,by='GO',all=T)
rownames(GO_AllLists)=GO_AllLists$GO
#GO_AllLists=GO_AllLists[,order(colnames(GO_AllLists))]
#GO_AllLists
#GO_AllLists_1=as.data.frame(t(data.frame(1:length(GO_AllLists_block_1[,c(1,3,4,10:12,15:16)]))))[-1,]
# colnames(GO_AllLists_block_1[,c(1,3,4,10:12,15:16)])
# colnames(GO_AllLists_1)
GO_AllLists_1=data.frame(GO=GO_AllLists$GO)
GO_AllLists_1$Category=apply(cbind(GO_AllLists$Category.x,GO_AllLists$Category.y),1,function(x){
  return(unique(na.omit(x)))
})
GO_AllLists_1$Description=apply(cbind(GO_AllLists$Description.x,GO_AllLists$Description.y),1,function(x){
  return(unique(na.omit(x)))
})
GO_AllLists_1$`#GeneInGO`=apply(cbind(GO_AllLists$`#GeneInGO.x`,GO_AllLists$`#GeneInGO.y`),1,function(x){
  return(unique(na.omit(x)))
})
GO_AllLists_1$`#GeneInGOAndHitList`=apply(cbind(GO_AllLists$`#GeneInGOAndHitList.x`,GO_AllLists$`#GeneInGOAndHitList.y`),1,function(x){
  return(sum(na.omit(x)))
})
GO_AllLists_1$`#GeneInHitList`=apply(cbind(GO_AllLists$`#GeneInHitList.x`,GO_AllLists$`#GeneInHitList.y`),1,function(x){
  return(sum(na.omit(x)))
})

GO_AllLists_1$GeneID=apply(cbind(GO_AllLists$GeneID.x,GO_AllLists$GeneID.y),1,function(x){
  return(paste(na.omit(x),collapse = '|'))
})
GO_AllLists_1$Hits=apply(cbind(GO_AllLists$Hits.x,GO_AllLists$Hits.y),1,function(x){
  return(paste(na.omit(x),collapse = '|'))
})
#GO_AllLists_1$GeneInGOAndHitList_1=GO_AllLists$`#GeneInGOAndHitList.x`; GO_AllLists_1$GeneInGOAndHitList_2=GO_AllLists$`#GeneInGOAndHitList.y`
return(GO_AllLists_1)
}
GO_AllLists_1_2=GO_AllLists_(GO_AllLists_block_1[,c(1,3,4,10:12,15:16)] ,GO_AllLists_block_2[,c(1,3,4,10:12,15:16)])
GO_AllLists_1_2_3=GO_AllLists_(GO_AllLists_1_2,GO_AllLists_block_3[,c(1,3,4,10:12,15:16)])
GO_AllLists_1_2_3_4=GO_AllLists_(GO_AllLists_1_2_3,GO_AllLists_block_4[,c(1,3,4,10:12,15:16)])
GO_AllLists_1_2_3_4_5=GO_AllLists_(GO_AllLists_1_2_3_4,GO_AllLists_block_5[,c(1,3,4,10:12,15:16)])
GO_AllLists_1_2_3_4_5_6=GO_AllLists_(GO_AllLists_1_2_3_4_5,GO_AllLists_block_6[,c(1,3,4,10:12,15:16)])
GO_AllLists_1_2_3_4_5_6_7=GO_AllLists_(GO_AllLists_1_2_3_4_5_6,GO_AllLists_block_7[,c(1,3,4,10:12,15:16)])
GO_AllLists_1_2_3_4_5_6_7_8=GO_AllLists_(GO_AllLists_1_2_3_4_5_6_7,GO_AllLists_block_8[,c(1,3,4,10:12,15:16)])

GO_AllLists=GO_AllLists_1_2_3_4_5_6_7_8
test=GO_AllLists[match(FINAL_GO$GO,GO_AllLists$GO),]
test$`#GeneInGOAndHitList` == test$`#GeneInGO`

##########################
table(Homo_sapiens.gene_info$type_of_gene)
nonproteins=Homo_sapiens.gene_info[Homo_sapiens.gene_info$type_of_gene == "other",-c(1,4,10,16)]
cat(file='fig5enrich/nonptrotein.txt',
    paste0(nonproteins$GeneID,collapse = ','),
    sep="\n",append = T)

assign(paste0('GO_AllLists_block_',9), 
       read_csv(paste0("fig5enrich/",list.files('fig5enrich','^block_9_*')[1],"/Enrichment_GO/GO_AllLists.csv"))
)
FINAL_GO$GO[match(FINAL_GO$GO,GO_AllLists_block_9$GO)]
GO_AllLists_1_2_3_4_5_6_7_8_9=GO_AllLists_(GO_AllLists_1_2_3_4_5_6_7_8,GO_AllLists_block_9[,c(1,3,4,10:12,15:16)])
GO_AllLists=GO_AllLists_1_2_3_4_5_6_7_8_9
test=GO_AllLists[match(FINAL_GO$GO,GO_AllLists$GO),]
test$`#GeneInGOAndHitList` == test$`#GeneInGO`

############################
table(Homo_sapiens.gene_info$type_of_gene)
nonproteins=Homo_sapiens.gene_info[Homo_sapiens.gene_info$type_of_gene == "pseudo",-c(1,4,10,16)]

##nonblock
nonblock_length=2500
total_length=length(nonproteins$GeneID)
for (i in 1:(as.integer(total_length/nonblock_length)+1)){
  if(i<= as.integer(total_length/nonblock_length)){
    assign(paste0('nonblock_',i),seq(1, total_length,  by=nonblock_length)[i]:(seq(1, total_length, by=nonblock_length)[i+1]-1))#i=1:7
  }else{
    assign(paste0('nonblock_',i),seq(1, total_length,  by=nonblock_length)[i]:total_length) #i=8
  }  
}

# # one file
# for(i in 1:(as.integer(total_length/nonblock_length)+1) ){
#   cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',0,'.txt'),
#       paste0(nonproteins$GeneID[get(paste0('nonblock_',i))],collapse = ' '),
#       #nonproteins$GeneID[get(paste0('nonblock_',i))],
#       sep="\n",append = T)  
# }
# Homo_sapiens.gene_nonblocks=t(read.delim('fig5enrich/Homo_sapiens.gene_info_0.txt',sep = ' ',header = F))
# write.csv(Homo_sapiens.gene_nonblocks,file = 'fig5enrich/Homo_sapiens.gene_nonblocks.csv',na = '')

## single file
for(i in 1:(as.integer(total_length/nonblock_length)+1) ){
  cat(file=paste0('fig5enrich/Homo_sapiens.nongene_info_',i,'.txt'), 
      'gene',
      sep="\n",append = F) 
  cat(file=paste0('fig5enrich/Homo_sapiens.nongene_info_',i,'.txt'),
      #paste0(nonproteins$GeneID[get(paste0('nonblock_',i))],collapse = ' '),
      nonproteins$GeneID[get(paste0('nonblock_',i))],
      sep="\n",append = T) 
  cat(file=paste0('fig5enrich/Homo_sapiens.nongene_info_',i,'.csv'), 
      'gene',
      sep="\n",append = F) 
  cat(file=paste0('fig5enrich/Homo_sapiens.nongene_info_',i,'.csv'),
      #paste0(nonproteins$GeneID[get(paste0('nonblock_',i))],collapse = ''),
      nonproteins$GeneID[get(paste0('nonblock_',i))],
      sep="\n",append = T) 
}

#http://metascape.org/gp/index.html#/main/step1

for (i in 1:4){
  assign(paste0('GO_AllLists_nonblock_',i), 
         read_csv(paste0("fig5enrich/",list.files('fig5enrich','^nonblock_*')[seq(1,8,2)][i],"/Enrichment_GO/GO_AllLists.csv"))
  )
}

GO_AllLists_non1_2=GO_AllLists_(GO_AllLists_nonblock_1[,c(1,3,4,10:12,15:16)] ,GO_AllLists_nonblock_2[,c(1,3,4,10:12,15:16)])
GO_AllLists_non1_2_3=GO_AllLists_(GO_AllLists_non1_2,GO_AllLists_nonblock_3[,c(1,3,4,10:12,15:16)])
GO_AllLists_non1_2_3_4=GO_AllLists_(GO_AllLists_non1_2_3,GO_AllLists_nonblock_4[,c(1,3,4,10:12,15:16)])

GO_AllLists_1_2_3_4_5_6_7_8_9_non=GO_AllLists_(GO_AllLists_1_2_3_4_5_6_7_8_9,GO_AllLists_non1_2_3_4)
GO_AllLists=GO_AllLists_1_2_3_4_5_6_7_8_9_non
test=GO_AllLists[match(FINAL_GO$GO,GO_AllLists$GO),]
test$`#GeneInGOAndHitList` - test$`#GeneInGO`


#######
 table(Homo_sapiens.gene_info$type_of_gene)
otherproteins=Homo_sapiens.gene_info[Homo_sapiens.gene_info$type_of_gene == "biological-region",-c(1,4,10,16)]
#########
# 
# table(Homo_sapiens.gene_info$type_of_gene)
# brproteins=Homo_sapiens.gene_info[Homo_sapiens.gene_info$type_of_gene == "biological-region",-c(1,4,10,16)]
# 
# ##brblock
# brblock_length=2500
# total_length=length(brproteins$GeneID)
# for (i in 1:(as.integer(total_length/brblock_length)+1)){
#   if(i<= as.integer(total_length/brblock_length)){
#     assign(paste0('brblock_',i),seq(1, total_length,  by=brblock_length)[i]:(seq(1, total_length, by=brblock_length)[i+1]-1))#i=1:7
#   }else{
#     assign(paste0('brblock_',i),seq(1, total_length,  by=brblock_length)[i]:total_length) #i=8
#   }  
# }
# #rm(list=ls(pattern = 'brblock_*'))
# 
# # # one file
# # for(i in 1:(as.integer(total_length/brblock_length)+1) ){
# #   cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',0,'.txt'),
# #       paste0(brproteins$GeneID[get(paste0('brblock_',i))],collapse = ' '),
# #       #brproteins$GeneID[get(paste0('brblock_',i))],
# #       sep="\n",append = T)  
# # }
# # Homo_sapiens.gene_brblocks=t(read.delim('fig5enrich/Homo_sapiens.gene_info_0.txt',sep = ' ',header = F))
# # write.csv(Homo_sapiens.gene_brblocks,file = 'fig5enrich/Homo_sapiens.gene_brblocks.csv',na = '')
# 
# ## single file
# for(i in 1:(as.integer(total_length/brblock_length)+1) ){
#   cat(file=paste0('fig5enrich/Homo_sapiens.brgene_info_',i,'.txt'), 
#       'gene',
#       sep="\n",append = F) 
#   cat(file=paste0('fig5enrich/Homo_sapiens.brgene_info_',i,'.txt'),
#       #paste0(brproteins$GeneID[get(paste0('brblock_',i))],collapse = ' '),
#       brproteins$GeneID[get(paste0('brblock_',i))],
#       sep="\n",append = T) 
#   cat(file=paste0('fig5enrich/Homo_sapiens.brgene_info_',i,'.csv'), 
#       'gene',
#       sep="\n",append = F) 
#   cat(file=paste0('fig5enrich/Homo_sapiens.brgene_info_',i,'.csv'),
#       #paste0(brproteins$GeneID[get(paste0('brblock_',i))],collapse = ''),
#       brproteins$GeneID[get(paste0('brblock_',i))],
#       sep="\n",append = T) 
# }

#########
# 
table(Homo_sapiens.gene_info$type_of_gene)
ncproteins=Homo_sapiens.gene_info[Homo_sapiens.gene_info$type_of_gene == "ncRNA",-c(1,4,10,16)]

##ncblock
ncblock_length=2500
total_length=length(ncproteins$GeneID)
for (i in 1:(as.integer(total_length/ncblock_length)+1)){
  if(i<= as.integer(total_length/ncblock_length)){
    assign(paste0('ncblock_',i),seq(1, total_length,  by=ncblock_length)[i]:(seq(1, total_length, by=ncblock_length)[i+1]-1))#i=1:7
  }else{
    assign(paste0('ncblock_',i),seq(1, total_length,  by=ncblock_length)[i]:total_length) #i=8
  }
}
#rm(list=ls(pattern = 'ncblock_*'))

# # one file
# for(i in 1:(as.integer(total_length/ncblock_length)+1) ){
#   cat(file=paste0('fig5enrich/Homo_sapiens.gene_info_',0,'.txt'),
#       paste0(ncproteins$GeneID[get(paste0('ncblock_',i))],collapse = ' '),
#       #ncproteins$GeneID[get(paste0('ncblock_',i))],
#       sep="\n",append = T)
# }
# Homo_sapiens.gene_ncblocks=t(read.delim('fig5enrich/Homo_sapiens.gene_info_0.txt',sep = ' ',header = F))
# write.csv(Homo_sapiens.gene_ncblocks,file = 'fig5enrich/Homo_sapiens.gene_ncblocks.csv',na = '')

## single file
for(i in 1:(as.integer(total_length/ncblock_length)+1) ){
  cat(file=paste0('fig5enrich/Homo_sapiens.ncgene_info_',i,'.txt'),
      'gene',
      sep="\n",append = F)
  cat(file=paste0('fig5enrich/Homo_sapiens.ncgene_info_',i,'.txt'),
      #paste0(ncproteins$GeneID[get(paste0('ncblock_',i))],collapse = ' '),
      ncproteins$GeneID[get(paste0('ncblock_',i))],
      sep="\n",append = T)
  cat(file=paste0('fig5enrich/Homo_sapiens.ncgene_info_',i,'.csv'),
      'gene',
      sep="\n",append = F)
  cat(file=paste0('fig5enrich/Homo_sapiens.ncgene_info_',i,'.csv'),
      #paste0(ncproteins$GeneID[get(paste0('ncblock_',i))],collapse = ''),
      ncproteins$GeneID[get(paste0('ncblock_',i))],
      sep="\n",append = T)
}


for (i in 1:2){
  assign(paste0('GO_AllLists_ncblock_',i), 
         read_csv(paste0("fig5enrich/",list.files('fig5enrich','^ncblock_*')[seq(1,4,2)][i],"/Enrichment_GO/GO_AllLists.csv"))
  )
}

GO_AllLists_nc1_2=GO_AllLists_(GO_AllLists_ncblock_1[,c(1,3,4,10:12,15:16)] ,GO_AllLists_ncblock_2[,c(1,3,4,10:12,15:16)])


GO_AllLists_1_2_3_4_5_6_7_8_9_nc=GO_AllLists_(GO_AllLists_1_2_3_4_5_6_7_8_9_non,GO_AllLists_nc1_2)
GO_AllLists=GO_AllLists_1_2_3_4_5_6_7_8_9_nc
test=GO_AllLists[match(FINAL_GO$GO,GO_AllLists$GO),]
test$`#GeneInGOAndHitList` - test$`#GeneInGO`

########################################
table(Homo_sapiens.gene_info$type_of_gene)
nopproteins=Homo_sapiens.gene_info[Homo_sapiens.gene_info$type_of_gene %in% c("rRNA",'scRNA','snoRNA','snRNA','tRNA'),-c(1,4,10,16)]
cat(file='fig5enrich/nopptrotein.txt',
    paste0(nopproteins$GeneID,collapse = ','),
    sep="\n",append = T)

assign(paste0('GO_AllLists_nopblock'), 
       read_csv(paste0("fig5enrich/nopblock_all.tt57slkyk/Enrichment_GO/GO_AllLists.csv"))
)

FINAL_GO$GO[match(FINAL_GO$GO,GO_AllLists_nopblock$GO)]
GO_AllLists_1_2_3_4_5_6_7_8_9_nop=GO_AllLists_(GO_AllLists_1_2_3_4_5_6_7_8_9_nc,GO_AllLists_nopblock[,c(1,3,4,10:12,15:16)])
GO_AllLists=GO_AllLists_1_2_3_4_5_6_7_8_9_nop
test=GO_AllLists[match(FINAL_GO$GO,GO_AllLists$GO),]
test$`#GeneInGOAndHitList` - test$`#GeneInGO`
#####################################

#save(GO_AllLists_1_2_3_4_5_6_7_8,GO_AllLists_1_2_3_4_5_6_7_8_9_nop,file = 'pathwayset.Rdata')
#save(test,file='top_pathwayset.Rdata')

##########
test %>% tidyr::separate_rows(GeneID,sep='[|]') %>% as.data.frame()->test.gmt.id
test %>% tidyr::separate_rows(Hits,sep='[|]') %>% as.data.frame()->test.gmt.symbol
#save(test.gmt.id,test.gmt.symbol,file = 'pathwayset.gmt.Rdata')
