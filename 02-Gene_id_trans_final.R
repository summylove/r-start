library(biomaRt)
#browseVignettes("") 
listMarts()
### BioMart databases can contain several datasets, for Ensembl
##the selected BioMart by using the function listDatasets().
ensembl=useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)
attributeslist=as.data.frame(attributes)

symbol <- getBM(attributes = c('entrezgene_accession','entrezgene_id','gene_biotype','hgnc_symbol'),
                mart = ensembl)
head(symbol)
dim(symbol) 

symbol2 <- symbol[-c(which(symbol$entrezgene_accession=="" & symbol$hgnc_symbol=="")),]  
dim(symbol2) 
#levels(factor(symbol2$gene_biotype))
symbol2=symbol2[which((symbol2$gene_biotype=="protein_coding")),]
dim(symbol2) 
symbol2=symbol2[!duplicated(symbol2$entrezgene_id),]
dim(symbol2) 
symbol2[which(symbol2$entrezgene_accession=="" ),
        'entrezgene_accession']=
  symbol2[which(symbol2$entrezgene_accession=="" ),'hgnc_symbol']
symbol2[which((symbol2$entrezgene_accession!=symbol2$hgnc_symbol)),1]=paste0(
  symbol2[which((symbol2$entrezgene_accession!=symbol2$hgnc_symbol)),1],"|",
  symbol2[which((symbol2$entrezgene_accession!=symbol2$hgnc_symbol)),4]
)
symbol2=symbol2[!duplicated(symbol2$entrezgene_accession),]
dim(symbol2) 

IDS=na.omit(symbol2[,c(1,2)])
dim(IDS)

save(IDS,file = "gencode.v31.annotation.ids.Rdata")


idtrans=function(genelist,f){
  # genelist=RMD_308_0
  load('gencode.v31.annotation.ids.Rdata')
  colnames(IDS)
  #colnames(IDS)
  if(f=='name'){
    print(paste0('mapping rate:',100*length(which(IDS$"entrezgene_accession" %in% genelist))/length(genelist),'%'))
    print('unmapped gene:')
    print(genelist[which(!( genelist %in% IDS$"entrezgene_accession"))])
    return(IDS[which(IDS$"entrezgene_accession" %in% genelist),"entrezgene_id"])
  }else{
    print(paste0('mapping rate:',100*length(which(IDS$entrezgene_id %in% genelist))/length(genelist),'%'))
    print('unmapped gene:')
    print(genelist[which(!( genelist %in% IDS$entrezgene_id))])
    return(IDS[which(IDS$entrezgene_id %in% genelist),"entrezgene_accession"])
  }
}
genelist <-idtrans(genelist,'name')
genelist <-idtrans(genelist,'id')


#### comparion ####
{
  
  symbol2[duplicated(symbol2$entrezgene_accession),]
  dim(symbol2) 
  unmap=symbol2[is.na(symbol2$entrezgene_id),'entrezgene_accession']
  unmap[!duplicated(unmap)]
  length(unmap)
}
{
 # BiocManager::install("EnsDb.Hsapiens.v86")  
  require(EnsDb.Hsapiens.v86)  
  keytypes(EnsDb.Hsapiens.v86)
  n=which(is.na(symbol2$entrezgene_id))
  
  
  V86= select(EnsDb.Hsapiens.v86,  
              keys =symbol2[n,'entrezgene_accession'] , 
              columns = "ENTREZID", 
              keytype = "SYMBOL"
  )
  colnames(V86)=c('entrezgene_accession','entrezgene_id')
  V0=symbol2[n,c('entrezgene_accession','entrezgene_id')]
  V0$entrezgene_id[(V0$entrezgene_accession %in% V86$entrezgene_accession)]=V86$entrezgene_id
  symbol2[n,'entrezgene_id']=V0$entrezgene_id
  
  unmap=symbol2[is.na(symbol2$entrezgene_id),'entrezgene_accession']
  unmap[!duplicated(unmap)]
  length(unmap)
} 
  
{
  ids=rtracklayer::import( "gencode.v31.annotation.gtf")
  idss=as.data.frame(mcols(ids))
  #[which(mcols(ids)$type=='gene'),])
  colnames(idss)
  head(idss)
  IDS=idss[,c(5,6,7,9)]

  IDS[1:5,]
  save(IDS,file = "gencode.v31.annotation.ids.Rdata")
}

{  
  require('rentrez')
  entrez_search()
  all_the_links <- entrez_link(dbfrom='gene', id=351, db="nuccore")
  entrez_fetch(all_the_links$links$gene_nuccore_refseqgene)
  
  require('annotate')
  browseVignettes('annotate')
  browseURL( entrezGeneByID(351))
  browseURL(entrezGeneQuery('APP'))

}