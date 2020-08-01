
#批处理需求：

library(xml2)
msigdb=read_xml('./msigdb_v7.1_files/msigdb_v7.1.xml')
genesetidtrans=data.frame()
for (i in 1:length(xml_contents(msigdb))){
  genesetid=data.frame(genesetID=xml_attrs(xml_child(msigdb, i))[["STANDARD_NAME"]],
                       geneset=xml_attrs(xml_child(msigdb, i))[["EXACT_SOURCE"]])
  genesetidtrans=rbind(genesetidtrans,genesetid)
}
#save(genesetidtrans,file = 'genesetidtrans.Rdata')

#批处理方案：
genesetidtrans=t(data.frame(row.names = c('STANDARD_NAME','CATEGORY_CODE','SYSTEMATIC_NAME',
                                          'EXACT_SOURCE','MEMBERS_EZID','MEMBERS_SYMBOLIZED')))
length(xml_contents(msigdb))
#msigdb=read_xml('./msigdb_v7.1_files/msigdb_v7.1.xml')
#save(msigdb,file = 'msigdb.Rdata')
as.character(xml_attrs(xml_child(msigdb, 1))[["STANDARD_NAME"]])
require(snowfall)
sfInit(parallel=TRUE, cpus=6)
sfLibrary(xml2)
sfExport('genesetidtrans')
genesetidtrans <- sfLapply(1:26860, function(x){#
  msigdb=read_xml('./msigdb_v7.1_files/msigdb_v7.1.xml')
  #msigdb=load(file = 'msigdb.Rdata')
  #as.character(xml_attrs(xml_child(msigdb, x))[["STANDARD_NAME"]])
                      c(as.character(xml_attrs(xml_child(msigdb, x))[["STANDARD_NAME"]]),
                       as.character(xml_attrs(xml_child(msigdb, x))[["CATEGORY_CODE"]]),
                       as.character(xml_attrs(xml_child(msigdb, x))[["SYSTEMATIC_NAME"]]),
                        as.character( xml_attrs(xml_child(msigdb, x))[["EXACT_SOURCE"]] ),
                       as.character( xml_attrs(xml_child(msigdb, x))[["MEMBERS_EZID"]] ),
                       as.character( xml_attrs(xml_child(msigdb, x))[["MEMBERS_SYMBOLIZED"]] )
                       )
}
)
sfStop()

genesetidtrans2=data.frame(STANDARD_NAME=as.character(unlist(unlist(rbind(genesetidtrans),recursive = F)[seq(1,26860*6,6)])),
                           CATEGORY_CODE=as.character(unlist(unlist(rbind(genesetidtrans),recursive = F)[seq(2,26860*6,6)])),
                           SYSTEMATIC_NAME=as.character(unlist(unlist(rbind(genesetidtrans),recursive = F)[seq(3,26860*6,6)])),
                           EXACT_SOURCE=as.character(unlist(unlist(rbind(genesetidtrans),recursive = F)[seq(4,26860*6,6)])),
                           MEMBERS_EZID=as.character(unlist(unlist(rbind(genesetidtrans),recursive = F)[seq(5,26860*6,6)])),
                           MEMBERS_SYMBOLIZED=as.character(unlist(unlist(rbind(genesetidtrans),recursive = F)[seq(6,26860*6,6)]))
                           )
strsplit(as.character(genesetidtrans2$STANDARD_NAME[1]),'[_]')[[1]]
genesetidtrans2$category=unlist(lapply(genesetidtrans2$STANDARD_NAME,function(x){
  unlist(strsplit(as.character(x),'_')[[1]])[1]
  }))
paste(unlist(strsplit(as.character(genesetidtrans2$STANDARD_NAME[1]),'_')[[1]])[-1],collapse = '_')
genesetidtrans2$term=unlist(lapply(genesetidtrans2$STANDARD_NAME,function(x){
  paste(unlist(strsplit(as.character(x),'_')[[1]])[-1],collapse = '_')
}))




#save(genesetidtrans,genesetidtrans2,file = "genesetidtrans.Rdata")
