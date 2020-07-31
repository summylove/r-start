
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
genesetidtrans=data.frame()
length(xml_contents(msigdb))
#save(msigdb,file = 'msigdb.Rdata')
require(snowfall)
sfInit(parallel=TRUE, cpus=6)
sfLibrary(xml2)
sfExport('genesetidtrans')
genesetidtrans <- sfLapply(1:26860, function(x){
  msigdb=read_xml('./msigdb_v7.1_files/msigdb_v7.1.xml')
  genesetidtrans[x,]=data.frame(genesetID=xml_attrs(xml_child(msigdb, x))[["STANDARD_NAME"]],
                       geneset=xml_attrs(xml_child(msigdb, x))[["EXACT_SOURCE"]])[1,]
}
)
sfStop()
save(genesetidtrans,file = "genesetidtrans.Rdata")