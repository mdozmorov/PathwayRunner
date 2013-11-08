library(biomaRt)

# Archived biomart version for mm9 
mart=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset="mmusculus_gene_ensembl")
# or the latest biomart version for mm
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# or the latest for hs
mart<-useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Information - lists of filters and attributes
listDatasets(useMart("ensembl")) # List available datasets
filters<-listFilters(mart) # Source IDs to convert
head(filters)
filters[grep("ensembl",filters[,1], ignore.case=T),] # List specific ones
attr<-listAttributes(mart) # Destination IDs
head(attr)
attr[grep("symbol",attr[,1]),] # List gene symbol-specific

IDs<-readLines("clipboard") # Read a list of Ensembl IDs
# Convert from Ensembl IDs to Gene Name and Description
genes<-getBM(attributes=c('ensembl_gene_id','external_gene_id','description','chromosome_name','start_position','end_position','strand'), 
             filters='ensembl_gene_id', values=IDs, mart=mart)#, uniqueRows=T)
# Convert from UCSC IDs to Gene Name and Description
genes<-getBM(attributes=c('ucsc','external_gene_id','description'), 
             filters='ucsc', values=IDs, mart=mart)#, uniqueRows=T)
# Convert Gene names IDs (hgnc_symbol, mgi_symbol)
genes<-getBM(attributes=c('wikigene_name','external_gene_id','description', 'ensembl_gene_id'), 
             filters='wikigene_name', values=IDs, mart=mart)#, uniqueRows=T)
# Convert entrezgene
genes<-getBM(attributes=c('external_gene_id','description'), filters='entrezgene', values="1", mart=mart, uniqueRows=T)
# Convert RefSEq IDs
genes<-getBM(attributes=c('refseq_mrna','external_gene_id','description'), 
             filters='refseq_mrna', values=IDs, mart=mart)#, uniqueRows=T)

write.table(genes,"clipboard-128",sep="\t")
write.table(genes,"F:/111.txt",sep="\t")