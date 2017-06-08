# BiomaRt: converting gene identifiers
library(biomaRt)
mart<-useMart("ensembl", dataset="hsapiens_gene_ensembl")

name2Entrez<-function(names){ # Gene name to Entrez
  as.vector(as.matrix(getBM(attributes='entrezgene', filters='hgnc_symbol', values=names, mart=mart, uniqueRows=T)))
}

uniprot2Entrez<-function(names){ # UniProt to Entrez
  as.vector(as.matrix(getBM(attributes='entrezgene', filters='uniprot_swissprot_accession', values=names, mart=mart, uniqueRows=T)))
}
# 146                      uniprot_sptrembl               UniProt/TrEMBL Accession(s) [e.g. A2MYD1]
# 147                     uniprot_swissprot              UniProt/Swissprot ID(s) [e.g. GAGE4_HUMAN]
# 148           uniprot_swissprot_accession            UniProt/Swissprot Accession(s) [e.g. Q13068]
# 150                      uniprot_genename                      UniProt Genename ID(s) [e.g. V4-4]
# 152      uniprot_genename_transcript_name Uniprot Genename Transcript Name ID(s) [e.g. SEPT1-202]

Entrez2name<-function(names){ # Entrez to gene name
  as.vector(as.matrix(getBM(attributes='wikigene_name', filters='entrezgene', values=names, mart=mart, uniqueRows=T)))
}