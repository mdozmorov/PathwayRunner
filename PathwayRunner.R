# Fooling around with reactome.db functionality
xx <- as.list(reactomeEXTID2PATHID) # maps Entrez Gene identiﬁers to Reactome pathway identiﬁers
xx[1]
xx<-as.list(reactomeGO2REACTOMEID) # maps GO identiﬁers to Reactome database identiﬁers
xx[1]
xx<-as.list(reactomePATHID2EXTID) # maps Reactome pathway identiﬁers to Entrez Gene identiﬁers
xx[1]
xx<-as.list(reactomePATHID2NAME) # maps Reactome pathway identiﬁers to pathway names used by Reactome for various pathways
xx[1]
xx<-as.list(reactomePATHNAME2ID) # maps Reactome pathway names to pathway identiﬁers used by Reactome for various pathways
xx[1]
xx<-as.list(reactomeREACTOMEID2GO) # maps Reactome database identiﬁers to Gene Ontoloty (GO) identiﬁers
xx[1]
reactomeMAPCOUNTS # provides the "map count" (i.e. the count of mapped keys) for each map in package reactome.db

#
# PathwayRunner computed enrichment of gene set(s) in all pathways using hypergeometric test
# Gene sets' pathway enrichment profiles can then be pairwise correlated, 
# and correlation coefficients clustered and visualized, showing relationships among gene sets
# by their pathway connectivity
# 
# --------------------
# PathwayRunner begins
# --------------------

# ToDo: Add conversion of UniProt IDs (or any gene IDs) to EntreZ IDs, used by reactome.db

library('reactome.db')
pathway.count<-length(names(as.list(reactomePATHID2EXTID))) # Total number of pathways (9,381)
pathid.exid<-as.list(reactomePATHID2EXTID) # Each element is a pathID, containing list of genes in that pathway. Used to test input gene set(s) against each pathway in this pathway list
total.gene.count<-length(names(as.list(reactomeEXTID2PATHID))) # Total number of genes (29,492)

# Create list of input gene set(s), to be tested for the enrichments against each pathway
test.count<-15 # For testing purposes, simply use genes from first test.count pathways
xx <- as.list(reactomePATHID2EXTID) # maps Entrez Gene identiﬁers to Reactome pathway identiﬁers
gene.list<-list() # Initialize list containing gene set(s) to be tested for enrichment
for (i in 1:test.count) { # Populate this list with genes 
  gene.list[[length(gene.list)+1]] <- list(xx[[i]]) # Incrementally append gene sets to the list
}

# Initialize empty matrix, which will contain enrichment results of testing gene set(s) [columns]
# vs. pathways [rows]
mtx.pval<-matrix(numeric(0), test.count, length(gene.list))

# Now, the actual enrichment analysis
for (i in 1:length(gene.list)){ # Test gene sets one by one
  for (j in 1:test.count) { # Against all pathways
    pathway.length<-length(pathid.exid[[j]]) # Number of genes in the current pathway
    pathway.genes<-as.character(pathid.exid[[j]]) # Gene IDs in the current pathway
    genes.count<-length(unlist(gene.list[[i]])) # Number of input genes in the current set
    genes.in.pathway<-length(intersect(unlist(gene.list[[i]]),pathway.genes)) # Number of input genes in the current set present in the current pathway
    # 2 x 2 contingency table for hypergeometric test.
    # Example from http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
    #
    #           | Genes in pathway                  | Genes out of pathway                  |
    # --------------------------------------------------------------------------------------------------------
    # Input     | genes.in.pathway                  | genes.count - genes.in.pathway        | genes.count
    # --------------------------------------------------------------------------------------------------------
    # Non-input | pathway.length - genes.in.pathway | total.gene.count - genes.count        | total.gene.count
    #           |                                   | - (pathway.length - genes.in.pathway) | - genes.count
    # --------------------------------------------------------------------------------------------------------
    #           | pathway.length                    | total.gene.count - pathway.length     | total.gene.count  
    #
    p.value<-fisher.test(matrix(c(genes.in.pathway, pathway.length - genes.in.pathway,
                                  genes.count - genes.in.pathway, total.gene.count - genes.count - (pathway.length - genes.in.pathway)), 2, 2), alternative='two.sided')
    if (p.value$p.value < .Machine$double.xmin) p.value$p.value <- .Machine$double.xmin # If p-value is smaller than the minimum finite number, set it to that number to avoid Inf, when -log10-transforming
    mtx.pval[j, i] <- round(-log10(p.value$p.value), 5) # Store -log10-transformed p-value
    if (p.value$estimate < 1) mtx.pval[j, i] <- -mtx.pval[j, i] # If odds ratio < 1, then we have depletion, add "-" to highlight this fact
  }
}

# ToDo: Add correlation, clustering and visualization of pathway connectivity
