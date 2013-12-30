# PathwayRunner computes enrichment of gene sets of interest in pre-defined gene sets using Fisher's exact test
# Example shows enrichment of 188 KEGG pathways (gene sets of interest) in 10,469 reactomeDb pathways (pre-defined gene sets)
# Prerequisite 1: Create gene.list, each element containing gene sets of interest. Can contain a single gene list. Example populates gene.list with gene sets defining KEGG pathways
# Prerequisite 2: Select pre-defined gene sets. Example uses gene sets defining reactomeDb pathways
source("utils.R") # functions for gene ID conversion

# Prerequisite 1: Create gene.list, each element containing gene sets of interest.
library(GSA)
geneset.obj<-GSA.read.gmt("data//msigdb.gmt") # Read all MSigDb gene sets
gene.list<-geneset.obj$genesets[grep("KEGG", geneset.obj$geneset.names)] # Get gene names for 186 KEGG pathways
gene.list<-lapply(gene.list, name2Entrez) # Convert gene names to EntrezIDs

# Prerequisite 2: Select pre-defined gene sets.
library('reactome.db')
pathway.count<-length(names(as.list(reactomePATHID2EXTID))) # Total number of pre-defined gene sets
pathid.exid<-as.list(reactomePATHID2EXTID) # Each element is a pathID, contains a list of genes in that pathway
total.gene.count<-length(names(as.list(reactomeEXTID2PATHID))) # Total number of genes

# mtx.pval will contain -log10-transformed enrichment p-values, with "-" added to indicate depletion
mtx.pval<-matrix(numeric(0), pathway.count, length(gene.list)) # Initialize empty [nxm] matrix, which will contain enrichment results of testing gene sets of interest [m] vs. pre-defined gene sets
colnames(mtx.pval)<-geneset.obj$geneset.names[grep("KEGG", geneset.obj$geneset.names)] # Set column names as the names of the gene sets of interest

Sys.time() # Starting time
# Now, the actual enrichment analysis
for (i in 1:length(gene.list)){ # Test gene sets of interest one by one
  for (j in 1:pathway.count) { # Against each pre-defined gene set
    pathway.length<-length(pathid.exid[[j]]) # Number of genes in the current pre-defined gene set
    pathway.genes<-as.character(pathid.exid[[j]]) # Gene IDs in the current  pre-defined gene set
    genes.count<-length(unlist(gene.list[[i]])) # Number of  genes in the current gene set of interest
    genes.in.pathway<-length(intersect(unlist(gene.list[[i]]), pathway.genes)) # Overlapping genes
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
Sys.time() # Ending time


# Clustering and visualization of pathway connectivity
mtx.cor<-cor(mtx.pval) #, use="pairwise.complete.obs")
library(gplots)
par(mar=c(10,6,6,5),oma=c(2,2,2,2)) #Make right and bottom margins larger
color<-colorRampPalette(c("blue","yellow"))
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"average" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
#Maximum distance seems to best separate correlations. Centroid linkage further sharpen the separation
h<-heatmap.2(as.matrix(mtx.cor),trace="none",density.info="none", col=color,distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)},cexRow=0.8,cexCol=0.8)

write.table(mtx.cor, "data//KEGG_vs_Reactome//KEGG_vs_Reactome.txt", sep="\t")

# # Experimenting with clustering/visualization with missing values
# # http://stat.ethz.ch/R-manual/R-devel/library/cluster/html/00Index.html
# library(cluster)
# plot(agnes(mtx.cor))
# library(ggplot2)
# x<-agnes(mtx.cor)
# y<-agnes(t(mtx.cor))
# # http://stats.stackexchange.com/questions/9050/how-to-display-a-matrix-of-correlations-with-missing-entries
# ggfluctuation(as.table(mtx.cor[x$order,y$order]), type="color") + labs(x="", y="") +
#   scale_fill_gradient(low = "blue",  high = "red")

# # Fooling around with reactome.db functionality
# xx <- as.list(reactomeEXTID2PATHID) # maps Entrez Gene identiﬁers to Reactome pathway identiﬁers
# xx[1]
# xx<-as.list(reactomeGO2REACTOMEID) # maps GO identiﬁers to Reactome database identiﬁers
# xx[1]
# xx<-as.list(reactomePATHID2EXTID) # maps Reactome pathway identiﬁers to Entrez Gene identiﬁers
# xx[1]
# xx<-as.list(reactomePATHID2NAME) # maps Reactome pathway identiﬁers to pathway names used by Reactome for various pathways
# xx[1]
# xx<-as.list(reactomePATHNAME2ID) # maps Reactome pathway names to pathway identiﬁers used by Reactome for various pathways
# xx[1]
# xx<-as.list(reactomeREACTOMEID2GO) # maps Reactome database identiﬁers to Gene Ontoloty (GO) identiﬁers
# xx[1]
# reactomeMAPCOUNTS # provides the "map count" (i.e. the count of mapped keys) for each map in package reactome.db
