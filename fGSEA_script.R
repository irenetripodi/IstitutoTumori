# 15/12/21
# 16/12/21

############################################
# Fast Gene Set Enrichment Analysis (fGESEA) 
############################################

# The package implements an algorithm for fast gene set enrichment analysis
#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("fgsea")

library(fgsea)

#install.packages("msigdbr")
library(msigdbr)
  
#install.packages("gplots")  # usefol to create plots (heatmap)
library(gplots)

#install.packages("RColorBrewer")  # tool to manage colors with R
library(RColorBrewer)

#----------------------------------------
# A brief inroduction to msigdbr package:
# Pathway analysis is a common task in genomics research and there are many available R-based software tools. 
# Depending on the tool, it may be necessary to import the pathways, translate genes to the appropriate species, convert between symbols and IDs, and format the resulting object.
# The Molecular Signatures Database (MSigDB) is one of the most widely used and comprehensive databases of gene sets for performing gene set enrichment analysis.
# Thus the msigdbr R package provides Molecular Signatures Database (MSigDB) gene sets typically used with the Gene Set Enrichment Analysis (GSEA) software.

# All gene sets in the database can be retrieved without specifying a collection/category.
# all_gene_sets = msigdbr(species = "Mus musculus")
# head(all_gene_sets) 
# There is a helper function to show the available species. 
# msigdbr_species()
# It is possible retrieve data for a specific collection, such as the hallmark gene sets.
# hallmark_gene_sets = msigdbr(species = "mouse", category = "H")
# head(hallmark_gene_sets)
#----------------------------------------

##############################
# 1) A20 IRES SPP1 Vs A20 IRES
##############################
  
data_1 <- read.table(file="A20_IRES_SPP1_Vs_A20_IRES.txt", sep="\t", header=T, quote="", dec=",")  # problem: the scientific notation e-10 are not considered 
#typeof(data_1) #it is a list
typeof(data_1$t) # the t column is composed by integers!
head(data_1$t)
  
data_1 <- data_1[order(data_1$t, decreasing = T),] # sort genes according to t-statistic and put in decreasing order.
# the other option is to sort genes according to LogFC
head(data_1$t)
  
# N.B.: the scientific notation e-10 is not considered, thus the order of the genes already performed is not correct!
data_1$t <- as.character(data_1$t) # convert the t column elements (that are integers) in characters 
typeof(data_1$t)
head(data_1$t)
  
data_1$t <- as.double(data_1$t) # convert the t column elements in double 
typeof(data_1$t)
head(data_1$t)
  
# Sort again the t-column of the data in decreasing order
data_1 <- data_1[order(data_1$t, decreasing = T),] 
head(data_1$t)
# now the order is correct!!

ranked_genes <- data_1$t # t-column ranked
head(ranked_genes)
names(ranked_genes) <- data_1$Symbol
head(names(ranked_genes))
  
hallmark <- msigdbr::msigdbr(species = "mouse", category = "H")
attributes(hallmark)
hallmark_list <- split(x = hallmark$gene_symbol, f = hallmark$gs_name)
head(hallmark_list)
  
?fgsea
gsea_1 <- fgsea(pathways = hallmark_list,
                  stats = ranked_genes,
                  nperm = 10000,
                  minSize = 15,
                  maxSize = 500)
  
# esxtract the NES and pvalue
NES_1 <- gsea_1$NES 
head(NES_1)
  
pval_1 <- gsea_1$pval
head(pval_1)
  
##############################
# 2) A20 IRES iOPN Vs A20 IRES
##############################

data_2 <- read.table(file="A20_IRES_iOPN_Vs_A20_IRES.txt", sep="\t", header=T, quote="",dec=",")  
#typeof(data_2) #it is a list
typeof(data_2$t) # the t column is composed by integers!
head(data_2$t)

data_2 <- data_2[order(data_2$t, decreasing = T),] # sort genes according to t-statistic and put in decreasing order.
# the other option is to sort genes according to LogFC
head(data_2$t)

# N.B.: the scientific notation e-10 is not considered, thus the order is not correct
data_2$t <- as.character(data_2$t) # convert the t column elements in characters 
typeof(data_2$t)
head(data_2$t)

data_2$t <- as.double(data_2$t)# convert the t column elements in double 
typeof(data_2$t)
head(data_2$t)

# Sort again the t-column in decreasing order
data_2 <- data_2[order(data_2$t, decreasing = T),] 
head(data_2$t)
# now the order is correct!!

ranked_genes <- data_2$t 
head(ranked_genes)
names(ranked_genes) <- data_2$Symbol
head(names(ranked_genes))

hallmark <- msigdbr::msigdbr(species = "mouse", category = "H")
attributes(hallmark)
hallmark_list <- split(x = hallmark$gene_symbol, f = hallmark$gs_name)
head(hallmark_list)

?fgsea
gsea_2 <- fgsea(pathways = hallmark_list,
            stats = ranked_genes,
            nperm = 10000,
            minSize = 15,
            maxSize = 500)

# esxtract the NES and pvalue
NES_2 <- gsea_2$NES 
head(NES_2)

pval_2 <- gsea_2$pval
head(pval_2)

###################################
# 3) A20 IRES SPP1 Vs A20 IRES iOPN
###################################

data_3 <- read.table(file="A20_IRES_SPP1_Vs_A20_IRES_iOPN.txt", sep="\t", header=T, quote="",dec=",")  
#typeof(data_3) #it is a list
typeof(data_3$t) # the t column is composed by integers!
head(data_3$t)

data_3 <- data_3[order(data_3$t, decreasing = T),] # sort genes according to t-statistic and put in decreasing order.
# the other option is to sort genes according to LogFC
head(data_3$t)

# N.B.: the scientific notation e-10 is not considered, thus the order is not correct
data_3$t <- as.character(data_3$t) # convert the t column elements in characters 
typeof(data_3$t)
head(data_3$t)

data_3$t <- as.double(data_3$t)# convert the t column elements in double 
typeof(data_3$t)
head(data_3$t)

# Sort again the t-column in decreasing order
data_3 <- data_3[order(data_3$t, decreasing = T),] 
head(data_3$t)
# now the order is correct!!

ranked_genes <- data_3$t 
head(ranked_genes)
names(ranked_genes) <- data_3$Symbol
head(names(ranked_genes))

hallmark <- msigdbr::msigdbr(species = "mouse", category = "H")
attributes(hallmark)
hallmark_list <- split(x = hallmark$gene_symbol, f = hallmark$gs_name)
head(hallmark_list)

?fgsea
gsea_3 <- fgsea(pathways = hallmark_list,
                stats = ranked_genes,
                nperm = 10000,
                minSize = 15,
                maxSize = 500)

# esxtract the NES and pvalue
NES_3 <- gsea_3$NES
head(NES_3)

hallmark_pathways_names <- gsea_3$pathway  # they are the names of the pahtways to be added in the heatmap, they are the same both for gsea_1, gsea_2 and gsea_3
head(hallmark_pathways_names)
View(hallmark_pathways_names)

pval_3 <- gsea_3$pval
head(pval_3)

#########
# HEATMAP
#########

# Create a matrix with the NES of the 3 groups of interest
all_NES <- cbind(NES_1, NES_2, NES_3)
head(all_NES)
typeof(all_NES) # it is a matrix composed by double

# rename the columns "NES_1", "NES_2" and "NES_3" with the names of the groups 
colnames(all_NES) <- c('A20 IRES SPP1 Vs A20 IRES','A20 IRES iOPN Vs A20 IRES','A20 IRES SPP1 Vs A20 IRES iOPN') 

# Create a matrix with the pvalues of the 3 groups of interest
all_pvalue <- cbind(pval_1, pval_2, pval_3)

# In the heatmap it is useful to display the values or "X" on all heatmap positions.
# cellnotes is the matrix of character strings which will be placed within each color cell, e.g. p-value symbols.

## come aveva fatto matteo:
## if p-value <= 0.05 (significant) leave the box of the heatmap empty, insted if is > 0.05 (not signficant) put a X in the box
## cellnotes <- ifelse(all_pvalue <= 0.05, "", "x") 

# N.B.: the heatmap should not contain rows that are not significant for all the 3 group (all the 3 boxes empty). 
# We want to remove the not significant rows and keep all the rows with a p_value <= 0.05 (significant)
all_pvalue_filtered <- all_pvalue[rowSums(all_pvalue[,1:3]<=0.05) != 0,] # 1:3 because there 3 columns
head(all_pvalue_filtered)
# filter the all_NES according to the p_values:
all_NES_filtered <- all_NES[rowSums(all_pvalue[,1:3]<=0.05) != 0,]  
head(all_NES_filtered)
# filter the names of the pathways according to the p_values:
hallmark_pathways_names_filtered <- hallmark_pathways_names[rowSums(all_pvalue[,1:3]<=0.05) != 0] # here we don't have to put the "," because it is not matrix but a vector of character of 1 column. 
head(hallmark_pathways_names_filtered)

# if p-value <= 0.05 (significant) put a X in the box, insted if is > 0.05 (not signficant) leave the box empty
cellnotes <- ifelse(all_pvalue_filtered <= 0.05, "x", "") 

mycol<- brewer.pal(11,"RdBu")
mycol <- mycol[11:1]

heatmap.2(all_NES_filtered,
          cellnote = cellnotes,
          notecol = "black",
          density.info="none",
          trace="none",
          col=mycol,
          Colv = F,
          Rowv = F,
          dendrogram="none",
          labRow = gsub("HALLMARK_","", hallmark_pathways_names_filtered), # substitute the first part of the string "HALLMARK_" with anything
          labCol = colnames(all_NES_filtered),
          cexRow=0.7,
          cexCol = 0.75,
          mar= c(10,18),
          key.xlab = "NES",
          key.title="NES",
          keysize = 1.1,
          lwid=c(0.5,2),
          main="HALLMARK",
          colsep=0:ncol(all_NES_filtered),
          rowsep=0:nrow(all_NES_filtered),
          sepcolor="black",
          sepwidth=c(0.000001,0.000001)
)

# Add the legend
## come aveva fatto matteo:
## legend("topright",legend = "X = pvalue > 0.05", bty="n",border=F,xpd=T)

legend("topright",legend = "X = pvalue <= 0.05",
       bty="n",border=F,xpd=T)


