# 13/12/21

# Osteopontin (OPN) is a protein encoded by SPP1 gene (secreted phosphoprotein 1). 
# AIM: Analyze the gene expression of SPP1 in ABC (Activated B Cells) and GCB (Germinal Center B Cells) in DLBCL (Diffuse Large B Cell Lymphoma)

# From GEO DATASETS (a db that stores curated gene expression DataSets) has been provided the data of the expression profiling by array
# file:///Users/irene/Library/Containers/com.apple.mail/Data/Library/Mail%20Downloads/0968FDC9-3AAC-413F-8182-981CA0283698/staudt%20DLBCL%20-%20GEO%20DataSets%20-%20NCBI.html

# The first paper of interest is the following:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53786
# Accession: GSE53786
# Platform: GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# Affymetrix Affymetrix HG-U133_Plus_2 Array annotation data (chip hgu133plus2) assembled using data from public repositories

# In addition we downloaded the Series Matrix file in txt format (GSE53786_series_matrix.txt). 
# It is a file with some rows of header and then a huge matrix composed as following: 
# - in the rows all the probeID information with the rows counts
# - in the columns all the information of the classification in ABC and GCB broups and the raw count

# Since in the rows there are probeIDs a further step is need: matching the probeIDs with the names of the genes.
# Then selecting only the probeID of SPP1 gene and drowing a boxplot of the expression profiling of the gene between the ABC and GCB class.

######################################################
# 1 # DLBCL cell-of-origin by gene expression in FFPET
######################################################
unmodified_datafame <- read.table("GSE53786_series_matrix.txt", header=F, sep = "\t")

input_dataframe <- read.table("GSE53786_series_matrix_modified.txt", header=FALSE, sep= "\t") # dataframe without the header, composed by just the huge matrix with probeIDs
#class(input_dataframe) # it is a dataframe

### 1A) Matching the probe id ("ID_REF" column) present in input_dataframe with the gene name of interest (SPP1)
# To do that it is necessary to istall Bioconductor

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("hgu133plus2.db") # HG-U133_Plus_2 affymatrix is the format present in GEODatasets

library(AnnotationDbi)
library(hgu133plus2.db)

column_probe_id <- input_dataframe[,1]
class(column_probe_id) # they are factors. Keys argument wants character vectors
column_probe_id <- as.character(column_probe_id) # to convert factors into character vectors
class(column_probe_id)

annotation <- AnnotationDbi::select(hgu133plus2.db, keys = column_probe_id, columns = c("PROBEID", "SYMBOL"), keytype = "PROBEID") 
# key = colonna delle ProbeID; columns = annotazioni che vuoi; keytype = tipo  di annotazione che dai in ingresso

# warning message: there are more the one probeID matching with the gene name
View(annotation)
probe_id_for_SPP1_1 = "1568574_x_at"
probe_id_for_SPP1_2 = "209875_s_at"
# These are the probeID matched with the gene of interest. It is normal that there are more that one probeID for just 1 gene. 

#------------------
# N.B. An alternative to do this step is to go to https://www.ncbi.nlm.nih.gov/geo/ and search for the GEO accession (GSE53786)
# In the platform section select the number of the platform (GPL570)  --> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
# Download full table that matches the probeID with the gene names
# Then load the table and put as input the name of the gene of interest and the results will be the probeIDs.
#------------------

### 1B) Selecting the header of the dataframe carrying the information related to ABC and GCB cells

header_dataframe <- unmodified_datafame[43,] # the numebr of the raw with this information is 43
class(header_dataframe)

### 1C) Selecting the row(s) with the probeID of interest related to the gene of interest

row_PROBEID_1568574_x_at = input_dataframe[(input_dataframe$V1 == probe_id_for_SPP1_1),]  
row_PROBEID_209875_s_at = input_dataframe[(input_dataframe$V1 == probe_id_for_SPP1_2),]  
View(row_PROBEID_209875_s_at)

### 1D) Converting the rows of the dataframes into vectors, useful for the boxplots
class = as.vector(t(header_dataframe)) # as.vector function with t() to transpose the data frame
class = class[-1] # removing the first row not meaningful
View(class)

# We want to split the rownames of "class" vector in order to keep only the ABC or GCB information
new_class <- c()  # create an empty vector that will contain only the new classes 
for (element in class) {
  group <- strsplit(element, ": ")[[1]][3]  # split the string according to ":" and select only the 3rd element of the string (the group ABC or GCB)
  new_class <- append(new_class, group)  # append the group (ABC or GCB) in the empty vector
}
head(new_class)

PROBEID_1568574_x_at_vector = as.vector(t(row_PROBEID_1568574_x_at))
PROBEID_1568574_x_at_vector = PROBEID_1568574_x_at_vector[-1]
View(PROBEID_1568574_x_at_vector)
typeof(PROBEID_1568574_x_at_vector) # the vector is a character
PROBEID_1568574_x_at_vector = as.numeric(PROBEID_1568574_x_at_vector)
typeof(PROBEID_1568574_x_at_vector)

PROBEID_209875_s_at_vector = as.vector(t(row_PROBEID_209875_s_at))
PROBEID_209875_s_at_vector = PROBEID_209875_s_at_vector[-1]
View(PROBEID_209875_s_at_vector)
PROBEID_209875_s_at_vector = as.numeric(PROBEID_209875_s_at_vector)

### 1E) Prepare the data for the boxplots of the SPP1 gene expression for the two groups ABC vs GCB
typeof(PROBEID_209875_s_at_vector)
typeof(new_class)

# Create the dataframe composed by the raw counts of the first probeID in the rows and by the groups (ABC and GCB) in the columns
data_1 <- data.frame(PROBEID_209875_s_at_vector, new_class) 
str(data_1)
str(data_1$PROBEID_209875_s_at_vector)
dim(data_1)

# Repeat the same process also for the other probeID
data_2 <- data.frame(PROBEID_1568574_x_at_vector, new_class)
str(data_2)
str(data_2$PROBEID_1568574_x_at_vector)
dim(data_2)

# In the columns of both dataframes are also other two classes (not only ABC and GCB, but also just DLBL and Unclassified DLBCL) 
# Thus removing the rows of the columns DLBL and Unclassified DLBCL
data_1 = subset(data_1, new_class=="ABC DLBCL" | new_class=="GCB DLBCL" ) # considering only the two groups ABC and GCB in the class
View(data_1)
dim(data_1)     

data_2 = subset(data_2, new_class=="ABC DLBCL" | new_class=="GCB DLBCL" ) # considering only the two groups ABC and GCB in the class
View(data_2)
dim(data_2)     

# sorting the dataframe by the vector new_class
# Putting first all the rows with GCB present in the columns and then ABC in the columns
data_1 = data_1[
  with(data_1, order(data_1$new_class, decreasing = T)),
  ]
View(data_1)
dim(data_1) 

data_2 = data_2[
  with(data_2, order(data_2$new_class, decreasing = T)),
  ]
View(data_2)
dim(data_2) 

# The dataframe$new_classes is composed by 4 factor (reffering to the 4 groups - ABC, GCB, DLBL, Unclassified DLBCL - present before) 
# There was a warning since we have deleting the classes DLBL and Unclassified DLBCL and thus we should have just 2 factors and not 4
is.factor(data_1$new_class)
nlevels(data_1$new_class) # 4 factors: it is not correct
data_1$new_class <- droplevels(data_1$new_class)
nlevels(data_1$new_class) # now 2 factors: now correct

str(data_2)
is.factor(data_2$new_class)
nlevels(data_2$new_class) # 4 factors: it is not correct
data_2$new_class <- droplevels(data_2$new_class)
nlevels(data_2$new_class) # now 2 fctors: now correct

### 1F) Performing the t test between the two groups, to calculate the p-values
# t.test() performs one and two sample t-tests on vectors of data
GCB_group_dataset1 <- subset(data_1, new_class == "GCB DLBCL") # Return subsets of vectors for only GCB class
ABC_group_dataset1 <- subset(data_1, new_class != "GCB DLBCL") # Return subsets of vectors for only ABC class
ttest_res_dataset1 <- t.test(GCB_group_dataset1$PROBEID_209875_s_at_vector, ABC_group_dataset1$PROBEID_209875_s_at_vector)
pvalue_dataset1 <- ttest_res_dataset1$p.value

GCB_group_dataset2 <- subset(data_2, new_class == "GCB DLBCL")
ABC_group_dataset2 <- subset(data_2, new_class != "GCB DLBCL")
ttest_res_dataset2 <- t.test(ABC_group_dataset2$PROBEID_1568574_x_at_vector, ABC_group_dataset2$PROBEID_1568574_x_at_vector)
pvalue_dataset2 <- ttest_res_dataset2$p.value

### 1G) Perform the boxplots of the SPP1 gene expression for the two groups ABC vs GCB
# Create two boxplots for the 2 probeIDs
boxplot(data_1$PROBEID_209875_s_at_vector ~ data_1$new_class, col = c("green", "red"), main = "BOXPLOT: SPP1 expression", sub = "GEO accession: GSE53786; probeID: 209875_s_at; pvalue = 0.904", 
        names = c("GCB DLBCL", "ABC DLBCL"), ylab = "", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(data_2$PROBEID_1568574_x_at_vector ~ data_2$new_class, col = c("green", "red"), main = "BOXPLOT: SPP1 expression", sub = "GEO accession: GSE53786; probeID: 1568574_x_at_vector; pvalue = 1",
        names = c("GCB DLBCL", "ABC DLBCL"), ylab = "", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)







