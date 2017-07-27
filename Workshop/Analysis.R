### First lets convert the individual salmon counts to a counts table
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("rjson","tximport"))

setwd("~/RNAseq_Workshop")
library("tximport")
library("rjson")

# Read in each of the quant.sf files produced by salmon
txi <- tximport(dir(path = "02-Salmon",full.names=T,recursive=T,pattern="quant.sf"),
                type = "salmon",txOut=TRUE,countsFromAbundance = "lengthScaledTPM")

# rename the columns to be our sample names
nms <- sub(".quant","",sapply(strsplit(dir(path="02-Salmon"),split="_"),"[[",1L))
colnames(txi$counts) = nms

# view the first 6 rows of the table
head(txi$counts)

# save the table
write.table(txi$counts,"counts_table.txt",sep="\t",quote=F)
