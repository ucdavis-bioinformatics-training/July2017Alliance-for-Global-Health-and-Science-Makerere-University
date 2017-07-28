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

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

library(edgeR)

# read in counts table
counts <- read.table("counts_table.txt",as.is=T)
head(counts)

# Create object for use in limma and edgeR analysis
d <- DGEList(counts)

# calculate normalization factors
d <- calcNormFactors(d)
d
dim(d)

# filter genes below X cpms in all samples 
cutoff <- 1
drop <- which(apply(cpm(d), 1, max) < cutoff)
if (length(d) > 1) d <- d[-drop,] 
dim(d) # number of genes left

# generate experiment information
samples <- colnames(d)
pheno <- sapply(strsplit(samples,split="\\."), function(x) paste(x[1:(length(x)-1)],collapse="."))
pdata <- data.frame(samples,pheno)
pdata

# plot MDS
pdf("MDS.pdf")
plotMDS(counts, col = as.numeric(factor(pdata$pheno)), labels = pdata$samples)
dev.off()

###################### Limma-voom analysis 
mm <- model.matrix(~0 + pheno, data=pdata) # specify model with no intercept for easier contrasts
colnames(mm) <- sub("pheno","",colnames(mm))
mm
pdf("voom.pdf")
y <- voom(d, mm, plot = T)
mtext(side = 3, line = 0.5, text = "1-Factor Model Without Intercept")
dev.off()

fit <- lmFit(y, mm)
head(coef(fit))

# Comparison between cultivars C and I5 at time 6
contr <- makeContrasts(hst1.1.mac1 - hst1,hst1.1.mac1 - mac1, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, coef=1, sort.by = "P", n = Inf)
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes

write.table(tmp2, file = "hst1vhst1.1.mac1.txt", row.names = F, sep = "\t", quote = F, na = "")
head(tmp2, 20)


