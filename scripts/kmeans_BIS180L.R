install.packages("MBCluster.Seq")
library("MBCluster.Seq")
data("Count")
head(Count)
dim(Count)

geneID <- 1:nrow(Count)
norm   <- rep(1, ncol(Count))
Treat  <- rep(1:4, 2)

mydata <- RNASeq.Data(Count, Normalize = norm, Treat, geneID)
str(mydata)

c0 <- KmeansPlus.RNASeq(data = mydata, nK = 10)$centers
c0
cls <- Cluster.RNASeq(data = mydata, model = "nbinom", centers = c0, method = "EM")$cluster
tree <- Hybrid.Tree(data = mydata, cluster = cls, model = "nbinom")
plotHybrid.Tree(merge = tree, cluster = cls, logFC = mydata$logFC, tree.title = NULL)

head(tree)


# this is the data workflow that the students did last week
# I will reproduce this analysis really quick so I can use the output
# for clustering

setwd("/Users/Cody_2/git.repos/BIS180L/data")
counts.data <- read.table("gh_internode_counts.tsv")
dim(counts.data)
counts.data <- counts.data[rownames(counts.data)!="*",]
counts.data[is.na(counts.data)] <- 0
counts.data <- counts.data[rowSums(counts.data > 10) >= 3,]


sample.description <- data.frame(
  sample=colnames(counts.data),

  #This next line searches for IMB211 or R500 in the colnames of counts.data and returns anything that matches
  #In this way we can extract the genotype info.
  gt=regmatches(colnames(counts.data),regexpr("R500|IMB211",colnames(counts.data))),

  #Now we use the same method to get the treatment
  trt=regmatches(colnames(counts.data),regexpr("NDP|DP",colnames(counts.data)))
  )
sample.description


# Now we can paste the trt and gt columns together to give a group identifier
sample.description$group <- paste(sample.description$gt,sample.description$trt,sep="_")

# set the reference treatment to "NDP"
sample.description$trt <- relevel(sample.description$trt,ref="NDP")

sample.description

library(edgeR)
dge.data <- DGEList(counts=counts.data, group=sample.description$group)
dim(dge.data) 
dge.data <- calcNormFactors(dge.data, method = "TMM")
dge.data$samples #
plotMDS(dge.data, method = "bcv") 
counts.data.normal <- cpm(dge.data)
counts.data.log <- log2(counts.data + 1)
design <- model.matrix(~gt+trt,data = sample.description)
rownames(design) <- sample.description$sample
design

#First the overall dispersion
dge.data <- estimateGLMCommonDisp(dge.data,design,verbose = TRUE)

#Then a trended dispersion based on count level
dge.data <- estimateGLMTrendedDisp(dge.data,design)

#And lastly we calculate the gene-wise dispersion, using the prior estimates to "squeeze" the dispersion towards the common dispersion.
dge.data <- estimateGLMTagwiseDisp(dge.data,design)

#We can examine this with a plot
plotBCV(dge.data)
fit <- glmFit(dge.data, design)
fit
gt.lrt <- glmLRT(fit,coef = "gtR500")
topTags(gt.lrt) # the top 10 most differentially expressed genes
summary(decideTestsDGE(gt.lrt,p.value=0.01)) #This uses the FDR.  0.05 would be OK also.

DEgene.gt <- topTags(gt.lrt,n = Inf)$table[topTags(gt.lrt,n = Inf)$table$FDR<0.01,]
write.csv(DEgene.gt,"DEgenes.gt.csv")
plotDE <- function(genes, dge, sample.description) {
  require(ggplot2)
  require(reshape2)
  tmp.data <- t(log2(cpm(dge[genes,])+1))
  tmp.data <- merge(tmp.data,sample.description,by.x="row.names",by.y="sample")
  tmp.data <- melt(tmp.data,value.name="log2_cpm",variable.name="gene")
  pl <- ggplot(tmp.data,aes(x=gt,y=log2_cpm,fill=trt))
  pl <- pl + facet_wrap( ~ gene)
  pl <- pl + ylab("log2(cpm)") + xlab("genotype")
  pl <- pl + geom_boxplot()
  pl + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
}

plotDE("Bra009785",dge.data,sample.description)
plotDE(rownames(DEgene.gt)[1:9],dge.data,sample.description)

design.interaction <- model.matrix(~gt*trt,data = sample.description)
rownames(design.interaction) <- sample.description$sample
design.interaction

#First the overall dispersion
interaction <- dge.data
interaction <- estimateGLMCommonDisp(interaction,design.interaction,verbose = TRUE)

#Then a trended dispersion based on count level
interaction <- estimateGLMTrendedDisp(interaction,design.interaction)

#And lastly we calculate the gene-wise dispersion, using the prior estimates to "squeeze" the dispersion towards the common dispersion.
interaction <- estimateGLMTagwiseDisp(interaction,design.interaction)

#We can examine this with a plot
plotBCV(interaction)
fit2 <- glmFit(interaction, design.interaction)
fit2
gt.lrt2 <- glmLRT(fit2,coef = "gtR500:trtDP")
topTags(gt.lrt2) # the top 10 most differentially expressed genes
summary(decideTestsDGE(gt.lrt2, p.value=0.05)) 


DEgene.GxE <- topTags(gt.lrt2,n = Inf)$table[topTags(gt.lrt2,n = Inf)$table$FDR<0.01,]
write.csv(DEgene.GxE,"DEgenes_GxE.csv")

head(DEgene.GxE)
dim(DEgene.GxE)

de_gene_names <- row.names(DEgene.GxE)
de_gene_names

dim(counts.data.normal)
head(counts.data.normal)

row.names(counts.data.normal) == de_gene_names

# subset normalized counts dataframe for clustering
GxE_counts <- as.data.frame(counts.data.normal[de_gene_names,])
dim(GxE_counts)
str(GxE_counts)

hclust(GxE_counts)
rownames(GxE_counts)

?hclust
hc <- hclust(dist(counts.data), "ave")
plot(hc)

?kmeans
(cl <- kmeans(counts.data, 10))
plot(counts.data, col = cl$cluster)

GxEnorm <- normalize(counts.data)
outsom <- som(GxEnorm, xdim = 5, ydim = 5)
plot(outsom)
str(outsom)


plot(outsom)
library("MBCluster.Seq")
counts.data
GeneID <- row.names(counts.data)
Norm   <- rep(1, ncol(counts.data))
Treat  <- rep(1:3, 4)
rna_out <- RNASeq.Data(counts.data, Normalize = NULL, Treat, GeneID)

c0 <- KmeansPlus.RNASeq(data = rna_out, nK = 10)$centers
cls <- Cluster.RNASeq(data = rna_out, model = "nbinom", centers = c0, method = "EM")$cluster
cls
tree <- Hybrid.Tree(data = rna_out, cluster = cls)
 
dim(counts.data)

corrmatrix <- cor(GxE_counts)
head(corrmatrix)


#############
#trouble shooting
#############
?prcomp
prcomp_counts <- prcomp(t(GxE_counts))
str(prcomp_counts)
plot(prcomp(GxE_counts))
head(prcomp_counts$rotation)


# load ggplot2
library(ggplot2)

# create data frame with scores
scores <- as.data.frame(prcomp_counts$rotation)[,c(1,2)]
head(scores)
rownames(scores)
rownames(clus)
str(scores)
str(clus)
plotting <- merge(clus, scores, by = "row.names")
str(plotting)
plotting$clus <- as.factor(plotting$clus)
clus
# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = clus)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 4, stat = "identity") 


d <- qplot(PC1, PC2, data=plotting, colour=clus)
d + scale_colour_brewer()


str(fit)
plot(fit$centers)
clus <- as.data.frame(fit$cluster)
str(clus)
names(clus) <- paste("cluster")
head(clus)
setwd("/Users/Cody_2/Box Sync/brassica_ein")
setwd("/Users/Cody_2/git.repos/BIS180L_web/data/")
genes2 <- read.table("GH_merged_v1.5_mapping.tsv", sep = "\t", header = TRUE)
head(genes2)




```{r, eval = FALSE}


#found the culpret libraries driving that weird response
genes <- read.table("voom_transform_brassica.csv", sep = ",", header = TRUE)
DE_genes <- read.table("DEgenes_GxE.csv", sep = ",")
DE_gene_names <- rownames(DE_genes)


head(genes)
lib_names <- names(genes)
lib_names
#drop these columns
# 38, 42, 46
genes <- genes[,-c(38,42,46)] 

GxE_counts <- as.data.frame(genes[DE_gene_names,])
head(GxE_counts)
genes_cor <- cor(t(GxE_counts))
dim(genes_cor)
min(genes_cor)
max(genes_cor)
mean(genes_cor)
median(genes_cor)


hist(genes_cor[upper.tri(genes_cor)])
genes_cor
head(genes_cor)[,1:10]
head(genes_cor)[1:10]
tail(genes_cor)[,1:10]
adjcent.z=abs(genes_cor) > 0.80
diag(adjcent.z)=0  ## genes do not connect themselves in the network
rownames(adjcent.z)=rownames(genes_cor)
colnames(adjcent.z)=colnames(genes_cor)
sum(adjcent.z)/2 ##  1251







#####yasu code
### fishers-z transformation (to make the sample correlation more comparable), 


n=48 ### this is the sample size 
z.s <- sqrt(n-3)*0.5*log((1+genes_cor)/(1-genes_cor))
summary(z.s[upper.tri(z.s)])
hist(z.s[upper.tri(z.s)])
head(z.s)
?qnorm



?edge.connectivity
?barabasi.game
g <- barabasi.game(20, m=1)
plot(g)
g2 <- barabasi.game(20, m=5)
edge.connectivity(g2)
plot(g2)

edge.connectivity(g, 100, 1)
edge.connectivity(g2, 100, 1)
edge.disjoint.paths(g2, 100, 1)

g <- erdos.renyi.game(50, 5/50)
g <- as.directed(g)
g <- induced.subgraph(g, subcomponent(g, 1))
graph.adhesion(g)
plot(g)



### cut off 
thre.z = qnorm(0.999999999999999995)  ## normal quanitle 
thre.z
adjcent.z=abs(z.s)>11  ## symmetric ajacency matrix: 1: there is an edge; 0 : there is no edge 
diag(adjcent.z)=0  ## genes do not connect themselves in the network
rownames(adjcent.z)=rownames(genes_cor)
colnames(adjcent.z)=colnames(genes_cor)
sum(adjcent.z)/2 ## 2155 edges
index=rowSums(adjcent.z)>0
weight.adjcent.z=adjcent.z[index,index]
head(weight.adjcent.z)

library(igraph)
g.temp=graph.adjacency(weight.adjcent.z, mode="undirected", diag=FALSE)
plot(g.temp)
community.leading=leading.eigenvector.community(g.temp)
community.leading
```

