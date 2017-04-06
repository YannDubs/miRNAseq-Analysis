# Packages
library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)
library(RColorBrewer)
library(limma)
library(gplots)
library(edgeR)
library(devtools)
library(ggbiplot)
library(DESeq)
library(geneplotter)
library(Rtsne)
library(NMF)
library(gridExtra)

set.seed(1234) # makes everything reproductible

meta_data <- read.table("../Data/raw_data/meta_data.txt", header = TRUE)
dataInitial <- read.table("../Data/raw_data/expn_matrix.txt", header = TRUE, row.names = 1)
spike_data <- read.table("../Data/raw_data/expn_matrix_spike.txt", header = TRUE, row.names = 1)

### NORMALIZATION ###

# finding sum of column and adding it as new row
spike_data <- rbind(spike_data, "Totals" = colSums(spike_data))
scaled_totals <- spike_data["Totals", ]
average_reads <- rowMeans(scaled_totals)

# NORMALIZING BY SPIKE-RNA
scaled_totals <- rbind(scaled_totals, "Scale" = scaled_totals[1, ] / average_reads)

# setting dim_names for upcoming matrix
dim_names <- list(row.names(dataInitial))
dim_names[2] <- list(colnames((dataInitial)))

# pre-making a matrix
normalized_data <- matrix(nrow = nrow(dataInitial), ncol = ncol(dataInitial), dimnames = dim_names)

# for loop which scales data by the spike miRNA and saves it as normalized_data
for(i in 1: ncol(dataInitial)){
    normalized_data[ , i] <- dataInitial[ , i] / (scaled_totals[2, i])
}
normalized_data <- as.data.frame(normalized_data)

total_data <- rbind(spike_data, "Total Spike" = colSums(spike_data), "Total Reads Norm" = colSums(normalized_data), "Total Reads Init" = colSums(dataInitial))

### KEEEP???
if (FALSE){
    total_cellular <- total_data["Total Reads Norm", ]
    average_reads <- rowMeans(total_cellular)
    cell_scales <- rbind(total_cellular, "Scale" = total_cellular[1, ] / average_reads) #dividing total spike values by average spike values
    dim_names <- list(row.names(dataInitial))       #setting dim_names for upcoming matrix
    dim_names[2] <- list(colnames((dataInitial)))
    normalized_data2 <- matrix(nrow = nrow(dataInitial), ncol = ncol(dataInitial), dimnames = dim_names) #pre-making a matrix
    for(i in 1: ncol(dataInitial)){   #for loop which scales data and saves it as normalized_data
        normalized_data2[ , i] <- normalized_data[ , i] / (cell_scales[2, i])
    }
    normalized_data <- as.data.frame(normalized_data2)
}
####

# DELETING DEAD CELLS
alive_indexes <- list() #initializing lists
k <- 1  #initializing counter
for(i in 1:ncol(total_data)){
    if(total_data["Total Reads Init", i] > 20000){
        alive_indexes[k] <- i
        k <- k + 1
    }
}

alive_data <- normalized_data[ , as.numeric(alive_indexes)]
alive_meta_data <- meta_data[as.numeric(alive_indexes), ]

# saving completely normalized data as n_data and n_meta_data
n_data <- alive_data
n_meta_data <- alive_meta_data

# makes group and design matrix
data_day <- as.character(n_meta_data$Population)
data_day <- recode(data_day, "HL60D0" = "Day 0", "HL60D1" = "Day 1", "HL60D3" = "Day 3", "HL60D7" = "Day 7")
data_day <- factor(data_day)

design <- model.matrix(~0 + n_meta_data$Population)
# deletes "n_meta_data$Population" from name of columns
colnames(design) <- gsub("n_meta_data\\$PopulationHL60", "", colnames(design))

# Keeps genes without all zeros miRNA
dataNZ <- n_data[which(rowSums(n_data) > 0),]

# making the contrast matrix for comparing each group
contrastMatrix <- makeContrasts(D1-D0,
D3-D0,
D7-D0,
D1-D3,
D1-D7,
D3-D7,
levels = c("D0","D1","D3","D7"))

# Removes unnecessary big data
remove(normalized_data)
remove(n_data)
remove(alive_data)
remove(alive_meta_data)
remove(dataInitial)

### DATA EXPLORATION ###

## PCA ##

#Log transform
log_normalized_dataNZ <- log(dataNZ)
log_normalized_dataNZ[log_normalized_dataNZ == "-Inf"] <- 0

data_pcaNZ <- prcomp(t(log_normalized_dataNZ), center = TRUE, scale. = TRUE)
ggbiplot(data_pcaNZ,
scale = 1,
obs.scale = 1,
varname.abbrev = FALSE,
var.axes = FALSE,
pc.biplot =TRUE,
circle = TRUE,
groups = data_day,
ellipse= TRUE) +
ggtitle("PCA with standardization") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

# will not be needed: don't keep big data
remove(data_pcaNZ)
remove(log_normalized_dataNZ)

# log transform
log_normalized_data <- log(dataNZ)
log_normalized_data[log_normalized_data == "-Inf"] <- 0

data_pca <- prcomp(t(log_normalized_data))
ggbiplot(data_pca,
scale = 1,
obs.scale = 1,
varname.abbrev = FALSE,
var.axes = FALSE,
pc.biplot =TRUE,
circle = TRUE,
groups = data_day,
ellipse= TRUE) +
ggtitle("PCA without standardization") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

# will not be needed: don't keep big data
remove(data_pca)
remove(log_normalized_data)

# arranging data based on total expression
indexTopFifty <- sort(rowSums(dataNZ), index=T, decreasing=TRUE)$ix[1:50]
topFifty <- dataNZ[indexTopFifty,]

# Log transform
log_topFifty <- log(topFifty)
log_topFifty[log_topFifty == "-Inf"] <- 0

# PCA
data_topFifty <- prcomp(t(log_topFifty))
ggbiplot(data_topFifty,
scale = 1,
obs.scale = 1,
varname.abbrev = FALSE,
var.axes = FALSE,
pc.biplot =TRUE,
circle = TRUE,
groups = data_day,
ellipse= TRUE) +
ggtitle("PCA top 50 genes") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

# will not be needed: don't keep big data
remove(data_topFifty)
remove(log_topFifty)
remove(topFifty)

## CORRELATION HEATMAP ##

# computing sample-sample correlations
normalized_data_cor <- cor(dataNZ)

cor_dat <- as.data.frame(normalized_data_cor)

labs <- list(Population = n_meta_data$Population)

NMF::aheatmap(cor_dat,
annCol = labs, annRow = labs, Rowv = NA, Colv = NA, labRow = NA, labCol = NA,
main = "Sample-Sample Correlation Matrix")

#will not be needed: don't keep big data
remove(normalized_data_cor)

## T-SNE ##

# ordering data by expn
alive_tsne_data <- cbind(dataNZ, "avg_expn" = rowMeans(dataNZ), "miRNA" = row.names(dataNZ))
alive_tsne_data <- dplyr::arrange(alive_tsne_data,desc(avg_expn))
alive_tsne_data <- t(alive_tsne_data)
colnames(alive_tsne_data) <- alive_tsne_data["miRNA", ]
#removes rows that were added for ordering
alive_tsne_data <- alive_tsne_data[!rownames(alive_tsne_data) %in% c("miRNA","avg_expn"), ]

# Runs tsne
tsne_out <- Rtsne(dist(alive_tsne_data))

# Show the objects in the 2D tsne representation
qplot(x=tsne_out$Y[,1], y=tsne_out$Y[,2],color=data_day) +
labs(colour = "Cell Type", x = "tsne1", y = "tsne2")+
ggtitle("t-SNE alive cells and non zero genes") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

#will not be needed: don't keep big data
remove(tsne_out)

# ordering data by expn
alive_tsne_data <- cbind(dataNZ, "avg_expn" = rowMeans(dataNZ), "miRNA" = row.names(dataNZ))
alive_tsne_data <- dplyr::arrange(alive_tsne_data,desc(avg_expn))
alive_tsne_data <- t(alive_tsne_data)
colnames(alive_tsne_data) <- alive_tsne_data["miRNA", ]
#removes rows that were added for ordering
alive_tsne_data <- alive_tsne_data[!rownames(alive_tsne_data) %in% c("miRNA","avg_expn"), ]

# Runs tsne
tsne_out <- Rtsne(dist(alive_tsne_data))

# Show the objects in the 2D tsne representation
qplot(x=tsne_out$Y[,1], y=tsne_out$Y[,2],color=data_day) +
labs(colour = "Cell Type", x = "tsne1", y = "tsne2")+
ggtitle("t-SNE alive cells and non zero genes") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

#will not be needed: don't keep big data
remove(tsne_out)

# tsne on Day 0/Day 1 cells only
indexDay01 <- which(data_day %in% c("Day 0","Day 1"))
# selecting only day 0 and day 1 cells
group_one_data <- alive_tsne_data[indexDay01, ]

# Run TSNE
tsne_out <- Rtsne(dist(group_one_data))

# Show the objects in the 2D tsne representation
qplot(x=tsne_out$Y[,1], y=tsne_out$Y[,2],color=data_day[indexDay01]) +
labs(colour = "Cell Type", x = "tsne1", y = "tsne2") +
ggtitle("t-SNE only day 0 and 1") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

#will not be needed: don't keep big data
remove(tsne_out)

tsne_out <- Rtsne(dist(group_one_data), dims = 3) # Run TSNE in 3 dimensions

plot_ly(x=tsne_out$Y[,1],
y=tsne_out$Y[,2],
z=tsne_out$Y[,3],
type="scatter3d",
color=data_day[indexDay01],
mode="markers") %>%
layout(title = 't-SNE day 0 and 1 3D')

#will not be needed: don't keep big data
remove(tsne_out)
remove(group_one_data)

#tsne on Day 3/Day 7 cells only
indexDay37 <- which(data_day %in% c("Day 3","Day 7"))
#selecting only day 3 and day 7 cells
group_two_data <- alive_tsne_data[indexDay37, ]

# Run TSNE
tsne_out <- Rtsne(dist(group_two_data))
# Show the objects in the 2D tsne representation
qplot(x=tsne_out$Y[,1], y=tsne_out$Y[,2],color=data_day[indexDay37]) +
labs(colour = "Cell Type", x = "tsne1", y = "tsne2") +
ggtitle("t-SNE only day 3 and 7") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

#will not be needed: don't keep big data
remove(tsne_out)
remove(alive_tsne_data)

tsne_out <- Rtsne(dist(group_two_data), dims = 3) # Run TSNE in 3 dimensions

plot_ly(x=tsne_out$Y[,1],
y=tsne_out$Y[,2],
z=tsne_out$Y[,3],
type="scatter3d",
color=data_day[indexDay37],
mode="markers") %>%
layout(title = 't-SNE day 3 and 7 3D')

#will not be needed: don't keep big data
remove(tsne_out)
remove(group_two_data)

### DATA ANALYSIS ###

## MEAN-VARIANCE ##

dge <- DGEList(counts=dataNZ, group = data_day)

# applies TMM normalization to dge
#dge <- calcNormFactors(dge)
data_voomed <- voom(dge,design,plot=TRUE)

#will not be needed: don't keep big data
remove(dge)

df <- data_voomed$E
dataf<- mutate(data.frame(t(df)), group = data_day) %>%
gather(mirna, exp, - group)

ggplot(dataf, aes(x=exp)) +
geom_histogram(aes(fill = group), colour = 'white')+
ggtitle("Distribution of counts with respect to the groups") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

#will not be needed: don't keep big data
remove(dataf)
remove(df)

# filtering out rows that have non zero counts in less than 30% of the samples
dataNZFilterd <-dataNZ[rowSums(dataNZ>=10)>=round(0.3*ncol(dataNZ)),]

dge <- DGEList(counts=dataNZFilterd, group = data_day)

#dge <- dge[isExpr,]
#dge <- calcNormFactors(dge)
data_voomed <- voom(dge,design,plot=TRUE)

#will not be needed: don't keep big data
remove(dge)
remove(dataNZFilterd)

data_voomed[1:2,1:10]

indexDay0 <- which(data_day == "Day 0")
indexDay1 <- which(data_day == "Day 1")
indexDay3 <- which(data_day == "Day 3")
indexDay7 <- which(data_day == "Day 7")
dataNZFilterdBis <- dataNZ[rowSums(dataNZ[,indexDay0]>=10)>=round(0.50*length(indexDay0)) |
rowSums(dataNZ[,indexDay1]>=10)>=round(0.50*length(indexDay1)) |
rowSums(dataNZ[,indexDay3]>=10)>=round(0.50*length(indexDay3)) |
rowSums(dataNZ[,indexDay7]>=10)>=round(0.50*length(indexDay7)) ,]

dge <- DGEList(counts=dataNZFilterdBis, group = data_day)

#dge <- calcNormFactors(dge)
data_voomed <- voom(dge,design,plot=TRUE)

df <- data_voomed$E
dataf<- mutate(data.frame(t(df)), group = data_day) %>%
gather(mirna, exp, - group)

ggplot(dataf, aes(x=exp)) +
geom_histogram(aes(fill = group), colour = 'white') +
ggtitle("Distribution of filtered counts with respect to the groups") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))


#will not be needed: don't keep big data
remove(dataf)
remove(df)

data_group <- as.character(data_day)
data_group[data_group == "Day 0" | data_group == "Day 1"] <- "Day01"
data_group[data_group == "Day 3" | data_group == "Day 7"] <- "Day37"

#makes a data frame with all information
df <- data_voomed$E
dataf4 <- data.frame(df)
dataf4$mirna <- as.factor(rownames(dataf4))
dataf4 <- gather(dataf4,samples, exp,-mirna);
dataf4$group <-  as.factor(data_group);
dataf4$day <-  data_day;
dataf4 <- mutate(dataf4, interGroup = paste(as.character(dataf4$mirna),data_group) )
dataf4 <- mutate(dataf4, interGroup2 = paste(as.character(dataf4$mirna),as.character(data_day)) )
dataf4$exp2 <- 2^dataf4$exp-0.5
randRNA <- unique(dataf4$mirna)[sample(1:length(unique(dataf4$mirna)), 12, replace=FALSE)]

ggplot(dataf4[dataf4$mirna %in% randRNA,], aes(x=exp, fill = group)) +
geom_density(alpha=0.5) +
ggtitle("Random density plot for interaction miRNA-group") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5)) +
theme(legend.position="none") + xlim(0,30) +
facet_wrap(~mirna,nrow = 4,ncol = 3)

#will not be needed: don't keep big data
remove(df)

ggplot(dataf4[dataf4$mirna %in% randRNA,], aes(x=exp2, fill = group)) +
geom_density(alpha=0.5) +
ggtitle("Random density plot for interaction miRNA-group: no log") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5)) +
theme(legend.position="none") + xlim(-500,10000) +
facet_wrap(~mirna,nrow = 5,ncol = 3)

### FIT LINEAR MODEL ###

# Convert population to a numeric variable.
Pop.num<-as.numeric(substr(n_meta_data$Population,6,6))
n_meta_data<-cbind(n_meta_data,Pop.num)

des.pop.num<-model.matrix(~Pop.num,n_meta_data)
v.pop.num<-voom(dge,des.pop.num, plot=FALSE)

# Next, fit linear model with day as a factor.
fit.pop.num <- lmFit(v.pop.num,des.pop.num)
efit.pop.num <- eBayes(fit.pop.num)

# store result
top.pop.num<-topTable(efit.pop.num,number=Inf)

# Look at table with top genes.
top.pop.num[1:10,]

fdr.pop.num<-nrow(subset(top.pop.num,adj.P.Val<0.01))

remove(v.pop.num)

hist(top.pop.num$adj.P.Val,breaks=29,main="Q-value Histogram - Day as numeric",xlab="Q-Value")
hist(top.pop.num$P.Value,breaks=29,main="P-value Histogram - Day as numeric",xlab="P-Value")

remove(top.pop.num)
des.pop.char<-model.matrix(~Population,n_meta_data)
v.pop.char<-voom(dge,des.pop.char, plot=FALSE)
fit.pop.char <- lmFit(v.pop.char,des.pop.char)
efit.pop.char <- eBayes(fit.pop.char)

top.pop.char<-topTable(efit.pop.char,number=Inf)
plotSA(efit.pop.char,main="Final Mean-Variance plot")

#will not be needed: don't keep big data
remove(fit)
remove(v.pop.char)

top.pop.char[1:10,]
fdr.pop.char<-nrow(subset(top.pop.char,adj.P.Val<0.01))
hist(top.pop.char$adj.P.Val,breaks=29,main="Q-value Histogram - Population as a factor",xlab="Q-Value")
hist(top.pop.char$P.Value,breaks=29,main="P-value Histogram - Population as a factor",xlab="P-Value")

remove(top.pop.char)

gsub("HL60D0","Day 0 or 1",n_meta_data$Population)
levels(n_meta_data$Population)
cluster<-character()
for(x in n_meta_data$Population){
    if(x == "HL60D0"){
        cluster<-append(cluster,"Day 0 or 1")
    }
    if(x == "HL60D1"){
        cluster<-append(cluster,"Day 0 or 1")
    }
    if(x == "HL60D3"){
        cluster<-append(cluster,"Day 3 or 7")
    }
    if(x == "HL60D7"){
        cluster<-append(cluster,"Day 3 or 7")
    }
}
n_meta_data<-cbind(n_meta_data,cluster) #make factor for cluster

des.cluster<-model.matrix(~cluster,n_meta_data)
v.cluster<-voom(dge,des.cluster, plot=FALSE)
fit.cluster <- lmFit(v.cluster,des.cluster)
efit.cluster <- eBayes(fit.cluster)
top.cluster<-topTable(efit.cluster,number=Inf)

fdr.cluster<-subset(top.cluster,adj.P.Val<0.01)
top.cluster[1:10,]

remove(v.cluster)
remove(efit.cluster)

hist(top.cluster$adj.P.Val,breaks=19,main="Q-value Histogram - Cluster as a factor",xlab="Q-Value")
hist(top.cluster$P.Value,breaks=19,main="P-value Histogram - Cluster as a factor",xlab="P-Value")

# Define plotting function
cluster1<-gsub("HL60D0","Day 0 or 1",n_meta_data$Population)
cluster1<-gsub("HL60D1","Day 0 or 1",cluster1)
cluster1<-gsub("HL60D3","Day 3 or 7",cluster1)
cluster1<-gsub("HL60D7","Day 3 or 7",cluster1)
meta2<-cbind(n_meta_data,cluster1)

#function which combines expression data with metaData
prepareData2 <- function(genes) {
    miniDat <-  subset(dataNZ, rownames(dataNZ) %in% genes)
    miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
    gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
    levels = genes))
    miniDat <- suppressWarnings(data.frame(meta2, miniDat))
    
    return(miniDat)
}

topGenes <- rownames(top.cluster[1:10,])
prepDatTop<-prepareData2(topGenes) #format data

ggplot(data=subset(prepDatTop,prepDatTop$gene==topGenes[1]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(topGenes[1]) +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))

remove(top.cluster)

fit <- lmFit(data_voomed, design)
fit <- contrasts.fit(fit, contrasts=contrastMatrix)
efit <- eBayes(fit)

dt <- decideTests(efit)
summary(dt)

nInversD1D3D7 <- nrow(dt[(dt[,4] * dt[,6]== -1), ])
nSameD1D3D7 <- nrow(dt[(dt[,4] * dt[,6]== 1), ])
nTotalD1D3D7 <- nrow(dt[dt[,4] != 0 | dt[,6] != 0,c(4,6)])
vennDiagram(dt[dt[,4] != 0 | dt[,6] != 0,c(4,6)], circle.col=c("turquoise", "salmon", "green","yellow"), main = "Venn diagramm differentially expressed D1-D3 and D3-7")

# Finding genes that are differntailly expressed in the opposite direction between d1-3 and d3-7
commonDEgenesInverseDay1337 <- which(dt[,4]*dt[,6] == -1)
names(commonDEgenesInverseDay1337)

# Finding genes that are differntailly expressed in the opposite direction between d1-3 and d3-7
commonDEgenesInverseDay1337not01 <- which(dt[,4]*dt[,6] == -1 & dt[,1] == 0)
names(commonDEgenesInverseDay1337not01)

randRNA <- unique(names(commonDEgenesInverseDay1337not01))[sample(1:length(unique(names(commonDEgenesInverseDay1337not01))), 4, replace=FALSE)]
prepDatTop<-prepareData2(randRNA) #format data

f1 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[1]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[1])

f2 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[2]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[2])#106b-5p

f3 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[3]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[3])

f4 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[4]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[4])

grid.arrange(f1,f2,f3,f4,nrow=2) #put plots together

# Finding genes that are not zero in at least one comparaison between the groups
commonDEgenesDay0137 <- which(dt[,2]!=0 & dt[,3]!=0 & dt[,4]!=0 & dt[,5]!=0)
names(commonDEgenesDay0137)

vennDiagram(dt[,2:5], circle.col=c("turquoise", "salmon", "green","yellow"), main = "Venn diagramm showing differences between clusters")

#List of genes different between the two clusters
topGenes <- names(commonDEgenesDay0137)
#topGenes<-c("hsa-let-7d-3p","hsa-let-7i-5p","hsa-miR-125a-5p","hsa-miR-132-3p","hsa-miR-143-3p","hsa-miR-145-5p","hsa-miR-181b-3p","hsa-miR-185-5p","hsa-miR-192-5p","hsa-miR-19b-3p","hsa-miR-223-3p","hsa-miR-23a-3p","hsa-miR-26a-5p","hsa-miR-30e-5p","hsa-miR-425-5p","hsa-miR-425-3p","hsa-miR-92a-3p")
prepDatTop<-prepareData2(topGenes) #format data

#make dataset for probe 125a - remove cell with top expression
p125<-subset(prepDatTop,prepDatTop$gene==topGenes[6])
p125<-p125[-which(p125$gExp==max(p125$gExp)),]

f1 <- ggplot(data=p125,aes(cluster1,gExp,colour=Population))+
geom_jitter()+
ggtitle(topGenes[6]) + ylim(-100,2000)

f2 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==topGenes[1]),aes(cluster1,gExp,colour=Population))+
geom_jitter()+
ggtitle(topGenes[16]) + ylim(-100,2000)

f3 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==topGenes[2]),aes(cluster1,gExp,colour=Population))+
geom_jitter()+
ggtitle(topGenes[5])+ ylim(-100,2000)#106b-5p

f4 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==topGenes[3]),aes(cluster1,gExp,colour=Population))+
geom_jitter()+
ggtitle(topGenes[1]) + ylim(-100,2000)

grid.arrange(f1,f2,f3,f4,nrow=2) #put plots together

# Finding genes that are differentially expressed between clusters but not whithin clusters
commonDEgenesBetweenNotWhithin <- which(dt[,2]!=0 & dt[,3]!=0 & dt[,4]!=0 & dt[,5]!=0 & dt[,1]==0 & dt[,6]==0)
names(commonDEgenesBetweenNotWhithin)

ggplot(dataf4[dataf4$mirna %in% names(commonDEgenesBetweenNotWhithin),], aes(x=exp, fill = group)) +
geom_density(alpha=0.5) +
ggtitle("Random point plot for interaction miRNA-group") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5)) +
theme(legend.position="none") + xlim(0,30) +
facet_wrap(~mirna,ncol = 3)

ggplot(dataf4[dataf4$mirna %in% names(commonDEgenesBetweenNotWhithin),], aes(x=exp2, fill = group)) +
geom_density(alpha=0.5) +
ggtitle("Random density plot for interaction miRNA-group: no log") +
theme_bw() +
theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5)) +
theme(legend.position="none") + xlim(-500,10000) +
facet_wrap(~mirna,ncol = 3)

randRNA <- unique(names(commonDEgenesBetweenNotWhithin))[sample(1:length(unique(names(commonDEgenesBetweenNotWhithin))), 4, replace=FALSE)]
prepDatTop<-prepareData2(randRNA) #format data

f1 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[1]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[1])

f2 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[2]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[2])#106b-5p

f3 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[3]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[3])

f4 <- ggplot(data=subset(prepDatTop,prepDatTop$gene==randRNA[4]),aes(cluster,gExp,colour=Population))+
geom_jitter()+
ggtitle(randRNA[4])

grid.arrange(f1,f2,f3,f4,nrow=2) #put plots together
