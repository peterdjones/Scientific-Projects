#Applying k-means clustering to RNA-seq data on gene expression of mouse endothelial cells
#in laminar vs oscillatory shear stress, with or with Cx40 treatment
#GSE118717 from GEO (https://www.ncbi.nlm.nih.gov/geo/gse118717).

#load required packages
library(ggplot2)
library(DESeq2)

#extract data (column one is gene names)
countsData <- read.delim("GSE118717_Raw_counts.txt")
trimData <- countsData[,2:13]

#create column data array for DESeq analysis
samples <- colnames(trimData)
conditions <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
colData <- cbind(samples, conditions)

#use DESeq to extract and normalize data
dds <- DESeqDataSetFromMatrix(countData = trimData,
                              colData = colData,
                              design= ~ conditions)
                          
dds <- DESeq(dds)
vst <- getVarianceStabilizedData(dds)
df <- vst
df <- df - rowMeans(df)
df <- as.data.frame(df)


#carry out kmeans, and place each gene into a cluster
res_km <- kmeans(df, 5, nstart = 10)
data_plot <- data.table(melt(data.table(class = as.factor(res_km$cluster), df)))
data_plot[, Time := rep(1:ncol(df), each = nrow(df))]
data_plot[, ID := rep(1:nrow(df), ncol(df))]
head(data_plot)

# prepare centroids to plot 
centers <- data.table(melt(res_km$centers))
setnames(centers, c("Var1", "Var2"), c("class", "Time"))
centers[, ID := class]
centers[, gr := as.numeric(as.factor(Time))]
head(centers)
head(data_plot)

# plot the results
ggplot(data_plot, aes(variable, value, group = ID)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.65) +
  
  geom_line(data = centers, aes(gr, value),
            color = "green", alpha = 0.80, size = 1.2) +
  xlab("Normalised Expression")+
  ylab("Sample")+
  theme(axis.text.x=element_text(angle=90))+
  theme_bw()

