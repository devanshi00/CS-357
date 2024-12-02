# example of experiments for Arrow head data

df <- read.table("../data/ArrowHead_TRAIN.csv", sep=",")
source("functions_lagrangian.R")
names(df)[1] = "classId"
classId <- as.factor(df$classId)
arrhBEATS = f.BEATS(sc = df)
arrhLag = f.Lagrangian(data = df[,-1], nw=ncol(arrhBEATS))
cosine_avg = f.cosineAvg(arrhLag, rows=nrow(arrhLag), nw=ncol(arrhBEATS))
seeds = c(1234, 4321, 1357, 2468, 5678, 9753, 8642, 6789, 2345, 3456)
arrhBEATS_clust = vector(mode="numeric", length=10)
arrhLag_clust = vector(mode="numeric", length=10)
arrhLag2_clust = vector(mode="numeric", length=10)
partition = 0.9
for(i in 1:10){
  set.seed(seeds[i])
  iS = createDataPartition(classId, p=partition)$Resample1
  sil = silhouette(cutree(hclust(dist(arrhBEATS[iS,])),length(levels(factor(classId)))),dist(arrhBEATS[iS,]))
  arrhBEATS_clust[i] = mean(sil[,3])
  arrhLag_clust[i] = f.Hclust(lag=arrhLag[iS,], nC=length(levels(factor(classId))))[[3]]
  arrhLag2_clust[i] = f.Hclust_Euc(cosine_avg=cosine_avg[iS], nC=length(levels(factor(classId))))[[3]]
}
mean(arrhBEATS_clust)
mean(arrhLag_clust)
mean(arrhLag2_clust)
save(arrhBEATS_clust, file = "../results/arrhBEATS_clust.Rdata")
save(arrhLag_clust, file = "../results/arrhLag_clust.Rdata")
save(arrhLag2_clust, file = "../results/arrhLag2_clust.Rdata")

arrhBEATS_kmeans = vector(mode="numeric", length=10)
arrhLag_kmeans = vector(mode="numeric", length=10)
arrhLag2_kmeans = vector(mode="numeric", length=10)
partition = 0.9
for(i in 1:10){
  set.seed(seeds[i])
  iS = createDataPartition(classId, p=partition)$Resample1
  arrhLag_kmeans[i] = f.kmeans(lag=arrhLag[iS,], nC=length(levels(factor(classId))))[[3]]
  arrhLag2_kmeans[i] = f.kmeans_Euc(cosine_avg=cosine_avg[iS], nC=length(levels(factor(classId))))[[3]]
}

mean(arrhLag_kmeans)
mean(arrhLag2_kmeans)
save(arrhLag_kmeans, file = "../results/arrhLag_kmeans.Rdata")
save(arrhLag2_kmeans, file = "../results/arrhLag2_kmeans.Rdata")

set.seed(1234)
partition = 0.1
iS = createDataPartition(classId, p=partition)$Resample1
w= ncol(arrhBEATS)
alpha = 10
arrhSAX = f.toSAX_sd(df, w, alpha)
arrhSAX_clust = vector(mode="numeric", length=10)
partition = 0.9
for(i in 1:10){
  set.seed(seeds[i])
  iS = createDataPartition(classId, p=partition)$Resample1
  sil = silhouette(cutree(hclust(dist(arrhSAX[iS,])),length(levels(factor(classId)))),dist(arrhSAX[iS,]))
  arrhSAX_clust[i] = mean(sil[,3])
}
mean(arrhSAX_clust)
save(arrhSAX_clust, file = "../results/arrhSAX_clust.Rdata")
