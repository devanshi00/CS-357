library("lsa")  # library to compute cosinus distance
library("party")
library("dtwclust")
library("cluster")
library("clValid")
library("caret")

# Lagrangian code and clustering
f.Lagrangian = function (data, nw){
  rows <- nrow(data)
  cols <- ncol(data)
  number_of_windows <- nw
  stepFloat <- cols / as.double(number_of_windows)
  step <- as.integer(ceiling(stepFloat))
  paa <- matrix(, nrow=rows, ncol=number_of_windows)
  for(i in 1:rows){
    section_start <- 1
    j <- 1
    sec <- data[i,]
    while(section_start <= cols-step){
      section <- subset(sec,select=section_start:(section_start+step-1))
      paa[i,j] <- rowMeans(section)
      section_start <- as.integer(j*stepFloat)
      j <- j+1
    }
  }
  
  lag <- matrix(, nrow=rows, ncol=number_of_windows)
  for(i in 1:rows){
    norm <- sqrt(sum(paa[i,]^2))
    for(j in 1:number_of_windows){
      lag[i,j] <- paa[i,j] / norm
    }
  }
  return(lag)
}

f.Hclust = function(lag, rows=nrow(lag), nC){
  matrixPrepared = matrix(, nrow = rows, ncol = rows)
  for(i in 1:rows){
    for(j in 1:rows){
      matrixPrepared[i,j] = cosine(lag[i,],lag[j,])
    }
  }
  dist2 = as.dist(1-matrixPrepared)
  number_of_cluster <- nC
  df1.hclust<- hclust(dist2, method="complete")
  plot(df1.hclust,hang=-1)
  df1.hclustk <- cutree(df1.hclust,k=number_of_cluster)
  sil <- silhouette(df1.hclustk, dist2)
  sil.summ <- summary(sil)
  return(list(df1.hclustk,table(df1.hclustk), sil.summ$si.summary[4], sil[,3]))
}

f.kmeans = function(lag, rows=nrow(lag), nC){
  matrixPrepared = matrix(, nrow = rows, ncol = rows)
  for(i in 1:rows){
    for(j in 1:rows){
      matrixPrepared[i,j] = cosine(lag[i,],lag[j,])
    }
  }
  dist2 = as.dist(1-matrixPrepared)
  number_of_cluster <- nC
  k= kmeans(dist2, centers=number_of_cluster)
  sil <- silhouette(k$cl, dist2)
  sil.summ <- summary(sil)
  return(list(k$cl,table(k$cl), sil.summ$si.summary[4], sil[,3]))
}

f.cosineAvg = function(lag, rows=nrow(lag), nw){
  avg <- vector(mode="numeric", length=nw)
  for(j in 1:nw){
    sum <- 0
    for(i in 1:rows){
      sum <- sum + lag[i,j]
    }
    avg[j] <- sum / rows
  }
  
  cosine_avg <- vector(mode="numeric", length=rows)
  for(i in 1:rows){
    cosine_avg[i] = cosine(lag[i,],avg)
  }
  return(cosine_avg)
}

f.Hclust_Euc = function(cosine_avg, nC){
  df2.hclust<- hclust(dist(cosine_avg,"euclidean"), method="complete")
  plot(df2.hclust,hang=-1)
  df2.hclustk <- cutree(df2.hclust,k=nC)
  sil <- silhouette(df2.hclustk, dist(cosine_avg, "euclidean"))
  sil.summ <- summary(sil)
  return(list(df2.hclustk,table(df2.hclustk), sil.summ$si.summary[4], sil[,3]))
}

f.kmeans_Euc = function(cosine_avg, nC){
  k= kmeans(cosine_avg, centers=nC)
  sil <- silhouette(k$cl, dist(cosine_avg,"euclidean"))
  sil.summ <- summary(sil)
  return(list(k$cl,table(k$cl), sil.summ$si.summary[4], sil[,3]))
}

# for BEATS and comparison
f.sliwi <- function(slide=8, window=64,size_data = n){
  iM <- list()
  nM <- as.integer((size_data-window)/slide + 1) # number of matrices
  for(i in 0:(nM-1)){
    index = (1+i*slide):(window+i*slide)
    iM[[i+1]] <- index
  }
  return(iM)
}

f.DCT <- function(A,U){
  B = U%*%A%*%t(U)
  return (B)
}

f.U <- function(n){
  res =NULL
  for (k in 0:(n - 1)) {
    res[[k + 1]] <- c(cos(pi/n * ((0:(n - 1)) + 0.5) * k))
  }
  res[[1]] =   res[[1]]*1/(sqrt(2))
  #res
  matrix(unlist(res), n,n, byrow = T)
}

f.roundingfactor = function(value){
  rf = 2
  x = abs(floor(log10(value))-4)
  if(x > 2){
    rf = x
  }
  return(rf)
}

q = c(16,12,14,14,18,24,49,72,11,12,13,17,22,35,64,92,10,14,16,22,37,55,78,95,16,19,24,29,56,64,
      87,98,24,26,40,51,68,81,103,112,40,58,57,87,109,104,121,100,51,60,69,80,103,113,120,103,61,55,
      56,62,77,92,101,99)
Z = matrix(q, 8, 8)

f.BEATS = function(sc,sl = 64,wd = 64, size =(ncol(sc)-1)){
  beatsData = NULL
  indices = f.sliwi(slide = sl, window = wd, size_data = size)
  for(i in 1:nrow(sc)){
    v = NULL
    a <- t(sc[i,-1])
    for(k in 1:length(indices)){
      matriz = matrix(a[indices[[k]]],sqrt(wd))
      v1 <- f.DCT(matriz,f.U(sqrt(wd)))
      x = Mod(eigen(round(v1/Z, f.roundingfactor(max(v1)))[1:4, 1:4])$values)
      sentence = !duplicated(x)
      if(sum(sentence)>= 3){
        x = x[sentence][1:3]
      }
      if(sum(sentence)<3){
        sentence[which(!sentence)[1]] = T
        x = x[sentence]
      }
      v = c(v, x)
    }
    beatsData <- rbind(beatsData, v)
  }
  rownames(beatsData) = 1:nrow(sc)
  return(beatsData)
}

# for SAX_sd and comparison

PAA_sd = function (x, w)
{
  if ((w - floor(w)) > 0) {
    stop("w (number of frames) must be an integer")
  }
  n <- length(x)
  if (w > n) {
    stop("cannot have more parts than the length of the series")
  }
  PAA <- rep(0, w)
  d = n/w
  breakpoints <- seq(0, n, d)
  for (i in 1:w) {
    init <- breakpoints[i] + 1
    end <- breakpoints[i + 1]
    frac_first <- ceiling(init) - init
    frac_end <- end - floor(end)
    interv = floor(init):ceiling(end)
    sec <- x[interv]
    if (frac_first > 0) {
      sec[1] = sec[1] * frac_first
    }
    if (frac_end > 0) {
      sec[length(sec)] = sec[length(sec)] * frac_end
    }
    PAA[i] = sd(sec)
  }
  PAA
}


f.toSAX_sd <- function(sc, w, alpha){
  saxSdData <- NULL
  for (i in 1:nrow(sc)) {
    a <- t(sc[i,-1])
    #x <- (a - mean(a)) /sd(a)
    #paax <- PAA(x, w) #generate PAA reductions
    paasd <- PAA_sd(a, w) #generate sd
    #SAXx <- convert.to.SAX.symbol( paax, alpha)
    saxSdData = rbind(saxSdData, paasd)
  }
  rownames(saxSdData) = 1:nrow(sc)
  return(as.data.frame(saxSdData))
}
