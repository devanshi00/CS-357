library("Rssa")

all = function(exp, nW){
  length = length(exp)
  stepFloat <- length / as.double(nW)
  step <- as.integer(ceiling(stepFloat))
  paa = vector(mode='numeric', length = nW)
  section_start <- 1
  j <- 1
  while(section_start <= length-step){
    section <- exp[section_start:(section_start+step-1)]
    paa[j] <- mean(section)
    section_start <- as.integer(j*stepFloat)
    j <- j+1
  }
  exp_lag = vector(mode='numeric', length = nW)
  norm <- sqrt(sum(paa^2))
  for(j in 1:57){
    exp_lag[j] <- paa[j] / norm
  }
  
  L = 1 + length(exp_lag) - 8 #8 is for number of elements in each row
  exp_ssa = ssa(exp_lag, L=L)  
  exp_ssa_recon = reconstruct(exp_ssa , groups=list(T=(1:4)))
  exp_after_smoothing = exp_ssa_recon$T
  window = 8
  slide = window/2
  rows = ceiling((length(exp_lag) - window+1)/(slide))
  #rows = as.integer(length(exp_lag)/window)
  exp_after_smoothing_mat = matrix(, nrow=rows, ncol=8) 
  start_sec = 0
  for(i in 1:rows){
    for(j in 1:window){
      exp_after_smoothing_mat[i,j] = exp_after_smoothing[start_sec+j]
    }
    start_sec = (start_sec + slide)
  }
  #construct ranges
  
  min = min(exp_after_smoothing_mat) - 0.01
  max = max(exp_after_smoothing_mat) + 0.01
  ranges <- list()
  step = as.double((max-min)/4)
  for(i in 1:4){
    ranges[[i]] = c(as.double(min), as.double(min+step))
    min = as.double(min + step)
  }
  
  #transform each element to its range
  
  exp_to_ranges = matrix(, ncol=8, nrow=rows)
  for(i in 1:rows){
    for(j in 1:8){
      for(k in 1:4){
        if(exp_after_smoothing_mat[i,j] >= ranges[[k]][1] & exp_after_smoothing_mat[i,j] < ranges[[k]][2]){
          exp_to_ranges[i,j] = k
        }
      }
    }
  }
  
  # number of each range happen in a matrix
  
  ranges_count = vector(mode="numeric", length=4)
  for(k in 1:4){
    for(i in 1:rows){
      for(j in 1:8){
        if(exp_to_ranges[i,j] == k){
          ranges_count[k] = ranges_count[k] + 1
        }
      }
    }
  }
  
  # probability of each element
  
  p = matrix(, ncol=8, nrow=rows)
  n = 8 * rows  # number of elements in matrix
  for(i in 1:rows){
    for(j in 1:8){
      for(k in 1:4){
        if(exp_to_ranges[i,j] == k){
          p[i,j] = ranges_count[k] / n
        }
      }
    }
  }
  
  # number of 2 ranges come together in rows
  ranges_joint_count = matrix(,ncol=4, nrow=4)
  for(k1 in 1:4){
    for(k2 in 1:4){
      r = 0
      for(i in 1:rows){
        see_k1 = FALSE
        see_k2 = FALSE
        for(j in 1:8){
          if(exp_to_ranges[i,j] == k1){
            see_k1 = TRUE
          }
          if(exp_to_ranges[i, j] == k2){
            see_k2 = TRUE
          }
        }
        if(see_k1==TRUE & see_k2==TRUE){
          r = r + 1
        }
      }
      ranges_joint_count[k1,k2] = r
    }
  }
  
  # joint probability(number of times two ranges come together in a row)
  # mutual entropy
  mutual_entropy = matrix(, ncol=rows, nrow=rows)
  gaus_x_y = matrix(, nrow=8, ncol=8)
  for(i in 1:rows){
    x = exp_to_ranges[i,]
    for(j in 1:rows){
      if(i != j){
        y = exp_to_ranges[j,]
        for(x_i in 1:8){
          for(y_j in 1:8){
            for(k1 in 1:4){
              for(k2 in 1:4){
                if(x[x_i] == k1 & y[y_j]==k2){
                  gaus_x_y[x_i,y_j] = (ranges_joint_count[k1,k2] / rows) 
                }
              }
            }
          }
        }
        mu = 0
        for(x_i in 1:8){
          for(y_j in 1:8){
            if(gaus_x_y[x_i,y_j] != 0){
              mu = mu + ( gaus_x_y[x_i,y_j] * log2( (gaus_x_y[x_i,y_j])/(p[i,x_i] * p[j,y_j]) ) )
            }
          }
        }
        mutual_entropy[i,j] = mu
      }
      else{
        mutual_entropy[i,j] = 0
      }
    }
  }
  
  return (list("lag" = exp_lag, "smooth_ver" = exp_after_smoothing, "smooth_matrix" = exp_after_smoothing_mat, "mutual_entropy" = mutual_entropy))
}

comparing = function(exp, mu){
  exp_cu = cusum(exp)
  v = vector(mode='numeric', length=nrow(mu)-1)
  for(i in 1:(nrow(mu)-1)){
    v[i] = mu[i,i+1]
  }
  return (list("cusum" = exp_cu, "mu_vec" = v, "max" = max(mu), "min" = min(mu[mu>0])))
}

find_fluc = function(comp){
  r = vector(mode='numeric', length=length(comp$mu_vec)-1)
  for(i in 1:length(r)){
    r[i] = abs(comp$mu_vec[i+1] - comp$mu_vec[i])
  }
  max1 = -1000000000
  max2 = -1000000000
  for(i in 1:length(r)){
    if(r[i] > max1){
      max1 = r[i]
      max1_index = i
    }
  } 
  for(i in 1:length(r)){
    if(r[i] > max2 & max1 != r[i]){
      max2 = r[i]
      max2_index = i
    }
  }
  return (list("range" = r, "max1_index" = max1_index, "max2_index" = max2_index ))
}

find_change_points = function(comp){
  it_is <- vector(mode='logical', length=length(comp$mu_vec) )
  it_is[1] = FALSE
  it_is[12] = FALSE
  for(i in 2:(length(comp$mu_vec)-1) ){
    if( ((comp$mu_vec[i+1] - comp$mu_vec[i]) < 0) & ((comp$mu_vec[i-1] - comp$mu_vec[i]) < 0) ){
      it_is[i] = TRUE
    }else if( ((comp$mu_vec[i+1] - comp$mu_vec[i]) > 0) & ((comp$mu_vec[i-1] - comp$mu_vec[i]) > 0) ){
      it_is[i] = TRUE
    }else{
      it_is[i] = FALSE
    }
  }
  return (it_is)
}

is_4_change_point = function(comp){
  it_is = FALSE
  if( ((comp$mu_vec[5] - comp$mu_vec[4]) < 0) & ((comp$mu_vec[3] - comp$mu_vec[4]) < 0) ){
    it_is = TRUE
  }else if( ((comp$mu_vec[5] - comp$mu_vec[4]) > 0) & ((comp$mu_vec[3] - comp$mu_vec[4]) > 0) ){
    it_is = TRUE
  }
  return (it_is)
}

is_7_change_point = function(comp){
  it_is = FALSE
  if( ((comp$mu_vec[8] - comp$mu_vec[7]) < 0) & ((comp$mu_vec[6] - comp$mu_vec[7]) < 0) ){
    it_is = TRUE
  }else if( ((comp$mu_vec[8] - comp$mu_vec[7]) > 0) & ((comp$mu_vec[6] - comp$mu_vec[7]) > 0) ){
    it_is = TRUE
  }
  return (it_is)
}

is_8_change_point = function(comp){
  it_is = FALSE
  if( ((comp$mu_vec[9] - comp$mu_vec[8]) < 0) & ((comp$mu_vec[7] - comp$mu_vec[8]) < 0) ){
    it_is = TRUE
  }else if( ((comp$mu_vec[9] - comp$mu_vec[8]) > 0) & ((comp$mu_vec[7] - comp$mu_vec[8]) > 0) ){
    it_is = TRUE
  }
  return (it_is)
}

is_6_change_point = function(comp){
  it_is = FALSE
  if( ((comp$mu_vec[7] - comp$mu_vec[6]) < 0) & ((comp$mu_vec[5] - comp$mu_vec[6]) < 0) ){
    it_is = TRUE
  }else if( ((comp$mu_vec[7] - comp$mu_vec[6]) > 0) & ((comp$mu_vec[5] - comp$mu_vec[6]) > 0) ){
    it_is = TRUE
  }
  return (it_is)
}

find_dif = function(fluc){
  max = -10000
  for(i in 1:length(fluc$range)){
    if(fluc$range[i] > max){
      max = fluc$range[i]
    }
  }
  min = 10000
  for(i in 1:length(fluc$range)){
    if(fluc$range[i] < min){
      min = fluc$range[i]
    }
  }
  return (list("max" = max, "min" = min, "dif" = ((max-min)/2) ))
}

find_mutual_dif = function(comp){
  max = -10000
  for(i in 1:length(comp$mu_vec)){
    if(comp$mu_vec[i] > max){
      max = comp$mu_vec[i]
    }
  }
  min = 10000
  for(i in 1:length(comp$mu_vec)){
    if(comp$mu_vec[i] < min){
      min = comp$mu_vec[i]
    }
  }
  return (list("max" = max, "min" = min, "dif" = ((max-min)/2) ))
}

if_dif_correct = function(fluc, dif){
  true_ones = NULL
  j = 1
  for(i in 1:length(fluc$range)){
    if(fluc$range[i] >= dif){
      true_ones[j] = i
      j = j + 1
    }
  }
  return (true_ones)
}

if_mutual_dif_correct = function(comp, dif){
  true_ones = NULL
  j = 1
  for( i in 1:(length(comp$mu_vec)-1) ){
    if( (abs(comp$mu_vec[i+1] - comp$mu_vec[i]) ) >= dif){
      true_ones[j] = i
      j = j + 1
    }
  }
  return (true_ones)
}

how_many_change_points = function(change_points){
  how_many_change_points = vector(mode='numeric', length=1000)
  for(i in 1:1000){
    for(j in 1:12){
      if(change_points[[i]][j] == TRUE){
        how_many_change_points[i] = how_many_change_points[i] + 1
      }
    }
  }
  return (how_many_change_points)
}

how_many_change_points_cusum = function(data_comp){
  how_many_change_points_cusum = vector(mode='numeric', length=1000)
  for(i in 1:1000){
    how_many_change_points_cusum[i] = length(data_comp[[i]]$cusum$violations$lower)
    if(length(data_comp[[i]]$cusum$violations$upper) != 0){
      for(j in 1:length(data_comp[[i]]$cusum$violations$upper)){
        same = FALSE
        if(length(data_comp[[i]]$cusum$violations$lower) != 0){
          for(k in 1:length(data_comp[[i]]$cusum$violations$lower)){
            if(data_comp[[i]]$cusum$violations$upper[j] == data_comp[[i]]$cusum$violations$lower[k]){
              same = TRUE
            }
          }
          if(same == FALSE){
            how_many_change_points_cusum[i] = how_many_change_points_cusum[i] + 1
          }
        }else{
            how_many_change_points_cusum[i] = how_many_change_points_cusum[i] + length(data_comp[[i]]$cusum$violations$upper)
        }
      }
    }
  }
  return(how_many_change_points_cusum)
}
