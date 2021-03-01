get.Agr.Metrics <- function(x, y){

  ### Calculate metrics of agreement as described in the paper:
  ### 
  ###   Duveiller, Fasbender & Meroni (2016) Revisiting the concept of 
  ###   a symmetric index of agreement for continuous datasets. 
  ###   Scientific Reports volume 6, Article number: 19401 
  ###   https://doi.org/10.1038/srep19401
  ###
  ### Please cite accordingly
  ### G. Duveiller 
    
  # size of the vectors (errors should be returned if not of the same size)
  n <- length(x)
  if(length(x)<=1){
    return(metrics = data.frame('L' = NA, 
                                'L.unsys' = NA, 
                                'f.sys' = NA,
                                'r' = NA, 
                                'b' = NA,
                                'a' = NA))}
  
  # Pearson moment-product coefficient of correlation
  r <- cor(x, y)
  
  # means of the vectors
  Xm <- mean(x); Ym <- mean(y)
  
  # covariance matrix
  S <- cov(cbind(x,y))
  
  # eigen decomposition
  E <- eigen(S)
  h <- cbind(x-Xm,y-Ym)%*%(E$vectors[,2])
  
  # kappa value to adjust when correlation is negative
  K <- sign(r<0) * 2*abs(sum((x - Xm)*(y - Ym)))
  
  # sum of squared differences (SSD)
  SSD <- sum((x - y)^2)
  # sum of squared (unsystematic) differences
  SSDu <- sum(2*(h)^2) 
  # sum of squared (systematic) differences (as diff of both)
  SSDs <- SSD - SSDu
  
  # the denominator, which maximizes the differences
  D <- sum((x - Xm)^2) + sum((y - Ym)^2) + n*(Xm - Ym)^2 + K
  
  # calc Lambda, the symmetric index of agreement
  L <- 1 - SSD/D
  
  # calc the unsystematic component of the index of agreement
  L.unsys <- 1 - SSDu/D
  
  # calc the fraction of systematic differences that are systematic
  f.sys <- SSDs/SSD
  
  b = (E$values[1] - S[1,1]) / S[1,2]
  a = Ym - b*Xm
  
  return(metrics = data.frame('L' = L, 
                              'L.unsys' = L.unsys,
                              'f.sys' = f.sys, 
                              'r' = r, 
                              'b' = b,
                              'a' = a))
}
