perma_predict_give_prob = function(
  # take in transformed train and test.
  # train and test matrix have row sample column gene.
  train_X,
  test_X,
  train_y,
  # test_y,
  perma_alpha,
  perma_tau
  # default K is K2.
) {

  X = as.matrix(train_X)
  # message(paste('dim(X) is',dim(X)))
  Xt = as.matrix(test_X)
  # message(paste('dim(Xt) is', dim(Xt)))
  
  
  temp = as.character(train_y)
  temp[temp == sort(levels(train_y))[1]] = 1
  temp[temp == sort(levels(train_y))[2]] = 2
  
  xclass = as.numeric(temp)  # 
  # message(table(xclass))
  # tclass = as.numeric(test_y)
  # message(table(tclass))
  
  # npts = length(tclass)
  
  ####
  # npts = length(tclass)   # number of testing points.
  npts = dim(Xt)[1]       # number of testing points.
  # npts = length(tclass)   
  nobs = dim(X)[1]        # number of training points
  # nobs = length(xclass)               # number of training points   286
  # n <- nobs + npts       # 286 + 13 = 299
  n = npts + nobs
  
  d1 <- rep(0, n * n)                    # initial distance matrix   299*299
  dim(d1) <- c(n, n)
  #
  
  Xt_X = rbind(Xt, X)  # Xt test in front, X in behind.
  # dim(Xt_X)
  
 
  # distance matrix dist()
  d1 = stats::dist(Xt_X, method = 'euclidean', diag = T, upper = T)
  d1 = as.matrix(d1)
  # print(paste('dim(d1) is', dim(d1)))
  
  train_median = median( d1[(npts+1):n, (npts+1):n] )
  d1 = d1/train_median
  
  ################# K2
  # message('K2 is on')
  Kmat <- K2(d1, tau = perma_tau)
  # message('K2 is done tempK2 is on')
  ###################################    it takes 4 sec
  tempK2 <- predict.self(Kmat, npts, nobs, class = xclass, alpha = perma_alpha) # 4 sec.
  #############################################################################################
  # tempK2
  # message('tempK2 is done.')
 
  # # probability of true class
  # predProbTrueK2 = rep(0, tNum); dim(predProbTrueK2) = c(tNum, 1)
  # 
  # predProbTrueK2[1:tNum, 1] <- tempK2$R3[cbind(seq(1:tNum), tclass)]
  
  ### output of predict.self, tempK2$R3: p1,...,p_k; RatioSum; Class: 1st largest; Class2: 2nd largest -- based on R1, R2, R3

  # probability of class 1 based on K2.
  predProbClass1K2 = rep(0, npts); dim(predProbClass1K2) = c(npts, 1)
  
  predProbClass1K2[1:npts, 1] <- tempK2$R3[, 1]  # R3
  # predProbClass1K2[1:npts, 1] <- tempK2$R2[, 1]    # R2
  
  # predicted class
  predClass = rep(0, npts); dim(predClass) = c(npts, 1)
  predClass[1:npts, 1] = tempK2$R3[, 4] # predicted class from perma.
  
  # rownames(predProbClass1K2) = "test_sample"
  # colnames(predProbClass1K2) = 'class_1_prob'
  
  
  
  return(data.frame(predProbClass1K2 = predProbClass1K2,
                    predClassPerma = predClass)) 
}

# Gaussian covariance
K2 = function(x, tau) {
  exp(-x^2/tau^2);        # tau square.
}

### R function to compute the predictive distribution for the classes
### Other functions needed: per.ratio
### Input: Kmat, npts, nobs, alpha (see per.ratio),
###        xclass -- class labels, from 1 to k
### Output: p1,...,p_k; RatioSum; Class: 1st largest; Class2: 2nd largest -- based on R1, R2, R3
predict.self = function(Kmat, npts, nobs, class, alpha) {
  # message('predict.self is on')
  n = npts + nobs;
  
  k = max(class);   # number of classes.   # class should be 1, 2, 3, 4, ..., can not be 0, 1, ....
  # k = length(unique(class));
  if(missing(alpha)) alpha = 1;
  
  # message(paste('npts is', npts, 'n is', n, 'k is', k))
  
  npts_k_3 = npts * k * 3
  # message(paste('npts_k_3 is', npts_k_3))
  # ratio = rep(0,npts*k*3);
  ratio = rep(0, npts_k_3)
  dim(ratio)=c(npts,k,3);
  
  for(r in 1:k) {
    nr = sum(class==r);
    index = as.logical(c(rep(1, npts), as.numeric(class==r)));
    temp = per.ratio(Kmat[index, index], npts, nr, alpha);
    ratio[,r,1] = temp$R1;
    ratio[,r,2] = temp$R2;
    ratio[,r,3] = temp$R3;
  }
  
  ans = rep(0,npts*(k+3)*3);
  dim(ans) = c(npts,k+3,3);
  Class = Class2 = rep(1,npts);     # initialize class to be all 1 at beginning.
  
  for(r in 1:3) {
    RatioSum = ratio[,,r] %*% rep(1, k);
    for(i in 1:npts) {
      Class[i] = sort.list(ratio[i,,r])[k];     # for class 1 and class 2, Class is the one whose ratio is biggest.
      Class2[i] = sort.list(ratio[i,,r])[k-1];  #                          Class2 is the one whose ratio is second biggest.
    }
    ans[,,r] = cbind(ratio[,,r]/(RatioSum %*% t(rep(1, k))), RatioSum, Class, Class2);
  }
  
  # directly give class, not only probability of being that class.           Class is the one whose ratio is biggest.
  #                                                                                  Class2 is the one whose ratio is second biggest. 
  #                                                                          Class here is in fact taking probability cutoff 0.5
  #                                                                          that if probability of being class 2 > 0.5, that sample goes
  #                                                                          to class 2.
  R1 = ans[,,1]; if(npts==1) R1 = t(R1); colnames(R1)=c(seq(1:k),"RatioSum","Class","Class2"); 
  R2 = ans[,,2]; if(npts==1) R2 = t(R2); colnames(R2)=c(seq(1:k),"RatioSum","Class","Class2");
  R3 = ans[,,3]; if(npts==1) R3 = t(R3); colnames(R3)=c(seq(1:k),"RatioSum","Class","Class2");
  
  list(R1=R1,R2=R2,R3=R3);
}

### R function to compute approximated permanent ratios based on R1,R2,R3
### Input: Kmat=K(c(tvec,x),c(tvec,x)) -- Matix with dim (npts+nobs)*(npts+nobs);
###        npts -- number of testing points;  nobs -- number of training points;
###        alpha -- parameter in per_alpha[K](X), (default value is 1 which leading to usual permanent)
### Output: R1, R2, R3 -- vectors per_alpha[K](X,t)/per_alpha[K](X) based on R1,R2,R3, dim: npts
per.ratio=function(Kmat, npts, nobs, alpha) {
  n = nobs + npts;
  if(missing(alpha)) alpha = 1;  # default alpha is 1.
  wt = diag(Kmat);
  Kmat = diag(1/sqrt(wt)) %*% Kmat %*% diag(1/sqrt(wt));
  Kxx = Kmat[(npts+1):n, (npts+1):n];
  Ktx = Kmat[1:npts, (npts+1):n];
  if(npts==1) Ktx = t(as.matrix(Ktx));
  kxx2m = Kxx^2;

  R1 = apply(Ktx^2, 1, sum) + alpha;           ## R1(t;X)
  R1 = R1 * wt[1:npts];
  
  temp = kxx2m; diag(temp) = alpha;
  R1i = apply(temp, 1, sum);                   ## R1(xi;X-i)
  
  temp = Kxx;  diag(temp) = alpha;
  R2 = ((Ktx %*% temp) * Ktx) %*% (1/R1i);     ## R2(t;X)
  R2 = as.vector(R2) + alpha;
  R2 = R2 * wt[1:npts];
  
  R1ij = matrix(0, nobs, nobs);                ## R1(xj;X-i-j)

  for(i in 1:nobs) {
    temp = kxx2m; temp[,i] = 0; diag(temp)=alpha;
    R1ij[i,] = apply(temp, 1, sum);
  }

  diag(R1ij) = R1i;

  R2i = rep(0, nobs);                          ## R2(xi;X-i)

  for(i in 1:nobs) {
    temp = Kxx; temp[,i] = 0; diag(temp)=alpha;
    temp = (temp %*% Kxx[i,]) * Kxx[i,] / R1ij[i,];
    temp[i] = alpha;
    R2i[i] = sum(temp);
  }
  
  temp = matrix(0, npts, nobs);
  
  for(i in 1:nobs) {
    temp1 = Kxx; temp1[i,] = 0; diag(temp1)=alpha;
    temp1 = Ktx %*% temp1;
    temp2 = Kxx[i,] / R1ij[i,]; temp2[i]= 0;
    temp[,i] = temp1 %*% temp2;
  }
  
  R3 = ((temp + alpha * Ktx) * Ktx) %*% (1/R2i);    ## R3(t;X)
  R3 = as.vector(R3) + alpha;
  R3 = R3 * wt[1:npts];
  
  list(R1=R1, R2=R2, R3=R3);
}

