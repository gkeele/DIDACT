# Creates a group-identity matrix given a vector of factors
incidence.matrix <- function(fact){
  m=diag(nlevels(fact))[fact,]
  colnames(m)=levels(fact)
  return(m)
}
# Creates constraint matrix
make.M <- function(X, n){
  j <- ncol(X)
  k <- (-1 + sqrt(j))*(j - 1)^(-3/2)
  m <- 1/sqrt(j - 1)
  c <- (n - 2)*k + m
  M <- diag(j - 1)
  M[M == 0] <- -k
  M <- rbind(M, rep(-m, j - 1))
  return(M)
}

# Update functions
update.beta <- function(X, y, res.var, prior.var){
  post.var <- chol2inv(chol(t(X) %*% X + res.var*chol2inv(chol(prior.var))))*res.var
  post.mean <- (1/res.var)*post.var %*% t(X) %*% y
  beta.vec <- c(mnormt::rmnorm(1, mean=post.mean, varcov=post.var))
  return(beta.vec)
}
update.sigma.2 <- function(X, y, prior.alpha, prior.beta, par.vec)
{
  post.alpha <- (nrow(X) + prior.alpha)/2
  post.beta <- (t(y - X %*% par.vec) %*% (y - X %*% par.vec) + prior.beta*prior.alpha)/2
  sigma.2 <- 1/rgamma(1, shape=post.alpha, rate=post.beta)
  return(sigma.2)
}
update.tau <- function(K, prior.alpha, prior.beta, par.vec)
{
  post.alpha <- (ncol(K) + prior.alpha)/2
  post.beta <- (t(par.vec) %*% chol2inv(chol(K)) %*% par.vec + prior.beta)/2
  tau <- 1/rgamma(1, shape=post.alpha, rate=post.beta)
  return(tau)
}

# Need to run this before evaluating utility
# Gibbs sampler for diallel data
#' @export
diallel.gibbs <- function(phenotype, sex, is.female=TRUE, mother.str, father.str, n.iter, burn.in, multi.chain=1, thin=1,
                          sigma.2.starter=5, taua.starter=2, taud.starter=2, tauo.starter=2, taue.starter=2,
                          strain.reorder=c(8, 7, 4, 6, 5, 1, 3, 2), use.constraint=T)
{
  # Defining strain columns and incidence matrices
  mother.str <- factor(mother.str, levels(mother.str)[strain.reorder])
  father.str <- factor(father.str, levels(father.str)[strain.reorder])
  mom.mat <- incidence.matrix(mother.str)
  pop.mat <- incidence.matrix(father.str)
  
  # For labeling
  strain.names <- colnames(mom.mat)

  # Labels
  add.names <- paste(strain.names, "add")
  dom.names <- paste(strain.names, "inbred")
  mat.names <- paste(strain.names, "mat")
  epi.names <- paste(rep(NA, (length(strain.names)*(length(strain.names)-1))/2), "epi")
  counter <- 1
  for(m in 1:length(strain.names)){
    for(n in 1:length(strain.names)){
      if(m > n){
        epi.names[counter] <- paste(strain.names[m], strain.names[n], sep=" ")
        counter <- counter + 1
      }
    }
  }
  # Setting up the data
  y <- as.vector(phenotype)
  n <- length(y)
  
  # By default: female = 1, male = 0
  if(!is.female){
    sex <- ifelse(sex, 1, 0)
  } 
  
  # Make X
  X <- cbind(rep(1, n), sex, ifelse(mother.str == father.str, 1, 0))
  
  # Setting up Z
  # Additive
  add.part <- mom.mat + pop.mat
  # Inbred
  inbred.part <- 1*(add.part == 2)
  # Maternal
  mat.part <- mom.mat - pop.mat
  # Epistatic
  mom.index <- apply(mom.mat, 1, function(x) which(x == 1))
  pop.index <- apply(pop.mat, 1, function(x) which(x == 1))
  combined.index <- cbind(mom.index, pop.index)
  sorted.index <- t(apply(combined.index, 1, function(x) sort(x, decreasing=T)))
  make.one <- paste(strain.names[sorted.index[,1]], strain.names[sorted.index[,2]], sep=" ")
  epi.part <- 1*(t(apply(matrix(make.one, nrow=length(make.one), ncol=1), 1, function(x) x == epi.names)))
  
  # Making constraint matrix
  if(use.constraint){
    M.add <- make.M(X=add.part, n=n)
    M.inbred <- make.M(X=inbred.part, n=n)
    M.mat <- make.M(X=mat.part, n=n)
    M.epi <- make.M(X=epi.part, n=n)
    
    # Transforming design matrices
    add.part <- add.part %*% M.add
    inbred.part <- inbred.part %*% M.inbred
    mat.part <- mat.part %*% M.mat
    epi.part <- epi.part %*% M.epi
  }

  # Overall design matrix
  X.all <- cbind(X, add.part, inbred.part, mat.part, epi.part)
  # Progress bar
  pb <- txtProgressBar(min=0, max=n.iter*multi.chain)
  # For multiple chains
  if (multi.chain > 1) {
    chain.list <- list()
  } 
  for (j in 1:multi.chain) {
    # Allocate memory
    p.mat <- matrix(0, n.iter, 60) # 3 fixed effects, 8 additive, 8 inbred, 8 maternal, 28 epistatic
    
    # Set hyperparameters
    ga.alpha <- 0.002
    ga.beta <- 2
    
    # Initialization
    fix.vec <- c(mean(phenotype[sex==0]), mean(phenotype[sex==1]), 0)
    rand.vec <- rnorm(ncol(X.all) - 3, 0, 20)
    beta.vec <- c(fix.vec, rand.vec)
    
    sigma.2 <- sigma.2.starter
    taua <- taua.starter
    taud <- taud.starter
    tauo <- tauo.starter
    taue <- taue.starter
    
    # Keep track of matrix rows
    counter <- 1
    
    # Updating parameters
    for(i in 1:((n.iter*thin)+burn.in)){
      # prior covariance matrix
      Sigma0 <- diag(c(rep(1000, 3), rep(taua, ncol(add.part)), 
                       rep(taud, ncol(inbred.part)), rep(tauo, ncol(mat.part)), 
                       rep(taue, ncol(epi.part))))
      
      # Update beta.vec
      beta.vec <- update.beta(X.all, y, sigma.2, Sigma0)
      
      # Update sigma.2
      sigma.2 <- update.sigma.2(X.all, y, ga.alpha, ga.beta, beta.vec)
      
      # Updating taua
      a.index <- 4:(4 + ncol(add.part) - 1)
      taua <- update.tau(diag(ncol(add.part)), ga.alpha, ga.beta, beta.vec[a.index])
      
      # Updating taud
      i.index <- (4 + ncol(add.part)):(4 + ncol(add.part) + ncol(inbred.part) - 1)
      taud <- update.tau(diag(ncol(inbred.part)), ga.alpha, ga.beta, beta.vec[i.index])
      
      # Updating tauo
      m.index <- (4 + ncol(add.part) + ncol(inbred.part)):(4 + ncol(add.part) + ncol(inbred.part) + ncol(mat.part) - 1)
      tauo <- update.tau(diag(ncol(mat.part)), ga.alpha, ga.beta, beta.vec[m.index])
      
      # Updating taue
      e.index <- (4 + ncol(add.part) + ncol(inbred.part) + ncol(mat.part)):(4 + ncol(add.part) + ncol(inbred.part) + ncol(mat.part) + ncol(epi.part) - 1)
      taue <- update.tau(diag(ncol(epi.part)), ga.alpha, ga.beta, beta.vec[e.index])
      
      # Keep values after burn-in
      if(i > burn.in){
        if((i-1-burn.in) %% thin == 0){
          if(use.constraint){
            p.mat[counter,] <- c(beta.vec[1:3], 
                                 M.add %*% beta.vec[a.index],
                                 M.inbred %*% beta.vec[i.index],
                                 M.mat %*% beta.vec[m.index],
                                 M.epi %*% beta.vec[e.index],
                                 sigma.2, taua, taud, tauo, taue)
          }
          else{
            p.mat[counter,] <- c(beta.vec, sigma.2, taua, taud, tauo, taue)
          }
          counter <- counter + 1
          # Makes a progress bar
          setTxtProgressBar(pb, counter)
        }
      }
    }
    colnames(p.mat) <- c("mu", "female", "inbred", add.names, dom.names, mat.names, epi.names, "sigma2", "tau2 add", "tau2 inbred", "tau2 mat", "tau2 epi")
    
    # If multiple chains, build up list
    if(multi.chain > 1){
      chain.list[[j]] <- coda::as.mcmc(p.mat)
    } 
  }
  
  # Return matrix or list of matrix
  if(multi.chain > 1){
    return(coda::as.mcmc.list(chain.list))
  }
  else{
    return(coda::as.mcmc(p.mat))
  }
}
