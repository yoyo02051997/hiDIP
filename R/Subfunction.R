#Data information
datainf <- function(data, datatype = "abundance"){
  
  if( datatype == 'abundance'){
    
    if (class(data) != "data.frame") data <- as.data.frame(data)
    name <- c("n", "S.obs", "C.hat", paste0("f",1:10))
    sub <- function(dat){
      out <- matrix(0, length(name), 1, dimnames = list(1:length(name), "value"))
      rownames(out) <- name
      out[1,1] <- as.integer(sum(dat))
      out[2,1] <-  round(c(sum(dat!=0, na.rm = T)), 0)
      out[4:13,1] <- sapply(1:10, function(x) sum(dat == x))
      f1 = sum(dat == 1)
      f2 = sum(dat == 2)
      n = sum(dat)
      out[3,1] <- if(f1 == 0 & f2 == 0) 0 else ifelse(f2>0, round(1 - (f1/n)*((n-1)*f1/((n-1)*f1+2*f2)), 3), round(1 - (f1/n)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)), 3)) 
      if(out[3,1]>=1) out[3,1] <- 0.999
      
      out
    }
    
    outmat <- apply(data, MARGIN = 2, function(x) sub(x))
    outmat[-3,] <- as.integer(outmat[-3,])  
    outmat = t(outmat)
    outmat <- data.frame(site = rownames(outmat), outmat)
    colnames(outmat) <- c("Assemblage",name)
    rownames(outmat) <- NULL
    
    res <- as.data.frame(outmat)
    
    
    
  } else {
    
    if (class(data) != "data.frame") data <- as.data.frame(data)
    name <- c("T", "U", "S.obs", "C.hat", paste0("Q",1:10))
    sub <- function(dat){
      out <- matrix(0, length(name), 1, dimnames = list(1:length(name), "value"))
      rownames(out) <- name
      out[1,1] <- as.integer(dat[1])
      dat <- dat[-1]
      out[2,1] <- as.integer(sum(dat))
      out[3,1] <-  round(c(sum(dat!=0, na.rm = T)), 0)
      out[5:14,1] <- sapply(1:10, function(x) sum(dat == x))
      f1 = sum(dat == 1)
      f2 = sum(dat == 2)
      n = sum(dat)
      out[4,1] <- if(f1 == 0 & f2 == 0) 0 else ifelse(f2>0, round(1 - (f1/n)*((out[1,1]-1)*f1/((out[1,1]-1)*f1+2*f2)), 3), round(1 - (f1/n)*((out[1,1]-1)*(f1-1)/((out[1,1]-1)*(f1-1)+2)), 3))
      if(out[4,1]>=1) out[4,1] <- 0.999
      
      out
    }
    
    outmat <- apply(data, MARGIN = 2, function(x) sub(x))
    outmat[-4,] <- round(outmat[-4,]) 
    outmat = t(outmat)
    outmat <- data.frame(site = rownames(outmat), outmat)
    colnames(outmat) <- c("Assemblage",name)
    rownames(outmat) <- NULL
    
    res <- as.data.frame(outmat)
    
  }
  
  return(res)
  
}
convToNewick <- function(tree){
  tree<-reorder.phylo(tree,"cladewise")
  n<-length(tree$tip)
  string<-vector(); string[1]<-"("; j<-2
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,2]<=n){
      string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
      if(!is.null(tree$edge.length)){
        string[j]<-paste(c(":",round(tree$edge.length[i],10)), collapse="")
        j<-j+1
      }
      v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
      while(length(v)>0&&k==v[length(v)]){
        string[j]<-")"; j<-j+1
        w<-which(tree$edge[,2]==tree$edge[k,1])
        if(!is.null(tree$edge.length)){
          string[j]<-paste(c(":",round(tree$edge.length[w],10)), collapse="")
          j<-j+1
        }
        v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
      }
      string[j]<-","; j<-j+1
    } else if(tree$edge[i,2]>=n){
      string[j]<-"("; j<-j+1
    }
  }
  if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)], ";")
  else string<-c(string[1:(length(string)-2)],";")
  string<-paste(string,collapse="")
  return(string)
}
datainfphy <- function(data, tree){
  tree <- newick2phylog(convToNewick(tree))
  name <- c("n", "S.obs", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")
  sub <- function(dat){
    tmp <- choose_data(dat, tree)
    #rtreephy <- newick2phylog(convToNewick(tree))
    out <- matrix(0, length(name), 1, dimnames = list(1:length(name), "value"))
    out[1, 1] <- sum(dat)
    out[2, 1] <- sum(dat > 0)
    f1 = sum(tmp[, 1] == 1)
    f2 = sum(tmp[, 1] == 2)
    node <- names(tree$parts)
    xx <- tmp[,1]
    names(xx) <- rownames(tmp)
    node2 = names(xx)[xx==2][names(xx)[xx==2] %in% node]
    if(length(node2)>0){
      rep2 <- sapply(node2, function(tx){
        sum(xx[tree$parts[[tx]]] %in% c(2,0)) == length(tree$parts[[tx]]) 
      })
      f2 <- f2 - sum(rep2)
    }
    node1 = names(xx)[xx==1][names(xx)[xx==1] %in% node]
    if(length(node1)>0){
      rep1 <- sapply(node1, function(tx){
        sum(xx[tree$parts[[tx]]] %in% c(1,0)) == length(tree$parts[[tx]]) 
      })
      f1 <- f1 - sum(rep1)
    }
    out[3, 1] <- f1
    out[4, 1] <- f2
    out[5, 1] <- round(sum(tmp[tmp[, 1] == 1, 2]), 2)
    out[6, 1] <- round(sum(tmp[tmp[, 1] == 2, 2]), 2)
    out[7, 1] <- sum(tmp[tmp[, 1]>0, 2])
    tmp[, 1] <- tmp[, 1]/tmp[nrow(tmp), 1]
    out[8, 1] <- sum(tmp[, 1]*tmp[, 2])
    # out[11, 1] <- -sum(tmp[tmp[,1]>0,2]*tmp[tmp[,1]>0,1]*log(tmp[tmp[,1]>0,1]))
    # out[12, 1] <- out[10, 1]*exp(out[11, 1]/out[10, 1])
    return(round(out,3))
  }
  outmat <- apply(data, 2,function(x){
    sub(x)
  })
  outmat = t(outmat)
  outmat <- data.frame(site = rownames(outmat), outmat)
  colnames(outmat) <- c("Assemblage",name)
  rownames(outmat) <- NULL
  outmat
}
plot.tree2 <- function(mat){
  #number of lower level must be large than or equal to the number of higher level
  t <- apply(mat,MARGIN = 1, function(x) length(unique(x)))
  if(sum((t[-1]-t[-length(t)])<0)>0) stop("number of lower level must be large than or equal to the number of higher level, please renew your structure matrix.")
  rownames(mat) <- paste0("level", 1:nrow(mat))
  colnames(mat) <- paste0("community", 1:ncol(mat))
  mat <- data.frame(t(mat), stringsAsFactors = F)
  m <- ncol(mat)
  mat$pathString <- apply(mat,1,paste,collapse="/") 
  population <- as.Node(mat)
  useRtreeList <- ToListExplicit(population, unname = TRUE)
  
  diagonalNetwork(useRtreeList,fontSize = 27, opacity = 10, linkColour = "#828282", nodeStroke = "#6495ED")
}
#Taxonomy
entropy_true <- function(p, q) {
  p <- p[p>0]
  sub <- function(q){
    if(q!=1) (1-sum(p^q))/(q-1)
    else -sum(p*log(p))
  }
  sapply(q, sub)
}

entropy_true_inc <- function(p, q){
  x = p[-1]
  U = sum(x)
  x <- x[x != 0]
  ai <- x/U
  Sub <- function(q){
    if (q == 1) -sum(ai*log(ai)) else qD_MLE(q,ai)
  }
  sapply(q, Sub)
}


entropy_est <- function(x, q){
  x = x[x > 0]
  n = sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  A = ifelse(f2>0, 2*f2/((n-1)*f1+2*f2), ifelse(f1>0, 2/((n-1)*(f1-1)+2), 1))
  
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)-1
    }
    else if(q==1){
      r <- 1:(n-1)
      a <- sum(x/n*(digamma(n)-digamma(x)))
      b <- ifelse(f1 == 0|A == 1, 0, f1/n*(1-A)^(1-n)*(-log(A)-sum((1-A)^r/r)))
      a+b
    }else if(abs(q-round(q)) == 0){
      a <- sum(exp(lchoose(x, q)-lchoose(n, q)))
      (1/(1-q))*(a-1)
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term <- qDSub(q, sort.data, n)
      
      a = sum(tab*term)
      b = ifelse(f1 == 0|A == 1, 0, f1/n*(1-A)^(1-n)*(A^(q-1)-sum(qDSub2(q, f1, A, n))))
      (1/(1-q))*(a+b-1)
    }
  }
  sapply(q, Sub)
}

max_est <- function(x, q){
  x = x[x > 0]
  n = sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  A = ifelse(f2>0, 2*f2/((n-1)*f1+2*f2), ifelse(f1>0, 2/((n-1)*(f1-1)+2), 1))
  
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      r <- 1:(n-1)
      a <- sum(x/n*(digamma(n)-digamma(x)))
      b <- ifelse(f1 == 0|A == 1, 0, f1/n*(1-A)^(1-n)*(-log(A)-sum((1-A)^r/r)))
      a+b
    }else if(abs(q-round(q)) == 0){
      a <- sum(exp(lchoose(x, q)-lchoose(n, q)))
      a
    }
    else {
      sort.data = sort(unique(x))
      tab = table(x)
      term <- qDSub(q, sort.data, n)
      
      r <- 0:(n-1)
      a = sum(tab*term)
      term2 <- qDSub2(q, f1, A, n)
      b = ifelse(f1 == 0|A == 1, 0, f1/n*(1-A)^(1-n)*(A^(q-1)-sum(term2)))
      a+b
    }
  }
  sapply(q, Sub)
}

hier.entropy <- function(dat, mat, q, weight = "size", type = 'mle', datatype = 'abundance'){
  rownames(dat) <- paste0("aelle", seq_len(nrow(dat)))
  population <- ncol(mat)
  H <- nrow(mat)
  if(datatype == 'incidence') {
    n <- sum(dat[-1, ])
    s <- nrow(dat[-1, ])
  } else{ 
    n <- sum(dat)
    s <- nrow(dat)
  }
  
  M <- vector("list", H)
  alpha.relative <- numeric(H)
  mat <- apply(mat, 2, as.character)
  index <- rev(lapply(seq_len(H), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x)))))
  if(type == 'est'){
    
    wij <- apply(dat, 2, function(x) sum(x)/n)
    
  } else {
    
    if(weight[1] == 'size' ) wij <- apply(dat, 2, function(x) sum(x)/n)
    if(weight[1] == 'equal') wij <- rep(1/ncol(dat), ncol(dat))
    
  }
  
  wij <- wij/sum(wij)
  w <- lapply(index, function(x) lapply(x, function(y) sum(wij[y])))
  if(datatype == 'abundance') pij <- apply(dat, 2, function(x) x/sum(x))
  if(datatype == 'incidence') pij <- apply(dat, 2, function(x){
    c(x[1], x[-1]/sum(x[-1]))
  } )
  if(datatype == 'incidence') weight.p <- pij*rbind(rep(1, population), t(replicate(s,wij))) else weight.p <- pij*t(replicate(s,wij))
  p <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
  xx <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
  x <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x]))))
  M <- lapply(index, function(x) lapply(x,length))
  FUN2 <- ifelse(datatype == 'abundance', entropy_true, entropy_true_inc)
  for(i in 1:H){
    if(type == 'est'){
      group <- do.call(rbind, xx[[i]])
      if(datatype == 'abundance'){
        alpha.relative[i] <- sum(unlist(w[[i]])*apply(group, 1, function(b) entropy_est(x = b, q)))
      } else {
        alpha.relative[i] <- sum(unlist(w[[i]])*apply(group, 1, function(b) entropy_est(x = b[-1], q)))
      }
    }else {
      group <- do.call(rbind, p[[i]])
      alpha.relative[i] <- sum(unlist(w[[i]])*apply(group, 1, function(b) FUN2(p = b, q)))
    }
  }
  if(q == 1) {
    D.alpha <- exp(alpha.relative)
  } else {
    D.alpha <- (alpha.relative*(1-q)+1)^(1/(1-q))
  }
  D.beta <- D.alpha[-1]/D.alpha[-length(D.alpha)]
  names(D.alpha) <- c(paste0("qD_alpha",seq_len(H-1)),"qD_gamma")
  names(alpha.relative) <- c(paste0("qH_alpha",seq_len(H-1)),"qH_gamma")
  all.pair <- combn(seq_len(H),2)
  if(H == 2){
    pair2 = pair <- all.pair
  } else {
    pair <- cbind(all.pair[,diff(all.pair) == 1], c(1,H))
    pair2 <- pair[ ,-ncol(pair)]
  } 
  temp <- apply(pair, 2, function(j){
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    num <- sapply(index[[j[1]]], function(a) seq_along(w.down)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    w.down.new <-w.down[num]
    ratio <- w.up/w.down.new
    if(q == 1) {
      Max <- -sum(w.up*log(ratio))
      Gamma.Max <- exp(Max)
    }
    else{
      if(type == 'est'){
        X <- do.call(rbind, xx[[j[1]]])
        if(datatype == 'abundance'){
          group.p <- apply(X, 1, function(x) {
            max_est(x, q)
          })
        } else {
          group.p <- apply(X, 1, function(x) {
            max_est(x[-1], q)
          })
        }
      } else {
        if(datatype == 'incidence') dat <- do.call(rbind,p[[j[[1]]]])[ ,-1]
        if(datatype == 'abundance') dat <- do.call(rbind,p[[j[[1]]]])
        if(q==0){
          group.p <- apply(dat, 1, function(N) sum(N>0))
        }else{
          group.p <-  rowSums(dat^q)
        }
      }
      Max <- sum(w.down.new*ratio*(1-ratio^(q-1))*group.p)/(q-1)
      Gamma.Max <- sum(w.up*(ratio^(q-1))*group.p)^(1/(1-q))
    }
    c(differential.relative = diff(alpha.relative[j])/Max, Gamma.Max = Gamma.Max)
  })
  differential.relative <- temp[1, ]
  if(q == 1) D.Beta_max <- temp[2,-length(D.alpha)] else D.Beta_max <- temp[2,-length(D.alpha)] / D.alpha[-length(D.alpha)]
  D.Diss <- Dissimilarity_measure(D.beta, D.Beta_max, q, pair2)
  names(differential.relative ) <- paste0("Delta ","level",apply(pair, 2,paste, collapse = ""))
  ifelse(H == 2, beta <- diff(alpha.relative), beta <- c(diff(alpha.relative), alpha.relative[length(alpha.relative)]-alpha.relative[1]))
  beta[beta < 0] <- 0
  differential.relative[differential.relative < 0] <- 0
  names(beta) <- paste0("qH_Beta ","level",apply(pair, 2,paste, collapse = ""))
  names(D.beta) <- paste0("qD_Beta ", seq_along(D.beta))
  names(D.Beta_max) <- paste0("qD_Beta_max ", seq_along(D.beta))
  result <- data.frame(relative = c(alpha.relative, beta, differential.relative, D.alpha, D.beta, D.Beta_max, D.Diss))
  out <- result
  out
}

Dissimilarity_measure <- function(D.beta, D.Beta_max, q, pair2){
  if(q!=1){
    sorensen <- (D.beta^(1-q)-1)/(D.Beta_max^(1-q)-1)
    Jaccard <- (D.beta^(q-1)-1)/(D.Beta_max^(q-1)-1)
    h <- (D.beta^(-1)-1)/(D.Beta_max^(-1)-1)
    turnover_rate <- (D.beta-1)/(D.Beta_max-1)
  } else {
    sorensen <- log(D.beta)/log(D.Beta_max)
    Jaccard <- sorensen
    h <- (D.beta^(-1)-1)/(D.Beta_max^(-1)-1)
    turnover_rate <- (D.beta-1)/(D.Beta_max-1)
  } 
  out <- c(sorensen, Jaccard, h, turnover_rate)
  names(out) <- c(paste0("1-CqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-UqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-SqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-VqN(", apply(pair2, 2, paste, collapse = '|'), ')'))
  out
}

Candf0Fun <- function(f1, f2, n) {
  if (f2 > 0) {
    C <- 1 - f1 / n * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))
    f0 <- (n - 1) / n * f1^2 / (2 * f2)
  } else if (f2 == 0 & f1 != 0) {
    C <- 1 - f1 / n * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
  } else {
    C <- 1
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
  }
  f0 <- ceiling(f0)
  return(c(C, f0))
}

Boots.population <- function(data, datatype){
  if(datatype == "abundance"){
    N = ncol(data)
    n = colSums(data)
    pool = rowSums(data)
    OBS = sum(pool > 0)
    data = data[pool > 0,]
    obs = sapply(1:N, function(k) sum(data[, k] > 0))
    F1 = sum(pool == 1)
    F2 = sum(pool == 2)
    F0 = round( ifelse(F2==0, F1*(F1-1)/2, F1^2/(2*F2)) )
    
    f1 = sapply(1:N, function(k) sum(data[,k] == 1))
    f2 = sapply(1:N, function(k) sum(data[,k] == 2))
    C =1-f1/n
    
    f0 = round(sapply(1:N, function(k) ifelse(f2[k] == 0, f1[k]*(f1[k]-1)/2, f1[k]^2/(2*f2[k]))))
    r.data = sapply(1:N, function(k) data[, k]/n[k])
    W = sapply(1:N, function(k) (1-C[k])/sum(r.data[, k]*(1-r.data[, k])^n[k]))
    
    ifelse(F0 > 0, boots.pop<-rbind(r.data,matrix(0, ncol=N, nrow=F0)), boots.pop <- r.data)
    
    for(i in 1:N)
    {
      if(f0[i]>0)
      {
        f0[i] = ifelse(f0[i]+obs[i] > OBS+F0, OBS+F0-obs[i], f0[i])
        boots.pop[ ,i][1:OBS] = boots.pop[ ,i][1:OBS]*(1-W[i]*(1-boots.pop[, i][1:OBS])^n[i])
        I = which(boots.pop[, i] == 0)
        II = sample(I, f0[i])
        boots.pop[II, i]=rep((1-C[i])/f0[i], f0[i])
      }
    }
  }else if(datatype == "incidence"){
    data <- as.matrix(data)
    t = data[1,]
    Tt = sum(t)
    data = data[-1,]
    N = ncol(data);
    pool = rowSums(data)
    OBS = length(pool[pool > 0])
    data = data[pool > 0,] 
    obs = sapply(1:N, function(k) length(data[, k][data[, k]>0]))
    Q1 = sum(pool == 1)
    Q2 = sum(pool == 2)
    Q0 = round(((Tt-1)/Tt)*ifelse(Q2 == 0, Q1*(Q1-1)/2, Q1^2/(2*Q2)))
    
    q1 = sapply(1:N, function(k) sum(data[,k] == 1))
    q2 = sapply(1:N, function(k) sum(data[,k] == 2))
    P1 = sapply(1:N, function(k) ifelse(q1[k]+q2[k] == 0, 0, 2*q2[k]/((t[k]-1)*q1[k]+2*q2[k])))
    
    q0 = round(sapply(1:N, function(k) ((t[k]-1)/t[k])*ifelse(q2[k] == 0, q1[k]*(q1[k]-1)/2, q1[k]^2/(2*q2[k]))))
    r.data = sapply(1:N, function(k) data[, k]/t[k])
    W = sapply(1:N, function(k) (1-P1[k])*(q1[k]/t[k])/sum(r.data[, k]*(1-r.data[, k])^t[k]))
    
    boots.pop = if(Q0 > 0) rbind(r.data,matrix(0,ncol=N,nrow=Q0)) else r.data 
    
    for(i in 1:N){
      if(q0[i]>0){
        q0[i] = ifelse(q0[i]+obs[i]>OBS+Q0, OBS+Q0-obs[i], q0[i])
        boots.pop[,i][1:OBS] = boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[, i][1:OBS])^t[i])
        I = which(boots.pop[,i] == 0)
        II = sample(I, q0[i])
        boots.pop[II, i] = rep((q1[i]/t[i])/q0[i], q0[i])
      }
    }
  }
  return(boots.pop)
}

Decomposition_Bootstrap <- function(data, mat, q, weight, nboot, conf = 0.95, type = "mle", datatype){
  if(is.null(nboot) | is.nan(nboot) | is.na(nboot) ) nboot <- 200
  m <- if(length(unique(unlist(mat[1,]))) == 1) nrow(mat)-1 else nrow(mat)
  N <- ncol(data)
  bootstrap <- function(i){
    boots.pop <- Boots.population(data, datatype)
    if (datatype == "incidence"){
      TT <- unlist(data[1, ])
      a <- sapply(1:ncol(boots.pop), function(i) {
        sapply(boots.pop[ ,i], function(p) rbinom(n = 1, size = TT[i], prob = p))
      })
      rbind(TT, a)
    }
    boots.data <- sapply(1:N, function(x) rmultinom(1, sum(data[, x]), boots.pop[, x]))
    
    
    out <- as.matrix(hier.entropy(boots.data, mat, q, weight, type = type, datatype = datatype))
    return(out)
  }
  
  if(nboot <= 0 | nboot == 1 ) {
    result <- sapply(1:3, bootstrap, simplify = 'array')
    result <- array(NA, dim = dim(result),dimnames = dimnames(result))
    est_mean <- apply(result, c(1:2), mean)
    confidence <- c((1-conf)/2, (1+conf)/2)
    est_LCL <- apply(result, c(1:2), function(x) quantile(x, confidence[1],na.rm = T))
    est_UCL <- apply(result, c(1:2), function(x) quantile(x, confidence[2],na.rm = T))
    sd <- apply(result, c(1:2), function(x) sd(x, na.rm = T))
    out <- array(c(est_mean, sd, est_LCL-est_mean, est_UCL-est_mean), dim = c(dim(est_mean), 4))
    dimnames(out) <- list(dimnames(sd)[[1]], dimnames(sd)[[2]], c('est_mean', 'sd', 'est_LCL', 'est_UCL'))
    out
  }else{
    result <- sapply(1:nboot, bootstrap, simplify = 'array')
    est_mean <- apply(result, c(1:2), mean)
    confidence <- c((1-conf)/2, (1+conf)/2)
    est_LCL <- apply(result, c(1:2), function(x) quantile(x, confidence[1],na.rm = T))
    est_UCL <- apply(result, c(1:2), function(x) quantile(x, confidence[2],na.rm = T))
    sd <- apply(result, c(1:2), function(x) sd(x, na.rm = T))
    out <- array(c(est_mean, sd, est_LCL-est_mean, est_UCL-est_mean), dim = c(dim(est_mean), 4))
    dimnames(out) <- list(dimnames(sd)[[1]], dimnames(sd)[[2]], c('est_mean', 'sd', 'est_LCL', 'est_UCL'))
    out
  }
}

hier.taxonomy_rel <- function(data, mat, q, weight = "size", nboot = 200, conf = 0.95, type = "mle", datatype){
  t <- apply(mat,MARGIN = 1, function(x) length(unique(x)))
  if(sum((t[-1]-t[-length(t)])<0)>0) stop("number of lower level must be large than or equal to the number of higher level, please renew your structure matrix.")
  if(ncol(data)!=ncol(mat)) stop('data does not match the structure matrix')
  if(nrow(mat)==1) stop('structure matrix has only one layer')
  if(nrow(data)<2) stop('please key in bigger than 2 number of species ')
  if(sum(round(data)!=data)>0) nboot=0
  if(length(unique(unlist(mat[1,])))!=1) stop("the highest level of structure matrix must have only one community")
  sub <- function(q){
    Boot <- Decomposition_Bootstrap(data, mat, q, weight, nboot, conf, type, datatype)
    est <- hier.entropy(dat = data, mat = mat, q = q, weight = weight, type = type, datatype = datatype)
    relative <- Boot[ ,1, ]
    LB <- est$relative-relative[ ,2]*qnorm((1+conf)/2)
    LB[LB<0] <- 0
    out_rel <- data.frame(est$relative, relative[ ,2], LB, est$relative+relative[ ,2]*qnorm((1+conf)/2))
    colnames(out_rel)  <- c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
    return(round(out_rel, 4))
  }
  output <- lapply(q, sub) %>% do.call(rbind, .)
  
  Order.q = rep(q, each = 24)
  
  Method = rep(rownames(output)[1:24], length(q))
  
  output = cbind(Method, Order.q, output, Decomposition = "relative")
  
  rownames(output) = NULL
  
  return(output)
}

entropy_true_abs <- function(p, Mi, q) {
  sub <- function(q){
    if(q!=1) {
      tmp <- replicate(ncol(p), unlist(Mi))
      pp <- p/tmp
      x <- c(tmp*t(apply(pp, 1, function(x) ifelse(x == 0, 0, x^q))))
      (sum(x))^(1/(1-q))/sum(unlist(Mi))
    } else {
      x <- c(p)
      x = x[x > 0]
      a <- -sum(x*log(x))
      c <- sum(sapply(1:nrow(p), function(i) sum(p[i, ])*log(Mi[[i]])))
      exp(a+c)/sum(unlist(Mi))
    }
  }
  sapply(q, sub)
}

entropy_est_abs <- function(xij, Mi, wij, q){
  x <- c(xij)
  x = x[x > 0]
  n = sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  A = ifelse(f2>0, 2*f2/((n-1)*f1+2*f2), ifelse(f1>0, 2/((n-1)*(f1-1)+2), 1))
  Sub <- function(q){
    if(q==1){
      
      a <- sum(apply(xij, 1, max_est, q=q) * unlist(wij))
      b <- sum( unlist(wij) * (  log(unlist(Mi)) - log(unlist(wij))) )
      # r <- 1:(n-1)
      # a <- sum(x/n*(digamma(n)-digamma(x)))
      # b <- ifelse(f1 == 0|A == 1, 0, f1/n*(1-A)^(1-n)*(-log(A)-sum((1-A)^r/r)))
      # c <- sum(sapply(1:nrow(xij), function(i) wij[[i]]*log(Mi[[i]])))
      
      
      exp(a+b) / sum(unlist(Mi))
    } else {
      tmp <- apply(xij, 1, max_est, q=q)
      tmp2 <- rowSums(xij)/n
      sum(unlist(Mi)^(1-q)*tmp2^q*tmp)^(1/(1-q)) / sum(unlist(Mi))
    }
  }
  #   if(q==0){
  #     sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)-1
  #   }
  #   else if(q==1){
  #     r <- 1:(n-1)
  #     a <- sum(x/n*(digamma(n)-digamma(x)))
  #     b <- ifelse(f1 == 0|A == 1, 0, f1/n*(1-A)^(1-n)*(-log(A)-sum((1-A)^r/r)))
  #     c <- sum(sapply(1:nrow(xij), function(i) wij[[i]]*log(Mi[[i]])))
  #     exp(a+b+c) / sum(unlist(Mi))
  #   }else if(abs(q-round(q)) == 0){
  #     a <- sum(exp(lchoose(x, q)-lchoose(n, q)))
  #     #ifelse(A==0,NA,A^(1/(1-q)))
  #     1/length(Mi)*a^(1/(1-q))
  #   }else {
  #     sort.data = sort(unique(x))
  #     tab = table(x)
  #     term = sapply(sort.data, function(z){
  #       k=0:(n-z)
  #       sum(choose(k-q, k)*exp(lchoose(n-k-1, z-1)-lchoose(n, z)))
  #     })
  #     r <- 0:(n-1)
  #     a = sum(tab*term)
  #     b = ifelse(f1 == 0|A == 1, 0, f1/n*(1-A)^(1-n)*(A^(q-1)-sum(choose(q-1, r)*(A-1)^r)))
  #     1/length(Mi)*(a+b)^(1/(1-q))
  #   }
  # }
  sapply(q, Sub)
}

hier.entropy_abs <- function(dat, mat, q, type = 'mle', datatype = 'abundance'){
  if(datatype == 'incidence') dat <- dat[-1, ]
  rownames(dat) <- paste0("aelle", seq_len(nrow(dat)))
  population <- ncol(mat)
  H <- nrow(mat)
  n <- sum(dat)
  s <-nrow(dat)
  M <- vector("list", H)
  alpha.relative <- numeric(H)
  mat <- apply(mat, 2, as.character)
  index <- rev(lapply(seq_len(H), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x)))))
  wij <- apply(dat, 2, function(x) sum(x)/n)
  wij <- wij/sum(wij)
  w <- lapply(index, function(x) lapply(x, function(y) sum(wij[y])))
  if(n == population) pij <- dat
  else pij <- apply(dat, 2, function(x) x/sum(x))
  weight.p <- pij*t(replicate(s,wij))
  p <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
  xx <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
  x <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x]))))
  M <- lapply(index, function(x) lapply(x,length))
  for(i in 1:H){
    if(type == 'est'){
      group <- do.call(rbind, xx[[i]])
      alpha.relative[i] <- entropy_est_abs(xij = group, Mi = M[[i]], wij = w[[i]], q)
    } else {
      group2 <- do.call(rbind, x[[i]])
      alpha.relative[i] <- entropy_true_abs(p = group2, Mi = M[[i]], q = q)
    }
  }
  
  if(q!=1) H.alpha <- (1-(population^(q-1)*(alpha.relative*population)^(1-q)))/(q-1) else H.alpha <- log(alpha.relative)
  names(alpha.relative) <- c(paste0("qD_alpha",seq_len(H-1)),"qD_gamma")
  names(H.alpha) <- c(paste0("qH_alpha",seq_len(H-1)),"qH_gamma")
  D.beta <- alpha.relative[-1]/alpha.relative[-length(alpha.relative)]
  all.pair <- combn(seq_len(H),2)
  if(H == 2){
    pair2 = pair <- all.pair
  } else {
    pair <- cbind(all.pair[,diff(all.pair) == 1], c(1,H))
    pair2 <- pair[ ,-ncol(pair)]
  } 
  
  differential.relative <- apply(pair, 2, function(j){
    up <- ifelse(q!=1, -diff(alpha.relative[j]^(1-q))*population^(1-q), diff(log(alpha.relative[j])))
    M1 <- unlist(M[[j[1]]])
    M2 <- unlist(M[[j[2]]])
    num <- sapply(index[[j[1]]], function(a) seq_along(M2)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    if(q == 1 ) {
      Max <- sum(unlist(w[[j[2]]])*log(M2))-sum(unlist(w[[j[1]]])*log(M1))
      Gamma.Max <- exp(Max)
    }else{
      if(type == 'est'){
        X <- do.call(rbind, xx[[j[1]]])
        tmp <- apply(X, 1, max_est, q = q)
        tmp2 <- unlist(w[[j[1]]])^q
        tmp3 <- unlist(sapply(M2, function(x) rep(x,x)))[1:length(M1)]
        Max <- sum((M1^(1-q)-tmp3^(1-q))*tmp2*tmp)
        Gamma.Max <- sum(tmp3^(1-q)*tmp2*tmp)^(1/(1-q)) / population
      } else {
        if(q==0){
          group.p <- apply(do.call(rbind,x[[j[[1]]]]),1,function(N) sum(N>0))
        }else{
          group.p <-  rowSums(do.call(rbind,x[[j[[1]]]])^q)
        }
        Max <- sum((M1^(1-q)-M2[num]^(1-q))*group.p)
        Gamma.Max <- sum(M2[num]^(1-q)*group.p)^(1/(1-q)) / population
      }
    }
    c(up = up, Max = Max, Gamma.Max = Gamma.Max)
  })
  differential <- differential.relative[1, ] / differential.relative[2, ]
  if(q == 1) D.Beta_max <- differential.relative[3,-length(alpha.relative)] else D.Beta_max <- differential.relative[3,-length(alpha.relative)] / alpha.relative[-length(alpha.relative)]
  D.Diss <- Dissimilarity_measure(D.beta, D.Beta_max, q, pair2)
  names(differential) <- paste0("Delta ","level(",apply(pair, 2,paste, collapse = ""), ")")
  ifelse(H == 2, beta <- diff(H.alpha), beta <- c(diff(H.alpha), H.alpha[length(H.alpha)]-H.alpha[1]))
  beta[beta < 0] <- 0
  names(beta) <- paste0("qH_Beta ","level",apply(pair, 2,paste, collapse = ""))
  names(D.beta) <- paste0("qD_Beta ", seq_along(D.beta))
  names(D.Beta_max) <- paste0("qD_Beta_max ", seq_along(D.beta))
  result <- data.frame(relative = c(H.alpha, beta, differential, alpha.relative, D.beta, D.Beta_max, D.Diss))
  out <- result
  out
}

Boots.population2 <- function(data, datatype){
  if(datatype == "abundance"){
    N = ncol(data)
    n = colSums(data)
    pool = rowSums(data)
    OBS = sum(pool > 0)
    data = data[pool > 0,]
    obs = sapply(1:N, function(k) sum(data[, k] > 0))
    F1 = sum(pool == 1)
    F2 = sum(pool == 2)
    F0 = round( ifelse(F2==0, F1*(F1-1)/2, F1^2/(2*F2)) )
    
    f1 = sapply(1:N, function(k) sum(data[,k] == 1))
    f2 = sapply(1:N, function(k) sum(data[,k] == 2))
    C =1-f1/n
    
    f0 = round(sapply(1:N, function(k) ifelse(f2[k] == 0, f1[k]*(f1[k]-1)/2, f1[k]^2/(2*f2[k]))))
    r.data = sapply(1:N, function(k) data[, k]/n[k])
    W = sapply(1:N, function(k) (1-C[k])/sum(r.data[, k]*(1-r.data[, k])^n[k]))
    
    ifelse(F0 > 0, boots.pop<-rbind(r.data,matrix(0, ncol=N, nrow=F0)), boots.pop <- r.data)
    
    for(i in 1:N)
    {
      if(f0[i]>0)
      {
        f0[i] = ifelse(f0[i]+obs[i] > OBS+F0, OBS+F0-obs[i], f0[i])
        boots.pop[ ,i][1:OBS] = boots.pop[ ,i][1:OBS]*(1-W[i]*(1-boots.pop[, i][1:OBS])^n[i])
        I = which(boots.pop[, i] == 0)
        II = sample(I, f0[i])
        boots.pop[II, i]=rep((1-C[i])/f0[i], f0[i])
      }
    }
  }else if(datatype == "incidence"){
    data <- as.matrix(data)
    t = data[1,]
    Tt = sum(t)
    data = data[-1,]
    N = ncol(data);
    pool = rowSums(data)
    OBS = length(pool[pool > 0])
    data = data[pool > 0,] 
    obs = sapply(1:N, function(k) length(data[, k][data[, k]>0]))
    Q1 = sum(pool == 1)
    Q2 = sum(pool == 2)
    Q0 = round(((Tt-1)/Tt)*ifelse(Q2 == 0, Q1*(Q1-1)/2, Q1^2/(2*Q2)))
    
    q1 = sapply(1:N, function(k) sum(data[,k] == 1))
    q2 = sapply(1:N, function(k) sum(data[,k] == 2))
    P1 = sapply(1:N, function(k) ifelse(q1[k]+q2[k] == 0, 0, 2*q2[k]/((t[k]-1)*q1[k]+2*q2[k])))
    
    q0 = round(sapply(1:N, function(k) ((t[k]-1)/t[k])*ifelse(q2[k] == 0, q1[k]*(q1[k]-1)/2, q1[k]^2/(2*q2[k]))))
    r.data = sapply(1:N, function(k) data[, k]/t[k])
    W = sapply(1:N, function(k) (1-P1[k])*(q1[k]/t[k])/sum(r.data[, k]*(1-r.data[, k])^t[k]))
    
    boots.pop = if(Q0 > 0) rbind(r.data,matrix(0,ncol=N,nrow=Q0)) else r.data 
    
    for(i in 1:N){
      if(q0[i]>0){
        q0[i] = ifelse(q0[i]+obs[i]>OBS+Q0, OBS+Q0-obs[i], q0[i])
        boots.pop[,i][1:OBS] = boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[, i][1:OBS])^t[i])
        I = which(boots.pop[,i] == 0)
        II = sample(I, q0[i])
        boots.pop[II, i] = rep((q1[i]/t[i])/q0[i], q0[i])
      }
    }
  }
  return(boots.pop)
}

Decomposition_Bootstrap2 <- function(data, mat, q, nboot, conf = 0.95, type = "mle", datatype){
  if(is.null(nboot) | is.nan(nboot) | is.na(nboot) ) nboot <- 200
  m <- if(length(unique(unlist(mat[1,]))) == 1) nrow(mat)-1 else nrow(mat)
  N <- ncol(data)
  bootstrap <- function(i){
    boots.pop <- Boots.population(data, datatype)
    if (datatype == "incidence"){
      TT <- unlist(data[1, ])
      a <- sapply(1:ncol(boots.pop), function(i) {
        sapply(boots.pop[ ,i], function(p) rbinom(n = 1, size = TT[i], prob = p))
      })
      rbind(TT, a)
    }
    boots.data <- sapply(1:N, function(x) rmultinom(1, sum(data[, x]), boots.pop[, x]))
    
    out <- as.matrix(hier.entropy_abs(boots.data, mat, q, type = type, datatype = datatype))
    return(out)
  }
  if(nboot <= 0 | nboot == 1 ) {
    result <- sapply(1:3, bootstrap, simplify = 'array')
    result <- array(NA, dim = dim(result),dimnames = dimnames(result))
    est_mean <- apply(result, c(1:2), mean)
    confidence <- c((1-conf)/2, (1+conf)/2)
    est_LCL <- apply(result, c(1:2), function(x) quantile(x, confidence[1],na.rm = T))
    est_UCL <- apply(result, c(1:2), function(x) quantile(x, confidence[2],na.rm = T))
    sd <- apply(result, c(1:2), function(x) sd(x, na.rm = T))
    out <- array(c(est_mean, sd, est_LCL-est_mean, est_UCL-est_mean), dim = c(dim(est_mean), 4))
    dimnames(out) <- list(dimnames(sd)[[1]], dimnames(sd)[[2]], c('est_mean', 'sd', 'est_LCL', 'est_UCL'))
    out
  }else{
    result <- sapply(1:nboot, bootstrap, simplify = 'array')
    est_mean <- apply(result, c(1:2), mean)
    confidence <- c((1-conf)/2, (1+conf)/2)
    est_LCL <- apply(result, c(1:2), function(x) quantile(x, confidence[1],na.rm = T))
    est_UCL <- apply(result, c(1:2), function(x) quantile(x, confidence[2],na.rm = T))
    sd <- apply(result, c(1:2), function(x) sd(x, na.rm = T))
    out <- array(c(est_mean, sd, est_LCL-est_mean, est_UCL-est_mean), dim = c(dim(est_mean), 4))
    dimnames(out) <- list(dimnames(sd)[[1]], dimnames(sd)[[2]], c('est_mean', 'sd', 'est_LCL', 'est_UCL'))
    out
  }
}

hier.taxonomy_abs <- function(data, mat, q, nboot = 200, conf = 0.95, type = "mle", datatype){
  t <- apply(mat,MARGIN = 1, function(x) length(unique(x)))
  if(sum((t[-1]-t[-length(t)])<0)>0) stop("number of lower level must be large than or equal to the number of higher level, please renew your structure matrix.")
  if(ncol(data)!=ncol(mat)) stop('data does not match the structure matrix')
  if(nrow(mat)==1) stop('structure matrix has only one layer')
  if(nrow(data)<2) stop('please key in bigger than 2 number of species ')
  if(sum(round(data)!=data)>0) nboot=0
  if(length(unique(unlist(mat[1,])))!=1) stop("the highest level of structure matrix must have only one community")
  sub <- function(q){
    Boot <- Decomposition_Bootstrap2(data, mat, q, nboot, conf, type, datatype)
    est <- hier.entropy_abs(dat = data, mat = mat, q = q, type = type, datatype = datatype)
    relative <- Boot[ ,1, ]
    LB <- est$relative-relative[ ,2]*qnorm((1+conf)/2)
    LB[LB<0] <- 0
    out_rel <- data.frame(est$relative, relative[ ,2], LB, est$relative+relative[ ,2]*qnorm((1+conf)/2))
    colnames(out_rel)  <- c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
    return(round(out_rel, 4))
  }
  output <- lapply(q, sub) %>% do.call(rbind, .)
  
  Order.q = rep(q, each = 24)
  
  Method = rep(rownames(output)[1:24], length(q))
  
  output = cbind(Method, Order.q, output, Decomposition = "absolute")
  
  rownames(output) = NULL
  
  return(output)
}

#Phylogeny
choose_data = function(data, tree){
  
  tmp <- data[names(tree$leaves)]  
  for(i in 1:length(tree$parts)){
    tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
    names(tmp)[length(tmp)] <- names(tree$parts)[i]
  }
  tmp <- data.frame('branch_abun' = tmp,"branch_length" = c(tree$leaves,tree$nodes))
  
  return(tmp)
}

TranMul = function(data, tree){
  rtree = newick2phylog(convToNewick(tree))
  data = cbind(data[names(rtree$leaves), ])
  rdata = apply(data, 1, sum)
  AlphaTmp = list()
  GammaTmp = choose_data(rdata, rtree)
  GammaTbar = sum(GammaTmp[, 1]*GammaTmp[, 2]) / sum(rdata)
  for(i in 1:ncol(data)){
    adata = aadata = data[, i]
    names(aadata) = rownames(data)
    names(adata) = tree$tip.label
    if(length(rtree$leaves)<=2){
      abun <- c(adata[names(adata)%in%rtree$tip.label], root = sum(adata))
      tmp = data.frame('branch_abun'=abun, "branch_length" = c(rtree$edge.length,0))
    } else{
      tmp = choose_data(aadata, rtree)
    }
    AlphaTmp[[i]] = tmp
  }
  output = list(Alpha=AlphaTmp, Gamma=GammaTmp)
  return(output)
}

mle.phy.q <- function(pp,LL,TT,q){
  LL <- LL[pp>0]
  pp <- pp[pp>0]
  if(q!=1){
    LL%*%((pp/TT)^q)
    #LL*((pp/T)^q)
  }else{
    -LL%*%(pp*log(pp))
    #-LL*((pp/T)*log(pp/T))
  }
}

delta = function(data, k, n, A){
  ans = sapply(1:length(k), function(i){
    if(k[i]<n){      
      data1 = data[data[,1]<=(n-k[i]) &  data[,1]>=1,]
      if( class(data1) == "numeric" ) data1 = t(as.matrix(data1))
      sum( data1[,2]*(data1[,1]/n)*exp(lchoose(n-data1[,1], k[i])-lchoose(n-1, k[i])) ) 
    }else{
      g1 = sum(data[data==1,2])
      g1*(1-A)^(k[i]-n+1)/n 
    }
  })
  return( ans )
}

est.phy.q <-  function(xx, LL, TT, n, q, rtreephy){ # proposed
  LL <- LL[names(xx)]
  tmp <- data.frame(abun = xx, length = LL)
  PD_obs <- sum(tmp[tmp[,1]>0,2])
  #f1 <- ifelse(datatype=='incidence', sum(rowSums(data)==1), sum(data==1) )
  f1 = sum(tmp[, 1]==1)
  f2 = sum(tmp[, 1]==2)
  g1 = sum(tmp[tmp[, 1]==1, 2])  
  g2 = sum(tmp[tmp[, 1]==2, 2]) 
  node <- names(rtreephy$parts)
  node2 = names(xx)[xx==2][names(xx)[xx==2] %in% node]
  if(length(node2)>0){
    rep2 <- sapply(node2, function(tx){
      sum(xx[rtreephy$parts[[tx]]] %in% c(2,0)) == length(rtreephy$parts[[tx]]) 
    })
    f2 <- f2 - sum(rep2)
  }
  node1 = names(xx)[xx==1][names(xx)[xx==1] %in% node]
  if(length(node1)>0){
    rep1 <- sapply(node1, function(tx){
      sum(xx[rtreephy$parts[[tx]]] %in% c(1,0)) == length(rtreephy$parts[[tx]]) 
    })
    f1 <- f1 - sum(rep1)
  }
  if(f2 > 0){
    A = 2*f2/((n-1)*f1+2*f2)
  }else if(f2 == 0 & f1 > 0){
    A = 2/((n-1)*(f1-1)+2)  
  }else{
    A = 1 
  }
  t_bar <- sum(xx*LL/n)
  # position0 <- which(q==0)
  # position1 <- which(q==1)
  # position_else = c(1:length(q))[-c(position1)]
  #position_else = c(1:length(q))[-c(position0,position1)]
  #ans = rep(0, length(q))
  if(q==0){
    ans = PD_obs+ifelse(g2>0, (n-1)/n*g1^2/(2*g2), (n-1)/n*g1*(f1-1)/2*(f2-1))
  }else if(q==1){
    q1 = sum(sapply(1:(n-1), function(r) {(1-A)^r/r} ))
    if(A < 1) h2 = (g1/n)*((1-A)^(-n+1))*(-log(A)-q1)
    if(A == 1) h2 = 0
    tmp2 = subset(tmp, tmp[,1]>=1 & tmp[,1]<=(n-1) )
    tmp2 = cbind(tmp2, sapply(tmp2[,1], function(x) sum( 1/(x:(n-1)))))
    h1 = sum(apply(tmp2, 1, prod))/n
    h = h1+h2
    #ans[position1] = t_bar*exp(h/t_bar)
    ans = h
  } else{
    r = 0 : (n-1)
    de = delta(tmp, r , n , A)
    a = sum( choose(q-1, r)*(-1)^r*de )/((t_bar)^q)
    if(A < 1) b = (g1*((1-A)^(1-n))/n)*(A^(q-1)-sum(choose(q-1, r)*(A-1)^r))/((t_bar)^q)
    if(A == 1) b = 0
    ans = a+b
  }
  return( ans )
}

Boots.pop = function(data, rtree, tmp){
  # if(datatype == "abundance"){
  N = ncol(data); n = colSums(data)
  pool=rowSums(data) ; OBS=length(pool)
  rtreephy <- newick2phylog(convToNewick(rtree))
  OBS_B <- dim(tmp[[1]])[1]
  obs <- colSums(data>0)
  TT <- sum(tmp[[1]][,1]/n[1]*tmp[[1]][,2])
  F1=sum(pool==1);F2=sum(pool==2)
  F0=ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2))*(sum(n)-1)/sum(n)  #pool assemblage f0 estimate
  F0_N <- round(F0)
  f1=sapply(1:N,function(k) sum(data[,k]==1))
  f2=sapply(1:N,function(k) sum(data[,k]==2))
  g1=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==1,2]) ))
  g2=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==2,2]) ))
  C <- ifelse(f2 == 0, C <- 1 - f1/n*(n-1)*(f1-1)/((n-1)*(f1-1) + 2), C <- 1 - f1/n*(n-1)*f1/((n-1)*f1 + 2*f2))
  f0 <- ifelse(f2 == 0, f0 <- f1*(f1-1)/2, f0 <- f1^2/(2*f2))*(n-1)/n
  f0_N <- round(f0)
  r.data=sapply(1:N,function(k) data[,k]/n[k])
  W=sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k]))
  g0 = sapply(1:N, function(k)  
    if((2*g2[k]*f1[k])>g1[k]*f2[k]) (n[k]-1)/n[k]*g1[k]^2/2/g2[k]
    else (n[k]-1)/n[k]*g1[k]*(f1[k]-1)/2/(f2[k]+1) )
  if(F0>0){boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=F0_N))
  }else{boots.pop=r.data}
  L = matrix(0, nrow=(OBS_B+F0_N), ncol=N)
  boots.pop2 = matrix(0, nrow=(OBS_B+F0_N), ncol=N) #obs branch abundance + undetected
  for(i in 1:N)
  { 
    if(f0_N[i]>0)
    {
      f0_N[i]=ifelse(f0_N[i]+obs[i]>OBS+F0_N, OBS+F0_N-obs[i],f0_N[i])
      boots.pop[,i][1:OBS] <- r.data[ ,i]*(1-W[i]*(1-r.data[ ,i])^n[i])   #species
      I=which(boots.pop[,i]==0) #ai+bi
      II=sample(I,f0_N[i])
      u.p <- (1-C[i])/f0_N[i]
      boots.pop[II,i]=rep(u.p,f0_N[i])
      da = boots.pop[1:OBS,i] #corrected observed relative species abundance
      names(da) = rownames(data)
      mat = choose_data(da, rtreephy) #corrected observed relative node abundance
      boots.pop2[,i]=c(mat[,1], boots.pop[,i][-(1:OBS)]) #add undetected species 
      F00 = sum(II > OBS) #not detect in pool assemblage 
      L[1:nrow(mat),i] = mat[,2] 
      if(F00>0 ){
        index = which(boots.pop2[,i] > 0)[which(boots.pop2[,i] > 0) > nrow(mat)]
        #un.sp <- rownames(data)[II[II<OBS]]
        #g0r <- g0[i]- sum(mat[un.sp,2])
        L[index, i] = g0[i]/F00
      }
    }else{
      L[seq_len(OBS_B), i] = tmp[[i]][ ,2] 
      boots.pop2[seq_len(OBS_B), i] = tmp[[i]][ ,1]
    }
  }
  if(F0_N==0){
    rownames(L) <- rownames(mat)
    rownames(boots.pop2) <- rownames(mat)
  } else{
    rownames(L) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
    rownames(boots.pop2) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
  }
  L[L>TT] <- TT
  return(list(p=boots.pop2,L=L,unseen=F0_N))
  # }
}

bootstrap.q.Beta = function(data, mat, rtree, tmp, q, nboot, wij, type, method){
  H <- nrow(mat)
  out = array(0, dim=c(8*H, length(q), nboot))
  pool <- rowSums(data)
  rtreephy <- newick2phylog(convToNewick(rtree))
  #if(datatype == "abundance"){
  n = colSums(data) ; N = ncol(data)
  pop = Boots.pop(data, rtree, tmp$Alpha)
  S = nrow(data)
  #B = length(c(phytree$leaves,phytree$parts))
  if(pop$unseen == 0) p = pop$p[1:S,]
  if(pop$unseen != 0) p = pop$p[c(1:S,tail(1:nrow(pop$p), pop$unseen)),]
  boot.data = array(0, dim = dim(p))
  rownames(boot.data) <- rownames(p)
  S = nrow(data)
  B = S+rtree$Nnode
  for(i in 1:nboot){
    L = pop$L
    if(pop$unseen == 0) p = pop$p[1:S,]
    if(pop$unseen != 0) p = pop$p[c(1:S,(B+1):nrow(pop$p)),]
    boot.data = array(0, dim = dim(p))
    for(j in 1:ncol(p)) boot.data[,j] = rmultinom(1,n[j],p[,j]) 
    rownames(boot.data) <- rownames(p)
    unseen = boot.data[-(1:S),]
    boot.data.obs <- boot.data[1:S,]
    #boot.data.obs <- boot.data.obs[rowSums(boot.data.obs)>0, ]
    #tip.boot <- names(rtreephy$leaves)[!names(rtreephy$leaves)%in%rownames(boot.data.obs)]
    #rtree.boot <- drop.tip(rtree, tip.boot)
    #rtreephy.boot <- newick2phylog(convToNewick(rtree.boot))
    boot.datatmp = apply(boot.data.obs, 2, function(x){
      #names(x) = names(rtreephy$leaves)
      tmp <- choose_data(x, rtreephy)
      abun <- tmp[ ,1]
      names(abun) <- rownames(tmp)
      return(abun)
    })  
    boot.datatmp = rbind(boot.datatmp, unseen)
    boot.gamma = rowSums(boot.datatmp)
    #boot.gamma = boot.gamma[boot.gamma]
    #boot.datatmp <- boot.datatmp[names(boot.gamma), ]
    boot.gamma <- boot.gamma[boot.gamma>0]
    boot.datatmp <- boot.datatmp[names(boot.gamma), ]
    L <- L[names(boot.gamma), ]
    L.gamma = rowSums(boot.datatmp[names(boot.gamma), ] * L) / boot.gamma
    boot.alpha <- apply(boot.datatmp, 2, function(s) data.frame(branch_abun = s, branch_length = L.gamma))
    boot.gamma <- data.frame(branch_abun = boot.gamma, branch_length = L.gamma)
    boot.tmp <- list(Gamma = boot.gamma, Alpha = boot.alpha)
    #sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, type))
    out[,,i] = sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, type))
  }
  #print(sum(is.infinite(apply(out, 3, sum))))
  #out[ , ,!is.infinite(apply(out, 3, sum))]
  #}
  return(out)
}

transconf = function(Bresult, est, conf){
  est.btse = sd(Bresult)
  est.LCL = est - qnorm(1-(1-conf)/2) * est.btse 
  est.UCL = est + qnorm(1-(1-conf)/2) * est.btse
  # if(any(est.LCL<0)) est.LCL[est.LCL<0] <- 0
  # if(any(est.UCL>1)) est.UCL[est.UCL>1] <- 1
  cbind(est = est, btse=est.btse, LCL=est.LCL, UCL = est.UCL)
}

phy.H.rel <- function(dat, mat, tmp, q, rtreephy, wij, type){
  n <- sum(dat)
  nsite <- ncol(mat)
  H <- nrow(mat)
  index <- lapply(rev(seq_len(H)), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x))))
  w <- lapply(index, function(x) sapply(x, function(y) sum(wij[y])))
  M <- lapply(index, function(x) sapply(x,length))
  ga <- tmp$Gamma[,1];gB <- tmp$Gamma[,2]
  names(gB) <- rownames(tmp$Gamma)
  gp=ga/n;TT=sum(gp*gB);
  aa <- sapply(tmp$Alpha, function(X){
    abun <- X[,1]
    names(abun) = rownames(X)
    return(abun)
  })
  B <-length(ga)
  pij <- sapply(1:nsite, function(x) {
    abun <- aa[,x]
    n <- sum(dat[,x])
    return(abun/n)
  }) #pk|im
  weight.p <- pij*t(replicate(B,wij))
  p <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
  x <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
  xtmp <- lapply(index, function(h) sapply(h, function(x) rowSums(cbind(aa[ ,x]))))
  if(type == "mle"){
    qD <- sapply(1:H, function(i) apply(p[[i]], 2, get("mle.phy.q"), LL = gB, TT = TT, q)) ##\sum{p^q} or -\sum{p*log(p)})
  } else if(type == "est"){
    qD <- sapply(1:H, function(i) sapply(seq_len(ncol(x[[i]])), function(j){
      n <- sum(x[[i]][,j])
      get("est.phy.q")(xx = xtmp[[i]][,j], LL = gB, TT = TT, n = n, q = q, rtreephy)
    }))
  }
  alpha.relative <- sapply(1:H, function(i){
    qD_est <- qD[[i]]
    if(q!=1){
      (w[[i]]%*%qD_est)^(1/(1-q))
    } else{
      exp(w[[i]]%*%(qD_est/TT+log(TT)))
    }
  })
  names(alpha.relative) <- c(paste0("qPD_alpha",seq_len(H-1)),"qPD_gamma")
  all.pair <- cbind(sapply(1:(H-1), function(i) cbind(i,i+1)), c(1,H))
  Ialpha.relative <- sapply(alpha.relative ,function(x){
    if(q!=1){
      alpha.relative.H.j <- (TT-x^(1-q)*TT^q)/(q-1)
    } else{
      alpha.relative.H.j <- (log(x)-log(TT))*TT
    } 
    alpha.relative.H.j}
  )
  names(Ialpha.relative) <- c(paste0("qI_alpha", seq_len(H-1)), "qI_gamma")
  Ibeta.relative <- apply(all.pair, 2, function(j){
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q)/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)
  })
  names(Ibeta.relative) <- c(paste0("qI_Beta ", "level(",apply(all.pair, 2,paste, collapse = "|"), ")"))
  differential.relative <-  apply(all.pair, 2, function(j){
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    num <- sapply(index[[j[1]]], function(a) seq_along(w.down)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    w.down.new <-w.down[num]
    ratio<- w.up/w.down.new
    est <- qD[[j[1]]]*TT^q
    Max <- ifelse(q!=1, sum(w.down.new*ratio*(1-ratio^(q-1))*est)/(q-1),
                  -TT*sum(w.up*log(ratio)))
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q)/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)/Max
    #sum(w.down*sapply(index[[j[2]]], function(g) sum(ratio[g]*(1-ratio[g]^(q-1))*group.p[g])))
  })
  names(differential.relative ) <- paste0("Delta ","level(",apply(all.pair, 2,paste, collapse = "|"), ")")
  beta <- apply(all.pair, 2, function(j) alpha.relative[j[2]]/alpha.relative[j[1]])
  betamax <- apply(all.pair, 2, function(j){
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    num <- sapply(index[[j[1]]], function(a) seq_along(w.down)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    w.down.new <-w.down[num]
    ratio<- w.up/w.down.new
    if(q==1){
      out <- exp(-sum(w.up*log(ratio)))
    } else{
      out <- (sum(w.down.new*(ratio)^q*qD[[j[[1]]]])/sum(w.down.new*(ratio)*qD[[j[[1]]]]))^(1/(1-q))
    }
    out
  })
  # gammamax <- betamax*apply(all.pair, 2, function(j){
  #   alpha.relative[j[1]]
  # })
  # gamma <- apply(all.pair, 2, function(j){''
  #   alpha.relative[j[2]]
  # })
  beta <-  beta[1:(H-1)]
  betamax <- betamax[1:(H-1)] 
  # gamma <-  gamma[1:(H-1)]
  # gammamax <- gammamax[1:(H-1)] 
  ifelse(q==1, C_1 <- log(beta)/log(betamax), C_1 <- (beta^(1-q)-1)/(betamax^(1-q)-1))
  ifelse(q==1, U_1 <- log(beta)/log(betamax), U_1 <- (beta^(q-1)-1)/(betamax^(q-1)-1))
  V_1 <- (beta-1)/(betamax-1)
  S_1 <- (beta^-1-1)/(betamax^-1-1)
  diff <- c("1-CqN" = C_1, "1-UqN" = U_1, "1-VqN" = V_1, "1-SqN" = S_1)
  names(diff) <- c(sapply(c("1-CqN", "1-UqN", "1-VqN", "1-SqN"), function(y) paste0(y,"(",apply(all.pair[,1:(H-1)], 2,paste, collapse = "|"), ")")))
  names(beta) <- paste0("qPD_Beta ",seq_len(H-1))
  names(betamax) <- paste0("qPD_Beta_max ",seq_len(H-1))
  out <-  c(Ialpha.relative ,Ibeta.relative, differential.relative, alpha.relative, beta, betamax, diff)
  #out <- data.frame(out)
  return(out)
}

phy.H.abs <- function(dat, mat, tmp, q, rtreephy, wij, type){
  n <- sum(dat)
  nsite <- ncol(mat)
  H <- nrow(mat)
  index <- lapply(rev(seq_len(H)), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x))))
  
  #rownames(aa) <- rownames(tmp$Gamma)
  
  # alpha.relative <- numeric(H)
  # alpha.absolute <- numeric(H)
  #wij <- colSums(dat)/n
  w <- lapply(index, function(x) sapply(x, function(y) sum(wij[y])))
  M <- lapply(index, function(x) sapply(x,length))
  ga <- tmp$Gamma[,1];gB <- tmp$Gamma[,2]
  names(gB) <- rownames(tmp$Gamma)
  gp=ga/n;TT=sum(gp*gB);
  aa <- sapply(tmp$Alpha, function(X){
    abun <- X[,1]
    names(abun) = rownames(X)
    return(abun)
  })
  B <-length(ga)
  pij <- sapply(1:nsite, function(x) {
    abun <- aa[,x]
    n <- sum(dat[,x])
    return(abun/n)
  }) #pk|im
  weight.p <- pij*t(replicate(B,wij))
  p <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
  #x <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
  nlist <- lapply(index, function(h) sapply(h, function(x) sum(dat[ ,x])))
  xtmp <- lapply(index, function(h) sapply(h, function(x) rowSums(cbind(aa[ ,x]))))
  if(type == "mle"){
    qD <- sapply(1:H, function(i) apply(p[[i]], 2, get("mle.phy.q"), LL = gB, TT = TT, q)) ##\sum{p^q} or -\sum{p*log(p)})
  } else if(type == "est"){
    qD <- sapply(1:H, function(i) sapply(seq_along(M[[i]]), function(j){
      n <- nlist[[i]][j]
      get("est.phy.q")(xtmp[[i]][,j], LL = gB, TT = TT, n = n, q = q, rtreephy)
    }))
  }
  alpha.relative <- sapply(1:H, function(i){
    qD_est <- qD[[i]]
    if(q!=1){
      (sum((nlist[[i]]/n)^q*M[[i]]^(1-q)*qD_est)^(1-q))/nsite
      #(w[[i]]%*%qD_est)^(1/(1-q))
    } else{
      exp(w[[i]]%*%(qD_est/TT+log(TT)))
    }
  })
  #names(alpha.absolute) <- c(paste0("alpha",seq_len(H-1)),"gamma")
  names(alpha.relative) <- c(paste0("qPD_alpha",seq_len(H-1)),"qPD_gamma")
  all.pair <- cbind(sapply(1:(H-1), function(i) cbind(i,i+1)), c(1,H))
  Ialpha.relative <- sapply(alpha.relative ,function(x){
    if(q!=1){
      alpha.relative.H.j <- (TT-x^(1-q)*TT^q*(nsite^(1-q)))/(q-1)
    } else{
      alpha.relative.H.j <- (log(x)-log(TT))*TT
    }
    alpha.relative.H.j})
  names(Ialpha.relative) <- c(paste0("qI_alpha", seq_len(H-1)), "qI_gamma")
  Ibeta.relative <- apply(all.pair, 2, function(j){
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q*(nsite^(1-q)))/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)
  })
  names(Ibeta.relative) <- c(paste0("qI_Beta ", "level(",apply(all.pair, 2,paste, collapse = "|"), ")"))
  differential.relative <-  apply(all.pair, 2, function(j){
    #qD = get(paste0("qD.", Q))
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    M1 <- M[[j[[1]]]]
    M2 <- M[[j[[2]]]]
    num <- sapply(index[[j[1]]], function(a) seq_along(M2)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    Max <- ifelse(q!=1, sum((M1^(1-q)-M2[num]^(1-q))*(nlist[[j[[1]]]]/n)^q*TT^q*qD[[j[[1]]]])/(q-1),
                  TT*(sum(w.down*log(M2))-sum(w.up*log(M1))))
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q*(nsite^(1-q)))/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)/Max
    #sum(w.down*sapply(index[[j[2]]], function(g) sum(ratio[g]*(1-ratio[g]^(q-1))*group.p[g])))
  })
  names(differential.relative ) <- paste0("Delta ","level(",apply(all.pair, 2,paste, collapse = "|"), ")")
  beta <- apply(all.pair, 2, function(j) alpha.relative[j[2]]/alpha.relative[j[1]])
  beta.max <- apply(all.pair, 2, function(j){
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    M1 <- M[[j[[1]]]]
    M2 <- M[[j[[2]]]]
    num <- sapply(index[[j[1]]], function(a) seq_along(M2)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    ifelse(q!=1,
           (sum(M2[num]^(1-q)*((nlist[[j[[1]]]]/n))^q*qD[[j[[1]]]])/sum(M1^(1-q)*((nlist[[j[[1]]]]/n))^q*qD[[j[[1]]]]))^(1/(1-q)),
           exp(sum(w.down*log(M2)))/exp(sum(w.up*log(M1))))
  })
  beta <-  beta[1:(H-1)]
  beta.max <- beta.max[1:(H-1)]
  ifelse(q==1, C_1 <- log(beta)/log(beta.max), C_1 <- (beta^(1-q)-1)/(beta.max^(1-q)-1))
  ifelse(q==1, U_1 <- log(beta)/log(beta.max), U_1 <- (beta^(q-1)-1)/(beta.max^(q-1)-1))
  V_1 <- (beta-1)/(beta.max-1)
  S_1 <- (beta^-1-1)/(beta.max^-1-1)
  diff <- c("1-C_qN" = C_1, "1-U_qN" = U_1, "1-V_qN" = V_1, "1-S_qN" = S_1)
  names(diff) <- c(sapply(c("1-CqN", "1-UqN", "1-VqN", "1-SqN"), function(y) paste0(y,"(",apply(all.pair[,1:(H-1)], 2,paste, collapse = "|"), ")")))
  names(beta) <- paste0("qPD_Beta ",seq_len(H-1))
  names(beta.max) <- paste0("qPD_Beta_max ",seq_len(H-1))
  #names(diff) <- c(sapply(c("1-C", "1-U", "1-V", "1-S"), function(y) withMathJax(paste0("$$",y,"_{",apply(all.pair, 2,paste, collapse = "|"), "}$$"))))
  #names(diff) <- c(sapply(c("1-C", "1-U", "1-V", "1-S"), function(y) paste0("$",y,"{",apply(all.pair, 2,paste, collapse = "|"), "}$")))
  #names(out) <- quote(names(out))
  #names(out) <- c(paste0("alpha of level ",seq_len(H-1)),"gamma", paste0("beta ",beta.pair[1, ]), paste0("differentiation ","level",apply(all.pair, 2,paste, collapse = "|")))
  out <-  c(Ialpha.relative ,Ibeta.relative, differential.relative, alpha.relative, beta, beta.max, diff)
  return(out)
}

#Functional
hier.func.abs <- function(dat, mat, dis, q, FDtype, FDtau, type, nboot, datatype){
  
  deltak1a <- function(k, n, x, a, zp){
    
    summa <- function(kk){
      
      index <- which(a >= 1 & a <= (n - kk))
      
      sum(x[index]/zp[index] * a[index]/n * exp(lchoose((n-a[index]), kk) - lchoose((n-1), kk)))
      
    }
    
    
    sapply(k, summa)
    
  }
  
  z.plus <- function(data, dis, FDtau, w){
    
    dis <- as.matrix(dis)
    dis[which(dis>FDtau,arr.ind = T)] <- FDtau
    
    
    mm <- round(apply(data, 1, function(x) {as.vector((1-dis/FDtau) %*% x )}))
    
    if (length(w) == 1){
      
      res <- rowSums( mm * as.matrix(replicate(ncol(data), w)))
      
    } else {
      
      res <- rowSums( mm * t(replicate(ncol(data), w)))
      
    }
    
    return(res)
    
    
  }
  
  x.plus <- function(data, w){
    
    
    if (length(w) == 1){
      
      res <- rowSums( t(data) * as.matrix(replicate(ncol(data), w)))
      
    } else {
      
      res <- rowSums( t(data) * t(replicate(ncol(data), w)))
      
    }
    
    return(res)
    
    
  }
  
  Boots.beta = function(data, dij, datatype){
    dij <- as.matrix(dij)
    #data is pooled matrix
    if(datatype == "abundance"){
      N = ncol(data); n = colSums(data)
      pool=rowSums(data) ; OBS=length(pool[pool>0])
      dij <- dij[pool>0, pool>0]
      data=data[pool>0,]
      obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]))
      F1=sum(pool==1);F2=sum(pool==2)
      F0=ceiling((sum(n)-1)/sum(n)*ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2)))
      f1=sapply(1:N,function(k) sum(data[,k]==1))
      f2=sapply(1:N,function(k) sum(data[,k]==2))
      C=sapply(1:N, function(k) {
        ifelse(f2[k]>0, 1-f1[k]*(n[k]-1)*f1[k]/n[k]/((n[k]-1)*f1[k]+2*f2[k]),
               1-f1[k]*(n[k]-1)*(f1[k]-1)/n[k]/((n[k]-1)*(f1[k]-1)+2))
      })
      # C=sapply(1:N,function(k) 1-f1[k]/n[k])
      f0=ceiling(sapply(1:N,function(k) ((n[k]-1)/n[k])*ifelse(f2[k]==0,f1[k]*(f1[k]-1)/2,f1[k]^2/(2*f2[k]))))
      r.data=sapply(1:N,function(k) data[,k]/n[k])
      W=sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k]))
      if(F0>0){boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=F0))
      }else{boots.pop=r.data}
      for(i in 1:N)
      {
        if(f0[i]>0)
        {
          f0[i]=ifelse(f0[i]+obs[i]>OBS+F0, OBS+F0-obs[i],f0[i])
          boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^n[i])   #
          I=which(boots.pop[,i]==0);II=sample(I,f0[i])
          boots.pop[II,i]=rep((1-C[i])/f0[i],f0[i])
        }
        #return(boots.pop)
      }
    }
    if(datatype == "incidence"){
      dat = data ; t = data[1, ] %>% as.numeric()
      data = data[-1, ]
      Tt=sum(t) ; N=ncol(data)
      pool=rowSums(data);OBS=length(pool[pool>0]);
      dij <- dij[pool>0, pool>0]
      data=data[pool>0,];
      obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]));
      Q1=sum(pool==1);Q2=sum(pool==2);
      Q0=ceiling(((Tt-1)/Tt)*ifelse(Q2==0,Q1*(Q1-1)/2,Q1^2/(2*Q2)));
      q1=sapply(1:N,function(k) sum(data[,k]==1));
      q2=sapply(1:N,function(k) sum(data[,k]==2));
      P1=sapply(1:N,function(k) ifelse(q1[k]+q2[k]==0,0,2*q2[k]/((t[k]-1)*q1[k]+2*q2[k])));##strange
      q0=ceiling(sapply(1:N,function(k) ((t[k]-1)/t[k])*ifelse(q2[k]==0,q1[k]*(q1[k]-1)/2,q1[k]^2/(2*q2[k]))));
      r.data=sapply(1:N,function(k) data[,k]/t[k]);
      W=sapply(1:N,function(k) (1-P1[k])*(q1[k]/t[k])/sum(r.data[,k]*(1-r.data[,k])^t[k]));##strange
      if(Q0>0){ boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=Q0))
      }else{boots.pop=r.data}
      for(i in 1:N)
      {
        if(q0[i]>0)
        {
          q0[i]=ifelse(q0[i]+obs[i]>OBS+Q0, OBS+Q0-obs[i],q0[i])
          boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^t[i])
          I=which(boots.pop[,i]==0);II=sample(I,q0[i])
          boots.pop[II,i]=rep((q1[i]/t[i])/q0[i],q0[i])
          
        }
      }
      #return(boots.pop)
    }
    
    X = rowSums(data)
    F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
    F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
    if (datatype=="abundance") {
      n_plus <- sum(n)
      F.0hat <- ifelse(F.2 > 0, ((n_plus-1)/n_plus) * (F.1^2/(2 * F.2)), ((n_plus-1)/n_plus)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((n_plus-2)* (n_plus-3)* (F11^2)/(4* n_plus* (n_plus-1)* F22)), ((n_plus-2)* (n_plus-3)* (F11*(F11-0.01))/(4 *n_plus * (n_plus-1))) )
      f0=F0
    } else if (datatype=="incidence") {
      T_plus <-sum(t)
      F.0hat <- ifelse(F.2 > 0, ((T_plus-1)/T_plus) * (F.1^2/(2 * F.2)), ((T_plus-1)/T_plus)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((T_plus-1)^2 * (F11^2)/(4* T_plus* T_plus* F22)), ((T_plus-1)* (T_plus-1)* (F11*(F11-0.01))/(4 *T_plus * T_plus)) )
      f0=Q0
    }
    
    if (f0==0) {
      d <- dij
    } else if (f0==1) {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0), length(X)*f0) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0, byrow = T)
      d00 = matrix(0, f0, f0)
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    } else {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0), length(X)*f0) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0, byrow = T)
      f0.num = (f0 * (f0-1) )/2
      random_d00 = as.vector(rmultinom(1, 1000, rep(1/f0.num, f0.num) ) )/1000
      d00 = matrix(0, f0, f0)
      d00[upper.tri(d00)] = (F00hat/2)*random_d00
      d00 <- pmax(d00, t(d00))###signmatrix
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    }
    
    
    return(list(pop=boots.pop, dij=d))
  }
  
  Dissimilarity_measure <- function(D.beta, D.Beta_max, q, pair2){
    if(q!=1){
      sorensen <- (D.beta^(1-q)-1)/(D.Beta_max^(1-q)-1)
      Jaccard <- (D.beta^(q-1)-1)/(D.Beta_max^(q-1)-1)
      h <- (D.beta^(-1)-1)/(D.Beta_max^(-1)-1)
      turnover_rate <- (D.beta-1)/(D.Beta_max-1)
    } else {
      sorensen <- log(D.beta)/log(D.Beta_max)
      Jaccard <- sorensen
      h <- (D.beta^(-1)-1)/(D.Beta_max^(-1)-1)
      turnover_rate <- (D.beta-1)/(D.Beta_max-1)
    } 
    out <- c(sorensen, Jaccard, h, turnover_rate)
    # names(out) <- c(paste0("1-CqN.", seq_along(D.beta)), paste0("1-UqN.", seq_along(D.beta)), paste0("1-SqN.", seq_along(D.beta)), paste0("1-VqN.", seq_along(D.beta)))
    names(out) <- c(paste0("1-CqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-UqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-SqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-VqN(", apply(pair2, 2, paste, collapse = '|'), ')'))
    out
  }
  
  max_est <- function(data, dij, FDtau, q, zz, all_data, de){
    
    dij = as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] = FDtau
    a = as.vector((1-dij/FDtau) %*% data )
    n = sum(data)
    a = round(a)
    data = data[a!=0]
    all_data = all_data[a != 0]
    zz = zz[a != 0]
    a = a[a!=0]
    
    data = data[zz!=0]
    all_data = all_data[zz != 0]
    a = a[zz!=0]
    zz = zz[zz != 0]
    
    f1.star = sum((data/a)[a==1]) ; f2.star = sum((data/a)[a==2])
    
    f1 = sum(a==1); f2 = sum(a==2)
    
    if(f2 > 0 & f1 > 0){
      a1 = 2*f2/((n-1)*f1+2*f2)
    }else if(f2 == 0 & f1 > 0 ){
      a1 = 2/((n-1)*(f1-1)+2)
    }else{
      a1 = 1
    }
    
    r = 0:(n-1)
    A = sum(sapply(r, function(k) choose(q-1,k) * (-1)^k) * de )
    
    if(f1.star!=0){
      
      B = ifelse(a1==1 , 0, (f1.star/n)*(1-a1)^(1-n)*(a1^(q-1)-sum(choose(q-1,r)*(a1-1)^r)))
      out = A + B
      
    } else {
      
      out = A
      
    }
    
    out
  }
  
  q1_est_FD <- function(data, dij, FDtau, q, zz, all_data, de){
    
    dij = as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] = FDtau
    a = as.vector((1-dij/FDtau) %*% data )
    n = sum(data)
    
    
    a = round(a)
    data = data[a!=0]
    all_data = all_data[a != 0]
    zz = zz[a != 0]
    a = a[a!=0]
    
    data = data[zz!=0]
    all_data = all_data[zz != 0]
    a = a[zz!=0]
    zz = zz[zz != 0]
    
    f1.star = sum((data/a)[a==1]) ; f2.star = sum((data/a)[a==2])
    
    f1 = sum(a==1); f2 = sum(a==2)
    
    
    if(f2 > 0 & f1 > 0){
      a.star = 2*f2/((n-1)*f1+2*f2)
    }else if(f2 == 0 & f1 > 0 ){
      a.star = 2/((n-1)*(f1-1)+2)
    }else{
      a.star = 1
    }
    
    
    A = sum(de[-1]/1:(n-1))
    
    B = ifelse(f1.star==0 | a.star==1, 0, (f1.star/n)*(1-a.star)^(1-n)*(-log(a.star) - sum(sapply(1:(n-1), function(x) (1-a.star)^x/x))) )
    
    out = A + B
    
    out
  }
  
  FD_est_abs <- function(data, dij, FDtau, q, Mi, wij, zz, all_data){ ## Diversity
    
    
    DE_maker2 <- function(data, dij, FDtau, zz, all_data){
      
      dij = as.matrix(dij)
      dij[which(dij>FDtau,arr.ind = T)] = FDtau
      a = as.vector((1-dij/FDtau) %*% data )
      n = sum(data)
      a = round(a)
      data = data[a!=0]
      all_data = all_data[a != 0]
      zz = zz[a != 0]
      a = a[a!=0]
      
      data = data[zz!=0]
      all_data = all_data[zz != 0]
      a = a[zz!=0]
      zz = zz[zz != 0]
      r = 0:(n-1)
      return(deltak1a(r,n,all_data,a, zz))
      
    }
    
    
    ind_de <- list()
    
    for(kk in 1 : nrow(data) ){
      
      
      ind_de[[kk]] <- DE_maker2(data[kk,], dij, FDtau, zz, all_data)
      
    }
    
    Sub1 <- function(q, de2){
      
      dij <- as.matrix(dij)
      dij[which(dij>FDtau,arr.ind = T)] <- FDtau
      
      
      if(q == 1){
        
        zzz <- zz
        all_datt <- all_data
        
        A =  sum( unlist(wij) * sapply(1:nrow(data), function(i){
          
          q1_est_FD(data[i, ], dij = dij, FDtau = FDtau, q = q, zz = zzz, all_data = all_datt, de = de2[[i]])
          
          
        }) )
        
        C1 <- sapply(1:nrow(data), function(i) wij[[i]]*log(Mi[[i]]))
        
        
        C2 <- sapply(1:nrow(data), function(i){
          
          max_est(data[i, ], dij = dij, FDtau = FDtau, q = q, zz = zzz, all_data = all_datt, de = de2[[i]])
          
          
        })
        
        
        B <- sum( (unlist(wij) * ( log(unlist(Mi)) - log(unlist(wij)) )) * C2 )
        
        out <- A + B
        
        res <- (exp(out) / sum(unlist(Mi)))
        res <- list('result' = res, 'max' = sum(C1 * C2))
        
      } else {
        
        # tmp <- apply(data, 1, max_est, dij = dij, FDtau = FDtau, q = q, zz = zz, all_data = all_data)
        
        tmp <- sapply(1:nrow(data), function(i){
          
          max_est(data[i, ], dij = dij, FDtau = FDtau, q = q, zz = zz, all_data = all_data, de = de2[[i]])
          
          
        })
        
        # tmp2 <- rowSums(matrix(data, nrow(data)))/n
        tmp2 <- unlist(wij)
        res <- sum(unlist(Mi)^(1-q)*tmp2^q*tmp)^(1/(1-q)) / sum(unlist(Mi))
        res <- list('result' = res, 'max' = tmp)
        
      }
      
      return(res)
      
      
    }
    
    lapply(q, Sub1, ind_de)
    
  }
  
  max_mle <- function(data, dij, FDtau, q, vv){
    
    dij <- as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] <- FDtau
    a <- as.vector((1 - dij/FDtau) %*% data )  
    data <- data[a!=0]
    v <- c(vv)[a != 0]
    a <- a[a!=0]
    #v <- data/a
    nplus <- sum(data)
    
    return( sum(v*(a/nplus)^q) )
  }
  
  FD_mle_abs <- function(data1, data2, dij, FDtau, q, wij, Mi, vv) {
    
    
    dij <- as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] <- FDtau
    n = sum(data2)
    
    
    Sub2 <- function(q){
      
      if(q == 1){
        
        vvv <- vv
        
        a <- as.vector(((1-dij/FDtau) %*% t(data2)))
        
        vv <- rep(vv, nrow(data2))
        v <- c(vv)[a != 0]
        a <- a[ a!= 0]
        
        
        A <- sum(-v*a/n*log(a/n))
        tmp2 <- rowSums(matrix(data2, nrow(data2)))/n
        # C1 <- sapply(1:nrow(data2), function(i) wij[[i]]*log(Mi[[i]]))
        C1 <- sapply(1:nrow(data2), function(i) tmp2[i]*log(Mi[[i]]))
        C2 <- apply(data1, 1, max_mle, dij = dij, FDtau = FDtau, q = q, vv = vvv)
        
        out = A + sum(C1 * C2)
        
        res <- exp(out) / sum(unlist(Mi))
        res <- list('result' = res, 'max' = sum(C1 * C2))
        
      } else {
        
        tmp <- apply(data1, 1, max_mle, dij = dij, FDtau = FDtau, q = q, vv = vv)
        # tmp2 <- rowSums(data)/n
        tmp2 <- rowSums(matrix(data2, nrow(data2)))/n
        # tmp2 <- unlist(wij)
        res <- sum(unlist(Mi)^(1-q)*tmp2^q*tmp)^(1/(1-q)) / sum(unlist(Mi))
        res <- list('result' = res, 'max' = tmp)
        
      }
      
      
    }
    
    lapply(q, Sub2)
    
  }
  
  Sub <- function(dat, mat, dis, q, FDtau, type, SD = F){
    
    vk <- function(pk_pool,dis, FDtau){
      dis[dis > FDtau] <- FDtau
      d <- dis/FDtau
      pk_pool/(pk_pool%*%(1-d))
      
    }
    
    
    dat <- as.matrix(dat)
    dis <- as.matrix(dis)
    population <- ncol(mat)
    H <- nrow(mat)
    mat[H,] <- paste0(mat[H,], seq(1,ncol(mat)))
    
    n_index <- which(rowSums(dat) > 0 )
    dat <- dat[n_index, ]
    dis <- dis[n_index, n_index]
    
    dat <- as.matrix(dat)
    dis <- as.matrix(dis)
    
    n <- sum(dat)
    s <- nrow(dat)
    
    M <- vector("list", H)
    alpha.relative <- list()
    
    check <- apply(dat, 2, function(x) sum(x))
    
    if( sum(check == 0) > 0 ){
      
      index <- which( check == 0)
      dat <- dat[, -index]
      mat <- mat[, -index]
      
      # if( !is.null(wij) ){
      #   
      #   wij <- wij[-index]
      #   
      # }
      
    }
    
    
    pij <- apply(dat, 2, function(x) x/sum(x))
    
    
    if(FDtau == "dmin"){
      
      mint <- min(dis[dis>0])
      FDtau <- mint
      
    }
    
    if(FDtau == "dmean"){
      
      tmp <- rowSums(pij)/sum(pij)
      meant <- c(tmp%*%dis%*%tmp)
      FDtau <- meant
      
    }
    
    
    if(FDtau == 'dmax'){
      
      maxt <- max(dis)
      FDtau <- maxt
      
    }
    
    
    mat <- apply(mat, 2, as.character)
    index <- rev(lapply(seq_len(H), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x)))))
    
    # if(wij[1] == 'size' ) wij <- apply(dat, 2, function(x) sum(x)/n)
    # if(wij[1] == 'equal') wij <- rep(1/ncol(dat), ncol(dat))
    
    wij <- apply(dat, 2, function(x) sum(x)/n)
    wij <- wij/sum(wij)
    
    
    w <- lapply(index, function(x) lapply(x, function(y) sum(wij[y])))
    pij <- apply(dat, 2, function(x) x/sum(x))
    
    weight.p <- pij*t(replicate(s,wij))
    
    p <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
    xx <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
    x <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x]))))
    pk_plpl <- do.call(rbind, p[[H]])
    #v_rel <- vk(pk_plpl, dis, FDtau)
    v_abs <- vk(rowSums(dat), dis, FDtau)
    M <- lapply(index, function(x) lapply(x,length))
    
    est.max <- list()
    mle.max <- list()
    
    for(i in 1:H){
      if(type == 'est'){
        group <- do.call(rbind, xx[[i]])
        zzl = round(z.plus(group, dis, FDtau, rep(1, length(unlist(w[[i]])))))
        xxx <- round(x.plus(group, rep(1, length(unlist(w[[i]])))))
        est.abs <- FD_est_abs(data = group, dij = dis, FDtau = FDtau, q = q, Mi = M[[i]], wij = w[[i]], zz = zzl, all_data = xxx)
        est.max[[i]] <- lapply(est.abs, function(x) x$max)
        alpha.relative[[i]] <- unlist(lapply(est.abs, function(x) x$result))
        
      } else {
        group1 <- do.call(rbind, p[[i]])
        group2 <- do.call(rbind, xx[[i]])
        mle.abs <- FD_mle_abs(data1 = group1, data2 = group2, dij = dis, FDtau = FDtau, q = q, wij = w[[i]], Mi = M[[i]], vv = v_abs)
        mle.max[[i]] <- lapply(mle.abs, function(x) x$max)
        alpha.relative[[i]] <- unlist(lapply(mle.abs, function(x) x$result))
      }
    }
    
    alpha.relative <- do.call(rbind, alpha.relative)
    
    
    H.alpha <- sapply(1:length(q), function(i){
      
      if(q[i] == 1){
        
        log(alpha.relative[,i])
        
      } else {
        
        (1-(population^(q[i]-1)*(alpha.relative[,i]*population)^(1-q[i])))/(q[i]-1)
        
      }
      
    })
    
    
    # if(q!=1) H.alpha <- (1-(population^(q-1)*(alpha.relative*population)^(1-q)))/(q-1) else H.alpha <- log(alpha.relative)
    
    rownames(alpha.relative) <- c(paste0("qFD_alpha",seq_len(H-1)),"qFD_gamma")
    rownames(H.alpha) <- c(paste0("qQ_alpha",seq_len(H-1)),"qQ_gamma")
    
    D.beta <- alpha.relative[-1,]/alpha.relative[-nrow(alpha.relative), ]
    D.beta <- matrix(D.beta, ncol = length(q))
    
    if( SD == F) D.beta[D.beta < 1] <- 1
    
    all.pair <- combn(seq_len(H),2)
    
    if(H == 2){
      pair2 = pair <- all.pair
    } else {
      pair <- cbind(all.pair[,diff(all.pair) == 1], c(1,H))
      pair2 <- pair[ ,-ncol(pair)]
    } 
    
    
    differential.relative <- lapply(1:ncol(pair), function(l){
      # up <- ifelse(q!=1, -diff(alpha.relative[j]^(1-q))*population^(1-q), diff(log(alpha.relative[j]))) #*#
      
      j <- pair[,l] 
      
      up <- c(sapply(1:length(q), function(i){
        
        if( q[i] == 1){
          
          diff(log(alpha.relative[j, i]))
          
        } else {
          
          -diff(alpha.relative[j, i]^(1-q[i]))*population^(1-q[i])
          
        }
        
      }))
      names(up) <- NULL
      
      M1 <- unlist(M[[j[1]]])
      M2 <- unlist(M[[j[2]]])
      num <- sapply(index[[j[1]]], function(a) seq_along(M2)[sapply(index[[j[2]]], function(b) all(a%in%b))])
      
      ss <- lapply(1:length(q), function(i){
        
        if( q[i] == 1 ) {
          
          if( type == 'est'){
            
            Max <- est.max[[j[2]]][[i]] - est.max[[j[1]]][[i]]
            Gamma.Max <- exp(Max)
            
            
          } else {
            
            Max <- mle.max[[j[2]]][[i]] - mle.max[[j[1]]][[i]]
            Gamma.Max <- exp(Max)
            
            
          }
          
          
        }else{
          if(type == 'est'){
            
            tmp <- est.max[[j[1]]][[i]]
            tmp2 <- unlist(w[[j[1]]])^q[i]
            tmp3 <- unlist(sapply(M2, function(x) rep(x,x)))[1:length(M1)]
            Max <- sum((M1^(1-q[i])-tmp3^(1-q[i]))*tmp2*tmp)
            Gamma.Max <- sum(tmp3^(1-q[i])*tmp2*tmp)^(1/(1-q[i])) / population
            
          } else {
            
            tmp <- mle.max[[j[1]]][[i]]
            tmp2 <- unlist(w[[j[1]]])^q[i]
            tmp3 <- unlist(sapply(M2, function(x) rep(x,x)))[1:length(M1)]
            Max <- sum((M1^(1-q[i])-tmp3^(1-q[i]))*tmp2*tmp)
            Gamma.Max <- sum(tmp3^(1-q[i])*tmp2*tmp)^(1/(1-q[i])) / population
            
            
          }
        }
        
        
        c(Max = Max, Gamma.Max = Gamma.Max)
        
      })
      
      return(rbind(up, do.call(cbind, ss)))
      
    })
    
    
    
    
    
    
    
    
    differential <- do.call(rbind, lapply(differential.relative, function(x) x[1,]/x[2,]))
    if( SD == F) differential[differential < 0] <- 0
    rownames(differential) <- paste0("Delta ","level(",apply(pair, 2,paste, collapse = "|"), ')')
    
    
    D.Beta_max <- sapply(1:length(q), function(i){
      
      if(q[i] == 1){
        
        do.call(rbind, lapply(differential.relative, function(x) x[3,i]))[-nrow(alpha.relative), ]
        
        
      } else {
        
        do.call(rbind, lapply(differential.relative, function(x) x[3,i]))[-nrow(alpha.relative), ] / alpha.relative[-nrow(alpha.relative), i]
        
      }
      
    })
    
    
    D.Beta_max <- matrix(D.Beta_max, ncol = length(q))
    
    
    rownames(D.Beta_max) <- paste0("qFD_Beta_max ", seq_len(nrow(D.beta)))
    
    D.Diss <- sapply(seq_along(q), function(i) Dissimilarity_measure(D.beta[,i], D.Beta_max[,i], q[i], pair2) )
    
    if( SD == F) D.Diss <- abs(D.Diss)
    
    if(H == 2){
      
      beta <- apply(H.alpha, 2, diff)
      
    } else {
      
      beta <- rbind(apply(H.alpha, 2, diff), H.alpha[nrow(H.alpha), ] - H.alpha[1, ])
    }
    
    
    if( SD == F) beta[beta < 0] <- 0
    beta <- matrix(beta, ncol = length(q))
    rownames(beta) <- paste0("qQ_Beta ","level(",apply(pair, 2,paste, collapse = "|"), ')')
    rownames(D.beta) <- paste0("qFD_Beta ", seq_len(nrow(D.beta)))
    
    # result <- data.frame(relative=c(alpha.relative, differential))
    TAU <- c('FDtau' = FDtau)
    
    result <- rbind(H.alpha, beta, differential, alpha.relative, D.beta, D.Beta_max, D.Diss)
    
    # result <- data.frame(relative = c(H.alpha, beta, differential, alpha.relative, D.beta, D.Beta_max, D.Diss, TAU))
    out <- result
    out
    
  }
  
  if(datatype == 'incidence'){
    
    su <- as.numeric(dat[1, ])
    dat <- dat[-1, ]
    
  }
  
  res <- Sub(dat, mat, dis, q, FDtau, type)
  
  AUC_func <- function(dat, mat, dis, q, type){
    tau1 <- seq(0.00001, 1, length.out = 30)
    ans <- sapply(q, function(s) {
      qFun = sapply(tau1, function(ta) Sub(dat, mat, dis, s, ta, type, SD = T))
      AUC <- apply(qFun, 1, function(i) {
        LA <- sum(i[seq_along(tau1[-1])]*diff(tau1))
        RA <- sum(i[-1]*diff(tau1))
        mean(c(LA, RA))
      })
      AUC
    })
    ans
  }
  if(FDtype == "AUC"){
    res.auc = AUC_func(dat, mat, dis, q, type)
    rownames(res.auc) = rownames(res)
  }
  if( nboot > 1 ){
    
    if( datatype == 'abundance'){
      
      n <- sum(dat)
      w <- colSums(dat)/n
      P <- Boots.beta(dat, dis, datatype)$pop
      P <- sapply(1:length(w), function(i) P[,i] * w[i])
      P <- P/sum(P)
      dij <- Boots.beta(dat, dis, datatype)$dij
      X <- replicate(nboot, matrix(rmultinom(1, n, c(P)), ncol = length(w)))
      
    }
    
    if( datatype == 'incidence' ){
      
      P <- Boots.beta(rbind(su, dat), dis, datatype)$pop
      dij <- Boots.beta(rbind(su, dat), dis, datatype)$dij
      X <- replicate(nboot, sapply(seq_along(su), function(i) rbinom(nrow(P),c(su[i]), P[,i])))
      
    }
    
    
    if( sum( dij > max(dis) ) > 0 ){
      
      
      dij[dij > max(dis)] <- (dij[dij > max(dis)]/max(dij[dij > max(dis)])) * max(dis)
      
      
    }
    
    dij <- replace(dij, dij == 0, 10^(-10))
    diag(dij) <- 0
    
    
    boot.est <- lapply(1:nboot, function(i){
      
      Sub(X[,,i], mat, dij, q, FDtau, type, SD = T)
      
    })
    
    est.se <- apply(array(as.numeric(unlist(boot.est)), dim = c(dim(res), nboot)), c(1,2), function(x) ifelse(is.na(sd(x[!is.na(x)])), 0, sd(x[!is.na(x)])))
    rownames(est.se) <- rownames(res)
    
    if(FDtype == "AUC"){
      boot.est.auc <- lapply(1:nboot, function(i){
        
        AUC_func(X[,,i], mat, dij, q, type)
        
      })
      est.se.auc <- apply(array(as.numeric(unlist(boot.est.auc)), dim = c(dim(res.auc), nboot)), c(1,2), function(x) ifelse(is.na(sd(x[!is.na(x)])), 0, sd(x[!is.na(x)])))
      rownames(est.se.auc) <- rownames(res)
    }
    
    
    
  } else {
    
    est.se <- matrix(NA, dim(res)[1], dim(res)[2])
    rownames(est.se) <- rownames(res)
    
    if(FDtype == "AUC"){
      est.se.auc <- matrix(NA, dim(res.auc)[1], dim(res.auc)[2])
      rownames(est.se.auc) <- rownames(res)
    }
  }
  
  if(FDtype == "AUC"){
    return( list('result' = res.auc, 'S.E.' = est.se.auc) )
  }else{
    return( list('result' = res, 'S.E.' = est.se) )
  }
  
}

hier.func.rel <- function(dat, mat, dis, q, FDtype, FDtau, weight, type, nboot, datatype){
  
  deltak1a <- function(k, n, x, a, zp){
    
    summa <- function(kk){
      
      index <- which(a >= 1 & a <= (n - kk))
      
      sum(x[index]/zp[index] * a[index]/n * exp(lchoose((n-a[index]), kk) - lchoose((n-1), kk)))
      
    }
    
    
    sapply(k, summa)
    
  }
  
  z.plus <- function(data, dis, FDtau, w){
    
    dis <- as.matrix(dis)
    dis[which(dis>FDtau,arr.ind = T)] <- FDtau
    
    
    mm <- round(apply(data, 1, function(x) {as.vector((1-dis/FDtau) %*% x )}))
    
    if (length(w) == 1){
      
      res <- rowSums( mm * as.matrix(replicate(ncol(data), w)))
      
    } else {
      
      res <- rowSums( mm * t(replicate(ncol(data), w)))
      
    }
    
    return(res)
    
    
  }
  
  x.plus <- function(data, w){
    
    
    if (length(w) == 1){
      
      res <- rowSums( t(data) * as.matrix(replicate(ncol(data), w)))
      
    } else {
      
      res <- rowSums( t(data) * t(replicate(ncol(data), w)))
      
    }
    
    return(res)
    
    
  }
  
  Dissimilarity_measure <- function(D.beta, D.Beta_max, q, pair2){
    if(q!=1){
      sorensen <- (D.beta^(1-q)-1)/(D.Beta_max^(1-q)-1)
      Jaccard <- (D.beta^(q-1)-1)/(D.Beta_max^(q-1)-1)
      h <- (D.beta^(-1)-1)/(D.Beta_max^(-1)-1)
      turnover_rate <- (D.beta-1)/(D.Beta_max-1)
    } else {
      sorensen <- log(D.beta)/log(D.Beta_max)
      Jaccard <- sorensen
      h <- (D.beta^(-1)-1)/(D.Beta_max^(-1)-1)
      turnover_rate <- (D.beta-1)/(D.Beta_max-1)
    } 
    out <- c(sorensen, Jaccard, h, turnover_rate)
    # names(out) <- c(paste0("1-CqN.", seq_along(D.beta)), paste0("1-UqN.", seq_along(D.beta)), paste0("1-SqN.", seq_along(D.beta)), paste0("1-VqN.", seq_along(D.beta)))
    names(out) <- c(paste0("1-CqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-UqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-SqN(", apply(pair2, 2, paste, collapse = '|'), ')'), paste0("1-VqN(", apply(pair2, 2, paste, collapse = '|'), ')'))
    out
  }
  
  FD_est_part <- function(data, dij, FDtau, q, zz, all_data, de){
    
    dij = as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] = FDtau
    a = as.vector((1-dij/FDtau) %*% data )
    n = sum(data)
    a = round(a)
    data = data[a!=0]
    all_data = all_data[a != 0]
    zz = zz[a != 0]
    a = a[a!=0]
    
    data = data[zz!=0]
    all_data = all_data[zz != 0]
    a = a[zz!=0]
    zz = zz[zz != 0]
    
    f1.star = sum((data/a)[a==1]) ; f2.star = sum((data/a)[a==2])
    
    f1 = sum(a==1); f2 = sum(a==2)
    
    if(f2 > 0 & f1 > 0){
      a1 = 2*f2/((n-1)*f1+2*f2)
    }else if(f2 == 0 & f1 > 0 ){
      a1 = 2/((n-1)*(f1-1)+2)
    }else{
      a1 = 1
    }
    
    if( q == 1 ){
      
      A = sum(de[-1]/1:(n-1))
      if(f1.star!=0){
        
        B = ifelse(f1.star==0|a1==1, 0, (f1.star/n)*(1-a1)^(1-n)*(-log(a1) - sum(sapply(1:(n-1), function(x) (1-a1)^x/x))) ) 
        out = A+B
      }else{
        
        out = A
      }
      
      
    } else {
      
      r = 0:(n-1)
      A = sum(sapply(r, function(k) choose(q-1,k) * (-1)^k) * de )
      
      if(f1.star!=0){
        
        B = ifelse(a1==1 , 0, (f1.star/n)*(1-a1)^(1-n)*(a1^(q-1)-sum(choose(q-1,r)*(a1-1)^r)))
        out = A + B
        
      } else {
        
        out = A
        
      }
      
    }
    
    out
  }
  
  max_est <- function(data, dij, FDtau, q, zz, all_data, de){
    
    dij = as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] = FDtau
    a = as.vector((1-dij/FDtau) %*% data )
    n = sum(data)
    a = round(a)
    data = data[a!=0]
    all_data = all_data[a != 0]
    zz = zz[a != 0]
    a = a[a!=0]
    
    data = data[zz!=0]
    all_data = all_data[zz != 0]
    a = a[zz!=0]
    zz = zz[zz != 0]
    
    f1.star = sum((data/a)[a==1]) ; f2.star = sum((data/a)[a==2])
    
    f1 = sum(a==1); f2 = sum(a==2)
    
    if(f2 > 0 & f1 > 0){
      a1 = 2*f2/((n-1)*f1+2*f2)
    }else if(f2 == 0 & f1 > 0 ){
      a1 = 2/((n-1)*(f1-1)+2)
    }else{
      a1 = 1
    }
    
    r = 0:(n-1)
    A = sum(sapply(r, function(k) choose(q-1,k) * (-1)^k) * de )
    
    if(f1.star!=0){
      
      B = ifelse(a1==1 , 0, (f1.star/n)*(1-a1)^(1-n)*(a1^(q-1)-sum(choose(q-1,r)*(a1-1)^r)))
      out = A + B
      
    } else {
      
      out = A
      
    }
    
    out
  }
  
  FD_est_rel <- function(data, dij, FDtau, q, wij, zz, all_data){ ## Diversity
    
    
    
    DE_maker2 <- function(data, dij, FDtau, zz, all_data){
      
      dij = as.matrix(dij)
      dij[which(dij>FDtau,arr.ind = T)] = FDtau
      a = as.vector((1-dij/FDtau) %*% data )
      n = sum(data)
      a = round(a)
      data = data[a!=0]
      all_data = all_data[a != 0]
      zz = zz[a != 0]
      a = a[a!=0]
      
      data = data[zz!=0]
      all_data = all_data[zz != 0]
      a = a[zz!=0]
      zz = zz[zz != 0]
      r = 0:(n-1)
      return(deltak1a(r,n,all_data,a, zz))
      
    }
    
    
    ind_de <- list()
    
    for(kk in 1 : nrow(data) ){
      
      
      ind_de[[kk]] <- DE_maker2(data[kk,], dij, FDtau, zz, all_data)
      
    }
    
    
    Sub1 <- function(q, de2){
      
      
      tmp <- sapply(1:nrow(data), function(i){
        
        FD_est_part(data[i, ], dij = dij, FDtau = FDtau, q = q, zz = zz, all_data = all_data, de = de2[[i]])
        
        
      }) 
      
      
      
      if( q == 1 ){
        
        tmp2 <- sapply(1:nrow(data), function(i){
          
          max_est(data[i, ], dij = dij, FDtau = FDtau, q = q, zz = zz, all_data = all_data, de = de2[[i]])
          
          
        }) 
        
        res <- exp(sum(unlist(wij) * tmp))
        res <- list('result' = res, 'max' = tmp2)
        
      } else {
        
        res <- sum(unlist(wij) * tmp)^(1/(1-q))
        res <- list('result' = res, 'max' = tmp)
        
      }
      
      return(res)
      
    }
    
    lapply(q, Sub1, ind_de)
    
  }
  
  Boots.beta = function(data, dij, datatype){
    dij <- as.matrix(dij)
    #data is pooled matrix
    if(datatype == "abundance"){
      N = ncol(data); n = colSums(data)
      pool=rowSums(data) ; OBS=length(pool[pool>0])
      dij <- dij[pool>0, pool>0]
      data=data[pool>0,]
      obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]))
      F1=sum(pool==1);F2=sum(pool==2)
      F0=ceiling((sum(n)-1)/sum(n)*ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2)))
      f1=sapply(1:N,function(k) sum(data[,k]==1))
      f2=sapply(1:N,function(k) sum(data[,k]==2))
      C=sapply(1:N, function(k) {
        ifelse(f2[k]>0, 1-f1[k]*(n[k]-1)*f1[k]/n[k]/((n[k]-1)*f1[k]+2*f2[k]),
               1-f1[k]*(n[k]-1)*(f1[k]-1)/n[k]/((n[k]-1)*(f1[k]-1)+2))
      })
      # C=sapply(1:N,function(k) 1-f1[k]/n[k])
      f0=ceiling(sapply(1:N,function(k) ((n[k]-1)/n[k])*ifelse(f2[k]==0,f1[k]*(f1[k]-1)/2,f1[k]^2/(2*f2[k]))))
      r.data=sapply(1:N,function(k) data[,k]/n[k])
      W=sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k]))
      if(F0>0){boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=F0))
      }else{boots.pop=r.data}
      for(i in 1:N)
      {
        if(f0[i]>0)
        {
          f0[i]=ifelse(f0[i]+obs[i]>OBS+F0, OBS+F0-obs[i],f0[i])
          boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^n[i])   #
          I=which(boots.pop[,i]==0);II=sample(I,f0[i])
          boots.pop[II,i]=rep((1-C[i])/f0[i],f0[i])
        }
        #return(boots.pop)
      }
    }
    if(datatype == "incidence"){
      dat = data ; t = data[1, ] %>% as.numeric()
      data = data[-1, ]
      Tt=sum(t) ; N=ncol(data)
      pool=rowSums(data);OBS=length(pool[pool>0]);
      dij <- dij[pool>0, pool>0]
      data=data[pool>0,];
      obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]));
      Q1=sum(pool==1);Q2=sum(pool==2);
      Q0=ceiling(((Tt-1)/Tt)*ifelse(Q2==0,Q1*(Q1-1)/2,Q1^2/(2*Q2)));
      q1=sapply(1:N,function(k) sum(data[,k]==1));
      q2=sapply(1:N,function(k) sum(data[,k]==2));
      P1=sapply(1:N,function(k) ifelse(q1[k]+q2[k]==0,0,2*q2[k]/((t[k]-1)*q1[k]+2*q2[k])));##strange
      q0=ceiling(sapply(1:N,function(k) ((t[k]-1)/t[k])*ifelse(q2[k]==0,q1[k]*(q1[k]-1)/2,q1[k]^2/(2*q2[k]))));
      r.data=sapply(1:N,function(k) data[,k]/t[k]);
      W=sapply(1:N,function(k) (1-P1[k])*(q1[k]/t[k])/sum(r.data[,k]*(1-r.data[,k])^t[k]));##strange
      if(Q0>0){ boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=Q0))
      }else{boots.pop=r.data}
      for(i in 1:N)
      {
        if(q0[i]>0)
        {
          q0[i]=ifelse(q0[i]+obs[i]>OBS+Q0, OBS+Q0-obs[i],q0[i])
          boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^t[i])
          I=which(boots.pop[,i]==0);II=sample(I,q0[i])
          boots.pop[II,i]=rep((q1[i]/t[i])/q0[i],q0[i])
          
        }
      }
      #return(boots.pop)
    }
    
    X = rowSums(data)
    F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
    F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
    if (datatype=="abundance") {
      n_plus <- sum(n)
      F.0hat <- ifelse(F.2 > 0, ((n_plus-1)/n_plus) * (F.1^2/(2 * F.2)), ((n_plus-1)/n_plus)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((n_plus-2)* (n_plus-3)* (F11^2)/(4* n_plus* (n_plus-1)* F22)), ((n_plus-2)* (n_plus-3)* (F11*(F11-0.01))/(4 *n_plus * (n_plus-1))) )
      f0=F0
    } else if (datatype=="incidence") {
      T_plus <-sum(t)
      F.0hat <- ifelse(F.2 > 0, ((T_plus-1)/T_plus) * (F.1^2/(2 * F.2)), ((T_plus-1)/T_plus)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((T_plus-1)^2 * (F11^2)/(4* T_plus* T_plus* F22)), ((T_plus-1)* (T_plus-1)* (F11*(F11-0.01))/(4 *T_plus * T_plus)) )
      f0=Q0
    }
    
    if (f0==0) {
      d <- dij
    } else if (f0==1) {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0), length(X)*f0) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0, byrow = T)
      d00 = matrix(0, f0, f0)
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    } else {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0), length(X)*f0) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0, byrow = T)
      f0.num = (f0 * (f0-1) )/2
      random_d00 = as.vector(rmultinom(1, 1000, rep(1/f0.num, f0.num) ) )/1000
      d00 = matrix(0, f0, f0)
      d00[upper.tri(d00)] = (F00hat/2)*random_d00
      d00 <- pmax(d00, t(d00))###signmatrix
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    }
    
    
    return(list(pop=boots.pop, dij=d))
  }
  
  max_mle <- function(data, dij, FDtau, q, vv){
    
    dij <- as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] <- FDtau
    a <- as.vector((1 - dij/FDtau) %*% data )  
    data <- data[a!=0]
    v <- c(vv)[a != 0]
    a <- a[a!=0]
    #v <- data/a
    nplus <- sum(data)
    
    return( sum(v*(a/nplus)^q) )
    
  }
  
  FD_mle_part <- function(data, dij, FDtau, q, vv){
    
    dij <- as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] <- FDtau
    a <- as.vector((1 - dij/FDtau) %*% data )  
    data <- data[a!=0]
    v <- c(vv)[a != 0]
    a <- a[a!=0]
    #v <- data/a
    nplus <- sum(data)
    
    if(q==1){
      
      res <- sum(-v*a/nplus*log(a/nplus))
      
    }else{
      
      res <- sum(v*(a/nplus)^q)
      
    }
    
    return( res )
    
  }
  
  FD_mle_rel <- function(data, dij, FDtau, q, wij, vv) {
    
    
    dij <- as.matrix(dij)
    dij[which(dij>FDtau,arr.ind = T)] <- FDtau
    
    
    
    Sub2 <- function(q){
      
      tmp <- sapply(1:nrow(data), function(i){
        
        FD_mle_part(data[i, ], dij = dij, FDtau = FDtau, q = q, vv = vv)
        
        
      }) 
      
      
      if( q == 1 ){
        
        tmp2 <- sapply(1:nrow(data), function(i){
          
          max_mle(data[i, ], dij = dij, FDtau = FDtau, q = q, vv = vv)
          
          
        }) 
        
        res <- exp(sum(unlist(wij) * tmp))
        res <- list('result' = res, 'max' = tmp2)
        
      } else {
        
        res <- sum(unlist(wij) * tmp)^(1/(1-q))
        res <- list('result' = res, 'max' = tmp)
        
      }
      
      
    }
    
    lapply(q, Sub2)
    
  }
  
  Sub <- function(dat, mat, dis, q, FDtau, wij, type, SD = F){
    
    vk <- function(pk_pool,dis, FDtau){
      dis[dis > FDtau] <- FDtau
      d <- dis/FDtau
      pk_pool/(pk_pool%*%(1-d))
      
    }
    
    
    dat <- as.matrix(dat)
    dis <- as.matrix(dis)
    population <- ncol(mat)
    H <- nrow(mat)
    mat[H,] <- paste0(mat[H,], seq(1,ncol(mat)))
    
    
    n_index <- which(rowSums(dat) > 0 )
    dat <- dat[n_index, ]
    dis <- dis[n_index, n_index]
    
    dat <- as.matrix(dat)
    dis <- as.matrix(dis)
    
    n <- sum(dat)
    s <- nrow(dat)
    
    M <- vector("list", H)
    alpha.relative <- list()
    
    check <- apply(dat, 2, function(x) sum(x))
    
    if( sum(check == 0) > 0 ){
      
      index <- which( check == 0)
      dat <- dat[, -index]
      mat <- mat[, -index]
      
      if( !is.character(wij) ){
        
        wij <- wij[-index]
        
      }
      
    }
    
    
    pij <- apply(dat, 2, function(x) x/sum(x))
    
    
    if(FDtau == "dmin"){
      
      mint <- min(dis[dis>0])
      FDtau <- mint
      
    }
    
    if(FDtau == "dmean"){
      
      tmp <- rowSums(pij)/sum(pij)
      meant <- c(tmp%*%dis%*%tmp)
      FDtau <- meant
      
    }
    
    
    if(FDtau == 'dmax'){
      
      maxt <- max(dis)
      FDtau <- maxt
      
    }
    
    
    mat <- apply(mat, 2, as.character)
    index <- rev(lapply(seq_len(H), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x)))))
    
    if(type == 'est'){
      
      wij <- apply(dat, 2, function(x) sum(x)/n)
      
    } else {
      
      if(wij[1] == 'size' ) wij <- apply(dat, 2, function(x) sum(x)/n)
      if(wij[1] == 'equal') wij <- rep(1/ncol(dat), ncol(dat))
      
    }
    
    wij <- wij/sum(wij)
    
    w <- lapply(index, function(x) lapply(x, function(y) sum(wij[y])))
    pij <- apply(dat, 2, function(x) x/sum(x))
    
    weight.p <- pij*t(replicate(s,wij))
    
    p <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
    xx <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
    x <- lapply(index, function(h) lapply(h, function(x) rowSums(as.matrix(weight.p[,x]))))
    pk_plpl <- do.call(rbind, p[[H]])
    v_rel <- vk(pk_plpl, dis, FDtau)
    M <- lapply(index, function(x) lapply(x,length))
    
    est.max <- list()
    mle.max <- list()
    
    for(i in 1:H){
      if(type == 'est'){
        group <- do.call(rbind, xx[[i]])
        zzl = round(z.plus(group, dis, FDtau, rep(1, length(unlist(w[[i]])))))
        xxx <- round(x.plus(group, rep(1, length(unlist(w[[i]])))))
        est.abs <- FD_est_rel(data = group, dij = dis, FDtau = FDtau, q = q, wij = w[[i]], zz = zzl, all_data = xxx)
        est.max[[i]] <- lapply(est.abs, function(x) x$max)
        alpha.relative[[i]] <- unlist(lapply(est.abs, function(x) x$result))
        
      } else {
        # group1 <- do.call(rbind, p[[i]])
        # group2 <- do.call(rbind, xx[[i]])
        group <- do.call(rbind, x[[i]])
        mle.abs <- FD_mle_rel(data = group, dij = dis, FDtau = FDtau, q = q, wij = w[[i]], vv = v_rel)
        mle.max[[i]] <- lapply(mle.abs, function(x) x$max)
        alpha.relative[[i]] <- unlist(lapply(mle.abs, function(x) x$result))
      }
    }
    
    alpha.relative <- do.call(rbind, alpha.relative)
    
    
    H.alpha <- sapply(1:length(q), function(i){
      
      if(q[i] == 1){
        
        log(alpha.relative[,i])
        
      } else {
        
        (1-(population^(q[i]-1)*(alpha.relative[,i]*population)^(1-q[i])))/(q[i]-1)
        
      }
      
    })
    
    
    # if(q!=1) H.alpha <- (1-(population^(q-1)*(alpha.relative*population)^(1-q)))/(q-1) else H.alpha <- log(alpha.relative)
    
    rownames(alpha.relative) <- c(paste0("qFD_alpha",seq_len(H-1)),"qFD_gamma")
    rownames(H.alpha) <- c(paste0("qQ_alpha",seq_len(H-1)),"qQ_gamma")
    
    D.beta <- alpha.relative[-1,]/alpha.relative[-nrow(alpha.relative), ]
    
    D.beta <- matrix(D.beta, ncol = length(q))
    if( SD == F) D.beta[D.beta < 1] <- 1
    
    all.pair <- combn(seq_len(H),2)
    
    if(H == 2){
      pair2 = pair <- all.pair
    } else {
      pair <- cbind(all.pair[,diff(all.pair) == 1], c(1,H))
      pair2 <- pair[ ,-ncol(pair)]
    } 
    
    
    
    differential.relative <- lapply(1:ncol(pair), function(l){
      # up <- ifelse(q!=1, -diff(alpha.relative[j]^(1-q))*population^(1-q), diff(log(alpha.relative[j]))) #*#
      
      j <- pair[,l] 
      
      up <- c(sapply(1:length(q), function(i){
        
        if( q[i] == 1){
          
          diff(log(alpha.relative[j, i]))
          
        } else {
          
          -diff(alpha.relative[j, i]^(1-q[i]))
          
        }
        
      }))
      
      names(up) <- NULL
      
      w.up <- unlist(w[[j[1]]])
      w.down <- unlist(w[[j[2]]])
      num <- sapply(index[[j[1]]], function(a) seq_along(w.down)[sapply(index[[j[2]]], function(b) all(a%in%b))])
      w.down.new <-w.down[num]
      ratio <- w.up/w.down.new
      
      ss <- lapply(1:length(q), function(i){
        
        if( q[i] == 1 ) {
          
          if( type == 'est'){
            
            Max <- -sum(w.up * log(ratio) * est.max[[j[1]]][[i]])
            Gamma.Max <- exp(Max)
            
            
          } else {
            
            Max <- -sum(w.up * log(ratio) * mle.max[[j[1]]][[i]])
            Gamma.Max <- exp(Max)
            
            
          }
          
          
        }else{
          
          if(type == 'est'){
            
            
            
            Max <- sum(w.down.new*ratio*(1-ratio^(q[i]-1))*est.max[[j[1]]][[i]])
            
            Gamma.Max <- sum(w.up*(ratio^(q[i]-1))*est.max[[j[1]]][[i]])^(1/(1-q[i]))
            
          } else {
            
            Max <- sum(w.down.new*ratio*(1-ratio^(q[i]-1))*mle.max[[j[1]]][[i]])
            
            Gamma.Max <- sum(w.up*(ratio^(q[i]-1))*mle.max[[j[1]]][[i]])^(1/(1-q[i]))
            
            
          }
        }
        
        
        c(Max = Max, Gamma.Max = Gamma.Max)
        
      })
      
      return(rbind(up, do.call(cbind, ss)))
      
    })
    
    
    
    
    differential <- do.call(rbind, lapply(differential.relative, function(x) x[1,]/x[2,]))
    if( SD == F) differential[differential < 0] <- 0
    rownames(differential) <- paste0("Delta ","level(",apply(pair, 2,paste, collapse = "|"), ')')
    
    
    D.Beta_max <- sapply(1:length(q), function(i){
      
      if(q[i] == 1){
        
        do.call(rbind, lapply(differential.relative, function(x) x[3,i]))[-nrow(alpha.relative), ]
        
        
      } else {
        
        do.call(rbind, lapply(differential.relative, function(x) x[3,i]))[-nrow(alpha.relative), ] / alpha.relative[-nrow(alpha.relative), i]
        
      }
      
    })
    
    D.Beta_max <- matrix(D.Beta_max, ncol = length(q))
    rownames(D.Beta_max) <- paste0("qFD_Beta_max ", seq_len(nrow(D.beta)))
    
    D.Diss <- sapply(seq_along(q), function(i) Dissimilarity_measure(D.beta[,i], D.Beta_max[,i], q[i], pair2) )
    
    if( SD == F) D.Diss <- abs(D.Diss)
    
    if(H == 2){
      
      beta <- apply(H.alpha, 2, diff)
      
    } else {
      
      beta <- rbind(apply(H.alpha, 2, diff), H.alpha[nrow(H.alpha), ] - H.alpha[1, ])
    }
    
    
    if( SD == F) beta[beta < 0] <- 0
    beta <- matrix(beta, ncol = length(q))
    rownames(beta) <- paste0("qQ_Beta ","level(",apply(pair, 2,paste, collapse = "|"), ')')
    rownames(D.beta) <- paste0("qFD_Beta ", seq_len(nrow(D.beta)))
    
    # result <- data.frame(relative=c(alpha.relative, differential))
    TAU <- c('FDtau' = FDtau)
    
    result <- rbind(H.alpha, beta, differential, alpha.relative, D.beta, D.Beta_max, D.Diss)
    
    # result <- data.frame(relative = c(H.alpha, beta, differential, alpha.relative, D.beta, D.Beta_max, D.Diss, TAU))
    out <- result
    out
    
  }
  
  
  if(datatype == 'incidence'){
    
    su <- as.numeric(dat[1, ])
    dat <- dat[-1, ]
    
  }
  
  res <- Sub(dat, mat, dis, q, FDtau, weight, type)
  
  AUC_func <- function(dat, mat, dis, q, weight, type){
    tau1 <- seq(0.00001, 1, length.out = 30)
    ans <- sapply(q, function(s) {
      qFun = sapply(tau1, function(ta) Sub(dat, mat, dis, s, ta, weight, type, SD = T))
      AUC <- apply(qFun, 1, function(i) {
        LA <- sum(i[seq_along(tau1[-1])]*diff(tau1))
        RA <- sum(i[-1]*diff(tau1))
        mean(c(LA, RA))
      })
      AUC
    })
    ans
  }
  
  if(FDtype == "AUC"){
    res.auc = AUC_func(dat, mat, dis, q, weight, type)
    rownames(res.auc) = rownames(res) 
  }
  if( nboot > 1 ){
    
    if( datatype == 'abundance'){
      
      n <- sum(dat)
      w <- colSums(dat)/n
      P <- Boots.beta(dat, dis, datatype)$pop
      P <- sapply(1:length(w), function(i) P[,i] * w[i])
      P <- P/sum(P)
      dij <- Boots.beta(dat, dis, datatype)$dij
      X <- replicate(nboot, matrix(rmultinom(1, n, c(P)), ncol = length(w)))
      
    }
    
    if( datatype == 'incidence' ){
      
      P <- Boots.beta(rbind(su, dat), dis, datatype)$pop
      dij <- Boots.beta(rbind(su, dat), dis, datatype)$dij
      X <- replicate(nboot, sapply(seq_along(su), function(i) rbinom(nrow(P),c(su[i]), P[,i])))
      
    }
    
    
    if( sum( dij > max(dis) ) > 0 ){
      
      
      dij[dij > max(dis)] <- (dij[dij > max(dis)]/max(dij[dij > max(dis)])) * max(dis)
      
      
    }
    
    dij <- replace(dij, dij == 0, 10^(-10))
    diag(dij) <- 0
    
    
    boot.est <- lapply(1:nboot, function(i){
      
      Sub(X[,,i], mat, dij, q, FDtau, weight, type, SD = T)
      
    })
    est.se <- apply(array(as.numeric(unlist(boot.est)), dim = c(dim(res), nboot)), c(1,2), function(x) ifelse(is.na(sd(x[!is.na(x)])), 0, sd(x[!is.na(x)])))
    rownames(est.se) <- rownames(res)
    
    if(FDtype == "AUC"){
      boot.est.auc <- lapply(1:nboot, function(i){
        
        AUC_func(X[,,i], mat, dij, q, weight, type)
        
      })
      est.se.auc <- apply(array(as.numeric(unlist(boot.est.auc)), dim = c(dim(res.auc), nboot)), c(1,2), function(x) ifelse(is.na(sd(x[!is.na(x)])), 0, sd(x[!is.na(x)])))
      rownames(est.se.auc) <- rownames(res)
    }
    
  } else {
    
    
    est.se <- matrix(NA, dim(res)[1], dim(res)[2])
    rownames(est.se) <- rownames(res)
    
    if(FDtype == "AUC"){
      est.se.auc <- matrix(NA, dim(res.auc)[1], dim(res.auc)[2])
      rownames(est.se.auc) <- rownames(res)
    }
  }
  
  if(FDtype == "AUC"){
    return( list('result' = res.auc, 'S.E.' = est.se.auc) )
  }else{
    return( list('result' = res, 'S.E.' = est.se) )
  }
  
  
}
