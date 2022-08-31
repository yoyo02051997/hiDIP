# DataInfo ----------------------
#' Exhibit basic data information
#'
#' \code{DataInfo}: exhibits basic data information
#'
#' @param data   \code{data.frames}
#' @param diversity selection of diversity type: \code{'TD'} = 'Taxonomic diversity', \code{'PD'} = 'Phylogenetic diversity', and \code{'FD'} = 'Functional diversity'.
#' @param datatype data type of input data: individual-based abundance data \code{(datatype = "abundance")}, sampling-unit-based incidence frequencies data \code{(datatype = "incidence")},
#' @return
#' a data.frame of basic data information incliuding sample size, observed species richness, sample coverage estimate, and the first ten abundance frequency counts.
#' @import dplyr
#' @import Rcpp
#' @import ggplot2
#' @import data.tree
#' @import ape
#' @import ade4
#' @import phytools
#' @import reshape2
#' @import networkD3
#' @import maps
#' @examples
#' ## Taxonomic diversity
#' data(macro)
#' DataInfo(data = macro, diversity = 'TD')
#'
#'
#' ## Phylogenetic diversity
#' data(macro)
#' data(macro_tree)
#' DataInfo(data = macro, diversity = 'PD', tree = macro_tree)
#'
#'
#' ## Functional diversity
#' data(macro)
#' DataInfo(data = macro, diversity = 'FD')
#' @export
DataInfo = function(data, diversity = "TD", datatype = "abundance", tree = NULL){
  if(diversity != "PD"){
    datainf(data = data, datatype = datatype)
  }else{
    datainfphy(data, tree)
  }
}


# hier.taxonomy -------------------------------------------------------------------
#' Decomposition of taxonomy diversity
#' @param data \code{data.frames}
#' @param mat hierarchical structure of data.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param weight weight for relative decomposition.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{20}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param type estimate type: estimate \code{(type = "est")}, empirical estimate \code{(type = "mle")}.Default is \code{"mle"}.
#' @param datatype data type of input data: individual-based abundance data \code{(datatype = "abundance")}, sampling-unit-based incidence frequencies data \code{(datatype = "incidence")}.
#' @param decomposition Relative decomposition: \code{(decomposition = "relative")}, Absolute decomposition: \code{(decomposition = "absolute")}.
#' @import dplyr
#' @import Rcpp
#' @import ggplot2
#' @import data.tree
#' @import ape
#' @import ade4
#' @import phytools
#' @import reshape2
#' @import networkD3
#' @import maps
#' @examples
#' ## Taxonomic diversity
#' data(macro)
#' data(macro_mat)
#' output1 = hier.taxonomy(data = macro, mat = macro_mat, q = seq(0, 2, 0.2))
#' output1
#' @export
hier.taxonomy = function(data, mat, q = seq(0, 2, 0.2), weight = "size", nboot = 20, 
                         conf = 0.95, type = "mle", datatype = "abundance", decomposition = "relative"){
  if(decomposition == "relative"){
    out = hier.taxonomy_rel(data, mat, q, weight = weight, nboot = nboot, conf = conf, type = type, datatype = datatype)
  }else{
    out = hier.taxonomy_abs(data, mat, q, nboot = nboot, conf = conf, type = type, datatype = datatype)
  }
  out
}
# hier.phylogeny -------------------------------------------------------------------
#' Decomposition of phylogeny diversity
#' @param data \code{data.frames}
#' @param mat hierarchical structure of data.
#' @param tree  a phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param weight weight for relative decomposition.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{20}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param type estimate type: estimate \code{(type = "est")}, empirical estimate \code{(type = "mle")}.Default is \code{"mle"}.
#' @param decomposition Relative decomposition: \code{(decomposition = "relative")}, Absolute decomposition: \code{(decomposition = "absolute")}.
#' @import dplyr
#' @import Rcpp
#' @import ggplot2
#' @import data.tree
#' @import ape
#' @import ade4
#' @import phytools
#' @import reshape2
#' @import networkD3
#' @import maps
#' @examples
#' ## Phylogeny diversity
#' data(macro)
#' data(macro_mat)
#' data(macro_tree)
#' output2 = hier.phylogeny(data = macro, mat = macro_mat, tree = macro_tree, q = seq(0, 2, 0.2))
#' output2
#' @export
hier.phylogeny <- function(data, mat, tree, q = seq(0, 2, 0.2), weight = "size", nboot = 20,
                           conf = 0.95, type = "mle", decomposition = "relative"){
  dat <- data[rowSums(data)>0, ]
  if(decomposition == "relative"){
    method = phy.H.rel
  }else{
    method = phy.H.abs
  }
  H <- nrow(mat)
  rtip <- tree$tip.label[!tree$tip.label %in% rownames(dat)]
  #tree$tip.label[!rtree$tip.label %in% rownames(dat)]
  rtree <- drop.tip(tree, rtip)
  tmp <- TranMul(dat, rtree)
  rtreephy <- newick2phylog(convToNewick(rtree))
  if(inherits(weight, "numeric")){
    wij <- weight
  } else if (weight == "size"){
    wij <- colSums(dat)/sum(dat)
  } else  {
    wij <- rep(1/ncol(dat), ncol(dat))
  } 
  est <- sapply(q, function(i) method(dat, mat, tmp, i, rtreephy, wij, type))
  if(nboot!=0){
    H <- nrow(mat)
    boot.est <- bootstrap.q.Beta(data = dat, mat, rtree = rtree, tmp = tmp, q = q, nboot = nboot, wij = wij, type = type, method)
    test <- boot.est[seq_len(H+1), , ]
    #is.infinite(sum(boot.est))
    #test <- boot.est[head(seq_len(dim(boot.est)[1]),H),1,]
    id <- apply(test, 1:2, function(x) {
      bb <- x
      q1 <- quantile(bb,0.25)
      q3 <- quantile(bb,0.75)
      q1 <- q1-1.5*(q3-q1)
      q3 <- q1+1.5*(q3-q1)
      which(bb >= q1 & bb <= q3)
    })
    index <- Reduce(function(x,y) {intersect(x,y)}, id)
    boot.est <- boot.est[ ,,index]
    #boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ][boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ]<0] <- 0
    #boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ][boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ]>1] <- 1
    #dim(boot.est)
    #diff.boot.est <- as.data.frame(apply(boot.est,3, tail, n = H))
    out = lapply(seq_along(q), function(i){x = sapply(seq_len(nrow(boot.est)), function(j) transconf(Bresult = boot.est[j,i,], est = est[j,i], conf)) %>% t()
    colnames(x) = c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
    rownames(x) = rownames(est)
    return(x)}) %>% do.call(rbind,.)
    
    Order.q = rep(q, each = 24)
    
    Method = rep(rownames(out)[1:24], length(q))
    
    rownames(out) = NULL
    
    out = cbind(Method, Order.q, as.data.frame(out), Decomposition = decomposition)
    
    out
    # 
    # CL = t(sapply(seq_len(nrow(boot.est)), function(j) transconf(matrix(boot.est[j,,], nrow=length(q)), est[j,], conf)))
    # rownames(CL) <- rownames(est)
    # colnames(CL) <- c(sapply(paste0("q=",q), function(k) paste(k,c("est", "bt.sd", "LCL", "UCL"), sep = "_")))
    # 
  }else{
    out <- lapply(seq_along(q), function(i){x = data.frame("Estimator" = est[,i], Bootstraps.e = NA, LB = NA, UB = NA) 
    colnames(x) = c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
    return(x)}) %>% do.call(rbind,.)
    
    Order.q = rep(q, each = 24)
    
    Method = rep(rownames(out)[1:24], length(q))
    
    rownames(out) = NULL
    
    out = cbind(Method, Order.q, as.data.frame(out), Decomposition = decomposition)
    out$'Bootstrap S.E.' = as.numeric(out$'Bootstrap S.E.')
    out$LCL = as.numeric(out$LCL)
    out$UCL = as.numeric(out$UCL)
    
    out
  }
  return(out)
}
# hier.functional -------------------------------------------------------------------
#' Decomposition of functional diversity
#' @param data \code{data.frames}
#' @param mat hierarchical structure of data.
#' @param dis species pairwise distance matrix for all species in the pooled assemblage.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param FDtype select FD type: \code{(FDtype = "tau_values")} for FD under specified threshold values, or \code{(FDtype = "AUC")} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is "tau_values".
#' @param FDtau a numerical vector between 0 and 1 specifying tau values (threshold levels).If NULL (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @param weight weight for relative decomposition.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{20}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param type estimate type: estimate \code{(type = "est")}, empirical estimate \code{(type = "mle")}.Default is \code{"mle"}.
#' @param decomposition Relative decomposition: \code{(decomposition = "relative")}, Absolute decomposition: \code{(decomposition = "absolute")}.
#' @import dplyr
#' @import Rcpp
#' @import ggplot2
#' @import data.tree
#' @import ape
#' @import ade4
#' @import phytools
#' @import reshape2
#' @import networkD3
#' @import maps
#' @examples
#' ## Functional diversity
#' data(macro)
#' data(macro_mat)
#' data(macro_dis)
#' output3 = hier.functional(data = macro, mat = macro_mat, dis = macro_dis, q = seq(0, 2, 0.2))
#' output3
#' @export
hier.functional <- function(data, mat, dis, q = seq(0, 2, 0.2), FDtype = "tau_values",
                            FDtau = NULL, weight = "size", nboot = 20, conf = 0.95, type = "mle", 
                            datatype = "abundance", decomposition = "relative"){
  
  if(is.null(FDtau)){
    if(datatype == "incidence"){
      data <- data[-1, ]
    }
    d = dis %>% as.matrix() 
    pij <- apply(data, 2, function(x) x/sum(x))
    tmp <- rowSums(pij)/sum(pij)
    FDtau <- c(tmp%*%d%*%tmp)
  } 
  if(decomposition == "relative"){
    if(FDtype != "AUC"){
      out = lapply(FDtau, function(x){
        qp = hier.func.rel(data, mat, dis, q, FDtype = FDtype, FDtau = x, weight, type, nboot, datatype)
        out <- lapply(seq_along(q), function(i) {
          
          nn <- abs(qnorm(1 - (1 - conf)/2))
          
          res <- cbind(qp$result[,i], qp$S.E.[,i], 
                       ifelse((qp$result[,i] -  nn * qp$S.E.[,i]) < 0, 0, qp$result[,i] -  nn * qp$S.E.[,i]), 
                       qp$result[,i] + nn * qp$S.E.[,i])
          
          colnames(res) <- c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
          
          res = as.data.frame(res)
          
          res$Tau = x
          
          res
        } ) %>% do.call(rbind, .)
        Order.q = rep(q, each = 24)
        
        Method = rep(rownames(out)[1:24], length(q))
        
        out = cbind(Method, Order.q, out, Decomposition = decomposition)
        
        rownames(out) = NULL
        
        out
      }) %>% do.call(rbind, .)
    }else{
      qp = hier.func.rel(data, mat, dis, q, FDtype = FDtype, FDtau = 0.5, weight, type, nboot, datatype)
      out <- lapply(seq_along(q), function(i) {
        
        nn <- abs(qnorm(1 - (1 - conf)/2))
        
        res <- cbind(qp$result[,i], qp$S.E.[,i], 
                     ifelse((qp$result[,i] -  nn * qp$S.E.[,i]) < 0, 0, qp$result[,i] -  nn * qp$S.E.[,i]), 
                     qp$result[,i] + nn * qp$S.E.[,i])
        
        colnames(res) <- c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
        
        res = as.data.frame(res)
        
        res
      } ) %>% do.call(rbind, .)
      Order.q = rep(q, each = 24)
      
      Method = rep(rownames(out)[1:24], length(q))
      
      out = cbind(Method, Order.q, out, Decomposition = decomposition)
      
      rownames(out) = NULL
      
      out
    }
  }else{
    if(FDtype != "AUC"){
      out = lapply(FDtau, function(x){
        qp = hier.func.abs(data, mat, dis, q, FDtype, x, type, nboot, datatype)
        out <- lapply(seq_along(q), function(i) {
          
          nn <- abs(qnorm(1 - (1 - conf)/2))
          
          res <- cbind(qp$result[,i], qp$S.E.[,i], 
                       ifelse((qp$result[,i] -  nn * qp$S.E.[,i]) < 0, 0, qp$result[,i] -  nn * qp$S.E.[,i]), 
                       qp$result[,i] + nn * qp$S.E.[,i])
          
          colnames(res) <- c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
          
          res = as.data.frame(res)
          
          res$Tau = x
          
          res
        } ) %>% do.call(rbind, .)
        Order.q = rep(q, each = 24)
        
        Method = rep(rownames(out)[1:24], length(q))
        
        out = cbind(Method, Order.q, out, Decomposition = decomposition)
        
        rownames(out) = NULL
        
        out
      }) %>% do.call(rbind, .)
    }else{
      qp = hier.func.abs(data, mat, dis, q, FDtype, FDtau = FDtau, type, nboot, datatype)
      out <- lapply(seq_along(q), function(i) {
        
        nn <- abs(qnorm(1 - (1 - conf)/2))
        
        res <- cbind(qp$result[,i], qp$S.E.[,i], 
                     ifelse((qp$result[,i] -  nn * qp$S.E.[,i]) < 0, 0, qp$result[,i] -  nn * qp$S.E.[,i]), 
                     qp$result[,i] + nn * qp$S.E.[,i])
        
        colnames(res) <- c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
        
        res = as.data.frame(res)
        
        res
      } ) %>% do.call(rbind, .)
      Order.q = rep(q, each = 24)
      
      Method = rep(rownames(out)[1:24], length(q))
      
      out = cbind(Method, Order.q, out, Decomposition = decomposition)
      
      rownames(out) = NULL
      
      out
    }
    
  }
  out
}
# gghier_taxonomy -------------------------------------------------------------------
#' ggplot2 extension for outcome from \code{hier_taxonomy}
#' @param outcome a list object computed by \code{hier_taxonomy}.
#' @param method \code{(method = 1)} diversity(alpha, gamma) based on Tsallis entropy (1988);\code{(method = 2)} beta diversity based on additive decomposition;
#' \code{(method = 3)} dissimilarity measure based on additive decomposition;\code{(method = 4)} diversity(alpha, gamma) based on Hill Number (1973);
#' \code{(method = 5)} beta diversity based on multiplicative  decomposition;\code{(method = 6)} dissimilarity measure based on multiplicative decomposition.
#' @examples
#' ## Taxonomy diversity
#' data(macro)
#' data(macro_mat)
#' output1 = hier.taxonomy(data = macro, mat = macro_mat, q = seq(0, 2, 0.2))
#' gghier_taxonomy(output1, method = 1)
#' @export
gghier_taxonomy = function(outcome, method = 1){
  if(method == 1){
    outcome = outcome[c(grep(c("qH_alpha"), outcome$Method),
                        grep(c("qH_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 2){
    outcome = outcome[grep("qH_Beta level", outcome$Method),]
  }else if(method == 3){
    outcome = outcome[grep("Delta level", outcome$Method),]
  }else if(method == 4){
    outcome = outcome[c(grep(c("qD_alpha"), outcome$Method),
                        grep(c("qD_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 5){
    outcome = outcome[grep("qD_Beta ", outcome$Method),]
  }else if(method == 6){
    outcome = outcome[grep("1-", outcome$Method),]
  }
  if(method != 6){
    out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                              fill = Method))+
      geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                              ymax = UCL, fill = Method), linetype = 0, 
                                          alpha = 0.2)
  }else{
    outcome$group = "1-C"
    outcome$group[grep("1-UqN", outcome$Method)] = "1-U"
    outcome$group[grep("1-SqN", outcome$Method)] = "1-S"
    outcome$group[grep("1-VqN", outcome$Method)] = "1-V"
    out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                              fill = Method))+
      facet_grid(group~.)+
      geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                              ymax = UCL, fill = Method), linetype = 0, 
                                          alpha = 0.2)
  }
  out = out  + theme_bw() + 
     theme(legend.position = "bottom",legend.box = "vertical", 
                                      legend.key.width = unit(1.2, "cm"), 
                                      legend.title = element_blank(), 
                                      legend.margin = margin(0, 0, 0, 0), 
                                      legend.box.margin = margin(-10, -10, -5, -10), 
                                      text = element_text(size = 16), 
                                      plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) + 
    guides(linetype = guide_legend(keywidth = 2.5))
  return(out)
}

# gghier_phylogeny -------------------------------------------------------------------
#' ggplot2 extension for outcome from \code{hier_phylogeny}
#' @param outcome a list object computed by \code{hier_phylogeny}.
#' @param method \code{(method = 1)} diversity(alpha, gamma) based on Tsallis entropy (1988);\code{(method = 2)} beta diversity based on additive decomposition;
#' \code{(method = 3)} dissimilarity measure based on additive decomposition;\code{(method = 4)} diversity(alpha, gamma) based on Hill Number (1973);
#' \code{(method = 5)} beta diversity based on multiplicative  decomposition;\code{(method = 6)} dissimilarity measure based on multiplicative decomposition.
#' @examples
#' ## Phylogeny diversity
#' data(macro)
#' data(macro_mat)
#' data(macro_tree)
#' output2 = hier.phylogeny(data = macro, mat = macro_mat, tree = macro_tree, q = seq(0, 2, 0.2))
#' gghier_phylogeny(output2, method = 1)
#' @export
gghier_phylogeny = function(outcome, method = 1){
  if(method == 1){
    outcome = outcome[c(grep(c("qI_alpha"), outcome$Method),
                        grep(c("qI_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 2){
    outcome = outcome[grep("qI_Beta level", outcome$Method),]
  }else if(method == 3){
    outcome = outcome[grep("Delta level", outcome$Method),]
  }else if(method == 4){
    outcome = outcome[c(grep(c("qPD_alpha"), outcome$Method),
                        grep(c("qPD_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 5){
    outcome = outcome[grep("qPD_Beta ", outcome$Method),]
  }else if(method == 6){
    outcome = outcome[grep("1-", outcome$Method),]
  }
  if(method != 6){
    out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                              fill = Method))+
      geom_line(size = 1.5) + 
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Method), linetype = 0, 
                  alpha = 0.2)
  }else{
    outcome$group = "1-C"
    outcome$group[grep("1-UqN", outcome$Method)] = "1-U"
    outcome$group[grep("1-SqN", outcome$Method)] = "1-S"
    outcome$group[grep("1-VqN", outcome$Method)] = "1-V"
    out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                              fill = Method))+
      facet_grid(group~.)+
      geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                              ymax = UCL, fill = Method), linetype = 0, 
                                          alpha = 0.2)
  }
  
  out = out  + theme_bw() + 
      theme(legend.position = "bottom",legend.box = "vertical", 
                                        legend.key.width = unit(1.2, "cm"), 
                                        legend.title = element_blank(), 
                                        legend.margin = margin(0, 0, 0, 0), 
                                        legend.box.margin = margin(-10, -10, -5, -10), 
                                        text = element_text(size = 16), 
                                        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) + 
    guides(linetype = guide_legend(keywidth = 2.5))
  return(out)
}
# gghier_functional -------------------------------------------------------------------
#' ggplot2 extension for outcome from \code{hier_functional}
#' @param outcome a list object computed by \code{hier_functional}.
#' @param method \code{(method = 1)} diversity(alpha, gamma) based on Tsallis entropy (1988);\code{(method = 2)} beta diversity based on additive decomposition;
#' \code{(method = 3)} dissimilarity measure based on additive decomposition;\code{(method = 4)} diversity(alpha, gamma) based on Hill Number (1973);
#' \code{(method = 5)} beta diversity based on multiplicative  decomposition;\code{(method = 6)} dissimilarity measure based on multiplicative decomposition.
#' @param profile a selection of profile versus to diversity.q-profile \code{(profile = "q")};tau-profil \code{(prifile = "tau")}
#' @examples
#' ## Functional diversity
#' data(macro)
#' data(macro_mat)
#' data(macro_dis)
#' output3 = hier.functional(data = macro, mat = macro_mat, dis = macro_dis, q = c(0, 1, 2), FDtau = c(0.2, 0.4, 0.6, 0.8, 1))
#' gghier_functional(output3, method = 1, profile = "q")
#' output4 = hier.functional(data = macro, mat = macro_mat, dis = macro_dis, q = c(0, 1, 2), FDtype = "AUC")
#' gghier_functional(output4, method = 1, profile = "q")
#' @export
gghier_functional = function(outcome, method = 1, profile = "q"){
  
  if (sum(colnames(outcome)[1:7] == c("Tau")) == 
      0) {
    class = "AUC"
  }else{
    class = "FD"
  }
  
  if(method == 1){
    outcome = outcome[c(grep(c("qQ_alpha"), outcome$Method),
                        grep(c("qQ_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 2){
    outcome = outcome[grep("qQ_Beta level", outcome$Method),]
  }else if(method == 3){
    outcome = outcome[grep("Delta level", outcome$Method),]
  }else if(method == 4){
    outcome = outcome[c(grep(c("qFD_alpha"), outcome$Method),
                        grep(c("qFD_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 5){
    outcome = outcome[grep("qFD_Beta ", outcome$Method),]
  }else if(method == 6){
    outcome = outcome[grep("1-", outcome$Method),]
  }
  if(profile == "q"){
    if(class != "AUC"){
      outcome$Tau = paste("Tau = ", round(outcome$Tau, 3), 
                          sep = "")
    }
    if(method != 6){
      out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                                fill = Method))+
        geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                                ymax = UCL, fill = Method), linetype = 0, 
                                            alpha = 0.2)
      if(class != "AUC"){
        out = out + facet_grid(. ~ Tau, scales = "free_y")
      }
    }else{
      outcome$group = "1-C"
      outcome$group[grep("1-UqN", outcome$Method)] = "1-U"
      outcome$group[grep("1-SqN", outcome$Method)] = "1-S"
      outcome$group[grep("1-VqN", outcome$Method)] = "1-V"
      out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                                fill = Method))+
        geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                                ymax = UCL, fill = Method), linetype = 0, 
                                            alpha = 0.2)
      if(class != "AUC"){
        out = out + facet_grid(group~Tau)
      }else{
        out = out + facet_grid(group~.)
      }
      
    }
  }else{
    outcome$Order.q = paste("Order q = ", outcome$Order.q, 
                            sep = "")
    if(method != 6){
      out = ggplot(outcome, aes(x = Tau, y = Estimator, colour = Method, 
                                fill = Method))+
        geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                                ymax = UCL, fill = Method), linetype = 0, 
                                            alpha = 0.2)+
        facet_grid(. ~ Order.q, scales = "free_y")
    }else{
      outcome$group = "1-C"
      outcome$group[grep("1-UqN", outcome$Method)] = "1-U"
      outcome$group[grep("1-SqN", outcome$Method)] = "1-S"
      outcome$group[grep("1-VqN", outcome$Method)] = "1-V"
      out = ggplot(outcome, aes(x = Tau, y = Estimator, colour = Method, 
                                fill = Method))+
        geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                                ymax = UCL, fill = Method), linetype = 0, 
                                            alpha = 0.2)+
        facet_grid(group~Order.q)
    }
  }
 
   
  out = out  + theme_bw() + 
      theme(legend.position = "bottom",legend.box = "vertical", 
                                        legend.key.width = unit(1.2, "cm"), 
                                        legend.title = element_blank(), 
                                        legend.margin = margin(0, 0, 0, 0), 
                                        legend.box.margin = margin(-10, -10, -5, -10), 
                                        text = element_text(size = 16), 
                                        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) + 
    guides(linetype = guide_legend(keywidth = 2.5))
  return(out)
}
