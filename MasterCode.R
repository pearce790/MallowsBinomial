#### Mallows Model Master Functions ####

#### Libraries ####
library(gtools)
library(nloptr)
library(lpSolve)

#### Mallows' and Mallows-Binomial Base Functions ####
d_R <- function(pi,pi0){
  ## Calculate the Kendall distance between two rankings, pi and pi0
  ## Inputs: pi is a (partial or complete) ranking of length R from a collection of J items, R<=J,
  ##         pi0 is a (complete) ranking of length J from a collection of J items.
  R <- length(pi)
  J <- length(pi0)
  if(length(setdiff(1:J,pi0))!=0){stop("pi0 must be a complete ranking")}
  if(any(pi>J,na.rm=T)){stop("pi cannot contain items not in pi0")}
  if(R>J){stop("R must be <= J")}
  dist <- 0
  for(r in 1:R){
    dist <- dist + (which(pi0 == pi[r]) - 1)
    pi0 <- setdiff(pi0,pi[r])
  }
  return(dist)
}
psi_R <- function(theta,J,R,log=FALSE){
  ## Calculate the normalizing constant of a (J,R) Mallows' model (J items; rankings of length R<=J)
  ## Inputs: theta = Mallows' scale parameter, R = length of partial rankings, 
  ##         J = number of total items, log = boolean to return log(psi_R).
  if(R>J){stop("R must be <= J")}
  if(R<=0 | J <= 0){stop("R and J must be positive integers")}
  if(theta<=0){stop("theta must be >=0")}
  psi_R <- prod((1-exp(-theta*(J-1:R+1)))/(1-exp(-theta)))
  if(log){psi_R <- log(psi_R)}
  return(psi_R)
}
dmm <- function(Pi,pi0,theta,log=FALSE){
  ## Calculate the joint likelihood of observations Pi from a Mallows' model with parameters pi0 and theta
  ## Inputs: Pi = I x R matrix of (partial) rankings, where I = number of reviewers and R = length of rankings,
  ##         pi0 = central ranking, theta = scale parameter, log = boolean to return loglikelihood.
  if(!is.matrix(Pi)){stop("Pi must be a matrix of (partial) rankings")}
  I <- nrow(Pi)
  R <- ncol(Pi)
  J <- length(pi0)
  if(log){
    return(-theta*sum(unlist(apply(Pi,1,function(pi){d_R(pi,pi0)})))-I*psi_R(theta,J,R,log=TRUE))
  }else{
    return(exp(-theta*sum(unlist(apply(Pi,1,function(pi){d_R(pi,pi0)}))))/(psi_R(theta,J,R,log=FALSE)^I))
  }
}
dmbm <- function(Pi,X,p,pi0,theta,M,log=FALSE){
  ## Calculate the joint likelihood of observations Pi and X from a Mallows-Binomial model with paramters p, theta.
  ## Inputs: Pi = I x R matrix of (partial) rankings, where I = number of reviewers and R = length of rankings,
  ##         X = I x J matrix of scores, where J = number of items, p = item quality vector of length J, 
  ##         theta = positive scale parameter, M = maximum integer score, log = boolean to return loglikelihood.
  if(length(p)!=ncol(X)){stop("p must equal ncol(X)")}
  if(log){
    return(dmm(Pi,pi0,theta,log)+sum(apply(X,1,function(x){dbinom(x,M,p,log)}),na.rm=T))
  }else{return(dmm(Pi,pi0,theta,log)*prod(apply(X,1,function(x){dbinom(x,M,p,log)}),na.rm=T))}
}
sample_mallows <- function(n,pi0,theta,R=length(pi0)){
  tmp <- function(pi0,theta){
    J <- length(pi0)
    Vj <- pi1 <- c()
    for(j in 1:(J-1)){
      probs <- exp(-theta*(0:(J-j)))/sum(exp(-theta*(0:(J-j))))
      Vj <- c(Vj,sample.int(J-j+1,size=1,prob=probs)-1)
    }
    Vj <- c(Vj,0)
    for(j in 1:J){
      pi1 <- c(pi1,pi0[Vj[j]+1])
      pi0 <- setdiff(pi0,pi0[Vj[j]+1])
    }
    return(pi1)
  }
  t(replicate(n,tmp(pi0,theta)))[,1:R]
}

#### Estimation Functions ####
phat_conditional <- function(X,M,order){
  ## Optimize for MLE of item qualities p conditional on (incomplete) ordering of the qualities
  ## Inputs: X = I x J matrix of scores, where I = number of reviewers and J = number of items,
  ##         M = maximum integer score, order = (incomplete) ordering of items.
  I <- nrow(X)
  J <- ncol(X)
  Xbar <- apply(X,2,function(x){mean(x,na.rm=TRUE)})
  opt_f <- function(p,Xbar,M,order,J){
    # quantity that must be optimized with respect to p
    return(-sum(log(p)*Xbar + log(1-p)*(M-Xbar),na.rm=TRUE))
  }
  opt_g <- function(p,Xbar,M,order,J){
    # contraints, which are met when all quantities returned are nonpositive
    return(c(-diff(p[order]),p[order[length(order)]]-p[setdiff(1:J,order)]))
  }
  res <- nloptr(x0=Xbar/M,
                eval_f = opt_f, eval_g_ineq = opt_g,
                lb = rep(0,J), ub = rep(1,J),
                Xbar = Xbar, M = M, order = order, J = J,
                opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_abs"=0.000001,maxeval=1000))
  res$solution
}
theta_conditional <- function(D,J,R){
  ## Optimize for MLE of scale parameter theta conditional on some value D.
  ## Inputs: D = average distance between the (partial) rankings and the central ranking,
  ##         J = number of items, R = size of (partial) rankings, order = (incomplete) ordering of items.
  if(D<0){stop("D cannot be less than 0")}
  if(D==0){return(1000)
  }else{
    return(optimize(function(theta){theta*D+psi_R(theta,J,R,log=TRUE)},
                    c(0,10^8))$minimum)
  }
}
getQ <- function(Pi,J){
  ## Calculate the Q matrix given a collection of (partial) rankings
  ## Inputs: Pi = I x R matrix of partial rankings, J = total number of items
  Q <- matrix(NA,nrow=J,ncol=J)
  for(i in 1:J){for(j in 1:J){
    Q[i,j] <- mean(apply(Pi,1,function(pi){
      if(i %in% pi & j %in% pi){return(which(pi==i)<which(pi==j))}
      if(i %in% pi & !(j %in% pi)){return(TRUE)}
      if(!(i %in% pi) & j %in% pi){return(FALSE)}
      if(!(i %in% pi) & !(j %in% pi)){return(FALSE)}
    }))
  }}
  return(Q)
}
totalcostheuristic_MB <- function(Q,X,M,order,J,R){
  ## Calculate the total cost heuristic given Q (sufficient statistics for Pi), X, and some (incomplete) order of items.
  ## Inputs: Q = a J x J empirical probability comparison matrix for items in the collection, X = I x J score matrix,
  ##         M = maximum integer score, order = (incomplete) ordering of items, R = length of partial rankings,
  ##         J = total number of items.
  
  #calculate conditionally-optimal p
  phat <- phat_conditional(X,M,order)
  
  #calculate conditionally-optimal theta
  S <- setdiff(1:J,order)
  D <- 0 
  if(length(S)>=2){D <- D + sum(apply(combinations(length(S),2,S),1,function(uv){min(Q[uv[1],uv[2]],Q[uv[2],uv[1]])}))}
  for(i in 1:length(order)){D <- D + sum(Q[setdiff(1:J,order[1:i]),order[i]])}
  thetahat <- theta_conditional(D,J,R)
  
  #calculate total cost heuristic
  Xbar <- apply(X,2,function(x){mean(x,na.rm=TRUE)})
  return(thetahat*D+psi_R(thetahat,J,R,log=TRUE)-sum(log(phat)*Xbar + log(1-phat)*(M-Xbar),na.rm=TRUE))
}
idx <- function(i,j,J){J*(i-1)+j}
lp_heuristic <- function(Q,I,J,order){
  coeffs <- Q*I
  S <- setdiff(1:J,order)
  coeffs <- coeffs[S,S]
  if(length(S)==1){
    pi0 <- c(order,S)
    dist <- sum(apply(Pi,1,function(rank){d_R(rank,pi0)}))
    return(list(pi0=pi0,dist=dist))
  }
  if(length(S)>=2){
    pairs <- combinations(length(S),2)
    pairs_constraints <- matrix(0,nrow=length(S)*(length(S)-1)/2,ncol=length(S)^2)
    for(row in 1:choose(length(S),2)){
      items <- pairs[row,]
      pairs_constraints[row,c(idx(items[1],items[2],length(S)),
                              idx(items[2],items[1],length(S)))] <- 1
    }
  }else{pairs_constraints <- matrix(0,nrow=0,ncol=length(S)^2)}
  if(length(S)>=3){
    trios <- combinations(length(S),3)
    trios_constraints <- matrix(0,nrow=nrow(trios),ncol=length(S)^2)
    for(row in 1:choose(length(S),3)){
      items <- trios[row,]
      trios_constraints[row,c(idx(items[1],items[2],length(S)),
                              idx(items[2],items[3],length(S)),
                              idx(items[3],items[1],length(S)))] <- 1
    }
  }else{trios_constraints <- matrix(0,nrow=0,ncol=length(S)^2)}
  
  sol <- lp(direction = "min",
            objective = -as.vector(t(coeffs)),
            const.mat = rbind(pairs_constraints,trios_constraints),
            const.dir = c(rep("==",nrow(pairs_constraints)),
                          rep(">=",nrow(trios_constraints))),
            const.rhs = 1)
  vals <- rep(0,J)
  vals[S] <- apply(matrix(sol$solution,byrow=T,ncol=length(S)),2,sum)
  vals[order] <- seq(-length(order),-1,length=length(order))
  vals <- round(vals)
  if(length(vals)-length(unique(vals))>0){
    cond <- TRUE
    Pi0 <- unique(t(apply(matrix(1:50),1,function(i){
      rank <- rank(vals,ties.method="random")
      unlist(lapply(1:J,function(j){which(rank==j)}))
    })))
    while(cond){
      Pi0 <- unique(rbind(Pi0,t(apply(matrix(1:50),1,function(i){
        rank <- rank(vals,ties.method="random")
        unlist(lapply(1:J,function(j){which(rank==j)}))
      }))))
      if(nrow(Pi0)==prod(factorial(table(vals)))){cond<-FALSE}
    }
    dists <- apply(Pi0,1,function(pi0){sum(unlist(apply(Pi,1,function(rank){d_R(rank,pi0)})))})
    return(list(pi0=Pi0[which.min(dists),],
                dist=dists[which.min(dists)]))
  }else{
    pi0 <- order(vals)
    dist <- sum(unlist(apply(Pi,1,function(rank){d_R(rank,pi0)})))
    return(list(pi0=pi0,dist=dist))
  }
}
totalcostheuristic_MB_LP <- function(Q,X,M,order,J,R,I){
  ## Calculate the total cost heuristic given Q (sufficient statistics for Pi), X, and some (incomplete) order of items.
  ## This functions uses a linear programming approach to find the heuristic
  ## Inputs: Q = a J x J empirical probability comparison matrix for items in the collection, X = I x J score matrix,
  ##         M = maximum integer score, order = (incomplete) ordering of items, R = length of partial rankings,
  ##         J = total number of items, I = total number of reviewers
  
  #calculate conditionally-optimal p
  phat <- phat_conditional(X,M,order)
  
  #calculate conditionally-optimal theta
  D <- lp_heuristic(Q,I,J,order)$dist/I
  thetahat <- theta_conditional(D,J,R)
  
  #calculate total cost heuristic
  Xbar <- apply(X,2,function(x){mean(x,na.rm=TRUE)})
  return(thetahat*D+psi_R(thetahat,J,R,log=TRUE)-sum(log(phat)*Xbar + log(1-phat)*(M-Xbar),na.rm=TRUE))
}
totalcostheuristic_M <- function(Q,order,J,R){
  ## Calculate the total cost heuristic given Q (sufficient statistics for Pi) and some (incomplete) order of items.
  ## Inputs: Q = a J x J empirical probability comparison matrix for items in the collection,
  ##         order = (incomplete) ordering of items, R = length of partial rankings, J = total number of items.

  #calculate conditionally-optimal theta
  S <- setdiff(1:J,order)
  D <- 0 
  if(length(S)>=2){D <- D + sum(apply(combinations(length(S),2,S),1,function(uv){min(Q[uv[1],uv[2]],Q[uv[2],uv[1]])}))}
  for(i in 1:length(order)){D <- D + sum(Q[setdiff(1:J,order[1:i]),order[i]])}
  thetahat <- theta_conditional(D,J,R)
  
  #calculate total cost heuristic
  return(thetahat*D+psi_R(thetahat,J,R,log=TRUE))
}
ASTAR_MB <- function(Pi,X,M,J,R){
  ## Find the MLE of a (J,R) Mallows-Binomial model
  ## Inputs: Pi = I x R matrix of rankings, X = I x J matrix of scores, 
  ##         M = maximum integer score, J = number of items in collection, R = length of partial ranking
  
  # initialize
  if(nrow(Pi) != nrow(X)){stop("Pi and X must have same number of rows")}
  I <- nrow(X)
  Q <- getQ(Pi,J)
  
  # set up and run through initial layer of tree
  nodes <- as.data.frame(matrix(c(1:J,rep(NA,J*(J-1))),nrow=J))
  names(nodes) <- c(paste0(1:(J-1)),"TCH")
  for(row in 1:J){
    order <- nodes[row,1:(J-1)]
    order <- as.numeric(order[!is.na(order)])
    nodes[row,J] <- totalcostheuristic_MB(Q,X,M,order,J,R)
  }
  num_nodes <- J
  
  # iteratively run through future layers until terminal node (i.e., the MLE) is found
  while(TRUE){
    #find current parent node
    curr_index <- which.min(nodes[,J])
    curr_parent <- nodes[curr_index,1:(J-1)]
    (curr_parent <- as.numeric(curr_parent[!is.na(curr_parent)]))
    
    #check if parent is terminal; if so, stop!
    if(length(curr_parent) == (J-1)){break()}
    
    #get queue (next level down from parent)
    queue <- as.data.frame(matrix(NA,nrow=J-length(curr_parent),ncol=J))
    names(queue) <- c(paste0(1:(J-1)),"TCH")
    queue[,1:length(curr_parent)] <- rep(curr_parent,each=J-length(curr_parent))
    queue[,length(curr_parent)+1] <- setdiff(1:J,curr_parent)
    
    #go through queue
    for(row in 1:nrow(queue)){
      order <- queue[row,1:(J-1)]
      order <- as.numeric(order[!is.na(order)])
      queue[row,J] <- totalcostheuristic_MB(Q,X,M,order,J,R)
    }
    nodes <- rbind(nodes,queue)
    num_nodes <- num_nodes + nrow(queue)
    
    #remove parent node from queue
    nodes <- nodes[-curr_index,]
    nodes
  }
  
  # output results
  pi0_mle = as.numeric(nodes[which.min(nodes[,J]),1:(J-1)])
  pi0_mle = c(pi0_mle,setdiff(1:J,pi0_mle))
  p_mle = phat_conditional(X,M,pi0_mle) + 
    seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(pi0_mle==j)}))] #additive term to correct ties
  theta_mle = theta_conditional(
    mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],pi0_mle)}))),J,R)
  return(list(pi0_mle = pi0_mle,p_mle=p_mle,theta_mle=theta_mle,
              num_nodes=num_nodes,nodes=nodes))
}
ASTAR_MB_LP <- function(Pi,X,M,J,R){
  ## Find the MLE of a (J,R) Mallows-Binomial model
  ## Inputs: Pi = I x R matrix of rankings, X = I x J matrix of scores, 
  ##         M = maximum integer score, J = number of items in collection, R = length of partial ranking
  
  # initialize
  if(nrow(Pi) != nrow(X)){stop("Pi and X must have same number of rows")}
  I <- nrow(X)
  Q <- getQ(Pi,J)
  
  # set up and run through initial layer of tree
  nodes <- as.data.frame(matrix(c(1:J,rep(NA,J*(J-1))),nrow=J))
  names(nodes) <- c(paste0(1:(J-1)),"TCH")
  for(row in 1:J){
    order <- nodes[row,1:(J-1)]
    order <- as.numeric(order[!is.na(order)])
    nodes[row,J] <- totalcostheuristic_MB_LP(Q,X,M,order,J,R,I)
  }
  num_nodes <- J
  
  # iteratively run through future layers until terminal node (i.e., the MLE) is found
  while(TRUE){
    #find current parent node
    curr_index <- which.min(nodes[,J])
    curr_parent <- nodes[curr_index,1:(J-1)]
    (curr_parent <- as.numeric(curr_parent[!is.na(curr_parent)]))
    
    #check if parent is terminal; if so, stop!
    if(length(curr_parent) == (J-1)){break()}
    
    #get queue (next level down from parent)
    queue <- as.data.frame(matrix(NA,nrow=J-length(curr_parent),ncol=J))
    names(queue) <- c(paste0(1:(J-1)),"TCH")
    queue[,1:length(curr_parent)] <- rep(curr_parent,each=J-length(curr_parent))
    queue[,length(curr_parent)+1] <- setdiff(1:J,curr_parent)
    
    #go through queue
    for(row in 1:nrow(queue)){
      order <- queue[row,1:(J-1)]
      order <- as.numeric(order[!is.na(order)])
      queue[row,J] <- totalcostheuristic_MB_LP(Q,X,M,order,J,R,I)
    }
    nodes <- rbind(nodes,queue)
    num_nodes <- num_nodes + nrow(queue)
    
    #remove parent node from queue
    nodes <- nodes[-curr_index,]
    nodes
  }
  
  # output results
  pi0_mle = as.numeric(nodes[which.min(nodes[,J]),1:(J-1)])
  pi0_mle = c(pi0_mle,setdiff(1:J,pi0_mle))
  p_mle = phat_conditional(X,M,pi0_mle) + 
    seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(pi0_mle==j)}))] #additive term to correct ties
  theta_mle = theta_conditional(
    mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],pi0_mle)}))),J,R)
  return(list(pi0_mle = pi0_mle,p_mle=p_mle,theta_mle=theta_mle,
              num_nodes=num_nodes,nodes=nodes))
}
ASTAR_M <- function(Pi,J,R){
  ## Find the MLE of a Mallows' model (no scores)
  ## Inputs: Pi = I x R matrix of rankings, J = number of items in collection, R = length of partial ranking
  
  # initialize
  I <- nrow(Pi)
  if(ncol(Pi) != R){stop("ncol(Pi) must equal R")}
  if(max(Pi,na.rm=T)>J){stop("Pi cannot contain elements greater than J")}
  Q <- getQ(Pi,J)
  
  # set up and run through initial layer of tree
  nodes <- as.data.frame(matrix(c(1:J,rep(NA,J*(J-1))),nrow=J))
  names(nodes) <- c(paste0(1:(J-1)),"TCH")
  for(row in 1:J){
    order <- nodes[row,1:(J-1)]
    order <- as.numeric(order[!is.na(order)])
    nodes[row,J] <- totalcostheuristic_M(Q,order,J,R)
  }
  num_nodes <- J
  
  # iteratively run through future layers until terminal node (i.e., the MLE) is found
  while(TRUE){
    #find current parent node
    curr_index <- which.min(nodes[,J])
    curr_parent <- nodes[curr_index,1:(J-1)]
    (curr_parent <- as.numeric(curr_parent[!is.na(curr_parent)]))
    
    #check if parent is terminal; if so, stop!
    if(length(curr_parent) == (J-1)){break()}
    
    #get queue (next level down from parent)
    queue <- as.data.frame(matrix(NA,nrow=J-length(curr_parent),ncol=J))
    names(queue) <- c(paste0(1:(J-1)),"TCH")
    queue[,1:length(curr_parent)] <- rep(curr_parent,each=J-length(curr_parent))
    queue[,length(curr_parent)+1] <- setdiff(1:J,curr_parent)
    
    #go through queue
    for(row in 1:nrow(queue)){
      order <- queue[row,1:(J-1)]
      order <- as.numeric(order[!is.na(order)])
      queue[row,J] <- totalcostheuristic_M(Q,order,J,R)
    }
    nodes <- rbind(nodes,queue)
    num_nodes <- num_nodes + nrow(queue)
    
    #remove parent node from queue
    nodes <- nodes[-curr_index,]
    nodes
  }
  
  # output results
  pi0_mle = as.numeric(nodes[which.min(nodes[,J]),1:(J-1)])
  pi0_mle = c(pi0_mle,setdiff(1:J,pi0_mle))
  unranked <- setdiff(1:J,as.numeric(Pi))
  pi0_mle = c(setdiff(pi0_mle,unranked),sample(unranked,length(unranked)))
  theta_mle = theta_conditional(
    mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],pi0_mle)}))),J,R)
  return(list(pi0_mle = pi0_mle,unranked=unranked,theta_mle=theta_mle,
              num_nodes=num_nodes,nodes=nodes))
}
FV_algorithm <- function(Pi,X,M,J,R,net=1){
  central_rank <- order(apply(getQ(Pi,J),2,sum),partial=apply(X,2,function(score){sum(score,na.rm=T)}))
  central_score<- order(apply(X,2,function(score){sum(score,na.rm=T)}),partial=apply(getQ(Pi,J),2,sum))
  
  diffs <- unique(matrix(c(central_rank,central_score),byrow=T,ncol=J,nrow=2))
  if(net>0){
    for(net_iter in 1:net){
      for(row in 1:nrow(diffs)){
        diffs <- rbind(diffs,matrix(unlist(lapply(1:(J-1),function(j){
          new <- diffs[row,]
          new[j:(j+1)] <- diffs[row,(j+1):j]
          new
        })),byrow=T,ncol=J))
      }
      diffs <- unique(diffs)
    }
  }
  
  
  costs <- apply(diffs,1,function(order){
    p <- phat_conditional(X,M,order) + 
      seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(order==j)}))] #additive term to correct ties
    theta <- theta_conditional(
      mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],order)}))),J,R)
    dmbm(Pi,X,round(p,6),order,theta,M,log=T)
  })
  num_nodes <- length(costs)
  
  
  pi0_FV <- diffs[which.max(costs),]
  
  p_FV <- phat_conditional(X,M,pi0_FV) + 
    seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(pi0_FV==j)}))] #additive term to correct ties
  theta_FV <- theta_conditional(
    mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],pi0_FV)}))),J,R)
  
  return(list(pi0_FV = pi0_FV,p_FV=p_FV,theta_FV=theta_FV,
              num_nodes = num_nodes))
}
greedy_algorithm <- function(Pi,X,M,J,R,net=0){
  ## Find the approximate MLE of a (J,R) Mallows-Binomial model
  ## Inputs: Pi = I x R matrix of rankings, X = I x J matrix of scores, 
  ##         M = maximum integer score, J = number of items in collection, R = length of partial ranking
  
  # initialize
  if(nrow(Pi) != nrow(X)){stop("Pi and X must have same number of rows")}
  I <- nrow(X)
  Q <- getQ(Pi,J)
  
  curr_ranking <- c()
  num_nodes <- 0
  while(length(curr_ranking)<(J-1)){
    S <- setdiff(1:J,curr_ranking)
    cost <- rep(NA,length(S))
    for(i in 1:length(S)){
      num_nodes <- num_nodes + 1
      s <- S[i]
      try_ranking <- c(curr_ranking,s)
      cost[i] <- totalcostheuristic_MB(Q,X,M,try_ranking,J,R)
    }
    curr_ranking <- c(curr_ranking,S[which.min(cost)])
  }
  if(net>0){
    diffs <- matrix(c(curr_ranking,setdiff(1:J,curr_ranking)),nrow=1)
    for(net_iter in 1:net){
      for(row in 1:nrow(diffs)){
        diffs <- rbind(diffs,matrix(unlist(lapply(1:(J-1),function(j){
          new <- diffs[row,]
          new[j:(j+1)] <- diffs[row,(j+1):j]
          new
        })),byrow=T,ncol=J))
      }
      diffs <- unique(diffs)
    }
    costs <- apply(unique(diffs),1,function(order){
      p <- phat_conditional(X,M,order) + 
        seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(order==j)}))] #additive term to correct ties
      theta <- theta_conditional(
        mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],order)}))),J,R)
      dmbm(Pi,X,round(p,6),order,theta,M,log=T)
    })
    num_nodes <- num_nodes + nrow(diffs)
    pi0_greedy <- diffs[which.max(costs),]
  }else {pi0_greedy <- c(curr_ranking,setdiff(1:J,curr_ranking))}
  
  
  
  p_greedy <- phat_conditional(X,M,pi0_greedy) + 
    seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(pi0_greedy==j)}))] #additive term to correct ties
  theta_greedy <- theta_conditional(
    mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],pi0_greedy)}))),J,R)
  
  return(list(pi0_greedy = pi0_greedy,p_greedy=p_greedy,theta_greedy=theta_greedy,
              num_nodes = num_nodes))
}
greedylocal_algorithm <- function(Pi,X,M,J,R){
  ## Find the approximate MLE of a (J,R) Mallows-Binomial model
  ## Inputs: Pi = I x R matrix of rankings, X = I x J matrix of scores, 
  ##         M = maximum integer score, J = number of items in collection, R = length of partial ranking
  
  # initialize
  if(nrow(Pi) != nrow(X)){stop("Pi and X must have same number of rows")}
  I <- nrow(X)
  Q <- getQ(Pi,J)
  
  num_nodes <- 0
  curr_ranking <- c()
  while(length(curr_ranking)<(J-1)){
    S <- setdiff(1:J,curr_ranking)
    cost <- rep(NA,length(S))
    for(i in 1:length(S)){
      num_nodes <- num_nodes + 1
      s <- S[i]
      try_ranking <- c(curr_ranking,s)
      cost[i] <- totalcostheuristic_MB(Q,X,M,try_ranking,J,R)
    }
    curr_ranking <- c(curr_ranking,S[which.min(cost)])
  }
  curr_ranking <- c(curr_ranking,setdiff(1:J,curr_ranking))
  curr_ranking_cost <- dmbm(Pi,X,round(phat_conditional(X,M,curr_ranking) + seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(curr_ranking==j)}))],6),
                            curr_ranking,theta_conditional(mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],curr_ranking)}))),J,R),
                            M,log=T)
  
  diffs <- matrix(unlist(lapply(1:(J-1),function(j){
    new <- curr_ranking
    new[j:(j+1)] <- curr_ranking[(j+1):j]
    new
  })),byrow=T,ncol=J)
  costs <- apply(diffs,1,function(order){
    p <- phat_conditional(X,M,order) + 
      seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(order==j)}))] #additive term to correct ties
    theta <- theta_conditional(
      mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],order)}))),J,R)
    dmbm(Pi,X,round(p,6),order,theta,M,log=T)
  })
  num_nodes <- num_nodes + length(costs)
  olddiffs <- matrix(NA,nrow=0,ncol=J)
  
  while(any(costs>curr_ranking_cost)){
    olddiffs <- rbind(olddiffs,diffs)
    curr_ranking <- diffs[which.max(costs),]
    curr_ranking_cost <- max(costs)
    
    diffs <- matrix(unlist(lapply(1:(J-1),function(j){
      new <- curr_ranking
      new[j:(j+1)] <- curr_ranking[(j+1):j]
      new
    })),byrow=T,ncol=J)
    diffs <- unique(diffs[!apply(diffs,1,function(diff){any(apply(olddiffs,1,function(diff2){all(diff==diff2)}))}),]) #remove orders already tested
    costs <- apply(diffs,1,function(order){
      p <- phat_conditional(X,M,order) + 
        seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(order==j)}))] #additive term to correct ties
      theta <- theta_conditional(
        mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],order)}))),J,R)
      dmbm(Pi,X,round(p,6),order,theta,M,log=T)
    })
    num_nodes <- num_nodes + length(costs)
  }

  pi0_greedy <- c(curr_ranking,setdiff(1:J,curr_ranking))
  p_greedy <- phat_conditional(X,M,pi0_greedy) + 
    seq(0,1e-8,length=J)[unlist(lapply(1:J,function(j){which(pi0_greedy==j)}))] #additive term to correct ties
  theta_greedy <- theta_conditional(
    mean(unlist(apply(Pi,1,function(pi){d_R(pi[!is.na(pi)],pi0_greedy)}))),J,R)
  
  return(list(pi0_greedylocal = pi0_greedy,p_greedylocal=p_greedy,theta_greedylocal=theta_greedy,
              num_nodes = num_nodes))
}


