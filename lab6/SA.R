load("bankdata.Rdata")
nclients <- dim(bankdata)[1]
plot(bankdata,main = "plot of data with result of combination4")

set.seed(12345)
first22_index <- sample(nrow(bankdata),22,replace = FALSE)
fullset <- 1:nclients
points(bankdata[first22_index,],col="red",pch=16)

crit <- function(dat, subs){
  s <- length(subs)
  dist <- matrix(rep(NA, nclients*s), ncol=s)
  for (i in 1:s){
    dist[, i] <- sqrt((dat[,1]-dat[subs[i],1])^2+(dat[,2]-dat[subs[i],2])^2)
  }
  sum(apply(dist, 1, min))
}

test <- crit(bankdata,first22_index)


proposal_distribution <- function(xt){
  current_index <- xt
  remove <- sample(length(current_index), 11)
  new <- sample(setdiff(fullset,current_index),11)
  new_index <- c(current_index[-remove],new)
  return(new_index)
}

simulated_annealing <- function(crit,start_index,initial_temp,stage,mj,schedule){
  xt <- start_index
  temp_j <- initial_temp
  values <- numeric(stage+1)
  for (j in 0:stage) {
    for (rep in c(1:3)) {
      for (i in 1:mj) {
        candidate <- proposal_distribution(xt)
        h <- exp((crit(bankdata,xt) - crit(bankdata,candidate))/temp_j)
        prob <- min(h,1)
        xt1 <- xt
        if(prob == 1 || runif(1) < prob){
          xt1 <- candidate
        }
        xt <- xt1
      }
    }
    if(schedule=="Geometric"){
      temp_j <- alpha*temp_j
    }
    else if(schedule=="Logarithmic"){
      temp_j <- temp_j/(1+log(1+j))
    }
    else if(schedule=="Linear"){
      temp_j <- max(temp_j - j * c, 1)
    }
    mj <- round(beta * mj)
    values[j+1] <- crit(bankdata, xt)
    plot(values,type="l")
  }
  return(xt)
}

alpha <- 0.9
beta <- 1.2
c <- 5


