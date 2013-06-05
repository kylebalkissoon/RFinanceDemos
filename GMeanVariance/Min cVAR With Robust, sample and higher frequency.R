#' Conditional Value At Risk portfolio allocation example
#' 
#'
#' This is a demonstration of using DEOptim, a package for Differential Evolution Optimization,
#' to optimize the wights of a group of ETFs to satisfy several different criteria 
#' (Minimize Robust conditional Value at Risk cVaR [MinCVaR_Rob], minimize conditional VaR covariance 
#' on daily returns[ MinCVaR_cov_daily], minimize Standard Deviation covariance [MinSD_cov], and 
#' Maximize Sharpe Ratio covariance [MaxSR_cov])
#' 
#' Inputs: A list of ETF symbols, and the benchmark Symbol
#' 

# debugging
#options(warn=2)
#options(error=dump.frames)

library(PerformanceAnalytics)
library(DEoptim)
library(quantmod)
library(robust)
library(doParallel)

# make a Cluster of 8 processors to run the optimization in parallel

cl <- makeCluster(8)
registerDoParallel(cl)

#
etflist <- c("TLT","SPY","XLB","XLV","XLP","XLY","XLE","XLF","XLI","XLK","XLU")
getSymbols(etflist,from="2003-12-01", to = "2013-05-01")
benchlist <- "SPY"
tickers <- etflist 

## Make matrix of Returns code stolen from vignette of deoptim
P <- NULL
seltickers <- NULL
for(ticker in tickers){
  tmp <- Ad(to.monthly(eval(parse(text=ticker))))
  if(is.null(P)){ timeP <- time(tmp) }
  if(any(time(tmp)!=timeP)) next
  else P <- cbind(P,as.numeric(tmp))
  seltickers <- c(seltickers,ticker)
}

P <- xts(P,order.by=timeP)
colnames(P) <- seltickers
#
# TODO: use simple returns, not log
#
R <- diff(log(P))
R <- R[-1,]
returnsDiscrete <- CalculateReturns(P, method='discrete')
returnsLog <- CalculateReturns(P, method='log')

numperiods = nrow(P)
numreturns = nrow(R)
numinstruments = ncol(P)

chart.CumReturns(R, legend.loc="topleft", 
                 main="Cumulative Monthly Returns",
                 colorset=rich10equal )
chart.Boxplot(R)

## Make Matrix of daily returns
Pday <- NULL
seltickers <- NULL
for(ticker in tickers){
  tmp <- Ad(to.daily(eval(parse(text=ticker))))
  if(is.null(Pday)){ timeP <- time(tmp) }
  if(any(time(tmp)!=timeP)) next
  else Pday <- cbind(Pday,as.numeric(tmp))
  seltickers <- c(seltickers,ticker)
}
Pday <- xts(Pday,order.by=timeP)
colnames(Pday) <- seltickers
#
# TODO: use simple returns, not log
#
Rday <- diff(log(Pday))
Rday <- Rday[-1,]

#Benchmarks
P2 <- NULL
seltickers <- NULL
for(ticker in benchlist){
  tmp <- Ad(to.monthly(eval(parse(text=ticker))))
  if(is.null(P2)){ timeP <- time(tmp) }
  if(any(time(tmp)!=timeP)) next
  else P2 <- cbind(P2,as.numeric(tmp))
  seltickers <- c(seltickers,ticker)
}
P2 <- xts(P2,order.by=timeP)
colnames(P2) <- seltickers
#
# TODO: use simple returns, not log
#
R2 <- diff(log(P2))
R2 <- R2[-1,]

initw <- rep(1/ncol(R),ncol(R))
objectivefun <- function(w){
  if(sum(w)==0){
		w <- w + 1e-2 
	}
	w <- w / sum(w)
  targ <- ES(weights=w,method="gaussian",portfolio_method="component",mu=mu,sigma=sigma)
  tmp1 <- targ$ES
  tmp2 <- max(targ$pct_contrib_ES-0.05,0)
  out <- tmp1 + 1e3 * tmp2
  return(out)
}
objectivefunsd <- function(w){
  if(sum(w)==0){
    w <- w + 1e-2 
  }
  w <- w / sum(w)
  targ <- StdDev(R=rollR, weights=w, portfolio_method="component", sigma=sigma)
  out <- targ$StdDev
  return(out)
}
objectivefunsr <- function(w){
  if(sum(w)==0){
    w <- w + 1e-2 
  }
  w <- w / sum(w)
  targ <- SharpeRatio(R=rollR, weights=w, FUN='StdDev', Rf=0, p=0.95, annualize=FALSE)
  rtarg <- -targ
  return(rtarg)
}

source("random_portfolios.R")
source("constraints.R")


#N <- ncol(R)
N <- numinstruments
minw <- 0
maxw <- 1
lower <- rep(minw,N)
upper <- rep(maxw,N)

eps <- 0.025

# set by to .001 for more precision
weight_seq <- generatesequence(min=minw,max=maxw,by=.001,rounding=3)
rpconstraint <- constraint( assets=N, min_sum=1-eps, max_sum=1+eps, min=lower, max=upper, weight_seq=weight_seq)

rp <- random_portfolios(rpconstraints=rpconstraint,permutations=N*10)
rp <- rp/rowSums(rp)
controlDE <- list(reltol=0.00001, steptol=150, itermax=2000, 
                  trace=250, NP=as.numeric(nrow(rp)),
                  initialpop=rp,strategy=6, c=0)
# in real use, remove the set.seed() line to change the RNG start
set.seed(1234)

preturn <- R

for (p in 1:numinstruments){
	preturn[,p] <- 0
}
optweights <- R
cnames <- colnames(preturn)
zeros <- c(0)

for ( i in 2:numinstruments){
  
  zeros <- cbind(zeros, 0) 
  } # add a zero row for the current month to preturn
colnames(zeros) <- cnames
nextdate <- tail(time(optweights)+1/frequency(optweights), n=1)
newrow <- as.xts(zeros, order.by=c(nextdate))
optweights <- rbind(R, newrow) # add a next period 0 row
preturn <- rbind(preturn, newrow)

for (z in 1:numinstruments){
	optweights[,z] <- 0
}

##Sample Cov Monthly
#
for (i in 2:numreturns){
  rollR = first(R,i)
  mu = colMeans(rollR)
  sigma = cov(rollR)
  weightvec = DEoptim(fn=objectivefun, lower=lower, upper=upper,
                       control=controlDE)
  # weightvec$optim$bestmem is a list 
  preturn[i+1,] = weightvec$optim$bestmem*R[i]
  optweights[i,] = weightvec$optim$bestmem
}
optweights2 <- optweights/rowSums(optweights)
portreturn_cov <- optweights2*R
portreturn_cov <- rowSums(portreturn_cov)
portreturn_cov <- xts(portreturn_cov,order.by=index(R))

colnames(portreturn_cov) <- "MinCVaR_cov"
OOSweights_cov <- weightvec$optim$bestmem/sum(weightvec$optim$bestmem)

twoassets <- merge.xts(portreturn_cov,R2)

##Sample COV Daily
for (i in 2:numreturns){
  rollR <- first(R,i)
  mu <- colMeans(rollR)
  rend <- endpoints(Rday)
  rollRday <- first(Rday,rend[i]+38)
  sigma <- cov(rollRday)
  sigma <- sigma*sqrt(23)
  weightvec <- DEoptim(fn=objectivefun, lower=lower,
                       upper=upper, control=controlDE)
  preturn[i+1,] <- weightvec$optim$bestmem*R[i]
  optweights[i+1,] <- weightvec$optim$bestmem
}
optweights2 <- optweights/rowSums(optweights)
portreturn_covd <- optweights2*R
portreturn_covd <- rowSums(portreturn_covd)
portreturn_covd <- xts(portreturn_covd, order.by=index(R))

colnames(portreturn_covd) <- "MinCVaR_cov_daily"
OOSweights_cov_d <- weightvec$optim$bestmem/sum(weightvec$optim$bestmem)

colofassets <- merge.xts(portreturn_covd, portreturn_cov, R2)

###ROBUST COVARIANCE 
for (i in 2:numreturns){
  rollR <- first(R,i)
  mu <- colMeans(rollR)
#rend = endpoints(Rday)
#rollRday = first(Rday,rend[i]+38)
  sigma <- covRob(rollRday)
  sigma <- sigma[[2]]
#sigma = sigma*sqrt(23)
  weightvec <- DEoptim(fn=objectivefun, lower=lower,
                       upper=upper, control=controlDE)  
  preturn[i+1,] <- weightvec$optim$bestmem*R[i]
  optweights[i+1,] <- weightvec$optim$bestmem
}

optweights2 <- optweights/rowSums(optweights)
portreturn_rob <- optweights2*R
portreturn_rob <- rowSums(portreturn_rob)
portreturn_rob <- xts(portreturn_rob,order.by=index(R))
colnames(portreturn_rob) <- "MinCVaR_Rob" 
colofassets <- merge.xts(colofassets,portreturn_rob)

for ( i in 2:numreturns){
  rollR <- first(R, i)
  mu <- colMeans(rollR)
  sigma <- cov(rollR)
  weightvec <- DEoptim(fn=objectivefunsd, lower=lower,
                       upper=upper, control=controlDE)
  preturn[i+1,] <- weightvec$optim$bestmem*R[i]
  optweights[i+1,] <- weightvec$optim$bestmem
}

optweights2 <- optweights/rowSums(optweights)
portreturn_minsd_cov <- optweights2*R
portreturn_minsd_cov <- rowSums(portreturn_minsd_cov)
portreturn_minsd_cov <- xts(portreturn_minsd_cov, order.by=index(R))

colnames(portreturn_minsd_cov) <- "MinSD_cov"
OOSweights_minsd_cov <- weightvec$optim$bestmem/sum(weightvec$optim$bestmem)

blob=merge.xts(portreturn_minsd_cov,R2)

### Max Sharpe

for (i in 2:numreturns){ 
  rollR <- first(R, i)
  mu <- colMeans(rollR)
  sigma <- cov(rollR)
  weightvec <- DEoptim(fn=objectivefunsr, lower=lower,
                       upper=upper, control=controlDE)
  preturn[i+1,] <- weightvec$optim$bestmem*R[i]
  optweights[i+1,] <- weightvec$optim$bestmem
}
portreturn_maxsr_cov <- optweights2*R
portreturn_maxsr_cov <- rowSums(portreturn_maxsr_cov)
portreturn_maxsr_cov <- xts(portreturn_maxsr_cov, order.by=index(R))

colnames(portreturn_maxsr_cov) <- "MaxSR_cov"
OOSweights_maxsr_cov <- weightvec$optim$bestmem/sum(weightvec$optim$bestmem)

colofassets <- merge.xts(portreturn_minsd_cov, 
                         R2, 
                         portreturn_rob, 
                         portreturn_covd, 
                         portreturn_maxsr_cov)
#colofassets2 <- colofassets['2007-12-31/2013-05-01']

chart.CumReturns(colofassets, wealth.index=TRUE, legend.loc="topleft")
