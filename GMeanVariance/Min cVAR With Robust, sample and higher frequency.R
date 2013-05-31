#' Conditional Value At Risk portfolio allocation example
#' 
#' ES(mu=colMeans(edhec),sigma=cov(edhec),m3=PerformanceAnalytics:::M3.MM(edhec),m4=PerformanceAnalytics:::M4.MM(edhec), weights=rep(1/ncol(edhec),ncol(edhec)))
#' 

# debugging
options(warn=2)
options(error=dump.frames)

library(PerformanceAnalytics)
library(DEoptim)
library(quantmod)
library(robust)
library(doParallel)

# make a Cluster of 8 processors to run the optimization in parallel

cl <- makeCluster(8)
registerDoParallel(cl)

#
etflist = c("TLT","SPY","XLB","XLV","XLP","XLY","XLE","XLF","XLI","XLK","XLU")
getSymbols(etflist,from="2003-12-01", to = "2013-05-01")
benchlist = "SPY"
tickers = etflist 

## Make matrix of Returns code stolen from vignette of deoptim
P<- NULL
seltickers<-NULL
for(ticker in tickers){
  tmp = Ad(to.monthly(eval(parse(text=ticker))))
  if(is.null(P)){timeP=time(tmp)}
  if(any(time(tmp)!=timeP)) next
  else P = cbind(P,as.numeric(tmp))
  seltickers = c(seltickers,ticker)
}

P = xts(P,order.by=timeP)
colnames(P) = seltickers
R = diff(log(P))
R = R[-1,]

chart.CumReturns(R, legend.loc="topleft", 
                 main="Cumulative Monthly Returns",
                 colorset=rich10equal )
chart.Boxplot(R)

## Make Matrix of daily returns
Pday<- NULL
seltickers<-NULL
for(ticker in tickers){
  tmp = Ad(to.daily(eval(parse(text=ticker))))
  if(is.null(Pday)){timeP=time(tmp)}
  if(any(time(tmp)!=timeP)) next
  else Pday = cbind(Pday,as.numeric(tmp))
  seltickers = c(seltickers,ticker)
}
Pday = xts(Pday,order.by=timeP)
colnames(Pday) = seltickers
Rday = diff(log(Pday))
Rday = Rday[-1,]

#Benchmarks
P2<- NULL
seltickers<-NULL
for(ticker in benchlist){
  tmp = Ad(to.monthly(eval(parse(text=ticker))))
  if(is.null(P2)){timeP=time(tmp)}
  if(any(time(tmp)!=timeP)) next
  else P2 = cbind(P2,as.numeric(tmp))
  seltickers = c(seltickers,ticker)
}
P2 = xts(P2,order.by=timeP)
colnames(P2) = seltickers
R2 = diff(log(P2))
R2 = R2[-1,]

initw = rep(1/ncol(R),ncol(R))
objectivefun = function(w){
  if(sum(w)==0){
		w = w + 1e-2 
	}
	w = w / sum(w)
  targ = ES(weights=w,method="gaussian",portfolio_method="component",mu=mu,sigma=sigma)
  tmp1 = targ$ES
  tmp2 = max(targ$pct_contrib_ES-0.05,0)
  out = tmp1 + 1e3 * tmp2
  return(out)
}

source("random_portfolios.R")
source("constraints.R")


N = ncol(R)
minw = 0
maxw = 1
lower = rep(minw,N)
upper = rep(maxw,N)

eps <- 0.025
# the weights are in tenths for testing
# set by to .001 for more precision
weight_seq<-generatesequence(min=minw,max=maxw,by=.1,rounding=3)
rpconstraint<-constraint( assets=N, min_sum=1-eps, max_sum=1+eps, min=lower, max=upper, weight_seq=weight_seq)

rp = random_portfolios(rpconstraints=rpconstraint,permutations=N*10)
rp <-rp/rowSums(rp)
controlDE = list(reltol=0.00001,steptol=150,itermax=2000,trace=250,NP=as.numeric(nrow(rp)),initialpop=rp,strategy=6,c=0)
set.seed(1234)

preturn = R
for (p in 1:ncol(preturn)){
	preturn[,p] = 0
}
optweights = R
cnames <- colnames(preturn)
zeros <- c(0)
for ( i in 2:ncol(optweights)){
  zeros <- cbind(zeros, 0) 
  }
colnames(zeros) <- cnames
nextdate <- tail(time(optweights)+1/frequency(optweights), n=1)
newrow <- as.xts(zeros, order.by=c(nextdate))
optweights=rbind(R, newrow)

for (z in 1:ncol(optweights)){
	optweights[,z] = 0
}

##Sample Cov Monthly
for (i in 2:length(R)){
  rollR = first(R,i)
  mu = colMeans(rollR)
  sigma = cov(rollR)
  weightvec = DEoptim(fn=objectivefun,lower=lower,upper=upper,control=controlDE)
  preturn[i+1,] = weightvec$optim$bestmem*R[i+1]
  optweights[i+1,] = weightvec$optim$bestmem
}
optweights2 = optweights/rowSums(optweights)
portreturn_cov = optweights2*R
portreturn_cov = rowSums(portreturn_cov)
portreturn_cov = xts(portreturn_cov,order.by=index(R))

colnames(portreturn_cov) = "MinCVaR_cov"
OOSweights_cov = weightvec$optim$bestmem/sum(weightvec$optim$bestmem)

twoassets = merge.xts(portreturn_cov,R2)

##Sample COV Daily

for (i in 2:length(R)){
  rollR = first(R,i)
  mu = colMeans(rollR)
  rend = endpoints(Rday)
  rollRday = first(Rday,rend[i]+38)
  sigma = cov(rollRday)
  sigma = sigma*sqrt(23)
  weightvec = DEoptim(fn=objectivefun,lower=lower,upper=upper,control=controlDE)
  preturn[i+1,] = weightvec$optim$bestmem*R[i+1]
  optweights[i+1,] = weightvec$optim$bestmem
}
optweights2 = optweights/rowSums(optweights)
portreturn_covd = optweights2*R
portreturn_covd = rowSums(portreturn_covd)
portreturn_covd = xts(portreturn_covd,order.by=index(R))

colnames(portreturn_covd) = "MinCVaR_cov_daily"
OOSweights_cov_d = weightvec$optim$bestmem/sum(weightvec$optim$bestmem)

colofassets = merge.xts(portreturn_covd,portreturn_cov,R2)

###ROBUST COVARIANCE 
for (i in 12:length(R)){
  rollR = first(R,i)
  mu = colMeans(rollR)
#rend = endpoints(Rday)
#rollRday = first(Rday,rend[i]+38)
  sigma =covRob(rollRday)
  sigma = sigma[[2]]
#sigma = sigma*sqrt(23)
  weightvec = DEoptim(fn=objectivefun,lower=lower,upper=upper,control=controlDE)  
  preturn[i+1,] = weightvec$optim$bestmem*R[i+1]
  optweights[i+1,] = weightvec$optim$bestmem
}

optweights2 = optweights/rowSums(optweights)
portreturn_rob = optweights2*R
portreturn_rob = rowSums(portreturn_rob)
portreturn_rob = xts(portreturn_rob,order.by=index(R))
colnames(portreturn_rob) = "MinCVaR_Rob" 
colofassets = merge.xts(colofassets,portreturn_rob)

chart.CumReturns(colofassets,wealth.index=TRUE,legend.loc="topleft")
