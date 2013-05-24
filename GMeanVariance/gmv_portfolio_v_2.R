library(PortfolioAnalytics)
library(PerformanceAnalytics)
library(DEoptim)
library(foreach)
library(iterators)
library(doParallel)
library(quantmod)

#symbol list stolen from deoptim vignette with some non yahoo symbols dropped

tickers = c( "VNO" , "VMC" , "WMT" , "WAG" , "DIS" , "WPO" , "WFC" , "WDC" ,
             "WY" , "WHR" , "WMB" , "WEC" , "XEL" , "XRX" , "XLNX" ,"ZION" ,"MMM" ,
             "ABT", "ADBE" , "AMD" , "AET" , "AFL" , "APD" , "ARG" ,"AA" , "AGN" ,
             "ALTR" , "MO" , "AEP" , "AXP" , "AIG" , "AMGN" , "APC" ,"ADI" , "AON" ,
             "APA", "AAPL" , "AMAT" ,"ADM" , "T" , "ADSK" , "ADP" , "AZO" , "AVY" ,
             "AVP", "BHI" , "BLL" , "BAC" , "BK" , "BCR" , "BAX" , "BBT" , "BDX" , "BMS" , "BBY" , 
             "BIG" , "HRB" , "BMC" , "BA" , "BMY" , "CA" , "COG" ,
             "CPB" , "CAH" , "CCL" , "CAT" , "CELG" , "CNP" , "CTL" , 
             "SCHW" , "CVX" , "CB" , "CI" ,"CINF" ,"CTAS" , "CSCO" , "C" , "CLF" ,
             "CLX", "CMS" , "KO" , "CCE" , "CL" , "CMCSA" ,"CMA" , "CSC" , "CAG" ,
             "COP" , "ED" , "GLW" , "COST" , "CVH" , "CSX" , "CMI" , "CVS" ,
             "DHR" , "DE")


getSymbols(tickers,from="2000-12-01", to = "2013-05-01")


## Make matrix of Returns code stolen from vignette of deoptim
P<- NULL
seltickers<-NULL
for(ticker in tickers){
  tmp = Cl(to.monthly(eval(parse(text=ticker))))
  if(is.null(P)){timeP=time(tmp)}
  if(any(time(tmp)!=timeP)) next
  else P = cbind(P,as.numeric(tmp))
  seltickers = c(seltickers,ticker)
}
P = xts(P,order.by=timeP)
colnames(P) = seltickers
R = diff(log(P))
R = R[-1,]



###### Calculate min variance portfolio

GMVconst = constraint(assets=colnames(R),min=rep(0.001,ncol(R)),max=rep(0.05,ncol(R)),min_sum=1,max_sum=1,
                      risk_aversion=1, 
                      weight_seq = generatesequence(by=.001))
GMVconst = add.objective(GMVconst,type="risk",name ="sd",enabled=TRUE,multiplier=0,risk_aversion=1)



GMVPortfolio <- optimize.portfolio.rebalancing(R, constraints=GMVconst,
                                               optimize_method="random", 
                                               trace=TRUE, 
                                               rebalance_on='months', 
                                               trailing_periods=NULL, 
                                               training_period=36,
                                               search_size=2000,
                                               verbose=FALSE,
                                               parallel=TRUE)