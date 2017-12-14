\name{optimPMG}
\alias{optimPMG}
\title{optimPMG}
\description{Estimates parameters of long-run and short-run relationships. Makes use of a "back-substitution" algorithm, as described by Pesaran, Shin and Smith (1999). Also estimates the information matrix as well as standard errors of estimations, as indicated in Equation 13 (Pesaran, Shin and Smith, 1999). Calculates Student's t-distribution type statistics, probs and confidence intervals. Also performs diagnostic tests of error terms, such as the Breusch-Godfrey autocorrelation test, the Goldfeld-Quandt heteroscedasticity test, the Conover nonparametric test of homogeneity of variance and the Jarque-Bera normality test}
\usage{optimPMG(dLL, maxIter, TetaStart, vecSR, vecLR, dataset, quantity, const)}
\arguments{
  \item{dLL}{a parameter indicating the convergence criterion; an optimization algorithm is stopped when an increase in concentrated log-likelihood function (Equation 8 in Pesaran, Shin and Smith (1999)) is less than dLL; the default value is  dLL=10^-10}
  \item{maxIter}{a maximum number of iterations; the default value is 200}
  \item{TetaStart}{a vector of first (initial) Teta values, from which the algorithm starts searching for parameters ensuring the maximization of log-likelihood function}
  \item{vecSR}{a list of vectors containing the column numbers of variables in short-run relationships for each group (alternatively a list of vectors containing the variables names instead of column numbers). In each vector of the list the first number must indicate dy (i.e., the dependant variable)}
  \item{vecLR}{a vector containing the column numbers of variables in long-run relationships (alternatively a vector containing the variables names instead of column numbers). The first number must indicate ly (i.e., the lagged dependant variable)}
  \item{dataset}{a panel data set in the form of stacked time series, containing variables of long-run and short-run relationships (i.e., including differentiated and lagged variables)}
  \item{quantity}{a vector of the number of time series observations in each group; in practice, it takes the form c(T1,...,Tn) since the PMG allows the numbers of time series observations to differ freely across groups (if the number of time series observations in each group is the same, then c(T,...,T) and T=T1=T2=...=Tn}
  \item{const}{logical. If TRUE (the default value), the intercept term is added to the model (i.e., to the short-run relationship)}
}
\details{Estimates parameters of long-run and short-run relationships. Also estimates the information matrix as well as standard errors of estimations, as indicated in Equation 13 (Pesaran, Shin and Smith, 1999). Calculates Student's t-distribution type statistics, probs and confidence intervals. Also performs diagnostic tests of error terms, such as the Breusch-Godfrey autocorrelation test, the Goldfeld-Quandt heteroscedasticity test and the Conover nonparametric test of homogeneity of variance and the Jarque-Bera normality test}
\value{
 \item{$LogL}{the concentrated log-likelihood function}
 \item{$dLogL}{the incresase of concentrated log-likelihood function in last iteration}
 \item{$i}{the number of iterations performed to achieve convergence}
 \item{$LR}{the estimated parameters of long-run relationships}
 \item{$SR}{the estimated parameters of short-run relationships}
 \item{$DiagTests}{results of diagnostic tests}
 \item{$residuals}{residuals}
}
\references{Pesaran, Shin and Smith (1999) <doi:10.1080/01621459.1999.10474156>}
\author{Lech Kujawski, Piotr Zientara}
\examples{
# first import DataExp, i=1...9, T1=T2=...T9=35
data(DataExp)
DataExp[1:5,]
# then prepare lags and diffs using LagPanel and DiffPanel
y10=data.frame(y10=DataExp[,1], row.names=row.names(DataExp))
cpi=data.frame(cpi=DataExp[,7], row.names=row.names(DataExp))
dy10=DiffPanel(variable=y10, quantity=rep(35,9))
dopeness=DiffPanel(variable=DataExp[,6], quantity=rep(35,9))
ly10=LagPanel(variable=y10, quantity=rep(35,9))
diip=DiffPanel(variable=DataExp[,11], quantity=rep(35,9))
dcrisk=DiffPanel(variable=DataExp[,9], quantity=rep(35,9))
ldcrisk=LagPanel(variable=dcrisk, quantity=rep(35,9))
dcpi=DiffPanel(variable=DataExp[,7], quantity=rep(35,9))
ddcpi=DiffPanel(variable=dcpi, quantity=rep(35,9))
ldebt=LagPanel(variable=DataExp[,4], quantity=rep(35,9))
# create homogenous preliminary dataset (containing "NA" as a result of DiffPanel, LagPanel)
dataPanel=cbind(y10, dy10, ly10, DataExp[,6], dopeness, diip,
DataExp[,11], ldcrisk, DataExp[,9], ddcpi, DataExp[,7])
dataPanel=data.frame(dataPanel)
names(dataPanel)=c("y10", "dy10", "ly10", "openess", "dopeness", "diip",
"iip", "ldcrisk", "crisk", "ddcpi", "cpi")
dataPanel[1:5,]
# prepare dataset and quantity for PMG or optimPMG functions using PanelNaOmit
dataPanel=PanelNaOmit(dataset=dataPanel, quantity=rep(35,9))
dataPanel$dataset[1:5,]
dataPanel$quantity
# optimPMG
OptimPmgExp=optimPMG(
dLL=10^-10,
maxIter=200,
TetaStart=rep(x=1, times=4), # note that length(TetaStart)=length(vecLR)-1
vecSR=list(SR1=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR2=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR3=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR4=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR5=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR6=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR7=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR8=c("dy10","dopeness","diip","ldcrisk","ddcpi"),
SR9=c("dy10","dopeness","diip","ldcrisk","ddcpi")),
vecLR=c("ly10","openess","iip","crisk","cpi"),
dataset=dataPanel$dataset,
quantity=dataPanel$quantity,
const=TRUE)
OptimPmgExp
}