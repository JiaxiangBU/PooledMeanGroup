\name{PanelNaOmit}
\alias{PanelNaOmit}
\title{PanelNaOmit}
\description{Prepares a panel data set for further calculations by eliminating "NA" and modifying quantity or a vector of the number of time series observations in each group}
\usage{PanelNaOmit(dataset, quantity)}
\arguments{
  \item{dataset}{a panel data set in the form of stacked time series, containing variables of long-run and short-run relationships (i.e., including differentiated and lagged variables from DiffPanel or LagPanel)}
  \item{quantity}{a vector of the number of time series observations in each group; in practice, it takes the form c(T1,...Tn) since the PMG allows the numbers of time series observations to differ freely across groups (if the number of time series observations in each group is the same, then c(T,...,T) and T=T1=T2=...=Tn}
}
\details{Eliminates "NA" and modifies quantity or a vector of the number of time series observations in each group}
\value{
 \item{$dataset}{panel data set for further calculations modified by eliminating "NA"}
 \item{$quantity}{modified vector of the number of time series observations in each group}
}
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
# create homogenous preliminary dataset (containing "NA") after DiffPanel, LagPanel
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
}