\name{DiffPanel}
\alias{DiffPanel}
\title{DiffPanel}
\description{Calculates first differences of a particular variable from a panel data set}
\usage{DiffPanel(variable, quantity)}
\arguments{
  \item{variable}{a particular variable from a panel data set in the form of stacked time series; in practice; a selected singular column from a panel data set}
  \item{quantity}{a vector of the number of time series observations in each group; in practice, it takes the form c(T1,...Tn) since the PMG allows the numbers of time series observations to differ freely across groups (if the number of time series observations in each group is the same, then c(T,...,T) and T=T1=T2=...=Tn}
}
\details{Calculates first differences of a particular variable from a panel data set in order to bring it to stationarity. Preserves the original dimension of time series observations in each group, completing data lost due to differentiating by inserting "NA"}
\value{First differences of a particular variable from a panel data set}
\author{Lech Kujawski, Piotr Zientara}
\examples{
# first import DataExp, i=1...9, T1=T2=...T9=35
data(DataExp)
DataExp[1:5,]
# then execute DiffPanel
y10=data.frame(y10=DataExp[,1], row.names=row.names(DataExp))
dy10=DiffPanel(variable=y10, quantity=rep(35,9))
diip=DiffPanel(variable=DataExp[,11], quantity=rep(35,9))
cbind(y10,dy10,diip)[1:5,]
}