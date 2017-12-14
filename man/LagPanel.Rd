\name{LagPanel}
\alias{LagPanel}
\title{LagPanel}
\description{Provides the first lag of a particular variable from a panel data set}
\usage{LagPanel(variable, quantity)}
\arguments{
  \item{variable}{a particular variable from a panel data set in the form of stacked time series; in practice; a selected singular column from a panel data set}
  \item{quantity}{a vector of the number of time series observations in each group; in practice, it takes the form c(T1,...Tn) since the PMG allows the numbers of time series observations to differ freely across groups (if the number of time series observations in each group is the same, then c(T,...,T) and T=T1=T2=...=Tn}
}
\details{Provides the first lag of a particular variable from a panel data set. Preserves the original dimension of time series observations in each group, completing data lost due to lagging by inserting "NA"}
\value{A lagged particular variable from a panel data set}
\author{Lech Kujawski, Piotr Zientara}
\examples{
# first import DataExp, i=1...9, T1=T2=...T9=35
data(DataExp)
DataExp[1:5,]
# then execute LagPanel
y10=data.frame(y10=DataExp[,1], row.names=row.names(DataExp))
ly10=LagPanel(variable=y10, quantity=rep(35,9))
ldebt=LagPanel(variable=DataExp[,4], quantity=rep(35,9))
cbind(y10,ly10,ldebt)[1:5,]
}