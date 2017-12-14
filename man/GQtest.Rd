\name{GQtest}
\alias{GQtest}
\title{GQtest}
\description{Tests for homoscedasticity}
\usage{GQtest(residuals, subsample, nep)}
\arguments{
  \item{residuals}{residuals for group i}
  \item{subsample}{the vector of c(s1i, s2i), where s1i=Ti/2 and s2i=s1i+1 (if Ti is EVEN) or where s1i=Ti/2-0.5 and s2i=s1i+1 (if Ti is UNEVEN)}
  \item{nep}{the number of estimated parameters for group i}
}
\details{Calculates F statistic with a prob}
\value{F statistic with a prob}
\author{Lech Kujawski, Piotr Zientara}
\examples{
# creating artificial variables
x1=rnorm(30,0,1)
x2=rnorm(30,0,1)
e=rnorm(30,0,0.2)
y=1+2*x1+3*x2+e
# any model
model=lm(y~x1+x2)
#BGtest
ExpGQtest=GQtest(residuals=resid(model), subsample=c(15,16), nep=3)
ExpGQtest
}