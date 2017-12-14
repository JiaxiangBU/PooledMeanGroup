\name{JBtest}
\alias{JBtest}
\title{JBtest}
\description{Tests for normality}
\usage{JBtest(residuals)}
\arguments{
  \item{residuals}{residuals for group i}
}
\details{Calculates chi-squared statistic with a prob}
\value{Chi-squared statistic with a prob}
\author{Lech Kujawski, Piotr Zientara}
\examples{
# creating artificial variables
x1=rnorm(30,0,1)
x2=rnorm(30,0,1)
e=rnorm(30,0,0.2)
y=1+2*x1+3*x2+e
# any model
model=lm(y~x1+x2)
#JBtest
ExpJBtest=JBtest(residuals=resid(model))
ExpJBtest
}