\name{BGtest}
\alias{BGtest}
\title{BGtest}
\description{Tests autocorrelation between the current and lagged residuals. The test is a joint test of the first P autocorrelations}
\usage{BGtest(residuals, explvariab, acor.ord)}
\arguments{
  \item{residuals}{residuals for group i}
  \item{explvariab}{explanatory variables (regressors) for group i}
  \item{acor.ord}{order of tested autocorrelations}
}
\details{Calculates statistics and probs of the Breusch-Godfrey autocorrelation test (with two variants: chi-squared and F}
\value{Chi-squared and F statistics with probs}
\author{Lech Kujawski, Piotr Zientara}
\examples{
# creating artificial variables
x1=rnorm(30,0,1)
x2=rnorm(30,0,1)
e=rnorm(30,0,0.2)
y=1+2*x1+3*x2+e
# any model
model=lm(y~x1+x2)
# BGtest
ExpBGtest=BGtest(residuals=resid(model), explvariab=cbind(x1,x2), acor.ord=4)
ExpBGtest
}