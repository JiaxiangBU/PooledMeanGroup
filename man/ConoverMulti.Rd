\name{ConoverMulti}
\alias{ConoverMulti}
\title{ConoverMulti}
\description{Tests for homoscedasticity among subsamples k within a particular group i (note that the Conover test is a non-parametric test)}
\usage{ConoverMulti(residuals, subsample)}
\arguments{
  \item{residuals}{residuals for group i}
  \item{subsample}{the vector of c(s1i, s2i,..., ski), where s1i+s2i+...+ski=Ti (i.e., the vector divides a particular group i into subsamples)}
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
# ConoverMulti
ExpConoverMulti=ConoverMulti(residuals=resid(model), subsample=c(10,10,10))
ExpConoverMulti
}