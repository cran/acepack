\name{ace}
\title{ace: Alternating Conditional Expectations}
\usage{
ace(x, y, wt, mon, lin, cat, circ, delrsq)

}
\arguments{
\item{x}{a matrix containing the independent variables.}
\item{y}{a vector containing the response variable.}
\item{wt}{an optional vector of weights.}
\item{mon}{an optional integer vector specifying which variables are to be
transformed by monotone transformations.
Positive values in \code{mon} refer to columns of the \code{x} matrix and zero
to the response variable.}
\item{lin}{an optional integer vector specifying which variables are to be
transformed by linear transformations.
Positive values in \code{lin} refer to columns of the \code{x} matrix and zero
to the response variable.}
\item{cat}{
an optional integer vector specifying which variables assume
categorical values.
Positive values in \code{cat} refer to columns of the \code{x}
matrix and zero to the response variable.}
\item{circ}{
an integer vector specifying which variables assume
circular (periodic) values.
Positive values in \code{circ} refer to columns of the \code{x}
matrix and zero to the response variable.}
\item{delrsq}{
termination threshold.  Iteration stops when R-squared
changes by less than \code{delrsq} in 3 consecutive iterations
(default 0.01).}
}

\bold{VALUE}
structure with the following components:
\item{x}{the input x matrix.}
\item{y}{the input y vector.}
\item{tx}{the transformed x values.}
\item{ty}{the transformed y values.}
\item{rsq}{the multiple R-squared value for the transformed values.}
\item{l}{not used in this version of ace}
\item{m}{not used in this version of ace}

\bold{DESCRIPTION }
Uses the alternating conditional expectations algorithm to find the transformations of y and x that maximise the proportion of variation in y explained by x.

\bold{REFERENCE}
 Breiman and Friedman, Journal of the American Statistical Association (September, 1985)
The R code is adapted from S code for avas() by Tibshirani, in the Statlib S archive; the FORTRAN is a double-precision version of FORTRAN code by Friedman and Spector in the Statlib general archive.

\examples{
TWOPI <- 8*atan(1)
x <- runif(200,0,TWOPI)
y <- exp(sin(x)+rnorm(200)/2)
a <- ace(x,y)
par(mfrow=c(3,1))
plot(a$y,a$ty)  # view the response transformation
plot(a$x,a$tx)  # view the carrier transformation
plot(a$tx,a$ty) # examine the linearity of the fitted model
}













