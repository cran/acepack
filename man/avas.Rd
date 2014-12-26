\name{avas}
\alias{avas}
\title{Additivity and variance stabilization for regression}
\usage{
avas(x, y, wt, mon, lin, cat, circ, delrsq,yspan)
}
\arguments{
    \item{x}{a matrix containing the independent variables.}
    \item{y}{a vector containing the response variable.}
    \item{wt}{an optional vector of weights.}
    \item{mon}{an optional integer vector specifying which variables are
	to be transformed by monotone transformations.  Positive values
	in \code{mon} refer to columns of the \code{x} matrix and zero
	to the response variable.}
    \item{lin}{an optional integer vector specifying which variables are
	to be transformed by linear transformations.  Positive values in
	\code{lin} refer to columns of the \code{x} matrix and zero to
	the response variable.}
    \item{cat}{an optional integer vector specifying which variables
	assume categorical values.  Positive values in \code{cat} refer
	to columns of the \code{x} matrix and zero to the response
	variable.}
    \item{circ}{an integer vector specifying which variables assume
	circular (periodic) values.  Positive values in \code{circ}
	refer to columns of the \code{x} matrix and zero to the response
	variable.}
    \item{delrsq}{termination threshold.  Iteration stops when R-squared
	changes by less than \code{delrsq} in 3 consecutive iterations
	(default 0.01).}
    \item{yspan}{Optional window size parameter for smoothing the
	variance.  Range is \eqn{[0,1]}.  Default is 0 (cross validated
	choice). .5 is a reasonable alternative to try.}
}
\value{
    A structure with the following components:
    \item{x}{the input x matrix.}
    \item{y}{the input y vector.}
    \item{tx}{the transformed x values.}
    \item{ty}{the transformed y values.}
    \item{rsq}{the multiple R-squared value for the transformed values.}
    \item{l}{not used in this version of avas}
    \item{m}{not used in this version of avas}
    \item{yspan}{span used for smoothing the variance}
    \item{iters}{iteration number and rsq for that iteration}
    \item{niters}{number of iterations used}
}
\references{
    Rob Tibshirani (1987),
    ``Estimating optimal transformations for regression''. 
    \emph{Journal of the American Statistical Association} \bold{83},
    394ff.
}
\examples{
TWOPI <- 8*atan(1)
x <- runif(200,0,TWOPI)
y <- exp(sin(x)+rnorm(200)/2)
a <- avas(x,y)
par(mfrow=c(3,1))
plot(a$y,a$ty)  # view the response transformation
plot(a$x,a$tx)  # view the carrier transformation
plot(a$tx,a$ty) # examine the linearity of the fitted model
}
\emph{File automatically converted from S(-PLUS) help format}
