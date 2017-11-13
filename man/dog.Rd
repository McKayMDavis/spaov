\name{dog}
\alias{dog}
\docType{data}
\title{Dog Diabetes Data
}
\description{
Data used to provide a simple example of spaov().
}
\usage{data("dog")}
\format{
  A data frame with 20 observations on the following 4 variables.
  \describe{
    \item{\code{Operations}}{a character vector}
    \item{\code{Dog}}{a numeric vector}
    \item{\code{Method}}{a character vector}
    \item{\code{LacticAcid}}{a numeric vector}
  }
}
\source{
Forbath, N., A. B. Kenshole, and G. Hetenyi, Jr. (1967), "Turnover of lactic acid in normal and diabetic dogs calculated by two tracer methods," \emph{Am. J. Physiol.}, v. 212,  pp. 1179-1183.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(dog)
str(dog)
}
