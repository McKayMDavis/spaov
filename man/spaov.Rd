\name{spaov}
\alias{spaov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Perform Split Plot Design ANOVA
}
\description{
The function parses factor columns and block as factors and creates an ANOVA table for the split plot ANOVA design and creates requirement check diagnostic plots (residuals vs fitted and QQ of residuals).
}
\usage{
spaov(resp, nestfct, crossfct, block, plots = TRUE)
}
\arguments{
  \item{resp}{\code{resp} is the response variable of the ANOVA model.}
  \item{nestfct}{\code{nestfct} is the factor to nest the block of the ANOVA model.}
  \item{crossfct}{\code{crossfct} is the factor to be crossed with the nest factor and the block of the ANOVA model.}
  \item{block}{\code{block} is the block of the ANOVA design.}
  \item{plots}{\code{plots} defaults to \code{TRUE} and displays requirements plots. When \code{FALSE} plots do not display.}
}
\author{McKay Davis}
\examples{
data("dog")
View(dog)
spaov(resp = dog$LacticAcid, nestfct = dog$Method, crossfct = dog$Operations, block = dog$Dog, plots = TRUE)
}
