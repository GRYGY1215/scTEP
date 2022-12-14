\name{clustering}
\alias{clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Using scDHA to conduct clustering with k from 5 to 10
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function conduct robust clustering using scDHA with k, number of cluster, set from 5 to 10.
}
\usage{
clustering(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
The data obtained using preprocessing function, a large list with exprression matrix, label, start cell type, all cell types.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
List with 6 clustering results with k, number of cluster, set from 5 to 10.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(dplyr)
library(foreach)
library(scDHA)
library(scTEP)
data('genesets')
data('goolam')

data = preprocessing(goolam)
data_fa = fa(data, genesets, data_org = 'mmu')
allCluster = clustering(data)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
