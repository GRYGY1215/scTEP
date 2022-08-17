\name{preprocessing}
\alias{preprocessing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Preprocessing the dataset.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Conduct preprocessing, including remove all zero columns and scale gene expression smaller than 100 by log transformation with 2 as base.
}
\usage{
preprocessing(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
A large SingleCellExperiment object contains the gene expression file, cell type, label.
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
A list contains the gene expression file, start cell type, cell type, label.
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
data('genesets')
data('goolam')

data = preprocessing(goolam)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

