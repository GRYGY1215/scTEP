\name{fa}
\alias{fa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Learn latent representation for intersection between dataset and pathways.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Load pathways acorrding to the orgnism of dataset, intersect the gene expression matrix with pathway by gene symbols, apply factor analysis for each intersection.
}
\usage{
fa(data, genesets, data_org = "hsa")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
The data obtained using preprocessing function, a large list with exprression matrix, label, start cell type, all cell types.
}
  \item{genesets}{
%%     ~~Describe \code{genesets} here~~
List of pathways, consist of gene symbols.
}
  \item{data_org}{
%%     ~~Describe \code{data_org} here~~
The orgnism of the dataset. Human sapiens or Mus musculus.
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
A matrix confists of latent representation of intersection between gene expression and pathways.
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
data_fa = fa(data, genesets, data_org = 'mmu')

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

