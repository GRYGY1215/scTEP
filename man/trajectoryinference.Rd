\name{trajectoryinference}
\alias{trajectoryinference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
conduct trajectory inferecnce
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Trajectory inference. First using the clustering results to calculate the pseudotime. Second, using MST algorithm to generate a MST. Third, sort the clusters in the MST according to the average pseudotime of outputs of first step. Then, we obtained a graph indicating the trajectory of dataset.
}
\usage{
trajectoryinference(data, data_org, fa, allCluster, ncores = 10L, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
the data contains expr
}
 \item{data_org}{
%%     ~~Describe \code{x} here~~
the orgnism of dataset
}
\item{fa}{
%%     ~~Describe \code{x} here~~
latent representation using factor analysis of intersect between data and pathway.
}
\item{allCluster}{
%%     ~~Describe \code{x} here~~
clustering results with k from 5 to 10
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
A list consist of follows.
(1) pseudotime: the pseudotime of cells.
(2) cluster: The clustering results.
(3) data_clus_cent: The center of each clusters, in reduced dimensional space.
(4) milestone_network: The milestone networks of the output graph.
(5) sorted_g: The igraph object of the output graph.
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
library(igraph)
data('genesets')
data('goolam')

data = preprocessing(goolam)
data_fa = fa(data, genesets, data_org = 'mmu')
allCluster = clustering(data)
out = trajectoryinference(data, 'mmu', data_fa, allCluster)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

