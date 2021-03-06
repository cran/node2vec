% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/node2vecR.R
\name{node2vecR}
\alias{node2vecR}
\title{Algorithmic Framework for Representational Learning on Graphs}
\usage{
node2vecR(
  data,
  p = NULL,
  q = NULL,
  directed = NULL,
  num_walks = NULL,
  walk_length = NULL,
  dim = NULL
)
}
\arguments{
\item{data}{input data for edges consisting of at least two columns of data and if there are weights,it must be in the third column.}

\item{p}{return parameter.Default to 1.}

\item{q}{in-out parameter.Default to 1.}

\item{directed}{the network is directed or undirected.Default to undirected.}

\item{num_walks}{number of walks per node.Default to 10.}

\item{walk_length}{number of nodes in each walk.Default to 80.}

\item{dim}{embedding dimensions.Default to 128.}
}
\value{
embedding results for each node
}
\description{
Algorithmic Framework for Representational Learning on Graphs
}
\examples{
#Parameters can be customized as needed
data(gene_edges)
use_data<-gene_edges
emb<-node2vecR(use_data,p=2,q=1,num_walks=5,walk_length=5,dim=10)





}
