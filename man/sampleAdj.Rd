% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleAdj.R
\name{sampleAdj}
\alias{sampleAdj}
\title{Sample clusters from adjacence matrix}
\usage{
sampleAdj(
  adj,
  size = 5,
  n = 5,
  maxTry = 1000,
  starTopo = 1,
  type = 1,
  seed = -1
)
}
\arguments{
\item{adj}{Adjazenzmatrix}

\item{size}{Size of clusters (number of links)}

\item{n}{Number of clusters}

\item{maxTry}{numbers of random tries to establish clusters}

\item{starTopo}{create star topology}

\item{type}{If set to 0, n refers to number of edges, choose
1 for number of nodes}

\item{seed}{if set to -1, no seed is used.}
}
\description{
Desc
}
