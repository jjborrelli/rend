% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_modelFUNCTIONS.R
\name{netHTML}
\alias{netHTML}
\title{Dynamic network visualization}
\usage{
netHTML(mat, dyn, path1 = getwd())
}
\arguments{
\item{mat}{initial adjacency matrix of the food web}

\item{dyn}{output of \code{CRsimulator}}

\item{path1}{where to save the html plot}
}
\value{
An html file of the changes to the food web during the simulation.
}

