% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{SecAct.heat.map}
\alias{SecAct.heat.map}
\title{Draw a heat map}
\usage{
SecAct.heat.map(
  fg.mat,
  title = NULL,
  colors = c("#03c383", "#aad962", "#fbbf45", "#ef6a32")
)
}
\arguments{
\item{fg.mat}{A matrix of values.}

\item{title}{The title for plot.}

\item{colors}{Colors.}
}
\value{
A ggplot2 object.
}
\description{
Draw a heat map of secreted proteins.
}
