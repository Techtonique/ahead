% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{direct_sampling}
\alias{direct_sampling}
\title{Direct sampling}
\usage{
direct_sampling(
  data = NULL,
  n = 100L,
  method = c("kde", "surrogate", "bootstrap"),
  kde = NULL,
  seed = NULL,
  ...
)
}
\arguments{
\item{data}{A numeric vector or matrix.}

\item{n}{The number of samples to draw.}

\item{method}{The method to use for sampling.}

\item{kde}{The kernel density estimate to use for sampling.}

\item{seed}{The seed to use for sampling.}

\item{...}{Additional arguments to pass to the density function.}
}
\description{
Direct sampling
}
