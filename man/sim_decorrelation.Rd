\name{sim_decorrelation}
\alias{sim_decorrelation}
\title{Decorrelation Parameter Sweep}
\description{
Sweeps the decorrelation parameter theta for both random forests and
ant colonies and returns empirical pairwise correlations.
}
\usage{
sim_decorrelation(theta_range = seq(0, 1, by = 0.05),
                  n_sites = 20, n_ants = 50, n_trees = 50,
                  n_features = 20, n_rep = 200)
}
\arguments{
  \item{theta_range}{Numeric vector of theta values.}
  \item{n_sites}{Number of candidate sites (ant colony).}
  \item{n_ants}{Number of ants.}
  \item{n_trees}{Number of trees.}
  \item{n_features}{Number of features (p).}
  \item{n_rep}{Monte Carlo replicates per theta.}
}
\value{
An S3 object of class \code{"decorrelation"} with \code{print()} and
\code{plot()} methods.
}
