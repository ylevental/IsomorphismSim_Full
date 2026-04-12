\name{sim_variance_decomp}
\alias{sim_variance_decomp}
\title{Variance Decomposition Simulation}
\description{
Simulates both an ant colony and a random forest across a range of
ensemble sizes and correlation levels, validating the variance
decomposition formula.
}
\usage{
sim_variance_decomp(N_range = c(5, 10, 25, 50, 100),
                    rho_range = seq(0, 1, by = 0.2),
                    sigma2 = 1, n_rep = 500)
}
\arguments{
  \item{N_range}{Integer vector of ensemble sizes.}
  \item{rho_range}{Numeric vector of pairwise correlation levels.}
  \item{sigma2}{Individual unit variance.}
  \item{n_rep}{Number of Monte Carlo replicates.}
}
\value{
An S3 object of class \code{"variance_decomp"} with \code{print()},
\code{summary()}, and \code{plot()} methods.
}
