\name{sim_margin_analysis}
\alias{sim_margin_analysis}
\title{Margin Analysis}
\description{
Computes margin distributions for boosting and ant colony quorum decisions.
}
\usage{
sim_margin_analysis(n_iterations = 200, n_ants = 100,
                    noise_level = 0.2)
}
\arguments{
  \item{n_iterations}{Number of boosting iterations.}
  \item{n_ants}{Number of ants.}
  \item{noise_level}{Noise level.}
}
\value{
An S3 object of class \code{"margin_analysis"} with \code{print()}
and \code{plot()} methods.
}
