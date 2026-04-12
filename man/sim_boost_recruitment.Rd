\name{sim_boost_recruitment}
\alias{sim_boost_recruitment}
\title{Boosting and Adaptive Recruitment Simulation}
\description{
Runs AdaBoost and ant colony adaptive recruitment in parallel.
}
\usage{
sim_boost_recruitment(n_iterations = 100, n_ants = 50,
                      n_sites = 10, evaporation_rate = 0.1,
                      noise_level = 0.15, n_rep = 100)
}
\arguments{
  \item{n_iterations}{Number of boosting iterations / recruitment waves.}
  \item{n_ants}{Number of ants.}
  \item{n_sites}{Number of candidate sites.}
  \item{evaporation_rate}{Pheromone evaporation rate.}
  \item{noise_level}{Noise level for data generation.}
  \item{n_rep}{Monte Carlo replicates for convergence curves.}
}
\value{
An S3 object of class \code{"boost_recruitment"} with elements
\code{$boost} and \code{$colony} (numeric learning curves),
and \code{print()} and \code{plot()} methods.
}
