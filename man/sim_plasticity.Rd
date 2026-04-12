\name{sim_plasticity}
\alias{sim_plasticity}
\title{Plasticity and Environmental Adaptation}
\description{
Demonstrates neural plasticity / colony adaptation correspondence,
with optional environmental shift.
}
\usage{
sim_plasticity(n_generations = 300, n_trails = 20,
               ltp_rate = 0.1, ltd_rate = 0.05,
               prune_threshold = 0.01, genesis_rate = 0.02,
               env_shift_at = NULL)
}
\arguments{
  \item{n_generations}{Number of generations.}
  \item{n_trails}{Number of trails.}
  \item{ltp_rate}{Long-term potentiation / reinforcement rate.}
  \item{ltd_rate}{Long-term depression / evaporation rate.}
  \item{prune_threshold}{Pruning threshold.}
  \item{genesis_rate}{New trail formation rate.}
  \item{env_shift_at}{Generation at which environment shifts (\code{NULL} = no shift).}
}
\value{
An S3 object of class \code{"plasticity_sim"} with \code{print()}
and \code{plot()} methods.
}
