\name{sim_colony_convergence}
\alias{sim_colony_convergence}
\title{Colony Convergence Simulation}
\description{
Tracks the colony's probability estimates over time, demonstrating
convergence to the best site via Thompson sampling and recruitment.
}
\usage{
sim_colony_convergence(n_ants = 50, n_sites = 5,
                       site_qualities = c(0.9, 0.7, 0.5, 0.3, 0.1),
                       n_steps = 200, p_explore = 0.3)
}
\arguments{
  \item{n_ants}{Number of ants.}
  \item{n_sites}{Number of candidate sites.}
  \item{site_qualities}{True quality vector.}
  \item{n_steps}{Number of time steps.}
  \item{p_explore}{Exploration probability.}
}
\value{
An S3 object of class \code{"colony_convergence"} with \code{print()}
and \code{plot()} methods.
}
