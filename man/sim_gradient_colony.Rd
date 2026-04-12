\name{sim_gradient_colony}
\alias{sim_gradient_colony}
\title{Gradient Descent and Generational Colony Learning}
\description{
Simulates generational pheromone evolution alongside neural network
training via SGD.
}
\usage{
sim_gradient_colony(n_generations = 50, n_ants = 30,
                    n_trails = 10, evaporation_rate = 0.05,
                    learning_rate = 0.05,
                    hidden_units = c(16, 8),
                    task = "classification", n_rep = 50)
}
\arguments{
  \item{n_generations}{Number of generations / epochs.}
  \item{n_ants}{Number of ants per generation.}
  \item{n_trails}{Number of trails / sites.}
  \item{evaporation_rate}{Pheromone evaporation rate.}
  \item{learning_rate}{Neural network learning rate.}
  \item{hidden_units}{Integer vector of hidden layer sizes.}
  \item{task}{\code{"classification"} or \code{"regression"}.}
  \item{n_rep}{Number of replicates for learning curves.}
}
\value{
An S3 object of class \code{"gradient_colony"} with elements
\code{$colony} and \code{$neural_net} (numeric vectors of normalized
performance), and \code{print()}, \code{summary()}, and \code{plot()} methods.
}
