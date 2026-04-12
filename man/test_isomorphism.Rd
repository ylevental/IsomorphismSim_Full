\name{test_isomorphism}
\alias{test_isomorphism}
\title{Test Isomorphism Between Two Learning Curves}
\description{
Compares two numeric vectors (e.g. colony fitness and NN loss) using
a Kolmogorov-Smirnov test, correlation analysis, or dynamic time warping.
}
\usage{
test_isomorphism(ant_result, ml_result,
                 method = c("ks", "correlation", "dtw"))
}
\arguments{
  \item{ant_result}{Numeric vector (colony learning curve).}
  \item{ml_result}{Numeric vector (ML learning curve).}
  \item{method}{One of \code{"ks"}, \code{"correlation"}, or \code{"dtw"}.}
}
\value{
An S3 object of class \code{"iso_test"} with \code{print()} and
\code{summary()} methods.
}
