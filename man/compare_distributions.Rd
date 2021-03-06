% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{compare_distributions}
\alias{compare_distributions}
\title{Compare Distributions}
\usage{
compare_distributions(data, limits = NULL, quantity = "Data",
  title = NULL)
}
\arguments{
\item{data}{data.frame or matrix; The first column is the empirical data, the others being the data from the theoretical distribution e.g. gene means or library sizes}

\item{quantity}{string; Name of thing being compared - used in the axis labels}

\item{limits;}{passed to scale_x_log10 "A numeric vector of length two providing limits of the scale.". Defaults to empirical range of the data.}
}
\value{
A ggplot2 object of pdf, cdf, qq, and pp plots.
}
\description{
\code{compare_distributions} Compare an empirical data distribution to fitted distributions using pdf, cdf, qq, and pp plots.
}
\examples{
fake_data <- data.frame(data1 = rlnorm(n = 500, meanlog = 0, sdlog = 1), # would usually empirical data
                        data2 = rlnorm(n = 500, meanlog = 0, sdlog = 1.5),
                        data3 = rlnorm(n = 500, meanlog = 0.5, sdlog = 1))
compare_distributions(fake_data)

}
