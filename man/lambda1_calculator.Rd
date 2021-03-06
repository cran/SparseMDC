% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lambda1_calculator.R
\name{lambda1_calculator}
\alias{lambda1_calculator}
\title{Lambda 1 Calcualtor}
\usage{
lambda1_calculator(dat_l, dim, nclust, nboot = 1000, alpha1 = 0.5,
  delta = 1e-07)
}
\arguments{
\item{dat_l}{list with D entries, each entry contains data d, p * n matrix.
This data should be centered and log-transformed.}

\item{dim}{Total number of conditions, D.}

\item{nclust}{Total number of clusters.}

\item{nboot}{The number of bootstrap repetitions used for estimating
lambda 1, the default value is 1000.}

\item{alpha1}{The quantile of the bootstrapped lambda 1 values to use,
range is (0,1). The default value is 0.5, the median of the calculated
lambda 1 values.}

\item{delta}{Small value term added to ensure existance, default value is
0.0000001.}
}
\value{
The estimated value of lambda1 for use in main SparseMDC
algorithm
}
\description{
Calculates the lambda 1 value for the SparseMDC algorithm. The lambda 1
value controls the number of marker genes selected for each
cluster in the output from SparseMDC. It is calculated as the value of
lambda 1 that results in no marker genes being selected when then are no
meaningful clusters present in the data. Please see the original manuscript
for further details.
}
\examples{
set.seed(10)
# Select small dataset for example
data_test <- data_biase[1:100,]
# Split data into condition A and B
data_A <- data_test[ , which(condition_biase == "A")]
data_B <- data_test[ , which(condition_biase == "B")]
data_C <- data_test[ , which(condition_biase == "C")]
# Store data as list
dat_l <- list(data_A, data_B, data_C)
# Pre-process the data
pdat <- pre_proc_data(dat_l, dim=3, norm = FALSE, log = TRUE,
center = TRUE)
lambda1 <- lambda1_calculator(pdat, dim = 3, nclust = 3 )

}
