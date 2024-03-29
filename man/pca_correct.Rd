% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalise.R
\name{pca_correct}
\alias{pca_correct}
\title{PCA-normalize}
\usage{
pca_correct(samples, num_comp, cumulative_variance)
}
\arguments{
\item{samples}{A data frame resulting from the \code{\link{normalize_reads}} function.}

\item{num_comp}{The number of PCA components to use in the normalization. If not specified, the function will use the cumulative variance to determine the number of components.}

\item{cumulative_variance}{The percentage of cumulative variance to be explained by the number of PCA components. If not specified, the function will use the \code{num_comp} parameter to determine the number of components.}
}
\value{
A normalized data frame.
}
\description{
This function performs PCA-based normalization on a normalized data frame that has undergone previous normalization steps such as read count normalization and GC-correction.
}
\references{
This function is based on a \href{https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com}{StackExchange answer}.
}
\seealso{
\code{\link{normalize_reads}}, \code{\link{gc_correct}}
}
