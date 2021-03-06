% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penPHcure.object.R
\name{penPHcure.object}
\alias{penPHcure.object}
\title{Penalized PH cure model object}
\arguments{
\item{AIC}{a list with elements containing the results of the selected model based on the Akaike information criterion (AIC). See Details.}

\item{BIC}{a list with elements containing the results of the selected model based on the Bayesian Information Criterion (BIC). See Details.}

\item{pen.type}{a character string indicating the type of penalty used, either \code{"SCAD"} or \code{"LASSO"}.}

\item{tuneGrid}{a data.frame containing the values of the AIC and BIC criteria for each combination of the tuning parameters.}

\item{pen.weights}{a list with elements named \code{CURE} and \code{SURV}, containing the penalty weights. For more details, see \code{\link{penPHcure}}.}

\item{N}{the sample size (number of individuals).}

\item{K}{the number of unique failure times.}

\item{isTies}{logical value: \code{TRUE} in case of tied event times.}

\item{censoring}{the proportion of censored individuals.}

\item{which.X}{character string indicating the method used to transform the covariates included in the incidence (cure) component from time-varying to time-invariant. See \code{\link{penPHcure}} for more details.}

\item{survform}{a formula object with all variables involved in the latency (survival) component of the model.}

\item{cureform}{a formula object with all variables involved in the incidence (survival) component of the model.}

\item{call}{object of class \code{call}, with all the specified arguments.}
}
\description{
This class of objects is returned by the function \code{\link{penPHcure}} when is called with the argument \code{pen.type = "SCAD" | "LASSO"}. Objects of this class have methods for the functions \code{summary} and \code{predict}.
}
\details{
The lists \code{AIC} and \code{BIC} contain the results of the selected model based on the Akaike information criterion (AIC) and Bayesian Information Criterion (BIC), respectively. They are composed by the following elements:
\itemize{
  \item \code{crit}: value of the minimized AIC/BIC criterion.
  \item \code{b}: a numeric vector with the estimated regression coefficients in the cure (incidence) component.
  \item \code{beta}: a numeric vector with the true estimated coefficients in the survival (latency) component.
  \item \code{cumhaz}: a numeric vector with the estimated cumulative baseline hazard function at the unique event times (reported in the \code{"names"} attribute).
  \item \code{tune_params}: a list with elements named \code{CURE} and \code{SURV} containing the selected tuning parameters, which minimize the AIC/BIC criterion.
}
}
\seealso{
\code{\link{penPHcure}}
}
