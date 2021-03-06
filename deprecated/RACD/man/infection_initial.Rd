% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ode_infection.R
\name{infection_initial}
\alias{infection_initial}
\title{Within-host Infection State ODE Model}
\usage{
infection_initial(a, zeta, psi, theta, ...)
}
\arguments{
\item{a}{age of individual}

\item{zeta}{individual biting heterogeneity}

\item{psi}{relative biting rate (used to calculate psi(a))}

\item{theta}{vector of model parameters (see \code{\link{RACD_Parameters}})}

\item{...}{additional parameters passed to numerical integrator \code{\link[deSolve]{ode}}}
}
\description{
ODE model for calculating the probability that an individual is in each state given their age and EIR heterogeneity attributes.
This function calls \code{\link{infection_ode}}
}
