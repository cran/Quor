% File Quor/man/conf.statement.Rd
% Part of the Quor package for R (http://www.R-project.org)
% Copyright (C) 2014, Adriano Polpo, Carlos A. de B. Pereira, Cassio P. de Campos.
% Distributed under FDL 1.3 or later
% http://www.gnu.org/licenses/fdl-1.3.html

\name{conf.statement}

\alias{conf.statement}
\alias{conf.statement.pooled}
\alias{print.conf.statement}

\title{Confidence Statement for Ordered Population Quantiles}
\description{
  Evaluates the Confidence Statement for Ordered Population Quantiles.
}
\usage{
 conf.statement(data,quantiles=NULL,ordering=NULL,verbose=TRUE,logscale=FALSE)
 conf.statement.pooled(data,quantiles=NULL,ordering=NULL,
                       verbose=TRUE,logscale=FALSE)
}
\arguments{
\item{data}{is a list with all groups to be test. Each element can be
            a vector with the elements, or a matrix, in which case each
            row will be considered as a different covariate to be tested.
            All the elements in the list must have the exact same dimension.}
\item{quantiles}{a vector of elements in [0,1] with length equal to the
                 number of groups. It tells us which quantile will be used
                 for each group. If null, then medians are compared by default.}
\item{ordering}{a matrix containing one permutation of 1:n per row, where n
                is the number of groups. It tells us which orderings are to be
                tested. If null, then it is assumed all orderings to be tested.}
\item{verbose}{if TRUE, display results on screen, FALSE otherwise.}
\item{logscale}{if TRUE, return confidence values in log scale.}
}
\details{
  The confidence statement is evaluated for the population quantiles consedering:
  The population median of the group 1 is lower than 
  the population median of the group 2, and this is lower than
  the population median of the group 3, and so...
  That is, M1 < M2 < ... < Mk, where
  Mj is the population median of the group j.

  The group j is the j-th vector in the data (a list object).

  If the conf.statement.pooled version is used, then each groups is checked
  against all the others, which are pooled together, so the computation always
  involves only two groups at a time (the target group and the pooled one).
}
\value{
  A list (with class "conf.statement") with components:
  \item{call}{command evaluated.}
  \item{confidence}{the (log) confidence statement level as a matrix permutations by variables.}
  \item{quantiles}{the quantiles of interested for each of the k groups.}
  \item{order}{the permutations of groups that were used.}
  \item{total.groups}{the total number of groups.}
  \item{total.covariates}{the total number of variables on which the confidence was computed.}
  \item{run.time}{Total time spent to run it.}
}
\examples{
  set.seed(42)
  # Statement: Q1 < Q2 < Q3, where Qj is the population
  #   median of the group j.
  data <- NULL
  data$x1 <- rnorm(10,0,1)
  data$x2 <- rnorm(10,1,1)
  data$x3 <- rnorm(10,2,1)
  conf.statement(data,ordering=as.integer(c(1,2,3)),verbose=FALSE)

  # Statement: Q2 < Q3 < Q1, where Qj is the population
  #   3rd quartile of the group j.
  conf.statement(data,ordering=as.integer(c(2,3,1)),
                 quantiles=c(0.75,0.75,0.75))

  # Evaluates statements Q1 < Q2 and Q2 < Q1 for the
  # data set gleason7, where Qj is the median of
  # population j
  data(gleason7)
  d <- list(x1 = gleason7[1:5,1], x2 = gleason7[,2])
  conf.statement(d,verbose=FALSE)
}
