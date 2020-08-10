\name{MAPCO}
\alias{MAPCO}
\alias{MAPCO_MLE}
\alias{data}

\title{A mixed model for Mendelian randomization in the presence of correlated pleiotropic effects}
\usage{MAPCO(by,bx,se_y,se_x,int_beta=0,int_sigma_alpha=NA,b0=0,tol=1e-8,n_iter=3000)}
\description{The package provides a Mendelian randomization approach for the inference of causal relationship between traits based on summary-level GWAS data, which takes account of the horizontal pleiotropy and correlation between genetic effects on traits. It requires neither the InSIDE assumption nor the conditions on the proportion of valid instrumental variables.
}

\arguments{
  \item{by}{A vector including the estimated genentic effects on the outcome.}
  \item{bx}{A vector including the estimated genentic effects on the exposure.}
  \item{se_y}{A vector including the standard errors of the estimated genentic effects on the outcome.}
  \item{se_x}{A vector including the standard errors of the estimated genentic effects on the exposure.}
  \item{int_beta}{The initial value of the causal effect for the EM algorithm. The default value is zero}
  \item{int_sigma_alpha}{The initial value of the standard deviation of the pleiotropy for the EM algorithm. If it is missing, then the
        standard deviation of bx would be used.}
  \item{b0}{The causal effect under the null hypothesis. The default value is zero.}
  \item{tol}{The tolenrence used in the EM algorithm. The default value is 1e-8.}
  \item{n_iter}{The maximum number of iterations in the EM algorithm. The default value is 3000.}
}

\value{
\item{beta}{The estiamted causal effect of the exposure on the outcome.}
\item{sigma_alpha}{The estiamted standard deviation of the pleiotropy.}
\item{rho}{The estiamted correlation coefficient between genetic effects on traits.}
\item{iter}{The number of iterations in the EM algorithm.}
\item{pval}{The p value for testing the null hypothesis for causal effect.}
\item{CI}{The 95 percent confidence interval for causal effect.}
}


\examples{
  by <- data[,1]
  bx <- data[,2]
  se_y <- data[,3]
  se_x <- data[,4]
  result <- MAPCO(by,bx,se_y,se_x)
  result
}



