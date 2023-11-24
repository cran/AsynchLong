\name{asynchTI}
\alias{asynchTI}
\title{
Regression Analysis for Time-Invariant Coefficients
}
\description{
Estimation of regression models for sparse asynchronous longitudinal
observations with time-invariant coefficients.}

\usage{
asynchTI(data.x, data.y, kType = "epan", lType = "identity", bw = NULL, 
         nCores = 1, verbose = TRUE, ...)
}
\arguments{
  \item{data.x}{
                A data.frame of covariates. 
                The structure of the data.frame must be 
                \{patient ID, time of measurement, measurement(s)\}. 
                Patient IDs must be of class integer or be able to be coerced to 
                class integer without loss of information.
                Missing values must be indicated as NA.
                All times will automatically be rescaled to [0,1].
  }

  \item{data.y}{
                A data.frame of response measurements. 
                The structure of the data.frame must be 
                \{patient ID, time of measurement, measurement\}. 
                Patient IDs must be of class integer or be able to be coerced to 
                class integer without loss of information.
                Missing values must be indicated as NA.
                All times will automatically be rescaled to [0,1].
  }

  \item{kType }{
                An object of class character indicating the type of 
                smoothing kernel to use in the estimating equation. 
                Must be one of \{"epan", "uniform", "gauss"\}, where
                "epan" is the Epanechnikov kernel and "gauss" is the
                Gaussian kernel.
  }

  \item{lType }{
                An object of class character indicating the type of link 
                function to use for the regression model. 
                Must be one of \{"identity","log","logistic"\}.
  }

  \item{bw    }{
                If provided, bw is an object of class numeric or a 
                numeric vector containing the bandwidths for which parameter
                estimates are to be obtained. 
                If NULL, an optimal bandwidth will be determined
                using an adaptive selection procedure. 
                The range of the bandwidth search space is taken
                to be \eqn{2*(Q3 - Q1)*n^{-0.7} to 2*(Q3 - Q1)*n^{-0.3}}{2*(Q3 - Q1)*n^(-0.7) to 2*(Q3 - Q1)*n^(-0.3)},
                where Q3 is the 0.75 quantile and Q1 is the 0.25 quantile
                of the pooled sample of measurement times for the covariate
                and response, and n is the number of patients.
                See original reference for details of the 
                selection procedure.
  }
  \item{nCores}{
                A numeric object. For auto-tune method, the number of cores
                to employ for calculation. If nCores > 1, the bandwidth
                search space will be distributed across the cores
                using parallel's parLapply.
  }
  \item{verbose }{ An object of class logical. TRUE results in screen prints.}
  \item{... }{Ignored. }
}
\details{
         For lType = "log" and lType = "logistic", parameter estimates are
         obtained by minimizing the estimating equation 
         using optim() with method="Nelder-Mead"; 
         all other arguments take their default values.

         For lType = "identity", parameter estimates are obtained using solve().
}
\value{
A list is returned. If bandwidths are provided, each element 
of the list is a matrix,
where the ith row corresponds to the ith bandwidth of argument ``bw" and the
columns correspond to the model parameters. If the bandwidth is
determined automatically, each element is a named vector calculated
at the optimal bandwidth.

  \item{betaHat }{The estimated model coefficients.}
  \item{stdErr  }{The standard error for each coefficient.}
  \item{zValue  }{The estimated z-value for each coefficient.}
  \item{pValue  }{The p-value for each coefficient.}

If the bandwidth is determined automatically, two additional list
elements are returned:

  \item{optBW   }{The estimated optimal bandwidth for each coefficient.}
  \item{minMSE  }{The mean squared error at the optimal bandwidth for
                  each coefficient.}
}

\references{
Cao, H., Zeng, D., and Fine, J. P. (2015)
Regression Analysis of sparse asynchronous longitudinal data.
Journal of the Royal Statistical Society: Series B, 77, 755-776.
}
\author{
Hongyuan Cao, Donglin Zeng, Jason P. Fine, and Shannon T. Holloway
}

\examples{

  data(asynchDataTI)

  res <- asynchTI(data.x = TI.x, 
                  data.y = TI.y,
                  bw = c(0.05, 0.03),
                  kType = "epan", 
                  lType = "identity")

}

