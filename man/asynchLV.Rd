\name{asynchLV}
\alias{asynchLV}
\title{
       Regression Analysis Using Last Value Carried Forward
}

\description{
             Estimation of regression models for sparse asynchronous 
             longitudinal observations using the last value carried forward 
             approach.
}

\usage{
       asynchLV(data.x, data.y, lType = "identity", ...)
}

\arguments{
  \item{data.x}{
                A data.frame of covariates. The structure of the data.frame must 
                be \{patient ID, time of measurement, measurement(s)\}. 
                Patient IDs must be of class integer or be able to be coerced 
                to class integer without loss of information.
                Missing values must be indicated as NA.
                All times will automatically be rescaled to [0,1].
  }

  \item{data.y}{
                A data.frame of response measurements. The structure of the 
                data.frame must be 
                \{patient ID, time of measurement, measurement\}. 
                Patient IDs must be of class integer or be able to be coerced 
                to integer without loss of information.
                Missing values must be indicated as NA.
                All times will automatically be rescaled to [0,1].
  }

  \item{lType }{
                An object of class character indicating the type of link 
                function to use for the regression model. 
                Must be one of \{"identity","log","logistic"\}.
  }

  \item{... }{Ignored. }

}

\details{
         For lType = "log" and lType = "logistic", parameter estimates are
         obtained by minimizing the estimating equation using R's optim()
         with method="Nelder-Mead"; 
         all other settings take their default values.

         For lType = "identity", parameter estimates are obtained use solve().
}
\value{
  A list is returned, the elements of which are named vectors:

  \item{betaHat }{The estimated model coefficients.}
  \item{stdErr  }{The standard error for each coefficient. }
  \item{zValue  }{The estimated z-value for each coefficient.}
  \item{pValue  }{The p-value for each coefficient.}

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

  res <- asynchLV(data.x = TI.x, 
                  data.y = TI.y,
                  lType = "identity")
  

}

