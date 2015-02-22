kernelFixed <- function(data.x, 
                        data.y,
                        bandwidth,
                        kType, 
                        lType,
                        time,
                        distanceFunction, ...){

  if( any(bandwidth < 1.5e-8) ) stop("All bandwidths must be > 0")

  nCov <- ncol(data.x) - 2L
  patientIDs <- sort(unique(data.y[,1L]))
  nPatients <- length(patientIDs)

  xIs <- list()
  for( i in 1L:nrow(data.y) ) {
    xIs[[i]] <- list()
    xIs[[i]]$v <- which( data.x[,1L] == data.y[i,1L] )
    xIs[[i]]$n <- length(xIs[[i]]$v) 
  }

  yIs <- list()
  for( i in 1L:nPatients ) {
    yIs[[i]] <- list()
    yIs[[i]]$v <- which( data.y[,1L] == patientIDs[i] )
    yIs[[i]]$n <- length(yIs[[i]]$v)
  }

  lbd <- length(bandwidth)

  bHat <- matrix(data = 0.0, 
                 nrow = lbd, 
                 ncol = nCov,
                 dimnames = list(NULL,colnames(data.x)[3L:(nCov+2L)]))

  sdVec <- bHat

  results <- matrix(data = 0.0,
                    nrow = nCov,
                    ncol = 4L,
                    dimnames = list(paste("beta",0L:{nCov-1L},sep=""),
                                    c("estimate","stdErr","z-value","p-value")))

  guess <- NULL

  for( bd in 1L:lbd ) {

    if(distanceFunction != "distanceLV") cat("Bandwidth: ", bandwidth[bd], "\n")

    bHat[bd,] <- betaHat(data.y = data.y,
                         data.x = data.x,
                         bandwidth = bandwidth[bd],
                         kType = kType,
                         lType = lType,
                         tt = time,
                         xIs = xIs,
                         yIs = yIs,
                         nPatients = nPatients,
                         guess = guess,
                         distanceFunction = distanceFunction)

    results[,1L] <- bHat[bd,]

    guess <- bHat[bd,]

    sdVec[bd,] <- SD(data.y = data.y,
                     data.x = data.x,
                     bandwidth = bandwidth[bd],
                     kType = kType,
                     lType = lType,
                     bHat = bHat[bd,],
                     tt = time,
                     xIs = xIs,
                     yIs = yIs,
                     nPatients = nPatients,
                     distanceFunction = distanceFunction)

    results[,2L] <- sdVec[bd,]
    results[,3L] <- bHat[bd,]/sdVec[bd,]
    results[,4L] <- 2.0*pnorm(-abs(results[,3L]))

    print(results)
    cat("\n")
  }

  zv <- bHat/sdVec
  pv <- 2.0*pnorm(-abs(results[,3L]))

  return( list( "betaHat" = bHat,
                "stdErr"  = sdVec,
                "zValue"  = zv,
                "pValue"  = pv ) )

}
