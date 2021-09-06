Adlas <- function(A,b,lambda,options=list())
{
   # Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes

   # This file is part of SLOPE Toolbox version 1.0.
   #
   #    The SLOPE Toolbox is free software: you can redistribute it
   #    and/or  modify it under the terms of the GNU General Public License
   #    as published by the Free Software Foundation, either version 3 of
   #    the License, or (at your option) any later version.
   #
   #    The SLOPE Toolbox is distributed in the hope that it will
   #    be useful, but WITHOUT ANY WARRANTY; without even the implied
   #    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   #    See the GNU General Public License for more details.
   #
   #    You should have received a copy of the GNU General Public License
   #    along with the SLOPE Toolbox. If not, see
   #    <http://www.gnu.org/licenses/>.

   # -------------------------------------------------------------
   # Start times
   # -------------------------------------------------------------
   t0 <- proc.time()[3]
   
   # -------------------------------------------------------------
   # Define function for retrieving option fields with defaults
   # -------------------------------------------------------------
   getDefaultField <- function(options,name,default)
   {  if (!is.null(options[[name]]))
      {  return(options[[name]])  }
      else
      {  return(default)  }
   }
   
   # -------------------------------------------------------------
   # Parse parameters
   # -------------------------------------------------------------
   iterations <- getDefaultField(options,"iterations", 10000)
   verbosity  <- getDefaultField(options,"verbosity" , 1)
   optimIter  <- getDefaultField(options,"optimIter" , 1)
   gradIter   <- getDefaultField(options,"gradIter"  , 20)
   tolInfeas  <- getDefaultField(options,"tolInfeas" , 1e-6)
   tolRelGap  <- getDefaultField(options,"tolRelGap" , 1e-6)
   xInit      <- getDefaultField(options,"xInit"     , vector())
   
   # -------------------------------------------------------------
   # Ensure that lambda is non-increasing
   # -------------------------------------------------------------
   n <- length(lambda)
   matrix(lambda,c(1,n))
   if ((n > 1) && any(lambda[2:n] > lambda[1:n-1]))
   {  stop("Lambda must be non-increasing.")
   }
   
   # -------------------------------------------------------------
   # Initialize
   # -------------------------------------------------------------
   # Get problem dimension
   n <- ncol(A)
   
   # Get initial lower bound on the Lipschitz constant
   rndKind  <- RNGkind()
   rndState <- tryCatch({.Random.seed},error=function(m) {} )
   set.seed(0,kind="Mersenne-Twister",normal.kind="Inversion")
   x <- matrix(rnorm(n),c(n,1))
   x <- x / sqrt(sum(x^2))
   x <- t(A) %*% (A %*% x)
   L <- sqrt(sum(x^2))
   if (!is.null(rndState))
      set.seed(rndState[-1],kind=rndKind[1],normal.kind=rndKind[2])
   
   # Set constants
   STATUS_RUNNING    <- 0
   STATUS_OPTIMAL    <- 1
   STATUS_ITERATIONS <- 2
   STATUS_MSG        <- c('Optimal','Iteration limit reached')
   
   # Initialize parameters and iterates
   if (length(xInit) == 0) xInit <- matrix(0,n,1)
   t       <- 1
   eta     <- 2
   lambda  <- matrix(lambda,nrow=length(lambda))
   b       <- matrix(b,nrow=length(b))
   x       <- xInit
   y       <- x
   Ax      <- A %*% x
   fPrev   <- Inf
   iter    <- 0
   status  <- STATUS_RUNNING
   Aprods  <- 2
   ATprods <- 1
   
   # Deal with Lasso case
   modeLasso <- (length(lambda) == 1)
   if (modeLasso)
      proxFunction <- function(x,lambda) { return(sign(x) * pmax(abs(x) - lambda,0)) }
   else
      proxFunction <- function(v1,v2) { return(proxSortedL1(v1,v2)) }
   
   if (verbosity > 0)
   {  printf <- function(...) invisible(cat(sprintf(...)))
      printf('%5s  %9s   %9s  %9s  %9s\n','Iter','||r||_2','Gap','Infeas.','Rel. gap')
   }
   
   # -------------------------------------------------------------
   # Main loop
   # -------------------------------------------------------------
   while (TRUE)
   {    
      # Compute the gradient at f(y)
      if ((iter %% gradIter) == 0) # Includes first iterations
      {  r <- (A %*% y) - b
         g <- t(A) %*% r
         f <- as.double(crossprod(r)) / 2
      }
      else
      {  r <- (Ax + ((tPrev - 1) / t) * (Ax - AxPrev)) - b
         g <- t(A) %*% r
         f <- as.double(crossprod(r)) / 2
      }
      
      # Increment iteration count
      iter <- iter + 1
      
      # Check optimality conditions
      if ((iter %% optimIter) == 0)
      {  # Compute 'dual', check infeasibility and gap
         if (modeLasso)
         { infeas    <- max(norm(g,'I')-lambda,0)
           objPrimal <- f + lambda * norm(y,'1')
           objDual   <- -f - as.double(crossprod(r,b))
         }
         else
         {  gs     <- sort(abs(g), decreasing=TRUE)
            ys     <- sort(abs(y), decreasing=TRUE)
            infeas <- max(max(cumsum(gs-lambda)),0)
            
            # Compute primal and dual objective
            objPrimal <-  f + as.double(crossprod(lambda,ys))
            objDual   <- -f - as.double(crossprod(r,b))
         }  
         
         # Format string
         if (verbosity > 0)
            str <- sprintf('   %9.2e  %9.2e  %9.2e',objPrimal - objDual, infeas/lambda[[1]], abs(objPrimal - objDual) / max(1,objPrimal))
         
         # Check primal-dual gap
         if ((abs(objPrimal - objDual)/max(1,objPrimal) < tolRelGap) &&
            (infeas < tolInfeas * lambda[[1]]))
            status <- STATUS_OPTIMAL
      }
      else
      {   str <- ''
      }
      
      if (verbosity > 0)
      {  if ((verbosity == 2) ||
             ((verbosity == 1) && ((iter %% optimIter) == 0)))
         {  printf('%5d  %9.2e%s\n', iter,f,str)
         }
      }
      
      # Stopping criteria
      if ((status == 0) && (iter >= iterations))
         status <- STATUS_ITERATIONS

      if (status != 0)
      {  if (verbosity > 0)
            printf('Exiting with status %d -- %s\n', status, STATUS_MSG[[status]])
         break
      }
      
      # Keep copies of previous values
      AxPrev <- Ax
      xPrev  <- x
      fPrev  <- f
      tPrev  <- t
      
      # Lipschitz search
      while (TRUE)
      {  # Compute prox mapping
         x <- proxFunction(y - (1/L)*g, lambda/L)
         d <- x - y
      
         Ax <- A %*% x
         r  <- Ax-b
         f  <- as.double(crossprod(r))/2
         q  <- fPrev + as.double(crossprod(d,g)) + (L/2)*as.double(crossprod(d))
                                 
         Aprods <- Aprods + 1
                                 
         if (q >= f*(1-1e-12))
            break
         else
            L <- L * eta
      } # Lipschitz search
      
      # Update
      t <- (1 + sqrt(1 + 4*t^2)) / 2
      y <- x + ((tPrev - 1) / t) * (x - xPrev)
   } # While (TRUE)
   
   
   # Information structure
   info <- c()
   info$runtime   <- proc.time()[3] - t0
   info$Aprods    <- Aprods + ceiling(iter / gradIter)
   info$ATprods   <- ATprods + iter
   info$objPrimal <- objPrimal
   info$objDual   <- objDual
   info$infeas    <- infeas
   info$status    <- status
   
   info$x         <- y
   info$L         <- L
   
   return(info)
} # Function Adlas
   
   
   
