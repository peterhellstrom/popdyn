#' @export
MatProp <- function(A){ # Check properties of matrix A.

if (Square(A)){
   s <- sqrt(floor(length(A)))
   {if (s > 1){
        Ap <- 1*(A>0)
        for (i in 1:s){Ap[i,i] <- 1}
        App <- Ap
        if (s > 2){for (i in 1:(s-2)){App <- 1*((App %*% Ap)>0)}}
        reducible <- (sum(App) < s*s)

        {if (!reducible){
           Ap <- 1*(A>0)
           App <- Ap
           for (i in 1:((s-1)^2)){App <- 1*((App %*% Ap)>0)}
           primitive <- !(sum(App) < s*s)
           }
           else{primitive <- FALSE}
        }
      }
      else{primitive <- (sum(A) > 0)
           reducible <- !primitive}
   }

   EV <- eigen(t(A),symmetric=FALSE)
   EVval <- EV$values
   EVvec <- EV$vectors
   j <- 1
   if (!(primitive)){
     repeat{if ((Im(EVval[j]) > pi/s) | (Re(EVval[j]) < 0)){j <- j+1}
              else{break}}
     }
   v <- abs(Re(EVvec[,j]))
   norm <- sum(v)
   v <- v/norm                         # right eigenvector (reproductive value)

   EV <- eigen(A,symmetric=FALSE)
   EVval <- EV$values
   EVvec <- EV$vectors
   j <- 1
   if (!(primitive)){
     repeat{if ((Im(EVval[j]) > pi/s) | (Re(EVval[j]) < 0)){j <- j+1}
              else{break}}
     }

   lambda <- abs(EVval[j])             # dominant eigenvalue
   w <- abs(Re(EVvec[,j]))
   norm <- sum(w)
   w <- w/norm                         # left eigenvector (stable structure)

   d <- 1
   repeat{
     if (d == s) {break}
       else{if (abs(1-abs(EVval[d+1])/lambda) < 1e-15){d <- d+1}
              else{break}}
     }                                  # period of matrix

   rho <- Inf                           # damping ratio
   if (d < s){rho <- lambda/abs(EVval[d+1])}

   S <- (v %*% t(w))/sum(w * v)         # sensitivity matrix
   E <- S*A/lambda                      # elasticity matrix

   MatProp <- list(EVval = EVval,       # eigenvalues
                   EVvec = EVvec,       # eigenvectors
                   S = S,               # sensitivity matrix (lambda1 on A)
                   E = E,               # elasticity matrix
                   reducible = reducible,
                   primitive = primitive,
                   d = d,               # period
                   w = w,               # stable structure
                   v = v,               # reproductive value
                   rho = rho,           # damping ratio
                   lambda1 = lambda     # dominant eigenvalue
                                      )
   }
MatProp
}
