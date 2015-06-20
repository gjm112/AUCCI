#' @title Qstat
#' 
#' 
#' @description Computes Q1 and Q2 from somewhere.  
#' 
#' @details A
#' 
#' @param data 
#' @param n.x
#' @param n.y
#' @param disease
#' @param marker
#' 
#' @return A data.frame containing Q1 and Q2  
#' 
#' @export Qstat


Qstat <- function(data, n.x, n.y, disease="disease", marker="marker") {
  x = "$"(data, marker)["$"(data, disease)==0]
  y = "$"(data, marker)["$"(data, disease)==1]
  Q2 <- Q1 <- 0
  for (i in 1:n.x) {
    Q1 = Q1 + ( sum(y > x[i]) + sum(y==x[i])/2 )^2
  }
  for (j in 1:n.y) {
    Q2 = Q2 + ( sum(y[j] > x) + sum(y[j]==x)/2 )^2
  }
  Q1 = Q1 /(n.x * n.y^2); Q2 = Q2 /(n.x^2 * n.y)
  return(data.frame(Q1 = Q1, Q2 = Q2))
}