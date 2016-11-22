#' lmp 
#'
#' Extracts the p-value from lm(). Source: Stephen Turner (http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html)
#' In short, lmp extracts the F-statistic from the lm output and calculates the p-value based on the F-distribution. 
#' @param modelobject lm output
#' 
#' @author Stephen Turner
#' 
#' @return A p-value for a linear regression as obtained through lm()
#' @keywords lm, p-value, fstatistic
#' @export
#' @examples 
#' set.seed(42)
#' n=20
#' d=data.frame(x1=rbinom(n,2,.5), x2=rbinom(n,2,.5))
#' d=transform(d, y=x1+x2+rnorm(n))
#' 
#' #fit the linear model
#' fit=lm(y ~ x1 + x2, data=d)
#' 
#' summary(fit) #shows that the F-test is 0.006641
#' names(summary(fit)) #can't access that p-value using this!
#' names(fit) # this doesn't work either
#' lmp(fit) # uses the above function to capture the F-test p-value.
#' #Source: http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}