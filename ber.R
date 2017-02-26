ber <- function(s1, s2)  {
  
  s1 = as.logical(as.integer(s1))
  
  s2 = as.logical(as.integer(s2))
  
  s3 = as.integer(xor(s1, s2))
  
  e = nnzero(s3)/length(s1)
  
  return(e)
  
}