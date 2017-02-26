aw <- function(s, snr) {
  
  s=c(s)
  x = noisew(1000,length(s)/1000,type = "gaussian",listen = FALSE,output = "array")
  
  np = sum(abs(x[1:length(x)])) ^ 2 / length(x)
  
  np = 10 * log10(np)
  
  p = sum(abs(s[1:length(s)])) ^ 2 / length(s)
  
  p = 10 * log10(p)
  
  n = snr * np
  
  np = p / n
  
  x=x*np
  
  r=s+t(x)
  
  return(t(r))
}