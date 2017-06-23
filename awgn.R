aw <- function(s, snr, f) {
  
  x = noisew(f,length(s)/f,type = "gaussian", output = "array")
  
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