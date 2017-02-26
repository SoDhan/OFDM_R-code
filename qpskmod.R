mod <- function(d)  {
  
  dnzr = 2 * d - 1
  
  spd = matrix(c(dnzr), nrow = 2, ncol = length(d) / 2)
  
  y = c(1)
  
  for (i in 1:length(d) / 2) {
    y1 = spd[1, i] * cos(2 * pi * (i/length(d)))
    
    y2 = spd[2, i] * sin(2 * pi * (i/length(d)))
    
    y = c(y, y1+y2)
    
  }
  y = matrix(y, nrow = 1)
  
  return(y)
  
}

