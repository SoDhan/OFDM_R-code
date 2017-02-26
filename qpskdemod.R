dmod <- function(sigrx)  {
  
   rxd=c(1)
   for (i in 1:length(sigrx) / 2) {
    zin = (cos(2 * pi * ((2*i)/(length(sigrx))))) * sigrx[2*(i+1)]
    
    zinint = (trapz((i/(length(sigrx))), zin)) * 2000
    
    if (zinint > 0)  {
      rxin = 1
      
    }
    else  {
      rxin = 0
      
    }
    
    zqd1 = (sin(2 * pi * ((2*i)/(length(sigrx))))) * sigrx[(2*(i+1)-1)]
    
    zqdint = (trapz((i/(length(sigrx))), zqd1)) * 2000
    
    if (zqdint > 0) {
      rxqd = 1
      
    }
    else  {
      rxqd = 0
      
    }
    rxd = c(rxd, rxin, rxqd)
    
  }
  return(rxd)
}