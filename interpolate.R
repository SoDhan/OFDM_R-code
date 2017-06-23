interpolate <- function(in1, locpil, n, method) {
  if (locpil[1] > 1) {
    slope = (in1[2] - in1[1]) / (locpil[2] - locpil[1])
    
    in1 = c(in1[1] - slope * (locpil[1] - 1), in1)
    
    locpil = c(1, locpil)
    
  }
  if (locpil[length(locpil)] < n)  {
    slope = (in1[length(in1)] - in1[length(in1) - 1]) / (locpil[length(locpil)] - locpil[length(locpil)] - 1)
    
    in1 = c(in1, in1[length(in1)] + slope * (n - locpil[length(locpil)]))
    
    locpil = c(locpil, n)
    
  }
  l = 1:n
  if (method == "linear")  {
    k = interp1(locpil ,Real(in1) ,l ,method = "linear")
    
    return(k)
    
  }
  else  {
    
    k = interp1(locpil,Real(in1) ,l ,method = "spline")
    
    return(k)
    
  }
}