mod <- function(d)  {
  data_NZR = 2 * d - 1
  s_p_data = matrix(data_NZR, 2, N/2)
  br = 1e6
  f = br
  T = 1 / br
  t = c(T/3, 2*T/3, T)
  y = c()
  y_in = c()
  y_qd = c()
  for (i in 1:(N/2)) {
    y1 = s_p_data[1, i] * cos(2 * pi * f * t)
    y2 = s_p_data[2, i] * sin(2 * pi * f * t)
    y_in = c(y_in, y1)
    y_qd = c(y_qd, y2)
    y = c(y, y1 + y2)
  }
  return(y)
}