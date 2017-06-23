dmod <- function(sigtx, N)  {
  br = 1e6
  f = br
  T = 1 / br
  t = c(T/3, 2*T/3, T)
  Rx_data = c()
  Rx_sig = sigtx
  for (i in 1:(N/2)) {
    la = (i - 1) * length(t) + 1
    ra = i * length(t)
    Z_in = Rx_sig[la:ra] * cos(2 * pi * f * t)
    Z_in_intg = (trapz(t, Z_in)) * (2 / T)
    Rx_in_data = ifelse(Re(Z_in_intg) > 0, 1, 0)
    Z_qd = Rx_sig[la:ra] * sin(2 * pi * f * t)
    Z_qd_intg = (trapz(t, Z_qd)) * (2 / T)
    Rx_qd_data = ifelse(Re(Z_qd_intg) > 0, 1, 0)
    Rx_data = c(Rx_data, Rx_in_data, Rx_qd_data)
  }
  return(Rx_data)
}