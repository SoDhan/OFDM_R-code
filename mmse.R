MMSE <- function(prx, ptx, n, pilfreq, hcir, snr)  {
  
  noiseVar = 10 ^ (snr * 0.1)
  
  np = n / pilfreq
  
  hls = ptx/prx
  
  k = 1:length(hcir) - 1
  
  hh = t(hcir) %*% hcir
  
  tmp = (Conj(hcir) * hcir) * k
  
  r = sum(tmp) / hh
  
  tmp=c(tmp)
  
  r2 = (t(tmp) %*% k) / hh
  
  t_rms = sqrt(r2 - r ^ 2)
  
  D = (complex(imaginary = 1) * 2 * pi * t_rms/ n)
  
  # Denomerator of Eq. (6.16) page 192
  k1 = matrix(rep((1:n - 1), np), ncol = np)
  
  k2 = matrix(rep((1:np - 1), n), nrow = n, byrow = TRUE)
  
  rf = 1 / (1 + (k1 - k2 * pilfreq) * D[1] )
  
  k3 = matrix(rep((1:np - 1), np), ncol = np)
  
  k4 = matrix(rep(1:np - 1, np), ncol = np, byrow = TRUE)
  
  rf2 = 1 / (1 +(k3 - k4) * pilfreq * D[1])
  
  Rhp = rf
  
  Rpp = rf2 + diag(1, length(hls)) / noiseVar
  
  k = t(Rhp %*% inv(Rpp) %*% hls)
  
  return(k)
  
}
