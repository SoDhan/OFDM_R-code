cplxfilt <- function (a,b,c) {
  library("signal")
  d <- complex(real = signal::filter(a, b, re(c)), imaginary = signal::filter(a, b, im(c)))
  return (d)
}