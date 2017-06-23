library("signal")

library("seewave")

library("Matrix")

library("pracma")

source("qpskmod.R")

source("qpskdemod.R")

source("mmse.R")

source("interpolate.R")

source("ber.R")

source("awgn.R")

N = 1024
pilotFrequency = 8
E = 0.1
Ncp = 64
L = 32
br=1e6;
f=br;
T=1/br;
t=c(T/3, 2*T/3, T);

D = as.integer(as.logical(as.integer(randn (1, N))))


D_Mod = mod(D)


D_Mod_serial = t(D_Mod)



PLoc = seq(1, N, pilotFrequency)


DLoc = setdiff(1:N, PLoc)


D_Mod_serial[PLoc] = E * D_Mod_serial[PLoc]

d_ifft = fft(D_Mod_serial, inverse = TRUE)




d_ifft_parallel = t(d_ifft)

CP_part = d_ifft_parallel[(nrow(d_ifft_parallel) - (Ncp - 1)):nrow(d_ifft_parallel),1:ncol(d_ifft_parallel)]


ofdm_cp = c(CP_part, d_ifft_parallel)



h = randn(L,1) + complex(imaginary = 1) * randn(L,1)

h = h / base::norm( h ,'2' )

H = fft(h)


d_channelled = complex(real = signal::filter(h, 1, Re(ofdm_cp)), imaginary = signal::filter(h, 1, Im(ofdm_cp)) )

channel_length = length(h)

H_power_dB = 10 * log10(abs(H %*% Conj(t(H))))

count = 0


snr_vector = seq(4, 44, 4)
r_NoCH=1

r_LS=1

r_MMSE=1


for (snr in snr_vector)  {
  SNR = snr
  
  
  count = count + 1
  
  
  print(c("step: ", count, " of: ", length(snr_vector)))
  
  
  ofdm_noisy_NoCH = c(aw(ofdm_cp, SNR, f))
  
  ofdm_noisy_with_chann = c(aw(d_channelled, SNR, f))
  
  
  ofdm_cp_removed_NoCH = ofdm_noisy_NoCH[(Ncp +1):(N+Ncp)]
  
  ofdm_cp_removed_with_chann = ofdm_noisy_with_chann[(Ncp +1):(N+Ncp)]
  
  d_parallel_fft_NoCH = fft(ofdm_cp_removed_NoCH, inverse = FALSE)
  
  d_parallel_fft_channel = fft(ofdm_cp_removed_with_chann, inverse = FALSE)
  
  TxP = D_Mod_serial[PLoc]

  RxP = d_parallel_fft_channel[PLoc]

  Hpilot_LS = RxP / TxP
  H_MMSE = MMSE(RxP, TxP, N, pilotFrequency, h, SNR)
  
  
  HData_LS = interpolate(t(Hpilot_LS), PLoc, N, "spline")
  
  
  
  HData_LS_parallel1 = HData_LS
  
  HData_MMSE_parallel1 = t(H_MMSE)
  
  
  
  
  d_received_NoCH = dmod(t(d_parallel_fft_NoCH), N)
  
  d_received_chann_LS = dmod(t(d_parallel_fft_channel) / HData_LS_parallel1, N)
  
  d_received_chann_MMSE = dmod(d_parallel_fft_channel / HData_MMSE_parallel1, N)
  
  D_no_pilots = D_Mod_serial[DLoc]
  
  Rec_d_NoCH = d_received_NoCH[DLoc]
  
  Rec_d_LS = d_received_chann_LS[DLoc]
  
  Rec_d_MMSE = d_received_chann_MMSE[DLoc]
  
  r_NoCH[count] = ber(D_no_pilots, Rec_d_NoCH)
  
  r_LS[count] = ber(D_no_pilots, Rec_d_LS)
  
  r_MMSE[count] = ber(D_no_pilots, Rec_d_MMSE)
  
}

semilogy(snr_vector, r_NoCH, type ='l', col="red")

par(new=TRUE)

semilogy(snr_vector, r_LS, type ='l', col="green")

par(new=TRUE)

semilogy(snr_vector, r_MMSE, type ='l', col="black")

H_power_esti_dB_LS = 10 * log10(abs(HData_LS_parallel1 %*% Conj(t(HData_LS_parallel1))))

H_power_esti_dB_MMSE = 10 * log10(abs(HData_MMSE_parallel1 %*% Conj(t(HData_MMSE_parallel1))))

plot(H_power_dB[seq(1, length(H_power_dB), 8)], type ='l', col="green")

par(new=TRUE)

plot(H_power_esti_dB_LS[seq(1, ncol(H_power_esti_dB_LS), 8)], type ='l', col="black")

par(new=TRUE)

plot(H_power_esti_dB_MMSE[1, seq(1, ncol(H_power_esti_dB_MMSE), 8)], type ='l', col="red")
###############THE END##################Thankyou For your patiance################