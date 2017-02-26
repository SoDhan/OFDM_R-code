library("seewave")

library("Matrix")

library("pracma")

source("qpskmod.R")


source("qpskdemod.R")


source("mmse.R")


source("interpolate.R")


source("ber.R")


source("awgn.R")


#give input

N = 128  #readline(prompt = "symbol length")


pilotFrequency = 16  #readline(prompt = "pilot frequency = ")


E = 0.1  #readline(prompt = "pilot energy = ")


Ncp = 16 #readline(prompt = " length of cp = ")


L = 128 #readline(prompt = "length of the channel (number of taps) = ")

## ----------------------------------------------------------------------
## data generation
D = as.integer(as.logical(as.integer(randn (1, N))))

## mapping (baseband modulation )
D_Mod_serial = mod(D)


## serial to parallel
#D_Mod_serial = t(D_Mod)


## specify Pilot & Date Locations
PLoc = seq(1, N, pilotFrequency)

# location of pilots
DLoc = setdiff(1:N, PLoc)

# location of data

## Pilot Insertion
D_Mod_serial[PLoc] = E * D_Mod_serial[PLoc]



###figure;
###imagesc(abs(D_Mod_serial ))
## inverse discret Fourier transform (IFFT)
#  Amplitude Modulation
d_ifft = fft(D_Mod_serial, inverse = TRUE)



## parallel to serail
d_ifft_parallel = t(d_ifft)



## Adding Cyclic Prefix
CP_part = d_ifft_parallel[(nrow(d_ifft_parallel) - (Ncp - 1)):nrow(d_ifft_parallel),1:ncol(d_ifft_parallel)]

# this is the Cyclic Prefix part to be appended.
ofdm_cp = c(CP_part, d_ifft_parallel)



## generating random channel
h = randn(L,1) + complex(imaginary = 1) * randn(L,1)

h = h / base::norm( h ,'2' )

# normalization
H = fft(h)

# Frequency-Domain Channel

d_channelled = t(filter(ofdm_cp, h, method = "convolution", sides = 1,circular = TRUE))

# channel effect
channel_length = length(h)

# True channel and its time-domain length
H_power_dB = 10 * log10(abs(H %*% Conj(t(H))))

# True channel power in dB

## add noise
count = 0


snr_vector = seq(1, 40, 4)
r_NoCH=1

r_LS=1

r_MMSE=1


for (snr in snr_vector)  {
  SNR = snr
  
  
  count = count + 1
  
  
  print(c("step: ", count, " of: ", length(snr_vector)))
  
  
  ofdm_noisy_NoCH = c(aw(ofdm_cp, SNR))
  
  ofdm_noisy_with_chann = c(aw(d_channelled, SNR))
  
  
  ## receiver
  #Remove Cyclic Prefix
  ofdm_cp_removed_NoCH = ofdm_noisy_NoCH[(Ncp +1):(N+Ncp)]
  
  ofdm_cp_removed_with_chann = ofdm_noisy_with_chann[(Ncp +1):(N+Ncp)]
  
  ## Discret Fourier transform (FFT)
  #  Amplitude Demodulation
  d_parallel_fft_NoCH = fft(ofdm_cp_removed_NoCH, inverse = FALSE)
  
  d_parallel_fft_channel = fft(ofdm_cp_removed_with_chann, inverse = FALSE)
  
  
  ## channel estimation
  # Extracting received pilots
  TxP = D_Mod_serial[PLoc]
  # trnasmitted pilots
  RxP = d_parallel_fft_channel[PLoc]
  # received pilots
  # Least-Square Estimation
  Hpilot_LS = RxP / TxP
  # LS channel estimation
  # MMSE Estimation:-
  H_MMSE = MMSE(RxP, TxP, N, pilotFrequency, h, SNR)
  
  # Interpolation p--->N
  HData_LS = interpolate(t(Hpilot_LS), PLoc, N, "spline")
  # Linear/Spline interpolation
  
  ## parallel to serial
  HData_LS_parallel1 = t(HData_LS)
  
  HData_MMSE_parallel1 = t(H_MMSE)
  
  
  
  ## demapping
  d_received_NoCH = dmod(t(d_parallel_fft_NoCH))
  # No Channel
  d_received_chann_LS = dmod(t(d_parallel_fft_channel) / HData_LS_parallel1)
  # LS channel estimation
  d_received_chann_MMSE = dmod(d_parallel_fft_channel / HData_MMSE_parallel1)
  # MMSE channel estimation
  ## Removing Pilots from received data and original data
  D_no_pilots = D[DLoc]
  # removing pilots from D
  Rec_d_NoCH = d_received_NoCH[DLoc]
  # removing pilots from d_received_NoCH
  Rec_d_LS = d_received_chann_LS[DLoc]
  # removing pilots from d_received_chann_LS
  Rec_d_MMSE = d_received_chann_MMSE[DLoc]
  # removing pilots from d_received_chann_MMSE
  ## Calculating BER
  r_NoCH[count] = ber(D_no_pilots, Rec_d_NoCH)
  
  r_LS[count] = ber(D_no_pilots, Rec_d_LS)
  
  r_MMSE[count] = ber(D_no_pilots, Rec_d_MMSE)
  
}

semilogy(snr_vector, r_NoCH,xlab="LS CE",ylab="MMSE CE",main ="orig. No Channel")

semilogy(snr_vector, r_LS)

semilogy(snr_vector, r_MMSE)

H_power_esti_dB_LS = 10 * log10(abs(HData_LS_parallel1 %*% Conj(t(HData_LS_parallel1))))
# Estimated channel power in dB
H_power_esti_dB_MMSE = 10 * log10(abs(HData_MMSE_parallel1 %*% Conj(t(HData_MMSE_parallel1))))
# Estimated channel power in dB
plot(H_power_dB[seq(1, length(H_power_dB), 8)],
     "b",
     xlab = "Time in samples",
     ylab = "Magnitude of coefficients",
     main = "ACTUAL AND ESTIMATED CHANNELS")

plot(
  H_power_esti_dB_LS[1, seq(1, ncol(H_power_esti_dB_LS), 8)],
  "b",
  xlab = "Time in samples",
  ylab = "Magnitude of coefficients",
  main = "ACTUAL AND ESTIMATED CHANNELS"
)

plot(
  H_power_esti_dB_MMSE[1, seq(1, ncol(H_power_esti_dB_MMSE), 8)],
  "b",
  xlab = "Time in samples",
  ylab = "Magnitude of coefficients",
  main = "ACTUAL AND ESTIMATED CHANNELS"
)

