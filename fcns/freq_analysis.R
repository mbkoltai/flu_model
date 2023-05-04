# Fourier analysis / wavelet analysis

# install.packages("biwavelet")
# wavelet analysis
library("biwavelet")

# Generate some sample data
x <- seq(0, 2*pi, length.out = 1024)
y <- sin(2*x) + 0.5*sin(4*x)

# Perform biwavelet analysis
bwt <- bwt(y)

# Plot the biwavelet coefficients
biplot(bwt, xlab = "Time", ylab = "Frequency")

# Fourier analysis
# Generate some sample data
x <- seq(0, 2*pi, length.out = 1024)
y <- sin(2*x) + 0.5*sin(4*x)
plot(x,y,type="l")

# Compute the Fourier coefficients
y_fft <- fft(y)
y_coeffs <- abs(y_fft)

# Plot the Fourier coefficients
plot(y_coeffs, type = "l", xlab = "Frequency", ylab = "Magnitude")

