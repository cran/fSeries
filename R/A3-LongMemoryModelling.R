# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA. 

# Copyrights (C)
# for this R-port: 
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTIONS:          FRACTIONAL GAUSSIAN NOISE:
#  fgnSim              Generates FGN Sequence after Durbin, Paxson, Beran
################################################################################


fgnSim = 
function(n = 1000, H = 0.7, 
method = c("beran", "durbin", "paxson"), mean = 0, std = 1)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a series of fractional Gaussian Noise
    
    # Arguments:
    #   n - length of desired series.
    #   H - the Hurst parameter.
    #   sigma - the standard deviation of the innovations used.
    
    # Details:
    #   FGN time series simulation. The FGN sequences use 
    #   functions going back to Taqqu et al., Paxson and 
    #   Beran (with Maechler's modifications from StatLib).

    # FUNCTION:
    
    # Settings:
    method = method[1]
    ans = NA
    
    # FROM DURBIN: 
    # Internal Function
    fgnSim.durbin = function(n, H, mean, std) {
        # Function to simulate a FGN  sequence, using the Durbin-Levinson 
        # coefficients, along the original C Function of Vadim Teverovsky
        # Original Cdurbin.C, double loop removed - now single loop in R!
        # This makes the function quite fast.
        n = n + 1
        h = H
        sigma = std
        normal = rnorm(n+1)
        sigma2 = sigma*sigma/2
        acov0 = 2*sigma2
        acov = vee = phi1 = phi2 = output = rep(0, n)
        I = 1:n
        acov = sigma2 *( (I+1)^(2*h) - 2*I^(2*h) + abs(I-1)^(2*h) )
        phi1[1] = acov[1]/acov0
        phi2[1] = phi1[1]
        vee0 = acov0
        vee[1] = vee0 * (1 - phi1[1]^2)
        output[1] = sqrt(vee0) * normal[1]      
        # Durbin-Levinson:
        for (i in 2:n){
            phi1[i] = acov[i]
            J = 1:(i-1)
            phi1[i] = phi1[i] - sum(phi2[J]*acov[i-J])
            phi1[i] = phi1[i]/vee[i-1]
            vee[i] = vee[i-1]*(1-phi1[i]^2)
            output[i] = sqrt(vee[i-1]) * normal[i]
            phi1[J] = phi2[J] - phi1[i]*phi2[i-J]
            output[i] = output[i] + sum(phi2[J] * output[i-J])
            phi2[1:i] =  phi1[1:i]}     
        # Return value:
        sigma*output[-1] + mean }
    
    # FROM PAXSON: 
    # Internal Function
    fgnSim.paxson = function(n, H, mean, std) {    
        # Description:
        #   Generates a FGN sequence by Paxson's FFT-Algorithm       
        # Details:
        #   This file contains a function for synthesizing approximate 
        #   fractional Gaussian noise.  
        #   * Note that the mean and variance of the returned points is 
        #     irrelevant, because any linear transformation applied to the 
        #     points preserves their correlational structure (and hence 
        #     their approximation to fractional Gaussian noise); and by 
        #     applying a linear transformation you can transform the 
        #     points to have any mean and variance you want.
        #   * If you're using the sample paths for simulating network 
        #     arrival counts, you'll want them to all be non-negative.  
        #     Hopefully you have some notion of the mean and variance 
        #     of your desired traffic, and can apply the corresponding
        #     transformation. If, after transforming, a fair number of 
        #     the points are still negative, then perhaps your traffic 
        #     is not well-modeled using fractional Gaussian noise.  
        #     You might instead consider using an exponential transformation.  
        #   * All of this is discussed in the technical report:
        #     Fast Approximation of Self-Similar Network Traffic,
        #     Vern Paxson, technical report LBL-36750/UC-405, April 1995.
        #     URL:ftp://ftp.ee.lbl.gov/papers/fast-approx-selfsim.ps.Z
        # Value:
        #   FGN vector of length n.
        # FUNCTION:
        # Internal Function: 
        FGN.B.est = function(lambda, H) {
            # Returns the estimate for B(lambda,H).
            d = -2*H - 1
            dprime = -2*H
            a = function(lambda,k) 2*k*pi+lambda
            b = function(lambda,k) 2*k*pi-lambda
            a1 = a(lambda,1); b1 = b(lambda,1)
            a2 = a(lambda,2); b2 = b(lambda,2)
            a3 = a(lambda,3); b3 = b(lambda,3)
            a4 = a(lambda,4); b4 = b(lambda,4)
            a1^d+b1^d+a2^d+b2^d+a3^d+b3^d + (a3^dprime+b3^dprime +
                a4^dprime+b4^ dprime)/(8*pi*H) }            
        # Internal Function:  
        FGN.spectrum = function(lambda, H) {
            # Returns an approximation of the power spectrum for 
            # fractional Gaussian noise at the given frequencies 
            # lambda and the given Hurst parameter H.
            2 * sin(pi*H) * gamma(2*H+1) * (1-cos(lambda)) * 
                (lambda^(-2*H-1) + FGN.B.est(lambda, H))}   
        # Returns a Fourier-generated sample path of a "self similar" process,
        # consisting of n points and Hurst parameter H (n should be even).
        n = n/2
        lambda = ((1:n)*pi)/n
        # Approximate ideal power spectrum.
        f = FGN.spectrum(lambda, H)
        # Adjust for estimating power
        # spectrum via periodogram.
        f = f * rexp(n)
        # Construct Corresponding Complex:
        # numbers with random phase.
        z = complex(modulus = sqrt(f), argument = 2*pi*runif(n))
        # Last element should have zero phase:
        z[n] = abs(z[n])
        # Expand z to correspond to a Fourier:
        # transform of a real-valued signal.
        zprime = c(0, z, Conj(rev(z)[-1]))
        # Inverse FFT gives sample path:
        z = Re(fft(zprime, inv = TRUE))     
        # Standardize:
        z = (z-mean(z))/sqrt(var(z))
        std*z + mean }
    
    # FROM BERAN: 
    # Internal Function
    fgnSim.beran = function(n, H, mean, std) { 
        # Description:
        #   Generates a FGN sequence by Beran's FFT-Algorithm  
        # Details:
        #   Uses the function simARMA0() from Beran's Long Memory 
        #   software package.
        # Value:
        #   FGN vector of length n.
        # FUNCTION:
        simFGN0 = function(n, H) {
            z = rnorm(2*n)
            zr = z[c(1:n)]
            zi = z[c((n+1):(2*n))]
            zic = -zi
            zi[1] = 0
            zr[1] = zr[1]*sqrt(2)
            zi[n] = 0
            zr[n] = zr[n]*sqrt(2)
            zr = c(zr[c(1:n)], zr[c((n-1):2)])
            zi = c(zi[c(1:n)], zic[c((n-1):2)])
            z = complex(real = zr,imaginary = zi)
            gksqrt = Re(gkFGN0(n, H))
            if (all(gksqrt > 0)) {
                gksqrt = sqrt(gksqrt)
                z = z*gksqrt
                z = fft(z, inverse = TRUE)
                z = 0.5*(n-1)**(-0.5)*z
                z = Re(z[c(1:n)]) } 
            else {
                gksqrt = 0*gksqrt
                cat("Re(gk)-vector not positive") }
            drop(z)}    
        gkFGN0 = function(n, H) {
            gammak = ckFGN0(n, H)
            ind = c(0:(n - 2), (n - 1), (n - 2):1)
            fft(c(gammak[ind+1]), inverse = TRUE) }
        ckFGN0 = function(n, H) {
            k = 0:(n-1)
            (abs(k-1)**(2*H)-2*abs(k)**(2*H)+abs(k+1)**(2*H))/2 }   
        # Generate Sequence:
        z = simFGN0(n = n, H = H)   
        # Standardize:
        # (z-mean(z))/sqrt(var(z))
        std*z + mean }
            
    # Generate Sequence:
    if (method == "beran")  
        ans = fgnSim.beran (n = n, H = H, mean = mean, std = std)
    if (method == "durbin") 
        ans = fgnSim.durbin(n = n, H = H, mean = mean, std = std)
    if (method == "paxson") 
        ans = fgnSim.paxson(n = n, H = H, mean = mean, std = std)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------
    
