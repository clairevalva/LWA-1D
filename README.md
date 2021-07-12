# LWA-1D
1D LWA code for Valva and Nakamura (2021) which can be found [here](https://doi.org/10.1029/2020JD034501)

==========================================

lwa1d.f90:  A Fortran 90 code to compute LWA in 1D traffic flow model per Valva and Nakamura (2021). This model uses Intel's MKL Vector Statistics Library to compute random numbers and netCDF to read forcing spectra data 'forcing.nc’.

forcing.nc:  A file containing the amplitude of the spectral coefficients for the sum III+IV+V in Eq. (6) of Valva and Nakamura (2021) for DJF 40-60S.  Based on 1979-2018 ERA-Interim 6-hourly data.  Variable name: ‘forcing’  Variable dimension: 41 x 61.  The first dimension spans zonal wavenumbers from -20 to 20.  The second dimension spans frequency from 0 to 1.491 CPD.  The negative wavenumber is equivalent to negative frequency. To match the observed variance of III+IV+V, multiply by 350.  (The mapping is performed in lwa1d.f90.)  

==========================================
