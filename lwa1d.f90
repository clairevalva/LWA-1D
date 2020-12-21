! ##### A Fortran 90 code to compute LWA in 1D traffic flow model per 
! ##### Valva and Nakamura (2021). This model uses intel's MKL Vector 
! ##### Statistics Library to compute random numbers and 
! ##### netCDF to read forcing spectra data 'forcing.nc'

include 'mkl_vsl.f90' 
       program MKL_VSL_GAUSSIAN
       USE MKL_VSL_TYPE 
       USE MKL_VSL
       USE NETCDF
      integer, parameter :: imax = 1001,mend=600001,mind = 180 
      integer, parameter :: mhalf = 1
      common /array/  aa(imax+1,4),cc(imax+1),a0(imax+1)
      common /brray/  f1(imax+1,3),f2(imax+1,3),f3(imax+1,3)
      common /crray/  d1(imax+1,3),d2(imax+1,3),d3(imax+1,3)
      common /drray/  aad(imax+1),bb(imax+1)
      real(kind=4) r1(1000),r2(1000),del(2460),ff(imax+1) ! 
      real(kind=4) s ! 
      real(kind=4) a, b, sigma !
      TYPE (VSL_STREAM_STATE) :: stream
      integer(kind=4) errcode 
      integer(kind=4) i,j 
      integer brng,method,seed,n
      real :: dt,dx,rkappa,pi,alpha,forcing(60,41),rpn(41,60)
      integer :: ncid, status,nDim,nVar,nAtt,uDimID,inq
      integer :: lonID,latID,vid2,varID

      isw = 0   ! 1: stationary wave forcing
                ! if 1, set mhalf = 72000 and mend = 672001
                ! otherwise mhalf = 1 and mhalf = 600001 
      gama = 15.          ! stationary wave amplitude (m/s)
      if(isw.ne.1) gama = 0.

      seed = 555    ! Randomization seed. Change every time to avoid identical
                    ! results
      ntape = 21    ! (output unit - change to avoid overwriting previous 
                    ! output) 

      pi = acos(-1.)
      dt = 120.
      dx = 28000. 
      alpha = 0.4         ! strength of nonlinearity
      uj = 40.            ! jet speed (m/s)
      a00 = 26.           ! mean LWA (m/s)
      eps = 0.86          ! strength of transient eddy forcing (obs: 1)
      amp = 350.*eps      ! = 350 is the default (same as obs)
      s00 = 3.e-5         ! constant LWA forcing (m/s**2)
      rkappa = 0.1*dx*dx/dt  ! diffusion coefficient
      write(6,*) 'rkappa =',rkappa

      tau = 10.*24.*3600.     ! damping time (sec) 
      dk = 2.*pi/28000000.   ! eddy forcing wavenumber increment (1/m)
      dw = 2.*pi/(179.*12.*3600.)   ! eddy forcing frequency increment (1/s)
!     n = 1000

!  ***** Read in eddy forcing spectra ******
       stat = nf90_open('forcing.nc',nf90_nowrite,ncid)
       stat = nf90_inquire(ncid,nDim,nVar,nAtt,uDimID)
       write(6,*) 'ndim,nvar,natt,uDimID =',nDim,nVar,nAtt,uDimID
       stat = nf90_inq_varid(ncid,"forcing",varID)
       write(6,*) 'Variable ID for forcing = ',varID
       stat = nf90_get_var(ncid,varID,forcing)
       stat = nf90_close(ncid)

      write(6,*) maxval(forcing),minval(forcing)

!  **** Calculate random phase *****
      brng=VSL_BRNG_MT19937 
      method = VSL_RNG_METHOD_UNIFORM_STD
      errcode=vslnewstream( stream, brng, seed )
      errcode=vsrnguniform( method, stream, 2460, del, 0.,1. ) 
      ii = 1
      do ik = 1,41
      do io = 1,60
       rpn(ik,io) = del(ii)
       ii = ii+1 
      enddo
      enddo

      rpn(:,:) = rpn(:,:)*2.*pi
      write(6,*) maxval(rpn),minval(rpn)

! ***** Initialization *****
      do i = 1,imax+1
        p = 2.*pi*float(i-1)/float(imax-1)
        a0(i) = gama*cos(2.*p)      ! STATIONARY WAVE
        aa(i,:) = a00
        if(isw.eq.1) aa(i,:) = 0.
      enddo

      f1(:,:) = 0.
      f2(:,:) = 0.
      f3(:,:) = 0.
      d1(:,:) = 0.
      d2(:,:) = 0.
      d3(:,:) = 0.

! ***** Main loop ****
        nfile = 0
      do m = 1,mend 
        t = dt*float(m-1)
        do i = 2,imax
         x = dx*float(i-1)
         f1(i,3) = (2.*alpha*a0(i+1)*aa(i+1,3) - 2.*alpha*  &
              a0(i-1)*aa(i-1,3))/(2.*dx) 
         if(isw.ne.1) f1(i,3) = 0.
         f2(i,3) = (alpha*aa(i+1,3)**2-alpha*     &
            aa(i-1,3)**2)/(2.*dx)
         f3(i,3) = -(uj*aa(i+1,3)-uj*      &
              aa(i-1,3))/(2.*dx)
         d1(i,3) = rkappa*(aa(i+1,1)+aa(i-1,1)-2.*aa(i,1))/(dx*dx)
         d2(i,3) = s00-aa(i,2)/tau 
         if(isw.ne.1) d2(i,3) = 0.
         

!! *** Transient eddy forcing calculations 
         d3(i,3) = 0.
      if(m.ge.mhalf) then
       do ik = 1,41
         rk = dk*(float(ik-1)-20.) 
        do io = 2,60
           w = dw*(float(io-1))
          if(ik.ne.21) then
           d3(i,3) = d3(i,3)+forcing(io,ik)*  &
                amp*cos(rk*x-w*t+rpn(ik,io)) 
          endif
        enddo
       enddo
       endif

!! *** 3rd-order Adams-Bashforth scheme ***
         
         f11 = (23.*f1(i,3)-16.*f1(i,2)+5.*f1(i,1))/12.
         f21 = (23.*f2(i,3)-16.*f2(i,2)+5.*f2(i,1))/12.
         f31 = (23.*f3(i,3)-16.*f3(i,2)+5.*f3(i,1))/12.
         d11 = (23.*d1(i,3)-16.*d1(i,2)+5.*d1(i,1))/12.
         d21 = (23.*d2(i,3)-16.*d2(i,2)+5.*d2(i,1))/12.
         d31 = (23.*d3(i,3)-16.*d3(i,2)+5.*d3(i,1))/12.
         aa(i,4) = aa(i,3)+dt*(f11+f21+f31+d11+d21+d31)
        enddo

!! *** Cyclic boundary condition ***
        aa(1,4) = aa(imax,4)
        aa(imax+1,4) = aa(2,4)

        write(6,*) 'End m =',m,m*120./(3600.*24.),&
               maxval(aa),minval(aa),nfile
  
!! *** Forward in time ***
        aa(:,1) = aa(:,2)
        aa(:,2) = aa(:,3)
        aa(:,3) = aa(:,4)
        f1(:,1) = f1(:,2)
        f1(:,2) = f1(:,3)
        f2(:,1) = f2(:,2)
        f2(:,2) = f2(:,3)
        f3(:,1) = f3(:,2)
        f3(:,2) = f3(:,3)
        d1(:,1) = d1(:,2)
        d1(:,2) = d1(:,3)
        d2(:,1) = d2(:,2)
        d2(:,2) = d2(:,3)
        d3(:,1) = d3(:,2)
        d3(:,2) = d3(:,3)

!! *** I/O (write sequentially to a binary file) ***
        if(mod(m,mind).eq.1.and.m.gt.mhalf) then
          aad(:) = aa(:,4)
          !aad(:) = d3(:,3)  ! forcing
          write(ntape) aad          
          nfile = nfile + 1
        endif
      enddo

      close(ntape)

      stop
      end 
