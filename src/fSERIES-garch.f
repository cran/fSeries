C ---------------------------------------------------------------------+
C  CONTENT: GARCH PACKAGE                                             ¦
C     garchsim          Process Simulation                             ¦
C     garchcvs          Conditional Variances                          ¦
C ---------------------------------------------------------------------+


C **********************************************************************
C garchsim:                                                            *
C Process Simulation                                                   *
C **********************************************************************


      subroutine garchsim(x, h, z, nt,
     &  omega, alpha, laga, na, beta, lagb, nb, h0)

C ----------------------------------------------------------------------
C    SIMULATION OF GARCH(p,q) TIME SERIES MODELS
C       *  x(t) = z(t) * h(t)**(1/2)         
C       *  h(t) = omega + A(L) * x(t)**2 + B(L) * h(t)                               
C     Arguments:
C       x(nt)          O  time series points to be calculated
C       h(nt)          O  garch volatilities
C       z(nt)          I  innovations for time series
C         nt           I  number of data points
C       omega          I  value of constant
C       alpha(na)      I  subset of alpha coefficients
C         laga(na)     I  lags belonging to alpha
C         na           I  number of alpha coefficients 
C                          the maximum lag is the order np
C       beta(nb)       I  subset of beta coefficients
C         lagb(nb)     I  lags belonging to beta
C         nb           I  number of beta coefficients
C                          the maximum lag is the order nq
C ----------------------------------------------------------------------


C
C>>> DECLARATIONS: 
      implicit double precision (a-h, o-z)
      dimension x(nt), h(nt), z(nt)
      dimension alpha(na), laga(na), beta(nb), lagb(nb)
C
C>>> SETTINGS:
      np=laga(na)
      nq=lagb(nb)
      maxpq=np
      if(nq.gt.np) maxpq=nq
C      
C>>> SET THE FIRST h(t) AND x(t):
      do i=1,maxpq
         h(i) = h0
         x(i) = z(i)
      enddo
C
C>>> CONTINUE BY ITERATION:
      do i=maxpq+1, nt
         h(i) = omega
         do j=1, na, 1
            h(i) = h(i) + alpha(j)*(x(i-laga(j)))**2 
         enddo
         do k=1, nq, 1
            h(i) = h(i) + beta(k)*h(i-lagb(k))
         enddo
         x(i)= dsqrt(h(i)) * z(i)
      enddo
C
      return
      end


C **********************************************************************
C  garchcvs                                                            *
C  Conditional Variances                                               *
C **********************************************************************


      subroutine garchcvs(x, h, z, nt, omega, alpha, laga, na, 
     &   beta, lagb, nb, h0)

C ----------------------------------------------------------------------
C    EVALUATION OF LOG-LIKELIHOOD FUNCTION FOR APARCH(p,q) MODELS
C       *  x(t) = z(t) * h(t)**(1/2)         
C       *  h(t) = omega + A(L) * x(t)**2 + B(L) * h(t)
C     Arguments:
C       x(nt)          I  time series points to be calculated
C       h(nt)          O  garch volatilities
C       z(nt)          I  innovations for time series
C         nt           I  number of data points
C       omega          I  value of constant
C       alpha(na)      I  subset of alpha coefficients
C         laga(na)     I  lags belonging to alpha and gamma
C         na           I  number of alpha coefficients 
C                           the maximum lag is the order np
C       beta(nb)       I  subset of beta coefficients
C         lagb(nb)     I  lags belonging to beta
C         nb           I  number of beta coefficients
C                           the maximum lag is the order nq
C ----------------------------------------------------------------------

C
C>>> DECLARATIONS:
      implicit  double precision (a-h, o-z)
      dimension x(nt), h(nt), z(nt)
      dimension alpha(na), laga(na), beta(nb), lagb(nb)
C
C>>> SETTINGS:
      np=laga(na)
      nq=lagb(nb)
      maxpq = np
      if(nq.gt.np) maxpq = nq

C>>> SET THE FIRST h(t):
      do i=1,maxpq
         h(i) = h0
	 z(i) = x(i)
      enddo
C
C>>> CONDITIONAL VARIANCES:
      do i=maxpq+1, nt
         h(i) = omega
         do j=1, na, 1
            h(i) = h(i) + alpha(j) * x(i-laga(j))**2
          enddo
          do k=1, nb, 1
             h(i) = h(i) + beta(k) * h(i-lagb(k))
          enddo
	  z(i) = x(i)/sqrt(h(i))
      enddo
C
      return
      end


C>>> LOG-LIKELIHOOD FUNCTION -- GAUSSIAN DENSITY:
C      f=0.0d0
C      do n=ncond+1,nt
C         hh=h(n)**deltainv
C         z=x(n)/hh
C         f=f-dlog(density(model,z, disparm)/hh)
C      enddo
C      fllh=f/(nt-ncond)
