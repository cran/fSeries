

C PART I: GARCH
C PART II: APARCH


C ##############################################################################


C ---------------------------------------------------------------------+
C  PART I: CONTENT: GARCH PACKAGE                                      ¦
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



C ##############################################################################



C ---------------------------------------------------------------------+
C  PART II: CONTENT: A-PARCH PACKAGE                                   ¦
C     aparchsim         Simulation of A-PARCH Processes                ¦
C     aparchllh         Evaluation of Log-Likelihood Function          ¦
C       gdensity          Normal Gaussian Densitz                      ¦
C       tdensity          Student t-distribution                       ¦
C       adensity          Symmetric alfa-stable distribution           ¦
C       gamma2, derfc2    Special Functions                            ¦
C ---------------------------------------------------------------------+


C **********************************************************************
C aparchsim:                                                           *
C Simulation of A-PARCH Processes                                      *
C **********************************************************************


      subroutine aparchsim( z, x, h, nt,
     &  omega, alfa, gamma, laga, na, beta, lagb, nb, delta)

C ----------------------------------------------------------------------
C    SIMULATION OF APARCH(p,q) TIME SERIES MODELS
C       *  x(t) = z(t) * h(t)**(1/delta)                                                                                   
C       *  h(t) = omega                              
C       *       + A(L) * [(|eps(t)|-gamma(t)*eps(t))**delta]   
C       *       + B(L) * h(t)                               
C     Arguments:
C       z(nt)          I  innovations for time series
C       x(nt)          O  time series points to be calculated
C       h(nt)          O  aparch volatilities
C         nt           I  number of data points
C       omega          I  value of constant
C       alfa(na)       I  subset of alfa coefficients
C         gamma(na)    I  assymetry coefficients
C         laga(na)     I  lags belonging to alpha and gamma
C         na           I  number of alfa coefficients 
C                          the maximum lag is the order np
C       beta(nb)       I  subset of beta coefficients
C         lagb(nb)     I  lags belonging to beta
C         nb           I  number of beta coefficients
C                          the maximum lag is the order nq
C       delta          I  value of recursion exponent
C ----------------------------------------------------------------------


C
C>>> DECLARATIONS: 
      implicit double precision (a-h, o-z)
      dimension z(nt), x(nt), h(nt)
      dimension alfa(na), gamma(na), laga(na), beta(nb), lagb(nb)
C
C>>> SETTINGS:
      np=laga(na)
      nq=lagb(nb)
      maxpq=np
      if(nq.gt.np) maxpq=nq
      deltainv=1.0d0/delta
C      
C>>> SET THE FIRST h(t) AND x(t):
      do n=1,maxpq
         h(n)=0.0d0
         x(n)=z(n)
      enddo
C
C>>> CONTINUE BY ITERATION:
      do n = maxpq+1, nt
         h(n)=omega
         do i = 1, na, 1
	        a = (abs(x(n-laga(i))) - gamma(i)*x(n-laga(i)))
            h(n) = h(n) + alfa(i) * a**delta
         enddo
         do j = 1, nb, 1
            h(n) = h(n) + beta(j) * h(n-lagb(j))
         enddo
         x(n)=h(n)**deltainv * z(n)
      enddo
C
      return
      end


C **********************************************************************
C  aparchllh                                                           *
C  Log-likelihood Function                                             *
C **********************************************************************

C        	result <- .Fortran("sumllh",
C    			as.integer(modsel),
C    			as.double(x),
C    			as.double(h),
C    			as.integer(nt),
C    			as.double(omega),
C    			as.double(alpha),
C			 	as.double(gamma),
C    			as.integer(laga),
C    			as.integer(na),
C				as.double(beta),
C    			as.integer(lagb),
C    			as.integer(nb),
C				as.double(delta),
C				as.double(disparm),
C    			as.integer(n.cond),
C    			as.double(0)
C    			)[[16]][
    			
       subroutine sumllh (x, h, nt,
     &  omega, alfa, gamma, laga, na, beta, lagb, nb, delta,
     &  ncond)   	
     
     C
C>>> DECLARATIONS:
      implicit  double precision (a-h, o-z)
      dimension x(nt), h(nt)
      dimension alfa(na), gamma(na), laga(na), beta(nb), lagb(nb)
C
C>>> ITERATE A-PARCH TIME SERIES PROCESS:
      do n = ncond+1, nt
         h(n) = omega
         do i = 1, na, 1
            a = (abs(x(n-laga(i))) - gamma(i)*x(n-laga(i)))           
            h(n) = h(n) + alfa(i) * abs(a)**delta
          enddo
          do j = 1, nb, 1
             h(n) = h(n) + beta(j) * h(n-lagb(j))
          enddo
      enddo    
C
      return
      end		
    			
    			
C ------------------------------------------------------------------------------



C **********************************************************************
C  aparchllh                                                           *
C  Log-likelihood Function                                             *
C **********************************************************************


      subroutine aparchllh (modsel, x, h, nt,
     &  omega, alfa, gamma, laga, na, beta, lagb, nb, delta, disparm,
     &  ncond, fllh)

C ----------------------------------------------------------------------
C    EVALUATION OF LOG-LIKELIHOOD FUNCTION FOR APARCH(p,q) MODELS
C       *  x(t) = z(t) * h(t)**(1/delta)                                                                                   
C       *  h(t) = omega                              
C       *       + A(L) * [(|eps(t)|-gamma(t)*eps(t))**delta]   
C       *       + B(L) * h(t)    
C     Arguments:
C       model          I  type of model assumption:    
C                           1: Gaussian innovations
C                           2: Student t-distributed innovations
C                           3: Symmetric alpha-stable innovations
C       x(nt)          I  time series points to be calculated
C       h(nt)          O  aparch volatilities
C         nt           I  number of data points
C       omega          I  value of constant
C       alfa(na)       I  subset of alfa coefficients
C         gamma(na)    I  assymetry coefficients
C         laga(na)     I  lags belonging to alpha and gamma
C         na           I  number of alfa coefficients 
C                           the maximum lag is the order np
C       beta(nb)       I  subset of beta coefficients
C         lagb(nb)     I  lags belonging to beta
C         nb           I  number of beta coefficients
C                           the maximum lag is the order nq
C       delta          I  value of volatility exponent
C       disparm        I  Distribution parameter
C                           model 2: df degress of freedom
C                           model 3: alpha stable index
C       ncond          I  conditining start value
C       fllh           O  value of log-likelihood function
C ----------------------------------------------------------------------

C
C>>> DECLARATIONS:
      implicit  double precision (a-h, o-z)
      dimension x(nt), h(nt)
      dimension alfa(na), gamma(na), laga(na), beta(nb), lagb(nb)
C
C>>> ITERATE A-PARCH TIME SERIES PROCESS:
      np=laga(na)
      nq=lagb(nb)
      maxpq = np
      if(nq.gt.np) maxpq = nq
      deltainv=1.0/delta
      do n = maxpq+1, nt
         h(n)=omega
         do i = 1, na, 1
            a = (abs(x(n-laga(i)))-gamma(i)*x(n-laga(i)))           
            h(n) = h(n) + alfa(i) * abs(a)**delta
          enddo
          do j = 1, nb, 1
             h(n) = h(n) + beta(j) * h(n-lagb(j))
          enddo
      enddo    
C
C>>> LOG-LIKELIHOOD FUNCTION -- GAUSSIAN DENSITY:
      f=0.0d0
      do n=ncond+1,nt
         hh=abs(h(n))**deltainv
         z=x(n)/hh
         f=f-dlog(density(modsel, z, disparm)/hh)
      enddo
      fllh=f/(nt-ncond)
C
      return
      end
      

C **********************************************************************
C density:                                                             *
C selector for density function                                        *
C **********************************************************************


      double precision function density(model, x, disparm)
      implicit double precision (a-h, o-z)
      density=0.0d0
      if (model.eq.1) density=gdensity(x)
      if (model.eq.2) density=tdensity(x, disparm)
      if (model.eq.3) density=adensity(x, disparm)
      return
      end


C **********************************************************************
C gdensity:                                                            *
C normal (Gaussian Distribution)                                       *
C **********************************************************************

    
      double precision function gdensity(x)
      implicit double precision (a-h, o-z)
      pi=4.0d0*datan(1.0d0)
      gdensity=(1.0d0/dsqrt(2.0d0*pi))*dexp(-x*x/2.0d0)
      return
      end


C **********************************************************************
C tdensity:                                                            *
C Student t-distribution function                                      *
C **********************************************************************


      double precision function tdensity(x, df)
      implicit double precision (a-h, o-z)    
      tdensity=1
      return
      end


C **********************************************************************
C adensity:                                                            *
C symmetric alfa-stable distribution function                          *
C **********************************************************************

      
      double precision function adensity(x, alpha)
      implicit double precision (a-h, o-z)
      integer first                          
      dimension ei(3),u(3)                  
      dimension s(4,19),r(19),q(0:5),p(0:5,0:19)
      dimension pd(0:4,0:19)                
      dimension alf2i(4)                     
      dimension znot(19),zn4(19),zn5(19),combo(0:5),zji(19,0:5)      
      SAVE first,pi,sqpi,a2,cpxp0,gpxp0,cpxpp0,gpxpp0
      SAVE cppp,gppp,znot,zn4,zn5,zji,ei,u,q
      SAVE a,ala,alf2,alf1,pialf,sp0,sppp0,xp0
      SAVE xpp0,xppp0,spzp1,rp0,rpp0,rppp0,rp1
      SAVE alf2i,r,b,c,p,pd,oldalf
      data combo/1.d0,5.d0,10.d0,10.d0,5.d0,1.d0/      
      data ((s(i,j),i=1,4),j=1,7)/
     1    1.85141 90959 d2,  -4.67693 32663 d2,    
     1    4.84247 20302 d2,  -1.76391 53404 d2,    
     1   -3.02365 52164 d2,   7.63519 31975 d2,    
     1   -7.85603 42101 d2,   2.84263 13374 d2,    
     1    4.40789 23600 d2,  -1.11811 38121 d3,    
     1    1.15483 11335 d3,  -4.19696 66223 d2,    
     1   -5.24481 42165 d2,   1.32244 87717 d3,    
     1   -1.35556 48053 d3,   4.88340 79950 d2,    
     1    5.35304 35018 d2,  -1.33745 70340 d3,    
     1    1.36601 40118 d3,  -4.92860 99583 d2,    
     1   -4.89889 57866 d2,   1.20914 18165 d3,    
     1   -1.22858 72257 d3,   4.40631 74114 d2,    
     1    3.29055 28742 d2,  -7.32117 67697 d2,    
     1    6.81836 41829 d2,  -2.28242 91084 d2/   
      data ((s(i,j),i=1,4),j=8,14)/
     1   -2.14954 02244 d2,   3.96949 06604 d2,    
     1   -3.36957 10692 d2,   1.09058 55709 d2,    
     1    2.11125 81866 d2,  -2.79211 07017 d2,    
     1    1.17179 66020 d2,   3.43946 64342 d0,    
     1   -2.64867 98043 d2,   1.19990 93707 d2,    
     1    2.10448 41328 d2,  -1.51108 81541 d2,    
     1    9.41057 84123 d2,  -1.72219 88478 d3,    
     1    1.40875 44698 d3,  -4.24725 11892 d2,    
     1   -2.19904 75933 d3,   4.26377 20422 d3,    
     1   -3.47239 81786 d3,   1.01743 73627 d3,    
     1    3.10474 90290 d3,  -5.42042 10990 d3,    
     1    4.22210 52925 d3,  -1.23459 71177 d3,    
     1   -5.14082 60668 d3,   1.10902 64364 d4,    
     1   -1.02703 37246 d4,   3.42434 49595 d3/   
      data ((s(i,j),i=1,4),j=15,19)/      
     1    1.12151 57876 d4,  -2.42435 29825 d4,    
     1    2.15360 57267 d4,  -6.84909 96103 d3,    
     1   -1.81206 31586 d4,   3.14301 32257 d4,    
     1   -2.41642 85641 d4,   6.91268 62826 d3,    
     1    1.73884 13126 d4,  -2.21083 97686 d4,    
     1    1.33979 99271 d4,  -3.12466 11987 d3,    
     1   -7.24357 75303 d3,   4.35453 99418 d3,    
     1    2.36161 55949 d2,  -7.65716 53073 d2,    
     1   -8.73767 25439 d3,   1.55108 52129 d4,    
     1   -1.37897 64138 d4,   4.63874 17712 d3 /  
c     
      xfun(z,alpha,a)=((1.0d0-z)**(-1.0d0/alpha)-1.0d0)/a
      zpfun(x,alpha,a)=alpha*a*(1+a*x)**(-alpha-1.0d0)  
      cfun(x)=0.5d0-datan(x)/pi            
      cden(x)=1.0d0/(pi*(1.0d0+x*x))             
      gfun(x)=0.5d0*derfc2(x/2.0d0)              
      gden(x)=exp(max(-160.d0,-x*x/4.d0))/(2.0d0*sqpi) 
c                                            
c  except on the first call,skip to 10.     
      if (first.eq.0)  then            
        pi=3.14159 26535 89793 d0            
        sqpi=sqrt(pi)                        
        a2=dsqrt(2.d0)-1                   
        cpxp0=1.0d0/pi                         
        gpxp0=1.0d0/(4.0d0*a2*sqpi)            
        cpxpp0=cpxp0*2.0d0                     
        gpxpp0=gpxp0*1.5d0                   
        cppp=cpxpp0*3.0d0-2.0d0/pi             
        gppp=gpxpp0*2.5d0-1.0d0/(32.0d0*sqpi*a2**3)
        do j=1,19                         
           znot(j)=dfloat(j)*.05d0              
           zn4(j)=(1-znot(j))**4              
           zn5(j)=(1-znot(j))*zn4(j)        
           do i=0,5                          
              zji(j,i)=combo(i)*(-znot(j))**(5-i)     
           enddo
        enddo
        do i=1,3                          
           ei(i)=i                              
           u(i)=1                               
        enddo                                  
        q(0)=0                               
        first=1
      endif
c                                            
c  if alpha is unchanged from last call,skip to 20  
      if (alpha.ne.oldalf) then
        a=2.0d0**(1.0d0/alpha)-1.0d0                 
        ala=alpha*a                        
        alf2=2.0d0-alpha                       
        alf1=alpha-1.0d0                       
        pialf=pi*alpha                     
        sp0=dgamma2(1.0d0/alpha)/pialf           
        sppp0=-dgamma2(3.0d0/alpha)/pialf        
        xp0=1.d0/ala                          
        xpp0=xp0*(1.0d0+alpha)/alpha         
        xppp0=xpp0*(1+2.0d0*alpha)/alpha 
        spzp1=(a**alpha)*dgamma2(alpha)*sin(pialf/2)/pi   
        rp0=-sp0*xp0+alf2*cpxp0+alf1*gpxp0 
        rpp0=-sp0*xpp0+alf2*cpxpp0+alf1*gpxpp0
        rppp0=-sp0*xppp0-sppp0*xp0**3
     &     +alf2*cppp+alf1*gppp                          
        rp1=-spzp1+alf2/pi               
        do i=1,4                         
           alf2i(i)=alf2**i-1               
        enddo
        do j=1,19                        
           r(j)=alf2*prod(s(1,j),alf2i,4)  
        enddo
        q(1)=rp0                             
        q(2)=rpp0/2                          
        q(3)=rppp0/6                         
        b=- prod(u,q(1),3)-prod(r,zn5,19)      
        c=rp1-prod(ei,q(1),3)-5.0d0*prod(r,zn4,19)
        q(4)=5*b-c                       
        q(5)=b-q(4)                        
        do i=0,5                         
           p(i,0)=q(i)                         
           do j=1,19                        
              p(i,j)=q(i)+prod(r,zji(1,i),j) 
           enddo
        enddo
        do i=1,5                         
           do j=0,19                        
              pd(i-1,j)=i*p(i,j)               
           enddo
        enddo
        oldalf=alpha                         
      endif
c                                            
c  now calculate new x:
      xa1=1.0d0+a*dabs(x)                   
      xa1a=xa1**(-alpha)                 
      z=1.0d0-xa1a                           
      zp=ala*xa1a/xa1                  
      x1=xfun(z,1.d0,1.d0)               
      x2=xfun(z,2.d0,a2)                 
      x1p=1/zpfun (x1,1.d0,1.d0)       
      x2p=1/zpfun (x2,2.d0,a2)         
      j=20*z                               
      j=min(j,19)                         
      rz=poly(p(0,j),z,5)               
      rpz=poly(pd(0,j),z,4)             
40    continue                               
      sfun=alf2*cfun(x1)+alf1*gfun(x2)+rz  
      sden=(alf2*cden(x1)*x1p+alf1*gden(x2)*x2p-rpz)*zp 
      if (x.lt.0.0d0) sfun=1.0d0-sfun            
      adensity=sden
c
      return                                 
      end                                    
c                                             
      function poly(a,x,k)               
c  computes k degree polynomial in x with coefficients a(j)
      implicit real*8 (a-h,o-z)           
      dimension a(0: k)                      
      poly=a(k)                            
      do 10 j=k-1,0,-1                   
10    poly=poly*x+a(j)                 
      return                                 
      end                                    
c                                             
      function prod(a,b,k)                 
c  computes inner product of two k-vectors a and b
      implicit real*8 (a-h,o-z)           
      dimension a(k),b(k)                   
      prod=0                               
      do 10 j=1,k                         
10    prod=prod+a(j)*b(j)              
      return                                 
      end                                    


C **********************************************************************
C dgamma2:                                                             *
C Gamma Function                                                       *
C **********************************************************************

   
      double precision function dgamma2(xx)
      implicit real*8 (a-h,o-z)
      dimension cof(6)
      data cof,stp  /76.18009173d0,-86.50532033d0,24.01409822d0,
     &-1.231739516d0,0.120858003d-2,-0.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*dlog(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      dgamma2=dexp(tmp+dlog(stp*ser))
      return
      end


C **********************************************************************
C derfc2:                                                              *
C Complementary Error Function                                         *
C **********************************************************************

     
      double precision function derfc2(x)
      double precision x,result
      call calerf(x,result,1)
      derfc2=result
      return
      end

      subroutine calerf(arg,result,jint)
      integer i,jint
      double precision
     1     a,arg,b,c,d,del,four,half,p,one,q,result,sixten,sqrpi,
     2     two,thresh,x,xbig,xden,xhuge,xinf,xmax,xneg,xnum,xsmall,
     3     y,ysq,zero
      dimension a(5),b(4),c(9),d(8),p(6),q(5)
C
C>>> Mathematical constants
      data four,one,half,two,zero/4.0d0,1.0d0,0.5d0,2.0d0,0.0d0/,
     1     sqrpi/5.6418958354775628695d-1/,thresh/0.46875d0/,
     2     sixten/16.0d0/
C
C>>> Machine-dependent Constants
      data xinf,xneg,xsmall/1.79d308,-26.628d0,1.11d-16/,
     1     xbig,xhuge,xmax/26.543d0,6.71d7,2.53d307/
C
C>>> Coefficients for approximation to erf in first interval
      data a/3.16112374387056560d00,1.13864154151050156d02,
     1       3.77485237685302021d02,3.20937758913846947d03,
     2       1.85777706184603153d-1/
      data b/2.36012909523441209d01,2.44024637934444173d02,
     1       1.28261652607737228d03,2.84423683343917062d03/
C
C>>> Coefficients for approximation to erfc in second interval
      data c/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/
      data d/1.57449261107098347D01,1.17693950891312499D02,
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/
C
C>>> Coefficients for approximation to erfc in third interval
      data p/3.05326634961232344d-1,3.60344899949804439d-1,
     1       1.25781726111229246d-1,1.60837851487422766d-2,
     2       6.58749161529837803d-4,1.63153871373020978d-2/
      data q/2.56852019228982242d00,1.87295284992346047d00,
     1       5.27905102951428412d-1,6.05183413124413191d-2,
     2       2.33520497626869185d-3/
      x=arg
      y=dabs(x)
      if (y.le.thresh) then
C
C>>> Evaluate erf for |x| <= 0.46875
        ysq=zero
        if (y.gt.xsmall) ysq=y*y
        xnum=a(5)*ysq
        xden=ysq
        do i=1, 3
          xnum=(xnum+a(i))*ysq
          xden=(xden+b(i))*ysq
        enddo    
        result=x*(xnum+a(4))/(xden+b(4))
        if (jint.ne.0) result=one-result
        if (jint.eq.2) result=dexp(ysq)*result
        return     
C
C>>> Evaluate erfc for 0.46875 <= |x| <= 4.0
      elseif (y.le.four) then
        xnum=c(9)*y
        xden=y
        do i=1, 7
          xnum=(xnum+c(i))*y
          xden=(xden+d(i))*y
        enddo     
        result=(xnum+c(8))/(xden+d(8))
        if (jint.ne.2) then
          ysq=aint(y*sixten)/sixten
          del=(y-ysq)*(y+ysq)
          result=dexp(-ysq*ysq)*dexp(-del)*result
        endif
C
C>>> Evaluate erfc for |x| > 4.0
      else
        result=zero
        if (y.ge.xbig) then
          if ((jint.ne.2).or.(y.ge.xmax)) goto 300
          if (y.ge.xhuge) then
            result=sqrpi/y
            goto 300
          endif
        endif
        ysq=one/(y*y)
        xnum=p(6)*ysq
        xden=ysq
        do i=1, 4
          xnum=(xnum+p(i))*ysq
          xden=(xden+q(i))*ysq
        enddo    
        result=ysq *(xnum+p(5))/(xden+q(5))
        result=(sqrpi-result)/y
        if (jint.ne.2) then
          ysq=aint(y*sixten)/sixten
          del=(y-ysq)*(y+ysq)
          result=dexp(-ysq*ysq)*dexp(-del)*result
        endif
      endif
C
C>>> Fix up for negative argument, erf, etc.
300   if (jint.eq.0) then
        result=(half-result)+half
        if (x.lt.zero) result=-result
      elseif (jint.eq.1) then
        if (x.lt.zero) result=two-result
      else
        if (x.lt.zero) then
          if (x.lt.xneg) then
            result=xinf
          else
            ysq=aint(x*sixten)/sixten
            del=(x-ysq)*(x+ysq)
            y=dexp(ysq*ysq)*dexp(del)
            result=(y+y)-result
          endif
        endif
      endif
      return
      end

      
C ##############################################################################
