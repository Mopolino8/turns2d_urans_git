module interpolators
  use input
  use grid
  use connectivity
  implicit none

  contains

!***********************************************************************
  subroutine spatial_interpolator_bicubic(q,nv)
!***********************************************************************
  integer :: nv
  real :: q(jmc,kmc,nv,ngrids)
  integer bcdim,ii,jj,iim,jjm,iip,jjp,kkp,id,j,k,n,im
  real djm1,dj0,djp1,dkm1,dk0,dkp1
  real w1,w2,w3,w4,w5,w6,w7,w8,w9,onefourth
  real,allocatable :: qbc(:,:)

  real a1,a2,a3,a4,a5,a6,a7,a8,a9
  logical :: loworder

  onefourth = 1./4

  bcdim=iieptr(ngrids)
  allocate(qbc(bcdim,nv))

!...LOOP THROUGH ALL THE MESHES AND COLLECT GLOBAL QBC

  do im = 1,ngrids
    do id=1,ndonor(im)

      ii=idonor(id,1,im)
      jj=idonor(id,2,im)

      iim=ii-1
      jjm=jj-1

      iip=ii+1
      jjp=jj+1

      djm1 = 1.+frac(id,1,im)
      dkm1 = 1.+frac(id,2,im)
      dj0  = frac(id,1,im)
      dk0  = frac(id,2,im)
      djp1 = 1.-frac(id,1,im)
      dkp1 = 1.-frac(id,2,im)

      w1   =    (dj0  * djp1 * dk0  * dkp1)*onefourth
      w2   = -2*(djm1 * djp1 * dk0  * dkp1)*onefourth
      w3   = -  (djm1 * dj0  * dk0  * dkp1)*onefourth
      w4   = -2*(dj0  * djp1 * dkm1 * dkp1)*onefourth
      w5   =  4*(djm1 * djp1 * dkm1 * dkp1)*onefourth
      w6   =  2*(djm1 * dj0  * dkm1 * dkp1)*onefourth
      w7   = -  (dj0  * djp1 * dkm1 * dk0 )*onefourth
      w8   =  2*(djm1 * djp1 * dkm1 * dk0 )*onefourth
      w9   =    (djm1 * dj0  * dkm1 * dk0 )*onefourth

      loworder = .false.
      lowordercheck: if (loworder) then

        a1 = dj0*dk0
        a2 = dk0 
        a3 = dj0*dk0
        a4 = dj0 
        a5 = 1.0
        a6 = dj0 
        a7 = dj0*dk0 
        a8 = dk0 
        a9 = dj0*dk0

        if (frac(id,1,im).gt.0) then
          w1 = 0.
          w4 = 0.
          w7 = 0.
          if (frac(id,2,im).gt.0) then
            w1 = 0.
            w2 = 0.
            w3 = 0.
            w5 = abs(djp1*dkp1)
            w6 = abs(a6*dkp1)
            w8 = abs(a8*djp1)
            w9 = abs(a9)
          else
            w7 = 0.
            w8 = 0.
            w9 = 0.
            w2 = abs(a2*djp1)
            w3 = abs(a3)
            w5 = abs(djp1*dkm1)
            w6 = abs(a6*dkm1)
          endif
        else
          w3 = 0.
          w6 = 0.
          w9 = 0.
          if (frac(id,2,im).gt.0) then
            w1 = 0.
            w2 = 0.
            w3 = 0.
            w4 = abs(a4*dkp1)
            w5 = abs(djm1*dkp1)
            w7 = abs(a7)
            w8 = abs(a8*djm1)
          else
            w7 = 0.
            w8 = 0.
            w9 = 0.
            w1 = abs(a1)
            w2 = abs(a2*djm1)
            w4 = abs(a4*dkm1)
            w5 = abs(djm1*dkm1)
          endif
        endif

      endif lowordercheck

!.....collect in global pointer qbc from pointer iisptr->iieptr

      do n=1,nv
        qbc(iisptr(im)-1+id,n) =                         &
                                 w1*q(iim,jjm,n,im)      &
                               + w2*q(ii ,jjm,n,im)      &
                               + w3*q(iip,jjm,n,im)      &
                               + w4*q(iim,jj ,n,im)      &
                               + w5*q(ii ,jj ,n,im)      &
                               + w6*q(iip,jj ,n,im)      &
                               + w7*q(iim,jjp,n,im)      &
                               + w8*q(ii ,jjp,n,im)      &
                               + w9*q(iip,jjp,n,im)
      enddo
    enddo

  enddo

!.....overwrite fringe points solution w/ donor global qbc solution
  do im = 1,ngrids
    do id=1,nfringe(im)
      j = ibc(id,im)
      ii = imesh(id,1,im)
      jj = imesh(id,2,im)

      do n = 1,nv
         q(ii,jj,n,im) = qbc(j,n)
      enddo

    enddo
  enddo
       
  deallocate(qbc)

  end subroutine spatial_interpolator_bicubic

!***********************************************************************
  subroutine spatial_interpolator_muscl(q,nv)
!***********************************************************************
  integer :: nv
  real :: q(jmc,kmc,nv,ngrids)

  integer :: j,k,jsf,ksf,jef,kef,jfine,kfine,n
  real :: f2i,f2i1,a1,a2,f3i,epsj,epsk

  epsj   =  (10./real(jmax(1)))**3
  epsk   =  (10./real(kmax(1)))**3

  do k = kbeg(1),kend(1)
    ksf = 2*k - kbeg(2)
    kef = 2*k - kbeg(2) + 1
    do j = jbeg(1),jend(1)
      jsf = 2*j - jbeg(2)
      jef = 2*j - jbeg(2) + 1

      do n = 1,nv
        do kfine = ksf,kef
          do jfine = jsf,jef
            q(jfine,kfine,n,2) = q(j,k,n,1)
          enddo
        enddo

        do kfine = ksf,kef
!          if (jsf.eq.jbeg(2).or.jef.eq.jend(2)) cycle
          f2i    = q(j+1,k,n,1) - q(j,k,n,1)
          f2i1   = q(j,k,n,1)   - q(j-1,k,n,1)

          a1     = 3.0*(f2i*f2i1+epsj)
          a2     = 2.0*(f2i-f2i1)**2 + a1
          f3i    = a1/a2
          f3i = 1.0

          q(jsf,kfine,n,2) = q(jsf,kfine,n,2) - &
            0.125*f3i*(q(j+1,k,n,1) - q(j-1,k,n,1))
          q(jef,kfine,n,2) = q(jef,kfine,n,2) + &
            0.125*f3i*(q(j+1,k,n,1) - q(j-1,k,n,1))
        enddo

        do jfine = jsf,jef
!          if (ksf.eq.kbeg(2).or.kef.eq.kend(2)) cycle
          f2i    = q(j,k+1,n,1) - q(j,k,n,1)
          f2i1   = q(j,k,n,1)   - q(j,k-1,n,1)
          a1     = 3.0*(f2i*f2i1+epsk)
          a2     = 2.0*(f2i-f2i1)**2 + a1
          f3i    = a1/a2
          f3i = 1.0

          q(jfine,ksf,n,2) = q(jfine,ksf,n,2) - &
            0.125*f3i*(q(j,k+1,n,1) - q(j,k-1,n,1))
          q(jfine,kef,n,2) = q(jfine,kef,n,2) + &
            0.125*f3i*(q(j,k+1,n,1) - q(j,k-1,n,1))
        enddo
      enddo
    enddo
  enddo

  end subroutine spatial_interpolator_muscl

!***********************************************************************
  subroutine fourier_interpolator(a,b,n,m)
!**********************************************************************
   integer :: n,m
   real :: a(n),b(m)

   integer :: i,nby2,mby2,m1
   complex, allocatable :: fhat(:),ghat(:)

   allocate(fhat(n),ghat(m))

!..forward fourier transform or DFT
   call ftf(a,fhat,n)

   do i=1,m
      ghat(i)=0
   enddo

   if (m.ge.n) then

!..fourier extension obtained by zero padding

     if (mod(n,2).eq.0) then
        nby2=n/2
        ghat(1:nby2)=fhat(1:nby2)
        ghat(nby2+1)=0.5*(fhat(nby2+1) + conjg(fhat(nby2+1)))
        ghat(m-nby2+2:m) = conjg(ghat(nby2:2:-1))
     else
        nby2=(n-1)/2
        ghat(1:nby2+1)=fhat(1:nby2+1)
        ghat(m-nby2+1:m) = conjg(ghat(nby2+1:2:-1))
     endif

   else

!.. fourier truncation (if m < n)

       if (mod(m,2).eq.0) then
          mby2=m/2
          ghat(1:mby2)=fhat(1:mby2)
          ghat(mby2+1)=fhat(mby2+1)+fhat(n-mby2+1)
          ghat(mby2+2:m) = conjg(ghat(mby2:2:-1))
       else
          mby2=(m-1)/2
          ghat(1:mby2+1)=fhat(1:mby2+1)
          ghat(mby2+2:m) = conjg(ghat(mby2+1:2:-1))
       endif
         
    endif

    ghat = ghat*m/n

!.. backward fourier transform or IDFT
    call ftb(ghat,b,m)

  end subroutine fourier_interpolator

!**********************************************************************
  subroutine ftf(u,uhat,N)
!..forward fourier transform
!*********************************************************************
    integer :: N
    real :: u(n)
    complex :: uhat(n)

    !local variables
    integer :: j,k
    real :: theta,pi

    pi = 4.*atan(1.)

    uhat = 0 
    do j = 0,N-1
      do k = 0,N-1
        theta = 2*pi*k*j/N
        uhat(j+1) = uhat(j+1) + u(k+1)*cmplx(cos(theta),-sin(theta))
      enddo
    enddo
   
 end subroutine ftf

!*********************************************************************
  subroutine ftb(uhat,u,N)
!..backward fourier transform
!*********************************************************************
    integer :: N
    real :: u(n)
    complex :: uhat(n)

    !local variables
    integer :: j,k
    real :: theta,pi

    pi = 4.*atan(1.)

    u = 0. 
    do j = 0,N-1
      do k = 0,N-1
        theta = 2*pi*k*j/N
        u(j+1) = u(j+1) + uhat(k+1)*cmplx(cos(theta),sin(theta))
      enddo
      u(j+1) = u(j+1)/N
    enddo
   
  end subroutine ftb

!*********************************************************************
  subroutine spline_interpolator(x,y,xk,yk,n,m)

!**********************************************************************
    integer :: n,m
    real :: x(n),y(n),xk(m),yk(m)
    real :: coef(2*n)

!..local variables

    real :: sigma,xout,yout
    integer :: i

!..extrapolate if outside the limit

    sigma=0.5
    call splico(x,y,coef,sigma,n)
      
    do i=1,m
         
      if (xk(i).lt.x(1)) then
         yk(i)=(xk(i)-x(2))*(y(2)-y(1))/(x(2)-x(1))+y(2)
      elseif (xk(i).gt.x(n)) then
         yk(i)=(xk(i)-x(n))*(y(n)-y(n-1))/(x(n)-x(n-1))+y(n)
      else
        xout=xk(i)
        call speval(x,y,coef,sigma,xout,yout,n)
        yk(i)=yout
      endif
       
    enddo
         
  end subroutine spline_interpolator
      
!**********************************************************************
      
  subroutine splico(x,y,coef,sigma,n)                               

!...this subroutine computes the spline coefficients necessary to fit an
!   exponential spline to the input abscissa x and ordinate y 
!
! argument list:
!     x(n)       input   input data abscissa
!     y(n)       input   input data ordinate
!     coef(2*n)  output  spline coeficients
!     sigma      input   spline tension factor
!     n          input   input array size
!
!**********************************************************************
    integer :: n
    real :: x(n),y(n),coef(n*2)
    real :: sigma
!                                                                            
!.. local variables
    integer :: i
    real :: delx1,delx2,delx12,dels,dx1,dx2
    real :: c1,c2,c3,slpp1,deln,delnm1,sigmap 
    real :: exps,sinhs,sinhin,diag1,diag2,diagin 
    real :: spdiag,slppn,delnn
   
!--> end points                                                                
!                                                                              
    delx1= x(2)-x(1)                                                  
    dx1= (y(2)-y(1))/delx1                                            
    delx2= x(3)-x(2)                                                  
    delx12= x(3)-x(1)                                                 
    c1= -(delx12+delx1)/delx12/delx1                                  
    c2= delx12/delx1/delx2                                            
    c3= -delx1/delx12/delx2                                           
    slpp1= c1*y(1)+c2*y(2)+c3*y(3)                                    
    deln= x(n)-x(n-1)                                                 
    delnm1= x(n-1)-x(n-2)                                             
    delnn= x(n)-x(n-2)                                                
    c1= (delnn+deln)/delnn/deln                                       
    c2= -delnn/deln/delnm1                                            
    c3= deln/delnn/delnm1                                             
    slppn= c3*y(n-2)+c2*y(n-1)+c1*y(n)                                
    sigmap= sigma*(n-1)/(x(n)-x(1))                                   
    dels= sigmap*delx1                                                
    exps= exp(dels)                                                   
    sinhs= (exps-1./exps)/2.                                          
    sinhin= 1./(delx1*sinhs)                                          
    diag1= sinhin*(dels*0.5*(exps+1./exps)-sinhs)                     
    diagin= 1./diag1                                                  
    coef(1)= diagin*(dx1-slpp1)                                       
    spdiag= sinhin*(sinhs-dels)                                       
    coef(n+1)= diagin*spdiag                                          
!                                                                       
!--> interior points                                                    
!                                                                       
    do i=2,n-1                                                        
      delx2= x(i+1)-x(i)                                              
      dx2= (y(i+1)-y(i))/delx2                                        
      dels= sigmap*delx2                                              
      exps= exp(dels)                                                 
      sinhs= (exps-1./exps)/2.                                        
      sinhin= 1./(delx2*sinhs)                                        
      diag2= sinhin*(dels*0.5*(exps+1./exps)-sinhs)                   
      diagin= 1./(diag1+diag2-spdiag*coef(n+i-1))                     
      coef(i)= diagin*(dx2-dx1-spdiag*coef(i-1))                      
      spdiag= sinhin*(sinhs-dels)                                     
      coef(i+n)= diagin*spdiag                                        
      dx1= dx2                                                        
      diag1= diag2                                                    
    enddo                                                             
    diagin= 1./(diag1-spdiag*coef(2*n-1))                             
    coef(n)= diagin*(slppn-dx2-spdiag*coef(n-1))                      
    do i=n-1,1,-1                                                     
      coef(i)= coef(i)-coef(i+n)*coef(i+1)                            
    enddo                                                             
                                                                        
  end subroutine splico

!**********************************************************************
                                                                        
  subroutine speval(x,y,coef,sigma,xout,yout,n)                     

!...this subroutine evaluates an exponential spline at the abscissa xout
!   and returns the value yout.
!
! argument list:
!     x(n)       input   input data abscissa
!     y(n)       input   input data ordinate
!     coef(2*n)  output  spline coeficients
!     sigma      input   spline tension factor
!     xout       input   abscissa at which spline evaluated
!     yout       input   spline value
!     n          input   input array size
!
!**********************************************************************

    integer :: n
    real :: x(n),y(n),coef(n*2)
    real :: sigma, xout, yout
!
!..   local variables
    integer :: idx,indm1,i
    real :: s,sigmap,del1,del2,dels,exps1,exps
    real :: sinhd1,sinhd2,sinhs

    idx = 0
    indm1 = 0

!**spline evaluation                                                           
    s= x(n)-x(1)                                                      
    sigmap= sigma*(n-1)/s                                             
    do i=1,n-1                                                        
      if(xout.ge.x(i) .and. xout.le.x(i+1)) then                      
        idx = i+1                                                    
        indm1= i                                                      
        goto 100                                                      
      endif                                                           
    enddo                                                             
 100 del1= xout-x(indm1)                                               
    del2= x(idx)-xout                                               
    dels= x(idx)-x(indm1)                                           
    exps1= exp(sigmap*del1)                                           
    sinhd1= 0.5*(exps1-1./exps1)                                      
    exps= exp(sigmap*del2)                                            
    sinhd2= 0.5*(exps-1./exps)                                        
    exps= exps1*exps                                                  
    sinhs= 0.5*(exps-1./exps)                                         
    yout= (coef(idx)*sinhd1+coef(indm1)*sinhd2)/sinhs +  &
          ((y(idx)-coef(idx))*del1+(y(indm1)-coef(indm1))*del2) &
          /dels      
!                                                                              
   end subroutine speval

!**********************************************************************

end module interpolators
