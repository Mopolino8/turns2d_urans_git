      subroutine calc_norm(a, b, norm_type, tmp_norm)
      real a, b, tmp_norm
      integer :: norm_type

      if(norm_type .eq. 2) then
         tmp_norm = (a - b)**2
      else
         print *, "The specified norm type is not defined"
         stop
      end if
      return
      end subroutine calc_norm


      subroutine calc_cf(jd, kd, x, y, xv, yv, q, xx, xy, yx, yy, cf)
      use params_global
      implicit none
      
      integer :: jd,kd
      real :: q(jd,kd,nq),x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax)
      real :: xx(jd,kd),xy(jd,kd)
      real :: yx(jd,kd),yy(jd,kd)
      real :: cf(jd)

      real :: zz1(jd)
      real :: fmm,rr,rho,u,v,um,vm,e,vsq,pp,pp1,cp,cn,cc,cmle,cpinv,uxi,vxi
      real :: ueta,veta,xix,xiy,etax,etay
      real :: tauw,sfdiv,cnv,ccv,cmlev,cpv,uinf2,amu
      real :: chord,alngth
      integer :: j, k, k2, k3, j1, j2, jm1, jp1, js, je

    
      k  = kbeg
      k2 = k+1
      k3 = k2+1
      j1 = jtail1
      j2 = jtail2

!     set normalization factor
      if (fsmach.ne.0) then
         fmm = fsmach
      elseif (fmtip.ne.0) then
         fmm = fmtip
      else
         fmm = 1.0
      endif

      
      alngth= 1.
      amu   = rinf*alngth/rey
      uinf2 = fmm**2
c     
      j = j1
      jp1 = j+1
      uxi = q(jp1,k,2)/q(jp1,k,1)-q(j,k,2)/q(j,k,1)
      vxi = q(jp1,k,3)/q(jp1,k,1)-q(j,k,3)/q(j,k,1)
      ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *     -.5*q(j,k3,2)/q(j,k3,1)
      veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *     -.5*q(j,k3,3)/q(j,k3,1)
      xix = xx(j,k)
      xiy = xy(j,k)
      etax = yx(j,k)
      etay = yy(j,k)
      
      zz1(j) = amu*((uxi*xiy+ueta*etay)-(vxi*xix+veta*etax))/(.5*rinf*uinf2)
      
      j = j2
      jm1 = j-1
      uxi = q(j,k,2)/q(j,k,1)-q(jm1,k,2)/q(jm1,k,1)
      vxi = q(j,k,3)/q(j,k,1)-q(jm1,k,3)/q(jm1,k,1)
      ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *     -.5*q(j,k3,2)/q(j,k3,1)
      veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *     -.5*q(j,k3,3)/q(j,k3,1)
      
      xix = xx(j,k)
      xiy = xy(j,k)
      etax = yx(j,k)
      etay = yy(j,k)
      
    
      zz1(j) = amu*((uxi*xiy+ueta*etay)-(vxi*xix+veta*etax))/(.5*rinf*uinf2)
      
      
c..   set limits
      
      js = j1+1
      je = j2-1
      do 110 j = js,je
         jp1 = j+1
         jm1 = j-1
         uxi=.5*(q(jp1,k,2)/q(jp1,k,1)-q(jm1,k,2)/q(jm1,k,1))
         vxi=.5*(q(jp1,k,3)/q(jp1,k,1)-q(jm1,k,3)/q(jm1,k,1))
         ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *        -.5*q(j,k3,2)/q(j,k3,1)
         veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *        -.5*q(j,k3,3)/q(j,k3,1)
         xix = xx(j,k)
         xiy = xy(j,k)
         
         etax = yx(j,k)
         etay = yy(j,k)

         zz1(j) = amu*((uxi*xiy+ueta*etay)-(vxi*xix+veta*etax))/(.5*rinf*uinf2)
 110  continue
      
      cf(:) = zz1(:)

      end subroutine calc_cf


      subroutine calc_cp(jd, kd, q, cp)
      use params_global
      implicit none
      integer :: jd, kd
      real :: q(jd, kd, nq), cp(jd)
      
      integer :: j
      integer :: j1, j2, k
      real :: fmm, rr, rho, u, v, e, vsq, pp, pp1, cp_
      
      j1 = jtail1
      j2 = jtail2
      k = kbeg
      
!     set normalization factor
      if (fsmach.ne.0) then
         fmm = fsmach
      elseif (fmtip.ne.0) then
         fmm = fmtip
      else
         fmm = 1.0
      endif
      
!     compute cp at grid points and store in an array 
      do j=jtail1,jtail2
         rr=1./q(j,k,1)
         rho=q(j,k,1)*q(j,k,nq)
         u=q(j,k,2)*rr
         v=q(j,k,3)*rr
         e=q(j,k,4)*q(j,k,nq)
         vsq=u*u+v*v
         pp=gm1*(e-rho*vsq/2.)
         pp1=pp/pinf
         cp_=2.*(pp1-1.)/(gamma*fmm**2)
         cp(j) = cp_

         rr=1./q(j,k-1,1)
         rho=q(j,k-1,1)*q(j,k-1,nq)
         u=q(j,k-1,2)*rr
         v=q(j,k-1,3)*rr
         e=q(j,k-1,4)*q(j,k-1,nq)
         vsq=u*u+v*v
         pp=gm1*(e-rho*vsq/2.)
         pp1=pp/pinf
         cp_=2.*(pp1-1.)/(gamma*fmm**2)
         cp(j) = 0.5*(cp(j) + cp_)
      end do
      end subroutine calc_cp
      
      subroutine objectivef(jd,kd,x,y,xv,yv,q,xx,xy,yx,yy,opt_obj,im)
      use params_global
      implicit none
      
      integer :: jd,kd,im
      real :: q(jd,kd,nq), x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real :: xx(jd,kd), xy(jd,kd)
      real :: yx(jd,kd), yy(jd,kd)
      real :: opt_obj
      
      real :: coeff(jd), tmp_norm
      integer :: j, k, np

      
      if (obj_ftype .eq. obj_ftype_cp) then
         call calc_cp(jd, kd, q, coeff)
      else
         print *, "The objective type is not defined"
         stop
      end if

      opt_obj = 0.0
      if (obj_pointwise) then
         opt_obj = coeff(obj_points(1)) -  obj_f(obj_points(1))
      else
         do np = 1, obj_npoints
            call calc_norm(coeff(obj_points(np)), obj_f(obj_points(np)), obj_normtype, tmp_norm)
            opt_obj = opt_obj + tmp_norm
         end do

         if(obj_if_reg) then
            do k=1,kd
               do j=1,jd
                  opt_obj = opt_obj + obj_reg_fac*(obj_coeff_prod(j,k,im) - 1.0)**2
               end do
            end do
         end if

      end if
      print *, "Objective function value: ", opt_obj
      write(747, *) opt_obj
      end subroutine      
