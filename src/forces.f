c**********************************************************************
      subroutine compute_forces(jd,kd,x,y,xv,yv,q,xx,xy,yx,yy,
     &  cfx_tot,cfy_tot,cm_tot,cl_tot,cd_tot,cpower_tot,im,fprefix)
c
c  this subroutine computes 2-d sectional forces, moment and power
c
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer jd,kd,im
      real q(jd,kd,nspec,nq),x(jd,kd,nspec),y(jd,kd,nspec)
      real xv(jmax,kmax,nspec),yv(jmax,kmax,nspec)
      real xx(jd,kd,nspec),xy(jd,kd,nspec)
      real yx(jd,kd,nspec),yy(jd,kd,nspec)
      real cfx_tot,cfy_tot,cm_tot,cl_tot,cd_tot,cpower_tot
      character*40, fprefix

      ! local variables

      integer nsp
      real x0,y0,cfx,cfy,cm0,cl,cd,cm,cpower
      real ca,sa
      character*40, forcefile, integer_string
      integer :: logf

c***  first executable statement
      pi  = 4.0 * atan( 1.0 )
      ca  = cos( pi * alfa / 180.0 )
      sa  = sin( pi * alfa / 180.0 )

      x0 = 0.25
      y0 = 0.00

      cfx_tot    = 0. 
      cfy_tot    = 0.
      cl_tot     = 0.
      cd_tot     = 0.
      cm_tot     = 0.
      cpower_tot = 0.

      write(integer_string,*) 10+im
      forcefile = trim(adjustl(fprefix))//'.'//trim(adjustl(integer_string))
      logf = 11
      open(unit=logf,file=forcefile,access='append',form='formatted')

!$OMP PARALLEL REDUCTION(+:cfx_tot,cfy_tot,cl_tot,cd_tot,cm_tot,cpower_tot) 
!$OMP& IF(NSPEC > 1)
!$OMP& PRIVATE(cfx,cfy,cm0,cl,cd,cm,cpower)
!$OMP DO ORDERED
      spectralloop: do nsp = 1,nspec
        call force2d(jd,kd,x(:,:,nsp),y(:,:,nsp),
     <       xv(:,:,nsp),yv(:,:,nsp),q(:,:,nsp,:),
     <       xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     <       cfx,cfy,cm0,cpower,nsp)

        cl = - cfx * sa + cfy * ca
        cd =   cfx * ca + cfy * sa
        cm =   cm0 - y0*cfx + x0*cfy

        cfx_tot    = cfx_tot    + cfx
        cfy_tot    = cfy_tot    + cfy
        !am cl_tot     = cl_tot     + cl
        !am cd_tot     = cd_tot     + cd
        !am cm_tot     = cm_tot     + cm
        !am cpower_tot = cpower_tot + cpower
        if(if_obj .and. if_objtot) then
          clobj(nsp) = cl
          cdobj(nsp) = cd
          cmobj(nsp) = cm

          cpower_tot = cpower_tot + cpower
          if(obj_ftype.eq.obj_ftype_cltot) then
            cl_tot     = cl_tot     + (cl-obj_f(nsp))**2
            cd_tot     = cd_tot     + cd**2
            cm_tot     = cm_tot     + cm**2
          elseif(obj_ftype.eq.obj_ftype_cdtot) then
            cl_tot     = cl_tot     + cl**2
            cd_tot     = cd_tot     + (cd-obj_f(nsp))**2
            cm_tot     = cm_tot     + cm**2
          elseif(obj_ftype.eq.obj_ftype_cmtot) then
            cl_tot     = cl_tot     + cl**2
            cd_tot     = cd_tot     + cd**2
            cm_tot     = cm_tot     + (cm-obj_f(nsp))**2
          else
            print '(A,I7,x,A)','OBJ TYPE: ',obj_ftype,' not supported '
            stop 
          end if
        else
          cl_tot     = cl_tot     + cl
          cd_tot     = cd_tot     + cd
          cm_tot     = cm_tot     + cm
          cpower_tot = cpower_tot + cpower
        endif

!$OMP ORDERED
        if (.not.timespectral) then
          if (iunst.eq.0) then
            write(logf,610) istep0,cl,cd,cm,cpower
          else
            write(logf,620) istep0,cl,cd,cm,cpower,theta_col
          endif
        else
          write(logf,630) istep0,nsp,cl,cd,cm,cpower
        endif
!$OMP END ORDERED

      enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL

!      if (timespectral) then
!        write(integer_string,*) 10+im
!        forcefile = 'fort_TS.'//trim(adjustl(integer_string))
!        open(unit=1000,file=forcefile,status='unknown',
!     &                                form='formatted')
!        write(1000,610) istep0,cfx_tot,cfy_tot,cm_tot,cpower_tot
!      endif

      close(logf)

610   FORMAT (i6,4(x,E18.10))
620   FORMAT (i6,5(x,E18.10))
630   FORMAT (2i6,4(x,E18.10))

      end subroutine compute_forces

c***********************************************************************
      subroutine force2d(jd,kd,x,y,xv,yv,q,xx,xy,yx,yy,
     <                   cfx,cfy,cm0,cpower,nsp)
c
c  this subroutine computes 2-d sectional cl, cd values.
c  modified version of subroutine clcd of arc2d.  vr. 01-15-91
c
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer jd,kd,nsp
      real q(jd,kd,nq),x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax)
      real xx(jd,kd),xy(jd,kd)
      real yx(jd,kd),yy(jd,kd)
      real cfx,cfy,cm0,cpower

      ! local variables
      
      real :: zz(jd),sxx(jd),syy(jd)
      real fmm,rr,rho,u,v,um,vm,e,vsq,pp,pp1,cp,cn,cc,cmle,cpinv,uxi,vxi
      real ueta,veta,xix,xiy,etaxxx,etax,etay
      real tauw,sfdiv,cnv,ccv,cmlev,cpv,uinf2,amu
      real chord,alngth

      integer k,k2,k3,j1,j2,j,jm1,jp1,jv,kv
      integer js,je
      !asitav
      character*60 int_str,filecp
      integer logcp,nobjf

c***  first executable statement

c..initialize constants and variables
c
c..set limits of integration
c     
      k  = kbeg
      k2 = k+1
      k3 = k2+1
      j1 = jtail1
      j2 = jtail2

      if (fsmach.ne.0) then
        fmm = fsmach
      elseif (fmtip.ne.0) then
        fmm = fmtip
      else
        fmm = 1.0
      endif
c
c..compute cp at grid points and store in an array 
c
      do 10 j=j1,j2
        rr=1./q(j,k,1)
        rho=q(j,k,1)*q(j,k,nq)
        u=q(j,k,2)*rr
        v=q(j,k,3)*rr
        e=q(j,k,4)*q(j,k,nq)
        vsq=u*u+v*v
        pp=gm1*(e-rho*vsq/2.)
        pp1=pp/pinf
        cp=2.*(pp1-1.)/(gamma*fmm**2)
        zz(j) = cp

        rr=1./q(j,k-1,1)
        rho=q(j,k-1,1)*q(j,k-1,nq)
        u=q(j,k-1,2)*rr
        v=q(j,k-1,3)*rr
        e=q(j,k-1,4)*q(j,k-1,nq)
        vsq=u*u+v*v
        pp=gm1*(e-rho*vsq/2.)
        pp1=pp/pinf
        cp=2.*(pp1-1.)/(gamma*fmm**2)
        zz(j) = 0.5*(zz(j) + cp)
10    continue

      !am !test overwrite with expt values (asitav)
      !am !----------------------------------------
      !am do j=j1,j2
      !am   nobjf = jd*(nsp-1)+j !jd per each TS
      !am   zz(j) = obj_f(nobjf)
      !am end do
      !am !----------------------------------------

      !if (iunst.eq.0.and..not.timespectral) then
      !if (mod(istep0,nrest).eq.0.or.istep0.eq.2) then
        write(int_str,"(I7)")nsp
        filecp='cp_TS'//trim(adjustl(int_str))//'.dat'
        !am open(unit=10,file='cp.dat',status='unknown')
        logcp = 90000+nsp
        open(unit=logcp,file=filecp,status='unknown')
        if(if_obj) then
          do j=j1,j2
            write(logcp,*) obj_cpx(j),zz(j),obj_cpy(j)
          enddo
        else
          do j=j1,j2
            write(logcp,*) x(j,k),zz(j),y(j,k)
          enddo
        end if
        close(logcp)
      !endif

c..compute normal force coefficient and chord directed force coeff
c..chord taken as one in all cases, modified later

      cn = 0.
      cc = 0.
      cmle = 0.
      cpinv = 0.
      do 11 j=j1,j2
        jv = j - nhalo
        kv = k - nhalo
        cn = cn - zz(j)*( xv(jv+1,kv) - xv(jv,kv))
        cc = cc + zz(j)*( yv(jv+1,kv) - yv(jv,kv))
        cmle = cmle + zz(j)*
     &     (   0.5*(xv(jv+1,kv)+xv(jv,kv))*(xv(jv+1,kv)-xv(jv,kv)) +
     &         0.5*(yv(jv+1,kv)+yv(jv,kv))*(yv(jv+1,kv)-yv(jv,kv)) )
        u=q(j,k,2)/q(j,k,1)
        v=q(j,k,3)/q(j,k,1)
        um=q(j,k-1,2)/q(j,k-1,1)
        vm=q(j,k-1,3)/q(j,k-1,1)
        cpinv = cpinv + 0.5*(u + um)*zz(j)*( yv(jv+1,kv) - yv(jv,kv))
     &                - 0.5*(v + vm)*zz(j)*( xv(jv+1,kv) - xv(jv,kv))
   11 continue

      cn = cn - half*cn
      cmle = cmle - half*cmle
      cc = cc + half*cc

c..viscous coefficent of friction calculation
c..re already has fsmach scaling

      cnv = 0.
      ccv = 0.
      cmlev = 0.
      cpv = 0.

      if(.not.invisc) then
        alngth= 1.
        amu   = rinf*alngth/rey
        uinf2 = fmm**2
c
        j = j1
        jp1 = j+1
        uxi = q(jp1,k,2)/q(jp1,k,1)-q(j,k,2)/q(j,k,1)
        vxi = q(jp1,k,3)/q(jp1,k,1)-q(j,k,3)/q(j,k,1)
        ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *       -.5*q(j,k3,2)/q(j,k3,1)
        veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *       -.5*q(j,k3,3)/q(j,k3,1)
        xix = xx(j,k)
        xiy = xy(j,k)
        etax = yx(j,k)
        etay = yy(j,k)
        sfdiv = 2./3*(uxi*xix + ueta*etax + vxi*xiy + veta*etay)
        sxx(j) = amu*(2.*(uxi*xix + ueta*etax) - sfdiv) 
        syy(j) = amu*(2.*(vxi*xiy + veta*etay) - sfdiv) 
        tauw= amu*((uxi*xiy+ueta*etay)+(vxi*xix+veta*etax))

        sxx(j) = sxx(j)/(.5*rinf*uinf2)
        syy(j) = syy(j)/(.5*rinf*uinf2)
        zz(j)= tauw/(.5*rinf*uinf2)

        j = j2
        jm1 = j-1
        uxi = q(j,k,2)/q(j,k,1)-q(jm1,k,2)/q(jm1,k,1)
        vxi = q(j,k,3)/q(j,k,1)-q(jm1,k,3)/q(jm1,k,1)
        ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *       -.5*q(j,k3,2)/q(j,k3,1)
        veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *       -.5*q(j,k3,3)/q(j,k3,1)
        xix = xx(j,k)
        xiy = xy(j,k)
        etax = yx(j,k)
        etay = yy(j,k)
        sfdiv = 2./3*(uxi*xix + ueta*etax + vxi*xiy + veta*etay)
        sxx(j) = amu*(2.*(uxi*xix + ueta*etax) - sfdiv) 
        syy(j) = amu*(2.*(vxi*xiy + veta*etay) - sfdiv) 
        tauw= amu*((uxi*xiy+ueta*etay)+(vxi*xix+veta*etax))

        sxx(j) = sxx(j)/(.5*rinf*uinf2)
        syy(j) = syy(j)/(.5*rinf*uinf2)
        zz(j)= tauw/(.5*rinf*uinf2)

c..set limits

        js = j1+1
        je = j2-1
        do 110 j = js,je
          jp1 = j+1
          jm1 = j-1
          uxi=.5*(q(jp1,k,2)/q(jp1,k,1)-q(jm1,k,2)/q(jm1,k,1))
          vxi=.5*(q(jp1,k,3)/q(jp1,k,1)-q(jm1,k,3)/q(jm1,k,1))
          ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *               -.5*q(j,k3,2)/q(j,k3,1)
          veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *               -.5*q(j,k3,3)/q(j,k3,1)
          xix = xx(j,k)
          xiy = xy(j,k)
          etax = yx(j,k)
          etay = yy(j,k)
          sfdiv = 2./3*(uxi*xix + ueta*etax + vxi*xiy + veta*etay)
          sxx(j) = amu*(2.*(uxi*xix + ueta*etax) - sfdiv) 
          syy(j) = amu*(2.*(vxi*xiy + veta*etay) - sfdiv) 
          tauw= amu*((uxi*xiy+ueta*etay)+(vxi*xix+veta*etax))

          sxx(j) = sxx(j)/(.5*rinf*uinf2)
          syy(j) = syy(j)/(.5*rinf*uinf2)
          zz(j)= tauw/(.5*rinf*uinf2)
110     continue

      if (iunst.eq.0.and..not.timespectral) then
        open(unit=10,file='cf.dat',status='unknown')
        do j=j1,j2
          write(10,*) x(j,k),zz(j)
        enddo
        close(10)
      endif
c
c..compute viscous normal and axial forces
c
        do 111 j=j1,j2
          jv = j - nhalo
          kv = k - nhalo
          ccv = ccv + zz(j) *( xv(jv+1,kv) - xv(jv,kv)) - 
     &                sxx(j)*( yv(jv+1,kv) - yv(jv,kv))
          cnv = cnv - zz(j) *( yv(jv+1,kv) - yv(jv,kv)) + 
     &                syy(j)*( xv(jv+1,kv) - xv(jv,kv))
          cmlev = cmlev + zz(j)*
     &    ( 0.5*(xv(jv+1,kv)+xv(jv,kv))*(yv(jv+1,kv) -yv(jv,kv))       +
     &      0.5*(yv(jv+1,kv)+yv(jv,kv))*(xv(jv+1,kv) -xv(jv,kv)) )     -
     &      0.5*(xv(jv+1,kv)+xv(jv,kv))*(xv(jv+1,kv) -xv(jv,kv))*syy(j)-
     &      0.5*(yv(jv+1,kv)+yv(jv,kv))*(yv(jv+1,kv) -yv(jv,kv))*sxx(j)
          u=q(j,k,2)/q(j,k,1)
          v=q(j,k,3)/q(j,k,1)
          um=q(j,k-1,2)/q(j,k-1,1)
          vm=q(j,k-1,3)/q(j,k-1,1)
          cpv = cpv + 0.5*(u + um)*zz(j) *( xv(jv+1,kv) - xv(jv,kv))
     &              - 0.5*(u + um)*sxx(j)*( yv(jv+1,kv) - yv(jv,kv))
     &              - 0.5*(v + vm)*zz(j) *( yv(jv+1,kv) - yv(jv,kv))
     &              + 0.5*(v + vm)*syy(j)*( xv(jv+1,kv) - xv(jv,kv))

111     continue

        cnv = cnv - half*cnv
        cmlev = cmlev - half*cmlev
        ccv = ccv + half*ccv
      endif

cc..modified for chord being non-unity
c
c      chord = sqrt((x(jtail1,k)-x(jle,k))**2
c     &       +     (y(jtail1,k)-y(jle,k))**2)
c      cl2d=cl2d/chord
c      clvv2d=clvv2d/chord
c      cd2d=cd2d/chord
c      cdvv2d=cdvv2d/chord
c      cm2d=cm2d/chord
c      cmvv2d=cmvv2d/chord

      cfx = cc + ccv
      cfy = cn + cnv
      cm0 =  cmle + cmlev
      cpower = (cpinv + cpv)/fmm

      return
      end

c*************************************************************************
