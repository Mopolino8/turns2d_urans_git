c***********************************************************************
      subroutine bc(q,x,y,xx,xy,yx,yy,ug,vg,im,jd,kd,bt)
c
c  explicitly update the mesh boundaries
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer im,jd,kd,j,k,js,je,ks,ke,idir,ib
      real q(jd,kd,nq),x(jd,kd),y(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd),bt(jd,kd)
      real ug(jd,kd),vg(jd,kd)

c***  first executable statement

C$AD II-LOOP
      do ib=1,nbc_all(im)
        js = jbcs_all(ib,im)
        je = jbce_all(ib,im)
        ks = kbcs_all(ib,im)
        ke = kbce_all(ib,im)
        if(js.lt.0) js = jmax + js + 1
        if(ks.lt.0) ks = kmax + ks + 1
        if(je.lt.0) je = jmax + je + 1
        if(ke.lt.0) ke = kmax + ke + 1

        if (js.gt.1) js = js + nhalo
        if (ks.gt.1) ks = ks + nhalo
        if (je.lt.jmax) je = je + nhalo - 1
        if (ke.lt.kmax) ke = ke + nhalo - 1
        if (je.eq.jmax) je = jd
        if (ke.eq.kmax) ke = kd

        idir = ibdir_all(ib,im)

c.. inviscid wall bc (k = 1 only) 
        if (ibtyp_all(ib,im).eq.3) then
          call bcwall(q,xx,xy,yx,yy,ug,vg,
     &         jd,kd,js,je,ks,ke,idir,.true.)

c.. viscous wall bc
        elseif (ibtyp_all(ib,im).eq.4) then
          call bcwall(q,xx,xy,yx,yy,ug,vg,
     &         jd,kd,js,je,ks,ke,idir,.false.)

c.. viscous or inviscid based on input 
        elseif (ibtyp_all(ib,im).eq.5) then
          call bcwall(q,xx,xy,yx,yy,ug,vg,
     &         jd,kd,js,je,ks,ke,idir,invisc)

c.. extrapolation bc
        elseif (ibtyp_all(ib,im).eq.10) then
          call bcextp(q,jd,kd,js,je,ks,ke,idir)

c.. symmetric bc
        elseif (ibtyp_all(ib,im).eq.11) then
          call bcsym(q,jd,kd,js,je,ks,ke,idir)

c.. periodic bc
        elseif (ibtyp_all(ib,im).eq.22) then
          call bcperiodic(q,jd,kd,js,je,ks,ke,idir)

c.. averaging bc for wake
        elseif (ibtyp_all(ib,im).eq.51) then
          call bcwake(q,x,y,jd,kd,js,je,ks,ke,idir)

c.. enforcing freestream bc
        elseif (ibtyp_all(ib,im).eq.46) then
          call bcout_inf(q,x,y,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)

c.. characteristic freesream bc
        elseif (ibtyp_all(ib,im).eq.47) then
!          if (.not.iprecon) then
!            call bcout(q,x,y,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
!          else
            call bcout_trkl(q,x,y,xx,xy,yx,yy,ug,vg,bt,
     &                       jd,kd,js,je,ks,ke,idir)
!          endif
        endif
      enddo
      return
      end

c***********************************************************************
      subroutine rbc(q,s,x,y,xx,xy,yx,yy,ug,vg,im,jd,kd,bt)
c
c  update residuals at mesh boundaries
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer im,jd,kd,j,k,js,je,ks,ke,idir,ib
      real s(jd,kd,nv),q(jd,kd,nq),x(jd,kd),y(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd),bt(jd,kd)
      real ug(jd,kd),vg(jd,kd)

c***  first executable statement

C$AD II-LOOP
      do ib=1,nbc_all(im)
        js = jbcs_all(ib,im)
        je = jbce_all(ib,im)
        ks = kbcs_all(ib,im)
        ke = kbce_all(ib,im)
        if(js.lt.0) js = jmax + js + 1
        if(ks.lt.0) ks = kmax + ks + 1
        if(je.lt.0) je = jmax + je + 1
        if(ke.lt.0) ke = kmax + ke + 1

        if (js.gt.1) js = js + nhalo
        if (ks.gt.1) ks = ks + nhalo
        if (je.lt.jmax) je = je + nhalo - 1
        if (ke.lt.kmax) ke = ke + nhalo - 1
        if (je.eq.jmax) je = jd
        if (ke.eq.kmax) ke = kd

        idir = ibdir_all(ib,im)

        if (ibtyp_all(ib,im).eq.3  .or.
     &      ibtyp_all(ib,im).eq.4  .or.
     &      ibtyp_all(ib,im).eq.5  .or.
     &      ibtyp_all(ib,im).eq.46 .or.
     &      ibtyp_all(ib,im).eq.47     ) then
          call zero_residuals(s,jd,kd,js,je,ks,ke,idir)

c.. inviscid wall bc (k = 1 only) 
        elseif (ibtyp_all(ib,im).eq.3) then
          call rbcwall(q,s,xx,xy,yx,yy,ug,vg,
     &         jd,kd,js,je,ks,ke,idir,.true.)

c.. viscous wall bc
        elseif (ibtyp_all(ib,im).eq.4) then
          call rbcwall(q,s,xx,xy,yx,yy,ug,vg,
     &         jd,kd,js,je,ks,ke,idir,.false.)

c.. viscous or inviscid based on input 
        elseif (ibtyp_all(ib,im).eq.5) then
          call rbcwall(q,s,xx,xy,yx,yy,ug,vg,
     &         jd,kd,js,je,ks,ke,idir,invisc)

c.. extrapolation bc
        elseif (ibtyp_all(ib,im).eq.10) then
          call rbcextp(q,s,jd,kd,js,je,ks,ke,idir)

c.. symmetric bc
        elseif (ibtyp_all(ib,im).eq.11) then
          call rbcsym(q,s,jd,kd,js,je,ks,ke,idir)

c.. periodic bc
        elseif (ibtyp_all(ib,im).eq.22) then
          call rbcperiodic(q,s,jd,kd,js,je,ks,ke,idir)

c.. averaging bc for wake
        elseif (ibtyp_all(ib,im).eq.51) then
          call rbcwake(q,s,jd,kd,js,je,ks,ke,idir)

c.. characteristic freesream bc
        elseif (ibtyp_all(ib,im).eq.47) then
!          if (.not.iprecon) then
!            call rbcout(q,s,x,y,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
!          else
            call rbcout_trkl(q,s,x,y,xx,xy,yx,yy,ug,vg,bt,
     &                       jd,kd,js,je,ks,ke,idir)
!          endif

        endif

      enddo

      return
      end

c***********************************************************************
      subroutine bcperiodic(q,jd,kd,js,je,ks,ke,idir)
c
c..periodic boundary for overlapping planes
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
      real scale
      
      iadd = sign(1,idir)

      pi = 4*atan(1.0)
      
      if(idir.eq.1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = j - js + 1
           jc = jd - 2*jj + jj1
 
           do k = ks,ke
             scale = q(jc,k,nq)/q(j,k,nq)
             q(j,k,1) = q(jc,k,1)*scale
             q(j,k,2) = q(jc,k,2)*scale
             q(j,k,3) = q(jc,k,3)*scale
             q(j,k,4) = q(jc,k,4)*scale
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j + 1
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             scale = q(jc,k,nq)/q(j,k,nq)
             q(j,k,1) = q(jc,k,1)*scale
             q(j,k,2) = q(jc,k,2)*scale
             q(j,k,3) = q(jc,k,3)*scale
             q(j,k,4) = q(jc,k,4)*scale
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks + 1
           kc = kd - 2*kk + kk1
 
           do j = js,je
             scale = q(j,kc,nq)/q(j,k,nq)
             q(j,k,1) = q(j,kc,1)*scale
             q(j,k,2) = q(j,kc,2)*scale
             q(j,k,3) = q(j,kc,3)*scale
             q(j,k,4) = q(j,kc,4)*scale
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k + 1
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             scale = q(j,kc,nq)/q(j,k,nq)
             q(j,k,1) = q(j,kc,1)*scale
             q(j,k,2) = q(j,kc,2)*scale
             q(j,k,3) = q(j,kc,3)*scale
             q(j,k,4) = q(j,kc,4)*scale
           enddo
        enddo

      endif

      return
      end

c***********************************************************************
      subroutine rbcperiodic(q,s,jd,kd,js,je,ks,ke,idir)
c
c..periodic boundary for overlapping planes
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real q(jd,kd,nq),s(jd,kd,nv)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
      real scale
      
      iadd = sign(1,idir)

      pi = 4*atan(1.0)
      
      if(idir.eq.1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = j - js + 1
           jc = jd - 2*jj + jj1
 
           do k = ks,ke
             scale = q(jc,k,nq)/q(j,k,nq)
             s(j,k,1) = s(jc,k,1)*scale
             s(j,k,2) = s(jc,k,2)*scale
             s(j,k,3) = s(jc,k,3)*scale
             s(j,k,4) = s(jc,k,4)*scale
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j + 1
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             scale = q(jc,k,nq)/q(j,k,nq)
             s(j,k,1) = s(jc,k,1)*scale
             s(j,k,2) = s(jc,k,2)*scale
             s(j,k,3) = s(jc,k,3)*scale
             s(j,k,4) = s(jc,k,4)*scale
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks + 1
           kc = kd - 2*kk + kk1
 
           do j = js,je
             scale = q(j,kc,nq)/q(j,k,nq)
             s(j,k,1) = s(j,kc,1)*scale
             s(j,k,2) = s(j,kc,2)*scale
             s(j,k,3) = s(j,kc,3)*scale
             s(j,k,4) = s(j,kc,4)*scale
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k + 1
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             scale = q(j,kc,nq)/q(j,k,nq)
             s(j,k,1) = s(j,kc,1)*scale
             s(j,k,2) = s(j,kc,2)*scale
             s(j,k,3) = s(j,kc,3)*scale
             s(j,k,4) = s(j,kc,4)*scale
           enddo
        enddo

      endif

      return
      end

c*********************************************************************
      subroutine bcout_inf(q,x,y,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
c
c     Outflow boundary condition 
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,js,je,ks,ke,idir
      real  q(jd,kd,nq)
      real  xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd),
     $      ug(jd,kd),vg(jd,kd),x(jd,kd),y(jd,kd)

      ! local variables
      real rj,rk
      integer j,k,iadd,iadir
      
c***  first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        do j  = js,je
          do k = ks,ke
            rj    = 1./q(j,k,nq)
            q(j,k,1) = rinf*rj
            q(j,k,2) = rinf*uinf*rj
            q(j,k,3) = rinf*vinf*rj
            q(j,k,4) = pinf/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                 q(j,k,3)**2)/q(j,k,1)
          enddo
        enddo
      elseif (iadir.eq.2) then
        do k  = ks,ke
          do j = js,je
            rk    = 1./q(j,k,nq)
            q(j,k,1) = rinf*rk
            q(j,k,2) = rinf*uinf*rk
            q(j,k,3) = rinf*vinf*rk
            q(j,k,4) = pinf/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                 q(j,k,3)**2)/q(j,k,1)
          enddo
        enddo
      endif

      return
      end

c***********************************************************************
      subroutine bcout(q,x,y,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
c
c     Outflow boundary condition 
c     Characteristic extrapolation based on Riemann invariants
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,js,je,ks,ke,idir
      real  q(jd,kd,nq)
      real  xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd),
     $      ug(jd,kd),vg(jd,kd),x(jd,kd),y(jd,kd)

      ! local variables
      real gm1i,gi,foso
      real xxn,xyn,yxn,yyn,snorm
      real rhoext,uext,vext,eext,pext
      real uinn,uexn,r1,r2,qn,cspe,c2,vel2,vel2ext,velcheck
      real stj,entro,u,v,press
      real rj,rjp1,rjp2,rk,rkp1,rkp2
      real rho1,rho2,u1,u2,v1,v2,p1,p2
      integer k,k1,k2,j,j1,j2,jc,kc,iadd,iadir
      
c***  first executable statement

      gm1i = 1./gm1
      gi   = 1./gamma

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then

        if (idir.eq.1) then
          j = je
        elseif (idir.eq.-1) then
          j = js
        endif
        j1 = j  + iadd
        j2 = j1 + iadd

        foso = 1.0
        do k = ks,ke
          rj    = 1./q(j,k,nq)
          rjp1  = 1./q(j1,k,1)
          rjp2  = 1./q(j2,k,1)
c..first-order extrapolation
          rhoext= (1.+foso)*q(j1,k,1)*q(j1,k,nq)-foso*q(j2,k,1)*q(j2,k,nq)
          uext  = (1.+foso)*q(j1,k,2)*rjp1-foso*q(j2,k,2)*rjp2
          vext  = (1.+foso)*q(j1,k,3)*rjp1-foso*q(j2,k,3)*rjp2
          eext  = (1.+foso)*q(j1,k,4)*q(j1,k,nq)-foso*q(j2,k,4)*q(j2,k,nq)
          pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2))
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(xx(j,k)**2+xy(j,k)**2)
          xxn = xx(j,k)*snorm
          xyn = xy(j,k)*snorm
c..calculate riemann invariants
          uinn = (uinf-ug(j,k))*xxn + (vinf-vg(j,k))*xyn
          r1 = uinn -2.*sqrt(gamma*pinf/rinf)*gm1i
          uexn = (uext-ug(j,k))*xxn + (vext-vg(j,k))*xyn
          r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
          qn = 0.5*(r1+r2)
          cspe = (r1-r2)*gm1*0.25
          c2 = cspe**2
c..is flow relatively subsonic or supersonic?
          vel2 = (uinf-ug(j,k))**2+(vinf-vg(j,k))**2
          vel2ext = (uext-ug(j,k))**2+(vext-vg(j,k))**2
          velcheck = abs(qn) !0.5*(vel2+vel2ext)
c..calculate contributions from interior and exterior
          if(qn.lt. 0) then
c..inflow boundary
            if(velcheck .lt. 1.0) then
c..fix three and extrapolate one (pressure)
              stj = qn - uinn
              entro = rinf**gamma/pinf
              u = uinf + stj*xxn
              v = vinf + stj*xyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rj
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..fix four
              q(j,k,1) = rinf*rj
              q(j,k,2) = rinf*uinf*rj
              q(j,k,3) = rinf*vinf*rj
              press = pinf
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          else
c..outflow boundary
            if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,and v
              stj = qn - uexn
              entro = rhoext**gamma/pext
              u = uext + stj*xxn
              v = vext + stj*xyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rj
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..extrapolate four
              q(j,k,1) = rhoext*rj
              q(j,k,2) = rhoext*uext*rj
              q(j,k,3) = rhoext*vext*rj
              press = pext
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          endif
        enddo

c..extrapolate everything to other halo cells

        foso = 0.0
        do jc = 1,je-js
          if (idir.eq.1) then
            j = je - jc
          elseif (idir.eq.-1) then
            j = js + jc
          endif
          j1 = j  + iadd
          j2 = j1 + iadd

          do k = ks,ke
            rho1 =  q(j1,k,1)*q(j1,k,nq)
            rho2 =  q(j2,k,1)*q(j2,k,nq)
            u1 = q(j1,k,2)/q(j1,k,1)  
            u2 = q(j2,k,2)/q(j2,k,1)  
            v1 = q(j1,k,3)/q(j1,k,1)  
            v2 = q(j2,k,3)/q(j2,k,1)  
            p1 = gm1*(q(j1,k,4) -0.5*(q(j1,k,2)**2+q(j1,k,3)**2)
     >                /q(j1,k,1))*q(j1,k,nq)
            p2 = gm1*(q(j2,k,4) -0.5*(q(j2,k,2)**2+q(j2,k,3)**2)
     >                /q(j2,k,1))*q(j2,k,nq)

            q(j,k,1) = ((1.+foso)*rho1 - foso*rho2 )/q(j,k,nq)
            q(j,k,2) = ((1.+foso)*u1 - foso*u2)*q(j,k,1)
            q(j,k,3) = ((1.+foso)*v1 - foso*v2)*q(j,k,1)
            press = (1.+foso)*p1 - foso*p2 
            q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     >                  +q(j,k,3)**2)/q(j,k,1)
          enddo
        enddo

      elseif (iadir.eq.2) then

        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif
        k1 = k  + iadd
        k2 = k1 + iadd

        foso = 1.0
        do j = js,je
          rk    = 1./q(j,k,nq)
          rkp1  = 1./q(j,k1,1)
          rkp2  = 1./q(j,k2,1)
c..first-order extrapolation
          rhoext= (1.+foso)*q(j,k1,1)*q(j,k1,nq)-foso*q(j,k2,1)*q(j,k2,nq)
          uext  = (1.+foso)*q(j,k1,2)*rkp1-foso*q(j,k2,2)*rkp2
          vext  = (1.+foso)*q(j,k1,3)*rkp1-foso*q(j,k2,3)*rkp2
          eext  = (1.+foso)*q(j,k1,4)*q(j,k1,nq)-foso*q(j,k2,4)*q(j,k2,nq)
          pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2))
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(yx(j,k)**2+yy(j,k)**2)
          yxn = yx(j,k)*snorm
          yyn = yy(j,k)*snorm
c..calculate riemann invariants
          uinn = (uinf-ug(j,k))*yxn + (vinf-vg(j,k))*yyn
          r1 = uinn -2.*sqrt(gamma*pinf/rinf)*gm1i
          uexn = (uext-ug(j,k))*yxn + (vext-vg(j,k))*yyn
          r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
          qn = 0.5*(r1+r2)
          cspe = (r2-r1)*gm1*0.25
          c2 = cspe**2
c..is flow relatively subsonic or supersonic?
          vel2 = (uinf-ug(j,k))**2+(vinf-vg(j,k))**2
          vel2ext = (uext-ug(j,k))**2+(vext-vg(j,k))**2
          velcheck = abs(qn) !0.5*(vel2+vel2ext)
c..calculate contributions from interior and exterior
          if(qn .lt. 0) then
c..inflow boundary
            if(velcheck .lt. 1.0) then
c..fix three and extrapolate one (pressure)
              stj = qn - uinn
              entro = rinf**gamma/pinf
              u = uinf + stj*yxn
              v = vinf + stj*yyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rk
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..fix four
              q(j,k,1) = rinf*rk
              q(j,k,2) = rinf*uinf*rk
              q(j,k,3) = rinf*vinf*rk
              press = pinf
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          else
c..outflow boundary
            if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,and v
              stj = qn - uexn
              entro = rhoext**gamma/pext
              u = uext + stj*yxn
              v = vext + stj*yyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rk
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..extrapolate four
              q(j,k,1) = rhoext*rk
              q(j,k,2) = rhoext*uext*rk
              q(j,k,3) = rhoext*vext*rk
              press = pext
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          endif
        enddo

c..extrapolate everything to other halo cells

        foso = 0.0
        do kc = 1,ke-ks
          if (idir.eq.2) then
            k = ke - kc
          elseif (idir.eq.-2) then
            k = ks + kc
          endif
          k1 = k  + iadd
          k2 = k1 + iadd

          do j = js,je
            rho1 =  q(j,k1,1)*q(j,k1,nq)
            rho2 =  q(j,k2,1)*q(j,k2,nq)
            u1 = q(j,k1,2)/q(j,k1,1)  
            u2 = q(j,k2,2)/q(j,k2,1)  
            v1 = q(j,k1,3)/q(j,k1,1)  
            v2 = q(j,k2,3)/q(j,k2,1)  
            p1 = gm1*(q(j,k1,4) -0.5*(q(j,k1,2)**2+q(j,k1,3)**2)
     >                /q(j,k1,1))*q(j,k1,nq)
            p2 = gm1*(q(j,k2,4) -0.5*(q(j,k2,2)**2+q(j,k2,3)**2)
     >                /q(j,k2,1))*q(j,k2,nq)

            q(j,k,1) = ((1.+foso)*rho1 - foso*rho2 )/q(j,k,nq)
            q(j,k,2) = ((1.+foso)*u1 - foso*u2)*q(j,k,1)
            q(j,k,3) = ((1.+foso)*v1 - foso*v2)*q(j,k,1)
            press = (1.+foso)*p1 - foso*p2 
            q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     >                  +q(j,k,3)**2)/q(j,k,1)
          enddo
        enddo

      endif

      return
      end

c***********************************************************************
      subroutine bcout_trkl(q,x,y,xx,xy,yx,yy,ug,vg,bt,jd,kd,
     <                      js,je,ks,ke,idir)
c
c     Outflow boundary condition with Turkel Preconditioning 
c     Characteristic extrapolation based on Riemann invariants
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,js,je,ks,ke,idir
      real  q(jd,kd,nq)
      real  xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd),
     $      ug(jd,kd),vg(jd,kd),x(jd,kd),y(jd,kd),bt(jd,kd)

      ! local variables
      real foso
      real rhoext,uext,vext,eext,pext
      real xxn,xyn,yxn,yyn,snorm
      real rho,rrho,uu,vv,e,uv2,cjkl,c2i,ge,qq
      real Z,R,S,Zi,bSq
      real rho1,rho2,u1,u2,v1,v2,p1,p2,press
      integer k,k1,k2,j,j1,j2,jc,kc,iadd,iadir
      real qext(4),qbar(4),delq(4),delw(4)
      
c***  first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then

        if (idir.eq.1) then
          j = je
        elseif (idir.eq.-1) then
          j = js
        endif
        j1 = j  + iadd
        j2 = j1 + iadd
c
        foso = 1.0
C$AD II-LOOP
        do k = ks,ke
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(xx(j,k)**2+xy(j,k)**2)
          xxn = xx(j,k)*snorm
          xyn = xy(j,k)*snorm
c..first-order extrapolation
          qext(1) = (1.+foso)*q(j1,k,1)*q(j1,k,nq)-foso*q(j2,k,1)*q(j2,k,nq)
          qext(2) = (1.+foso)*q(j1,k,2)*q(j1,k,nq)-foso*q(j2,k,2)*q(j2,k,nq)
          qext(3) = (1.+foso)*q(j1,k,3)*q(j1,k,nq)-foso*q(j2,k,3)*q(j2,k,nq)
          qext(4) = (1.+foso)*q(j1,k,4)*q(j1,k,nq)-foso*q(j2,k,4)*q(j2,k,nq)

          qbar(1) = 0.5*(qext(1)+rinf)
          qbar(2) = 0.5*(qext(2)+rinf*uinf)
          qbar(3) = 0.5*(qext(3)+rinf*vinf)
          qbar(4) = 0.5*(qext(4)+einf)
       
          delq(1) = 0.5*(rinf-qext(1))
          delq(2) = 0.5*(rinf*uinf-qext(2))
          delq(3) = 0.5*(rinf*vinf-qext(3))
          delq(4) = 0.5*(einf-qext(4))

          rho = qbar(1)
          rrho = 1./rho
          uu = qbar(2)*rrho
          vv = qbar(3)*rrho
          e = qbar(4)*rrho
          uv2 = 0.5*(uu*uu+vv*vv)
          cjkl = sqrt(ggm1*(e-uv2))
          c2i = 1./(cjkl*cjkl)
          ge = gamma*e - gm1*uv2
          qq = (uu - ug(j,k))*xxn + (vv - vg(j,k))*xyn

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          Z = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
          Zi = 1./Z
          R = 0.5*((1.-bSq)*qq+Z)
          S = 0.5*((1.-bSq)*qq-Z)

c..multiplying by sign(lamda_pa)*inv(X_pa)
          delw(1) = (1.-uv2*gm1*c2i)*sign(1.,qq)*delq(1)
          delw(1) = delw(1) + uu*c2i*gm1*sign(1.,qq)*delq(2)
          delw(1) = delw(1) + vv*c2i*gm1*sign(1.,qq)*delq(3)
          delw(1) = delw(1) - gm1*c2i*sign(1.,qq)*delq(4)
       
          delw(2) = (vv*xxn-uu*xyn)*sign(1.,qq)*delq(1)
          delw(2) = delw(2) + xyn*sign(1.,qq)*delq(2)
          delw(2) = delw(2) - xxn*sign(1.,qq)*delq(3)
       
          delw(3) = -Zi*(uv2*gm1*c2i*S+bSq*qq)*sign(1.,R+bSq*qq)*delq(1)
          delw(3) = delw(3) + (xxn*bSq+S*uu*gm1*c2i)*Zi*
     c                        sign(1.,R+bSq*qq)*delq(2)
          delw(3) = delw(3) + (xyn*bSq+S*vv*gm1*c2i)*Zi*
     c                        sign(1.,R+qq*bSq)*delq(3)
          delw(3) = delw(3) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(4)
       
          delw(4) = Zi*(uv2*gm1*c2i*R+qq*bSq)*sign(1.,S+bSq*qq)*delq(1)
          delw(4) = delw(4) - (xxn*bSq+R*uu*gm1*c2i)*Zi*
     c                        sign(1.,S+bSq*qq)*delq(2)
          delw(4) = delw(4) - (xyn*bSq+R*vv*gm1*c2i)*Zi*
     c                        sign(1.,S+qq*bSq)*delq(3)
          delw(4) = delw(4) + R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(4)
       
c.. multiplying by (X_pa)
          delq(1) = delw(1)+delw(3)+delw(4)
          delq(2) = uu*delw(1)+xyn*delw(2)+(uu+xxn*R/bSq)*delw(3)
          delq(2) = delq(2) + (uu+xxn*S/bSq)*delw(4)
          delq(3) = vv*delw(1)-xxn*delw(2)+(vv+xyn*R/bSq)*delw(3)
          delq(3) = delq(3) + (vv+xyn*S/bSq)*delw(4)
          delq(4) = uv2*delw(1)+(uu*xyn-vv*xxn)*delw(2)
          delq(4) = delq(4) + (ge+R*qq/bSq)*delw(3)
          delq(4) = delq(4) + (ge+S*qq/bSq)*delw(4)

c.. qbar + delq
          q(j,k,1) = (qbar(1)-delq(1))/q(j,k,nq)
          q(j,k,2) = (qbar(2)-delq(2))/q(j,k,nq)
          q(j,k,3) = (qbar(3)-delq(3))/q(j,k,nq)
          q(j,k,4) = (qbar(4)-delq(4))/q(j,k,nq)

        enddo

c..extrapolate everything to other halo cells

        foso = 0.0
C$AD II-LOOP
        do jc = 1,je-js
          if (idir.eq.1) then
            j = je - jc
          elseif (idir.eq.-1) then
            j = js + jc
          endif
          j1 = j  + iadd
          j2 = j1 + iadd

C$AD II-LOOP
          do k = ks,ke
            rho1 =  q(j1,k,1)*q(j1,k,nq)
            rho2 =  q(j2,k,1)*q(j2,k,nq)
            u1 = q(j1,k,2)/q(j1,k,1)  
            u2 = q(j2,k,2)/q(j2,k,1)  
            v1 = q(j1,k,3)/q(j1,k,1)  
            v2 = q(j2,k,3)/q(j2,k,1)  
            p1 = gm1*(q(j1,k,4) -0.5*(q(j1,k,2)**2+q(j1,k,3)**2)
     >                /q(j1,k,1))*q(j1,k,nq)
            p2 = gm1*(q(j2,k,4) -0.5*(q(j2,k,2)**2+q(j2,k,3)**2)
     >                /q(j2,k,1))*q(j2,k,nq)

            q(j,k,1) = ((1.+foso)*rho1 - foso*rho2 )/q(j,k,nq)
            q(j,k,2) = ((1.+foso)*u1 - foso*u2)*q(j,k,1)
            q(j,k,3) = ((1.+foso)*v1 - foso*v2)*q(j,k,1)
            press = (1.+foso)*p1 - foso*p2 
            q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     >                  +q(j,k,3)**2)/q(j,k,1)
          enddo
        enddo

      elseif(iadir.eq.2) then

        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif
        k1 = k  + iadd
        k2 = k1 + iadd
c
        foso = 1.0
C$AD II-LOOP
        do j = js,je
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(yx(j,k)**2+yy(j,k)**2)
          yxn = yx(j,k)*snorm
          yyn = yy(j,k)*snorm
c..first-order extrapolation
          qext(1) = (1.+foso)*q(j,k1,1)*q(j,k1,nq)-foso*q(j,k2,1)*q(j,k2,nq)
          qext(2) = (1.+foso)*q(j,k1,2)*q(j,k1,nq)-foso*q(j,k2,2)*q(j,k2,nq)
          qext(3) = (1.+foso)*q(j,k1,3)*q(j,k1,nq)-foso*q(j,k2,3)*q(j,k2,nq)
          qext(4) = (1.+foso)*q(j,k1,4)*q(j,k1,nq)-foso*q(j,k2,4)*q(j,k2,nq)

          qbar(1) = 0.5*(qext(1)+rinf)
          qbar(2) = 0.5*(qext(2)+rinf*uinf)
          qbar(3) = 0.5*(qext(3)+rinf*vinf)
          qbar(4) = 0.5*(qext(4)+einf)
       
          delq(1) = 0.5*(rinf-qext(1))
          delq(2) = 0.5*(rinf*uinf-qext(2))
          delq(3) = 0.5*(rinf*vinf-qext(3))
          delq(4) = 0.5*(einf-qext(4))

          rho = qbar(1)
          rrho = 1./rho
          uu = qbar(2)*rrho
          vv = qbar(3)*rrho
          e = qbar(4)*rrho
          uv2 = 0.5*(uu*uu+vv*vv)
          cjkl = sqrt(ggm1*(e-uv2))
          c2i = 1./(cjkl*cjkl)
          ge = gamma*e - gm1*uv2
          qq = (uu-ug(j,k))*yxn + (vv-vg(j,k))*yyn

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          Z = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
          Zi = 1./Z
          R = 0.5*((1.-bSq)*qq+Z)
          S = 0.5*((1.-bSq)*qq-Z)

c..multiplying by sign(lamda_pa)*inv(X_pa)
          delw(1) = (1.-uv2*gm1*c2i)*sign(1.,qq)*delq(1)
          delw(1) = delw(1) + uu*c2i*gm1*sign(1.,qq)*delq(2)
          delw(1) = delw(1) + vv*c2i*gm1*sign(1.,qq)*delq(3)
          delw(1) = delw(1) - gm1*c2i*sign(1.,qq)*delq(4)
       
          delw(2) = (vv*yxn-uu*yyn)*sign(1.,qq)*delq(1)
          delw(2) = delw(2) + yyn*sign(1.,qq)*delq(2)
          delw(2) = delw(2) - yxn*sign(1.,qq)*delq(3)
       
          delw(3) = -Zi*(uv2*gm1*c2i*S+bSq*qq)*sign(1.,R+bSq*qq)*delq(1)
          delw(3) = delw(3) + (yxn*bSq+S*uu*gm1*c2i)*Zi*
     c                        sign(1.,R+bSq*qq)*delq(2)
          delw(3) = delw(3) + (yyn*bSq+S*vv*gm1*c2i)*Zi*
     c                        sign(1.,R+qq*bSq)*delq(3)
          delw(3) = delw(3) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(4)
       
          delw(4) = Zi*(uv2*gm1*c2i*R+qq*bSq)*sign(1.,S+bSq*qq)*delq(1)
          delw(4) = delw(4) - (yxn*bSq+R*uu*gm1*c2i)*Zi*
     c                        sign(1.,S+bSq*qq)*delq(2)
          delw(4) = delw(4) - (yyn*bSq+R*vv*gm1*c2i)*Zi*
     c                        sign(1.,S+qq*bSq)*delq(3)
          delw(4) = delw(4) + R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(4)
       
c.. multiplying by (X_pa)
          delq(1) = delw(1)+delw(3)+delw(4)
          delq(2) = uu*delw(1)+yyn*delw(2)+(uu+yxn*R/bSq)*delw(3)
          delq(2) = delq(2) + (uu+yxn*S/bSq)*delw(4)
          delq(3) = vv*delw(1)-yxn*delw(2)+(vv+yyn*R/bSq)*delw(3)
          delq(3) = delq(3) + (vv+yyn*S/bSq)*delw(4)
          delq(4) = uv2*delw(1)+(uu*yyn-vv*yxn)*delw(2)
          delq(4) = delq(4) + (ge+R*qq/bSq)*delw(3)
          delq(4) = delq(4) + (ge+S*qq/bSq)*delw(4)

c.. qbar - delq
          q(j,k,1) = (qbar(1)-delq(1))/q(j,k,nq)
          q(j,k,2) = (qbar(2)-delq(2))/q(j,k,nq)
          q(j,k,3) = (qbar(3)-delq(3))/q(j,k,nq)
          q(j,k,4) = (qbar(4)-delq(4))/q(j,k,nq)

        enddo

c..extrapolate everything to other halo cells

        foso = 0.0
C$AD II-LOOP
        do kc = 1,ke-ks
          if (idir.eq.2) then
            k = ke - kc
          elseif (idir.eq.-2) then
            k = ks + kc
          endif
          k1 = k  + iadd
          k2 = k1 + iadd

C$AD II-LOOP
          do j = js,je
            rho1 =  q(j,k1,1)*q(j,k1,nq)
            rho2 =  q(j,k2,1)*q(j,k2,nq)
            u1 = q(j,k1,2)/q(j,k1,1)  
            u2 = q(j,k2,2)/q(j,k2,1)  
            v1 = q(j,k1,3)/q(j,k1,1)  
            v2 = q(j,k2,3)/q(j,k2,1)  
            p1 = gm1*(q(j,k1,4) -0.5*(q(j,k1,2)**2+q(j,k1,3)**2)
     >                /q(j,k1,1))*q(j,k1,nq)
            p2 = gm1*(q(j,k2,4) -0.5*(q(j,k2,2)**2+q(j,k2,3)**2)
     >                /q(j,k2,1))*q(j,k2,nq)

            q(j,k,1) = ((1.+foso)*rho1 - foso*rho2 )/q(j,k,nq)
            q(j,k,2) = ((1.+foso)*u1 - foso*u2)*q(j,k,1)
            q(j,k,3) = ((1.+foso)*v1 - foso*v2)*q(j,k,1)
            press = (1.+foso)*p1 - foso*p2 
            q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     >                  +q(j,k,3)**2)/q(j,k,1)
          enddo
        enddo

      endif

      return
      end

c*********************************************************************
      subroutine rbcout_trkl(q,dq,x,y,xx,xy,yx,yy,ug,vg,bt,jd,kd,
     <                      js,je,ks,ke,idir)
c
c     Outflow boundary condition with Turkel Preconditioning 
c     Characteristic extrapolation based on Riemann invariants
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,js,je,ks,ke,idir
      real  q(jd,kd,nq),dq(jd,kd,nv)
      real  xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd),
     $      ug(jd,kd),vg(jd,kd),x(jd,kd),y(jd,kd),bt(jd,kd)

      ! local variables
      real foso
      real rhoext,uext,vext,eext,pext
      real xxn,xyn,yxn,yyn,snorm
      real rho,rrho,uu,vv,e,uv2,cjkl,c2i,ge,qq
      real Z,R,S,Zi,bSq
      real u1,u2,v1,v2
      real drho1,drho2,du1,du2,dv1,dv2,dp1,dp2,dpress
      integer k,k1,k2,j,j1,j2,jc,kc,iadd,iadir
      real qext(4),qbar(4)
      real sext(4),sbar(4),dels(4),delw(4)
      
c***  first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then

        if (idir.eq.1) then
          j = je
        elseif (idir.eq.-1) then
          j = js
        endif
        j1 = j  + iadd
        j2 = j1 + iadd
c
        foso = 1.0
C$AD II-LOOP
        do k = ks,ke
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(xx(j,k)**2+xy(j,k)**2)
          xxn = xx(j,k)*snorm
          xyn = xy(j,k)*snorm
c..first-order extrapolation
          qext(1) = (1.+foso)*q(j1,k,1)*q(j1,k,nq)-foso*q(j2,k,1)*q(j2,k,nq)
          qext(2) = (1.+foso)*q(j1,k,2)*q(j1,k,nq)-foso*q(j2,k,2)*q(j2,k,nq)
          qext(3) = (1.+foso)*q(j1,k,3)*q(j1,k,nq)-foso*q(j2,k,3)*q(j2,k,nq)
          qext(4) = (1.+foso)*q(j1,k,4)*q(j1,k,nq)-foso*q(j2,k,4)*q(j2,k,nq)

          sext(1) = (1.+foso)*dq(j1,k,1)*q(j1,k,nq)-foso*dq(j2,k,1)*q(j2,k,nq)
          sext(2) = (1.+foso)*dq(j1,k,2)*q(j1,k,nq)-foso*dq(j2,k,2)*q(j2,k,nq)
          sext(3) = (1.+foso)*dq(j1,k,3)*q(j1,k,nq)-foso*dq(j2,k,3)*q(j2,k,nq)
          sext(4) = (1.+foso)*dq(j1,k,4)*q(j1,k,nq)-foso*dq(j2,k,4)*q(j2,k,nq)

          qbar(1) = 0.5*(qext(1)+rinf)
          qbar(2) = 0.5*(qext(2)+rinf*uinf)
          qbar(3) = 0.5*(qext(3)+rinf*vinf)
          qbar(4) = 0.5*(qext(4)+einf)

          sbar(1) = 0.5*sext(1)
          sbar(2) = 0.5*sext(2)
          sbar(3) = 0.5*sext(3)
          sbar(4) = 0.5*sext(4)
       
          dels(1) = -0.5*sext(1) 
          dels(2) = -0.5*sext(2)
          dels(3) = -0.5*sext(3)
          dels(4) = -0.5*sext(4)

          rho = qbar(1)
          rrho = 1./rho
          uu = qbar(2)*rrho
          vv = qbar(3)*rrho
          e = qbar(4)*rrho
          uv2 = 0.5*(uu*uu+vv*vv)
          cjkl = sqrt(ggm1*(e-uv2))
          c2i = 1./(cjkl*cjkl)
          ge = gamma*e - gm1*uv2
          qq = (uu - ug(j,k))*xxn + (vv - vg(j,k))*xyn

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          Z = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
          Zi = 1./Z
          R = 0.5*((1.-bSq)*qq+Z)
          S = 0.5*((1.-bSq)*qq-Z)

c..multiplying by sign(lamda_pa)*inv(X_pa)
          delw(1) = (1.-uv2*gm1*c2i)*sign(1.,qq)*dels(1)
          delw(1) = delw(1) + uu*c2i*gm1*sign(1.,qq)*dels(2)
          delw(1) = delw(1) + vv*c2i*gm1*sign(1.,qq)*dels(3)
          delw(1) = delw(1) - gm1*c2i*sign(1.,qq)*dels(4)
       
          delw(2) = (vv*xxn-uu*xyn)*sign(1.,qq)*dels(1)
          delw(2) = delw(2) + xyn*sign(1.,qq)*dels(2)
          delw(2) = delw(2) - xxn*sign(1.,qq)*dels(3)
       
          delw(3) = -Zi*(uv2*gm1*c2i*S+bSq*qq)*sign(1.,R+bSq*qq)*dels(1)
          delw(3) = delw(3) + (xxn*bSq+S*uu*gm1*c2i)*Zi*
     c                        sign(1.,R+bSq*qq)*dels(2)
          delw(3) = delw(3) + (xyn*bSq+S*vv*gm1*c2i)*Zi*
     c                        sign(1.,R+qq*bSq)*dels(3)
          delw(3) = delw(3) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*dels(4)
       
          delw(4) = Zi*(uv2*gm1*c2i*R+qq*bSq)*sign(1.,S+bSq*qq)*dels(1)
          delw(4) = delw(4) - (xxn*bSq+R*uu*gm1*c2i)*Zi*
     c                        sign(1.,S+bSq*qq)*dels(2)
          delw(4) = delw(4) - (xyn*bSq+R*vv*gm1*c2i)*Zi*
     c                        sign(1.,S+qq*bSq)*dels(3)
          delw(4) = delw(4) + R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*dels(4)
       
c.. multiplying by (X_pa)
          dels(1) = delw(1)+delw(3)+delw(4)
          dels(2) = uu*delw(1)+xyn*delw(2)+(uu+xxn*R/bSq)*delw(3)
          dels(2) = dels(2) + (uu+xxn*S/bSq)*delw(4)
          dels(3) = vv*delw(1)-xxn*delw(2)+(vv+xyn*R/bSq)*delw(3)
          dels(3) = dels(3) + (vv+xyn*S/bSq)*delw(4)
          dels(4) = uv2*delw(1)+(uu*xyn-vv*xxn)*delw(2)
          dels(4) = dels(4) + (ge+R*qq/bSq)*delw(3)
          dels(4) = dels(4) + (ge+S*qq/bSq)*delw(4)

c.. sbar + dels
          dq(j,k,1) = (sbar(1)-dels(1))/q(j,k,nq)
          dq(j,k,2) = (sbar(2)-dels(2))/q(j,k,nq)
          dq(j,k,3) = (sbar(3)-dels(3))/q(j,k,nq)
          dq(j,k,4) = (sbar(4)-dels(4))/q(j,k,nq)

        enddo

c..extrapolate everything to other halo cells

        foso = 0.0
C$AD II-LOOP
        do jc = 1,je-js
          if (idir.eq.1) then
            j = je - jc
          elseif (idir.eq.-1) then
            j = js + jc
          endif
          j1 = j  + iadd
          j2 = j1 + iadd

C$AD II-LOOP
          do k = ks,ke
            drho1 =  dq(j1,k,1)*q(j1,k,nq)
            drho2 =  dq(j2,k,1)*q(j2,k,nq)

            u1 = q(j1,k,2)/q(j1,k,1)  
            u2 = q(j2,k,2)/q(j2,k,1)  
            du1 = dq(j1,k,2)/q(j1,k,1) - q(j1,k,2)*dq(j1,k,1)/q(j1,k,1)**2
            du2 = dq(j2,k,2)/q(j2,k,1) - q(j2,k,2)*dq(j2,k,1)/q(j2,k,1)**2

            v1 = q(j1,k,3)/q(j1,k,1)  
            v2 = q(j2,k,3)/q(j2,k,1)  
            dv1 = dq(j1,k,3)/q(j1,k,1) - q(j1,k,3)*dq(j1,k,1)/q(j1,k,1)**2
            dv2 = dq(j2,k,3)/q(j2,k,1) - q(j2,k,3)*dq(j2,k,1)/q(j2,k,1)**2

            dp1 = gm1*(dq(j1,k,4) + 
     >            0.5*(q(j1,k,2)**2+q(j1,k,3)**2)*dq(j1,k,1)/q(j1,k,1)**2 -
     >            (q(j1,k,2)*dq(j1,k,2) + q(j1,k,3)*dq(j1,k,3))/q(j1,k,1))*
     >            q(j1,k,nq)
            dp2 = gm1*(dq(j2,k,4) + 
     >            0.5*(q(j2,k,2)**2+q(j2,k,3)**2)*dq(j2,k,1)/q(j2,k,1)**2 -
     >            (q(j2,k,2)*dq(j2,k,2) + q(j2,k,3)*dq(j2,k,3))/q(j2,k,1))*
     >            q(j2,k,nq)

            dq(j,k,1) = ((1.+foso)*drho1 - foso*drho2 )/q(j,k,nq)
            dq(j,k,2) = ((1.+foso)*du1 - foso*du2)*q(j,k,1) + 
     >                 ((1.+foso)*u1 - foso*u2)*dq(j,k,1)
            dq(j,k,3) = ((1.+foso)*dv1 - foso*dv2)*q(j,k,1) + 
     >                 ((1.+foso)*v1 - foso*v2)*dq(j,k,1)
            dpress = (1.+foso)*dp1 - foso*dp2 
            dq(j,k,4) = dpress/(gm1*q(j,k,nq)) -
     >            0.5*(q(j,k,2)**2+q(j,k,3)**2)*dq(j,k,1)/q(j,k,1)**2 +
     >            (q(j,k,2)*dq(j,k,2) + q(j,k,3)*dq(j,k,3))/q(j,k,1)
          enddo
        enddo

      elseif(iadir.eq.2) then

        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif
        k1 = k  + iadd
        k2 = k1 + iadd
c
        foso = 1.0
C$AD II-LOOP
        do j = js,je
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(yx(j,k)**2+yy(j,k)**2)
          yxn = yx(j,k)*snorm
          yyn = yy(j,k)*snorm
c..first-order extrapolation
          qext(1) = (1.+foso)*q(j,k1,1)*q(j,k1,nq)-foso*q(j,k2,1)*q(j,k2,nq)
          qext(2) = (1.+foso)*q(j,k1,2)*q(j,k1,nq)-foso*q(j,k2,2)*q(j,k2,nq)
          qext(3) = (1.+foso)*q(j,k1,3)*q(j,k1,nq)-foso*q(j,k2,3)*q(j,k2,nq)
          qext(4) = (1.+foso)*q(j,k1,4)*q(j,k1,nq)-foso*q(j,k2,4)*q(j,k2,nq)

          sext(1) = (1.+foso)*dq(j,k1,1)*q(j,k1,nq)-foso*dq(j,k2,1)*q(j,k2,nq)
          sext(2) = (1.+foso)*dq(j,k1,2)*q(j,k1,nq)-foso*dq(j,k2,2)*q(j,k2,nq)
          sext(3) = (1.+foso)*dq(j,k1,3)*q(j,k1,nq)-foso*dq(j,k2,3)*q(j,k2,nq)
          sext(4) = (1.+foso)*dq(j,k1,4)*q(j,k1,nq)-foso*dq(j,k2,4)*q(j,k2,nq)

          qbar(1) = 0.5*(qext(1)+rinf)
          qbar(2) = 0.5*(qext(2)+rinf*uinf)
          qbar(3) = 0.5*(qext(3)+rinf*vinf)
          qbar(4) = 0.5*(qext(4)+einf)

          sbar(1) = 0.5*sext(1)
          sbar(2) = 0.5*sext(2)
          sbar(3) = 0.5*sext(3)
          sbar(4) = 0.5*sext(4)
       
          dels(1) = -0.5*sext(1)
          dels(2) = -0.5*sext(2)
          dels(3) = -0.5*sext(3)
          dels(4) = -0.5*sext(4)

          rho = qbar(1)
          rrho = 1./rho
          uu = qbar(2)*rrho
          vv = qbar(3)*rrho
          e = qbar(4)*rrho
          uv2 = 0.5*(uu*uu+vv*vv)
          cjkl = sqrt(ggm1*(e-uv2))
          c2i = 1./(cjkl*cjkl)
          ge = gamma*e - gm1*uv2
          qq = (uu-ug(j,k))*yxn + (vv-vg(j,k))*yyn

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          Z = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
          Zi = 1./Z
          R = 0.5*((1.-bSq)*qq+Z)
          S = 0.5*((1.-bSq)*qq-Z)

c..multiplying by sign(lamda_pa)*inv(X_pa)
          delw(1) = (1.-uv2*gm1*c2i)*sign(1.,qq)*dels(1)
          delw(1) = delw(1) + uu*c2i*gm1*sign(1.,qq)*dels(2)
          delw(1) = delw(1) + vv*c2i*gm1*sign(1.,qq)*dels(3)
          delw(1) = delw(1) - gm1*c2i*sign(1.,qq)*dels(4)
       
          delw(2) = (vv*yxn-uu*yyn)*sign(1.,qq)*dels(1)
          delw(2) = delw(2) + yyn*sign(1.,qq)*dels(2)
          delw(2) = delw(2) - yxn*sign(1.,qq)*dels(3)
       
          delw(3) = -Zi*(uv2*gm1*c2i*S+bSq*qq)*sign(1.,R+bSq*qq)*dels(1)
          delw(3) = delw(3) + (yxn*bSq+S*uu*gm1*c2i)*Zi*
     c                        sign(1.,R+bSq*qq)*dels(2)
          delw(3) = delw(3) + (yyn*bSq+S*vv*gm1*c2i)*Zi*
     c                        sign(1.,R+qq*bSq)*dels(3)
          delw(3) = delw(3) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*dels(4)
       
          delw(4) = Zi*(uv2*gm1*c2i*R+qq*bSq)*sign(1.,S+bSq*qq)*dels(1)
          delw(4) = delw(4) - (yxn*bSq+R*uu*gm1*c2i)*Zi*
     c                        sign(1.,S+bSq*qq)*dels(2)
          delw(4) = delw(4) - (yyn*bSq+R*vv*gm1*c2i)*Zi*
     c                        sign(1.,S+qq*bSq)*dels(3)
          delw(4) = delw(4) + R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*dels(4)
       
c.. multiplying by (X_pa)
          dels(1) = delw(1)+delw(3)+delw(4)
          dels(2) = uu*delw(1)+yyn*delw(2)+(uu+yxn*R/bSq)*delw(3)
          dels(2) = dels(2) + (uu+yxn*S/bSq)*delw(4)
          dels(3) = vv*delw(1)-yxn*delw(2)+(vv+yyn*R/bSq)*delw(3)
          dels(3) = dels(3) + (vv+yyn*S/bSq)*delw(4)
          dels(4) = uv2*delw(1)+(uu*yyn-vv*yxn)*delw(2)
          dels(4) = dels(4) + (ge+R*qq/bSq)*delw(3)
          dels(4) = dels(4) + (ge+S*qq/bSq)*delw(4)

c.. qbar - dels
          dq(j,k,1) = (sbar(1)-dels(1))/q(j,k,nq)
          dq(j,k,2) = (sbar(2)-dels(2))/q(j,k,nq)
          dq(j,k,3) = (sbar(3)-dels(3))/q(j,k,nq)
          dq(j,k,4) = (sbar(4)-dels(4))/q(j,k,nq)

        enddo

c..extrapolate everything to other halo cells

        foso = 0.0
C$AD II-LOOP
        do kc = 1,ke-ks
          if (idir.eq.2) then
            k = ke - kc
          elseif (idir.eq.-2) then
            k = ks + kc
          endif
          k1 = k  + iadd
          k2 = k1 + iadd

C$AD II-LOOP
          do j = js,je
            drho1 =  dq(j,k1,1)*q(j,k1,nq)
            drho2 =  dq(j,k2,1)*q(j,k2,nq)

            u1 = q(j,k1,2)/q(j,k1,1)  
            u2 = q(j,k2,2)/q(j,k2,1)  
            du1 = dq(j,k1,2)/q(j,k1,1) - q(j,k1,2)*dq(j,k1,1)/q(j,k1,1)**2
            du2 = dq(j,k2,2)/q(j,k2,1) - q(j,k2,2)*dq(j,k2,1)/q(j,k2,1)**2

            v1 = q(j,k1,3)/q(j,k1,1)  
            v2 = q(j,k2,3)/q(j,k2,1)  
            dv1 = dq(j,k1,3)/q(j,k1,1) - q(j,k1,3)*dq(j,k1,1)/q(j,k1,1)**2
            dv2 = dq(j,k2,3)/q(j,k2,1) - q(j,k2,3)*dq(j,k2,1)/q(j,k2,1)**2

            dp1 = gm1*(dq(j,k1,4) + 
     >            0.5*(q(j,k1,2)**2+q(j,k1,3)**2)*dq(j,k1,1)/q(j,k1,1)**2 -
     >            (q(j,k1,2)*dq(j,k1,2) + q(j,k1,3)*dq(j,k1,3))/q(j,k1,1))*
     >            q(j,k1,nq)
            dp2 = gm1*(dq(j,k2,4) + 
     >            0.5*(q(j,k2,2)**2+q(j,k2,3)**2)*dq(j,k2,1)/q(j,k2,1)**2 -
     >            (q(j,k2,2)*dq(j,k2,2) + q(j,k2,3)*dq(j,k2,3))/q(j,k2,1))*
     >            q(j,k2,nq)

            dq(j,k,1) = ((1.+foso)*drho1 - foso*drho2 )/q(j,k,nq)
            dq(j,k,2) = ((1.+foso)*du1 - foso*du2)*q(j,k,1) + 
     >                 ((1.+foso)*u1 - foso*u2)*dq(j,k,1)
            dq(j,k,3) = ((1.+foso)*dv1 - foso*dv2)*q(j,k,1) + 
     >                 ((1.+foso)*v1 - foso*v2)*dq(j,k,1)
            dpress = (1.+foso)*dp1 - foso*dp2 
            dq(j,k,4) = dpress/(gm1*q(j,k,nq)) -
     >            0.5*(q(j,k,2)**2+q(j,k,3)**2)*dq(j,k,1)/q(j,k,1)**2 +
     >            (q(j,k,2)*dq(j,k,2) + q(j,k,3)*dq(j,k,3))/q(j,k,1)
          enddo
        enddo

      endif

      return
      end

c***********************************************************************
      subroutine bctany(q,uwall,vwall,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
c
c  inviscid : tangency b.c. at a eta equal constant surface
c  viscous  : no slip condition
c
c  presently, set up for surface at k = kreq
c             direction of freestream = idir   
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer jd,kd,ks,ke,idir
      real q(jd,kd,nq)
      real uwall(mdim), vwall(mdim)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)

      ! local variables
      real :: u(mdim), v(mdim)
      real :: scale(mdim)

      integer js,je,k,k1,k2,j,iadd,iadir
      real foso,t,u1,u2,v1
      
c***  first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bctany'
      elseif(iadir.eq.2) then
        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif
        k1 = k + iadd
        k2 = k1 + iadd
        foso = 1.0

        if( invisc ) then

c..for inviscid:
c  extrapolate contravariant velocities to the surface
c  solve for the surface cartesian velocity components

cjdb..note: really we are extrapolating the physical values!

C$AD II-LOOP
          do 11 j = js,je
            u1 = (q(j,k1,2)*xx(j,k)+q(j,k1,3)*xy(j,k))/q(j,k1,1)
            u2 = (q(j,k2,2)*xx(j,k)+q(j,k2,3)*xy(j,k))/q(j,k2,1)
            u(j)  = (1.+foso)*u1 - foso*u2
            v1 = (q(j,k1,2)*yx(j,k)+q(j,k1,3)*yy(j,k))/q(j,k1,1)
            v(j)  = ug(j,k)*yx(j,k)   + vg(j,k)*yy(j,k)   + 
     &              ug(j,k1)*yx(j,k1) + vg(j,k1)*yy(j,k1) - v1
   11     continue
        else
c
c..for viscous flows the no slip condition at the surface requires
c  setting the contravariant velocities u = 0., & v = 0. for
c  no suction or injection at the surface. this satisfies tangency
c  at the walls  
c..grs
c  solve for the surface cartesian velocity components
c
C$AD II-LOOP
          do 55 j = js, je
            u1 = (q(j,k1,2)*xx(j,k)+q(j,k1,3)*xy(j,k))/q(j,k1,1)
            v1 = (q(j,k1,2)*yx(j,k)+q(j,k1,3)*yy(j,k))/q(j,k1,1)
            u(j) = ug(j,k)*xx(j,k)   + vg(j,k)*xy(j,k)   +
     &             ug(j,k1)*xx(j,k1) + vg(j,k1)*xy(j,k1) - u1
            v(j) = ug(j,k)*yx(j,k)   + vg(j,k)*yy(j,k)   +
     &             ug(j,k1)*yx(j,k1) + vg(j,k1)*yy(j,k1) - v1
   55     continue
        endif
        call uv(js,je,k,xx,xy,yx,yy,u,v,uwall,vwall,jd,kd)

      endif
c
      return
      end

c***********************************************************************
      subroutine rbctany(q,s,duwall,dvwall,xx,xy,yx,yy,ug,vg,
     &                   jd,kd,js,je,ks,ke,idir)
c
c  inviscid : tangency b.c. at a eta equal constant surface
c  viscous  : no slip condition
c
c  presently, set up for surface at k = kreq
c             direction of freestream = idir   
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer jd,kd,ks,ke,idir
      real q(jd,kd,nq),s(jd,kd,nv)
      real duwall(mdim), dvwall(mdim)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)

      ! local variables
      real :: du(mdim), dv(mdim)
      real :: scale(mdim)

      integer js,je,k,k1,k2,j,iadd,iadir
      real foso,t,du1,du2,dv1
      
c***  first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bctany'
      elseif(iadir.eq.2) then
        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif
        k1 = k + iadd
        k2 = k1 + iadd
        foso = 1.0

        if( invisc ) then

c..for inviscid:
c  extrapolate contravariant velocities to the surface
c  solve for the surface cartesian velocity components

cjdb..note: really we are extrapolating the physical values!

C$AD II-LOOP
          do 11 j = js,je
            du1 = (s(j,k1,2)*xx(j,k)+s(j,k1,3)*xy(j,k))/q(j,k1,1) - 
     &            (q(j,k1,2)*xx(j,k)+q(j,k1,3)*xy(j,k))*s(j,k1,1)/
     &             q(j,k1,1)**2
            du2 = (s(j,k2,2)*xx(j,k)+s(j,k2,3)*xy(j,k))/q(j,k2,1) - 
     &            (q(j,k2,2)*xx(j,k)+q(j,k2,3)*xy(j,k))*s(j,k2,1)/
     &             q(j,k2,1)**2
            du(j) = (1.+foso)*du1 - foso*du2
            dv1 = (s(j,k1,2)*yx(j,k)+s(j,k1,3)*yy(j,k))/q(j,k1,1) - 
     &            (q(j,k1,2)*yx(j,k)+q(j,k1,3)*yy(j,k))*s(j,k1,1)/
     &             q(j,k1,1)**2
            dv(j) = -dv1 !need to add contribution from dug and dvg
   11     continue
        else
c
c..for viscous flows the no slip condition at the surface requires
c  setting the contravariant velocities u = 0., & v = 0. for
c  no suction or injection at the surface. this satisfies tangency
c  at the walls  
c..grs
c  solve for the surface cartesian velocity components
c
C$AD II-LOOP
          do 55 j = js, je
            du1 = (s(j,k1,2)*xx(j,k)+s(j,k1,3)*xy(j,k))/q(j,k1,1) - 
     &            (q(j,k1,2)*xx(j,k)+q(j,k1,3)*xy(j,k))*s(j,k1,1)/
     &             q(j,k1,1)**2
            dv1 = (s(j,k1,2)*yx(j,k)+s(j,k1,3)*yy(j,k))/q(j,k1,1) - 
     &            (q(j,k1,2)*yx(j,k)+q(j,k1,3)*yy(j,k))*s(j,k1,1)/
     &             q(j,k1,1)**2
            du(j) = -du1 
            dv(j) = -dv1 !need to add contribution from dug and dvg
   55     continue
        endif
        call ruv(js,je,k,xx,xy,yx,yy,du,dv,duwall,dvwall,jd,kd)

      endif
c
      return
      end

c***********************************************************************
      subroutine bcwake(q,x,y,jd,kd,js,je,ks,ke,idir)
c
c  for complete cgrid   
c  treatment for in the trailing wake at k = 1
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real q(jd,kd,nq)
      real x(jd,kd), y(jd,kd)

      ! local variables
      
      integer k,k1,j,jj,kc,iadd,iadir
      real scale1,scale2
      
c**   first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwake'
      elseif(iadir.eq.2) then

        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
            k1 = ke + kc
          elseif (idir.eq.-2) then
            k = ks + kc - 1
            k1 = ks - kc
          endif

c..in the wake, average points above and below

          do j = js, je
            jj = jd - j + 1
            scale1= q(jj,k1,nq)/q(j,k,nq)
            scale2= q(j,k1,nq)/q(jj,k,nq)
            q(j,k,1)  = q(jj,k1,1)*scale1
            q(j,k,2)  = q(jj,k1,2)*scale1
            q(j,k,3)  = q(jj,k1,3)*scale1
            q(j,k,4)  = q(jj,k1,4)*scale1
            q(jj,k,1) = q(j,k1,1)*scale2
            q(jj,k,2) = q(j,k1,2)*scale2
            q(jj,k,3) = q(j,k1,3)*scale2
            q(jj,k,4) = q(j,k1,4)*scale2
          enddo
       enddo
      
      endif
c
      return
      end

c***********************************************************************
      subroutine rbcwake(q,s,jd,kd,js,je,ks,ke,idir)
c
c  for complete cgrid   
c  treatment for in the trailing wake at k = 1
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real q(jd,kd,nq),s(jd,kd,nv)
      real x(jd,kd), y(jd,kd)

      ! local variables
      
      integer k,k1,j,jj,kc,iadd,iadir
      real scale1,scale2
      
c**   first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in rbcwake'
      elseif(iadir.eq.2) then

        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
            k1 = ke + kc
          elseif (idir.eq.-2) then
            k = ks + kc - 1
            k1 = ks - kc
          endif

c..in the wake, average points above and below

          do 30 j = js, je
            jj = jd - j + 1
            scale1= q(jj,k1,nq)/q(j,k,nq)
            scale2= q(j,k1,nq)/q(jj,k,nq)
            s(j,k,1)  = s(jj,k1,1)*scale1
            s(j,k,2)  = s(jj,k1,2)*scale1
            s(j,k,3)  = s(jj,k1,3)*scale1
            s(j,k,4)  = s(jj,k1,4)*scale1
            s(jj,k,1) = s(j,k1,1)*scale2
            s(jj,k,2) = s(j,k1,2)*scale2
            s(jj,k,3) = s(j,k1,3)*scale2
            s(jj,k,4) = s(j,k1,4)*scale2
   30     continue
       enddo
      
      endif
c
      return
      end

c***********************************************************************
      subroutine zero_residuals(s,jd,kd,js,je,ks,ke,idir)
c
c  zero the residual 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real s(jd,kd,nv)

      ! local variables
      
      integer j,k,n,iadd,iadir
      
c**   first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      do j = js,je
        do k = ks,ke
          do n = 1,4
            s(j,k,n) = 0.0
          enddo
        enddo
      enddo
 
      return
      end

c***********************************************************************
      subroutine bcwall(q,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,
     &                          idir,invsc)
c
c  generic solid wall boundary conditions.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      logical invsc
      real q(jd,kd,nq)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real ug(jd,kd),vg(jd,kd)
      
      ! local variables
      real :: uwall(mdim),vwall(mdim)

      integer k,k1,k2,k3,kc,ihigh,j,iadd,iadir
      real foso,rj,us,vs,ue,ve,t,scal,rscal,ajacinv
      real rho1,rho2,u1,u2,v1,v2,p1,p2,press
      logical tmp

c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

c...setting wind tunnel walls to be inviscid
      tmp = invisc
      !am invisc = invsc

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwall'
      elseif(iadir.eq.2) then

        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif

        k1 = k  + iadd
        k2 = k1 + iadd

        if (invsc) then
          foso = 1.0
        else
          foso = 0.0
        endif

c..compute surface velocities

        call bctany(q,uwall,vwall,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)

c..extrapolate pressure and density to the surface

        do 11 j = js,je
          rho1 =  q(j,k1,1)*q(j,k1,nq)
          rho2 =  q(j,k2,1)*q(j,k2,nq)

          p1 = gm1*(q(j,k1,4) -0.5*(q(j,k1,2)**2+q(j,k1,3)**2)
     >              /q(j,k1,1))*q(j,k1,nq)
          p2 = gm1*(q(j,k2,4) -0.5*(q(j,k2,2)**2+q(j,k2,3)**2)
     >              /q(j,k2,1))*q(j,k2,nq)

          q(j,k,1) = ((1.+foso)*rho1 - foso*rho2 )/q(j,k,nq)
          q(j,k,2) = uwall(j)*q(j,k,1)
          q(j,k,3) = vwall(j)*q(j,k,1)
          press = (1.+foso)*p1 - foso*p2 
          q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     >                +q(j,k,3)**2)/q(j,k,1)
 11     continue

c..extrapolate everything to other halo cells

        do kc = 1,ke-ks
          if (idir.eq.2) then
            k = ke - kc
          elseif (idir.eq.-2) then
            k = ks + kc
          endif
          k1 = k  + iadd
          k2 = k1 + iadd

          foso = 1.0

          do j = js,je
            rho1 =  q(j,k1,1)*q(j,k1,nq)
            rho2 =  q(j,k2,1)*q(j,k2,nq)
            u1 = q(j,k1,2)/q(j,k1,1)  
            u2 = q(j,k2,2)/q(j,k2,1)  
            v1 = q(j,k1,3)/q(j,k1,1)  
            v2 = q(j,k2,3)/q(j,k2,1)  
            p1 = gm1*(q(j,k1,4) -0.5*(q(j,k1,2)**2+q(j,k1,3)**2)
     >                /q(j,k1,1))*q(j,k1,nq)
            p2 = gm1*(q(j,k2,4) -0.5*(q(j,k2,2)**2+q(j,k2,3)**2)
     >                /q(j,k2,1))*q(j,k2,nq)

            q(j,k,1) = ((1.+foso)*rho1 - foso*rho2 )/q(j,k,nq)
            q(j,k,2) = ((1.+foso)*u1 - foso*u2)*q(j,k,1)
            q(j,k,3) = ((1.+foso)*v1 - foso*v2)*q(j,k,1)
            press = (1.+foso)*p1 - foso*p2 
            q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     >                  +q(j,k,3)**2)/q(j,k,1)
          enddo
        enddo

c..slowly turn on the wall boundary condition 
c  over 30 time steps
        t = float(istep0)/30.
        if(t.gt.1..or.iread.gt.0) t=1.
        scal = (10.-15.*t+6.*t*t)*t**3
        rscal = 1.-scal
        do k=ks,ke
          do j=js,je
            ajacinv = 1./q(j,k,nq)
            q(j,k,1) = q(j,k,1)*scal+rscal*rinf*ajacinv
            q(j,k,2) = q(j,k,2)*scal+rscal*rinf*uinf*ajacinv
            q(j,k,3) = q(j,k,3)*scal+rscal*rinf*vinf*ajacinv
            q(j,k,4) = q(j,k,4)*scal+rscal*einf*ajacinv
          enddo
        enddo
      
      endif

c...resetting the inviscid to original value
      !am invisc = tmp

      return
      end

c***********************************************************************
      subroutine rbcwall(q,s,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,
     &                          idir,invsc)
c
c  generic solid wall boundary conditions.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      logical invsc
      real q(jd,kd,nq),s(jd,kd,nv)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real ug(jd,kd),vg(jd,kd)
      
      ! local variables
      real :: uwall(mdim),vwall(mdim)
      real :: duwall(mdim),dvwall(mdim)

      integer k,k1,k2,k3,kc,ihigh,j,iadd,iadir
      real foso,rj,us,vs,ue,ve,t,scal,rscal,ajacinv
      real u1,u2,v1,v2
      real drho1,drho2,du1,du2,dv1,dv2,dp1,dp2,dpress
      logical tmp

c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

c...setting wind tunnel walls to be inviscid
      tmp = invisc
      !am invisc = invsc

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwall'
      elseif(iadir.eq.2) then

        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif

        k1 = k  + iadd
        k2 = k1 + iadd

        if (invsc) then
          foso = 1.0
        else
          foso = 0.0
        endif

c..compute surface velocities

        call bctany(q,uwall,vwall,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
        call rbctany(q,s,duwall,dvwall,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)

c..extrapolate pressure and density to the surface

        do 11 j = js,je
          drho1 =  s(j,k1,1)*q(j,k1,nq)
          drho2 =  s(j,k2,1)*q(j,k2,nq)

          dp1 = gm1*(s(j,k1,4) + 
     >          0.5*(q(j,k1,2)**2+q(j,k1,3)**2)*s(j,k1,1)/q(j,k1,1)**2 -
     >          (q(j,k1,2)*s(j,k1,2) + q(j,k1,3)*s(j,k1,3))/q(j,k1,1))*
     >          q(j,k1,nq)
          dp2 = gm1*(s(j,k2,4) + 
     >          0.5*(q(j,k2,2)**2+q(j,k2,3)**2)*s(j,k2,1)/q(j,k2,1)**2 -
     >          (q(j,k2,2)*s(j,k2,2) + q(j,k2,3)*s(j,k2,3))/q(j,k2,1))*
     >          q(j,k2,nq)

          s(j,k,1) = ((1.+foso)*drho1 - foso*drho2 )/q(j,k,nq)
          s(j,k,2) = duwall(j)*q(j,k,1) + uwall(j)*s(j,k,1)
          s(j,k,3) = dvwall(j)*q(j,k,1) + vwall(j)*s(j,k,1)
          dpress = (1.+foso)*dp1 - foso*dp2 
          s(j,k,4) = dpress/(gm1*q(j,k,nq)) -
     >          0.5*(q(j,k,2)**2+q(j,k,3)**2)*s(j,k,1)/q(j,k,1)**2 +
     >          (q(j,k,2)*s(j,k,2) + q(j,k,3)*s(j,k,3))/q(j,k,1)
 11     continue

c..extrapolate everything to other halo cells

        do kc = 1,ke-ks
          if (idir.eq.2) then
            k = ke - kc
          elseif (idir.eq.-2) then
            k = ks + kc
          endif
          k1 = k  + iadd
          k2 = k1 + iadd

          foso = 1.0

          do j = js,je
            drho1 =  s(j,k1,1)*q(j,k1,nq)
            drho2 =  s(j,k2,1)*q(j,k2,nq)

            u1 = q(j,k1,2)/q(j,k1,1)  
            u2 = q(j,k2,2)/q(j,k2,1)  
            du1 = s(j,k1,2)/q(j,k1,1) - q(j,k1,2)*s(j,k1,1)/q(j,k1,1)**2
            du2 = s(j,k2,2)/q(j,k2,1) - q(j,k2,2)*s(j,k2,1)/q(j,k2,1)**2

            v1 = q(j,k1,3)/q(j,k1,1)  
            v2 = q(j,k2,3)/q(j,k2,1)  
            dv1 = s(j,k1,3)/q(j,k1,1) - q(j,k1,3)*s(j,k1,1)/q(j,k1,1)**2
            dv2 = s(j,k2,3)/q(j,k2,1) - q(j,k2,3)*s(j,k2,1)/q(j,k2,1)**2

            dp1 = gm1*(s(j,k1,4) + 
     >            0.5*(q(j,k1,2)**2+q(j,k1,3)**2)*s(j,k1,1)/q(j,k1,1)**2 -
     >            (q(j,k1,2)*s(j,k1,2) + q(j,k1,3)*s(j,k1,3))/q(j,k1,1))*
     >            q(j,k1,nq)
            dp2 = gm1*(s(j,k2,4) + 
     >            0.5*(q(j,k2,2)**2+q(j,k2,3)**2)*s(j,k2,1)/q(j,k2,1)**2 -
     >            (q(j,k2,2)*s(j,k2,2) + q(j,k2,3)*s(j,k2,3))/q(j,k2,1))*
     >            q(j,k2,nq)

            s(j,k,1) = ((1.+foso)*drho1 - foso*drho2 )/q(j,k,nq)
            s(j,k,2) = ((1.+foso)*du1 - foso*du2)*q(j,k,1) + 
     >                 ((1.+foso)*u1 - foso*u2)*s(j,k,1)
            s(j,k,3) = ((1.+foso)*dv1 - foso*dv2)*q(j,k,1) + 
     >                 ((1.+foso)*v1 - foso*v2)*s(j,k,1)
            dpress = (1.+foso)*dp1 - foso*dp2 
            s(j,k,4) = dpress/(gm1*q(j,k,nq)) -
     >            0.5*(q(j,k,2)**2+q(j,k,3)**2)*s(j,k,1)/q(j,k,1)**2 +
     >            (q(j,k,2)*s(j,k,2) + q(j,k,3)*s(j,k,3))/q(j,k,1)
          enddo
        enddo

c..slowly turn on the wall boundary condition 
c  over 30 time steps
        t = float(istep0)/30.
        if(t.gt.1..or.iread.gt.0) t=1.
        scal = (10.-15.*t+6.*t*t)*t**3
        rscal = 1.-scal
        do k=ks,ke
          do j=js,je
            ajacinv = 1./q(j,k,nq)
            s(j,k,1) = s(j,k,1)*scal
            s(j,k,2) = s(j,k,2)*scal
            s(j,k,3) = s(j,k,3)*scal
            s(j,k,4) = s(j,k,4)*scal
          enddo
        enddo
      
      endif

c...resetting the inviscid to original value
      !am invisc = tmp

      return
      end

c***********************************************************************
      subroutine bcextp(q,jd,kd,js,je,ks,ke,idir)
c
c..symmetric bc
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,j2,k2,jc,kc
      integer iadd,iadir
      real foso,scale1,scale2

      foso = 0.0
c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      if (iadir.eq.1) then
        do jc = 1,je-js+1
          if (idir.eq.1) then
            j = je - jc + 1
          elseif (idir.eq.-1) then
            j = js + jc - 1
          endif
          j1 = j + iadd
          j2 = j1 + iadd
          do k = ks,ke
            scale1 = q(j1,k,nq)/q(j,k,nq)
            scale2 = q(j2,k,nq)/q(j,k,nq)
            q(j,k,1) = (1.+foso)*q(j1,k,1)*scale1 - foso*q(j2,k,1)*scale2
            q(j,k,2) = (1.+foso)*q(j1,k,2)*scale1 - foso*q(j2,k,2)*scale2
            q(j,k,3) = (1.+foso)*q(j1,k,3)*scale1 - foso*q(j2,k,3)*scale2
            q(j,k,4) = (1.+foso)*q(j1,k,4)*scale1 - foso*q(j2,k,4)*scale2
          enddo
        enddo
      elseif (iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
          elseif (idir.eq.-2) then
            k = ks + kc - 1
          endif
          k1 = k + iadd
          k2 = k1 + iadd
          do j = js,je
            scale1 = q(j,k1,nq)/q(j,k,nq)
            scale2 = q(j,k2,nq)/q(j,k,nq)
            q(j,k,1) = (1.+foso)*q(j,k1,1)*scale1 - foso*q(j,k2,1)*scale2
            q(j,k,2) = (1.+foso)*q(j,k1,2)*scale1 - foso*q(j,k2,2)*scale2
            q(j,k,3) = (1.+foso)*q(j,k1,3)*scale1 - foso*q(j,k2,3)*scale2
            q(j,k,4) = (1.+foso)*q(j,k1,4)*scale1 - foso*q(j,k2,4)*scale2
          enddo
        enddo
      endif
c
      return
      end

c***********************************************************************
      subroutine rbcextp(q,s,jd,kd,js,je,ks,ke,idir)
c
c..symmetric bc
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real q(jd,kd,nq),s(jd,kd,nv)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,j2,k2,jc,kc
      integer iadd,iadir
      real foso,scale1,scale2

      foso = 0.0
c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      if (iadir.eq.1) then
        do jc = 1,je-js+1
          if (idir.eq.1) then
            j = je - jc + 1
          elseif (idir.eq.-1) then
            j = js + jc - 1
          endif
          j1 = j + iadd
          j2 = j1 + iadd
          do k = ks,ke
            scale1 = q(j1,k,nq)/q(j,k,nq)
            scale2 = q(j2,k,nq)/q(j,k,nq)
            s(j,k,1) = (1.+foso)*s(j1,k,1)*scale1 - foso*s(j2,k,1)*scale2
            s(j,k,2) = (1.+foso)*s(j1,k,2)*scale1 - foso*s(j2,k,2)*scale2
            s(j,k,3) = (1.+foso)*s(j1,k,3)*scale1 - foso*s(j2,k,3)*scale2
            s(j,k,4) = (1.+foso)*s(j1,k,4)*scale1 - foso*s(j2,k,4)*scale2
          enddo
        enddo
      elseif (iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
          elseif (idir.eq.-2) then
            k = ks + kc - 1
          endif
          k1 = k + iadd
          k2 = k1 + iadd
          do j = js,je
            scale1 = q(j,k1,nq)/q(j,k,nq)
            scale2 = q(j,k2,nq)/q(j,k,nq)
            s(j,k,1) = (1.+foso)*s(j,k1,1)*scale1 - foso*s(j,k2,1)*scale2
            s(j,k,2) = (1.+foso)*s(j,k1,2)*scale1 - foso*s(j,k2,2)*scale2
            s(j,k,3) = (1.+foso)*s(j,k1,3)*scale1 - foso*s(j,k2,3)*scale2
            s(j,k,4) = (1.+foso)*s(j,k1,4)*scale1 - foso*s(j,k2,4)*scale2
          enddo
        enddo
      endif
c
      return
      end

c***********************************************************************
      subroutine bcsym(q,jd,kd,js,je,ks,ke,idir)
c
c..symmetric bc
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,jc,kc
      integer iadd,iadir
      real scale

c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      if (iadir.eq.1) then
        do jc = 1,je-js+1
          if (idir.eq.1) then
            j = je - jc + 1
            j1 = je + jc
          elseif (idir.eq.-1) then
            j = js + jc - 1
            j1 = js - jc
          endif
          do k = ks,ke
            scale = q(j1,k,nq)/q(j,k,nq)
            q(j,k,1) = q(j1,k,1)*scale
            q(j,k,2) = -q(j1,k,2)*scale
            q(j,k,3) = q(j1,k,3)*scale
            q(j,k,4) = q(j1,k,4)*scale
          enddo
        enddo
      elseif (iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
            k1 = ke + kc
          elseif (idir.eq.-2) then
            k = ks + kc - 1
            k1 = ks - kc
          endif
          do j = js,je
            scale = q(j,k1,nq)/q(j,k,nq)
            q(j,k,1) = q(j,k1,1)*scale
            q(j,k,2) = q(j,k1,2)*scale
            q(j,k,3) = -q(j,k1,3)*scale
            q(j,k,4) = q(j,k1,4)*scale
          enddo
        enddo
      endif
c
      return
      end

c***********************************************************************
      subroutine rbcsym(q,s,jd,kd,js,je,ks,ke,idir)
c
c..symmetric bc
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real q(jd,kd,nq),s(jd,kd,nv)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,jc,kc
      integer iadd,iadir
      real scale

c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      if (iadir.eq.1) then
        do jc = 1,je-js+1
          if (idir.eq.1) then
            j = je - jc + 1
            j1 = je + jc
          elseif (idir.eq.-1) then
            j = js + jc - 1
            j1 = js - jc
          endif
          do k = ks,ke
            scale = q(j1,k,nq)/q(j,k,nq)
            s(j,k,1) = s(j1,k,1)*scale
            s(j,k,2) = -s(j1,k,2)*scale
            s(j,k,3) = s(j1,k,3)*scale
            s(j,k,4) = s(j1,k,4)*scale
          enddo
        enddo
      elseif (iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
            k1 = ke + kc
          elseif (idir.eq.-2) then
            k = ks + kc - 1
            k1 = ks - kc
          endif
          do j = js,je
            scale = q(j,k1,nq)/q(j,k,nq)
            s(j,k,1) = s(j,k1,1)*scale
            s(j,k,2) = s(j,k1,2)*scale
            s(j,k,3) = -s(j,k1,3)*scale
            s(j,k,4) = s(j,k1,4)*scale
          enddo
        enddo
      endif
c
      return
      end

c***********************************************************************
      subroutine uv(js,je,k,xx,xy,yx,yy,u,v,uwall,vwall,jd,kd)
c
c  solve for the cartesian momemtum components (ru,rv) from
c  the contravariant velocity components (u,v) along lines of j.
c     
c  note: consistent with metfv     
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real u(mdim), v(mdim)
      real uwall(mdim),vwall(mdim)
      ! local variables
      integer js,je,k,j
      real t11,t12,t21,t22,vol

c***  first executable statement

      do 101 j = js,je
c
        t11   = (yy(j,k))
        t12   =-(xy(j,k))
c
        t21   =-(yx(j,k))
        t22   = (xx(j,k))
c
        vol=1./(t22*t11-t12*t21)
c
c..the following are the physical plane velocities u & v
c
        uwall(j) = (t11*u(j) +t12*v(j))*vol
        vwall(j) = (t21*u(j) +t22*v(j))*vol

  101 continue

      return
      end
      
c***********************************************************************
      subroutine ruv(js,je,k,xx,xy,yx,yy,du,dv,duwall,dvwall,jd,kd)
c
c  solve for the cartesian momemtum components (ru,rv) from
c  the contravariant velocity components (u,v) along lines of j.
c     
c  note: consistent with metfv     
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real du(mdim), dv(mdim)
      real duwall(mdim),dvwall(mdim)
      ! local variables
      integer js,je,k,j
      real t11,t12,t21,t22,vol

c***  first executable statement

      do 101 j = js,je
c
        t11   = (yy(j,k))
        t12   =-(xy(j,k))
c
        t21   =-(yx(j,k))
        t22   = (xx(j,k))
c
        vol=1./(t22*t11-t12*t21)
c
c..the following are the physical plane velocities u & v
c
        duwall(j) = (t11*du(j) +t12*dv(j))*vol
        dvwall(j) = (t21*du(j) +t22*dv(j))*vol

  101 continue

      return
      end
      
c***********************************************************************
