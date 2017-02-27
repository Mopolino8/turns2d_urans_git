c***********************************************************************
      subroutine komegasst(q,s,turmu,x,y,xv,yv,xx,xy,yx,yy,ug,vg,jd,kd,
     >		tscale,iblank,im)
c  turbulent eddy viscosity. model is one equation spalart-
c  allmaras. ref{ aiaa 92-0439}.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer im,jd,kd
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd), tscale(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd) 
      integer iblank(jd,kd)

      ! local variables
      real,allocatable :: bjmat(:,:,:),sjmat(:,:,:)
      real,allocatable :: aj(:,:,:),cj(:,:,:)
      real,allocatable :: ak(:,:,:),ck(:,:,:)
      real,allocatable :: rho(:,:),u(:,:),v(:,:)
      real,allocatable :: vmul(:,:),tauxx(:,:),tauxy(:,:),tauyy(:,:)
      real,allocatable :: ux(:,:),uy(:,:),vx(:,:),vy(:,:)
      real,allocatable :: cdkw(:,:),F1(:,:),F2(:,:),sn(:,:)
      real,allocatable :: vort(:,:),strain(:,:)
      
      integer k,j,n,jloc,kloc
      real relfac,restke,restomega,resmax,oat,tscal,dtpseudo_turb
      real tkelim,tomegalim

      allocate(bjmat(jd,kd,2),sjmat(jd,kd,2))
      allocate(aj(jd,kd,2),cj(jd,kd,2))
      allocate(ak(jd,kd,2),ck(jd,kd,2))
      allocate(rho(jd,kd),u(jd,kd),v(jd,kd))
      allocate(vmul(jd,kd),tauxx(jd,kd),tauxy(jd,kd),tauyy(jd,kd))
      allocate(ux(jd,kd),uy(jd,kd),vx(jd,kd),vy(jd,kd))
      allocate(cdkw(jd,kd),F1(jd,kd),F2(jd,kd),sn(jd,kd))
      allocate(vort(jd,kd),strain(jd,kd))

c***  first executable statement

      tkelim = 1e-20
      tomegalim = 1e-20
c...laminar co-efficient of viscosity calculation
      
      call lamvis(q,vmul,jd,kd)

c...compute velocities, turbulent stress and velocity gradients

      do  k = 1,kd
      do  j = 1,jd
        rho(j,k) = q(j,k,1)*q(j,k,nq)
        u(j,k)  = q(j,k,2)/q(j,k,1)
        v(j,k)  = q(j,k,3)/q(j,k,1)
        sn(j,k) = 1.e10
      enddo
      enddo

c...apply boundary condition

      call komegabc(q,vmul,rho,x,y,xv,yv,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)

c...compute distance to wall

      if(bodyflag(im)) call kwdist(sn,x,y,xv,yv,jd,kd)

c...compute vorticity, strain, turbulent stress and cross-diffusion

      call kwvortic(vort,strain,u,v,xx,xy,yx,yy,jd,kd,jbeg,jend,kbeg,kend)
      call turbstress(q,turmu,tauxx,tauxy,tauyy,u,v,ux,uy,vx,vy,
     &                xx,xy,yx,yy,jd,kd,jbeg,jend,kbeg,kend)
      call crossdiff(q,rho,cdkw,xx,xy,yx,yy,jd,kd,jbeg,jend,kbeg,kend)

c...compute blending function F1
 
      call blend(q,rho,vmul,cdkw,sn,F1,F2,jd,kd,jbeg,jend,kbeg,kend)

c...compute eddy viscosity

      call findturm(q,rho,strain,F2,turmu,jd,kd)

c...compute intermittancy factor, if transition model used

      if (itrans.eq.1) then
        call sst_gammatheta(q,s,turmu,sn,vort,xx,xy,yx,yy,ug,vg,
     >                        jd,kd,tscale,iblank,im)
      endif

c...compute rhs and lhs

      call komegarhslhs(q,s,vmul,turmu,tauxx,tauxy,tauyy,cdkw,ux,uy,vx,vy,
     &                  F1,rho,u,v,ug,vg,xx,xy,yx,yy,aj,bjmat,cj,ak,ck,
     &                  sjmat,jd,kd,jbeg,jend,kbeg,kend)

c...invert using DDADI

c..set time-accuracy
      oat = 1.0
      if (ntac.eq.-2) oat = 0.5
      if (ntac.ge.2 .and. istep.gt.1) oat = 2./3.
      if (ntac.eq.3 .and. istep.gt.2) oat = 6./11.

      restke=0.0
      restomega=0.0
      do k = kbeg,kend
      do j = jbeg,jend

        tscal = tscale(j,k)
        if(timeac.eq.1) then
          dtpseudo_turb=10.0
          tscal = max(iblank(j,k),0)*( 1.0 + 0.002*sqrt(q(j,k,nq)))
     <                       /(1.+sqrt(q(j,k,nq)))
          tscal = tscal*dtpseudo_turb
          tscal = tscal/(1.+tscal/h/oat)
        endif

        do n = 1,2
          aj(j,k,n) = aj(j,k,n)*tscal
          cj(j,k,n) = cj(j,k,n)*tscal
          ak(j,k,n) = ak(j,k,n)*tscal
          ck(j,k,n) = ck(j,k,n)*tscal
          bjmat(j,k,n) = 1.+ bjmat(j,k,n)*tscal
          sjmat(j,k,n) = sjmat(j,k,n)*tscal
          !write(12,'(i4,3f12.5)'),n,aj(j,k,n),bjmat(j,k,n),cj(j,k,n)
          !write(13,'(i4,3f12.5)'),n,ak(j,k,n),bjmat(j,k,n),ck(j,k,n)
        enddo
      enddo
      enddo

      call lsolvej2(aj,bjmat,cj,sjmat,jd,kd,jbeg,jend,kbeg,kend)

      do k = kbeg,kend
      do j = jbeg,jend
        sjmat(j,k,1) = sjmat(j,k,1)*bjmat(j,k,1)
        sjmat(j,k,2) = sjmat(j,k,2)*bjmat(j,k,2)
      enddo
      enddo

      call lsolvek2(ak,bjmat,ck,sjmat,jd,kd,jbeg,jend,kbeg,kend)

      relfac = 1.
      resmax = 0.
      do k = kbeg,kend
      do j = jbeg,jend
        s(j,k,nmv+1) = relfac*sjmat(j,k,1)
        s(j,k,nmv+2) = relfac*sjmat(j,k,2)

        q(j,k,nmv+1) = q(j,k,nmv+1) + s(j,k,nmv+1)*max(iblank(j,k),0)
        q(j,k,nmv+2) = q(j,k,nmv+2) + s(j,k,nmv+2)*max(iblank(j,k),0)

        q(j,k,nmv+1) = max(q(j,k,nmv+1),tkelim*rho(j,k))
        q(j,k,nmv+2) = max(q(j,k,nmv+2),tomegalim*rho(j,k))

        restke = restke + s(j,k,nmv+1)**2
        restomega = restomega + s(j,k,nmv+2)**2
        if (max(restke,restomega).gt.resmax) then
          jloc = j
          kloc = k
          resmax = max(restke,restomega)
        endif
      enddo
      enddo

      restke = sqrt(restke/jmax/kmax)
      restomega = sqrt(restomega/jmax/kmax)
      resmax = sqrt(resmax/jmax/kmax)
!$OMP ORDERED
      if( mod(istep,npnorm).eq.0) then
        if (.not.timespectral.or.itn.eq.1) then
           write(1233+im,*),istep0,restke,restomega,resmax
        endif
      endif
!$OMP END ORDERED
c
c...apply boundary condition again

      call komegabc(q,vmul,rho,x,y,xv,yv,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)

c...compute eddy viscosity

      call findturm(q,rho,strain,F2,turmu,jd,kd)

      return
      end

c*************************************************************
      subroutine findturm(q,rho,strain,F2,turmu,jd,kd)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq),rho(jd,kd),strain(jd,kd),F2(jd,kd),turmu(jd,kd)
      ! local variables
      integer j,k
      real a1,rei
 
      data a1 /0.31/
      
      rei = 1./rey
 
      do k   = 1,kd
      do j   = 1,jd
        turmu(j,k) = a1*rho(j,k)*q(j,k,nmv+1)/
     &                max(a1*q(j,k,nmv+2),F2(j,k)*rho(j,k)*strain(j,k)*rei)
        turmu(j,k) = max(turmu(j,k),1e-20)
      enddo
      enddo

      return
      end

c*************************************************************
        subroutine kwdist(sn,x,y,xv,yv,jd,kd)
c
c     calculate distance to solid surface
c     not rigorous now
c     
c*************************************************************
        use params_global
c*************************************************************
        implicit none
c*************************************************************
        integer jd,kd
        real sn(jd,kd)
        real x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax)

        ! local variables
        integer j,k,jv

c***  first executable statement
       
        k=1
        do j=jtail1,jtail2
        sn(j,k)=0.0
        enddo
        
        do j=jtail1,jtail2
        do k=2,kd
          jv = j - nhalo
          sn(j,k)=sqrt((x(j,k)-0.5*(xv(jv,1)+xv(jv+1,1)))**2 +
     &                 (y(j,k)-0.5*(yv(jv,1)+yv(jv+1,1)))**2)
        enddo
        enddo

        do j=1,jtail1-1
        do k=1,kd
          jv = jtail1 - nhalo
          sn(j,k)=sqrt((x(j,k)-xv(jv,1))**2 + (y(j,k)-yv(jv,1))**2)
!          sn(j,k)=1.e10
        enddo
        enddo

        do j=jtail2+1,jd
        do k=1,kd
          jv = jtail1 - nhalo
          sn(j,k)=sqrt((x(j,k)-xv(jv,1))**2 + (y(j,k)-yv(jv,1))**2)
!          sn(j,k)=1.e10
        enddo
        enddo

        return
        end

c*************************************************************
        subroutine kwvortic(vort,strain,u,v,xx,xy,yx,yy,jd,kd,js,je,ks,ke)
c*************************************************************
        use params_global
c*************************************************************
        implicit none
c*************************************************************
        integer jd,kd,js,je,ks,ke
        real vort(jd,kd),strain(jd,kd)
        real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
        real u(jd,kd),v(jd,kd)

        ! local variables
        integer j,k,jp,jm1,kp,km1
        real ux,vx,uy,vy
        real usi,vsi,ueta,veta
        real sxx,sxy,syy

c**   first executable statement

        do j = js,je
          jp = j + 1
          jm1 = j - 1
          do k = ks,ke
            kp = k + 1
            km1 = k - 1

            usi  = 0.5*(u(jp,k)-u(jm1,k))
            vsi  = 0.5*(v(jp,k)-v(jm1,k))
            
            ueta = 0.5*(u(j,kp)-u(j,km1))
            veta = 0.5*(v(j,kp)-v(j,km1))
              
            ux = xx(j,k)*usi + yx(j,k)*ueta
            uy = xy(j,k)*usi + yy(j,k)*ueta

            vx = xx(j,k)*vsi + yx(j,k)*veta
            vy = xy(j,k)*vsi + yy(j,k)*veta

            sxx = ux
            sxy = 0.5*(uy + vx)
            syy = vy

            vort(j,k) = abs(uy - vx)
            strain(j,k) = sqrt(2.*(sxx*sxx + 2.*sxy*sxy + syy*syy))
          enddo
        enddo

        return
        end

c*************************************************************
      subroutine turbstress(q,turmu,tauxx,tauxy,tauyy,
     &                   u,v,ux,uy,vx,vy,xx,xy,yx,yy,jd,kd,js,je,ks,ke)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq),turmu(jd,kd)
      real tauxx(jd,kd),tauxy(jd,kd),tauyy(jd,kd)
      real u(jd,kd),v(jd,kd)
      real ux(jd,kd),uy(jd,kd),vx(jd,kd),vy(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)

! local variables
      integer j,k,jm1,km1,jp,kp
      real usi,vsi,ueta,veta,sfdiv,rei

      rei = 1./rey
      do j = js,je 
        jp = j + 1
        jm1 = j - 1
        do k = ks,ke
          kp = k + 1
          km1 = k - 1

          usi  = 0.5*(u(jp,k)-u(jm1,k))
          vsi  = 0.5*(v(jp,k)-v(jm1,k))
          
          ueta = 0.5*(u(j,kp)-u(j,km1))
          veta = 0.5*(v(j,kp)-v(j,km1))
            
          ux(j,k) = xx(j,k)*usi + yx(j,k)*ueta
          uy(j,k) = xy(j,k)*usi + yy(j,k)*ueta

          vx(j,k) = xx(j,k)*vsi + yx(j,k)*veta
          vy(j,k) = xy(j,k)*vsi + yy(j,k)*veta

          sfdiv = 2./3*(ux(j,k) + vy(j,k))

          tauxx(j,k) = turmu(j,k)* rei * ( 2.0 * ux(j,k) - sfdiv )
     &                            -2./3*q(j,k,nmv+1)
          tauyy(j,k) = turmu(j,k)* rei * ( 2.0 * vy(j,k) - sfdiv )
     &                            -2./3*q(j,k,nmv+1)
          tauxy(j,k) = turmu(j,k)* rei * ( uy(j,k) + vx(j,k) )

        enddo
      enddo

      return
      end

c***********************************************************************
      subroutine crossdiff(q,rho,cdkw,xx,xy,yx,yy,jd,kd,js,je,ks,ke)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq),rho(jd,kd),cdkw(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)

      ! local variables
      integer j,k,jp,kp,jm1,km1
      real ksi,keta,omegasi,omegaeta,kx,ky,omegax,omegay

      include 'komega.h'

      do j = js,je
        jp = j+1
        jm1 = j-1
        do k = ks,ke
          kp = k+1
          km1 = k-1
          ksi = 0.5*(q(jp,k,nmv+1)/rho(jp,k)-q(jm1,k,nmv+1)/rho(jm1,k))
          keta = 0.5*(q(j,kp,nmv+1)/rho(j,kp)-q(j,km1,nmv+1)/rho(j,km1))

          omegasi = 0.5*(q(jp,k,nmv+2)/rho(jp,k)-q(jm1,k,nmv+2)/rho(jm1,k))
          omegaeta = 0.5*(q(j,kp,nmv+2)/rho(j,kp)-q(j,km1,nmv+2)/rho(j,km1))

          kx = ksi*xx(j,k) + keta*yx(j,k)
          ky = ksi*xy(j,k) + keta*yy(j,k)
          
          omegax = omegasi*xx(j,k) + omegaeta*yx(j,k)
          omegay = omegasi*xy(j,k) + omegaeta*yy(j,k)

          cdkw(j,k) = kx*omegax + ky*omegay
          cdkw(j,k) = cdkw(j,k)*2*sigmaw2*rho(j,k)*rho(j,k)/q(j,k,nmv+2)
        enddo
      enddo

      return
      end

c***********************************************************************
      subroutine blend(q,rho,vmul,cdkw,sn,F1,F2,jd,kd,js,je,ks,ke)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq),rho(jd,kd),vmul(jd,kd),cdkw(jd,kd),sn(jd,kd)
      real F1(jd,kd),F2(jd,kd)

      ! local variables
      integer j,k
      real a1,a2,a3,rei,arg1,arg2,cd
      real Ry,F3

      include 'komega.h'

      rei = 1./rey
      do j = js,je
      do k = ks,ke
        a1 = rei*sqrt(q(j,k,nmv+1)*rho(j,k))/(0.09*q(j,k,nmv+2)*sn(j,k))
        a2 = 500*rei*rei*vmul(j,k)/(q(j,k,nmv+2)*sn(j,k)**2)
        cd = max(cdkw(j,k),1e-10)
        a3 = 4*q(j,k,nmv+1)*sigmaw2/(cd*sn(j,k)**2)
        arg1 = min(max(a1,a2),a3)
        arg2 = max(2*a1,a2)
        F1(j,k) = tanh(arg1**4) 
        F2(j,k) = tanh(arg2**2) 
      enddo
      enddo

c...Modified blending function, F1 when transition model is used

      if (itrans.eq.1) then
        do j = js,je
        do k = ks,ke
          Ry = sn(j,k)*sqrt(rho(j,k)*q(j,k,nmv+1))*rey/vmul(j,k)
          F3 = exp(-(Ry/120)**8)
          F1(j,k) = max(F1(j,k),F3)
        enddo
        enddo
      endif

      return
      end

c***********************************************************************
      subroutine komegarhslhs(q,s,vmul,turmu,tauxx,tauxy,tauyy,cdkw,
     &                        ux,uy,vx,vy,F1,rho,u,v,ug,vg,xx,xy,yx,yy,
     &                        aj,bjmat,cj,ak,ck,sjmat,jd,kd,js,je,ks,ke)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq), s(jd,kd,nv),rho(jd,kd)
      real vmul(jd,kd), turmu(jd,kd)
      real tauxx(jd,kd),tauxy(jd,kd),tauyy(jd,kd),cdkw(jd,kd)
      real ux(jd,kd),uy(jd,kd),vx(jd,kd),vy(jd,kd)
      real F1(jd,kd)
      real u(jd,kd),v(jd,kd),ug(jd,kd),vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real aj(jd,kd,2),bjmat(jd,kd,2),cj(jd,kd,2)
      real ak(jd,kd,2),ck(jd,kd,2),sjmat(jd,kd,2)

      ! local variables
      integer j,k,n,jp,kp,jm1,km1
      real uu,vv,fwd,bck
      real,allocatable :: up(:),um(:),vp(:),vm(:)

      allocate(up(mdim),um(mdim),vp(mdim),vm(mdim))

c***  first executable statement

      do k = 1,kd
      do j = 1,jd
        do n = 1,2
          aj(j,k,n)=0.0
          bjmat(j,k,n)=0.0
          cj(j,k,n)=0.0
          ak(j,k,n)=0.0
          ck(j,k,n)=0.0
          sjmat(j,k,n)=s(j,k,n+nmv)
        enddo
      enddo
      enddo

c...lhs contribution from the convection term

      do k=ks,ke
        do j=js-1,je+1
          uu=xx(j,k)*(u(j,k)-ug(j,k))+xy(j,k)*(v(j,k)-vg(j,k))
          up(j) = 0.5*(uu+abs(uu))
          um(j) = 0.5*(uu-abs(uu))
        enddo

        do j = js,je
          jp = j + 1
          jm1 = j - 1

          if(up(j).gt.1.0e-12) then
            fwd=1.0
          else
            fwd=0.0
          endif

          if(um(j).lt.-1.0e-12) then
            bck=1.0
          else
            bck=0.0
          endif

          sjmat(j,k,1) = sjmat(j,k,1) - up(j)*(q(j,k,nmv+1) - q(jm1,k,nmv+1)) - 
     &                          um(j)*(q(jp,k,nmv+1) - q(j,k,nmv+1))
          sjmat(j,k,2) = sjmat(j,k,2) - up(j)*(q(j,k,nmv+2) - q(jm1,k,nmv+2)) - 
     &                          um(j)*(q(jp,k,nmv+2) - q(j,k,nmv+2))
          do n = 1,2
            aj(j,k,n) = aj(j,k,n) - up(j)!fwd*(up(jm1)+um(jm1))
            cj(j,k,n) = cj(j,k,n) + um(j)!bck*(up(jp)+um(jp))
            bjmat(j,k,n) = bjmat(j,k,n) + up(j)- um(j)
          enddo
        enddo
      enddo

      do j=js,je
        do k=ks-1,ke+1
          vv = yx(j,k)*(u(j,k)-ug(j,k))+yy(j,k)*(v(j,k)-vg(j,k))
          vp(k) = 0.5*(vv+abs(vv))
          vm(k) = 0.5*(vv-abs(vv))
        enddo

        do k = ks,ke
          kp = k + 1
          km1 = k - 1

          if(vp(k).gt.1.0e-12) then
            fwd=1.0
          else
            fwd=0.0
          endif

          if(vm(k).lt.-1.0e-12) then
            bck=1.0
          else
            bck=0.0
          endif

          sjmat(j,k,1) = sjmat(j,k,1) - vp(k)*(q(j,k,nmv+1) - q(j,km1,nmv+1)) -
     &                          vm(k)*(q(j,kp,nmv+1) - q(j,k,nmv+1))
          sjmat(j,k,2) = sjmat(j,k,2) - vp(k)*(q(j,k,nmv+2) - q(j,km1,nmv+2)) -
     &                          vm(k)*(q(j,kp,nmv+2) - q(j,k,nmv+2))
          do n = 1,2
            ak(j,k,n) = ak(j,k,n) - vp(k) !fwd*(vp(km1)+vm(km1))
            ck(j,k,n) = ck(j,k,n) + vm(k) !bck*(vp(kp)+vm(kp))
            bjmat(j,k,n) = bjmat(j,k,n) + vp(k) - vm(k)
          enddo
        enddo
      enddo
   
c...rhs and lhs contribution from the source term

      call komegasource(q,rho,turmu,tauxx,tauxy,tauyy,cdkw,ux,uy,vx,vy,
     &                  F1,bjmat,sjmat,jd,kd,js,je,ks,ke)


c...rhs contribution from the diffusion term

      call komegadiffus(q,rho,vmul,turmu,F1,xx,xy,yx,yy,
     &                  aj,bjmat,cj,ak,ck,sjmat,jd,kd,js,je,ks,ke)

      return
      end

c***********************************************************************
      subroutine komegadiffus(q,rho,vmul,turmu,F1,xx,xy,yx,yy,
     & aj,bjmat,cj,ak,ck,sjmat,jd,kd,js,je,ks,ke)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq),rho(jd,kd)
      real vmul(jd,kd),turmu(jd,kd)
      real F1(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real aj(jd,kd,2),bjmat(jd,kd,2),cj(jd,kd,2)
      real ak(jd,kd,2),ck(jd,kd,2),sjmat(jd,kd,2)
      
    ! local variables
      integer j,k

      real,allocatable :: chp(:,:),dnuhp(:,:)
      real,allocatable :: dmxh(:,:),dmyh(:,:)

      real rei,vnulh,vnuh,dxp,dxm,dcp,dcm,ax,cx,dyp,dym,ay,cy
      real sigmak,sigmaw

      include 'komega.h'

      allocate(chp(jd,kd),dnuhp(jd,kd))
      allocate(dmxh(jd,kd),dmyh(jd,kd))

c***  first executable statement

      rei = 1./rey

c...compute j direction differences.

c.....compute half-point co-efficients

      do j=js-1,je
      do k=ks-1,ke+1
        dmxh(j,k) = 0.5*(xx(j,k)+xx(j+1,k))
        dmyh(j,k) = 0.5*(xy(j,k)+xy(j+1,k))
      enddo
      enddo

c.....compute j-direction fluxes for k-equation

      do j=js-1,je
      do k=ks-1,ke+1
        vnulh    = 0.5*(vmul(j,k)+vmul(j+1,k))
        vnuh     = 0.5*(turmu(j,k)+turmu(j+1,k))
        sigmak   = F1(j,k)*sigmak1 + (1.-F1(j,k))*sigmak2
        chp(j,k) = rei*(vnulh+sigmak*vnuh)       
        dnuhp(j,k) = q(j+1,k,nmv+1)/rho(j+1,k)-q(j,k,nmv+1)/rho(j,k)
      enddo
      enddo

      do j=js,je
      do k=ks-1,ke+1

        dxp=dmxh(j,k)*xx(j,k)+dmyh(j,k)*xy(j,k)
        dxm=dmxh(j-1,k)*xx(j,k)+dmyh(j-1,k)*xy(j,k)

c.....enforce positivity (as suggested by overflow)
          
        dcp    = dxp*(chp(j,k))
        dcm    = dxm*(chp(j-1,k))
        ax=0.0
        cx=0.0
        if(k.ne.ke+1.and.k.ne.ks-1) then
          ax       = max(dcm,0.0)
          cx       = max(dcp,0.0)
        endif

c.....compute fluxes.
        sjmat(j,k,1)=sjmat(j,k,1)-ax*dnuhp(j-1,k)+cx*dnuhp(j,k)

c.....jacobian terms

        aj(j,k,1) = aj(j,k,1) - ax/rho(j-1,k)
        cj(j,k,1) = cj(j,k,1) - cx/rho(j+1,k)
        bjmat(j,k,1) = bjmat(j,k,1)+ (ax + cx)/rho(j,k)

      enddo
      enddo

c.....boundary terms in jacobians

      do  k = 1,kd
        aj(je+1,k,1) = 0.0
        cj(js-1,k,1) = 0.0
      enddo

c.....compute j-direction fluxes for omega-equation

      do j=js-1,je
      do k=ks-1,ke+1
        vnulh     = 0.5*(vmul(j,k)+vmul(j+1,k))
        vnuh     = 0.5*(turmu(j,k)+turmu(j+1,k))
        sigmaw   = F1(j,k)*sigmaw1 + (1.-F1(j,k))*sigmaw2
        chp(j,k) = rei*(vnulh+sigmaw*vnuh)       
        dnuhp(j,k) =q(j+1,k,nmv+2)/rho(j+1,k)-q(j,k,nmv+2)/rho(j,k)
      enddo
      enddo

      do j=js,je
      do k=ks-1,ke+1

        dxp=dmxh(j,k)*xx(j,k)+dmyh(j,k)*xy(j,k)
        dxm=dmxh(j-1,k)*xx(j,k)+dmyh(j-1,k)*xy(j,k)

c.....enforce positivity (as suggested by overflow)
          
        dcp    = dxp*(chp(j,k))
        dcm    = dxm*(chp(j-1,k))
        ax=0.0
        cx=0.0

        if(k.ne.ke+1.and.k.ne.ks-1) then
          ax       = max(dcm,0.0)
          cx       = max(dcp,0.0)
        endif

c.....compute fluxes.
        sjmat(j,k,2)=sjmat(j,k,2)-ax*dnuhp(j-1,k)+cx*dnuhp(j,k)

c.....jacobian terms

        aj(j,k,2) = aj(j,k,2) - ax/rho(j-1,k)
        cj(j,k,2) = cj(j,k,2) - cx/rho(j+1,k)
        bjmat(j,k,2) = bjmat(j,k,2)+ (ax + cx)/rho(j,k)

      enddo
      enddo

c.....boundary terms in jacobians

      do  k = 1,kd
        aj(je+1,k,2) = 0.0
        cj(js-1,k,2) = 0.0
      enddo


c...compute k direction differences.

c.....compute half-point co-efficients

      do k=ks-1,ke
      do j=js-1,je+1
        dmxh(j,k) = 0.5*(yx(j,k)+yx(j,k+1))
        dmyh(j,k) = 0.5*(yy(j,k)+yy(j,k+1))
      enddo
      enddo

c.....compute k-direction fluxes for k-equation

      do k=ks-1,ke
      do j=js-1,je+1
        vnulh     = 0.5*(vmul(j,k)+vmul(j,k+1))
        vnuh     =  0.5*(turmu(j,k)+turmu(j,k+1))
        sigmak   = F1(j,k)*sigmak1 + (1.-F1(j,k))*sigmak2
        chp(j,k) =rei*(vnulh+sigmak*vnuh)       
        dnuhp(j,k) =q(j,k+1,nmv+1)/rho(j,k+1)-q(j,k,nmv+1)/rho(j,k)
      enddo
      enddo

      do k=ks,ke
      do j=js-1,je+1

        dyp=dmxh(j,k)*yx(j,k)+dmyh(j,k)*yy(j,k)
        dym=dmxh(j,k-1)*yx(j,k)+dmyh(j,k-1)*yy(j,k)

c.....enforce positivity (as suggested by overflow)
       
        dcp    = dyp*(chp(j,k))
        dcm    = dym*(chp(j,k-1))

        ay=0.0
        cy=0.0
        if(j.ne.js-1.and.j.ne.je+1) then
          ay       = max(dcm,0.0)
          cy       = max(dcp,0.0)
        endif

c.....compute fluxes.
        sjmat(j,k,1)=sjmat(j,k,1)-ay*dnuhp(j,k-1)+cy*dnuhp(j,k)

c.....jacobian terms

          ak(j,k,1) = ak(j,k,1) - ay/rho(j,k-1)
          ck(j,k,1) = ck(j,k,1) - cy/rho(j,k+1)
          bjmat(j,k,1) = bjmat(j,k,1)+ (ay + cy)/rho(j,k)

      enddo
      enddo

c.....boundary terms in jacobians

      do  j = 1,kd
        ak(j,ke+1,1) = 0.0
        ck(j,ks-1,1) = 0.0
      enddo

c.....compute k-direction fluxes for omega-equation

      do k=ks-1,ke
      do j=js-1,je+1
        vnulh     = 0.5*(vmul(j,k)+vmul(j,k+1))
        vnuh     =  0.5*(turmu(j,k)+turmu(j,k+1))
        sigmaw   = F1(j,k)*sigmaw1 + (1.-F1(j,k))*sigmaw2
        chp(j,k) =rei*(vnulh+sigmaw*vnuh)       
        dnuhp(j,k) =q(j,k+1,nmv+2)/rho(j,k+1)-q(j,k,nmv+2)/rho(j,k)
      enddo
      enddo

      do k=ks,ke
      do j=js-1,je+1

        dyp=dmxh(j,k)*yx(j,k)+dmyh(j,k)*yy(j,k)
        dym=dmxh(j,k-1)*yx(j,k)+dmyh(j,k-1)*yy(j,k)

c.....enforce positivity (as suggested by overflow)
       
        dcp    = dyp*(chp(j,k))
        dcm    = dym*(chp(j,k-1))

        ay=0.0
        cy=0.0
        if(j.ne.js-1.and.j.ne.je+1) then
          ay       = max(dcm,0.0)
          cy       = max(dcp,0.0)
        endif

c.....compute fluxes.
        sjmat(j,k,2)=sjmat(j,k,2)-ay*dnuhp(j,k-1)+cy*dnuhp(j,k)

c.....jacobian terms

        ak(j,k,2) = ak(j,k,2) - ay/rho(j,k-1)
        ck(j,k,2) = ck(j,k,2) - cy/rho(j,k+1)
        bjmat(j,k,2) = bjmat(j,k,2)+ (ay + cy)/rho(j,k)

      enddo
      enddo

c.....boundary terms in jacobians

      do  j = 1,kd
        ak(j,ke+1,2) = 0.0
        ck(j,ks-1,2) = 0.0
      enddo

      return
      end

c***********************************************************************
        subroutine komegasource(q,rho,turmu,tauxx,tauxy,tauyy,cdkw,ux,uy,vx,vy,
     &                          F1,bjmat,sjmat,jd,kd,js,je,ks,ke)
c***********************************************************************
        use params_global
c***********************************************************************
        implicit none
        integer jd,kd,js,je,ks,ke
        real q(jd,kd,nq),rho(jd,kd),turmu(jd,kd)
        real tauxx(jd,kd),tauxy(jd,kd),tauyy(jd,kd),cdkw(jd,kd)
        real ux(jd,kd),uy(jd,kd),vx(jd,kd),vy(jd,kd)
        real F1(jd,kd)
        real bjmat(jd,kd,2),sjmat(jd,kd,2)
        
        ! local variables
        integer j,k
        real rei,prod,kprod,kdes,omegades,omegacrossdiff,ksust,omegasust
        real prodlim,gammap,beta
        
        include 'komega.h'

c**   first executable statement

        rei = 1./rey

        do k=ks,ke
        do j=js,je
          prod = tauxx(j,k)*ux(j,k) + tauxy(j,k)*(uy(j,k) + vx(j,k)) + 
     &           tauyy(j,k)*vy(j,k)
          kdes = betastar*rey*q(j,k,nmv+2)*q(j,k,nmv+1)/rho(j,k)
          ksust = betastar*rey*tkeinf*tomegainf*rho(j,k)
          prodlim = 10.*kdes
          prod = min(prod,prodlim)
          kprod = prod

c...modification for transition model

          if (itrans.eq.1) then
            kprod = q(j,k,nmv+nturb+1)/rho(j,k)*prod
            kdes = min(max(q(j,k,nmv+nturb+1)/rho(j,k),0.1),1.0)*kdes
            ksust = min(max(q(j,k,nmv+nturb+1)/rho(j,k),0.1),1.0)*ksust
          endif

          beta = F1(j,k)*beta1 + (1-F1(j,k))*beta2 
          omegades = beta*rey*q(j,k,nmv+2)*q(j,k,nmv+2)/rho(j,k)
          omegasust = beta*rey*tomegainf**2*rho(j,k)
          omegacrossdiff = (1-F1(j,k))*cdkw(j,k)*rei
          sjmat(j,k,1) = sjmat(j,k,1) + kprod - kdes + ksust
          gammap = F1(j,k)*gamma1 + (1-F1(j,k))*gamma2 
          sjmat(j,k,2) = sjmat(j,k,2) + gammap*rho(j,k)/turmu(j,k)*prod   
     &                        - omegades + omegacrossdiff + omegasust
          bjmat(j,k,1) = bjmat(j,k,1) + kdes/q(j,k,nmv+1) 
          bjmat(j,k,2) = bjmat(j,k,2) + (2*omegades+abs(omegacrossdiff))/q(j,k,nmv+2) 
        enddo
        enddo

        return
        end

c*************************************************************
        subroutine komegabc(q,vmul,rho,x,y,xv,yv,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)
c*************************************************************
        use params_global
c*************************************************************
        implicit none
c*************************************************************
        integer jd,kd
        real q(jd,kd,nq),vmul(jd,kd),rho(jd,kd)
        real x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax)
        real u(jd,kd),v(jd,kd),ug(jd,kd),vg(jd,kd)
        real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
        integer im

c..   local variables
        integer js,je,ks,ke,idir
        integer ib

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

c.. inviscid wind tunnel wall bc
          if (ibtyp_all(ib,im).eq.3.or.ibtyp_all(ib,im).eq.4) then
            call komegabc_extpt(q,jd,kd,js,je,ks,ke,idir)

c.. wall bc at l = 1 (only interior portion of wall)
          elseif (ibtyp_all(ib,im).eq.5) then
            call komegabc_wall(q,vmul,x,y,xv,yv,jd,kd,js,je,ks,ke,idir)

c.. symmetric bc
          elseif (ibtyp_all(ib,im).eq.10.or.ibtyp_all(ib,im).eq.11) then
            call komegabc_sym(q,jd,kd,js,je,ks,ke,idir)

c.. periodic bc
          elseif (ibtyp_all(ib,im).eq.22) then
            call komegabc_periodic(q,jd,kd,js,je,ks,ke,idir)

c.. averaging bc for wake
          elseif (ibtyp_all(ib,im).eq.51) then
            call komegabc_wake(q,jd,kd,js,je,ks,ke,idir)

c.. freesream enforcing bc
          elseif (ibtyp_all(ib,im).eq.46) then
            call komegabc_inf(q,rho,jd,kd,js,je,ks,ke,idir)

c.. freesream characteristic bc
          elseif (ibtyp_all(ib,im).eq.47) then
            call komegabc_out(q,rho,u,v,ug,vg,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)

          endif
        enddo

        return
        end

c*************************************************************
      subroutine komegabc_wall(q,vmul,x,y,xv,yv,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq),vmul(jd,kd)
      real x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,j2,k2,jv,kv,jc,kc,iadd,iadir
      real rei,dely2,t,scal,qwall5,qwall6,vmulwall,omegawall,foso

      include 'komega.h'

      foso = 1.0
      iadd = sign(1,idir)
      iadir = abs(idir)

      t = (float(istep0)-1.)/30.
      if(t.gt.1..or.iread.gt.0) t=1.
      scal = (10.-15.*t+6.*t*t)*t**3

      rei = 1./rey

      if(iadir.eq.1) then
        if (idir.eq.1) then
          j = je
        elseif (idir.eq.-1) then
          j = js
        endif
        j1 = j + iadd
        jv = j1 - nhalo
        do k=ks,ke
          kv = k - nhalo
          dely2 = (x(j1,k)-0.5*(xv(jv,kv)+xv(jv,kv+1)))**2 +
     &            (y(j1,k)-0.5*(yv(jv,kv)+yv(jv,kv+1)))**2
          vmulwall = 0.5*(vmul(j1,k) + vmul(j,k))
          omegawall = 60.*vmulwall*rei*rei/(beta1*dely2)
          qwall5 = -q(j1,k,nmv+1) 
          qwall6 = 2*omegawall - 6.*vmul(j1,k)*rei*rei/(beta1*dely2)
          q(j,k,nmv+1) = scal*qwall5 + (1.-scal)*q(j,k,nmv+1)
          q(j,k,nmv+2) = scal*qwall6 + (1.-scal)*q(j,k,nmv+2)
        enddo

        do jc = 1,je-js
          if (idir.eq.1) then
            j = je - jc
          elseif (idir.eq.-1) then
            j = js + jc
          endif
          j1 = j  + iadd
          j2 = j1 + iadd
          do j = js,je
            q(j,k,nmv+1) = (1.+foso)*q(j1,k,nmv+1) - foso*q(j2,k,nmv+1) 
            q(j,k,nmv+2) = (1.+foso)*q(j1,k,nmv+2) - foso*q(j2,k,nmv+2) 
          enddo
        enddo

      elseif(iadir.eq.2) then
        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif
        k1 = k + iadd
        kv = k1 - nhalo
        do j=js,je
          jv = j - nhalo
          dely2 = (x(j,k1)-0.5*(xv(jv,kv)+xv(jv+1,kv)))**2 +
     &            (y(j,k1)-0.5*(yv(jv,kv)+yv(jv+1,kv)))**2
          vmulwall = 0.5*(vmul(j,k1) + vmul(j,k1))
          omegawall = 60.*vmulwall*rei*rei/(beta1*dely2)
          qwall5 = -q(j,k1,nmv+1) 
          qwall6 = 2*omegawall - 6.*vmul(j,k1)*rei*rei/(beta1*dely2)
          q(j,k,nmv+1) = scal*qwall5 + (1.-scal)*q(j,k,nmv+1)
          q(j,k,nmv+2) = scal*qwall6 + (1.-scal)*q(j,k,nmv+2)
        enddo

        do kc = 1,ke-ks
          if (idir.eq.2) then
            k = ke - kc
          elseif (idir.eq.-2) then
            k = ks + kc
          endif
          k1 = k  + iadd
          k2 = k1 + iadd
          do j = js,je
            q(j,k,nmv+1) = (1.+foso)*q(j,k1,nmv+1) - foso*q(j,k2,nmv+1) 
            q(j,k,nmv+2) = (1.+foso)*q(j,k1,nmv+2) - foso*q(j,k2,nmv+2) 
          enddo
        enddo

      endif

      return
      end

c*************************************************************
      subroutine komegabc_extpt(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,j2,k2,jc,kc,iadd,iadir
      real foso

      foso = 1.0
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        do jc = 1,je-js+1
          if (idir.eq.1) then
            j = je - jc + 1
          elseif (idir.eq.-1) then
            j = js + jc - 1
          endif
          j1 = j + iadd
          j2 = j1 + iadd
          do k=ks,ke
            q(j,k,nmv+1) = (1.+foso)*q(j1,k,nmv+1) - foso*q(j2,k,nmv+1)
            q(j,k,nmv+2) = (1.+foso)*q(j1,k,nmv+2) - foso*q(j2,k,nmv+2)
          enddo
        enddo
      elseif(iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
          elseif (idir.eq.-2) then
            k = ks + kc - 1
          endif
          k1 = k + iadd
          k2 = k1 + iadd
          do j=js,je
            q(j,k,nmv+1) = (1.+foso)*q(j,k1,nmv+1) - foso*q(j,k2,nmv+1)
            q(j,k,nmv+2) = (1.+foso)*q(j,k1,nmv+2) - foso*q(j,k2,nmv+2)
          enddo
        enddo
      endif

      return
      end

c*************************************************************
      subroutine komegabc_sym(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,jc,kc,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        do jc = 1,je-js+1
          if (idir.eq.1) then
            j = je - jc + 1
            j1 = je + jc
          elseif (idir.eq.-1) then
            j = js + jc - 1
            j1 = js - jc
          endif
          do k=ks,ke
            q(j,k,nmv+1) = q(j1,k,nmv+1)
            q(j,k,nmv+2) = q(j1,k,nmv+2)
          enddo
        enddo
      elseif(iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
            k1 = ke + kc
          elseif (idir.eq.-2) then
            k = ks + kc - 1
            k1 = ks - kc
          endif
          do j=js,je
            q(j,k,nmv+1) = q(j,k1,nmv+1)
            q(j,k,nmv+2) = q(j,k1,nmv+2)
          enddo
        enddo
      endif

      return
      end

c*************************************************************
      subroutine komegabc_periodic(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(idir.eq.1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = j - js + 1
           jc = jd - 2*jj + jj1
 
           do k = ks,ke
             q(j,k,nmv+1) = q(jc,k,nmv+1)
             q(j,k,nmv+2) = q(jc,k,nmv+2)
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j + 1
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             q(j,k,nmv+1) = q(jc,k,nmv+1)
             q(j,k,nmv+2) = q(jc,k,nmv+2)
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks + 1
           kc = kd - 2*kk + kk1
 
           do j = js,je
             q(j,k,nmv+1) = q(j,kc,nmv+1)
             q(j,k,nmv+2) = q(j,kc,nmv+2)
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k + 1
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             q(j,k,nmv+1) = q(j,kc,nmv+1)
             q(j,k,nmv+2) = q(j,kc,nmv+2)
           enddo
        enddo

      endif

      return
      end

c*************************************************************
      subroutine komegabc_wake(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jj,k1,kc,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in komegabc_wake'
      elseif(iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
            k1 = ke + kc
          elseif (idir.eq.-2) then
            k = ks + kc - 1
            k1 = ks - kc
          endif

          do j=js,je
            jj = jd - j + 1
            q(j,k,nmv+1)  = q(jj,k1,nmv+1)
            q(j,k,nmv+2)  = q(jj,k1,nmv+2)
            q(jj,k,nmv+1) = q(j,k1,nmv+1)
            q(jj,k,nmv+2) = q(j,k1,nmv+2)
          enddo
        enddo

      endif

      return
      end

c*************************************************************
      subroutine komegabc_inf(q,rho,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq),rho(jd,kd)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        do j = js,je
          do k = ks,ke
            q(j,k,nmv+1)=tkeinf*rho(j,k)
            q(j,k,nmv+2)=tomegainf*rho(j,k)
          enddo
        enddo
      elseif(iadir.eq.2) then
        do k = ks,ke
          do j=js,je
            q(j,k,nmv+1)=tkeinf*rho(j,k)
            q(j,k,nmv+2)=tomegainf*rho(j,k)
          enddo
        enddo
      endif

      return
      end

c*************************************************************
      subroutine komegabc_out(q,rho,u,v,ug,vg,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq),rho(jd,kd),u(jd,kd),v(jd,kd)
      real ug(jd,kd),vg(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,jc,kc,iadd,iadir
      real uu

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        if (idir.eq.1) then
          j = je
        elseif (idir.eq.-1) then
          j = js
        endif
        j1 = j  + iadd

        do k = ks,ke
          q(j,k,nmv+1)=tkeinf*rho(j,k)
          q(j,k,nmv+2)=tomegainf*rho(j,k)
          uu = (u(j,k)-ug(j,k)+u(j1,k)-ug(j1,k))*xx(j,k)
          uu = uu+(v(j,k)-vg(j,k)+v(j1,k)-vg(j1,k))*xy(j,k)
          uu = uu*iadd
          if(uu.lt.0.) then
            q(j,k,nmv+1) = q(j1,k,nmv+1)
            q(j,k,nmv+2) = q(j1,k,nmv+2)
          endif
        enddo

c..extrapolate everything to other halo cells

        do jc = 1,je-js
          if (idir.eq.1) then
            j = je - jc
          elseif (idir.eq.-1) then
            j = js + jc
          endif
          j1 = j  + iadd

          do k = ks,ke
            q(j,k,nmv+1) = q(j1,k,nmv+1)
            q(j,k,nmv+2) = q(j1,k,nmv+2)
          enddo

        enddo

      elseif(iadir.eq.2) then
        if (idir.eq.2) then
          k = ke
        elseif (idir.eq.-2) then
          k = ks
        endif
        k1 = k  + iadd

        do j=js,je
          q(j,k,nmv+1)=tkeinf*rho(j,k)
          q(j,k,nmv+2)=tomegainf*rho(j,k)
          uu=   (u(j,k)-ug(j,k)+u(j,k1)-ug(j,k1))*yx(j,k)
          uu=uu+(v(j,k)-vg(j,k)+v(j,k1)-vg(j,k1))*yy(j,k)
          uu=uu*iadd
          if(uu.lt.0.) then
            q(j,k,nmv+1)=q(j,k1,nmv+1)
            q(j,k,nmv+2)=q(j,k1,nmv+2)
          endif
        enddo

c..extrapolate everything to other halo cells

        do kc = 1,ke-ks
          if (idir.eq.2) then
            k = ke - kc
          elseif (idir.eq.-2) then
            k = ks + kc
          endif
          k1 = k  + iadd

          do j = js,je
            q(j,k,nmv+1) = q(j,k1,nmv+1)
            q(j,k,nmv+2) = q(j,k1,nmv+2)
          enddo

        enddo

      endif

      return
      end

c*************************************************************
        subroutine lsolvej2(a,b,c,s,jd,kd,js,je,ks,ke)
c*************************************************************
        use params_global
c*************************************************************
        
        integer jd,kd,js,je,ks,ke

        real a(jd,kd,2),b(jd,kd,2),c(jd,kd,2)
        real s(jd,kd,2)

        ! local variables

        integer j,k,n
        real bb

c***  first executable statement
 
        do n = 1,2
          do k = ks,ke
            bb       = 1./b(js,k,n)
            s(js,k,n)  = s(js,k,n)*bb
            c(js,k,n)  = c(js,k,n)*bb
          enddo
 
          do  j = js+1,je
          do  k = ks,ke
            bb      = 1./(b(j,k,n) - a(j,k,n)*c(j-1,k,n))
            s(j,k,n)  = (s(j,k,n) - a(j,k,n)*s(j-1,k,n))*bb
            c(j,k,n)  = c(j,k,n)*bb
          enddo
          enddo

          do j = je-1,js,-1
          do k = ks,ke
            s(j,k,n)    = s(j,k,n) - c(j,k,n)*s(j+1,k,n)
          enddo
          enddo
        enddo

        return
        end

c*************************************************************
        subroutine lsolvek2(a,b,c,s,jd,kd,js,je,ks,ke)
c*************************************************************
        use params_global
c*************************************************************
        implicit none
c*************************************************************
        integer jd,kd,js,je,ks,ke
        real a(jd,kd,2),b(jd,kd,2),c(jd,kd,2)
        real s(jd,kd,2)

        ! local variables

        integer j,k,n
        real bb

c***  first executable statement

        do n = 1,2
          do j = js,je
            bb       = 1./b(j,ks,n)
            s(j,ks,n)  = s(j,ks,n)*bb
            c(j,ks,n)  = c(j,ks,n)*bb
          enddo
 
          do  k = ks+1,ke
          do  j = js,je
            bb      = 1./(b(j,k,n) - a(j,k,n)*c(j,k-1,n))
            s(j,k,n)  = (s(j,k,n) - a(j,k,n)*s(j,k-1,n))*bb
            c(j,k,n)  = c(j,k,n)*bb
          enddo
          enddo

          do k = ke-1,ks,-1
          do j = js,je
            s(j,k,n)    = s(j,k,n) - c(j,k,n)*s(j,k+1,n)
          enddo
          enddo
        enddo

        return
        end

c*************************************************************
