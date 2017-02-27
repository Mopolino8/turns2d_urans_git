c***********************************************************************
      subroutine sst_gammatheta(q,s,turmu,sn,vort,xx,xy,yx,yy,ug,vg,
     >                          jd,kd,tscale,iblank,im)
c  correlation-based transition model - Langtry and Menter. 
c  Ref{ AIAA Journal, Vol. 47, No. 12}.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer im,jd,kd
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd)
      real sn(jd,kd), vort(jd,kd), tscale(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd) 
      integer iblank(jd,kd)

      ! local variables
      real,allocatable :: itmc(:,:),ret(:,:)
      real,allocatable :: bjmat(:,:,:),sjmat(:,:,:)
      real,allocatable :: aj(:,:,:),cj(:,:,:)
      real,allocatable :: ak(:,:,:),ck(:,:,:)
      real,allocatable :: itmcsep(:,:),rho(:,:),u(:,:),v(:,:),vmag(:,:)
      real,allocatable :: vmul(:,:),strain(:,:),du_ds(:,:)
      
      integer k,j,n,jloc,kloc
      real relfac,resitmc,resret,resmax,oat,tscal,dtpseudo_turb
      real itmclim,retlim

      allocate(itmc(jd,kd),ret(jd,kd))
      allocate(bjmat(jd,kd,2),sjmat(jd,kd,2))
      allocate(aj(jd,kd,2),cj(jd,kd,2))
      allocate(ak(jd,kd,2),ck(jd,kd,2))
      allocate(itmcsep(jd,kd),rho(jd,kd),u(jd,kd),v(jd,kd),vmag(jd,kd))
      allocate(vmul(jd,kd),strain(jd,kd),du_ds(jd,kd))

c***  first executable statement

      itmclim = 1e-20
      retlim = 1e-20
c...laminar co-efficient of viscosity calculation
      
      call lamvis(q,vmul,jd,kd)

      do  k = 1,kmax
      do  j = 1,jmax
        rho(j,k) = q(j,k,1)*q(j,k,nq)
        u(j,k)  = q(j,k,2)/q(j,k,1)
        v(j,k)  = q(j,k,3)/q(j,k,1)
        vmag(j,k) = sqrt(u(j,k)*u(j,k) + v(j,k)*v(j,k))
      enddo
      enddo

c...apply boundary condition

      call gammathetabc(q,rho,u,v,ug,vg,xx,xy,yx,yy,im)

c...initiate local working variable for intermittancy and Re-theta

      do  k = 1,kmax
      do  j = 1,jmax
        itmc(j,k) = q(j,k,7)
        ret(j,k) = q(j,k,8)
      enddo
      enddo

c...compute strain rate

      call calc_strain_du_ds(strain,du_ds,u,v,vmag,xx,xy,yx,yy,jd,kd)

c...compute rhs and lhs

      call gammathetarhslhs(q,s,itmc,ret,itmcsep,vmul,turmu,vort,strain,
     &                      du_ds,rho,u,v,vmag,ug,vg,xx,xy,yx,yy,
     &                      sn,aj,bjmat,cj,ak,ck,sjmat,jd,kd)

c...invert using DDADI

c..set time-accuracy
      oat = 1.0
      if (ntac.eq.-2) oat = 0.5
      if (ntac.ge.2 .and. istep.gt.1) oat = 2./3.
      if (ntac.eq.3 .and. istep.gt.2) oat = 6./11.

      resitmc=0.0
      resret=0.0
      do k = 1,kmax
      do j = 1,jmax

        tscal = tscale(j,k)
        if(timeac.eq.1.and.1.eq.0) then
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
          sjmat(j,k,n)  = sjmat(j,k,n)*tscal
          !write(12,'(i4,3f12.5)'),n,aj(j,k,n),bjmat(j,k,n),cj(j,k,n)
          !write(13,'(i4,3f12.5)'),n,ak(j,k,n),bjmat(j,k,n),ck(j,k,n)
        enddo
      enddo
      enddo

      call lsolvej2(aj,bjmat,cj,sjmat,jd,kd)

      do k = 1,kmax
      do j = 1,jmax
        sjmat(j,k,1) = sjmat(j,k,1)*bjmat(j,k,1)
        sjmat(j,k,2) = sjmat(j,k,2)*bjmat(j,k,2)
      enddo
      enddo

      call lsolvek2(ak,bjmat,ck,sjmat,jd,kd)

      relfac = 1.
      resmax = 0.
      do k = 1,kmax
      do j = 1,jmax
        sjmat(j,k,1) = relfac*sjmat(j,k,1)
        sjmat(j,k,2) = relfac*sjmat(j,k,2)

        itmc(j,k) = itmc(j,k) + sjmat(j,k,1)*max(iblank(j,k),0)
        ret(j,k) = ret(j,k) + sjmat(j,k,2)*max(iblank(j,k),0)

        itmc(j,k) = max(itmc(j,k),itmclim*rho(j,k))
        ret(j,k) = max(ret(j,k),retlim*rho(j,k))

c....limit intermittancy to one
        itmc(j,k) = min(itmc(j,k),rho(j,k))

        resitmc = resitmc + sjmat(j,k,1)**2
        resret = resret + sjmat(j,k,2)**2
        if (max(resitmc,resret).gt.resmax) then
          jloc = j
          kloc = k
          resmax = max(resitmc,resret)
        endif

c...Modification in intermittancy for separation-induced transition
        itmc(j,k) = max(itmc(j,k),rho(j,k)*itmcsep(j,k))

      enddo
      enddo

      resitmc = sqrt(resitmc/jmax/kmax)
      resret = sqrt(resret/jmax/kmax)
      resmax = sqrt(resmax/jmax/kmax)
      if( mod(istep,npnorm).eq.0) then
         write(1333+im,*),istep0,resitmc,resret,resmax
      endif
c
!..update the global variables

      do  k = 1,kmax
      do  j = 1,jmax
        q(j,k,7) = itmc(j,k)
        q(j,k,8) = ret(j,k)
      enddo
      enddo

c...apply boundary condition again

      call gammathetabc(q,rho,u,v,ug,vg,xx,xy,yx,yy,im)

      return
      end

c*************************************************************
      subroutine calc_strain_du_ds(strain,du_ds,u,v,vmag,xx,xy,yx,yy,jd,kd)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real strain(jd,kd),du_ds(jd,kd)
      real u(jd,kd),v(jd,kd),vmag(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)

! local variables
      integer j,k,jm1,km1,jp,kp
      real usi,vsi,ueta,veta
      real ux,uy,vx,vy
      real sxx,sxy,syy

      do j = 2,jd-1 
        jp = j + 1
        jm1 = j - 1
        do k = 2,kd-1
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

          strain(j,k) = sqrt(2.*(sxx*sxx + 2.*sxy*sxy + syy*syy))
          du_ds(j,k) = (u(j,k)*u(j,k)*sxx + 2.*u(j,k)*v(j,k)*sxy + 
     &                  v(j,k)*v(j,k)*syy)/(vmag(j,k)*vmag(j,k))
        enddo
      enddo

      return
      end

c***********************************************************************
      subroutine gammathetarhslhs(q,s,itmc,ret,itmcsep,vmul,turmu,vort,
     &                   strain,du_ds,rho,u,v,vmag,ug,vg,xx,xy,yx,yy,sn,
     &                   aj,bjmat,cj,ak,ck,sjmat,jd,kd)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real q(jd,kd,nq),s(jd,kd,nv),itmc(jd,kd),ret(jd,kd),itmcsep(jd,kd)
      real rho(jd,kd), vmul(jd,kd), turmu(jd,kd)
      real vort(jd,kd),strain(jd,kd),du_ds(jd,kd)
      real u(jd,kd),v(jd,kd),vmag(jd,kd),ug(jd,kd),vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real sn(jd,kd)
      real aj(jd,kd,2),bjmat(jd,kd,2),cj(jd,kd,2)
      real ak(jd,kd,2),ck(jd,kd,2),sjmat(jd,kd,2)

      ! local variables
      integer j,k,n,jp,kp,jm1,km1
      real uu,vv,fwd,bck
      real,allocatable :: up(:),um(:),vp(:),vm(:)

      allocate(up(mdim),um(mdim),vp(mdim),vm(mdim))

c***  first executable statement

      do k = 1,kmax
      do j = 1,jmax
        do n = 1,2
          aj(j,k,n)=0.0
          bjmat(j,k,n)=0.0
          cj(j,k,n)=0.0
          ak(j,k,n)=0.0
          ck(j,k,n)=0.0
          sjmat(j,k,n) = s(j,k,n+nmv+2)
        enddo
      enddo
      enddo

c...lhs contribution from the convection term

      do k=2,kmax-1
        do j=1,jmax
          uu=xx(j,k)*(u(j,k)-ug(j,k))+xy(j,k)*(v(j,k)-vg(j,k))
          up(j) = 0.5*(uu+abs(uu))
          um(j) = 0.5*(uu-abs(uu))
        enddo

        do j = 2,jmax-1
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

          sjmat(j,k,1) = sjmat(j,k,1) - up(j)*(itmc(j,k) - itmc(jm1,k)) -
     &                                  um(j)*(itmc(jp,k) - itmc(j,k))
          sjmat(j,k,2) = sjmat(j,k,2) - up(j)*(ret(j,k) - ret(jm1,k)) - 
     &                                  um(j)*(ret(jp,k) - ret(j,k))
          do n = 1,2
            aj(j,k,n) = aj(j,k,n) - up(j)!fwd*(up(jm1)+um(jm1))
            cj(j,k,n) = cj(j,k,n) + um(j)!bck*(up(jp)+um(jp))
            bjmat(j,k,n) = bjmat(j,k,n) + up(j)- um(j)
          enddo
        enddo
      enddo

      do j=2,jmax-1
        do k=1,kmax
          vv = yx(j,k)*(u(j,k)-ug(j,k))+yy(j,k)*(v(j,k)-vg(j,k))
          vp(k) = 0.5*(vv+abs(vv))
          vm(k) = 0.5*(vv-abs(vv))
        enddo

        do k = 2,kmax-1
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

          sjmat(j,k,1) = sjmat(j,k,1) - vp(k)*(itmc(j,k) - itmc(j,km1)) -
     &                                  vm(k)*(itmc(j,kp) - itmc(j,k))
          sjmat(j,k,2) = sjmat(j,k,2) - vp(k)*(ret(j,k) - ret(j,km1)) -
     &                                  vm(k)*(ret(j,kp) - ret(j,k))
          do n = 1,2
            ak(j,k,n) = ak(j,k,n) - vp(k) !fwd*(vp(km1)+vm(km1))
            ck(j,k,n) = ck(j,k,n) + vm(k) !bck*(vp(kp)+vm(kp))
            bjmat(j,k,n) = bjmat(j,k,n) + vp(k) - vm(k)
          enddo
        enddo
      enddo
   
c...rhs and lhs contribution from the source term

      call gammathetasource(q,itmc,ret,itmcsep,rho,vmul,vort,strain,du_ds,vmag,sn,
     &                  bjmat,sjmat,jd,kd)


c...rhs contribution from the diffusion term

      call gammathetadiffus(itmc,ret,rho,vmul,turmu,
     &                      xx,xy,yx,yy,aj,bjmat,cj,ak,ck,sjmat,jd,kd)

      return
      end

c***********************************************************************
      subroutine gammathetadiffus(itmc,ret,rho,vmul,turmu,
     &                            xx,xy,yx,yy,aj,bjmat,cj,ak,ck,sjmat,jd,kd)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
      integer jd,kd
      real itmc(jd,kd),ret(jd,kd),rho(jd,kd)
      real vmul(jd,kd),turmu(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real aj(jd,kd,2),bjmat(jd,kd,2),cj(jd,kd,2)
      real ak(jd,kd,2),ck(jd,kd,2),sjmat(jd,kd,2)
      
    ! local variables
      integer j,k

      real,allocatable :: chp(:,:),dnuhp(:,:)
      real,allocatable :: dmxh(:,:),dmyh(:,:)

      real rei,vnulh,vnuh,dxp,dxm,dcp,dcm,ax,cx,dyp,dym,ay,cy

      include 'gammatheta.h'

      allocate(chp(jd,kd),dnuhp(jd,kd))
      allocate(dmxh(jd,kd),dmyh(jd,kd))

c***  first executable statement

        rei = 1./rey
c       compute j direction differences.

c       compute half-point co-efficients

        do k=1,kmax
        do j=1,jmax-1
          dmxh(j,k) = 0.5*(xx(j,k)+xx(j+1,k))
          dmyh(j,k) = 0.5*(xy(j,k)+xy(j+1,k))
        enddo
        enddo

c     j-direction diffusion for k-equation

        do k=1,kmax
        do j=1,jmax-1
          vnulh     = 0.5*(vmul(j,k)+vmul(j+1,k))
          vnuh     = 0.5*(turmu(j,k)+turmu(j+1,k))
          chp(j,k) = rei*(vnulh+vnuh/sigmaf)       
          dnuhp(j,k) =itmc(j+1,k)/rho(j+1,k)-itmc(j,k)/rho(j,k)
        enddo
        enddo

        do k=1,kmax
        do j=2,jmax-1

          dxp=dmxh(j,k)*xx(j,k)+dmyh(j,k)*xy(j,k)
          dxm=dmxh(j-1,k)*xx(j,k)+dmyh(j-1,k)*xy(j,k)

c         enforce positivity (as suggested by overflow)
          
          dcp    = dxp*(chp(j,k))
          dcm    = dxm*(chp(j-1,k))
          ax=0.0
          cx=0.0

          if(k.ne.kmax.and.k.ne.1) then
            ax       = max(dcm,0.0)
            cx       = max(dcp,0.0)
          endif
c          compute fluxes.

          sjmat(j,k,1)=sjmat(j,k,1)-ax*dnuhp(j-1,k)+cx*dnuhp(j,k)

c          jacobian terms

          aj(j,k,1) = aj(j,k,1) - ax/rho(j-1,k)
          cj(j,k,1) = cj(j,k,1) - cx/rho(j+1,k)
          bjmat(j,k,1) = bjmat(j,k,1)+ (ax + cx)/rho(j,k)

        enddo
        enddo

c      j-direction diffusion for omega-equation

        do k=1,kmax
        do j=1,jmax-1
          vnulh     = 0.5*(vmul(j,k)+vmul(j+1,k))
          vnuh     = 0.5*(turmu(j,k)+turmu(j+1,k))
          chp(j,k) = rei*sigmat*(vnulh+vnuh)       
          dnuhp(j,k) =ret(j+1,k)/rho(j+1,k)-ret(j,k)/rho(j,k)
        enddo
        enddo

        do k=1,kmax
        do j=2,jmax-1

          dxp=dmxh(j,k)*xx(j,k)+dmyh(j,k)*xy(j,k)
          dxm=dmxh(j-1,k)*xx(j,k)+dmyh(j-1,k)*xy(j,k)

c         enforce positivity (as suggested by overflow)
          
          dcp    = dxp*(chp(j,k))
          dcm    = dxm*(chp(j-1,k))
          ax=0.0
          cx=0.0

          if(k.ne.kmax.and.k.ne.1) then
            ax       = max(dcm,0.0)
            cx       = max(dcp,0.0)
          endif
c          compute fluxes.

          sjmat(j,k,2)=sjmat(j,k,2)-ax*dnuhp(j-1,k)+cx*dnuhp(j,k)

c          jacobian terms

          aj(j,k,2) = aj(j,k,2) - ax/rho(j-1,k)
          cj(j,k,2) = cj(j,k,2) - cx/rho(j+1,k)
          bjmat(j,k,2) = bjmat(j,k,2)+ (ax + cx)/rho(j,k)

        enddo
        enddo

c       compute k direction differences.

c       compute half-point co-efficients

        do k=1,kmax-1
        do j=1,jmax
          dmxh(j,k) = 0.5*(yx(j,k)+yx(j,k+1))
          dmyh(j,k) = 0.5*(yy(j,k)+yy(j,k+1))
       enddo
       enddo

c     k-direction diffusion for k-equation

        do k=1,kmax-1
        do j=1,jmax
          vnulh     = 0.5*(vmul(j,k)+vmul(j,k+1))
          vnuh     =  0.5*(turmu(j,k)+turmu(j,k+1))
          chp(j,k) =rei*(vnulh+vnuh/sigmaf)       
          dnuhp(j,k) =itmc(j,k+1)/rho(j,k+1)-itmc(j,k)/rho(j,k)
       enddo
       enddo

       do k=2,kmax-1
       do j=1,jmax

         dyp=dmxh(j,k)*yx(j,k)+dmyh(j,k)*yy(j,k)
         dym=dmxh(j,k-1)*yx(j,k)+dmyh(j,k-1)*yy(j,k)

c      enforce positivity (as suggested by overflow)
       
          dcp    = dyp*(chp(j,k))
          dcm    = dym*(chp(j,k-1))

          ay=0.0
          cy=0.0
          if(j.ne.1.and.j.ne.jmax) then
            ay       = max(dcm,0.0)
            cy       = max(dcp,0.0)
          endif

c        compute fluxes.

          sjmat(j,k,1)=sjmat(j,k,1)-ay*dnuhp(j,k-1)+cy*dnuhp(j,k)

c        jacobian terms

          ak(j,k,1) = ak(j,k,1) - ay/rho(j,k-1)
          ck(j,k,1) = ck(j,k,1) - cy/rho(j,k+1)
          bjmat(j,k,1) = bjmat(j,k,1)+ (ay + cy)/rho(j,k)

        enddo
        enddo

c     k-direction diffusion for omega-equation

        do k=1,kmax-1
        do j=1,jmax
          vnulh     = 0.5*(vmul(j,k)+vmul(j,k+1))
          vnuh     =  0.5*(turmu(j,k)+turmu(j,k+1))
          chp(j,k) =rei*sigmat*(vnulh+vnuh)       
          dnuhp(j,k) =ret(j,k+1)/rho(j,k+1)-ret(j,k)/rho(j,k)
       enddo
       enddo

       do k=2,kmax-1
       do j=1,jmax

         dyp=dmxh(j,k)*yx(j,k)+dmyh(j,k)*yy(j,k)
         dym=dmxh(j,k-1)*yx(j,k)+dmyh(j,k-1)*yy(j,k)

c      enforce positivity (as suggested by overflow)
       
          dcp    = dyp*(chp(j,k))
          dcm    = dym*(chp(j,k-1))

          ay=0.0
          cy=0.0
          if(j.ne.1.and.j.ne.jmax) then
            ay       = max(dcm,0.0)
            cy       = max(dcp,0.0)
          endif

c        compute fluxes.

          sjmat(j,k,2)=sjmat(j,k,2)-ay*dnuhp(j,k-1)+cy*dnuhp(j,k)

c        jacobian terms

          ak(j,k,2) = ak(j,k,2) - ay/rho(j,k-1)
          ck(j,k,2) = ck(j,k,2) - cy/rho(j,k+1)
          bjmat(j,k,2) = bjmat(j,k,2)+ (ay + cy)/rho(j,k)

        enddo
        enddo

      return
      end

c***********************************************************************
      subroutine gammathetasource(q,itmc,ret,itmcsep,rho,vmul,
     &                            vort,strain,du_ds,vmag,
     &                            sn,bjmat,sjmat,jd,kd)
c***********************************************************************
        use params_global
c***********************************************************************
        implicit none
        integer jd,kd
        real q(jd,kd,nq), itmc(jd,kd), ret(jd,kd), itmcsep(jd,kd) 
        real rho(jd,kd), vort(jd,kd),strain(jd,kd),du_ds(jd,kd),vmag(jd,kd)
        real vmul(jd,kd),sn(jd,kd)
        real bjmat(jd,kd,2),sjmat(jd,kd,2)
        
        ! local variables
        integer j,k,iter
        real rei
        real rey_t,rey_tc,re_theta,re_theta_lim
        real rey_t_lim, rey_t_flim1, rey_t_flim2, rey_t_flim3
        real flen,fsublayer,rw,re_v,re_omega,r_t
        real f_onset,f_onset1,f_onset2,f_onset3
        real f_turb,f_lambda,f_wake,f_theta,f_reattach
        real tu,theta,theta_bl,lambda
        real time_scale,delta,delta_bl
        real prod, des
        real var1
        
        include 'gammatheta.h'

c**   first executable statement

        rei = 1./rey
        rey_t_lim = 1870.
        rey_t_flim1 = 400.
        rey_t_flim2 = 596.
        rey_t_flim3 = 1200.
 
        re_theta_lim = 20.*rei

        do k=2,kmax-1
        do j=2,jmax-1

c..SOURCE TERM FOR INTERMITTANCY EQUATION

!...Langtry & Menter's Correlations

!....use dimensional value of re_theta_t in correlation expressions

          rey_t = ret(j,k)/rho(j,k)*rey

!....re_theta_critical correlation

          if(rey_t.le.rey_t_lim) then
            rey_tc = rey_t - (3.96035 - 0.0120656  *rey_t + 
     &                                  868.23e-6  *rey_t**2 -
     &                                  696.506e-9 *rey_t**3 + 
     &                                  174.105e-12*rey_t**4 )
          else
            rey_tc = rey_t - ( 593.11 + 0.482*(rey_t - 1870.0) )
          endif

!....f_length correlation

          if(rey_t.lt.rey_t_flim1) then
            flen = 39.8189 - 0.011927*rey_t - 1.32567e-4*rey_t**2
          elseif(rey_t.ge.rey_t_flim1 .and. rey_t.lt.rey_t_flim2) then
            flen = 263.404 - 1.23939*rey_t + 0.00194548*rey_t**2 - 
     &                 1.01695e-6*rey_t**3
          elseif(rey_t.ge.rey_t_flim2 .and. rey_t.lt.rey_t_flim3) then
            flen = 0.5 - 3e-4*(rey_t - 596.)
          else
            flen = 0.3188
          endif

!....viscous sublayer modification to f_length

          rw = rey*rey*sn(j,k)**2*q(j,k,6)/(500.*vmul(j,k))
          fsublayer = exp(-(2.5*rw)**2)
          flen = flen*(1. - fsublayer) + 40.*fsublayer

!....vorticity reynolds number
          re_v = rey*rho(j,k)*sn(j,k)**2*strain(j,k)/vmul(j,k)

!....f_onset controls transition onset location

          r_t = rho(j,k)*q(j,k,5)/(q(j,k,6)*vmul(j,k))
          f_onset1 = re_v/(2.193*rey_tc)
          f_onset2 = min(max(f_onset1,f_onset1**4),2.)
          f_onset3 = max(1. - (0.4*r_t)**3,0.)
          f_onset = max(f_onset2 - f_onset3, 0.)

!....compute f_turb

          f_turb = exp(-(0.25*r_t)**4)

          prod = flen*c_a1*strain(j,k)*sqrt(f_onset*rho(j,k)*itmc(j,k))
          prod = prod*(1. - c_e1*itmc(j,k)/rho(j,k))

          des = c_a2*itmc(j,k)*vort(j,k)*f_turb
          des = des*(c_e2*itmc(j,k)/rho(j,k) - 1.)

          var1 = 1.5*flen*c_a1*c_e1*strain(j,k)*
     &           sqrt(f_onset*itmc(j,k)/rho(j,k)) + c_a2*vort(j,k)*f_turb
c          var1 = var1 - 0.5*flen*c_e1*strain(j,k)*
c     &           sqrt(f_onset*rho(j,k)/itmc(j,k))
          bjmat(j,k,1) = bjmat(j,k,1) + var1
          sjmat(j,k,1) = sjmat(j,k,1) + prod - des

c..SOURCE TERM FOR RE_THETA EQUATION

!...limit minimum value of turbulent intensity to 0.027 (as recommended by langtry)
          tu = max(100.*sqrt(2./3*q(j,k,5)/rho(j,k))/vmag(j,k),0.027)

!...iterate for value of transition momentum thickness reynolds number, re_theta
!...lambda is the pressure gradient paramter

          f_lambda = 1. ! initial value of lambda taken as f_lambda(lambda=0)

          do iter = 1,10  ! 10 iterations should be sufficient for convergence
            if(tu.le.1.3) then
              re_theta = f_lambda * (1173.51-589.428*tu+0.2196/(tu*tu))
            else
              re_theta = 331.5 * f_lambda*(tu-0.5658)**(-0.671)
            endif
            re_theta = max(re_theta*rei,re_theta_lim)
            if (fmtip.eq.0) then
              theta = re_theta*vmul(j,k)/(rho(j,k)*fsmach)
            else
              theta = re_theta*vmul(j,k)/(rho(j,k)*fmtip)
            endif
            lambda = rey*rho(j,k)*theta*theta*du_ds(j,k)/vmul(j,k)
            lambda = min(max(-0.1,lambda),0.1)

            if(lambda.le.0.) then
              f_lambda = 1. - (-12.986*lambda - 123.66*lambda**2 - 
     &                         405.689*lambda**3)*exp(-(2./3*tu)**1.5)
            else
              f_lambda = 1. + 0.275*(1.-exp(-35.*lambda))*exp(-2.*tu)
            endif
          enddo

!...calculate blending function f_theta

          time_scale = 500.0*vmul(j,k)/(rho(j,k)*vmag(j,k)*vmag(j,k))
          theta_bl = ret(j,k)*vmul(j,k)/(rho(j,k)*rho(j,k)*vmag(j,k))
          delta_bl = 7.5*theta_bl
          delta = 50.0*vort(j,k)*sn(j,k)*delta_bl/vmag(j,k) + 1.e-20
          re_omega = q(j,k,6)*sn(j,k)*sn(j,k)*rey*rey/vmul(j,k)

          f_wake = exp(-(1e-5*re_omega)**2)
               
          var1 = (itmc(j,k)/rho(j,k) - 1./c_e2)/(1. - 1./c_e2)
          var1 = 1. - var1*var1
          f_theta = min(max(f_wake*exp(-(sn(j,k)/delta)**4),var1),1.)

          var1 = c_theta*(1.-f_theta)*rey/time_scale

          bjmat(j,k,2) = bjmat(j,k,2) + var1
          sjmat(j,k,2) = sjmat(j,k,2) + var1*(rho(j,k)*re_theta - ret(j,k))

!...Modification in intermittancy for separation-induced transition

          f_reattach = exp(-(0.05*r_t)**4)
          itmcsep(j,k) = s1*max(0.,re_v/(3.235*rey_tc)-1.)*f_reattach
          itmcsep(j,k) = min(itmcsep(j,k),2.0)*f_theta

        enddo
        enddo


        return
        end

c*************************************************************
        subroutine gammathetabc(q,rho,u,v,ug,vg,xx,xy,yx,yy,im)
c*************************************************************
        use params_global
c*************************************************************
        implicit none
c*************************************************************
        real q(jmax,kmax,nq),rho(jmax,kmax)
        real u(jmax,kmax),v(jmax,kmax),ug(jmax,kmax),vg(jmax,kmax)
        real xx(jmax,kmax),xy(jmax,kmax),yx(jmax,kmax),yy(jmax,kmax)
        integer im

c..   local variables
        integer js,je,ks,ke,idir
        integer ib

        do ib=1,nbc_all(im)
          js = jbcs_all(ib,im)
          je = jbce_all(ib,im)
          ks = kbcs_all(ib,im)
          ke = kbce_all(ib,im)
          if(js.lt.0) js = jmax+js+1
          if(ks.lt.0) ks = kmax+ks+1
          if(je.lt.0) je = jmax+je+1
          if(ke.lt.0) ke = kmax+ke+1
          idir = ibdir_all(ib,im)

c.. inviscid wind tunnel wall bc
          if (ibtyp_all(ib,im).eq.3.or.ibtyp_all(ib,im).eq.4) then
            call gammathetabc_extpt(q,js,je,ks,ke,idir)

c.. wall bc at l = 1 (only interior portion of wall)
          elseif (ibtyp_all(ib,im).eq.5) then
            call gammathetabc_wall(q,js,je,ks,ke,idir)

c.. symmetric bc
          elseif (ibtyp_all(ib,im).eq.10.or.ibtyp_all(ib,im).eq.11) then
            call gammathetabc_sym(q,js,je,ks,ke,idir)

c.. periodic bc
          elseif (ibtyp_all(ib,im).eq.22) then
            call gammathetabc_periodic(q,js,je,ks,ke,idir)

c.. averaging bc for wake
          elseif (ibtyp_all(ib,im).eq.51) then
            call gammathetabc_wake(q,js,je,ks,ke,idir)

c.. averaging bc for wake of O-grid
          elseif (ibtyp_all(ib,im).eq.52) then
            call gammathetabc_wake_ogrid(q,js,je,ks,ke,idir)

c.. freesream bc
          elseif (ibtyp_all(ib,im).eq.47) then
            call gammathetabc_out(q,rho,u,v,ug,vg,xx,xy,yx,yy,js,je,ks,ke,idir)

          endif
        enddo

        return
        end

c*************************************************************
      subroutine gammathetabc_wall(q,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        j = js
        j1 = j + iadd
        do k=ks,ke
          q(j,k,7) = q(j1,k,7)
          q(j,k,8) = q(j1,k,8)
        enddo
      elseif(iadir.eq.2) then
        k = ks
        k1 = k + iadd
        do j=js,je
          q(j,k,7) = q(j,k1,7)
          q(j,k,8) = q(j,k1,8)
        enddo
      endif

      return
      end

c*************************************************************
      subroutine gammathetabc_extpt(q,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        j = js
        j1 = j + iadd
        do k=ks,ke
          q(j,k,7) = q(j1,k,7)
          q(j,k,8) = q(j1,k,8)
        enddo
      elseif(iadir.eq.2) then
        k = ks
        k1 = k + iadd
        do j=js,je
          q(j,k,7) = q(j,k1,7)
          q(j,k,8) = q(j,k1,8)
        enddo
      endif

      return
      end

c*************************************************************
      subroutine gammathetabc_sym(q,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        j = js
        j1 = j + iadd
        do k=ks,ke
          q(j,k,7) = q(j1,k,7)
          q(j,k,8) = q(j1,k,8)
        enddo
      elseif(iadir.eq.2) then
        k = ks
        k1 = k + iadd
        do j=js,je
          q(j,k,7) = q(j,k1,7)
          q(j,k,8) = q(j,k1,8)
        enddo
      endif

      return
      end

c*************************************************************
      subroutine gammathetabc_periodic(q,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(idir.eq.1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = j - js
           jc = jmax - 2*jj + jj1
 
           do k = ks,ke
             q(j,k,7) = q(jc,k,7)
             q(j,k,8) = q(jc,k,8)
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             q(j,k,7) = q(jc,k,7)
             q(j,k,8) = q(jc,k,8)
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks
           kc = kmax - 2*kk + kk1
 
           do j = js,je
             q(j,k,7) = q(j,kc,7)
             q(j,k,8) = q(j,kc,8)
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             q(j,k,7) = q(j,kc,7)
             q(j,k,8) = q(j,kc,8)
           enddo
        enddo

      endif

      return
      end

c*************************************************************
      subroutine gammathetabc_wake(q,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jj,j1,k1,iadd,iadir
      real qav1,qav2

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in gammathetabc_wake'
      elseif(iadir.eq.2) then
        k  = ks
        k1 = k + iadd
        do j=js,je
          jj = jmax - j + 1
          qav1 = 0.5*(q(j,k1,7)+q(jj,k1,7))
          qav2 = 0.5*(q(j,k1,8)+q(jj,k1,8))
          q(j,k,7)  = qav1
          q(j,k,8)  = qav2
          q(jj,k,7) = qav1
          q(jj,k,8) = qav2
        enddo
      endif

      return
      end

c*************************************************************
      subroutine gammathetabc_wake_ogrid(q,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jj,j1,iadd,iadir
      real qav1,qav2

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        j  = js
        j1 = j + iadd
        jj = jmax - j1 + 1
        do k=ks,ke
          qav1 = 0.5*(q(j1,k,7)+q(jj,k,7))
          qav2 = 0.5*(q(j1,k,8)+q(jj,k,8))
          q(j,k,7) = qav1
          q(j,k,8) = qav2
        enddo
      elseif(iadir.eq.2) then
      print*,'idir = ',idir,' is not implemented in 
     c                              gammathetabc_wake_ogrid'
      endif

      return
      end

c*************************************************************
      subroutine gammathetabc_out(q,rho,u,v,ug,vg,xx,xy,yx,yy,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      real q(jmax,kmax,nq),rho(jmax,kmax),u(jmax,kmax),v(jmax,kmax)
      real ug(jmax,kmax),vg(jmax,kmax)
      real xx(jmax,kmax),xy(jmax,kmax),yx(jmax,kmax),yy(jmax,kmax)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,iadd,iadir
      real uu

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        j  = js
        j1 = j + iadd
        do k=ks,ke
          q(j,k,7)=itmcinf*rho(j,k)
          q(j,k,8)=retinf*rho(j,k)
          uu=   (u(j,k)-ug(j,k)+u(j1,k)-ug(j1,k))*xx(j,k)
          uu=uu+(v(j,k)-vg(j,k)+v(j1,k)-vg(j1,k))*xy(j,k)
          uu=uu*iadd
          if(uu.lt.0.) then
            q(j,k,7)=q(j1,k,7)
            q(j,k,8)=q(j1,k,8)
          endif
        enddo
      elseif(iadir.eq.2) then
        k  = ks
        k1 = k + iadd
        do j=js,je
          q(j,k,7)=itmcinf*rho(j,k)
          q(j,k,8)=retinf*rho(j,k)
          uu=   (u(j,k)-ug(j,k)+u(j,k1)-ug(j,k1))*yx(j,k)
          uu=uu+(v(j,k)-vg(j,k)+v(j,k1)-vg(j,k1))*yy(j,k)
          uu=uu*iadd
          if(uu.lt.0.) then
            q(j,k,7)=q(j,k1,7)
            q(j,k,8)=q(j,k1,8)
          endif
        enddo
      endif

      return
      end

c*************************************************************
