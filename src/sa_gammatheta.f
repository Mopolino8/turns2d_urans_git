      module sa_transition_model
      
      use params_global

      implicit none

      private
      public :: sa_gammatheta

      integer,allocatable,private :: check_alloc(:)

      real :: flen_global, alpha_global

      contains

c***********************************************************************
      subroutine init_gammatheta(q,jd,kd)
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)

      real :: tu,re_theta

      if(allocated(check_alloc)) deallocate(check_alloc)

      allocate(check_alloc(1))

      itmcinf = 1.
      flen_global = 12.0
      alpha_global = 0.85

      tu = tuinf
      if(tu.le.1.3) then
        re_theta = (1173.51-589.428*tu+0.2196/(tu*tu))
      else
        re_theta = 331.5*(tu-0.5658)**(-0.671)
      endif
      retinf = max(re_theta,20.)/rey

      if (iread.eq.0) then
         q(:,:,nmv+nturb+1) = itmcinf 
         q(:,:,nmv+nturb+2) = retinf 
      endif
      
      end subroutine init_gammatheta

c***********************************************************************
      subroutine sa_gammatheta(q,s,turmu,sn,vort,xx,xy,yx,yy,ug,vg,
     >                          jd,kd,tscale,iblank,im)
c  correlation-based transition model - Langtry and Menter. 
c  Ref{ AIAA Journal, Vol. 47, No. 12}.
c
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

      if(.not.allocated(check_alloc).or.istep.eq.1) then
        call init_gammatheta(q,jd,kd)
      endif

      itmclim = 1e-20
      retlim = 1e-20
c...laminar co-efficient of viscosity calculation
      
      call lamvis(q,vmul,jd,kd)

      do  k = 1,kd
      do  j = 1,jd
        rho(j,k) = q(j,k,1)*q(j,k,nq)
        u(j,k)  = q(j,k,2)/q(j,k,1)
        v(j,k)  = q(j,k,3)/q(j,k,1)
        vmag(j,k) = sqrt(u(j,k)*u(j,k) + v(j,k)*v(j,k))
        strain(j,k) = 0.0
        du_ds(j,k)  = 0.0
      enddo
      enddo

c...apply boundary condition

      call gammathetabc(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)

c...initiate local working variable for intermittancy and Re-theta

      do  k = 1,kd
      do  j = 1,jd
        itmc(j,k) = q(j,k,nmv+nturb+1)
        ret(j,k) = q(j,k,nmv+nturb+2)
      enddo
      enddo

c...compute strain rate

      call calc_strain_du_ds(strain,du_ds,u,v,vmag,xx,xy,yx,yy,
     &                       jd,kd,jbeg,jend,kbeg,kend)

c...compute rhs and lhs

      call gammathetarhslhs(q,s,itmc,ret,itmcsep,vmul,turmu,vort,strain,
     &                      du_ds,rho,u,v,vmag,ug,vg,xx,xy,yx,yy,
     &                      sn,aj,bjmat,cj,ak,ck,sjmat,
     &                      jd,kd,jbeg,jend,kbeg,kend)

c...invert using DDADI

c..set time-accuracy
      oat = 1.0
      if (ntac.eq.-2) oat = 0.5
      if (ntac.ge.2 .and. istep.gt.1) oat = 2./3.
      if (ntac.eq.3 .and. istep.gt.2) oat = 6./11.

      resitmc=0.0
      resret=0.0
      do k = kbeg,kend
      do j = jbeg,jend

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
        sjmat(j,k,1) = relfac*sjmat(j,k,1)
        sjmat(j,k,2) = relfac*sjmat(j,k,2)

        itmc(j,k) = itmc(j,k) + sjmat(j,k,1)*max(iblank(j,k),0)
        ret(j,k) = ret(j,k) + sjmat(j,k,2)*max(iblank(j,k),0)

        itmc(j,k) = max(itmc(j,k),itmclim)
        ret(j,k) = max(ret(j,k),retlim)

c....limit intermittancy to one
        itmc(j,k) = min(itmc(j,k),1.)

        resitmc = resitmc + sjmat(j,k,1)**2
        resret = resret + sjmat(j,k,2)**2
        if (max(resitmc,resret).gt.resmax) then
          jloc = j
          kloc = k
          resmax = max(resitmc,resret)
        endif

c...Modification in intermittancy for separation-induced transition
        itmc(j,k) = max(itmc(j,k),itmcsep(j,k))

      enddo
      enddo

      resitmc = sqrt(resitmc/jd/kd)
      resret = sqrt(resret/jd/kd)
      resmax = sqrt(resmax/jd/kd)
      if( mod(istep,npnorm).eq.0) then
         write(1333+im,*),istep0,resitmc,resret,resmax
      endif
c
!..update the global variables

      do  k = kbeg,kend
      do  j = jbeg,jend
        q(j,k,nmv+nturb+1) = itmc(j,k)
        q(j,k,nmv+nturb+2) = ret(j,k)
      enddo
      enddo

c...apply boundary condition again

      call gammathetabc(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)

      end subroutine sa_gammatheta

c*************************************************************
      subroutine calc_strain_du_ds(strain,du_ds,u,v,vmag,xx,xy,yx,yy,
     &                             jd,kd,js,je,ks,ke)
c*************************************************************
      integer jd,kd,js,je,ks,ke
      real strain(jd,kd),du_ds(jd,kd)
      real u(jd,kd),v(jd,kd),vmag(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)

! local variables
      integer j,k,jm1,km1,jp,kp
      real usi,vsi,ueta,veta
      real ux,uy,vx,vy
      real sxx,sxy,syy

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

          strain(j,k) = sqrt(2.*(sxx*sxx + 2.*sxy*sxy + syy*syy))
          du_ds(j,k) = (u(j,k)*u(j,k)*sxx + 2.*u(j,k)*v(j,k)*sxy + 
     &                  v(j,k)*v(j,k)*syy)/(vmag(j,k)*vmag(j,k))
        enddo
      enddo

      end subroutine calc_strain_du_ds

c***********************************************************************
      subroutine gammathetarhslhs(q,s,itmc,ret,itmcsep,vmul,turmu,vort,
     &                   strain,du_ds,rho,u,v,vmag,ug,vg,xx,xy,yx,yy,sn,
     &                   aj,bjmat,cj,ak,ck,sjmat,jd,kd,js,je,ks,ke)
c***********************************************************************

      integer jd,kd,js,je,ks,ke
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

      do k = 1,kd
      do j = 1,jd
        do n = 1,2
          aj(j,k,n)=0.0
          bjmat(j,k,n)=0.0
          cj(j,k,n)=0.0
          ak(j,k,n)=0.0
          ck(j,k,n)=0.0
          sjmat(j,k,n) = s(j,k,n+nmv+nturb)
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

      call gammathetasource(q,itmc,ret,itmcsep,rho,vmul,turmu,
     &                  vort,strain,du_ds,vmag,sn,
     &                  bjmat,sjmat,jd,kd,js,je,ks,ke)


c...rhs contribution from the diffusion term

      call gammathetadiffus(itmc,ret,vmul,turmu,xx,xy,yx,yy,
     &                 aj,bjmat,cj,ak,ck,sjmat,jd,kd,js,je,ks,ke)

      end subroutine gammathetarhslhs

c***********************************************************************
      subroutine gammathetadiffus(itmc,ret,vmul,turmu,xx,xy,yx,yy,
     &                 aj,bjmat,cj,ak,ck,sjmat,jd,kd,js,je,ks,ke)
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      real itmc(jd,kd),ret(jd,kd)
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

c...compute j direction differences.

c.....compute half-point co-efficients

      do j=js-1,je
      do k=ks-1,ke+1
        dmxh(j,k) = 0.5*(xx(j,k)+xx(j+1,k))
        dmyh(j,k) = 0.5*(xy(j,k)+xy(j+1,k))
      enddo
      enddo

c.....compute j-direction fluxes for gamma-equation

      do j=js-1,je
      do k=ks-1,ke+1
        vnulh     = 0.5*(vmul(j,k)+vmul(j+1,k))
        vnuh     = 0.5*(turmu(j,k)+turmu(j+1,k))
        chp(j,k) = rei*(vnulh+vnuh/sigmaf)       
        dnuhp(j,k) =itmc(j+1,k)-itmc(j,k)
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

        aj(j,k,1) = aj(j,k,1) - ax
        cj(j,k,1) = cj(j,k,1) - cx
        bjmat(j,k,1) = bjmat(j,k,1)+ (ax + cx)

      enddo
      enddo

c.....boundary terms in jacobians

      do  k = 1,kd
        aj(je+1,k,1) = 0.0
        cj(js-1,k,1) = 0.0
      enddo

c.....compute j-direction fluxes for retheta-equation

      do j=js-1,je
      do k=ks-1,ke+1
        vnulh     = 0.5*(vmul(j,k)+vmul(j+1,k))
        vnuh     = 0.5*(turmu(j,k)+turmu(j+1,k))
        chp(j,k) = rei*sigmat*(vnulh+vnuh)       
        dnuhp(j,k) =ret(j+1,k)-ret(j,k)
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

        aj(j,k,2) = aj(j,k,2) - ax
        cj(j,k,2) = cj(j,k,2) - cx
        bjmat(j,k,2) = bjmat(j,k,2)+ (ax + cx)

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

c.....compute k-direction fluxes for gamma-equation

      do k=ks-1,ke
      do j=js-1,je+1
        vnulh     = 0.5*(vmul(j,k)+vmul(j,k+1))
        vnuh     =  0.5*(turmu(j,k)+turmu(j,k+1))
        chp(j,k) =rei*(vnulh+vnuh/sigmaf)       
        dnuhp(j,k) =itmc(j,k+1)-itmc(j,k)
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

        ak(j,k,1) = ak(j,k,1) - ay
        ck(j,k,1) = ck(j,k,1) - cy
        bjmat(j,k,1) = bjmat(j,k,1)+ (ay + cy)

      enddo
      enddo

c.....boundary terms in jacobians

      do  j = 1,kd
        ak(j,ke+1,1) = 0.0
        ck(j,ks-1,1) = 0.0
      enddo

c.....compute k-direction fluxes for retheta-equation

      do k=ks-1,ke
      do j=js-1,je+1
        vnulh     = 0.5*(vmul(j,k)+vmul(j,k+1))
        vnuh     =  0.5*(turmu(j,k)+turmu(j,k+1))
        chp(j,k) =rei*sigmat*(vnulh+vnuh)       
        dnuhp(j,k) =ret(j,k+1)-ret(j,k)
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

        ak(j,k,2) = ak(j,k,2) - ay
        ck(j,k,2) = ck(j,k,2) - cy
        bjmat(j,k,2) = bjmat(j,k,2)+ (ay + cy)

      enddo
      enddo

c.....boundary terms in jacobians

      do  j = 1,kd
        ak(j,ke+1,2) = 0.0
        ck(j,ks-1,2) = 0.0
      enddo


      end subroutine gammathetadiffus

c***********************************************************************
      subroutine gammathetasource(q,itmc,ret,itmcsep,rho,vmul,turmu,
     &                            vort,strain,du_ds,vmag,
     &                            sn,bjmat,sjmat,jd,kd,js,je,ks,ke)
c***********************************************************************
        integer jd,kd,js,je,ks,ke
        real q(jd,kd,nq), itmc(jd,kd), ret(jd,kd), itmcsep(jd,kd) 
        real rho(jd,kd), vort(jd,kd),strain(jd,kd),du_ds(jd,kd),vmag(jd,kd)
        real vmul(jd,kd),turmu(jd,kd),sn(jd,kd)
        real bjmat(jd,kd,2),sjmat(jd,kd,2)
        
        ! local variables
        integer j,k,iter
        real rei
        real rey_t,rey_tc,re_theta,re_theta_lim
        real rey_t_lim, rey_t_flim1, rey_t_flim2, rey_t_flim3
        real flen,fsublayer,re_v,r_t
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

        do k=ks,ke
        do j=js,je

c..SOURCE TERM FOR INTERMITTANCY EQUATION

!...Langtry & Menter's Correlations

!....use dimensional value of re_theta_t in correlation expressions

          rey_t = ret(j,k)*rey

!....re_theta_critical correlation (using linear relationship for now)

          rey_tc = alpha_global*rey_t

!....f_length correlation (use simple correlation for now)

          flen = flen_global/rey_t

!....vorticity reynolds number
          re_v = rey*rho(j,k)*sn(j,k)**2*strain(j,k)/vmul(j,k)

!....f_onset controls transition onset location

          r_t = turmu(j,k)/vmul(j,k)
          f_onset1 = re_v/(2.193*rey_tc)
          f_onset2 = min(max(f_onset1,f_onset1**4),2.)
          f_onset3 = max(1. - (0.4*r_t)**3,0.)
          f_onset = max(f_onset2 - f_onset3, 0.)

!....compute f_turb

          f_turb = exp(-(0.25*r_t)**4)

          prod = flen*c_a1*strain(j,k)*sqrt(f_onset*itmc(j,k))
          prod = prod*(1. - c_e1*itmc(j,k))

          des = c_a2*itmc(j,k)*vort(j,k)*f_turb
          des = des*(c_e2*itmc(j,k) - 1.)

          var1 = 1.5*flen*c_a1*c_e1*strain(j,k)*
     &           sqrt(f_onset*itmc(j,k)) + c_a2*vort(j,k)*f_turb
          bjmat(j,k,1) = bjmat(j,k,1) + var1
          sjmat(j,k,1) = sjmat(j,k,1) + prod - des

c..SOURCE TERM FOR RE_THETA EQUATION

!...limit minimum value of turbulent intensity to 0.027 (as recommended by langtry)
          tu = tuinf 

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
          theta_bl = ret(j,k)*vmul(j,k)/(rho(j,k)*vmag(j,k))
          delta_bl = 7.5*theta_bl
          delta = 50.0*vort(j,k)*sn(j,k)*delta_bl/vmag(j,k) + 1.e-20

          f_wake = 1. 
               
          var1 = (itmc(j,k) - 1./c_e2)/(1. - 1./c_e2)
          var1 = 1. - var1*var1
          f_theta = min(max(f_wake*exp(-(sn(j,k)/delta)**4),var1),1.)

          var1 = c_theta*(1.-f_theta)*rey/time_scale

          bjmat(j,k,2) = bjmat(j,k,2) + var1
          sjmat(j,k,2) = sjmat(j,k,2) + var1*(re_theta - ret(j,k))

!...Modification in intermittancy for separation-induced transition

          f_reattach = exp(-(0.05*r_t)**4)
          itmcsep(j,k) = s1*max(0.,re_v/(3.235*rey_tc)-1.)*f_reattach
          itmcsep(j,k) = min(itmcsep(j,k),2.0)*f_theta

        enddo
        enddo

        end subroutine gammathetasource

c*************************************************************
        subroutine gammathetabc(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)
c*************************************************************
        integer jd,kd
        real q(jd,kd,nq)
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

c.. wall bc (extrapolation)
          if (ibtyp_all(ib,im).eq.3.or.ibtyp_all(ib,im).eq.4.or.
     &        ibtyp_all(ib,im).eq.5) then
            call gammathetabc_extpt(q,jd,kd,js,je,ks,ke,idir)

c.. symmetric bc
          elseif (ibtyp_all(ib,im).eq.10.or.ibtyp_all(ib,im).eq.11) then
            call gammathetabc_sym(q,jd,kd,js,je,ks,ke,idir)

c.. periodic bc
          elseif (ibtyp_all(ib,im).eq.22) then
            call gammathetabc_periodic(q,jd,kd,js,je,ks,ke,idir)

c.. averaging bc for wake
          elseif (ibtyp_all(ib,im).eq.51) then
            call gammathetabc_wake(q,jd,kd,js,je,ks,ke,idir)

c.. freesream enforcing bc
          elseif (ibtyp_all(ib,im).eq.46) then
            call gammathetabc_inf(q,jd,kd,js,je,ks,ke,idir)

c.. freesream bc
          elseif (ibtyp_all(ib,im).eq.47) then
            call gammathetabc_out(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)

          endif
        enddo

        end subroutine gammathetabc

c*************************************************************
      subroutine gammathetabc_extpt(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,j2,k2,jc,kc,iadd,iadir
      real foso

      foso = 0.0
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
            q(j,k,nmv+nturb+1) = (1.+foso)*q(j1,k,nmv+nturb+1) - foso*q(j2,k,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = (1.+foso)*q(j1,k,nmv+nturb+2) - foso*q(j2,k,nmv+nturb+2)
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
            q(j,k,nmv+nturb+1) = (1.+foso)*q(j,k1,nmv+nturb+1) - foso*q(j,k2,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = (1.+foso)*q(j,k1,nmv+nturb+2) - foso*q(j,k2,nmv+nturb+2)
          enddo
        enddo
      endif

      end subroutine gammathetabc_extpt

c*************************************************************
      subroutine gammathetabc_sym(q,jd,kd,js,je,ks,ke,idir)
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
            q(j,k,nmv+nturb+1) = q(j1,k,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = q(j1,k,nmv+nturb+2)
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
            q(j,k,nmv+nturb+1) = q(j,k1,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = q(j,k1,nmv+nturb+2)
          enddo
        enddo
      endif

      end subroutine gammathetabc_sym

c*************************************************************
      subroutine gammathetabc_periodic(q,jd,kd,js,je,ks,ke,idir)
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
             q(j,k,nmv+nturb+1) = q(jc,k,nmv+nturb+1)
             q(j,k,nmv+nturb+2) = q(jc,k,nmv+nturb+2)
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j + 1
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             q(j,k,nmv+nturb+1) = q(jc,k,nmv+nturb+1)
             q(j,k,nmv+nturb+2) = q(jc,k,nmv+nturb+2)
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks + 1
           kc = kd - 2*kk + kk1
 
           do j = js,je
             q(j,k,nmv+nturb+1) = q(j,kc,nmv+nturb+1)
             q(j,k,nmv+nturb+2) = q(j,kc,nmv+nturb+2)
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k + 1
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             q(j,k,nmv+nturb+1) = q(j,kc,nmv+nturb+1)
             q(j,k,nmv+nturb+2) = q(j,kc,nmv+nturb+2)
           enddo
        enddo

      endif
      end subroutine gammathetabc_periodic

c*************************************************************
      subroutine gammathetabc_wake(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jj,k1,kc,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in gammathetabc_wake'
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
            q(j,k,nmv+nturb+1)  = q(jj,k1,nmv+nturb+1)
            q(j,k,nmv+nturb+2)  = q(jj,k1,nmv+nturb+2)
            q(jj,k,nmv+nturb+1) = q(j,k1,nmv+nturb+1)
            q(jj,k,nmv+nturb+2) = q(j,k1,nmv+nturb+2)
          enddo
        enddo

      endif

      end subroutine gammathetabc_wake

c*************************************************************
      subroutine gammathetabc_inf(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        do j = js,je
          do k = ks,ke
            q(j,k,nmv+nturb+1)=itmcinf
            q(j,k,nmv+nturb+2)=retinf
          enddo
        enddo
      elseif(iadir.eq.2) then
        do k = ks,ke
          do j = js,je
            q(j,k,nmv+nturb+1)=itmcinf
            q(j,k,nmv+nturb+2)=retinf
          enddo
        enddo
      endif

      end subroutine gammathetabc_inf

c*************************************************************
      subroutine gammathetabc_out(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq),u(jd,kd),v(jd,kd)
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
          q(j,k,nmv+nturb+1) = itmcinf
          q(j,k,nmv+nturb+2) = retinf
          uu = (u(j,k)-ug(j,k)+u(j1,k)-ug(j1,k))*xx(j,k)
          uu = uu+(v(j,k)-vg(j,k)+v(j1,k)-vg(j1,k))*xy(j,k)
          uu = uu*iadd
          if(uu.lt.0.) then
            q(j,k,nmv+nturb+1) = q(j1,k,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = q(j1,k,nmv+nturb+2)
          endif
        enddo

c..extrapolate to other halo cells

        do jc = 1,je-js
          if (idir.eq.1) then
            j = je - jc
          elseif (idir.eq.-1) then
            j = js + jc
          endif
          j1 = j  + iadd

          do k = ks,ke
            q(j,k,nmv+nturb+1) = q(j1,k,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = q(j1,k,nmv+nturb+2)
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
          q(j,k,nmv+nturb+1) = itmcinf
          q(j,k,nmv+nturb+2) = retinf
          uu=   (u(j,k)-ug(j,k)+u(j,k1)-ug(j,k1))*yx(j,k)
          uu=uu+(v(j,k)-vg(j,k)+v(j,k1)-vg(j,k1))*yy(j,k)
          uu=uu*iadd
          if(uu.lt.0.) then
            q(j,k,nmv+nturb+1) = q(j,k1,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = q(j,k1,nmv+nturb+2)
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
            q(j,k,nmv+nturb+1) = q(j,k1,nmv+nturb+1)
            q(j,k,nmv+nturb+2) = q(j,k1,nmv+nturb+2)
          enddo

        enddo

      endif

      end subroutine gammathetabc_out

c*************************************************************
      end module sa_transition_model
