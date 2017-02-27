c***********************************************************************
      subroutine vmu_sa(q,s,turmu,x,y,xv,yv,xx,xy,yx,yy,ug,vg,jd,kd,
     >                  tscale,iblank,im,nsp)
c  turbulent eddy viscosity. model is one equation spalart-
c  allmaras. ref{ aiaa 92-0439}.
c
c***********************************************************************
      use params_global
      use sa_transition_model
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer im,jd,kd,nsp
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd), tscale(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd) 
      integer iblank(jd,kd)

      ! local variables
      real,allocatable :: sjmat(:,:),bjmat(:,:)
      real,allocatable :: aj(:,:),bj(:,:),cj(:,:)
      real,allocatable :: ak(:,:),bk(:,:),ck(:,:)
      real,allocatable :: vmul(:,:),vort(:,:),sn(:,:)
      real,allocatable :: u(:,:),v(:,:)
      real,allocatable :: s1(:,:),s2(:,:)
      real,allocatable :: trip(:,:)
      
      integer k,j,jloc,kloc,niter,ihuang,logv,nhuang
      real relfac,vnulim,dmaxtur,dl2norm
      real resl2,tscal,dtpseudo_turb

      allocate(sjmat(jd,kd),bjmat(jd,kd))
      allocate(aj(jd,kd),bj(jd,kd),cj(jd,kd))
      allocate(ak(jd,kd),bk(jd,kd),ck(jd,kd))
      allocate(vmul(jd,kd),vort(jd,kd),sn(jd,kd))
      allocate(u(jd,kd),v(jd,kd))
      allocate(s1(jd,kd),s2(jd,kd))
      allocate(trip(jd,kd))
c***  first executable statement

      relfac=1.0
      vnulim=1.0e-20
      nturiter=1
      nhuang=1


      if(obj_samod .and.(mod(istep0,nrest).eq.0.or.istep0.eq.2)) then
         call store_coeff_prod(x, y, obj_coeff_prod(1:jd,1:kd,im), jd, kd,nsp,im)
      end if

c...laminar co-efficient of viscosity calculation
      
      call lamvis(q,vmul,jd,kd)

c...for compatibility with turbulence model

      do  k = 1,kd
        do  j = 1,jd
          u(j,k)  = q(j,k,2)/q(j,k,1)
          v(j,k)  = q(j,k,3)/q(j,k,1)
          sjmat(j,k) = 0.0
          vort(j,k) = 0.0
          sn(j,k) = 1.e10

          aj(j,k) = 0.
          bj(j,k) = 0.
          cj(j,k) = 0.
          ak(j,k) = 0.
          bk(j,k) = 0.
          ck(j,k) = 0.
          bjmat(j,k) = 0.
        enddo
      enddo

      call turbc(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)

c...compute vorticity & distance function
       
      if(bodyflag(im)) then 
        call vortic(vort,u,v,xx,xy,yx,yy,jd,kd,jbeg,jend,kbeg,kend)
        call dist(sn,x,y,xv,yv,jd,kd)
      endif

      call c_turm(q,turmu,vmul,jd,kd)
      call ftrip(q,vort,sn,trip,u,v,jd,kd)

      do 30 niter=1,nturiter

        if (itrans.eq.1) then
          call sa_gammatheta(q,s,turmu,sn,vort,xx,xy,yx,yy,ug,vg,
     >                         jd,kd,tscale,iblank,im)
        endif

        call rhslhs(q,s,sn,vort,u,v,vmul,trip,sjmat,xx,xy,yx,yy,
     &              ug,vg,aj,bj,cj,ak,bk,ck,jd,kd,jbeg,jend,kbeg,kend,
     &              nsp,im)

c   scale by time-step

        dl2norm=0.0

        do k = kbeg,kend
        do j = jbeg,jend

          dl2norm=dl2norm+sjmat(j,k)*sjmat(j,k)
          tscal = tscale(j,k)
          
          if(cfltur.gt.1.and.idual_turb.eq.1) then
            dtpseudo_turb=cfltur/(bj(j,k)+bk(j,k))
            tscal=dtpseudo_turb*float(iblank(j,k))
            tscal=tscal/(1.+tscal/h)
          endif

          if(timeac.eq.1) then
            dtpseudo_turb=1.0
            tscal = max(iblank(j,k),0)*( 1.0 + 0.002*sqrt(q(j,k,nq)))
     <                         /(1.+sqrt(q(j,k,nq)))
            tscal = tscal*dtpseudo_turb
            tscal = tscal/(1.+tscal/h)
          endif

          sjmat(j,k)  = sjmat(j,k)*tscal
          aj(j,k) = aj(j,k)*tscal
          cj(j,k) = cj(j,k)*tscal
          ak(j,k) = ak(j,k)*tscal
          ck(j,k) = ck(j,k)*tscal
          bjmat(j,k) = 1.+(bj(j,k)+bk(j,k))*tscal
          s1(j,k) = 0.
          s2(j,k) = 0.

        enddo
        enddo

        do ihuang=1,nhuang

          if(ihuang.gt.1) then
            do k = kbeg,kend
            do j = jbeg,jend
              sjmat(j,k)  = bjmat(j,k)*(sjmat(j,k)-s1(j,k))
            enddo
            enddo
          endif

          call lsolvej(aj,bjmat,cj,sjmat,jd,kd,jbeg,jend,kbeg,kend)

          do k = kbeg,kend
          do j = jbeg,jend
            s1(j,k) = s1(j,k)+sjmat(j,k)
            sjmat(j,k)  = (s1(j,k)-s2(j,k))*bjmat(j,k)
          enddo
          enddo

          call lsolvek(ak,bjmat,ck,sjmat,jd,kd,jbeg,jend,kbeg,kend)

          do k = kbeg,kend
          do j = jbeg,jend
            s2(j,k)  = s2(j,k)+sjmat(j,k)
            sjmat(j,k)  = s2(j,k)
          enddo
          enddo

        enddo

        tscheck: if (.not.timespectral.or.itn.eq.1) then
          resl2=0.0
          do k = kbeg,kend
          do j = jbeg,jend
            sjmat(j,k) = relfac*sjmat(j,k)
            resl2=resl2+sjmat(j,k)*sjmat(j,k)
          enddo
          enddo
          resl2=sqrt(resl2/jd/kd)
          dl2norm=sqrt(dl2norm/jd/kd)

          if( mod(istep0,npnorm).eq.0 ) then
!$OMP ORDERED
            if(niter.eq.nturiter) write(1233+im,*)istep0,resl2,dl2norm
!$OMP END ORDERED
          endif
        endif tscheck
        
        do k = kbeg,kend
        do j = jbeg,jend
          q(j,k,nmv+1) = max( (q(j,k,nmv+1) + sjmat(j,k)*max(iblank(j,k),0)),vnulim)
        enddo
        enddo

        call turbc(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)

30      continue

      call c_turm(q,turmu,vmul,jd,kd)

      return
      end

c*************************************************************
      subroutine lamvis(q,vmul,jd,kd)
c
c     calculate laminar viscosity
c
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************

      integer jd,kd
      real q(jd,kd,nq),vmul(jd,kd)
        
! local variables
        
      integer k,j
      real gkpr,prtr,dre,c2b,c2bp,ra,uvel,vvel
      real ei,tt

c***  first executable statement
        
      gkpr = gamma/pr
      prtr = pr/0.9
      dre  = .5/rey
      c2b  =198.6/tinf
      c2bp = c2b +1.

C$AD II-LOOP
      do 10 k = 1,kd
C$AD II-LOOP
      do 10 j = 1,jd
         ra    = 1./q(j,k,1)
         uvel  = q(j,k,2)*ra
         vvel  = q(j,k,3)*ra
         ei    =   q(j,k,4)*ra-0.5*(uvel**2+vvel**2)
         tt    = ggm1*ei
         vmul(j,k)  = c2bp*tt*sqrt(tt)/(c2b + tt)
10    continue

      return
      end

c*************************************************************
      subroutine dist(sn,x,y,xv,yv,jd,kd)
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
     &               (y(j,k)-0.5*(yv(jv,1)+yv(jv+1,1)))**2)
      enddo
      enddo


      do j=1,jtail1-1
      do k=1,kd
        jv = jtail1 - nhalo
        sn(j,k)=sqrt((x(j,k)-xv(jv,1))**2 + (y(j,k)-yv(jv,1))**2)
c        sn(j,k)=1.e10
      enddo
      enddo

      do j=jtail2+1,jd
      do k=1,kd
        jv = jtail1 - nhalo
        sn(j,k)=sqrt((x(j,k)-xv(jv,1))**2 + (y(j,k)-yv(jv,1))**2)
c        sn(j,k)=1.e10
      enddo
      enddo

      return
      end

c*************************************************************
      subroutine vortic(vort,u,v,xx,xy,yx,yy,jd,kd,js,je,ks,ke)
c*************************************************************
        use params_global
c*************************************************************
        implicit none
c*************************************************************
        integer jd,kd,js,je,ks,ke
        real vort(jd,kd)
        real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
        real u(jd,kd),v(jd,kd)

        ! local variables
        integer j,k
        real dx2,dy2,ux,vx,uy,vy,tx

c**   first executable statement

        dx2     = 0.5
        dy2     = 0.5
      
C$AD II-LOOP
        do 10  k = ks,ke
C$AD II-LOOP
        do 10  j = js,je

          ux = (u(j+1,k) - u(j-1,k)) * dx2
          vx = (v(j+1,k) - v(j-1,k)) * dx2
 
          uy = (u(j,k+1) - u(j,k-1)) * dy2
          vy = (v(j,k+1) - v(j,k-1)) * dy2

          tx = xy(j,k)*ux -xx(j,k)*vx +yy(j,k)*uy -yx(j,k)*vy
          vort(j,k) = abs(tx)

10      continue

      return
      end

c*************************************************************
      subroutine ftrip(q,vort,sn,trip,u,v,jd,kd)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************

      integer jd,kd
      real q(jd,kd,nq),u(jd,kd),v(jd,kd)
      real vort(jd,kd),sn(jd,kd)
      real trip(jd,kd)

      include 'sadata.h'

      ! local variables
      integer j1,j2,j,k
      real dx2,dy2
      
c***  first executable statement
      
      j1  = jtail1
      j2  = jtail2
      dx2 = 0.5
      dy2 = 0.5

C$AD II-LOOP
      do 10  k   = 1,kd
C$AD II-LOOP
      do 10  j   = 1,jd
        trip(j,k)=0.0
10    continue

      return
      end

c*************************************************************
      subroutine c_turm(q,turmu,vmul,jd,kd)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      real turmu(jd,kd),vmul(jd,kd)

! local variables
      integer j,k
      real chi,vnul,fv1

      include 'sadata.h'

C$AD II-LOOP
      do 10  k   = 1,kd
C$AD II-LOOP
      do 10  j   = 1,jd
        vnul = vmul(j,k)/(q(j,k,1)*q(j,k,nq))
        chi  = q(j,k,nmv+1)/vnul
        fv1 = chi**3/(chi**3+cv1**3)
        turmu(j,k)=q(j,k,1)*q(j,k,nq)*q(j,k,nmv+1)*fv1
10    continue

      return
      end

c***********************************************************************
      subroutine rhslhs(q,s,sn,vort,u,v,vmul,trip,sjmat,xx,xy,yx,yy,
     &                  ug,vg,aj,bj,cj,ak,bk,ck,jd,kd,js,je,ks,ke,nsp,im)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,nsp,im
      real q(jd,kd,nq), s(jd,kd,nv), sjmat(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real aj(jd,kd),bj(jd,kd),cj(jd,kd)
      real ak(jd,kd),bk(jd,kd),ck(jd,kd)
      real vmul(jd,kd),vort(jd,kd),sn(jd,kd)
      real u(jd,kd),v(jd,kd)
      real trip(jd,kd)

      ! local variables
      integer j,k

c***  first executable statement

C$AD II-LOOP
      do k = 1,kd
C$AD II-LOOP
      do j = 1,jd
        aj(j,k)=0.0
        bj(j,k)=0.0
        cj(j,k)=0.0
        ak(j,k)=0.0
        bk(j,k)=0.0
        ck(j,k)=0.0
        sjmat(j,k) = s(j,k,nmv+1)
      enddo
      enddo


      call convec(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
     &            ug,vg,aj,bj,cj,ak,bk,ck,jd,kd,js,je,ks,ke)

      call diffus(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
     &            aj,bj,cj,ak,bk,ck,jd,kd,js,je,ks,ke)

      call source(q,sn,u,v,sjmat,vmul,vort,xx,xy,yx,yy,
     &            aj,bj,cj,ak,bk,ck,jd,kd,js,je,ks,ke,nsp,im)

      return
      end

c***********************************************************************
      subroutine convec(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
     &                  ug,vg,aj,bj,cj,ak,bk,ck,jd,kd,js,je,ks,ke)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq), sjmat(jd,kd)
      real ug(jd,kd),vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real aj(jd,kd),bj(jd,kd),cj(jd,kd)
      real ak(jd,kd),bk(jd,kd),ck(jd,kd)

      real vmul(jd,kd),sn(jd,kd)
      real u(jd,kd),v(jd,kd)

      !local variables
      real :: um(jd,kd),up(jd,kd),vm(jd,kd),vp(jd,kd)

      integer j,k
      real uu,vv,fwd,bck

c***  first executable statement

      do 10 k=ks-1,ke+1
      do 10 j=js-1,je+1
        uu=xx(j,k)*(u(j,k)-ug(j,k))+xy(j,k)*(v(j,k)-vg(j,k))
        vv=yx(j,k)*(u(j,k)-ug(j,k))+yy(j,k)*(v(j,k)-vg(j,k))
        up(j,k)=0.5*(uu+abs(uu))
        um(j,k)=0.5*(uu-abs(uu))
        vp(j,k)=0.5*(vv+abs(vv))
        vm(j,k)=0.5*(vv-abs(vv))
 10   continue

      do 20 k=ks,ke
      do 20 j=js,je
        sjmat(j,k)=sjmat(j,k)-up(j,k)*(q(j,k,nmv+1)-q(j-1,k,nmv+1))
     &                       -um(j,k)*(q(j+1,k,nmv+1)-q(j,k,nmv+1))
        sjmat(j,k)=sjmat(j,k)-vp(j,k)*(q(j,k,nmv+1)-q(j,k-1,nmv+1))
     &                       -vm(j,k)*(q(j,k+1,nmv+1)-q(j,k,nmv+1))
 20   continue

c  j-direction jacobian terms

      do 30 k=ks-1,ke+1
      do 30 j=js,je

        if(up(j+1,k).gt.1.0e-12) then
          fwd=1.0
        else
          fwd=0.0
        endif

        if(um(j-1,k).lt.-1.0e-12) then
          bck=1.0
        else
          bck=0.0
        endif

        aj(j+1,k)=aj(j+1,k)-fwd*(up(j,k)+um(j,k))
        bj(j,k)=  bj(j,k)+ up(j,k)-um(j,k)
        cj(j-1,k)=cj(j-1,k)+bck*(up(j,k)+um(j,k))
30    continue

c  k-direction jacobian terms

      do 40 k=ks,ke
      do 40 j=js-1,je+1

        if(vp(j,k+1).gt.1.0e-12) then
          fwd=1.0
        else
          fwd=0.0
        endif

        if(vm(j,k-1).lt.-1.0e-12) then
          bck=1.0
        else
          bck=0.0
        endif

        ak(j,k+1)=ak(j,k+1)-fwd*(vp(j,k)+vm(j,k))
        bk(j,k)=  bk(j,k)+ vp(j,k)-vm(j,k)
        ck(j,k-1)=ck(j,k-1)+bck*(vp(j,k)+vm(j,k))
40    continue

      return
      end

c***********************************************************************
      subroutine diffus(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
     &                  aj,bj,cj,ak,bk,ck,jd,kd,js,je,ks,ke)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq), sjmat(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real aj(jd,kd),bj(jd,kd),cj(jd,kd)
      real ak(jd,kd),bk(jd,kd),ck(jd,kd)
      
      real vmul(jd,kd),sn(jd,kd)
      real u(jd,kd),v(jd,kd)
      
      include 'sadata.h'

      ! local variables
      real :: chp(jd,kd),dnuhp(jd,kd)
      real :: dmxh(jd,kd),dmyh(jd,kd)
        
      integer j,k
      real vnulh,vnuh,dxp,dxm,c2,dcp,dcm,ax,cx,dyp,dym,ay,cy

c***  first executable statement
      
c...compute j direction differences.

c.....compute half-point co-efficients

      do 10 k=ks-1,ke + 1
      do 10 j=js-1,je

        dmxh(j,k) = 0.5*(xx(j,k)+xx(j+1,k))
        dmyh(j,k) = 0.5*(xy(j,k)+xy(j+1,k))
        vnulh     = 0.5*(vmul(j,k)/q(j,k,1)/q(j,k,nq)+
     &                   vmul(j+1,k)/q(j+1,k,1)/q(j+1,k,nq))
        vnuh = 0.5*(q(j,k,nmv+1)+q(j+1,k,nmv+1))
        chp(j,k) = (1.+cb2)/sigma/rey*(vnulh+vnuh)
        dnuhp(j,k) = q(j+1,k,nmv+1)-q(j,k,nmv+1)

10    continue

c.....compute j direction fluxes

      do 20 j = js,je
      do 20 k = ks-1,ke + 1
      
        dxp = dmxh(j,k)*xx(j,k)   + dmyh(j,k)*xy(j,k)
        dxm = dmxh(j-1,k)*xx(j,k) + dmyh(j-1,k)*xy(j,k)
        
        c2 = cb2/sigma/rey
        c2 = c2*(vmul(j,k)/q(j,k,1)/q(j,k,nq)+q(j,k,nmv+1))

c.....enforce positivity (as suggested by overflow)

        dcp = dxp*(chp(j,k)-c2)
        dcm = dxm*(chp(j-1,k)-c2)
        ax = 0.0
        cx = 0.0
        if (k.ne.ke+1.and.k.ne.ks-1) then
            ax = max(dcm,0.0)
            cx = max(dcp,0.0)
        endif

c.....compute fluxes.
        sjmat(j,k) = sjmat(j,k)-ax*dnuhp(j-1,k)+cx*dnuhp(j,k)

c.....jacobian terms

        aj(j,k) = aj(j,k) - ax
        cj(j,k) = cj(j,k) - cx
        bj(j,k) = bj(j,k) + ax + cx

20    continue

c.....boundary terms in jacobians

      do 30 k = 1,kd
        aj(je+1,k) = 0.0
        cj(js-1,k) = 0.0
30    continue

c...compute k direction differences.

c.....compute half-point co-efficients

      do 40 k=ks-1,ke
      do 40 j=js-1,je+1

        dmxh(j,k) = 0.5*(yx(j,k)+yx(j,k+1))
        dmyh(j,k) = 0.5*(yy(j,k)+yy(j,k+1))
        vnulh = 0.5*(vmul(j,k)/q(j,k,1)/q(j,k,nq)+
     &               vmul(j,k+1)/q(j,k+1,1)/q(j,k+1,nq))
        vnuh = 0.5*(q(j,k,nmv+1)+q(j,k+1,nmv+1))
        chp(j,k) = (1.+cb2)/sigma/rey*(vnulh+vnuh)
        dnuhp(j,k) = q(j,k+1,nmv+1)-q(j,k,nmv+1)

40    continue

c.....compute k direction fluxes

      do 50 k = ks,ke
      do 50 j = js-1,je+1

        dyp = dmxh(j,k)*yx(j,k)   + dmyh(j,k)*yy(j,k)
        dym = dmxh(j,k-1)*yx(j,k) + dmyh(j,k-1)*yy(j,k)
        
        c2 = cb2/sigma/rey
        c2 = c2*(vmul(j,k)/q(j,k,1)/q(j,k,nq)+q(j,k,nmv+1))

c.....enforce positivity (as suggested by overflow)
       
        dcp = dyp*(chp(j,k)-c2)
        dcm = dym*(chp(j,k-1)-c2)
       
        ay = 0.0
        cy = 0.0
        if (j.ne.js-1.and.j.ne.je+1) then
          ay = max(dcm,0.0)
          cy = max(dcp,0.0)
        endif

c.....compute fluxes.

        sjmat(j,k) = sjmat(j,k)-ay*dnuhp(j,k-1)+cy*dnuhp(j,k)

c.....jacobian terms

        ak(j,k) = ak(j,k) - ay
        ck(j,k) = ck(j,k) - cy
        bk(j,k) = bk(j,k) + ay + cy

50    continue

c.....boundary terms in jacobians

      do 60 j=1,jd
        ak(j,ke+1) = 0.0
        ck(j,ks-1) = 0.0
60    continue

      return
      end

c***********************************************************************
      subroutine source(q,sn,u,v,sjmat,vmul,vort,xx,xy,yx,yy,
     &                  aj,bj,cj,ak,bk,ck,jd,kd,js,je,ks,ke,nsp,im)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
      integer jd,kd,js,je,ks,ke,nsp,im
      real q(jd,kd,nq), sjmat(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real aj(jd,kd),bj(jd,kd),cj(jd,kd)
      real ak(jd,kd),bk(jd,kd),ck(jd,kd)
      
      real vmul(jd,kd),sn(jd,kd),vort(jd,kd)
      real u(jd,kd),v(jd,kd)
      
      include 'sadata.h'

      ! local variables
      real :: chp(jd,kd),dnuhp(jd,kd)
      real :: dmxh(jd,kd),dmyh(jd,kd)
      real :: prod_store(jd, kd), dest_store(jd, kd)
      real rcv2,d2min,rmax,stilim,fturf
      real vnul,chi,d,fv1,fv2,fv3,ft2,dchi
      real dfv1,dfv2,dfv3,dft2,d2,stilda,r,g,fw
      real dstild,dr,dg,dfw,pro,des,prod,dest,dpro,ddes
      real tk1,tk2
      integer j,k,isour

c**   first executable statement
      prod_store(:,:) = 0.0
      dest_store(:,:) = 0.0

      rcv2=(1./5.)
      d2min=1.e-12
      rmax=10.0
      stilim=1.e-12
      fturf=0.0
      isour = 0

C$AD II-LOOP
      do 10 k=ks,ke
C$AD II-LOOP
      do 10 j=js,je
      vnul = vmul(j,k)/q(j,k,1)/q(j,k,nq)
      chi  = q(j,k,nmv+1)/vnul
      d    = sn(j,k)
      chi  = max(chi,1.e-12)
      fv1  = (chi**3)/(chi**3+cv1**3)

      if(isour.eq.0) then
        fv2  = 1.-chi/(1.+chi*fv1)
        fv3  = 1.0
        ft2  = fturf*ct3*exp(-1.*ct4*chi*chi)
        dchi = 1./vnul
        dfv1 = (3.*cv1**3)*(chi**2)*dchi*(1./(chi**3+cv1**3))**2
        dfv2 = (fv2-1.)/chi/vnul+(1.-fv2)*(1-fv2)*(dfv1+fv1/chi/vnul)
        dfv3 = 0.0
      else
        fv2  = 1./(1.+chi*rcv2)
        fv2  = fv2**3
        fv3  = (1.+chi*fv1)*(1.-fv2)/chi
        ft2  = fturf*ct3*exp(-1.*ct4*chi*chi)
        dchi = 1./vnul
        dfv1 = (3.*cv1**3)*(chi**2)*dchi*(1./(chi**3+cv1**3))**2
        dfv2 = (-3.*rcv2)*fv2*dchi/(1.+chi*rcv2)
        dfv3 = ( (chi*dfv1 + dchi*fv1)*(1.-fv2) 
     &       - (1.+chi*fv1)*dfv2- dchi*fv3)/chi 
      endif

      dft2 = (-2.*ct4)*chi*dchi*ft2
      d2   = max(d**2,d2min)

c....for new definition of s_{tilda}, refer aiaa-95-0312

      stilda = vort(j,k)*fv3 + vnul/(d2*akt*akt*rey)*chi*fv2
      stilda = max(stilda,stilim)

      r  = q(j,k,nmv+1)/(d2*akt*akt*rey)/stilda
      r  = min(r,rmax)
      g  = r + cw2*(r**6 - r)
      fw = (1.+cw3**6)/(g**6+cw3**6)
      fw = g*(fw**(1./6.))

       dstild = vort(j,k)*dfv3+vnul*(dchi*fv2+chi*dfv2)/(d2*akt*akt*rey)
       dr  = vnul*(dchi-chi*dstild/stilda)/stilda/(d2*akt*akt*rey)
       dg  = dr*(1.+cw2*(6.*r**5.-1.))
       dfw = ((1.+cw3**6)/(g**6+cw3**6))**(1./6.)
       dfw =  dfw*dg*(1.- g**6/(g**6+cw3**6))
       
       pro  = cb1*stilda*(1.-ft2)
       des  = (cw1*fw-cb1/(akt*akt)*ft2)/d2/rey

       prod = cb1*stilda*(1.-ft2)*q(j,k,nmv+1)
       dest = (cw1*fw-cb1/(akt*akt)*ft2)*q(j,k,nmv+1)*q(j,k,nmv+1)/d2/rey
       if(obj_samod) then
          prod_store(j,k) = prod
          dest_store(j,k) = dest

          prod = prod*obj_coeff_prod(j,k,im)
       end if



       dpro = pro*dstild/stilda - cb1*stilda*dft2
  
       ddes = (cw1*dfw-cb1/(akt*akt)*dft2)/d2/rey*vnul*chi
       ddes = ddes + des
       ddes = ddes*q(j,k,nmv+1)
       dpro = dpro*q(j,k,nmv+1)

c...modification for transition model

        if (itrans.eq.1) then
          prod = q(j,k,nmv+nturb+1)*prod
          dest = min(max(q(j,k,nmv+nturb+1),0.1),1.0)*dest
        endif

        sjmat(j,k) = sjmat(j,k) + prod-dest
        tk1=max(des*q(j,k,nmv+1)-pro,0.0)
        tk2=max(ddes-dpro,0.0)
        bk(j,k) = bk(j,k)+tk1+tk2
c       bj(j,k) = bj(j,k)+tk1+tk2

10    continue


      if(obj_samod .and.(mod(istep0,nrest).eq.0.or.istep0.eq.2)) then
         call store_saproduction(prod_store, jd, kd,nsp,im)
         !amm call store_sadestruction(dest_store, jd, kd,nsp,im)
      end if


      return
      end

c*************************************************************
      subroutine turbc(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,im)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      real u(jd,kd),v(jd,kd),ug(jd,kd),vg(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      integer im

c..   local variables
      integer js,je,ks,ke,idir
      integer j,k,ib

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

c.. inviscid wind tunnel wall bc
        if (ibtyp_all(ib,im).eq.3.or.ibtyp_all(ib,im).eq.4) then
          call turbc_extpt(q,jd,kd,js,je,ks,ke,idir)

c.. wall bc at l = 1 (only interior portion of wall)
        elseif (ibtyp_all(ib,im).eq.5) then
          call turbc_wall(q,jd,kd,js,je,ks,ke,idir)

c.. symmetric bc
        elseif (ibtyp_all(ib,im).eq.10.or.ibtyp_all(ib,im).eq.11) then
          call turbc_sym(q,jd,kd,js,je,ks,ke,idir)

c.. periodic bc
        elseif (ibtyp_all(ib,im).eq.22) then
          call turbc_periodic(q,jd,kd,js,je,ks,ke,idir)

c.. averaging bc for wake
        elseif (ibtyp_all(ib,im).eq.51) then
          call turbc_wake(q,jd,kd,js,je,ks,ke,idir)

c.. freesream enforcing bc
        elseif (ibtyp_all(ib,im).eq.46) then
          call turbc_inf(q,jd,kd,js,je,ks,ke,idir)

c.. freesream characteristic bc
        elseif (ibtyp_all(ib,im).eq.47) then
          call turbc_out(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)

        endif
      enddo

      return
      end

c*************************************************************
      subroutine turbc_wall(q,jd,kd,js,je,ks,ke,idir)
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
      real t,foso

      foso = 1.0
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        if (idir.eq.1) then
          j = je
        elseif (idir.eq.-1) then
          j = js
        endif
        j1 = j + iadd
        do k=ks,ke
          q(j,k,nmv+1) = -q(j1,k,nmv+1)
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
          q(j,k,nmv+1) = -q(j,k1,nmv+1)
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
          enddo
        enddo

      endif

      return
      end

c*************************************************************
      subroutine turbc_extpt(q,jd,kd,js,je,ks,ke,idir)
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
          enddo
        enddo
      endif

      return
      end

c*************************************************************
      subroutine turbc_sym(q,jd,kd,js,je,ks,ke,idir)
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
          enddo
        enddo
      endif

      return
      end

c*************************************************************
      subroutine turbc_periodic(q,jd,kd,js,je,ks,ke,idir)
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
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j + 1
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             q(j,k,nmv+1) = q(jc,k,nmv+1)
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks + 1
           kc = kd - 2*kk + kk1
 
           do j = js,je
             q(j,k,nmv+1) = q(j,kc,nmv+1)
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k + 1
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             q(j,k,nmv+1) = q(j,kc,nmv+1)
           enddo
        enddo

      endif

      return
      end

c*************************************************************
      subroutine turbc_wake(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jj,j1,k1,kc,iadd,iadir

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in turbc_wake'
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
            q(jj,k,nmv+1) = q(j,k1,nmv+1)
          enddo
        enddo

      endif

      return
      end

c*************************************************************
      subroutine turbc_inf(q,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
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
            q(j,k,nmv+1)=vnuinf
          enddo
        enddo
      elseif(iadir.eq.2) then
        do k = ks,ke
          do j = js,je
            q(j,k,nmv+1)=vnuinf
          enddo
        enddo
      endif

      return
      end

c*************************************************************
      subroutine turbc_out(q,u,v,ug,vg,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      real u(jd,kd),v(jd,kd),ug(jd,kd),vg(jd,kd)
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
          q(j,k,nmv+1) = vnuinf
          uu = (u(j,k)-ug(j,k)+u(j1,k)-ug(j1,k))*xx(j,k)
          uu = uu+(v(j,k)-vg(j,k)+v(j1,k)-vg(j1,k))*xy(j,k)
          uu = uu*iadd
          if(uu.lt.0.) q(j,k,nmv+1) = q(j1,k,nmv+1)
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
            q(j,k,nmv+1) = q(j1,k,nmv+1)
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
          q(j,k,nmv+1) = vnuinf
          uu=   (u(j,k)-ug(j,k)+u(j,k1)-ug(j,k1))*yx(j,k)
          uu=uu+(v(j,k)-vg(j,k)+v(j,k1)-vg(j,k1))*yy(j,k)
          uu=uu*iadd
          if(uu.lt.0.) q(j,k,nmv+1) = q(j,k1,nmv+1)
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
          enddo

        enddo

      endif

      return
      end

c*************************************************************
      subroutine lsolvej(a,b,c,s,jd,kd,js,je,ks,ke)
c*************************************************************
      use params_global
c*************************************************************
        
      integer jd,kd,js,je,ks,ke

      real a(jd,kd),b(jd,kd),c(jd,kd)
      real s(jd,kd)

c.. local variables

        integer j,k
        real bb

c***  first executable statement

      do k = ks,ke
        bb      = 1./b(js,k)
        s(js,k) = s(js,k)*bb
        c(js,k) = c(js,k)*bb
      enddo
 
      do  j = js+1,je
      do  k = ks,ke
         bb     = 1./(b(j,k) - a(j,k)*c(j-1,k))
         s(j,k) = (s(j,k) - a(j,k)*s(j-1,k))*bb
         c(j,k) = c(j,k)*bb
      enddo
      enddo

      do j = je-1,js,-1
      do k = ks,ke
         s(j,k) = s(j,k) - c(j,k)*s(j+1,k)
      enddo
      enddo

      return
      end

c*************************************************************
      subroutine lsolvek(a,b,c,s,jd,kd,js,je,ks,ke)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd,js,je,ks,ke
      real a(jd,kd),b(jd,kd),c(jd,kd)
      real s(jd,kd)

c.. local variables

      integer j,k
      real bb

c***  first executable statement

      do j = js,je
        bb       = 1./b(j,ks)
        s(j,ks)  = s(j,ks)*bb
        c(j,ks)  = c(j,ks)*bb
      enddo
 
      do  k = ks+1,ke
      do  j = js,je
        bb     = 1./(b(j,k) - a(j,k)*c(j,k-1))
        s(j,k) = (s(j,k) - a(j,k)*s(j,k-1))*bb
        c(j,k) = c(j,k)*bb
      enddo
      enddo

      do k = ke-1,ks,-1
      do j = js,je
         s(j,k) = s(j,k) - c(j,k)*s(j,k+1)
      enddo
      enddo
     
      return
      end



c***********************************************************************
      subroutine store_saproduction(prod_store,jd,kd,nsp,im)
c*************************************************************
      use params_global, only: nspec,num_grids
c*************************************************************
c     
c     write the production term out
c     Anand

c***********************************************************************
      implicit none
      integer jd,kd,nsp,im
      real :: prod_store(jd,kd)

      ! local variables
      integer j,k,logprod
      character*60 :: int_str,gfile

      do j=1,jd
         do k=1,kd
             if (isnan(prod_store(j,k))) then
                print *, j, k, prod_store(j,k)
                stop
             end if
         end do
      end do

      if(nspec.le.1) then
        gfile='production.dat'
        logprod = 22000+nspec
      else
        write(int_str,'(i7)')nsp
        gfile='production_TS'//trim(adjustl(int_str))//'.dat'
        logprod = 22000+nsp
      end if
      
      !am open(2211,file='production.dat',form='unformatted')
      if(im.eq.1) open(logprod,file=gfile,form='unformatted')
      !write(logprod) jd, kd
      write(logprod) ((prod_store(j,k),j=1,jd),k=1,kd)
      if(im.eq.num_grids) close(logprod)
      return

      end


c***********************************************************************
      subroutine store_sadestruction(dest_store,jd,kd,nsp,im)
c*************************************************************
      use params_global, only: nspec,num_grids
c*************************************************************
c     
c     write the destruction term out
c     Anand

c***********************************************************************
      implicit none
      integer jd,kd,nsp,im
      real :: dest_store(jd,kd)

      ! local variables
      integer j,k,logdest
      character*60 :: int_str,gfile

      do j=1,jd
         do k=1,kd
             if (isnan(dest_store(j,k))) then
                print *, j, k, dest_store(j,k)
                stop
             end if
         end do
      end do

      if(nspec.le.1) then
        gfile='destruction.dat'
        logdest = 23000+nspec
      else
        write(int_str,'(i7)')nsp
        gfile='destruction_TS'//trim(adjustl(int_str))//'.dat'
        logdest = 23000+nsp
      end if

      !am open(2211,file='destruction.dat',form='unformatted')
      if(im.eq.1) open(logdest,file=gfile,form='unformatted')
      !write(logdest) jd, kd
      write(logdest) ((dest_store(j,k),j=1,jd),k=1,kd)
      if(im.eq.num_grids) close(logdest)
      return

      end

c******************************************************************
      subroutine store_coeff_prod(x, y, coeff_prod, jd, kd,nsp,im)
c******************************************************************
      use params_global, only: num_grids
c******************************************************************
      implicit none
      integer jd,kd,j,k,nsp,im,logbeta
      real :: x(jd,kd), y(jd,kd), coeff_prod(jd,kd)
      
      logbeta = 99000+nsp
      if(nsp.eq.1) then
        if(im.eq.1) then
          open( logbeta, file = 'coeffprod.plt', form = 'formatted')
          write(logbeta, *) "TITLE = ""Coeff Prod"""
          write(logbeta, *) "VARIABLES = ""X"" ""Y"" ""beta"""
        endif
        write(logbeta,*) "ZONE"
        write(logbeta,*) "I = ",jd,", J = ",kd
        write(logbeta,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
        write(logbeta,'(F22.13)') ((x(j,k),j=1,jd),k=1,kd)
        write(logbeta,'(F22.13)') ((y(j,k),j=1,jd),k=1,kd)
        write(logbeta,*) ((coeff_prod(j,k),j=1,jd),k=1,kd)
        if(im.eq.num_grids) close(logbeta)
      endif !nsp
      end subroutine store_coeff_prod
c******************************************************************


