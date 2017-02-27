c***********************************************************************
      subroutine vmu_sa_ad(q,qb,s,sb,turmu,x,y,xv,yv,xx,xy,yx,yy,ug,vg,jd,kd,
     >                  tscale,iblank,im,isolve)
c  turbulent eddy viscosity. model is one equation spalart-
c  allmaras. ref{ aiaa 92-0439}.
c
c***********************************************************************
      use params_global
      use sa_transition_model
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer im,jd,kd,isolve
      real q(jd,kd,nq), qb(jd,kd,nq), s(jd,kd,nv), sb(jd,kd,nv)
      real turmu(jd,kd), tscale(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd) 
      integer iblank(jd,kd)

      ! local variables
      real,allocatable :: sjmat(:,:),sjmatb(:,:),bjmat(:,:)
      real,allocatable :: aj(:,:),bj(:,:),cj(:,:)
      real,allocatable :: ak(:,:),bk(:,:),ck(:,:)
      real,allocatable :: vmul(:,:),vmulb(:,:)
      real,allocatable :: vort(:,:),vortb(:,:),sn(:,:)
      real,allocatable :: u(:,:),ub(:,:),v(:,:),vb(:,:)
      real,allocatable :: s1(:,:),s2(:,:)
      real,allocatable :: trip(:,:)
      
      integer k,j,jloc,kloc,niter,ihuang,logv,nhuang
      real relfac,vnulim,dmaxtur,dl2norm
      real resl2,tscal,dtpseudo_turb
      real tempb,tempb0

      allocate(sjmat(jd,kd),sjmatb(jd,kd),bjmat(jd,kd))
      allocate(aj(jd,kd),bj(jd,kd),cj(jd,kd))
      allocate(ak(jd,kd),bk(jd,kd),ck(jd,kd))
      allocate(vmul(jd,kd),vmulb(jd,kd))
      allocate(vort(jd,kd),vortb(jd,kd),sn(jd,kd))
      allocate(u(jd,kd),ub(jd,kd),v(jd,kd),vb(jd,kd))
      allocate(s1(jd,kd),s2(jd,kd))
      allocate(trip(jd,kd))
c***  first executable statement

      relfac=1.0
      vnulim=1.0e-20
      nturiter=1
      nhuang=1

c...for compatibility with turbulence model

      do  k = 1,kd
        do  j = 1,jd
          u(j,k)  = q(j,k,2)/q(j,k,1)
          v(j,k)  = q(j,k,3)/q(j,k,1)
          sjmat(j,k) = 0.0
          vort(j,k) = 0.0
          sn(j,k) = 1.e10
          sjmatb(j,k) = sb(j,k,5)
          vortb(j,k) = 0.0
          vmulb(j,k) = 0.0
          ub(j,k) = 0.0
          vb(j,k) = 0.0

          aj(j,k) = 0.
          bj(j,k) = 0.
          cj(j,k) = 0.
          ak(j,k) = 0.
          bk(j,k) = 0.
          ck(j,k) = 0.
          bjmat(j,k) = 0.
        enddo
      enddo

c...laminar co-efficient of viscosity calculation
      
      call lamvis(q,vmul,jd,kd)

c...compute vorticity & distance function
       
      if(bodyflag(im)) then 
        call vortic(vort,u,v,xx,xy,yx,yy,jd,kd,jbeg,jend,kbeg,kend)
        call dist(sn,x,y,xv,yv,jd,kd)
      endif

      call c_turm(q,turmu,vmul,jd,kd)
      call ftrip(q,vort,sn,trip,u,v,jd,kd)

c   compute LHS      

      call convec_lhs(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
     &            ug,vg,aj,bj,cj,ak,bk,ck,jd,kd,jbeg,jend,kbeg,kend)

      call diffus_lhs(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
     &            aj,bj,cj,ak,bk,ck,jd,kd,jbeg,jend,kbeg,kend)

      call source_lhs(q,sn,u,v,sjmat,vmul,vort,xx,xy,yx,yy,
     &            aj,bj,cj,ak,bk,ck,jd,kd,jbeg,jend,kbeg,kend)


      do 30 niter=1,nturiter

        call RHSLHS_BQ(q, qb, s, sn, vort, vortb, u, ub, v, vb, 
     +                 vmul, vmulb, trip, sjmat, sjmatb, xx, xy, yx,
     +                 yy, ug, vg, aj, bj, cj, ak, bk, ck, jd, kd,
     +                 jbeg, jend, kbeg, kend)

        if(bodyflag(im)) then 
          call VORTIC_BQ(vort, vortb, u, ub, v, vb, xx, xy, yx, yy, jd,
     +                   kd, jbeg, jend, kbeg, kend)
        endif

        DO k=kd,1,-1
          DO j=jd,1,-1
            tempb = vb(j, k)/q(j, k, 1)
            qb(j, k, 3) = qb(j, k, 3) + tempb
            qb(j, k, 1) = qb(j, k, 1) - q(j, k, 3)*tempb/q(j, k, 1)
            vb(j, k) = 0.0
            tempb0 = ub(j, k)/q(j, k, 1)
            qb(j, k, 2) = qb(j, k, 2) + tempb0
            qb(j, k, 1) = qb(j, k, 1) - q(j, k, 2)*tempb0/q(j, k, 1)
            ub(j, k) = 0.0
          ENDDO
        ENDDO

        call LAMVIS_BQ(q, qb, vmul, vmulb, jd, kd)

        call TURBC_BQ(q, qb, u, v, ug, vg, xx, xy, yx, yy, jd, kd, 
     +                im)

        if (isolve==0) return
c   scale by time-step

        dl2norm=0.0

        do k = kbeg,kend
        do j = jbeg,jend

          dl2norm=dl2norm+qb(j,k,5)*qb(j,k,5)
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

          qb(j,k,5)  = qb(j,k,5)*tscal
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
              qb(j,k,5)  = bjmat(j,k)*(qb(j,k,5)-s1(j,k))
            enddo
            enddo
          endif

          call lsolvej(aj,bjmat,cj,qb(:,:,5),jd,kd,jbeg,jend,kbeg,kend)

          do k = kbeg,kend
          do j = jbeg,jend
            s1(j,k) = s1(j,k)+qb(j,k,5)
            qb(j,k,5)  = (s1(j,k)-s2(j,k))*bjmat(j,k)
          enddo
          enddo

          call lsolvek(ak,bjmat,ck,qb(:,:,5),jd,kd,jbeg,jend,kbeg,kend)

          do k = kbeg,kend
          do j = jbeg,jend
            s2(j,k)  = s2(j,k)+qb(j,k,5)
            qb(j,k,5)  = s2(j,k)
          enddo
          enddo

        enddo

        tscheck: if (.not.timespectral.or.itn.eq.1) then
          resl2=0.0
          do k = kbeg,kend
          do j = jbeg,jend
            qb(j,k,5) = relfac*qb(j,k,5)
            resl2=resl2+qb(j,k,5)*qb(j,k,5)
          enddo
          enddo
          resl2=sqrt(resl2/jd/kd)
          dl2norm=sqrt(dl2norm/jd/kd)

          if( mod(istep0,npnorm).eq.0 ) then
!$OMP ORDERED
            if(niter.eq.nturiter) write(11233+im,*)istep0,resl2,dl2norm
!$OMP END ORDERED
          endif
        endif tscheck
        
        do k = kbeg,kend
        do j = jbeg,jend
          sb(j,k,5) = (sb(j,k,5) + qb(j,k,5)*max(iblank(j,k),0))
        enddo
        enddo
        print*,'im,sum(qb): ',im,sum(qb(:,:,5))

30      continue

      return
      end

c***********************************************************************
      subroutine convec_lhs(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
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

c  j-direction jacobian terms

      do 30 k=ks-1,ke+1
      do 30 j=js-1,je+1

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

        cj(j,k)=cj(j,k)-fwd*(up(j,k)+um(j,k))
        bj(j,k)=  bj(j,k)+ up(j,k)-um(j,k)
        aj(j,k)=aj(j,k)+bck*(up(j,k)+um(j,k))
30    continue

c  k-direction jacobian terms

      do 40 k=ks-1,ke+1
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

        ck(j,k)=ck(j,k)-fwd*(vp(j,k)+vm(j,k))
        bk(j,k)=  bk(j,k)+ vp(j,k)-vm(j,k)
        ak(j,k)=ak(j,k)+bck*(vp(j,k)+vm(j,k))
40    continue

      return
      end

c***********************************************************************
      subroutine diffus_lhs(q,sn,u,v,sjmat,vmul,xx,xy,yx,yy,
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
        vnuh = 0.5*(q(j,k,5)+q(j+1,k,5))
        chp(j,k) = (1.+cb2)/sigma/rey*(vnulh+vnuh)
        dnuhp(j,k) = q(j+1,k,5)-q(j,k,5)

10    continue

c.....compute j direction fluxes

      do 20 k = ks-1,ke + 1
      do 20 j = js,je
      
        dxp = dmxh(j,k)*xx(j,k)   + dmyh(j,k)*xy(j,k)
        dxm = dmxh(j-1,k)*xx(j,k) + dmyh(j-1,k)*xy(j,k)
        
        c2 = cb2/sigma/rey
        c2 = c2*(vmul(j,k)/q(j,k,1)/q(j,k,nq)+q(j,k,5))

c.....enforce positivity (as suggested by overflow)

        dcp = dxp*(chp(j,k)-c2)
        dcm = dxm*(chp(j-1,k)-c2)
        ax = 0.0
        cx = 0.0
        if (k.ne.ke+1.and.k.ne.ks-1) then
            ax = max(dcm,0.0)
            cx = max(dcp,0.0)
        endif

c.....jacobian terms

        cj(j-1,k) = cj(j-1,k) - ax
        aj(j+1,k) = aj(j+1,k) - cx
        bj(j,k) = bj(j,k) + ax + cx

20    continue

c.....boundary terms in jacobians

      do 30 k = 1,kd
        aj(js-1,k) = 0.0
        aj(js,k) = 0.0
        cj(je,k) = 0.0
        cj(je+1,k) = 0.0
30    continue

c...compute k direction differences.

c.....compute half-point co-efficients

      do 40 k=ks-1,ke
      do 40 j=js-1,je+1

        dmxh(j,k) = 0.5*(yx(j,k)+yx(j,k+1))
        dmyh(j,k) = 0.5*(yy(j,k)+yy(j,k+1))
        vnulh = 0.5*(vmul(j,k)/q(j,k,1)/q(j,k,nq)+
     &               vmul(j,k+1)/q(j,k+1,1)/q(j,k+1,nq))
        vnuh = 0.5*(q(j,k,5)+q(j,k+1,5))
        chp(j,k) = (1.+cb2)/sigma/rey*(vnulh+vnuh)
        dnuhp(j,k) = q(j,k+1,5)-q(j,k,5)

40    continue

c.....compute j direction fluxes

      do 50 j = js-1,je+1
      do 50 k = ks,ke

        dyp = dmxh(j,k)*yx(j,k)   + dmyh(j,k)*yy(j,k)
        dym = dmxh(j,k-1)*yx(j,k) + dmyh(j,k-1)*yy(j,k)
        
        c2 = cb2/sigma/rey
        c2 = c2*(vmul(j,k)/q(j,k,1)/q(j,k,nq)+q(j,k,5))

c.....enforce positivity (as suggested by overflow)
       
        dcp = dyp*(chp(j,k)-c2)
        dcm = dym*(chp(j,k-1)-c2)
       
        ay = 0.0
        cy = 0.0
        if (j.ne.js-1.and.j.ne.je+1) then
          ay = max(dcm,0.0)
          cy = max(dcp,0.0)
        endif

c.....jacobian terms

        ck(j,k-1) = ck(j,k-1) - ay
        ak(j,k+1) = ak(j,k+1) - cy
        bk(j,k) = bk(j,k) + ay + cy

50    continue

c.....boundary terms in jacobians

      do 60 j=1,jd
        ak(j,ks-1) = 0.0
        ak(j,ks) = 0.0
        ck(j,ke) = 0.0
        ck(j,ke+1) = 0.0
60    continue

      return
      end

c***********************************************************************
      subroutine source_lhs(q,sn,u,v,sjmat,vmul,vort,xx,xy,yx,yy,
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
      
      real vmul(jd,kd),sn(jd,kd),vort(jd,kd)
      real u(jd,kd),v(jd,kd)
      
      include 'sadata.h'

      ! local variables
      real :: chp(jd,kd),dnuhp(jd,kd)
      real :: dmxh(jd,kd),dmyh(jd,kd)
      
      real rcv2,d2min,rmax,stilim,fturf
      real vnul,chi,d,fv1,fv2,fv3,ft2,dchi
      real dfv1,dfv2,dfv3,dft2,d2,stilda,r,g,fw
      real dstild,dr,dg,dfw,pro,des,prod,dest,dpro,ddes
      real tk1,tk2
      integer j,k,isour

c**   first executable statement

      rcv2=(1./5.)
      d2min=1.e-12
      rmax=10.0
      stilim=1.e-12
      fturf=0.0
      isour = 0

      do 10 k=ks,ke
      do 10 j=js,je
      vnul = vmul(j,k)/q(j,k,1)/q(j,k,nq)
      chi  = q(j,k,5)/vnul
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

      r  = q(j,k,5)/(d2*akt*akt*rey)/stilda
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

       prod = cb1*stilda*(1.-ft2)*q(j,k,5)
       if(obj_samod) then
          prod = prod*obj_coeff_prod(j,k)
       end if

       dest = (cw1*fw-cb1/(akt*akt)*ft2)*q(j,k,5)*q(j,k,5)/d2/rey

       dpro = pro*dstild/stilda - cb1*stilda*dft2
  
       ddes = (cw1*dfw-cb1/(akt*akt)*dft2)/d2/rey*vnul*chi
       ddes = ddes + des
       ddes = ddes*q(j,k,5)
       dpro = dpro*q(j,k,5)

c...modification for transition model

        if (itrans.eq.1) then
          prod = q(j,k,5)*prod
          dest = min(max(q(j,k,5),0.1),1.0)*dest
        endif

        tk1=max(des*q(j,k,5)-pro,0.0)
        tk2=max(ddes-dpro,0.0)
        bk(j,k) = bk(j,k)+tk1+tk2
c       bj(j,k) = bj(j,k)+tk1+tk2

10    continue

      return
      end

c*************************************************************
