c***********************************************************************
      subroutine rhsup( q,s,xx,xy,yx,yy,x,y,xv,yv,xold,yold,
     &     xole,yole,iblank,ugv,vgv,jd,kd,nei,im,bt)
c
c  muscl approach:
c  qt = 0                       1st-order  ; |irhsy| = 1
c  qt = 0.25
c    th = 0     upwind-biased   2nd-order  ; |irhsy| = 2
c    th = 1/3   upwind-biased   3rd-order  ; |irhsy| = 3
c
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,nei,im
      real q(jd,kd,nq), s(jd,kd,nv)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real ugv(jmax,kmax), vgv(jmax,kmax)
      real bt(jd,kd)
      real xold(jmax,kmax),yold(jmax,kmax)
      real xole(jmax,kmax),yole(jmax,kmax)
      integer iblank(jd,kd)

      ! local variables
      real :: ppp(jd,kd)
      real :: tj(mdim), xa(mdim), ya(mdim)
      real :: f(mdim,nmv),ql(mdim,nmv),qr(mdim,nmv)
      real :: fmin(nmv),fmax(nmv)
      integer :: iba(mdim)
      real :: a(mdim),b(mdim),c(mdim),fin(mdim),fout(mdim)
      real :: bbt(mdim)


      integer irhs,ilima,k,j,jv,kv,jj
      real th,qt,eps,epsj,epsk,rhoi,dx2,dy2,temp
      integer ibmin,ibmax

c**   first executable statement

      irhs   = iabs(irhsy)
      ilima = abs(ilim)
c     limter = 1
c     if( irhs .eq. 2 .and. limter.eq.1 ) limter = 2
      th     = real( irhs - 2 )/real( irhs )
      qt     = 0.25
      if( irhs .eq. 1 ) qt = 0.0
      eps    = 1.e-6
      epsj   =  (10./real(jmax))**3
      epsk   =  (10./real(kmax))**3
c
C$AD II-LOOP
      do k = 1, kd
C$AD II-LOOP
      do j = 1, jd
        ppp(j,k) = gm1*( q(j,k,4)-0.5*(q(j,k,2)**2+q(j,k,3)**2)
     <                                                      /q(j,k,1) )
      enddo
      enddo
c
c..xi fluxes
c..nonconservative variables
c     
C$AD II-LOOP
      do 13 k = kbeg,kend
        kv = k - nhalo + 1
c
C$AD II-LOOP
        do 15 j = 1,jd
          rhoi   = 1.0/q(j,k,1)
          f(j,1) = q(j,k,1)*q(j,k,nq)
          f(j,2) = q(j,k,2)*rhoi
          f(j,3) = q(j,k,3)*rhoi
          f(j,4) = ppp(j,k)*q(j,k,nq)
          iba(j) = abs(iblank(j,k))
          bbt(j) = 1.!bt(j,k)
   15   continue
         
c..at boundaries
        ibmin = 1
        ibmax = 1
c..limit
        if(ilima.eq.0)then
          call iflux(f,ql,qr,1,jd,jd-1,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
        endif
        if(ilima.eq.1)then
         call muscld_new(f,ql,qr,1,jd,jd-1,th,qt,epsj,fmin,fmax,ibmin,ibmax,iba)
c         call muscld(f,ql,qr,1,jd,jd-1,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
        if(ilima.eq.2)then
          call quad(f,ql,qr,1,jd,jd-1,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
        if(ilim.eq.4)then
          call weno(f,ql,qr,1,jd,jd-1,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
        endif
c
c..metric terms
c.do finite-volume like on finer mesh

        if(nei.eq.0) then
C$AD II-LOOP
         do j = jbeg-1,jend
          jv = j - nhalo + 1
          xa(j) = (yold(jv,kv)-yold(jv,kv-1))
          ya(j) =-(xold(jv,kv)-xold(jv,kv-1))
         enddo
        else
C$AD II-LOOP
         do j = jbeg-1,jend
          jv = j - nhalo + 1
          xa(j) = (yv(jv,kv)-yv(jv,kv-1))
          ya(j) =-(xv(jv,kv)-xv(jv,kv-1))
         enddo
        endif

	if(iunst.gt.0) then
C$AD II-LOOP
        do j = jbeg-1,jend
          jv = j - nhalo + 1
          dx1 = xv(jv,kv)-xold(jv,kv-1)
          dy1 = yv(jv,kv)-yold(jv,kv-1)
          dx2 = xold(jv,kv)-xv(jv,kv-1)
          dy2 = yold(jv,kv)-yv(jv,kv-1)
          tj(j) = -0.5*( dx1*dy2 - dx2*dy1 )/dt
        enddo
        if(ntac.eq.2.and.istep.gt.1) then
C$AD II-LOOP
         do j = jbeg-1,jend
          jv = j - nhalo + 1
          dx1 = xold(jv,kv)-xole(jv,kv-1)
          dy1 = yold(jv,kv)-yole(jv,kv-1)
          dx2 = xole(jv,kv)-xold(jv,kv-1)
          dy2 = yole(jv,kv)-yold(jv,kv-1)
          temp = -0.5*( dx1*dy2 - dx2*dy1 )/dt
          tj(j)= 1.5*tj(j)-0.5*temp
         enddo
        endif
        elseif (timespectral) then
C$AD II-LOOP
          do j = jbeg-1,jend
            jv = j - nhalo + 1
            tj(j) = - 0.5*( ugv(jv,kv) + ugv(jv,kv-1) )*xa(j)
     &              - 0.5*( vgv(jv,kv) + vgv(jv,kv-1) )*ya(j)
          enddo
 	else
 	do j = jbeg-1,jend
 	tj(j)=0.0
 	enddo
 	endif

c..compute the generalized numerical flux in roe!'s upwinding
c
        if (.not. iprecon) then
           call roeflx( f,ql,qr,xa,ya,tj,jbeg-1,jend )
        else
           call roetrklflx( f,ql,qr,xa,ya,tj,jbeg-1,jend,bbt)
        endif
c
C$AD II-LOOP
        do 12 j = jbeg,jend
          s(j,k,1) = s(j,k,1) - ( f(j,1) - f(j-1,1) )
          s(j,k,2) = s(j,k,2) - ( f(j,2) - f(j-1,2) )
          s(j,k,3) = s(j,k,3) - ( f(j,3) - f(j-1,3) )
          s(j,k,4) = s(j,k,4) - ( f(j,4) - f(j-1,4) )
   12   continue

   13 continue
c
c..eta fluxes
c
C$AD II-LOOP
      do 23 j = jbeg,jend
        jv = j - nhalo + 1
c
C$AD II-LOOP
        do 25 k = 1,kd
          rhoi   = 1.0/q(j,k,1)
          f(k,1) = q(j,k,1)*q(j,k,nq)
          f(k,2) = q(j,k,2)*rhoi
          f(k,3) = q(j,k,3)*rhoi
          f(k,4) = ppp(j,k)*q(j,k,nq)
          iba(k) = abs(iblank(j,k))
          bbt(k) = 1.!bt(j,k)
   25   continue

c..at boundaries
        ibmin = 1
        ibmax = 1
c..limit
        if(ilima.eq.0)then
          call iflux(f,ql,qr,1,kd,kd-1,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
        endif
        if(ilima.eq.1)then
         call muscld_new(f,ql,qr,1,kd,kd-1,th,qt,epsk,fmin,fmax,ibmin,ibmax,iba)
c          call muscld(f,ql,qr,1,kd,kd-1,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
        if(ilima.eq.2)then
          call quad(f,ql,qr,1,kd,kd-1,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
	if(ilim.eq.4)then
          call weno(f,ql,qr,1,kd,kd-1,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
	endif

c..metric terms
c.do finite-volume like on finer mesh

        if(nei.eq.0) then
C$AD II-LOOP
         do k = kbeg-1,kend
          kv = k - nhalo + 1
          xa(k) =-(yold(jv,kv)-yold(jv-1,kv))
          ya(k) = (xold(jv,kv)-xold(jv-1,kv))
         enddo
        else
C$AD II-LOOP
         do k = kbeg-1,kend
          kv = k - nhalo + 1
          xa(k) =-(yv(jv,kv)-yv(jv-1,kv))
          ya(k) = (xv(jv,kv)-xv(jv-1,kv))
         enddo
        endif

 	if(iunst.gt.0) then
C$AD II-LOOP
        do k = kbeg-1,kend
          kv = k - nhalo + 1
          dx1 = xv(jv,kv)-xold(jv-1,kv)
          dy1 = yv(jv,kv)-yold(jv-1,kv)
          dx2 = xold(jv,kv)-xv(jv-1,kv)
          dy2 = yold(jv,kv)-yv(jv-1,kv)
          tj(k)  = 0.5*( dx1*dy2 - dx2*dy1 )/dt
        enddo
        if(ntac.eq.2.and.istep.gt.1) then
C$AD II-LOOP
         do k = kbeg-1,kend
          kv = k - nhalo + 1
          dx1 = xold(jv,kv)-xole(jv-1,kv)
          dy1 = yold(jv,kv)-yole(jv-1,kv)
          dx2 = xole(jv,kv)-xold(jv-1,kv)
          dy2 = yole(jv,kv)-yold(jv-1,kv)
          temp  = 0.5*( dx1*dy2 - dx2*dy1 )/dt
          tj(k)= 1.5*tj(k)-0.5*temp
         enddo
        endif
        elseif (timespectral) then
C$AD II-LOOP
          do k = kbeg-1,kend
            kv = k - nhalo + 1
            tj(k) = - 0.5*( ugv(jv,kv) + ugv(jv-1,kv) )*xa(k)
     &              - 0.5*( vgv(jv,kv) + vgv(jv-1,kv) )*ya(k)
          enddo
 	else
C$AD II-LOOP
 	do k = kbeg-1,kend
 	tj(k)=0.
 	enddo
 	endif
c
c..compute the generalized numerical flux in roe!'s upwinding
c
        if (.not. iprecon) then
           call roeflx( f,ql,qr,xa,ya,tj,kbeg-1,kend )
        else
           call roetrklflx( f,ql,qr,xa,ya,tj,kbeg-1,kend,bbt)
        endif
c
C$AD II-LOOP
        do 22 k = kbeg,kend
          s(j,k,1) = s(j,k,1) - ( f(k,1) - f(k-1,1) )
          s(j,k,2) = s(j,k,2) - ( f(k,2) - f(k-1,2) )
          s(j,k,3) = s(j,k,3) - ( f(k,3) - f(k-1,3) )
          s(j,k,4) = s(j,k,4) - ( f(k,4) - f(k-1,4) )
   22   continue

   23 continue
c
      return
      end

c***********************************************************************
      subroutine roeflx( f,ql,qr,xa,ya,tj,is,ie )
c
c  compute the generalized numerical flux in roe!'s upwinding
c  by s.o.                          
c  mod by jdb to incorporate smoother entropy check
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer is,ie
      real f(mdim,nmv)
      real tj(mdim), xa(mdim), ya(mdim)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      real eps,rlft,ulft,vlft,plft
      real rlfti,rulft,rvlft,uvl,elft,hlft,clft
      real rrht,urht,vrht,prht
      real rrhti,rurht,rvrht,uvr,erht,hrht,crht
      real tklft,tomegalft,tkrht,tomegarht
      real rat,rati,rav,uav,vav,hav,uv,cav,tkav,tomegaav
      real aq1,aq2,aq3,aq4,aq5,aq6,ri1,ri2,ri3,rr2,rr,r0,r1,r2,r3
      real uu,c2,c2i,auu,aupc,aumc,uulft,uurht,upclft,upcrht
      real umclft,umcrht,dauu,dauus,daupc,daumc,daumcs,rcav,aquu
      real daupcs,c2ih,ruuav,b1,b2,b3,b4,b5,b6,b7,b8,b9,aj
      real plar,eplft,eprht,fssub

      integer i,i1

c***  first executable statement

      eps = 1.e-6
C$AD II-LOOP
      do 11 i = is,ie
c
       i1    = i + 1
       rlft = ql(i,1)
       ulft = ql(i,2)
       vlft = ql(i,3)
       plft = ql(i,4)
       rlfti = 1.0/rlft
       rulft = rlft*ulft
       rvlft = rlft*vlft
       uvl = 0.5*( ulft*ulft + vlft*vlft )
       elft = plft/gm1 + rlft*uvl
       hlft = ( elft + plft )*rlfti
       clft = sqrt( gm1*( hlft - uvl ) )
c
       rrht = qr(i1,1)
       urht = qr(i1,2)
       vrht = qr(i1,3)
       prht = qr(i1,4)
       rrhti = 1.0/rrht
       rurht = rrht*urht
       rvrht = rrht*vrht
       uvr = 0.5*( urht*urht + vrht*vrht )
       erht = prht/gm1 + rrht*uvr
       hrht = ( erht + prht )*rrhti
       crht = sqrt( gm1*( hrht - uvr ) )
c
       rat  = sqrt( rrht*rlfti )
       rati = 1.0/( rat + 1. )
       rav  =   rat*rlft
       uav  = ( rat*urht + ulft )*rati
       vav  = ( rat*vrht + vlft )*rati
       hav  = ( rat*hrht + hlft )*rati
       uv   = 0.5*( uav*uav + vav*vav )
       cav  = sqrt( gm1*( hav - uv ) )
c
       aq1  = rrht - rlft
       aq2  = urht - ulft
       aq3  = vrht - vlft
       aq4  = prht - plft
c
       ri1 = xa(i)
       ri2 = ya(i)
       ri3 = tj(i)
       rr2 = ri1*ri1 + ri2*ri2
       rr  = sqrt( rr2 )
       r0  = 1.0 / rr
       r1  = ri1*r0
       r2  = ri2*r0
       r3  = ri3*r0
c
       uu  = r1*uav + r2*vav + r3
       c2  = cav*cav
       c2i = 1.0/c2
c
       auu   = abs( uu    )
       aupc  = abs( uu+cav )
       aumc  = abs( uu-cav )
c     
       uulft = r1*ulft + r2*vlft + r3
       uurht = r1*urht + r2*vrht + r3
       upclft= uulft + clft
       upcrht= uurht + crht
       umclft= uulft - clft
       umcrht= uurht - crht
c
       dauu = 4.*(uurht-uulft)+eps
       dauus = amax1(dauu,0.0)
ccray       auu = cvmgt(auu**2/dauu+0.25*dauu,auu,auu.le.0.5*dauus)
       if( auu.le.0.5*dauus ) then
         auu = auu**2/dauu+0.25*dauu
       end if
       daupc = 4.*(upcrht-upclft)+eps
       daupcs = amax1(daupc,0.0)
ccray       aupc = cvmgt(aupc**2/daupc+0.25*daupc,aupc,aupc.le.0.5*daupcs)
       if( aupc.le.0.5*daupcs ) then
         aupc = aupc**2/daupc+0.25*daupc
       end if
       daumc = 4.*(umcrht-umclft)+eps
       daumcs = amax1(daumc,0.0)
ccray       aumc = cvmgt(aumc**2/daumc+0.25*daumc,aumc,aumc.le.0.5*daumcs)
       if( aumc.le.0.5*daumcs ) then
         aumc = aumc**2/daumc+0.25*daumc
       end if
c     
       rcav = rav*cav
       aquu = uurht - uulft
       c2ih = 0.5*c2i
       ruuav= auu*rav
       b1   = auu*( aq1 - c2i*aq4 )
       b2   = c2ih*aupc*( aq4 + rcav*aquu )
       b3   = c2ih*aumc*( aq4 - rcav*aquu )
       b4   = b1 + b2 + b3
       b5   = cav*( b2 - b3 )
       b6   = ruuav*( aq2 - r1*aquu )
       b7   = ruuav*( aq3 - r2*aquu )
c
       aq1 = b4
       aq2 = uav*b4 + r1*b5 + b6
       aq3 = vav*b4 + r2*b5 + b7
       aq4 = hav*b4 + ( uu-r3 )*b5 + uav*b6 + vav*b7 - c2*b1/gm1
c
       aj    = 0.5*rr
       plar  = plft + prht
       eplft = elft + plft
       eprht = erht + prht
       fssub = rr*r3
       fssub = 0.0
       f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
       f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 )
     <                                          -fssub*rinf*uinf
       f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
     <                                          -fssub*rinf*vinf
       f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
     <                                          -fssub*einf

c      fssub = rr*(r1*uinf + r2*vinf + r3)
c      f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
c      f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 ) 
c    <                            -fssub*rinf*uinf-rr*r1*pinf
c      f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
c    <                            -fssub*rinf*vinf-rr*r2*pinf
c      f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
c    <                            -fssub*(einf+pinf)+rr*r3*pinf
   11 continue
c
      return
      end


c***********************************************************************
      subroutine roetrklflx( f,ql,qr,xa,ya,tj,is,ie,b)
c
c  compute the generalized numerical flux in roe!'s upwinding
c  by s.o. and uses turkel preconditioning                         
c  mod by jdb to incorporate smoother entropy check
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer is,ie
      real f(mdim,nmv)
      real tj(mdim), xa(mdim), ya(mdim),b(mdim)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      real eps,rlft,ulft,vlft,plft
      real rlfti,rulft,rvlft,uvl,elft,hlft,clft
      real rrht,urht,vrht,prht
      real rrhti,rurht,rvrht,uvr,erht,hrht,crht
      real tklft,tomegalft,tkrht,tomegarht
      real rat,rati,rav,uav,vav,hav,uv,cav,tkav,tomegaav
      real aq1,aq2,aq3,aq4,aq5,aq6,ri1,ri2,ri3,rr2,rr,r0,r1,r2,r3
      real uumxt,uu,c2,c2i,auu,aupc,aumc,uulft,uurht,upclft,upcrht
      real umclft,umcrht,dauu,dauus,daupc,daumc,daumcs,rcav,aquu
      real daupcs,c2ih,ruuav,b1,b2,b3,b4,b5,b6,b7,b8,b9,aj
      real plar,eplft,eprht,fssub
      real R,S,X,bSq

      integer i,i1

c***  first executable statement

      eps = 1.e-6
C$AD II-LOOP
      do 11 i = is,ie
c
       i1    = i + 1
       rlft = ql(i,1)
       ulft = ql(i,2)
       vlft = ql(i,3)
       plft = ql(i,4)
       rlfti = 1.0/rlft
       rulft = rlft*ulft
       rvlft = rlft*vlft
       uvl = 0.5*( ulft*ulft + vlft*vlft )
       elft = plft/gm1 + rlft*uvl
       hlft = ( elft + plft )*rlfti
       clft = sqrt( gm1*( hlft - uvl ) )
c
       rrht = qr(i1,1)
       urht = qr(i1,2)
       vrht = qr(i1,3)
       prht = qr(i1,4)
       rrhti = 1.0/rrht
       rurht = rrht*urht
       rvrht = rrht*vrht
       uvr = 0.5*( urht*urht + vrht*vrht )
       erht = prht/gm1 + rrht*uvr
       hrht = ( erht + prht )*rrhti
       crht = sqrt( gm1*( hrht - uvr ) )
c
       rat  = sqrt( rrht*rlfti )
       rati = 1.0/( rat + 1. )
       rav  =   rat*rlft
       uav  = ( rat*urht + ulft )*rati
       vav  = ( rat*vrht + vlft )*rati
       hav  = ( rat*hrht + hlft )*rati
       uv   = 0.5*( uav*uav + vav*vav )
       cav  = sqrt( gm1*( hav - uv ) )
c
       aq1  = rrht - rlft
       aq2  = urht - ulft
       aq3  = vrht - vlft
       aq4  = prht - plft
c
       ri1 = xa(i)
       ri2 = ya(i)
       ri3 = tj(i)
       rr2 = ri1*ri1 + ri2*ri2
       rr  = sqrt( rr2 )
       r0  = 1.0 / rr
       r1  = ri1*r0
       r2  = ri2*r0
       r3  = ri3*r0
c
       uumxt  = r1*uav + r2*vav
       uu  = uumxt + r3
       c2  = cav*cav
       c2i = 1.0/c2
c
       bSq = Mp**2/(b(i)-Mp**2*(b(i)-1))

       X = sqrt( (1.-bSq)*uu*(1.-bSq)*uu+4.*bSq*c2 )
       auu   = abs( uu    )
       aupc  = 0.5*abs( (1.+bSq)*uu + X )
       aumc  = 0.5*abs( (1.+bSq)*uu - X )
c     
       uulft = r1*ulft + r2*vlft + r3
       uurht = r1*urht + r2*vrht + r3
       X = sqrt( (1.-bSq)*uulft*((1.-bSq)*uulft)+4.*bSq*clft*clft )
       upclft= 0.5*( (1.+bSq)*uulft + X )
       umclft= 0.5*( (1.+bSq)*uulft - X )
       X = sqrt( (1.-bSq)*uurht*((1.-bSq)*uurht)+4.*bSq*crht*crht )
       upcrht= 0.5*( (1.+bSq)*uurht + X )
       umcrht= 0.5*( (1.+bSq)*uurht - X )
c
       dauu = 4.*(uurht-uulft)+eps
       dauus = amax1(dauu,0.0)
ccray       auu = cvmgt(auu**2/dauu+0.25*dauu,auu,auu.le.0.5*dauus)
       if( auu.le.0.5*dauus ) then
         auu = auu**2/dauu+0.25*dauu
       end if
c
       daupc = 4.*(upcrht-upclft)+eps
       daupcs = amax1(daupc,0.0)
ccray       aupc = cvmgt(aupc**2/daupc+0.25*daupc,aupc,aupc.le.0.5*daupcs)
       if( aupc.le.0.5*daupcs ) then
         aupc = aupc**2/daupc+0.25*daupc
       end if
c
       daumc = 4.*(umcrht-umclft)+eps
       daumcs = amax1(daumc,0.0)
ccray       aumc = cvmgt(aumc**2/daumc+0.25*daumc,aumc,aumc.le.0.5*daumcs)
       if( aumc.le.0.5*daumcs ) then
         aumc = aumc**2/daumc+0.25*daumc
       end if
c     
       rcav = rav*cav
       aquu = uurht - uulft
       c2ih = 0.5*c2i
       ruuav= auu*rav
       X = sqrt( (1.-bSq)*uu*((1.-bSq)*uu)+4.*bSq*c2 )
       R     = 0.5*( (1.-bSq)*uu + X)
       S     = 0.5*( (1.-bSq)*uu - X )
       b1   = auu*( aq1 - c2i*aq4 )
       b2   = aupc*( aq4/R + rav*aquu )/X
       b3   = aumc*(-aq4/S - rav*aquu )/X
       b4   = b1 + b2 + b3
       b5   = ( R*b2 + S*b3 )
       b6   = ruuav*( aq2 - r1*aquu )
       b7   = ruuav*( aq3 - r2*aquu )
c
       aq1 = b4
       aq2 = uav*b4 + r1*b5 + b6
       aq3 = vav*b4 + r2*b5 + b7
       aq4 = hav*b4 + uumxt*b5 + uav*b6 + vav*b7 - c2*b1/gm1
c
       aj    = 0.5*rr
       plar  = plft + prht
       eplft = elft + plft
       eprht = erht + prht
       fssub = rr*r3
       fssub = 0.0
       f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
       f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 )
     <                                          -fssub*rinf*uinf
       f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
     <                                          -fssub*rinf*vinf
       f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
     <                                          -fssub*einf

c      fssub = rr*(r1*uinf + r2*vinf + r3)
c      f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
c      f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 ) 
c    <                            -fssub*rinf*uinf-rr*r1*pinf
c      f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
c    <                            -fssub*rinf*vinf-rr*r2*pinf
c      f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
c    <                            -fssub*(einf+pinf)+rr*r3*pinf
   11 continue
c
      return
      end

c*************************************************************************
      subroutine iflux(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  ibmin,ibmax,iba)
c
c 2nd order weno scheme with Van-Leer limiter
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************

      integer is,ie,im,ibmin,ibmax,iba(mdim)
      real th,qt,eps
      real f(mdim,nmv),fmin(nmv),fmax(nmv)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      integer i,n
      real s1,s2,slope,lim_vanleer

c..this is just 1st order upwind

      if(qt.eq.0.0)then
C$AD II-LOOP
        do n=1,nmv
C$AD II-LOOP
          do i=is,ie
            ql(i,n)=f(i,n)
            qr(i,n)=f(i,n)
          enddo
        enddo
      else
C$AD II-LOOP
        do n=1,nmv
C$AD II-LOOP
          do i = is+1,ie-1
            s1     = f(i,n)-f(i-1,n)
            s2     = f(i+1,n)-f(i,n)
            !slope  = lim_vanleer(s1,s2)
            slope  = (s1+s2)/2
            ql(i,n)= f(i,n) + 0.5*slope
            qr(i,n)= f(i,n) - 0.5*slope
          enddo

          if(ibmin.eq.2) then
            s1   = f(is,n)-fmin(n)
            s2   = f(is+1,n)-f(is,n)
            slope  = lim_vanleer(s1,s2)
            ql(is,n)= f(is,n) + 0.5*slope
            qr(is,n)= f(is,n) - 0.5*slope
          else
            ql(is,n)= f(is,n)
            qr(is,n)= f(is,n)
          endif

          if(ibmax.eq.2) then
            s1   = f(ie,n)-f(ie-1,n)
            s2   = fmax(n)-f(ie,n)
            slope  = lim_vanleer(s1,s2)
            ql(ie,n)= f(ie,n) + 0.5*slope
            qr(ie,n)= f(ie,n) - 0.5*slope
          else
            ql(ie,n)= f(ie,n)
            qr(ie,n)= f(ie,n)
          endif
        enddo
      endif

      return
      end

c*************************************************************************
      real function lim_minmod(x,y)
c*************************************************************************
      real :: x,y

      if(sign(1.0,x).ne.sign(1.0,y)) then
        lim_minmod = 0.0
      elseif(abs(x).lt.abs(y)) then
        lim_minmod = x
      else
        lim_minmod = y
      endif

      end function lim_minmod

c*************************************************************************
      real function lim_vanleer(x,y)
c*************************************************************************
      real :: x,y
      real :: a,b
      real :: lim_minmod
      a = 0.5*(x+y)
      b = 2.0*lim_minmod(x,y)

      lim_vanleer = lim_minmod(a,b)

      end function lim_vanleer

c*************************************************************************
      subroutine muscld(f,ql,qr,
     <                     is,imax,im,th,qt,eps,fmin,fmax,ibmin,ibmax)
c
c  muscl interpolation for higher order accuracy
c  differentiable limiter for 3rd order accuracy
c
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************

      real f(mdim,nmv)
      real ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)
      ! local variables
      real :: f2(mdim,nmv)
      integer is,imax,im,ibmin,ibmax,n,i
      real th,qt,eps,thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt

c***  first executable statement

      if(qt.eq.0.0)then
C$AD II-LOOP
        do 20 n=1,nmv
C$AD II-LOOP
        do 20 i=is,imax
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   20   continue
        return
      else
        thm = 1.0-th
        thp = 1.0+th
C$AD II-LOOP
        do 1 n=1,nmv
C$AD II-LOOP
          do 11 i=is,im
            f2(i,n) = f(i+1,n) - f(i,n)
 11       continue
C$AD II-LOOP
          do 12 i=is+1,im
           f2i    = f2(i,n)
           f2i1   = f2(i-1,n)
           a1     = 3.0*f2i*f2i1
           a2     = 2.0*(f2i-f2i1)**2 + a1
           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
           f3i    = ( a1 +epsf )/( a2 +epsf )
           f3qt   = qt*f3i
           ql(i,n)= f(i,n)+f3qt*( thm*f2i1 + thp*f2i )
           qr(i,n)= f(i,n)-f3qt*( thp*f2i1 + thm*f2i )
 12     continue
c..first-order at boundary
           qr(is  ,n)= f(is,n)
           ql(is  ,n)= f(is,n)
c..central at boundary?
c         qr(is,n) = 0.5*(f(is,n)+f(is+1,n))
c         ql(is,n) = 0.5*(f(is,n)+f(is+1,n))
          if(ibmin.eq.2) then
            f2i    = f(is+1,n)-f(is,n)
            f2i1   = f(is,n)-fmin(n)
            a1     = 3.0*f2i*f2i1
            a2     = 2.0*(f2i-f2i1)**2 + a1
            epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
c	    print*,'fix the epsilon - look at turns3d'
            f3i    = ( a1 +epsf )/( a2 +epsf )
            f3qt   = qt*f3i
            ql(is,n)= f(is,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(is,n)= f(is,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
c..first-order at boundary
           ql(imax,n)= f(imax,n)
           qr(imax,n)= f(imax,n)
c..central at boundary?
c         ql(imax,n) = 0.5*(f(imax,n)+f(im,n))
c         qr(imax,n) = 0.5*(f(imax,n)+f(im,n))
          if(ibmax.eq.2) then
            f2i    = fmax(n)-f(imax,n)
            f2i1   = f(imax,n)-f(im,n)
            a1     = 3.0*f2i*f2i1
            a2     = 2.0*(f2i-f2i1)**2 + a1
            epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
            f3i    = ( a1 +epsf )/( a2 +epsf )
            f3qt   = qt*f3i
            ql(imax,n)= f(imax,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(imax,n)= f(imax,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
    1   continue
      endif
c
      return
      end

c************************************************************************
      subroutine muscld_new(f,ql,qr,is,imax,im,th,qt,eps,
     <     fmin,fmax,ibmin,ibmax,iba)
c
c  muscl interpolation for higher order accuracy
c  differentiable limiter for 3rd order accuracy
c
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************
      real f(mdim,nmv)
      real ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)
      integer iba(mdim)
      
      ! local variables
      real :: f2(mdim,nmv)
      integer is,imax,im,ibmin,ibmax,n,i
      real thm,thp,f2i,f2i1,a1,a2,f3i,f3qt
      real th,qt,eps

      if(qt.eq.0.0)then
C$AD II-LOOP
        do 20 n=1,nmv
C$AD II-LOOP
        do 20 i=is,imax
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   20   continue
        return
      else
        thm = 1.0-th
        thp = 1.0+th
C$AD II-LOOP
        do 1 n=1,nmv
C$AD II-LOOP
          do 11 i=is,im
            f2(i,n) = f(i+1,n) - f(i,n)
 11       continue
C$AD II-LOOP
          do 12 i=is+1,im
           f2i    = f2(i,n)
           f2i1   = f2(i-1,n)
           a1     = 3.0*(f2i*f2i1+eps)
           a2     = 2.0*(f2i-f2i1)**2 + a1
c           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
           f3i    = a1/a2
           f3qt   = qt*f3i
     >           *iba(i)*iba(i+1)*iba(i-1)
           ql(i,n)= f(i,n)+f3qt*( thm*f2i1 + thp*f2i )
           qr(i,n)= f(i,n)-f3qt*( thp*f2i1 + thm*f2i )
 12     continue
c..first-order at boundary
           qr(is  ,n)= f(is,n)
           ql(is  ,n)= f(is,n)
c..central at boundary?
c         qr(is,n) = 0.5*(f(is,n)+f(is+1,n))
c         ql(is,n) = 0.5*(f(is,n)+f(is+1,n))
          if(ibmin.eq.2) then
            f2i    = f(is+1,n)-f(is,n)
            f2i1   = f(is,n)-fmin(n)
            a1     = 3.0*(f2i*f2i1+eps)
            a2     = 2.0*(f2i-f2i1)**2 + a1
c           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
            f3i    = a1/a2
            f3qt   = qt*f3i
            ql(is,n)= f(is,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(is,n)= f(is,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
c..first-order at boundary
           ql(imax,n)= f(imax,n)
           qr(imax,n)= f(imax,n)
c..central at boundary?
c         ql(imax,n) = 0.5*(f(imax,n)+f(im,n))
c         qr(imax,n) = 0.5*(f(imax,n)+f(im,n))
          if(ibmax.eq.2) then
            f2i    = fmax(n)-f(imax,n)
            f2i1   = f(imax,n)-f(im,n)
            a1     = 3.0*(f2i*f2i1+eps)
            a2     = 2.0*(f2i-f2i1)**2 + a1
c           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
            f3i    = a1/a2
            f3qt   = qt*f3i
            ql(imax,n)= f(imax,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(imax,n)= f(imax,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
    1   continue
      endif
c
      return
      end


c************************************************************************
      subroutine weno3(f,ql,qr,is,imax,im,th,qt,eps,
     <     fmin,fmax,ibmin,ibmax,iba)
c
c  3rd order WENO interpolation for higher order accuracy
c
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************
      real f(mdim,nmv)
      real ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)
      integer iba(mdim)
      
      ! local variables
      real :: f2(mdim,nmv)
      integer is,imax,im,ibmin,ibmax,n,i
      real f2i,f2i1,a1,a2,w01,w02,w11,w12,f3
      real th,qt,eps,epsw


      epsw = 1e-6

      if(qt.eq.0.0)then
C$AD II-LOOP
        do 20 n=1,nmv
C$AD II-LOOP
        do 20 i=is,imax
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   20   continue
        return
      else
C$AD II-LOOP
        do 1 n=1,nmv
C$AD II-LOOP
          do 11 i=is,im
            f2(i,n) = f(i+1,n) - f(i,n)
 11       continue
C$AD II-LOOP
          do 12 i=is+1,im
           f2i    = f2(i,n)
           f2i1   = f2(i-1,n)
           a1     = 1./(f2i1**2+epsw)
           a2     = 1./(f2i**2+epsw)
           w01     = a1/(a1+2*a2) 
           w02     = 1. - w01 
           w11     = 2*a1/(2*a1+a2) 
           w12     = 1. - w11 
           f3   = 0.5*iba(i)*iba(i+1)*iba(i-1)
           ql(i,n)= f(i,n)+f3*( w01*f2i1 + w02*f2i )
           qr(i,n)= f(i,n)-f3*( w11*f2i1 + w12*f2i )
 12     continue
c..first-order at boundary
           qr(is  ,n)= f(is,n)
           ql(is  ,n)= f(is,n)
c..central at boundary?
          if(ibmin.eq.2) then
            f2i    = f(is+1,n)-f(is,n)
            f2i1   = f(is,n)-fmin(n)
            a1     = 1./(f2i1**2+epsw)
            a2     = 1./(f2i**2+epsw)
            w01     = a1/(a1+2*a2) 
            w02     = 1. - w01 
            w11     = 2*a1/(2*a1+a2) 
            w12     = 1. - w11 
            f3   = 0.5
            ql(is,n)= f(is,n)+f3*( w01*f2i1 + w02*f2i )
            qr(is,n)= f(is,n)-f3*( w11*f2i1 + w12*f2i )
          endif
c..first-order at boundary
           ql(imax,n)= f(imax,n)
           qr(imax,n)= f(imax,n)
c..central at boundary?
          if(ibmax.eq.2) then
            f2i    = fmax(n)-f(imax,n)
            f2i1   = f(imax,n)-f(im,n)
            a1     = 1./(f2i1**2+epsw)
            a2     = 1./(f2i**2+epsw)
            w01     = a1/(a1+2*a2) 
            w02     = 1. - w01 
            w11     = 2*a1/(2*a1+a2) 
            w12     = 1. - w11 
            f3   = 0.5
            ql(imax,n)= f(imax,n)+f3*( w01*f2i1 + w02*f2i )
            qr(imax,n)= f(imax,n)-f3*( w11*f2i1 + w12*f2i )
          endif
    1   continue
      endif
c
      return
      end

c**************************************************************************
      subroutine weno(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  ibmin,ibmax,iba)
c
c 5th order weno scheme
c note: qr(is) and ql(ie) are never used.
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************

      integer is,ie,im,ibmin,ibmax,iba(mdim)
      real th,qt,eps,at1
      real f(mdim,nmv),fmin(nmv),fmax(nmv)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      real :: f0(mdim),f1(mdim),f2(mdim),f2m(mdim)
      real :: slope(mdim)

      real dre,duj,vjp1,vjm1,ejp1,ammd,fl,fr
      integer i,n
      real f0bmin,f1bmin,f2bmin,f0bmax,f1bmax,f2bmax,at,s1,ati1,t1
      real weno5,blankv
c*************************************************************************
      ammd(fl,fr) = 0.5*(sign(1.,fl)+sign(1.,fr))*amin1(abs(fl),abs(fr))
c*************************************************************************

c..this is just 1st order upwind
      if(qt.eq.0.0)then
C$AD II-LOOP
        do n=1,nmv
C$AD II-LOOP
          do i=is,ie
            ql(i,n)=f(i,n)
            qr(i,n)=f(i,n)
          enddo
        enddo
        return
      else
c
c 5th order weno scheme follows
c
C$AD II-LOOP
        do 1 n=1,nmv
c
c..let's load up a few difference arrays
c
C$AD II-LOOP
          do i=is,ie
            f0(i) = f(i,n)
          enddo
c..1st difference at i+1/2
C$AD II-LOOP
          do i=is,ie-1
            f1(i) = f0(i+1) - f0(i)
          enddo
c..2nd difference at i
C$AD II-LOOP
          do i=is+1,ie-1
            f2(i)  = f1(i) - f1(i-1)
          enddo
c..extrapolate at the boundaries
          f2(is) = 2.*f2(is+1)-f2(is+2)
          f2(ie) = 2.*f2(ie-1)-f2(ie-2)
c..modify at boundaries, if needed
          if(ibmin.eq.2) then
            f0bmin = fmin(n)
            f1bmin = f0(is) - f0bmin
            f2(is) = f1(is) - f1bmin
            f2bmin = 2.*f2(is)-f2(is+1)
          endif
          if(ibmax.eq.2) then
            f0bmax = fmax(n)
            f1bmax = f0bmax - f0(ie)
            f2(ie) = f1bmax - f1(ie-1)
            f2bmax = 2.*f2(ie)-f2(ie-1)
          endif
c..limit 2nd difference to i+1/2
C$AD II-LOOP
          do i = is,ie-1
            f2m(i) = ammd(f2(i),f2(i+1))
          enddo
c
c..now combine everything to get ql and qr
c
c
c..first-order at boundary?
c
          at    = f1(is) - 0.5*f2m(is)
          slope(i) = ammd(at,2.*f1(is))
          ql(is,n) = f0(is) + 0.5*slope(i)
c
          if(ibmin.eq.2) then
            i = is
            s1    = ammd(f1(i),f1bmin)
            at    = f1(i)  - 0.5*f2m(i)
            ati1  = f1bmin + 0.5*ammd(f2(i),f2bmin)
            t1    = ammd(at,ati1)
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            ql(i,n) = f0(i) + 0.5*slope(i)
          endif
c
          at1   = f1(ie-1) + 0.5*f2m(ie-1)
          slope(i) = ammd(at1,2.*f1(ie-1))
          qr(ie,n) = f0(ie) - 0.5*slope(i)
c
          if(ibmax.eq.2) then
            i = ie
            s1    = ammd(f1bmax,f1(i-1))
            at    = f1bmax  - 0.5*ammd(f2(i-1),f2bmax)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
            t1    = ammd(at,ati1)
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            qr(i,n) = f0(i) - 0.5*slope(i)
          endif
c
c..sonic-a near the boundary?
c
C$AD II-LOOP
          do i=is+1,ie-1,ie-is-2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope(i)
            qr(i,n) = f0(i) - 0.5*slope(i)
          enddo
c
c..suresh at interior
c
C$AD II-LOOP
          do i=is+2,ie-2
          ql(i,n) = weno5(f0(i-2),f0(i-1),f0(i),f0(i+1),f0(i+2))
          qr(i,n) = weno5(f0(i+2),f0(i+1),f0(i),f0(i-1),f0(i-2))
          blankv  = iba(i)*iba(i-1)*iba(i+1)*iba(i-2)*iba(i+2)
	  ql(i,n) = ql(i,n)*blankv+(1.-blankv)*f0(i)
	  qr(i,n) = qr(i,n)*blankv+(1.-blankv)*f0(i)
          enddo
c
c..sonic-a near the boundary?
c
          i=is+2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope(i)
c
c..sonic-a near the boundary?
c
          i=ie-1
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            qr(i,n) = f0(i) - 0.5*slope(i)
    1   continue
      endif
c
      return
      end
c*************************************************************************
      function weno5(a,b,c,d,e)
c
c
c*************************************************************************
      implicit none
c*************************************************************************
      real weno5
      real a,b,c,d,e
      real b1,b2,epsw,djm1,ejm1,dj,ej,djp1,ejp1
      real dis0,dis1,dis2,q30,q31,q32,d01,d02,a1ba0,a2ba0
      real w0,w1,w2
      
      b1 = 13./12.
      b2 = 1./6.
      epsw = 1.e-6
      djm1 = a-2.*b+c
      ejm1 = a-4.*b+3.*c
      dj   = b-2.*c+d
      ej   = b-d
      djp1 = c-2.*d+e
      ejp1 = 3.*c-4.*d+e
      dis0 = b1*djm1*djm1+0.25*ejm1*ejm1+epsw
      dis1 = b1*dj*dj+0.25*ej*ej+epsw
      dis2 = b1*djp1*djp1+0.25*ejp1*ejp1+epsw
      q30 = 2.*a-7.*b+11.*c
      q31 = -b+5.*c+2.*d
      q32 = 2.*c+5.*d-e
      d01 = dis0/dis1
      d02 = dis0/dis2
      a1ba0 = 6.*d01
      a2ba0 = 3.*d02
      w0 = 1./(1.+a1ba0+a2ba0)
      w1 = a1ba0*w0
      w2 = 1.-w0-w1
      weno5 = b2*(w0*q30+w1*q31+w2*q32)
      return
      end

c*************************************************************************
      subroutine quad(f,ql,qr,
     <                    is,ie,im,th,qt,eps,fmin,fmax,ibmin,ibmax)
c
c..piecewise quadratic reconstruction with 6th order compact evaluation of
c  nodal derivatives and new monotonicity-preserving constraint.
c
c**************************************************************************
      use params_global
c**************************************************************************
      implicit none
c**************************************************************************
      integer is,ie,im,ibmin,ibmax
      real th,qt,eps
      real f(mdim,nmv),ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)

      ! local variables
      real :: d1(mdim),d2(mdim)
      real :: con1(mdim),con2(mdim),con3(mdim)
      real :: s1(mdim),sf(mdim),sc(mdim),sb(mdim)

      real amedian,ammd
      real b,c,d,fl,fr
      integer n,i,is1,is2,ie2,ie1
      real f1i,f1i1,f1i2,f1i3,f2i1,f2i2,f2i3,f2i,f2mi
      real f1bmin,f2bmin,f2mi1,f2mi2,at,sl,ati1,t1,slopes2
      real f1i4,f1bmax,f2bmax,f2mi3,at1,slopee2,sm,cm,si,ci
      real sr,cr,cl
c**************************************************************************
      ammd(fl,fr) = 0.5*(sign(1.,fl)+sign(1.,fr))*amin1(abs(fl),abs(fr))
      amedian(b,c,d)=b+0.5*(sign(1.,c-b)+sign(1.,d-b))*
     &     amin1(abs(c-b),abs(d-b))
c**************************************************************************

c***  first executable statement

      if(qt.eq.0.)then
C$AD II-LOOP
        do 10 n=1,nmv
C$AD II-LOOP
        do 10 i=is,ie
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   10   continue
        return
      else
        is1=is+1
        is2=is+2
        ie2=ie-2
        ie1=ie-1
c
C$AD II-LOOP
        do 1 n=1,nmv
c
c..slope:
c
        i=is
        f1i=f(i+1,n)-f(i,n)
        f1i1=f(i+2,n)-f(i+1,n)
        f1i2=f(i+3,n)-f(i+2,n)
        f1i3=f(i+4,n)-f(i+3,n)
        f2i1=f1i1-f1i
        f2i2=f1i2-f1i1
        f2i3=f1i3-f1i2
        f2i=2.*f2i1-f2i2
        if(ibmin.eq.2) then
          f1bmin=f(i,n)-fmin(n)
          f2i=f1i-f1bmin
          f2bmin=2.*f2i-f2i1
        endif
        f2mi=ammd(f2i,f2i1)
        f2mi1=ammd(f2i1,f2i2)
        f2mi2=ammd(f2i2,f2i3)
        at=f1i-0.5*f2mi
        d1(i)=ammd(at,2.*f1i)
        if(ibmin.eq.2) then
          sl=ammd(f1i,f1bmin)
          at=f1i-0.5*f2mi
          ati1=f1bmin+0.5*ammd(f2i,f2bmin)
          t1=ammd(at,ati1)
          d1(i)=sign(1.,t1)*amin1(0.5*abs(at+ati1),
     &          amax1(2.*abs(sl),abs(t1)))
        endif
c
c..sonica scheme at the boundary:
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=is1
        at=f1i1-0.5*f2mi1
        ati1=f1i+0.5*f2mi
        sl=ammd(f1i1,f1i)
        t1=ammd(at,ati1)
        d1(i)=sign(1.,t1)*
     &        amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=is2
        at=f1i2-0.5*f2mi2
        ati1=f1i1+0.5*f2mi1
        sl=ammd(f1i2,f1i1)
        t1=ammd(at,ati1)
        slopes2=sign(1.,t1)*
     &          amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        con1(i)=0.25
        con2(i)=1.
        con3(i)=0.75*(f(i+1,n)-f(i-1,n))
        con3(i)=con3(i)-con1(i)*d1(i-1)
c
C$AD II-LOOP
        do i=is+3,ie-3
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(f(i+2,n)+28.*f(i+1,n)-28.*f(i-1,n)-f(i-2,n))/36.
        enddo
c
        i=ie
        f1i1=f(i,n)-f(i-1,n)
        f1i2=f(i-1,n)-f(i-2,n)
        f1i3=f(i-2,n)-f(i-3,n)
        f1i4=f(i-3,n)-f(i-4,n)
        f2i1=f1i1-f1i2
        f2i2=f1i2-f1i3
        f2i3=f1i3-f1i4
        f2i=2.*f2i1-f2i2
        if(ibmax.eq.2) then
          f1bmax=fmax(n)-f(i,n)
          f2i=f1bmax-f1i1
          f2bmax=2.*f2i-f2i1
        endif
        f2mi1=ammd(f2i1,f2i)
        f2mi2=ammd(f2i2,f2i1)
        f2mi3=ammd(f2i3,f2i2)
        at1=f1i1+0.5*f2mi1
        d1(i)=ammd(at1,2.*f1i1)
        if(ibmax.eq.2) then
          sl=ammd(f1bmax,f1i1)
          at=f1bmax-0.5*ammd(f2i1,f2bmax)
          ati1=f1i1+0.5*f2mi1
          t1=ammd(at,ati1)
          d1(i)=sign(1.,t1)*amin1(0.5*abs(at+ati1),
     &          amax1(2.*abs(sl),abs(t1)))
        endif
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=ie1
        at=f1i1-0.5*f2mi1
        ati1=f1i2+0.5*f2mi2
        sl=ammd(f1i1,f1i2)
        t1=ammd(at,ati1)
        d1(i)=sign(1.,t1)*
     &        amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=ie2
        at=f1i2-0.5*f2mi2
        ati1=f1i3+0.5*f2mi3
        sl=ammd(f1i2,f1i3)
        t1=ammd(at,ati1)
        slopee2=sign(1.,t1)*
     &          amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        con1(i)=0.25
        con2(i)=1.
        con3(i)=0.75*(f(i+1,n)-f(i-1,n))
        con3(i)=con3(i)-con1(i)*d1(i+1)
        call tridag(con1,con2,con1,con3,d1,is2,ie2)
c
c..curvature:
c
        i=is
        d2(i)=0.
        i=is1
        d2(i)=0.
        i=is2
        con1(i)=0.1
        con2(i)=1.
        con3(i)=1.2*(f(i+1,n)-2.*f(i,n)+f(i-1,n))
        con3(i)=con3(i)-con1(i)*d2(i-1)
C$AD II-LOOP
        do i=is+3,ie-3
          con1(i)=2./11.
          con2(i)=1.
          con3(i)=(3.*f(i+2,n)+48.*f(i+1,n)-102.*f(i,n)+48.*f(i-1,n)
     &             +3.*f(i-2,n))/44.
        enddo
        i=ie
        d2(i)=0.
        i=ie1
        d2(i)=0.
        i=ie2
        con1(i)=0.1
        con2(i)=1.
        con3(i)=1.2*(f(i+1,n)-2.*f(i,n)+f(i-1,n))
        con3(i)=con3(i)-con1(i)*d2(i+1)
        call tridag(con1,con2,con1,con3,d2,is2,ie2)
c
c..correct the sign of the slope and curvature:
c
C$AD II-LOOP
        do i=is1,ie
          s1(i)=f(i,n)-f(i-1,n)
        enddo
C$AD II-LOOP
        do i=is,ie2
          sf(i)=(-f(i+2,n)+4.*f(i+1,n)-3.*f(i,n))/2.
        enddo
C$AD II-LOOP
        do i=is1,ie1
          sc(i)=(f(i+1,n)-f(i-1,n))/2.
        enddo
C$AD II-LOOP
        do i=is2,ie
          sb(i)=(f(i-2,n)-4.*f(i-1,n)+3.*f(i,n))/2.
        enddo
c
C$AD II-LOOP
        do i=is2,ie2
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=f(i+2,n)-2.*f(i+1,n)+f(i,n)
          if(sm.eq.sc(i)) cm=f(i+1,n)-2.*f(i,n)+f(i-1,n)
          if(sm.eq.sb(i)) cm=f(i,n)-2.*f(i-1,n)+f(i-2,n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=f(i+1,n)-2.*f(i,n)+f(i-1,n)
        enddo
c
C$AD II-LOOP
        do i=is2,ie2
          si=d1(i)
          ci=d2(i)
c..impose monotonicity-preserving constraints:
c..method 1:
c
c         if(abs(d1(i)).gt.2.0*abs(s1(i+1)).or.abs(d1(i)).gt.
c    &       2.0*abs(s1(i))) then
c           if(abs(d1(i)).le.abs(d1(i+1))) then
c             sr=d1(i)
c             cr=s1(i+1)-d1(i)
c           else
c             sr=2.*s1(i+1)-d1(i+1)
c             cr=-2.*s1(i+1)+2.*d1(i+1)
c           endif
c           if(abs(d1(i)).le.abs(d1(i-1))) then
c             sl=d1(i)
c             cl=-s1(i)+d1(i)
c           else
c             sl=2.*s1(i)-d1(i-1)
c             cl=2.*s1(i)-2.*d1(i-1)
c           endif
c           si=ammd(sl,sr)
c           ci=ammd(cl,cr)
c         endif
c
c..method 2:
c
          if(abs(d1(i)).gt.2.0*abs(s1(i+1)).and.abs(d1(i)).gt.
     &       2.0*abs(s1(i))) then
            if(abs(d1(i)).le.abs(d1(i+1))) then
              sr=d1(i)
              cr=2.*(s1(i+1)-d1(i))
            else
              sr=2.*s1(i+1)-d1(i+1)
              cr=2.*(-s1(i+1)+d1(i+1))
            endif
            if(abs(d1(i)).le.abs(d1(i-1))) then
              sl=d1(i)
              cl=2.*(-s1(i)+d1(i))
            else
              sl=2.*s1(i)-d1(i-1)
              cl=2.*(s1(i)-d1(i-1))
            endif
            si=ammd(sl,sr)
            ci=ammd(cl,cr)
          else
            if(abs(d1(i)).gt.2.0*abs(s1(i+1))) then
              if(abs(d1(i)).le.abs(d1(i+1))) then
                sr=d1(i)
                cr=2.*(s1(i+1)-d1(i))
              else
                sr=2.*s1(i+1)-d1(i+1)
                cr=2.*(-s1(i+1)+d1(i+1))
              endif
              si=ammd(d1(i),sr)
              ci=ammd(d2(i),cr)
            endif
            if(abs(d1(i)).gt.2.0*abs(s1(i))) then
              if(abs(d1(i)).le.abs(d1(i-1))) then
                sl=d1(i)
                cl=2.*(-s1(i)+d1(i))
              else
                sl=2.*s1(i)-d1(i-1)
                cl=2.*(s1(i)-d1(i-1))
              endif
              si=ammd(sl,d1(i))
              ci=ammd(cl,d2(i))
            endif
          endif
c..piecewise quadratic reconstruction:
          ql(i,n)=f(i,n)+0.5*si+ci/12.
          qr(i,n)=f(i,n)-0.5*si+ci/12.
        enddo
c
c..consistency with the sonica scheme:
c
        qr(is2,n)=f(is2,n)-0.5*slopes2
        ql(ie2,n)=f(ie2,n)+0.5*slopee2
 1      continue
      endif
c
      return
      end

c***********************************************************************
      subroutine tridag(a,b,c,f,z,ni,nl)

c     tri diagonal matrix inversion
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,ni,nl
      parameter (jd=1001)
      real a(jd),b(jd),c(jd),f(jd),z(jd)

      ! local variables
      real :: w(jd),g(jd)

      integer nipl,j,nd,j1
      real d,rd

c
      w(ni)=c(ni)/b(ni)
      g(ni)=f(ni)/b(ni)
      nipl=ni+1
      do 10 j=nipl,nl
      d=b(j)-a(j)*w(j-1)
      rd=1.0/d
      w(j)=c(j)*rd
      g(j)=(f(j)-a(j)*g(j-1))*rd
 10   continue
      z(nl)=g(nl)
      nd=nl-ni
      do 20 j1=1,nd
      j=nl-j1
      z(j)=g(j)-w(j)*z(j+1)
 20   continue
      return
      end

c***********************************************************************
