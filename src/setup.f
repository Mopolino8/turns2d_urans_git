c***********************************************************************
      subroutine determine_size(jgmx,kgmx,nspec,nq,nv,ipointer,nmesh,
     &     igrd,igrdv,igrdb,iqdim,isdim,mdim,nadd)
c
c     Determines the sizes of all arrays for allocation purposes
c     Also sets pointers
c
c     igrd  -> coordinates,metrics  -> ipointer(:,1)
c     iqdim -> q-variables          -> ipointer(:,2)
c     isdim -> flux                 -> ipointer(:,3)
c     igrdv -> cell vertex coord    -> ipointer(:,4)
c     igrdb -> refined mesh         -> ipointer(:,5)
c***********************************************************************
      implicit none
c***********************************************************************

      integer nq,nv,nspec,nmesh,igrd,igrdv,igrdb,iqdim,isdim,mdim,nadd
      integer jgmx(nmesh),kgmx(nmesh)
      integer ipointer(nmesh,5)

      ! local variables
      
      integer i

c**   first executable statement
      do i=1,5
         ipointer(1,i)=1
      enddo

      do i=2,nmesh
         ipointer(i,1)=ipointer(i-1,1) + jgmx(i-1)*kgmx(i-1)*nspec
         ipointer(i,2)=ipointer(i-1,2) + jgmx(i-1)*kgmx(i-1)*nq*nspec
         ipointer(i,3)=ipointer(i-1,3) + jgmx(i-1)*kgmx(i-1)*nv*nspec
         ipointer(i,4)=ipointer(i-1,4) + (jgmx(i-1)-nadd)*(kgmx(i-1)-nadd)*nspec
         ipointer(i,5)=ipointer(i-1,5) + 4*jgmx(i-1)*kgmx(i-1)*nspec
      enddo

      i=nmesh+1
      igrd  = ipointer(i-1,1) + jgmx(i-1)*kgmx(i-1)*nspec
      iqdim = ipointer(i-1,2) + jgmx(i-1)*kgmx(i-1)*nq*nspec
      isdim = ipointer(i-1,3) + jgmx(i-1)*kgmx(i-1)*nv*nspec
      igrdv = ipointer(i-1,4) + (jgmx(i-1)-nadd)*(kgmx(i-1)-nadd)*nspec
      igrdb = ipointer(i-1,5) + 4*jgmx(i-1)*kgmx(i-1)*nspec

      mdim=0
      do i=1,nmesh
         if (mdim.lt.jgmx(i)) mdim=jgmx(i)
         if (mdim.lt.kgmx(i)) mdim=kgmx(i)
      enddo
      
      return
      end

c***********************************************************************
      subroutine set_pointers_globals(im,ipointer,ig,igq,
     &     igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
c
c     Sets the pointers of the current mesh to the global size parameters
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer im,ig,igq,igs,igv,igb,jd,kd,nmesh
      integer ipointer(nmesh,5),jgmx(nmesh),kgmx(nmesh)
      
c***  first executable statement

      ig=ipointer(im,1)
      igq=ipointer(im,2)
      igs=ipointer(im,3)
      igv=ipointer(im,4)
      igb=ipointer(im,5)

      jd=jgmx(im)
      kd=kgmx(im)

      jmax = jd - nadd
      kmax = kd - nadd

      jbeg = nhalo + 1
      kbeg = nhalo + 1
      jend = jd-nhalo
      kend = kd-nhalo

      jtail1 = jtail(im)
      jtail1 = jtail1 + nhalo
C
      if (half.eq.0) then
        jle = (jd+1)/2
        jtail2 = jd - jtail1 + 1
      else
        jle = jd - 1
        jtail2 = jd -  1
      endif

      !amm rotf = irot_dir_all(im)*rf !asitav
      if(rf.ne.0) then
        rotf = irot_dir_all(im)*rf 
      elseif(rf_tef.ne.0 .and. motiontype.eq.3) then
        rotf = irot_dir_all(im)*rf_tef
      end if

      xac = xac_all(im)
      yac = yac_all(im)

      return
      end

c**********************************************************************
      subroutine find_cellcenter(x,y,xv,yv,jd,kd,im)
c.. find cell center
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,im
      real x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax)

! local variables
      integer j,k,jv,kv,nh 

      integer idir,ib,js,je,ks,ke
      integer jc,jj,jj1,kc,kk,kk1,k1
      integer iadir,iadd

      x = 0.; y = 0.

      do j = jbeg,jend
        do k = kbeg,kend
          jv = j - nhalo + 1
          kv = k - nhalo + 1
          x(j,k) = 0.25*(xv(jv,kv)+xv(jv-1,kv)+xv(jv,kv-1)+xv(jv-1,kv-1))
          y(j,k) = 0.25*(yv(jv,kv)+yv(jv-1,kv)+yv(jv,kv-1)+yv(jv-1,kv-1))
        enddo
      enddo

!...boundary values (need not be accurate - should not be used anywhere) 

      j = jbeg - 1
      do k = kbeg,kend
        jv = j - nhalo + 1
        kv = k - nhalo + 1
        x(j,k) = x(j+1,k) - 
     &      0.5*( xv(jv+1,kv) - xv(jv,kv) + xv(jv+1,kv-1) - xv(jv,kv-1))
        y(j,k) = y(j+1,k) -
     &      0.5*( yv(jv+1,kv) - yv(jv,kv) + yv(jv+1,kv-1) - yv(jv,kv-1))
      enddo

      j = jend + 1
      do k = kbeg,kend
        jv = j - nhalo + 1
        kv = k - nhalo + 1
        x(j,k) = x(j-1,k) - 
     &      0.5*( xv(jv-2,kv) - xv(jv-1,kv) + xv(jv-2,kv-1) - xv(jv-1,kv-1))
        y(j,k) = y(j-1,k) -
     &      0.5*( yv(jv-2,kv) - yv(jv-1,kv) + yv(jv-2,kv-1) - yv(jv-1,kv-1))
      enddo
 
      k = kbeg - 1
      do j = jbeg,jend
        jv = j - nhalo + 1
        kv = k - nhalo + 1
        x(j,k) = x(j,k+1) -
     &      0.5*( xv(jv,kv+1) - xv(jv,kv) + xv(jv-1,kv+1) - xv(jv-1,kv))
        y(j,k) = y(j,k+1) -
     &      0.5*( yv(jv,kv+1) - yv(jv,kv) + yv(jv-1,kv+1) - yv(jv-1,kv))
      enddo

      k = kend + 1
      do j = jbeg,jend
        jv = j - nhalo + 1
        kv = k - nhalo + 1
        x(j,k) = x(j,k-1) -
     &      0.5*( xv(jv,kv-2) - xv(jv,kv-1) + xv(jv-1,kv-2) - xv(jv-1,kv-1))
        y(j,k) = y(j,k-1) -
     &      0.5*( yv(jv,kv-2) - yv(jv,kv-1) + yv(jv-1,kv-2) - yv(jv-1,kv-1))
      enddo

      j = jbeg - 1
      k = kbeg - 1
      x(j,k) = 2*x(j+1,k) - x(j+2,k)
      y(j,k) = 2*y(j+1,k) - y(j+2,k)

      k = kend + 1
      x(j,k) = 2*x(j+1,k) - x(j+2,k)
      y(j,k) = 2*y(j+1,k) - y(j+2,k)

      j = jend + 1
      k = kbeg - 1
      x(j,k) = 2*x(j,k+1) - x(j,k+2)
      y(j,k) = 2*y(j,k+1) - y(j,k+2)

      k = kend + 1
      x(j,k) = 2*x(j,k-1) - x(j,k-2)
      y(j,k) = 2*y(j,k-1) - y(j,k-2)

      do nh = 2,nhalo
        j = jbeg - nh
        do k = kbeg - nh + 1,kend + nh - 1
          x(j,k) = 2*x(j+1,k) - x(j+2,k)
          y(j,k) = 2*y(j+1,k) - y(j+2,k)
        enddo

        j = jend + nh
        do k = kbeg - nh + 1,kend + nh - 1
          x(j,k) = 2*x(j-1,k) - x(j-2,k)
          y(j,k) = 2*y(j-1,k) - y(j-2,k)
        enddo
 
        k = kbeg - nh
        do j = jbeg - nh + 1,jend + nh - 1
          x(j,k) = 2*x(j,k+1) - x(j,k+2)
          y(j,k) = 2*y(j,k+1) - y(j,k+2)
        enddo

        k = kend + nh
        do j = jbeg - nh + 1,jend + nh - 1
          x(j,k) = 2*x(j,k-1) - x(j,k-2)
          y(j,k) = 2*y(j,k-1) - y(j,k-2)
        enddo

        j = jbeg - nh
        k = kbeg - nh
        x(j,k) = 2*x(j+1,k) - x(j+2,k)
        y(j,k) = 2*y(j+1,k) - y(j+2,k)

        k = kend + nh
        x(j,k) = 2*x(j+1,k) - x(j+2,k)
        y(j,k) = 2*y(j+1,k) - y(j+2,k)

        j = jend + nh
        k = kbeg - nh
        x(j,k) = 2*x(j,k+1) - x(j,k+2)
        y(j,k) = 2*y(j,k+1) - y(j,k+2)

        k = kend + nh
        x(j,k) = 2*x(j,k-1) - x(j,k-2)
        y(j,k) = 2*y(j,k-1) - y(j,k-2)

      enddo

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
        iadir = abs(idir)
        iadd = sign(1,idir)

        if (ibtyp_all(ib,im).eq.22) then

          if(idir.eq.1) then
            jj  = je - js + 1
            do j = js,je
               jj1 = j - js + 1
               jc = jd - 2*jj + jj1
 
               do k = ks,ke
                 x(j,k) = x(jc,k)
                 y(j,k) = y(jc,k)
               enddo
            enddo

          elseif(idir.eq.-1) then
            jj  = je - js + 1
            do j = js,je
               jj1 = je - j + 1
               jc = 1 + 2*jj - jj1
 
               do k = ks,ke
                 x(j,k) = x(jc,k)
                 y(j,k) = y(jc,k)
               enddo
            enddo

          elseif(idir.eq.2) then
            kk  = ke - ks + 1
            do k = ks,ke
               kk1 = k - ks + 1
               kc = kd - 2*kk + kk1
 
               do j = js,je
                 x(j,k) = x(j,kc)
                 y(j,k) = y(j,kc)
               enddo
            enddo

          elseif(idir.eq.-2) then
            kk  = ke - ks + 1
            do k = ks,ke
               kk1 = ke - k + 1
               kc = 1 + 2*kk - kk1 
 
               do j = js,je
                 x(j,k) = x(j,kc)
                 y(j,k) = y(j,kc)
               enddo
            enddo

          endif

        elseif (ibtyp_all(ib,im).eq.51) then

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

              do j = js, je
                jj = jd - j + 1
                x(j,k)  = x(jj,k1)
                y(j,k)  = y(jj,k1)
                x(jj,k)  = x(j,k1)
                y(jj,k)  = y(j,k1)
              enddo
           enddo
          
          endif
        endif
      enddo

      end subroutine find_cellcenter

c***********************************************************************
      subroutine find_coord_TS(xv,yv,xv1,yv1,tspec,jd,kd)
c
c***********************************************************************
      use params_global
      implicit none

      integer :: jd,kd
      real :: xv(jmax,kmax),yv(jmax,kmax),xv1(jmax,kmax),yv1(jmax,kmax)

c...local variables

      integer :: j,k
      real :: theta,cs,ss
      real :: tspec
 
      pi = 4.0*atan(1.0)

      if (motiontype.eq.1) then
!am?        theta=angmax*pi/180.*sin(rotf*tspec)
        theta=angmax*pi/180.*(1.-cos(rotf*tspec))
        theta_col = alfa+theta*180./pi !asitav
        print*,'tspec,theta_col: ',tspec,theta_col
        ss = sin(theta)
        cs = cos(theta)
c
        do k = 1,kmax
          do j = 1,jmax
            xv(j,k) = (xv1(j,k)-xac)*cs + (yv1(j,k)-yac)*ss + xac
            yv(j,k) = (yv1(j,k)-yac)*cs - (xv1(j,k)-xac)*ss + yac
          enddo
        enddo
      elseif (motiontype.eq.2) then
        cs = cos(rotf*tspec)
        ss = sin(rotf*tspec)

        do k = 1,kmax
          do j = 1,jmax
            xv(j,k) = (xv1(j,k)-xac)*cs - (yv1(j,k)-yac)*ss + xac
            yv(j,k) = (yv1(j,k)-yac)*cs + (xv1(j,k)-xac)*ss + yac
          enddo
        enddo
      elseif (motiontype.eq.3) then
        call move_tef_ts(xv,yv,xv1,yv1,tspec)
      endif

      end subroutine find_coord_TS

c***********************************************************************
      subroutine find_coord_TS_orig(xv,yv,xv1,yv1,tspec,jd,kd)
c
c***********************************************************************
      use params_global
      implicit none

      integer :: jd,kd
      real :: xv(jmax,kmax),yv(jmax,kmax),xv1(jmax,kmax),yv1(jmax,kmax)

c...local variables

      integer :: j,k
      real :: theta,cs,ss
      real :: tspec
 
      pi = 4.0*atan(1.0)

      if (motiontype.eq.1) then
        theta=angmax*pi/180.*sin(rotf*tspec)
        ss = sin(theta)
        cs = cos(theta)
c
        do k = 1,kmax
          do j = 1,jmax
            xv(j,k) = (xv1(j,k)-xac)*cs + (yv1(j,k)-yac)*ss + xac
            yv(j,k) = (yv1(j,k)-yac)*cs - (xv1(j,k)-xac)*ss + yac
          enddo
        enddo
      elseif (motiontype.eq.2) then
        cs = cos(rotf*tspec)
        ss = sin(rotf*tspec)

        do k = 1,kmax
          do j = 1,jmax
            xv(j,k) = (xv1(j,k)-xac)*cs - (yv1(j,k)-yac)*ss + xac
            yv(j,k) = (yv1(j,k)-yac)*cs + (xv1(j,k)-xac)*ss + yac
          enddo
        enddo
      elseif (motiontype.eq.3) then
        call move_tef_ts(xv,yv,xv1,yv1,tspec)
      endif

      end subroutine find_coord_TS_orig
c***********************************************************************
      subroutine metfv( q,x,y,xv,yv,xx,xy,yx,yy,jd,kd,im)
c  finite volume formulation         11/4/87  s.o.
c  compute the metrics and the jacobian for computational space
c  uniform computational space, deltas = 1 < averaged metrics >
c  modified to incorporate refined mesh 1996
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,im
 
      real q(jd,kd,nq), x(jd,kd),y(jd,kd)
      real xv(jmax,kmax), yv(jmax,kmax)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables

      integer,allocatable :: jjp(:),jjr(:),kkp(:),kkr(:)
      real fac1,fac2,dx2,dy2,qj,a3
      integer j,k,jc,kc,j1,k1,nh,nneg

      integer idir,ib,js,je,ks,ke
      integer jj,jj1,kk,kk1
      integer iadir,iadd

      allocate(jjp(jd),jjr(jd),kkp(kd),kkr(kd))

c***  first executable statement

      fac1    = 0.5
      fac2    = 0.5

      do 1  j = 1,jd
       jjp(j) = j + 1
       jjr(j) = j - 1
    1 continue
      jjp(jd)=jd
      jjr(1   )=1
      do 2  k = 1,kd
        kkp(k) = k + 1
        kkr(k) = k - 1
    2 continue
      kkp(kd)=kd
      kkr(1   )=1

c..find the volume near j,k

      do k=2,kmax
      do j=2,jmax
          jc = j + nhalo - 1 
          kc = k + nhalo - 1 
          dx1 = xv(j,k)-xv(j-1,k-1)
          dy1 = yv(j,k)-yv(j-1,k-1)
          dx2 = xv(j-1,k)-xv(j,k-1)
          dy2 = yv(j-1,k)-yv(j,k-1)
          a3 = 0.5*( dx1*dy2 -dx2*dy1 )
          q(jc,kc,nq) = 1./a3
      enddo
      enddo

c..xi derivatives
c
      do j=2,jmax
        do k=2,kmax
          jc = j + nhalo - 1 
          kc = k + nhalo - 1 
          xx(jc,kc) = 0.5*(yv(j,k) - yv(j,k-1) + yv(j-1,k) - yv(j-1,k-1))
          xy(jc,kc) =-0.5*(xv(j,k) - xv(j,k-1) + xv(j-1,k) - xv(j-1,k-1))
        enddo
      enddo
c
c..eta derivatives
c
      do k=2,kmax
        do j=2,jmax
          jc = j + nhalo - 1 
          kc = k + nhalo - 1 
          yx(jc,kc) =-0.5*(yv(j,k) - yv(j-1,k) + yv(j,k-1) - yv(j-1,k-1))
          yy(jc,kc) = 0.5*(xv(j,k) - xv(j-1,k) + xv(j,k-1) - xv(j-1,k-1))
        enddo
      enddo
      do k = kbeg,kend
        do j = jbeg,jend
          qj     = q(j,k,nq)
          xx(j,k)= qj*(xx(j,k))
          xy(j,k)= qj*(xy(j,k))
          yx(j,k)= qj*(yx(j,k))
          yy(j,k)= qj*(yy(j,k))
        enddo
      enddo
c
c..boundary values
c
      do nh = 1,nhalo
        j = jbeg - nh
        j1 = jbeg + nh - 1
        do k = kbeg - nh + 1,kend + nh - 1
          q(j,k,nq) = q(j1,k,nq)
          xx(j,k) = xx(j1,k)
          xy(j,k) = xy(j1,k)
          yx(j,k) = yx(j1,k)
          yy(j,k) = yy(j1,k)
        enddo

        j = jend + nh
        j1 = jend - nh + 1
        do k = kbeg - nh + 1,kend + nh - 1
          q(j,k,nq) = q(j1,k,nq)
          xx(j,k) = xx(j1,k)
          xy(j,k) = xy(j1,k)
          yx(j,k) = yx(j1,k)
          yy(j,k) = yy(j1,k)
        enddo
 
        k = kbeg - nh
        k1 = kbeg + nh - 1
        do j = jbeg - nh,jend + nh
          q(j,k,nq) = q(j,k1,nq)
          xx(j,k) = xx(j,k1)
          xy(j,k) = xy(j,k1)
          yx(j,k) = yx(j,k1)
          yy(j,k) = yy(j,k1)
        enddo

        k = kend + nh
        k1 = kend - nh + 1
        do j = jbeg - nh,jend + nh
          q(j,k,nq) = q(j,k1,nq)
          xx(j,k) = xx(j,k1)
          xy(j,k) = xy(j,k1)
          yx(j,k) = yx(j,k1)
          yy(j,k) = yy(j,k1)
        enddo
      enddo

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
        iadir = abs(idir)
        iadd = sign(1,idir)

        if (ibtyp_all(ib,im).eq.22) then

          if(idir.eq.1) then
            jj  = je - js + 1
            do j = js,je
               jj1 = j - js + 1
               jc = jd - 2*jj + jj1
 
               do k = ks,ke
                 q(j,k,nq) = q(jc,k,nq)
                 xx(j,k) = xx(jc,k)
                 xy(j,k) = xy(jc,k)
                 yx(j,k) = yx(jc,k)
                 yy(j,k) = yy(jc,k)
               enddo
            enddo

          elseif(idir.eq.-1) then
            jj  = je - js + 1
            do j = js,je
               jj1 = je - j + 1
               jc = 1 + 2*jj - jj1
 
               do k = ks,ke
                 q(j,k,nq) = q(jc,k,nq)
                 xx(j,k) = xx(jc,k)
                 xy(j,k) = xy(jc,k)
                 yx(j,k) = yx(jc,k)
                 yy(j,k) = yy(jc,k)
               enddo
            enddo

          elseif(idir.eq.2) then
            kk  = ke - ks + 1
            do k = ks,ke
               kk1 = k - ks + 1
               kc = kd - 2*kk + kk1
 
               do j = js,je
                 q(j,k,nq) = q(j,kc,nq)
                 xx(j,k) = xx(j,kc)
                 xy(j,k) = xy(j,kc)
                 yx(j,k) = yx(j,kc)
                 yy(j,k) = yy(j,kc)
               enddo
            enddo

          elseif(idir.eq.-2) then
            kk  = ke - ks + 1
            do k = ks,ke
               kk1 = ke - k + 1
               kc = 1 + 2*kk - kk1 
 
               do j = js,je
                 q(j,k,nq) = q(j,kc,nq)
                 xx(j,k) = xx(j,kc)
                 xy(j,k) = xy(j,kc)
                 yx(j,k) = yx(j,kc)
                 yy(j,k) = yy(j,kc)
               enddo
            enddo

          endif

        elseif (ibtyp_all(ib,im).eq.51) then

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

              do j = js, je
                jj = jd - j + 1
                q(j,k,nq) = q(jj,k1,nq)
                xx(j,k)  = -xx(jj,k1)
                xy(j,k)  = -xy(jj,k1)
                yx(j,k)  = -yx(jj,k1)
                yy(j,k)  = -yy(jj,k1)

                q(jj,k,nq) = q(j,k1,nq)
                xx(jj,k)  = -xx(j,k1)
                xy(jj,k)  = -xy(j,k1)
                yx(jj,k)  = -yx(j,k1)
                yy(jj,k)  = -yy(j,k1)
              enddo
            enddo

          endif

        endif
      enddo
c
c..check for negative jacobians
c
      nneg = 0
        do 910 k = 1, kd
          do 920 j = 1, jd
ccray            nneg = cvmgt(nneg+1,nneg,q(j,k,nq).le.0.)
            if( q(j,k,nq).le.0.0 ) then
              nneg = nneg+1
            end if
 920      continue
 910    continue
c
      if(nneg .ne. 0) then
        write(6,*) nneg, ' negative jacobians in block'
        do 74 k = 1,kd
        do 74 j = 1,jd
          if( q(j,k,nq).le.0.0 ) then
            write(6,603) q(j,k,nq), j, k
          end if
   74   continue
      endif
c
  603 format( ' ',10x,'negative jacobian = ',1p,e10.3,1x,'at j,k =',
     $                 2i5,5x)
c
      return
      end

c***********************************************************************
      subroutine find_timemetrics_TS(x,y,xv,yv,tspec,ug,vg,ugv,vgv,jd,kd)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax),tspec
      real ug(jd,kd), vg(jd,kd),ugv(jmax,kmax),vgv(jmax,kmax)

c...local variables

      integer :: j,k

      pi = 4.0*atan(1.0)

      if (motiontype.eq.1) then
!am?        do k = 1,kd
!am?        do j = 1,jd
!am?          ug(j,k) =  angmax*pi/180.*rotf*(y(j,k)-yac)*cos(rotf*tspec)
!am?          vg(j,k) = -angmax*pi/180.*rotf*(x(j,k)-xac)*cos(rotf*tspec)
!am?        enddo
!am?        enddo
!am?        do k = 1,kmax
!am?        do j = 1,jmax
!am?          ugv(j,k) =  angmax*pi/180.*rotf*(yv(j,k)-yac)*cos(rotf*tspec)
!am?          vgv(j,k) = -angmax*pi/180.*rotf*(xv(j,k)-xac)*cos(rotf*tspec)
!am?        enddo
!am?        enddo

        do k = 1,kd
        do j = 1,jd
          ug(j,k) =  angmax*pi/180.*rotf*(y(j,k)-yac)*sin(rotf*tspec)
          vg(j,k) = -angmax*pi/180.*rotf*(x(j,k)-xac)*sin(rotf*tspec)
        enddo
        enddo
        do k = 1,kmax
        do j = 1,jmax
          ugv(j,k) =  angmax*pi/180.*rotf*(yv(j,k)-yac)*sin(rotf*tspec)
          vgv(j,k) = -angmax*pi/180.*rotf*(xv(j,k)-xac)*sin(rotf*tspec)
        enddo
        enddo
      elseif (motiontype.eq.2) then
        do k = 1,kd
        do j = 1,jd
          ug(j,k) = - rotf*(y(j,k)-yac)
          vg(j,k) =   rotf*(x(j,k)-xac)
        enddo
        enddo
        do k = 1,kmax
        do j = 1,jmax
          ugv(j,k) = - rotf*(yv(j,k)-yac)
          vgv(j,k) =   rotf*(xv(j,k)-xac)
        enddo
        enddo
      elseif(motiontype.eq.3) then
        call gettimemetrics_tef_TS(x,y,xv,yv,ug,vg,ugv,vgv,jd,kd,tspec)
      endif

      end subroutine find_timemetrics_TS

c***********************************************************************
      subroutine find_timemetrics_TS_orig(x,y,xv,yv,tspec,ug,vg,ugv,vgv,jd,kd)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax),tspec
      real ug(jd,kd), vg(jd,kd),ugv(jmax,kmax),vgv(jmax,kmax)

c...local variables

      integer :: j,k

      pi = 4.0*atan(1.0)

      if (motiontype.eq.1) then
        do k = 1,kd
        do j = 1,jd
          ug(j,k) =  angmax*pi/180.*rotf*(y(j,k)-yac)*cos(rotf*tspec)
          vg(j,k) = -angmax*pi/180.*rotf*(x(j,k)-xac)*cos(rotf*tspec)
        enddo
        enddo
        do k = 1,kmax
        do j = 1,jmax
          ugv(j,k) =  angmax*pi/180.*rotf*(yv(j,k)-yac)*cos(rotf*tspec)
          vgv(j,k) = -angmax*pi/180.*rotf*(xv(j,k)-xac)*cos(rotf*tspec)
        enddo
        enddo
      elseif (motiontype.eq.2) then
        do k = 1,kd
        do j = 1,jd
          ug(j,k) = - rotf*(y(j,k)-yac)
          vg(j,k) =   rotf*(x(j,k)-xac)
        enddo
        enddo
        do k = 1,kmax
        do j = 1,jmax
          ugv(j,k) = - rotf*(yv(j,k)-yac)
          vgv(j,k) =   rotf*(xv(j,k)-xac)
        enddo
        enddo
      elseif(motiontype.eq.3) then
        call gettimemetrics_tef_TS(x,y,xv,yv,ug,vg,ugv,vgv,jd,kd,tspec)
      endif

      end subroutine find_timemetrics_TS_orig

c***********************************************************************
      subroutine computeTScoefs(N,Ds)
c***********************************************************************
        implicit none
c***********************************************************************
        integer :: N
        real :: Ds(N,N)

!local variables
        integer :: i,j
        real :: pi
      
        pi = 4.*atan(1.)
      
        if (mod(N,2).eq.0) then
          do i = 1,N
            do j = 1,N
              if (i.eq.j) then
                Ds(i,j) = 0.
              else
                Ds(i,j) = 0.5*(-1)**(i-j)/tan(pi*(i-j)/N)
              endif
            enddo
          enddo
        else
          do i = 1,N
            do j = 1,N
              if (i.eq.j) then
                Ds(i,j) = 0.
              else
                Ds(i,j) = 0.5*(-1)**(i-j)/sin(pi*(i-j)/N)
              endif
            enddo
          enddo
        endif
      
      end subroutine computeTScoefs

c**********************************************************************
      subroutine read_inputs(jgmx,kgmx,nmesh,bcfile)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer nmesh
      integer jgmx(nmesh),kgmx(nmesh)
      character*40,bcfile

      ! local variables
      real a,alphamean,cs,ss,tu
      integer istop,i,i1,ii,im,ihar
      logical file_exists

C..Parameter for boundary conditions
      integer :: irot_dir
      integer :: nbc,ibtyp(25),ibdir(25)
      integer :: jbcs(25),kbcs(25)
      integer :: jbce(25),kbce(25)
      real    :: xcen,ycen
      
      namelist/bcinp/ nbc,irot_dir,xcen,ycen,ibtyp,ibdir,jbcs,jbce,kbcs,kbce

      namelist/flapinp/ pxc,dela,rf_tef,
     &                 nharmflap,ampFlap,phiFlap,theta_f0,theta_finit

c***  first executable statement

c..read in the input and write it back out

      read (5, inputs)
      write (6, inputs)

      inquire(file=bcfile, exist=file_exists)
      
      if (file_exists) then
        open(unit=10,file=bcfile,status='unknown')
        do im = 1,num_grids
          irot_dir = 1
          xcen = 0.; ycen = 0.
          read(10,bcinp)
          write(6,bcinp)
          nbc_all(im) = nbc
          irot_dir_all(im) = irot_dir
          xac_all(im) = xcen
          yac_all(im) = ycen
          ibtyp_all(:,im) = ibtyp(:)
          ibdir_all(:,im) = ibdir(:)
          jbcs_all(:,im) = jbcs(:)
          jbce_all(:,im) = jbce(:)
          kbcs_all(:,im) = kbcs(:)
          kbce_all(:,im) = kbce(:)
        enddo
        close(10)
      else
!        write(6,*) 'USING DEFAULT BC'
!        write(6,*) '    - FIRST MESH --> C-type Airfoil Mesh'
!        write(6,*) '    - SECOND MESH --> Rectangular Wind-tunnel Mesh '
        write(6,*),"Error: need BC file -> ",bcfile
        stop
      endif 

c..In case of timespectral calculation

      if (timespectral) then
        iunst = 0
      else
        nspec = 1
      endif

      pi=4*atan(1.0)

!asitav
!----------------------------------------------
c..read in the TEF inputs and write it back out
      if (iteflap.eq.1) then
        rewind(5)
        read(5,flapinp)
        write (6, flapinp)

        dela=dela*pi/180
        theta_f0 = -theta_f0*pi/180.0
        theta_finit = -theta_finit*pi/180.0
        do ihar=1,nharmflap
           ampFlap(ihar)=-ampFlap(ihar)*pi/180.
           phiFlap(ihar)= phiFlap(ihar)*pi/180.
        enddo
      endif
!----------------------------------------------

c..calculate a few things based on inputs

      if (invisc) then
        rmue = 0.0
        lamin = .true.
      endif
      if (lamin) iturb = -1

c..write summary of inputs (and check validity!)

      istop = 0
      write(6,*) ' '
      write(6,*) 'here is your input:'

      write(6,*) 'restart'
      if(iread.eq.0) then
        write(6,*) '  initial run w/o restart for this case '
      elseif(iread.eq.1) then
        write(6,*) '  restart run for this case '
      else
        write(6,*) '**invalid iread, must be 0 or 1**'
        istop = istop+1
      endif

      write(6,*) '  this run will last until ',nsteps,' time steps'
      write(6,*) 'flow parameters'
      write(6,*) '  free stream mach number is ',fsmach
      if(alfa.ne.0.) write(6,*) '  flow at ',
     <               alfa,' degrees angle of attack'
      if(invisc) then
        write(6,*) '  inviscid flow'
      else
        write(6,*) '  viscous flow'
        write(6,*) '    reynolds number = ',rey
        if(lamin) then
          write(6,*) '      laminar flow'
        else
          write(6,*) '      turbulent flow'
          if (iturb.eq.0) then
            write(6,*) '     Baldwin-Lomax turbulence model'
          elseif (iturb.eq.1) then
            write(6,*) '     Spalart-Allmaras turbulence model'
            if (itrans.eq.1) then
              write(6,*) '     Modified Langter-Mentry transition model'
            endif
          elseif (iturb.eq.2) then
            write(6,*) '     k-omega SST turbulence model'
            if (itrans.eq.1) then
              write(6,*) '       Langter-Mentry transition model'
            endif
          endif
        endif
      endif

      write(6,*) 'time info'
      if(iunst.eq.0) then
        write(6,*) '  steady flow '
      elseif(iunst.eq.1) then
        write(6,*) '  unsteady flow (pitching oscillation)'
      elseif(iunst.eq.2) then
        write(6,*) '  unsteady flow (pitching ramp)'
      elseif(iunst.eq.3) then
        write(6,*) '  unsteady flow (prescribed pitching)' 
      elseif(iunst.eq.4) then
        write(6,*) '  unsteady flow (pitching oscillation about 1/4c)'
      elseif(iunst.eq.5) then
        write(6,*) '  unsteady flow (pitching oscillation about 1/4c
     c               with deformation)'
      else
        write(6,*) '**invalid iunst, must be between 0 and 5**'
        istop = istop+1
      endif
      if(ntac.eq.1) then
        write(6,*) '  1st order in time '
      elseif(ntac.eq.2) then
        write(6,*) '  2nd order in time '
      elseif(ntac.eq.3) then
        write(6,*) '  3rd order in time '
      else
        write(6,*) '**invalid ntac, must be 1 or 2 or 3**'
        istop = istop+1
      endif
      if(itnmax.eq.1) then
        write(6,*) '  no use of newton iterations '
      elseif(itnmax.gt.1) then
        write(6,*) '  using newton iterations with ',
     <        itnmax,' iterations'
      else
        write(6,*) '**invalid itnmax, must be 0 or greater**'
        istop = istop+1
      endif
      if(iunst.gt.0 .and. dt.gt.0.10) then
        write(6,*) '**invalid dt for unsteady flow**'
        istop=istop+1
      elseif(dt.gt.0.) then
        write(6,*) '  dt is ',dt
      else
        write(6,*) '**invalid dt, must be greater than 0'
        istop=istop+1
      endif
      if(timeac.ne.1.0 .and. iunst.gt.0) then
        write(6,*) '**invalid timeac for unsteady flow**'
        istop=istop+1
      elseif(timeac.le.1.0 .and. timeac.ge.0.0) then
        write(6,*) '  timeac for jacobian scaling is ',timeac
      else
        write(6,*) '**invalid timeac, must be between 0 and 1**'
        istop=istop+1
      endif
c     
      write(6,*) 'algorithm stuff'
      if(epse.gt. 0.0 .and. epse.lt. 1.0) then
        write(6,*) '  dissipation for ilu3d is ',epse
      else
        write(6,*) '**invalid epse, must be between 0 and 1**'
        istop=istop+1
      endif
      if(irhsy.eq.-1) then
        write(6,*) '  1st order in space '
      elseif(irhsy.eq.-2) then
        write(6,*) '  2nd order in space '
      elseif(irhsy.eq.-3) then
        write(6,*) '  3rd order in space '
      else
        write(6,*) '**invalid irhsy, must be between -1 and -3**'
        istop=istop+1
      endif
      if(ilim.ge.1) then
        write(6,*) '  limiting in both directions'
      elseif(ilim.le.-1) then
        write(6,*) '  limiting in only the j-direction'
      elseif(ilim.eq.0) then
        write(6,*) '  no limiting'
      endif
      if(abs(ilim).eq.0) then
        write(6,*) '     using muscl scheme'
      elseif(abs(ilim).eq.1) then
        write(6,*) '     using muscl with korens differentiable limiter'
      elseif(abs(ilim).eq.2) then
        write(6,*) '     using muscl with van a. differentiable limiter'
      elseif(abs(ilim).eq.3) then
        write(6,*) '     using muscl with c-o minmod scheme'            
      elseif(abs(ilim).eq.4) then
        write(6,*) '     using sonic-a scheme of hyunh et al.'          
      elseif(abs(ilim).eq.5) then
        write(6,*) '     using sonic extension of c-o minmod scheme'    
      elseif(abs(ilim).eq.7) then
        write(6,*) '     using cubic interpolation with no limiting'    
      elseif(abs(ilim).eq.8) then
        write(6,*) '     using quartic interpolation with no limiting'  
      elseif(abs(ilim).eq.9) then
        write(6,*) '     using quadratic reconstruction (pade) '        
      endif
      if(iread.eq.0) then
        write(6,*) '  totime = ',totime
      else
        write(6,*) '  totime = ',totime,' will be reset from q-file'
      endif
c
      if(iunst.eq.1) then
        write(6,*) 'reduced frequency is ',rf
        write(6,*) 'amplitude of pitching oscillation is (deg) ',angmax
      endif
c     
      if(jint+kint.gt.2) then
        write(6,*) 'coarsening of grid for this run, use with caution'
        if(jint.ne.1)
     <     write(6,*) ' use every ',jint,' points in j-direction'
        if(kint.ne.1)
     <     write(6,*) ' use every ',kint,' points in k-direction'
      endif
c
      if(istop.gt.0) write(6,*) 'note: ',istop,' errors in input******'
      write(6,*) ' '
c
c...setting pseudo-time step cycle
c
      ii = itnmax + 1
      do i = 1,itnmax
        if (dtpseudo(i).eq.0) then
          ii = i
          exit 
        endif
      enddo
      if (ii.eq.1) then
        dtpseudo = 1.
      else
        do i = ii,itnmax
          i1 = mod(i,ii-1)
          if (i1.eq.0) i1 = ii-1
          dtpseudo(i) = dtpseudo(i1)
        enddo
      endif

c...setting freestream conditions

      gm1   = gamma -1.0
      ggm1  = gamma*gm1
      cs    = cos( pi*alfa/180.0 )
      ss    = sin( pi*alfa/180.0 )
      theta_col=alfa
c
      einf  = 1.0/ggm1 +0.5*fsmach**2
      pinf  = 1.0/gamma
      rinf  = 1.0
      ainf  = gamma*pinf/rinf
      htinf = gamma*einf/rinf - 0.5*gm1*fsmach**2

      nq = 5
      nv = 4
      nmv = 4

      if (iturb.eq.1) then
        nturb = 1
        nq = nq + nturb
        nv = nv + nturb
      elseif (iturb.eq.2) then
        nturb = 2
        nq = nq + nturb
        nv = nv + nturb
      endif

      if (iturb.eq.1.or.iturb.eq.2) then
        if (itrans.eq.1) then
          nq = nq + 2
          nv = nv + 2
        endif
      endif
c
      uinf  = fsmach*cs
      vinf  = fsmach*ss
c
      if (fsmach.ne.0) then
        rey = rey/fsmach
      else
        rey = rey/fmtip
      endif

      if (motiontype.eq.1.and.iunst.eq.2) then
        iunst = 1
        write(6,*) "Warning: iunst = 2 not supported for pitching motion, 
     &reset to 1"
      endif

      if (motiontype.eq.2) rf = fmtip/rartio

      if(dt.lt.0.0 .and. iunst.ge.1) then
         !amm if (motiontype.eq.1.or.motiontype.eq.2) then
         if (motiontype.eq.1.or.motiontype.eq.2.or.motiontype.eq.3) then
           if (rf.ne.0) then
             dt = abs(dt)*pi/abs(rf)/180.
           !asitav tef
           !----------------
           elseif (rf_tef.ne.0..and.motiontype.eq.3) then     
             dt = abs(dt)*pi/abs(rf_tef)/180.
           !----------------
           else
             stop "Error: negative time inputted with zero frequency"
           endif
         endif
      endif

      h  = dt
      hd = .5*dt
c
      if (itrans.eq.0) then
        tkeinf = 1e-6*fsmach**2 
        tomegainf = 5*fsmach/rey 
      else
        !tu = 3
        !tkeinf = (tu/100*fsmach)**2*1.5
        !tomegainf = tkeinf/12
     
        !tu = 0.1
        !tkeinf = (tu/100*fsmach)**2*1.5
        !tomegainf = tkeinf/10

        tu = tuinf
        tkeinf = (tu/100*fsmach)**2*1.5
        tomegainf = tkeinf/vnuinf
      endif

	 nmesh=num_grids

!      if (iread.eq.0) then
         if(num_grids.gt.1) read(1) nmesh
         read(1) (jgmx(i),kgmx(i),i=1,nmesh)
         write(6,*) (jgmx(i),kgmx(i),i=1,nmesh)
!      else
!
!         if(num_grids.gt.1) read(4) nmesh
!         read(4) (jgmx(i),kgmx(i),i=1,nmesh)
!
!      endif

      nadd = 2*nhalo - 1
      jgmx = jgmx + nadd
      kgmx = kgmx + nadd

      return
      end

      subroutine obj_setup(jd,kd)
      use params_global
      implicit none
      integer :: jd, kd
      
      integer j, k, ib, np,unitcp,nsp,nobjf
      integer maxj, maxk
      logical file_exists
      real tmp
      character(60) :: int_str,cpfile

      namelist/objinp/obj_ftype,obj_pointwise,obj_npoints,obj_points,obj_normtype,obj_samod,obj_if_reg,obj_reg_fac

      inquire(file='obj.inp', exist=file_exists)
      
      if (file_exists) then
        open(unit=10,file='obj.inp',status='unknown')
        read(10,objinp)
        write(6,objinp)
        close(10)
      else
         print *, "obj.inp not found"
         stop
      endif 

      
      if(obj_npoints .gt. 0) then
         print *, "Using specified points in obj_setup"
      else
         print *, "Setting points from jtail1 to jtail2 in 
     &obj_setup and reseting obj_npoints"
         obj_npoints = jtail2 - jtail1 +  1
         do np = 1, obj_npoints
            obj_points(np) = jtail1 + np - 1
         end do
      end if
      
      
      
      if(obj_ftype .eq. obj_ftype_cp) then
         print *, "Benchmark type is coefficient of pressure"

         allocate(obj_f(jd*nspec))
         obj_f = 0.0
         !for plotting purpose
         allocate(obj_cpx(jd),obj_cpy(jd))

!$OMP PARALLEL REDUCTION(+:obj_f) 
!$OMP& IF(NSPEC > 1)
!$OMP& PRIVATE(tmp,unitcp,cpfile,int_str,j,nobjf)
!$OMP DO ORDERED
         spectralloop: do nsp=1,nspec
           unitcp = 221100+nsp

           write(int_str,'(I7)')nsp        
           cpfile = 'cp_benchmark_TS'//trim(adjustl(int_str))//'.dat'
           print*,'nsp,cpfile: ',nsp,trim(adjustl(cpfile))
           !am inquire(file="cp_benchmark.dat", exist=file_exists)
           inquire(file=cpfile, exist=file_exists)
           if(file_exists) then
              open(unitcp, file=cpfile, form="formatted")
              do j=jtail1,jtail2
                nobjf = jd*(nsp-1)+j
                !am read(unitcp,*) tmp, tmp, obj_f(nobjf)
                read(unitcp,*) obj_cpx(j), obj_cpy(j), obj_f(nobjf)
              end do
              close(unitcp)
           else
              print *, 'The cp benchmark does not exist, 
     >                  setting it to zero'
              obj_f(:) = 0.0
           end if
        enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL
!
      elseif(obj_ftype .eq. obj_ftype_cltot) then
         if_objtot = .TRUE.
         print *, "Benchmark type is total coefficient of lift"
         allocate(obj_f(nspec))
         allocate(clobj(nspec),cdobj(nspec),cmobj(nspec)) !used in forces.f to save airloads
         inquire(file="cl_benchmark.dat", exist=file_exists)
         if(file_exists) then
            open(2211, file="cl_benchmark.dat", form="formatted")
            !am read(2211,*) ((obj_f(j),tmp), j=1,nspec) !format: cl(nsp),nsp
            do j=1,nspec
              read(2211,*) obj_f(j),tmp !format: cl(nsp),nsp
            end do
            close(2211)
         else
            print *, 'The cl benchmark does not exist, 
     &setting it to zero'
            obj_f(:) = 0.0
         end if
!
      else
         print *, "Benchmark type is not defined"
         stop
      end if
      
      end subroutine obj_setup


      subroutine reset_qb(jd,kd,qb)
      use params_global
      implicit none
      
      integer :: jd,kd
      real :: qb(jd,kd,nq)
      
      qb(:,:,:) = 0.0
      end subroutine reset_qb

c**********************************************************************
      subroutine read_coeff_prod(coeff_prod, jd, kd)
      implicit none
      integer jd, kd
      real :: coeff_prod(jd, kd)
      integer j, k
      open(2211,file='beta.opt',form='formatted')
      read(2211, *) ((coeff_prod(j,k),j=1,jd),k=1,kd)
      close(2211)
      return
      end
