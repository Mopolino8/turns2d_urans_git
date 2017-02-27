c***********************************************************************
      subroutine preilu2d(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     >            tscale,bt)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd),bt(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables
      real,allocatable :: d(:,:),st(:,:,:),tmp(:)
      real,allocatable :: a(:,:),b(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)
      integer,allocatable :: ms(:),me(:)
      integer j,k,m,i,k1,jj
      integer itgs,j1,n
      real eps2,epsv,dj,uu,vv,uv2,cjkl,rj1,rj2
      real qq1,qqx,rr1,rk1,rk2,qq2,qqy,rr2,vnu,svt
      real ri1,ri2,qq,cc,sp1,sm2,chkx,spec,s1,s2,bSq,X,Y
      real a2,a4,chky,sv,di,s3,s4
      real scale1,scale2,sav1,sav2,sav3,sav4,scale,term

      allocate(d(jd,kd),st(jd,kd,4),tmp(4))
      allocate(a(mdim,4),b(mdim,4))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))
      allocate(ms(mdim*2),me(mdim*2))
c***  first executable statement

      eps2 = epse*2.
      epsv = 1. + eps2

c..set-up for hyper-plane loop

      do 1 m=js+ks,je+ke
        i     = max(js,m-ke)
        ms(m) = i
        me(m) = min(je,m-js)
    1 continue
c
c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..setup d, the diagonal term and and store some arrays
c
      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
          uu  = q(j,k,2)
          vv  = q(j,k,3)
          uv2 = 0.5*( uu*uu + vv*vv )
          cjkl = sqrt( ggm1*( q(j,k,4) - uv2 ) )
c
          rj1 = xx(j,k)
          rj2 = xy(j,k)
          qq1 = rj1*uu + rj2*vv
          qqx = abs( -rj1*ug(j,k)-rj2*vg(j,k) + qq1 )
          rr1 = sqrt( rj1**2 + rj2**2 )
          rk1 = yx(j,k)
          rk2 = yy(j,k)
          qq2 = rk1*uu + rk2*vv
          qqy = abs( -rk1*ug(j,k)-rk2*vg(j,k) + qq2 )
          rr2 = sqrt( rk1**2 + rk2**2 )
          vnu = 2.0*rr2*rr2*(rmue+turmu(j,k))/(rey*q(j,k,1))
          svt = tscale(j,k)
          ge(j,k) = gamma*q(j,k,4) - gm1*uv2

          uv(j,k) = uv2
          qx(j,k) = qq1
          cx(j,k) = cjkl*rr1
          qy(j,k) = qq2
          cy(j,k) = cjkl*rr2
          vn(j,k) = vnu

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          X = sqrt(((1.-bSq)*qqx)*((1.-bSq)*qqx)
     c       +4.*bSq*cx(j,k)*cx(j,k))
          term = 0.5*((bSq+1.)*qqx+X)
          Y = sqrt(((1.-bSq)*qqy)*((1.-bSq)*qqy)
     c       +4.*bSq*cy(j,k)*cy(j,k))
          term = term + 0.5*((bSq+1.)*qqy+Y)
          d(j,k) = 1.0/(1.0+svt*(term*epsv+vnu))

  111 continue

c..apply point Jacobi update for the boundary residuals
c
      j = js-1
      do k = ks,ke
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

      j = je+1
      do k = ks,ke
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

      k = ks-1
      do j = js,je
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

      k = ke+1
      do j = js,je
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

c..store RHS in separate variable
c
      do k = ks-1, ke+1
      do j = js-1, je+1
        st(j,k,1) = s(j,k,1)
        st(j,k,2) = s(j,k,2)
        st(j,k,3) = s(j,k,3)
        st(j,k,4) = s(j,k,4)
      enddo
      enddo

      do 1000 itgs=1,1 
c
c..forward sweep
c..loop on hyper-plane
c
      do 120 m = js+ks,je+ke
c
c..setup a contribution in j-direction
c
      do 121 j = ms(m),me(m)
        k  = m-j
        j1 = j-1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        rr1 = sqrt( ri1**2 + ri2**2 )
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq

        bSq = Mp**2/(bt(j1,k)-Mp**2*(bt(j1,k)-1))

        X = sqrt(((1.-bSq)*qqx)*((1.-bSq)*qqx)
     c     +4.*bSq*cc*cc)

        chkx = 0.5 + sign(0.5,(bSq+1.)*qqx+X)
        sp1 = 0.5*((bSq+1.)*abs(qqx)+X)
        sm2 = eps2*sp1
        spec= chkx*sp1 + sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        a(j,1) = a4                  + chkx*qqx*s1
        a(j,2) = ri1*a2 + uu*a4      + chkx*qqx*s2
        a(j,3) = ri2*a2 + vv*a4      + chkx*qqx*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 + chkx*qqx*s4

        do n=1,4
           tmp(n)=a(j,n)
        enddo
        
        a(j,1) =uv(j1,k)*tmp(1)-uu*tmp(2)-vv*tmp(3)+tmp(4)
        a(j,2) =uu*a(j,1)
        a(j,3) =vv*a(j,1)
        a(j,4) =ge(j1,k)*a(j,1)

        cjkl= sqrt(ggm1*(q(j1,k,4)-uv(j1,k)))
        term=gm1*(bSq-1.)/(cjkl*cjkl)
        do n=1,4
           a(j,n)=a(j,n)*term
           a(j,n)=a(j,n)+tmp(n)
        enddo

        a(j,1) = a(j,1) + spec*s1
        a(j,2) = a(j,2) + spec*s2
        a(j,3) = a(j,3) + spec*s3
        a(j,4) = a(j,4) + spec*s4
         
  121 continue
c
c..setup b contribution in k-direction
c
      do 122 j = ms(m),me(m)
        k  = m-j
        k1 = k-1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        rr2 = sqrt( ri1**2 + ri2**2 )
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq

        bSq = Mp**2/(bt(j,k1)-Mp**2*(bt(j,k1)-1))

        Y = sqrt(((1.-bSq)*qqy)*((1.-bSq)*qqy)
     c     +4.*bSq*cc*cc)

        chky = 0.5 + sign(0.5,(bSq+1.)*qqy+Y)
        sp1 = 0.5*((bSq+1.)*abs(qqy)+Y)
        sm2 = eps2*sp1
        spec= chky*sp1 + vnu + sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        b(j,1) = a4                  + chky*qqy*s1
        b(j,2) = ri1*a2 + uu*a4      + chky*qqy*s2
        b(j,3) = ri2*a2 + vv*a4      + chky*qqy*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 + chky*qqy*s4

        do n=1,4
           tmp(n)=b(j,n)
        enddo
        
        b(j,1) =uv(j,k1)*tmp(1)-uu*tmp(2)-vv*tmp(3)+tmp(4)
        b(j,2) =uu*b(j,1)
        b(j,3) =vv*b(j,1)
        b(j,4) =ge(j,k1)*b(j,1)

        cjkl= sqrt(ggm1*(q(j,k1,4)-uv(j,k1)))
        term=gm1*(bSq-1.)/(cjkl*cjkl)
        do n=1,4
           b(j,n)=b(j,n)*term
           b(j,n)=b(j,n)+tmp(n)
        enddo

        b(j,1) = b(j,1) + spec*s1
        b(j,2) = b(j,2) + spec*s2
        b(j,3) = b(j,3) + spec*s3
        b(j,4) = b(j,4) + spec*s4
         
  122 continue
c
c..bi-diagonal inversion
c     
      if(iunst.eq.2) then
        do 123 j = ms(m),me(m)
          k  = m-j
          svt = tscale(j,k)
          sv = svt*0.5
          di = d(j,k)
          s(j,k,1) = ( st(j,k,1) + sv*( a(j,1)+b(j,1) ) )*di
          s2 = ( st(j,k,2) + sv*( a(j,2)+b(j,2) ) )*di
          s3 = ( st(j,k,3) + sv*( a(j,3)+b(j,3) ) )*di
          s4 = 1.0 + (2.*di*sv*rotf)**2
          s(j,k,2) = ( s2 + di*2.*sv*rotf*s3)/s4
          s(j,k,3) = ( s3 - di*2.*sv*rotf*s2)/s4
          s(j,k,4) = ( st(j,k,4) + sv*( a(j,4)+b(j,4) ) )*di
  123   continue
      else

cdir$ ivdep
        do 124 j = ms(m),me(m)
          k  = m-j
          svt = tscale(j,k)
          sv = svt*0.5
          di = d(j,k)
          s(j,k,1) = ( st(j,k,1) + sv*( a(j,1)+b(j,1) ) )*di
          s(j,k,2) = ( st(j,k,2) + sv*( a(j,2)+b(j,2) ) )*di
          s(j,k,3) = ( st(j,k,3) + sv*( a(j,3)+b(j,3) ) )*di
          s(j,k,4) = ( st(j,k,4) + sv*( a(j,4)+b(j,4) ) )*di
 124    continue
      endif
c
  120 continue
c
c..backward sweep
c..loop on hyper-plane
c
      do 220 m = je+ke,js+ks,-1
c
c..setup a contribution in j-direction
c
      do 221 j = ms(m),me(m)
        k  = m-j
        j1 = j+1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        rr1 = sqrt( ri1**2 + ri2**2 )
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq

        bSq = Mp**2/(bt(j1,k)-Mp**2*(bt(j1,k)-1))

        X = sqrt(((1.-bSq)*qqx)*((1.-bSq)*qqx)
     c     +4.*bSq*cc*cc)

        chkx = 0.5 - sign(0.5,(bSq+1.)*qqx-X)
        sp1 = 0.5*((bSq+1.)*abs(qqx)+X)
        sm2 = eps2*sp1
        spec= chkx*sp1 + sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        a(j,1) = a4                  + chkx*qqx*s1
        a(j,2) = ri1*a2 + uu*a4      + chkx*qqx*s2
        a(j,3) = ri2*a2 + vv*a4      + chkx*qqx*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 + chkx*qqx*s4

        do n=1,4
           tmp(n)=a(j,n)
        enddo
        
        a(j,1) =uv(j1,k)*tmp(1)-uu*tmp(2)-vv*tmp(3)+tmp(4)
        a(j,2) =uu*a(j,1)
        a(j,3) =vv*a(j,1)
        a(j,4) =ge(j1,k)*a(j,1)

        cjkl= sqrt(ggm1*(q(j1,k,4)-uv(j1,k)))
        term=gm1*(bSq-1.)/(cjkl*cjkl)
        do n=1,4
           a(j,n)=a(j,n)*term
           a(j,n)=a(j,n)+tmp(n)
        enddo

        a(j,1) = a(j,1) - spec*s1
        a(j,2) = a(j,2) - spec*s2
        a(j,3) = a(j,3) - spec*s3
        a(j,4) = a(j,4) - spec*s4
         
  221 continue
c
c..setup b contribution in k-direction
c
      do 222 j = ms(m),me(m)
        k  = m-j
        k1 = k+1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        rr2 = sqrt( ri1**2 + ri2**2 )
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq

        bSq = Mp**2/(bt(j,k1)-Mp**2*(bt(j,k1)-1))

        Y = sqrt(((1.-bSq)*qqy)*((1.-bSq)*qqy)
     c     +4.*bSq*cc*cc)

        chky = 0.5 - sign(0.5,(bSq+1.)*qqy-Y)
        sp1 = 0.5*((bSq+1.)*abs(qqy)+Y)
        sm2 = eps2*sp1
        spec= chky*sp1 + vnu + sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        b(j,1) = a4                  + chky*qqy*s1
        b(j,2) = ri1*a2 + uu*a4      + chky*qqy*s2
        b(j,3) = ri2*a2 + vv*a4      + chky*qqy*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 + chky*qqy*s4

        do n=1,4
           tmp(n)=b(j,n)
        enddo
        
        b(j,1) =uv(j,k1)*tmp(1)-uu*tmp(2)-vv*tmp(3)+tmp(4)
        b(j,2) =uu*b(j,1)
        b(j,3) =vv*b(j,1)
        b(j,4) =ge(j,k1)*b(j,1)

        cjkl= sqrt(ggm1*(q(j,k1,4)-uv(j,k1)))
        term=gm1*(bSq-1.)/(cjkl*cjkl)
        do n=1,4
           b(j,n)=b(j,n)*term
           b(j,n)=b(j,n)+tmp(n)
        enddo

        b(j,1) = b(j,1) - spec*s1
        b(j,2) = b(j,2) - spec*s2
        b(j,3) = b(j,3) - spec*s3
        b(j,4) = b(j,4) - spec*s4
         
  222 continue
c
c..bi-diagonal inversion
c     
cdir$ ivdep
      do 223 j = ms(m),me(m)
        k  = m-j
        svt = tscale(j,k)
        sv = svt*0.5
        di = d(j,k)
        s(j,k,1) = s(j,k,1) - sv*( a(j,1)+b(j,1) )*di
        s(j,k,2) = s(j,k,2) - sv*( a(j,2)+b(j,2) )*di
        s(j,k,3) = s(j,k,3) - sv*( a(j,3)+b(j,3) )*di
        s(j,k,4) = s(j,k,4) - sv*( a(j,4)+b(j,4) )*di
  223 continue
c
  220 continue
 1000 continue
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c***********************************************************************
      subroutine ilu2d(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     >                 tscale)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables
      real,allocatable :: d(:,:),st(:,:,:),tmp(:)
      real,allocatable :: a(:,:),c(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)
      real,allocatable :: b(:,:)
      integer,allocatable :: ms(:),me(:)
      integer j,k,m,i,k1,jj
      integer itgs,j1
      real eps2,epsv,dj,uu,vv,uv2,cjkl,rj1,rj2
      real qq1,qqx,rr1,rk1,rk2,qq2,qqy,rr2,vnu,svt
      real ri1,ri2,qq,cc,sp1,sm2,chkx,spec,s1,s2
      real a2,a4,chky,sv,di,s3,s4
      real scale1,scale2,sav1,sav2,sav3,sav4,scale

      allocate(d(jd,kd),st(jd,kd,4),tmp(4))
      allocate(a(mdim,4),c(mdim,4))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))
      allocate(b(mdim,4))
      allocate(ms(mdim*2),me(mdim*2))

c***  first executable statement

      eps2 = epse*2.
      epsv = 1. + eps2

c..set-up for hyper-plane loop

      do 1 m=js+ks,je+ke
        i     = max(js,m-ke)
        ms(m) = i
        me(m) = min(je,m-js)
    1 continue

c
c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..setup d, the diagonal term and and store some arrays
c
      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
          uu  = q(j,k,2)
          vv  = q(j,k,3)
          uv2 = 0.5*( uu*uu + vv*vv )
          cjkl = sqrt( ggm1*( q(j,k,4) - uv2 ) )
c
          rj1 = xx(j,k)
          rj2 = xy(j,k)
          qq1 = rj1*uu + rj2*vv
          qqx = abs( -rj1*ug(j,k)-rj2*vg(j,k) + qq1 )
          rr1 = sqrt( rj1**2 + rj2**2 )
          rk1 = yx(j,k)
          rk2 = yy(j,k)
          qq2 = rk1*uu + rk2*vv
          qqy = abs( -rk1*ug(j,k)-rk2*vg(j,k) + qq2 )
          rr2 = sqrt( rk1**2 + rk2**2 )
          vnu = 2.0*rr2*rr2*(rmue+turmu(j,k))/(rey*q(j,k,1))
          svt = tscale(j,k)
          d(j,k)=1.0/
     &     ( 1.0+svt*((qqx+qqy +cjkl*(rr1+rr2))*epsv +vnu) )
c
          uv(j,k) = uv2
          qx(j,k) = qq1
          cx(j,k) = cjkl*rr1
          qy(j,k) = qq2
          cy(j,k) = cjkl*rr2
          vn(j,k) = vnu
          ge(j,k) = gamma*q(j,k,4) - gm1*uv2
  111 continue

c..apply point Jacobi update for the boundary residuals
c
      j = js-1
      do k = ks,ke
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

      j = je+1
      do k = ks,ke
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

      k = ks-1
      do j = js,je
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

      k = ke+1
      do j = js,je
        s(j,k,1) = s(j,k,1)*d(j,k)
        s(j,k,2) = s(j,k,2)*d(j,k)
        s(j,k,3) = s(j,k,3)*d(j,k)
        s(j,k,4) = s(j,k,4)*d(j,k)
      enddo

c..store RHS in separate variable
c
      do k = ks-1, ke+1
      do j = js-1, je+1
        st(j,k,1) = s(j,k,1)
        st(j,k,2) = s(j,k,2)
        st(j,k,3) = s(j,k,3)
        st(j,k,4) = s(j,k,4)
      enddo
      enddo
 
      do 1000 itgs=1,1 
c
c..forward sweep
c..loop on hyper-plane
c
      do 120 m = js+ks,je+ke
c
c..setup a contribution in j-direction
c
      do 121 j = ms(m),me(m)
        k  = m-j
        j1 = j-1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 + sign( 0.5, qqx+cc )
        spec= chkx*( qqx + sp1 ) + sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        a(j,1) = a4                  +spec*s1
        a(j,2) = ri1*a2 + uu*a4      +spec*s2
        a(j,3) = ri2*a2 + vv*a4      +spec*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 +spec*s4
  121 continue
c
c..setup b contribution in k-direction
c
      do 122 j = ms(m),me(m)
        k  = m-j
        k1 = k-1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq
        sp1 = abs(qqy) + cc
        sm2 = eps2*sp1
        chky= 0.5 + sign( 0.5, qqy+cc )
        spec= chky*( qqy + sp1 ) + vnu + sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        b(j,1) = a4                  +spec*s1
        b(j,2) = ri1*a2 + uu*a4      +spec*s2
        b(j,3) = ri2*a2 + vv*a4      +spec*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 +spec*s4
  122 continue
c
c..bi-diagonal inversion
c     
      if(iunst.eq.2) then
        do 123 j = ms(m),me(m)
          k  = m-j
          svt = tscale(j,k)
          sv = svt*0.5
          di = d(j,k)
          s(j,k,1) = ( st(j,k,1) + sv*( a(j,1)+b(j,1) ) )*di
          s2 = ( st(j,k,2) + sv*( a(j,2)+b(j,2) ) )*di
          s3 = ( st(j,k,3) + sv*( a(j,3)+b(j,3) ) )*di
          s4 = 1.0 + (2.*di*sv*rotf)**2
          s(j,k,2) = ( s2 + di*2.*sv*rotf*s3)/s4
          s(j,k,3) = ( s3 - di*2.*sv*rotf*s2)/s4
          s(j,k,4) = ( st(j,k,4) + sv*( a(j,4)+b(j,4) ) )*di
  123   continue
      else

cdir$ ivdep
        do 124 j = ms(m),me(m)
          k  = m-j
          svt = tscale(j,k)
          sv = svt*0.5
          di = d(j,k)
          s(j,k,1) = ( st(j,k,1) + sv*( a(j,1)+b(j,1) ) )*di
          s(j,k,2) = ( st(j,k,2) + sv*( a(j,2)+b(j,2) ) )*di
          s(j,k,3) = ( st(j,k,3) + sv*( a(j,3)+b(j,3) ) )*di
          s(j,k,4) = ( st(j,k,4) + sv*( a(j,4)+b(j,4) ) )*di
 124    continue
      endif
c
  120 continue
c
c..backward sweep
c..loop on hyper-plane
c
      do 220 m = je+ke,js+ks,-1
c
c..setup a contribution in j-direction
c
      do 221 j = ms(m),me(m)
        k  = m-j
        j1 = j+1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 - sign( 0.5, qqx-cc )
        spec= chkx*( qqx - sp1 ) - sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        a(j,1) = a4                  +spec*s1
        a(j,2) = ri1*a2 + uu*a4      +spec*s2
        a(j,3) = ri2*a2 + vv*a4      +spec*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 +spec*s4
  221 continue
c
c..setup b contribution in k-direction
c
      do 222 j = ms(m),me(m)
        k  = m-j
        k1 = k+1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq
        sp1 = abs(qqy) + cc
        sm2 = eps2*sp1
        chky= 0.5 - sign( 0.5, qqy-cc )
        spec= chky*( qqy - sp1 ) - vnu - sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        b(j,1) = a4                  +spec*s1
        b(j,2) = ri1*a2 + uu*a4      +spec*s2
        b(j,3) = ri2*a2 + vv*a4      +spec*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 +spec*s4
  222 continue
c
c..bi-diagonal inversion
c     
cdir$ ivdep
      do 223 j = ms(m),me(m)
        k  = m-j
        svt = tscale(j,k)
        sv = svt*0.5
        di = d(j,k)
        s(j,k,1) = s(j,k,1) - sv*( a(j,1)+b(j,1) )*di
        s(j,k,2) = s(j,k,2) - sv*( a(j,2)+b(j,2) )*di
        s(j,k,3) = s(j,k,3) - sv*( a(j,3)+b(j,3) )*di
        s(j,k,4) = s(j,k,4) - sv*( a(j,4)+b(j,4) )*di
  223 continue
c
  220 continue
 1000 continue
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c***********************************************************************
      subroutine arc2d(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     c                                                iblank,tscale,bt)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      integer n
      integer iblank(jd,kd)
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd),bt(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables

      real,allocatable :: p(:,:),d2px(:,:),d2py(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)

      integer j,k
      real dj,uu,vv,uv2,cjkl,rj1,rj2
      real qq1,rr1,rk1,rk2,qq2,rr2,vnu
      real c2i
      real a1,a2,a3,a4,a5,a6

      allocate(p(jd,kd),d2px(jd,kd),d2py(jd,kd))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))

c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..multiply by T_si_inverse

      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
         uu  = q(j,k,2)
         vv  = q(j,k,3)
         uv2 = 0.5*( uu*uu + vv*vv )
         cjkl   = sqrt( ggm1*( q(j,k,4) - uv2 ) )
         c2i    = 1.0/(cjkl*cjkl)
         p(j,k) = ( q(j,k,4) - uv2 )*gm1*q(j,k,1)
c
         rj1 = xx(j,k)
         rj2 = xy(j,k)
         qq1 = rj1*uu + rj2*vv
         qx(j,k) = qq1 - ( rj1*ug(j,k) + rj2*vg(j,k) )
         rr1 = sqrt( rj1**2 + rj2**2 )
         rj1 = rj1/rr1
         rj2 = rj2/rr1

         rk1 = yx(j,k)
         rk2 = yy(j,k)
         qq2 = rk1*uu + rk2*vv
         qy(j,k) = qq2 - ( rk1*ug(j,k) + rk2*vg(j,k) )
         rr2 = sqrt( rk1**2 + rk2**2 )
         rk1 = rk1/rr2
         rk2 = rk2/rr2

         vnu = (rmue+turmu(j,k))/(rey*q(j,k,1))

         uv(j,k) = uv2
         cx(j,k) = cjkl*rr1
         cy(j,k) = cjkl*rr2
         vn(j,k) = vnu
         ge(j,k) = gamma*q(j,k,4) - gm1*uv2

         a1 = s(j,k,2)*uu + s(j,k,3)*vv - s(j,k,4)
         a1 = a1*gm1*c2i
         a1 = a1+s(j,k,1)*( 1.0 - uv2*gm1*c2i )

         a2 =(rj1*vv-rj2*uu)*s(j,k,1)+rj2*s(j,k,2)-rj1*s(j,k,3)
  
         a3 = uv2*s(j,k,1)-s(j,k,2)*uu-s(j,k,3)*vv+s(j,k,4)
         a3 = a3*gm1*c2i

         a4 = qq1*s(j,k,1)/rr1-rj1*s(j,k,2)-rj2*s(j,k,3)
         a4 = a4

         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) = 0.5*(a3-a4/cjkl)
         s(j,k,4) = 0.5*(a3+a4/cjkl)

  111 continue

      do 112 k = ks-1,ke+1
      do 112 j = js,je
         d2px(j,k) = abs( p(j+1,k) - 2.*p(j,k) + p(j-1,k) )
         d2px(j,k) = d2px(j,k)/abs( p(j+1,k) + 2.*p(j,k) + p(j-1,k) )
  112 continue

      j=js-1
      do 113 k = ks-1,ke+1
         d2px(j,k) = d2px(j+1,k)
         d2px(j,k) = 0
  113 continue

      j=je+1
      do 114 k = ks-1,ke+1
         d2px(j,k) = d2px(j-1,k)
         d2px(j,k) = 0
  114 continue

c   finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,1,xx,xy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qx,cx,js-1,je+1,ks-1,ke+1,1,xx,
     >                  xy,vn,jd,kd,iblank,tscale,bt)
      endif
       

c multiply by N_inverse, s(j,k,1) is unchanged

      do 121 k = ks-1,ke+1
      do 121 j = js-1,je+1

         cjkl = sqrt( ggm1*( q(j,k,4) - uv(j,k) ) )
        
         rj1  = xx(j,k)
         rj2  = xy(j,k)
         rr1  = sqrt(rj1*rj1+rj2*rj2)
        
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         rr2  = sqrt(rk1*rk1+rk2*rk2)

         a1 = (rj1*rk1+rj2*rk2)/rr1/rr2
         a2 = (rj1*rk2-rj2*rk1)/rr1/rr2
         a3 = cjkl

         a4 = a1*s(j,k,2)+(cjkl*s(j,k,3)-cjkl*s(j,k,4))*a2
         a5 =-0.5*(a2*s(j,k,2)/cjkl-(a1+1)*s(j,k,3)+(a1-1)*s(j,k,4))
         a6 = 0.5*(a2*s(j,k,2)/cjkl-(a1-1)*s(j,k,3)+(a1+1)*s(j,k,4))


         s(j,k,2)=a4
         s(j,k,3)=a5
         s(j,k,4)=a6
  121 continue

      do 122 k = ks,ke
      do 122 j = js-1,je+1
         d2py(j,k) = abs( p(j,k+1) - 2.*p(j,k) + p(j,k-1) ) 
         d2py(j,k) = d2py(j,k)/abs( p(j,k+1) + 2.*p(j,k) + p(j,k-1) )
  122 continue

      k=ks-1
      do 123 j = js-1,je+1
         d2py(j,k) = d2py(j,k+1)
         d2py(j,k) = 0
  123 continue

      k=ke+1
      do 124 j = js-1,je+1
         d2py(j,k) = d2py(j,k-1)
         d2py(j,k) = 0
  124 continue

c...finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,2,yx,yy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qy,cy,ks-1,ke+1,js-1,je+1,2,yx,
     >                  yy,vn,jd,kd,iblank,tscale,bt)
      endif
  

c...multiply by T_eta
      do 131 k = ks-1,ke+1
      do 131 j = js-1,je+1
         uu   = q(j,k,2)
         vv   = q(j,k,3)
         cjkl = sqrt(ggm1*(q(j,k,4)-uv(j,k)))
         c2i  = 1.0/(cjkl*cjkl)
         
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         qq2  = qy(j,k) + ( rk1*ug(j,k) + rk2*vg(j,k) )
         rr2  = sqrt(rk1*rk1+rk2*rk2)
         rk1  = rk1/rr2
         rk2  = rk2/rr2
         
         a1 = s(j,k,1)+s(j,k,3)+s(j,k,4)

         a2 = uu*s(j,k,1)+rk2*s(j,k,2)+(uu+rk1*cjkl)*s(j,k,3)
         a2 = a2 + (uu-rk1*cjkl)*s(j,k,4)
 
         a3 = vv*s(j,k,1)-rk1*s(j,k,2)+(vv+rk2*cjkl)*s(j,k,3)
         a3 = a3 + (vv-rk2*cjkl)*s(j,k,4)

         a4 = uv(j,k)*s(j,k,1)+(rk2*uu-rk1*vv)*s(j,k,2)
         a4 = a4 + (ge(j,k)+cjkl*qq2/rr2)*s(j,k,3)
         a4 = a4 + (ge(j,k)-cjkl*qq2/rr2)*s(j,k,4)
         
         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) = a3
         s(j,k,4) = a4
  131 continue

      if (iunst.eq.2) then
         do k=ks-1,ke+1
            do j=js-1,je+1
              a2 = s(j,k,2)
              a3 = s(j,k,3)
              a4 = 1.0 + (tscale(j,k)*rotf)**2
              s(j,k,2) = ( a2 + tscale(j,k)*rotf*a3)/a4
              s(j,k,3) = ( a3 - tscale(j,k)*rotf*a2)/a4
            enddo
         enddo
      endif
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c***********************************************************************
      subroutine arc2d_precon(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     c                                                 iblank,tscale,bt)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      integer n
      integer iblank(jd,kd)
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd),bt(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables

      real,allocatable :: p(:,:),d2px(:,:),d2py(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)

      integer j,k
      real dj,uu,vv,uv2,cjkl,rj1,rj2
      real qqmxt,qq1,rr1,rk1,rk2,qq2,rr2,vnu
      real c2i,X,Y,Z,X1,Y1,Z1,bSq
      real a1,a2,a3,a4,a5,a6

      allocate(p(jd,kd),d2px(jd,kd),d2py(jd,kd))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))

c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..multiply by T_si_inverse

      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
         uu  = q(j,k,2)
         vv  = q(j,k,3)
         uv2 = 0.5*( uu*uu + vv*vv )
         cjkl   = sqrt( ggm1*( q(j,k,4) - uv2 ) )
         c2i    = 1.0/(cjkl*cjkl)
         p(j,k) = ( q(j,k,4) - uv2 )*gm1*q(j,k,1)
c
         rj1 = xx(j,k)
         rj2 = xy(j,k)
         qq1 = rj1*( uu - ug(j,k) ) + rj2*( vv - vg(j,k) )
         rr1 = sqrt( rj1**2 + rj2**2 )
         rj1 = rj1/rr1
         rj2 = rj2/rr1

         rk1 = yx(j,k)
         rk2 = yy(j,k)
         qq2 = rk1*( uu - ug(j,k) ) + rk2*( vv - vg(j,k) )
         rr2 = sqrt( rk1**2 + rk2**2 )
         rk1 = rk1/rr2
         rk2 = rk2/rr2

         vnu = (rmue+turmu(j,k))/(rey*q(j,k,1))

         uv(j,k) = uv2
         qx(j,k) = qq1
         cx(j,k) = cjkl*rr1
         qy(j,k) = qq2
         cy(j,k) = cjkl*rr2
         vn(j,k) = vnu
         ge(j,k) = gamma*q(j,k,4) - gm1*uv2

         bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

         X = sqrt(((1.-bSq)*qq1)*((1.-bSq)*qq1)
     c       +4.*bSq*cx(j,k)*cx(j,k))/rr1
         Y = 0.5*((1.0 - bSq)*qq1/rr1+X)
         Z = 0.5*((1.0 - bSq)*qq1/rr1-X)

         qqmxt = qq1 + ug(j,k)*rj1*rr1 + vg(j,k)*rj2*rr1

         a1 = s(j,k,2)*uu + s(j,k,3)*vv - s(j,k,4)
         a1 = a1*gm1*c2i
         a1 = a1+s(j,k,1)*( 1.0 - uv2*gm1*c2i )

         a2 =(rj1*vv-rj2*uu)*s(j,k,1)+rj2*s(j,k,2)-rj1*s(j,k,3)
  
         a3 = uv2*s(j,k,1)-s(j,k,2)*uu-s(j,k,3)*vv+s(j,k,4)
         a3 = a3*gm1*c2i

         a4 = qqmxt*s(j,k,1)/rr1-rj1*s(j,k,2)-rj2*s(j,k,3)
         a4 = a4*bSq

         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) =-(Z*a3+a4)/X
         s(j,k,4) = (Y*a3+a4)/X

  111 continue

      do 112 k = ks-1,ke+1
      do 112 j = js,je
         d2px(j,k) = abs( p(j+1,k) - 2.*p(j,k) + p(j-1,k) )
         d2px(j,k) = d2px(j,k)/abs( p(j+1,k) + 2.*p(j,k) + p(j-1,k) )
  112 continue

      j=js-1
      do 113 k = ks-1,ke+1
         d2px(j,k) = d2px(j+1,k)
         d2px(j,k) = 0
  113 continue

      j=je+1
      do 114 k = ks-1,ke+1
         d2px(j,k) = d2px(j-1,k)
         d2px(j,k) = 0
  114 continue

c   finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,1,xx,xy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qx,cx,js-1,je+1,ks-1,ke+1,1,xx,
     >                  xy,vn,jd,kd,iblank,tscale,bt)
      endif
       

c multiply by N_inverse, s(j,k,1) is unchanged

      do 121 k = ks-1,ke+1
      do 121 j = js-1,je+1

         cjkl = sqrt( ggm1*( q(j,k,4) - uv(j,k) ) )
        
         rj1  = xx(j,k)
         rj2  = xy(j,k)
         rr1  = sqrt(rj1*rj1+rj2*rj2)
        
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         rr2  = sqrt(rk1*rk1+rk2*rk2)

         bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

         X = sqrt((( 1.0 - bSq )*qx(j,k))*(( 1.0 - bSq )*qx(j,k))
     c      +4.*bSq*cx(j,k)*cx(j,k))/rr1
         Y = 0.5*(( 1.0 - bSq )*qx(j,k)/rr1+X)
         Z = 0.5*(( 1.0 - bSq )*qx(j,k)/rr1-X)

         X1 = sqrt((( 1.0 - bSq )*qy(j,k))*(( 1.0 - bSq )*qy(j,k))
     c       +4.*bSq*cy(j,k)*cy(j,k))/rr2
         Y1 = 0.5*(( 1.0 - bSq )*qy(j,k)/rr2+X1)
         Z1 = 0.5*(( 1.0 - bSq )*qy(j,k)/rr2-X1)

         a1 = (rj1*rk1+rj2*rk2)/rr1/rr2
         a2 = (rj1*rk2-rj2*rk1)/rr1/rr2
         a3 = cjkl

         a4 = a1*s(j,k,2)+(Y*s(j,k,3)+Z*s(j,k,4))*a2/bSq
         a5 =-(a2*bSq*s(j,k,2)-(Y*a1-Z1)*s(j,k,3)-(Z*a1-Z1)*s(j,k,4))/X1
         a6 = (a2*bSq*s(j,k,2)-(Y*a1-Y1)*s(j,k,3)-(Z*a1-Y1)*s(j,k,4))/X1


         s(j,k,2)=a4
         s(j,k,3)=a5
         s(j,k,4)=a6
  121 continue

      do 122 k = ks,ke
      do 122 j = js-1,je+1
         d2py(j,k) = abs( p(j,k+1) - 2.*p(j,k) + p(j,k-1) ) 
         d2py(j,k) = d2py(j,k)/abs( p(j,k+1) + 2.*p(j,k) + p(j,k-1) )
  122 continue

      k=ks-1
      do 123 j = js-1,je+1
         d2py(j,k) = d2py(j,k+1)
         d2py(j,k) = 0
  123 continue

      k=ke+1
      do 124 j = js-1,je+1
         d2py(j,k) = d2py(j,k-1)
         d2py(j,k) = 0
  124 continue

c...finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,2,yx,yy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qy,cy,ks-1,ke+1,js-1,je+1,2,yx,
     >                  yy,vn,jd,kd,iblank,tscale,bt)
      endif
  

c...multiply by T_eta
      do 131 k = ks-1,ke+1
      do 131 j = js-1,je+1
         uu   = q(j,k,2)
         vv   = q(j,k,3)
         cjkl = sqrt(ggm1*(q(j,k,4)-uv(j,k)))
         c2i  = 1.0/(cjkl*cjkl)
         
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         qq2  = qy(j,k)
         rr2  = sqrt(rk1*rk1+rk2*rk2)
         rk1  = rk1/rr2
         rk2  = rk2/rr2
         
         bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

         X = sqrt(((1.-bSq)*qq2)*((1.-bSq)*qq2)
     c       +4.*bSq*cy(j,k)*cy(j,k))/rr2
         Y = 0.5*((1.-bSq)*qq2/rr2+X)
         Z = 0.5*((1.-bSq)*qq2/rr2-X)
         
         qqmxt = qq2 + ug(j,k)*rk1*rr2 + vg(j,k)*rk2*rr2

         a1 = s(j,k,1)+s(j,k,3)+s(j,k,4)

         a2 = uu*s(j,k,1)+rk2*s(j,k,2)+(uu+rk1*Y/bSq)*s(j,k,3)
         a2 = a2 + (uu+rk1*Z/bSq)*s(j,k,4)
 
         a3 = vv*s(j,k,1)-rk1*s(j,k,2)+(vv+rk2*Y/bSq)*s(j,k,3)
         a3 = a3 + (vv+rk2*Z/bSq)*s(j,k,4)

         a4 = uv(j,k)*s(j,k,1)+(rk2*uu-rk1*vv)*s(j,k,2)
         a4 = a4 + (ge(j,k)+Y*qqmxt/rr2/bSq)*s(j,k,3)
         a4 = a4 + (ge(j,k)+Z*qqmxt/rr2/bSq)*s(j,k,4)
         
         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) = a3
         s(j,k,4) = a4
  131 continue

      if (iunst.eq.2) then
         do k=ks-1,ke+1
            do j=js-1,je+1
              a2 = s(j,k,2)
              a3 = s(j,k,3)
              a4 = 1.0 + (tscale(j,k)*rotf)**2
              s(j,k,2) = ( a2 + tscale(j,k)*rotf*a3)/a4
              s(j,k,3) = ( a3 - tscale(j,k)*rotf*a2)/a4
            enddo
         enddo
      endif
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c**********************************************************************
      subroutine lhsinv_up(q,s,qn,cn,ms,me,ns,ne,idir,xn,yn,vnu,
     >                     jd,kd,iblank,tscale,bt)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer ms,me,ns,ne,idir,jd,kd
      integer iblank(jd,kd)
      real q(jd,kd,nq),s(jd,kd,nv),tscale(jd,kd),bt(jd,kd)
      real qn(jd,kd),cn(jd,kd),xn(jd,kd),yn(jd,kd),vnu(jd,kd)

      !local variables
      integer,allocatable :: iblnk(:)
      real,allocatable :: jcbn(:),tscal(:)
      real,allocatable :: vistrm1(:),vistrm2(:),vistrm3(:),g(:)
      real,allocatable :: a(:),b(:),c(:),d(:),e(:),f(:)
      real,allocatable :: a1(:),b1(:),c1(:),d1(:),e1(:),f1(:)
      real,allocatable :: diag_plus(:,:),diag_minus(:,:)

      integer i,m,n
      real svt,dis2,dis4,bSq,X
      real c2,c2m,c4m2,c4m,c4,c4p
      real eig1,eig2,eig3,epsval,eps,fac

      allocate(iblnk(me))
      allocate(jcbn(me),tscal(me))
      allocate(vistrm1(me),vistrm2(me),vistrm3(me),g(me))
      allocate(a(me),b(me),c(me),d(me),e(me),f(me))
      allocate(a1(me),b1(me),c1(me),d1(me),e1(me),f1(me))
      allocate(diag_plus(me,nq),diag_minus(me,nq))

c..some initialization

      epsval  = 0.05
      fac     = 1.05

      do n = ns,ne
        if (idir.eq.1) then
          do m = ms,me
            tscal(m) = tscale(m,n)
            jcbn(m)  = q(m,n,nq)
            iblnk(m) = max(iblank(m,n),0)
            bSq = Mp**2/(bt(m,n)-Mp**2*(bt(m,n)-1))
            X = sqrt(((1.-bSq)*qn(m,n))*((1.-bSq)*qn(m,n))
     c         +4.*bSq*cn(m,n)*cn(m,n))
            eig1 = qn(m,n)
            eig2 = 0.5*((bSq+1.)*qn(m,n)+X)
            eig3 = 0.5*((bSq+1.)*qn(m,n)-X)
            eps  = epsval*sqrt( xn(m,n)**2 + yn(m,n)**2 )
            diag_plus(m,1)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
            diag_plus(m,2)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
            diag_plus(m,3)  = 0.5*(eig2 + fac*sqrt(eig2**2 + eps**2))
            diag_plus(m,4)  = 0.5*(eig3 + fac*sqrt(eig3**2 + eps**2))
            diag_minus(m,1) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
            diag_minus(m,2) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
            diag_minus(m,3) = 0.5*(eig2 - fac*sqrt(eig2**2 + eps**2))
            diag_minus(m,4) = 0.5*(eig3 - fac*sqrt(eig3**2 + eps**2))
          enddo

          do m = ms+1,me-1
            vistrm1(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
            vistrm1(m) = vistrm1(m) + (xn(m+1,n)*xn(m+1,n)
     c                              + yn(m+1,n)*yn(m+1,n))/q(m+1,n,nq)
            vistrm1(m) = vistrm1(m)*q(m,n,nq)*vnu(m,n)

            vistrm3(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
            vistrm3(m) = vistrm3(m) + (xn(m-1,n)*xn(m-1,n)
     c                              + yn(m-1,n)*yn(m-1,n))/q(m-1,n,nq)
            vistrm3(m)= vistrm3(m)*q(m,n,nq)*vnu(m,n)

            vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
          enddo
          m = ms
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0
          m = me
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0

        elseif (idir.eq.2) then

          do m = ms,me
            tscal(m) = tscale(n,m)
            jcbn(m)  = q(n,m,nq)
            iblnk(m) = max(iblank(n,m),0)
            bSq = Mp**2/(bt(n,m)-Mp**2*(bt(n,m)-1))
            X= sqrt(((1.-bSq)*qn(n,m))*((1.-bSq)*qn(n,m))
     c         +4.*bSq*cn(n,m)*cn(n,m))
            eig1 = qn(n,m)
            eig2 = 0.5*((bSq+1.)*qn(n,m)+X)
            eig3 = 0.5*((bSq+1.)*qn(n,m)-X)
            eps  = epsval*sqrt( xn(n,m)**2 + yn(n,m)**2 )
            diag_plus(m,1)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
            diag_plus(m,2)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
            diag_plus(m,3)  = 0.5*(eig2 + fac*sqrt(eig2**2 + eps**2))
            diag_plus(m,4)  = 0.5*(eig3 + fac*sqrt(eig3**2 + eps**2))
            diag_minus(m,1) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
            diag_minus(m,2) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
            diag_minus(m,3) = 0.5*(eig2 - fac*sqrt(eig2**2 + eps**2))
            diag_minus(m,4) = 0.5*(eig3 - fac*sqrt(eig3**2 + eps**2))
          enddo

          do m = ms+1,me-1
            vistrm1(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
            vistrm1(m) = vistrm1(m) + (xn(n,m+1)*xn(n,m+1)
     c                              + yn(n,m+1)*yn(n,m+1))/q(n,m+1,nq)
            vistrm1(m) = vistrm1(m)*q(n,m,nq)*vnu(n,m)
         
            vistrm3(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
            vistrm3(m) = vistrm3(m) + (xn(n,m-1)*xn(n,m-1)
     c                              + yn(n,m-1)*yn(n,m-1))/q(n,m-1,nq)
            vistrm3(m) = vistrm3(m)*q(n,m,nq)*vnu(n,m)

            vistrm2(m) = 0.5*(vistrm1(m)+vistrm3(m))
          enddo
          m = ms
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0
          m = me
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0
 
        endif
 
c
        do i = 1,nmv
          do m=ms,me
            a1(m) = 0.
            b1(m) = 0.
            c1(m) = 0.
            d1(m) = 0.
            e1(m) = 0.
            f1(m) = 0.

            if (idir.eq.1) then
               f1(m) = s(m,n,i)
            elseif (idir.eq.2) then
               f1(m) = s(n,m,i)
            endif

c        euler terms
            b1(m) = b1(m) - diag_plus(m,i)
            d1(m) = d1(m) + diag_minus(m,i)
            c1(m) = c1(m) + diag_plus(m,i) - diag_minus(m,i)

c        viscous terms
            b1(m) = b1(m) - 0.5*vistrm1(m)
            c1(m) = c1(m) + vistrm2(m)
            d1(m) = d1(m) - 0.5*vistrm3(m)

c        multiply by delta t
            if(m.le.me-1) then
              b1(m) = b1(m)*tscal(m+1)
            else
              b1(m) = b1(m)*tscal(m)
            endif
            if(m.ge.ms+1) then
              d1(m) = d1(m)*tscal(m-1)
            else
              d1(m) = d1(m)*tscal(m)
            endif
            c1(m) = c1(m)*tscal(m) + 1.0
          enddo

          do m = ms+1,me-1
            a(m-ms) = a1(m) 
            b(m-ms) = b1(m) 
            c(m-ms) = c1(m) 
            d(m-ms) = d1(m) 
            e(m-ms) = e1(m) 
            f(m-ms) = f1(m) 
          enddo
          f(1) = f(1) - b1(ms)*f1(ms)/c1(ms)
          f(me-ms-1) = f(me-ms-1) - d1(me)*f1(me)/c1(me)
          d(1) = 0
          e(1) = 0
          e(2) = 0
          b(me-ms-1) = 0
          a(me-ms-1) = 0
          a(me-ms-2) = 0

          call pentadag(a,b,c,d,e,f,me-ms-1)

          do m=ms+1,me-1
            if (idir.eq.1) then
              s(m,n,i)=f(m-ms)
            elseif (idir.eq.2) then
              s(n,m,i)=f(m-ms)
            endif
          enddo
        enddo
      enddo

      return
      end

c**********************************************************************
      subroutine lhsinv(q,s,qn,cn,d2p,ms,me,ns,ne,idir,xn,yn,vnu,jd,kd,
     >                                                iblank,tscale,bt)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer ms,me,ns,ne,idir,jd,kd
      integer iblank(jd,kd)
      real q(jd,kd,nq),s(jd,kd,nv),tscale(jd,kd),bt(jd,kd)
      real qn(jd,kd),cn(jd,kd),d2p(jd,kd),xn(jd,kd),yn(jd,kd),vnu(jd,kd)

      !local variables
      integer,allocatable :: iblnk(:) 
      real,allocatable :: jcbn(:),tscal(:)
      real,allocatable :: vistrm1(:),vistrm2(:),vistrm3(:),g(:)
      real,allocatable :: a(:),b(:),c(:),d(:),e(:),f(:)
      real,allocatable :: a1(:),b1(:),c1(:),d1(:),e1(:),f1(:)
      real,allocatable :: atmp(:),btmp(:),ctmp(:),dtmp(:),etmp(:)
      real,allocatable :: diag(:,:)

      integer i,m,n
      real svt,dis2,dis4,X,bSq
      real c2,c2m,c4m2,c4m,c4,c4p

      allocate(iblnk(me))
      allocate(jcbn(me),tscal(me))
      allocate(vistrm1(me),vistrm2(me),vistrm3(me),g(me))
      allocate(a(me),b(me),c(me),d(me),e(me),f(me))
      allocate(a1(me),b1(me),c1(me),d1(me),e1(me),f1(me))
      allocate(atmp(me),btmp(me),ctmp(me),dtmp(me),etmp(me))
      allocate(diag(me,nq))

c..some initialization

      dis2 = 10.0
      dis4 = 0.1

      do n = ns,ne
        if (idir.eq.1) then
          do m = ms,me
            jcbn(m)  = q(m,n,nq)
            tscal(m) = tscale(m,n)
            g(m)     = d2p(m,n)
            iblnk(m) = max(iblank(m,n),0)
            bSq = Mp**2/(bt(m,n)-Mp**2*(bt(m,n)-1))
            X= sqrt(((1.-bSq)*qn(m,n))*((1.-bSq)*qn(m,n))
     c         +4.*bSq*cn(m,n)*cn(m,n))
            diag(m,1) = qn(m,n)
            diag(m,2) = qn(m,n)
            diag(m,3) = 0.5*((bSq+1.)*qn(m,n)+X)
            diag(m,4) = 0.5*((bSq+1.)*qn(m,n)-X)
c           spectral radius divided by jacobian
            diag(m,nq)= 0.5*((bSq+1.)*abs(qn(m,n))+X)/q(m,n,nq)
          enddo

          do m = ms+1,me-1
            vistrm1(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
            vistrm1(m) = vistrm1(m) + (xn(m+1,n)*xn(m+1,n)
     c                              +yn(m+1,n)*yn(m+1,n))/q(m+1,n,nq)
            vistrm1(m) = vistrm1(m)*q(m,n,nq)*vnu(m,n)

            vistrm3(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
            vistrm3(m) = vistrm3(m) + (xn(m-1,n)*xn(m-1,n)
     c                              +yn(m-1,n)*yn(m-1,n))/q(m-1,n,nq)
            vistrm3(m)= vistrm3(m)*q(m,n,nq)*vnu(m,n)

            vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
          enddo
          m = ms
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0
          m = me
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0

        elseif (idir.eq.2) then

          do m = ms,me
            jcbn(m) = q(n,m,nq)
            tscal(m) = tscale(n,m)
            g(m)= d2p(n,m)
            iblnk(m) = max(iblank(n,m),0)
            bSq = Mp**2/(bt(n,m)-Mp**2*(bt(n,m)-1))
            X= sqrt(((1.-bSq)*qn(n,m))*((1.-bSq)*qn(n,m))
     c         +4.*bSq*cn(n,m)*cn(n,m))
            diag(m,1)= qn(n,m)
            diag(m,2)= qn(n,m)
            diag(m,3)= 0.5*((bSq+1.)*qn(n,m)+X)
            diag(m,4)= 0.5*((bSq+1.)*qn(n,m)-X)
            diag(m,nq)= 0.5*((bSq+1.)*abs(qn(n,m))+X)/q(n,m,nq)
          enddo

          do m = ms+1,me-1
            vistrm1(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
            vistrm1(m) = vistrm1(m) + (xn(n,m+1)*xn(n,m+1)
     c                              +yn(n,m+1)*yn(n,m+1))/q(n,m+1,nq)
            vistrm1(m) = vistrm1(m)*q(n,m,nq)*vnu(n,m)
        
            vistrm3(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
            vistrm3(m) = vistrm3(m) + (xn(n,m-1)*xn(n,m-1)
     c                              +yn(n,m-1)*yn(n,m-1))/q(n,m-1,nq)
            vistrm3(m) = vistrm3(m)*q(n,m,nq)*vnu(n,m)

            vistrm2(m) = 0.5*(vistrm1(m)+vistrm3(m))
          enddo
          m = ms
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0
          m = me
          vistrm1(m) = 0
          vistrm2(m) = 0
          vistrm3(m) = 0
 
        endif

 
        do m = ms+3,me-3
c         second order dissipation
          c2   = 0.5*((g(m)+g(m+1)))*dis2
          c2m  = 0.5*((g(m)+g(m-1)))*dis2

          c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
          c4   = (dis4-min(dis4,1.*c2))
          c4m  = (dis4-min(dis4,1.*c2m))
          c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))

          atmp(m) = (diag(m+1,nq)+diag(m+2,nq))*c4p
          btmp(m) = -((diag(m+2,nq)+diag(m+1,nq))*c4p+
     c             (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
          ctmp(m) = +((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c             (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
          dtmp(m) = -((diag(m-2,nq)+diag(m-1,nq))*c4m2+
     c             (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
          etmp(m) = (diag(m-1,nq)+diag(m-2,nq))*c4m2
        enddo 

        m    = ms+2
        c2   = 0.5*((g(m)+g(m+1)))*dis2
        c2m  = 0.5*((g(m)+g(m-1)))*dis2
        c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
        c4   = (dis4-1.*min(dis4,c2))
        c4m  = (dis4-min(dis4,1.*c2m))
        c4m2 = (dis4-min(dis4,1.*dis2*(g(m-1))))
        atmp(m) = (diag(m+1,nq)+diag(m+2,nq))*c4p
        btmp(m) = -((diag(m+2,nq)+diag(m+1,nq))*c4p+
     c           (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
        ctmp(m) = ((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c          (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
        dtmp(m) = -((2.*diag(m-1,nq))*c4m2+
     c           (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
        
        m    = ms+1
        c2   = 0.5*((g(m)+g(m+1)))*dis2
        c2m  = ((g(m)))*dis2
        c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
        c4   = (dis4-min(dis4,1.*c2))
        c4m  = (dis4-min(dis4,1.*c2m))
        atmp(m) = (diag(m+1,nq)+diag(m+2,nq))*c4p
        btmp(m) = -((diag(m+2,nq)+diag(m+1,nq))*c4p+
     c           (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
        ctmp(m) = (2.*(diag(m,nq))*(2.*c4m+c2m)+
     c           (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
        m    = me-2
        c2   = 0.5*((g(m)+g(m+1)))*dis2
        c2m  = 0.5*((g(m)+g(m-1)))*dis2
        c4p  = (dis4-min(dis4,1.*dis2*(g(m+1))))
        c4   = (dis4-min(dis4,1.*c2))
        c4m  = (dis4-min(dis4,1.*c2m))
        c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))
        btmp(m) = -(2.*(diag(m+1,nq))*c4p+
     c           (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
        ctmp(m) = ((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c          (diag(m+1,nq)+diag(m,nq))*(2.*c4+c2))
        dtmp(m) = -((diag(m-2,nq)+diag(m-1,nq))*c4m2+
     c           (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
        etmp(m) = (diag(m-1,nq)+diag(m-2,nq))*c4m2
        
        m    = me-1  
        c2   = ((g(m)))*dis2
        c2m  = 0.5*((g(m)+g(m-1)))*dis2
        c4   = (dis4-min(dis4,1.*c2))
        c4m  = (dis4-min(dis4,1.*c2m))
        c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))
        ctmp(m) = +((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c           2.*(diag(m,nq))*(2.*c4+c2))
        dtmp(m) = -((diag(m-2,nq)+diag(m-1,nq))*c4m2+
     c           (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
        etmp(m) = (diag(m-1,nq)+diag(m-2,nq))*c4m2

        do m = ms,me
          atmp(m) = atmp(m)*jcbn(m)
          btmp(m) = btmp(m)*jcbn(m)
          ctmp(m) = ctmp(m)*jcbn(m)
          dtmp(m) = dtmp(m)*jcbn(m)
          etmp(m) = etmp(m)*jcbn(m)
        enddo

        do i = 1,nmv
          do m = ms,me
            a1(m) = atmp(m)
            b1(m) = btmp(m)
            c1(m) = ctmp(m)
            d1(m) = dtmp(m)
            e1(m) = etmp(m)
            if (idir.eq.1) then
               f1(m) = s(m,n,i)
            elseif (idir.eq.2) then
               f1(m) = s(n,m,i)
            endif

c        euler terms
            b1(m) = b1(m) - diag(m,i)*0.5
            d1(m) = d1(m) + diag(m,i)*0.5

c        viscous terms
            b1(m) = b1(m) - 0.5*vistrm1(m)
            c1(m) = c1(m) + vistrm2(m)
            d1(m) = d1(m) - 0.5*vistrm3(m)

c        multiply by delta t
            if(m.le.me-2) then
              a1(m) = a1(m)*tscal(m+2)
            else
              a1(m) = a1(m)*tscal(m)
            endif
            if(m.le.me-1) then
              b1(m) = b1(m)*tscal(m+1)
            else
              b1(m) = b1(m)*tscal(m)
            endif
            if(m.ge.ms+1) then
              d1(m) = d1(m)*tscal(m-1)
            else
              d1(m) = d1(m)*tscal(m)
            endif
            if(m.ge.ms+2) then
              e1(m) = e1(m)*tscal(m-2)
            else
              e1(m) = e1(m)*tscal(m)
            endif

c         add identity to matrix
            c1(m) = c1(m)*tscal(m) + 1.0
          enddo

          do m = ms+1,me-1
            a(m-ms) = a1(m) 
            b(m-ms) = b1(m) 
            c(m-ms) = c1(m) 
            d(m-ms) = d1(m) 
            e(m-ms) = e1(m) 
            f(m-ms) = f1(m) 
          enddo
          f(1) = f(1) - b1(ms)*f1(ms)/c1(ms)
          f(me-ms-1) = f(me-ms-1) - d1(me)*f1(me)/c1(me)
          d(1) = 0
          e(1) = 0
          e(2) = 0
          b(me-ms-1) = 0
          a(me-ms-1) = 0
          a(me-ms-2) = 0

          call pentadag(a,b,c,d,e,f,me-ms-1)

          do m=ms+1,me-1
            if (idir.eq.1) then
              s(m,n,i)=f(m-ms)
            elseif (idir.eq.2) then
              s(n,m,i)=f(m-ms)
            endif
          enddo
        enddo
      enddo 

      return
      end

c**********************************************************************
      subroutine pentadag(a,b,c,d,e,f,N)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer N
      real a(N),b(N),c(N),d(N),e(N),f(N)
 
      !local variables
      integer i
      real d0,d1,d2

c     making diagonal 1.
      i = 1
      d0 = 1.0/c(i)
      d(i+1) = d(i+1)*d0
      e(i+2) = e(i+2)*d0
      f(i)   = f(i)*d0
      
c     making diagonal 1. and lower diagonal 0.
      i = 2
      d1 = b(i-1)
      d0 = 1.0/(c(i)-d1*d(i))
      f(i)   = (f(i)-d1*f(i-1))*d0
      d(i+1) = (d(i+1)-d1*e(i+1))*d0
      e(i+2) = e(i+2)*d0
      
c     making diaongal 1. and lower diagonals 0.
      do i = 3,N
         d1 = a(i-2)
         d2 = (b(i-1) - a(i-2)*d(i-1))
         d0 = 1./(c(i)-d2*d(i)-d1*e(i))
         f(i) = (f(i)-d2*f(i-1)-d1*f(i-2))*d0
         if (i.le.N-1) then
            d(i+1) = (d(i+1)-d2*e(i+1))*d0
            if (i.le.n-2) then
               e(i+2) = (e(i+2))*d0
            endif
         endif
      enddo 

c     backward sweep
      i  = N-1
      d1 = d(i+1)
      f(i) = f(i) - d1*f(i+1)
      do i = N-2,1,-1
         d1 = d(i+1)
         d2 = e(i+2)
         f(i) = f(i) - d1*f(i+1) - d2*f(i+2)
      enddo

      return
      end
c**********************************************************************

