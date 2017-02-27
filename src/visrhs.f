c***********************************************************************
      subroutine visrhs(turmu,q,s,xx,xy,yx,yy,ug,vg,
     >                  tscale,iblank,jd,kd)
c
c  compute the viscous rhs 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd), tscale(jd,kd)
      real ug(jd,kd),vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      integer iblank(jd,kd)
c
      ! local variables
      real :: f2(mdim),f3(mdim),f4(mdim)
      real :: f2v(mdim),f3v(mdim),f4v(mdim),u(mdim),v(mdim),
     <         e(mdim)
      real :: a1v(mdim),a2v(mdim),a4v(mdim),a5v(mdim)
      real :: b1v(mdim),b2v(mdim),b4v(mdim),b5v(mdim),d1v(mdim)
      real :: d2v(mdim),d4v(mdim),d5v(mdim),d7v(mdim)


      real gkpr,prtr,dr,c2b,c2bp,t4,t5,t6,t7,t8
      integer j,k,km1,k1,jm1,j1,je,ke
      real ra,tt,vmue,turm,vnu,gkap,rj,b1,b2,b5,b4
      real du,dv,dei,rjm,rjp,d1,d2,d4,d5,d7
      real ujp1,ujm1,dvj
      real ejm1,dej,rk,rkm,rkp,a1,a2,a4,a5
      real ukp1,ukm1,duk
      real vkp1,vkm1,dvk,ekp1,ekm1,dek
      real dre,duj,vjp1,vjm1,ejp1

c*** first executable statement

      gkpr = gamma/pr
      prtr = pr/0.9
      dre  = .5/rey
      c2b  =198.6/tinf
      c2bp = c2b +1.

C$AD II-LOOP
      do j = jbeg,jend
C$AD II-LOOP
        do k = 1,kd
          ra    = 1./q(j,k,1)
          u(k)  = q(j,k,2)*ra
          v(k)  = q(j,k,3)*ra
          e(k)  = q(j,k,4)*ra-.5*(u(k)**2+v(k)**2)
          tt    = ggm1*e(k)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,nq)
          rjm    = 1./q(j-1,k,nq)
          rjp    = 1./q(j+1,k,nq)

          b4v(k) = (yx(j,k)**2+yy(j,k)**2)*rj
          b1v(k) = (b4v(k)+yx(j,k)**2/3.*rj)*vnu*dre
          b2v(k) = (b4v(k)+yy(j,k)**2/3.*rj)*vnu*dre
          b5v(k) = (yx(j,k)*yy(j,k)/3.*rj)*vnu*dre
          b4v(k) = b4v(k)*gkpr*gkap*dre

          d1v(k) = (4./3.*xx(j-1,k)*yx(j-1,k)+xy(j-1,k)*yy(j-1,k))*rjm
          d1v(k) = d1v(k)+(4./3.*xx(j+1,k)*yx(j+1,k)+xy(j+1,k)*
     1    yy(j+1,k))*rjp
          d1v(k) = vnu*dre*d1v(k)

          d2v(k) = (4./3.*xy(j-1,k)*yy(j-1,k)+xx(j-1,k)*yx(j-1,k))*rjm
          d2v(k) = d2v(k)+(4./3.*xy(j+1,k)*yy(j+1,k)+xx(j+1,k)*
     1    yx(j+1,k))*rjp
          d2v(k) = vnu*dre*d2v(k)

          d4v(k) = (xy(j-1,k)*yy(j-1,k)+xx(j-1,k)*yx(j-1,k))*rjm
          d4v(k) = d4v(k)+(xy(j+1,k)*yy(j+1,k)+xx(j+1,k)*
     1    yx(j+1,k))*rjp
          d4v(k) = d4v(k)*gkpr*gkap*dre

          d5v(k) = (-2./3.*xy(j-1,k)*yx(j-1,k)+xx(j-1,k)*yy(j-1,k))*rjm
          d5v(k) = d5v(k)+(-2./3.*xy(j+1,k)*yx(j+1,k)+xx(j+1,k)*
     1    yy(j+1,k))*rjp
          d5v(k) = vnu*dre*d5v(k)

          d7v(k) = (-2./3.*xx(j-1,k)*yy(j-1,k)+xy(j-1,k)*yx(j-1,k))*rjm
          d7v(k) = d7v(k)+(-2./3.*xx(j+1,k)*yy(j+1,k)+xy(j+1,k)*
     1    yx(j+1,k))*rjp
          d7v(k) = vnu*dre*d7v(k)
        enddo

C$AD II-LOOP
        do  k = 1,kd-1
          k1    = k+1

          b1    = b1v(k1)+b1v(k)
          b2    = b2v(k1)+b2v(k)
          b5    = b5v(k1)+b5v(k)
          b4    = b4v(k1)+b4v(k)
          d1    = d1v(k)
          d2    = d2v(k)
          d4    = d4v(k)
          d5    = d5v(k)
          d7    = d7v(k)
          ujp1  = q(j+1,k,2)/q(j+1,k,1)
          ujm1  = q(j-1,k,2)/q(j-1,k,1)
          duj   = 0.5*(ujp1-ujm1)

          vjp1  = q(j+1,k,3)/q(j+1,k,1)
          vjm1  = q(j-1,k,3)/q(j-1,k,1)
          dvj   = 0.5*(vjp1-vjm1)

          ejp1  = q(j+1,k,4)/q(j+1,k,1)-.5*(ujp1**2+vjp1**2)
          ejm1  = q(j-1,k,4)/q(j-1,k,1)-.5*(ujm1**2+vjm1**2)
          dej   = 0.5*(ejp1-ejm1)

          t4    = u(k1)*b1v(k1)+u(k)*b1v(k)
          t5    = v(k1)*b2v(k1)+v(k)*b2v(k)
          t6    = u(k1)*b5v(k1)+u(k)*b5v(k)
          t7    = v(k1)*b5v(k1)+v(k)*b5v(k)
          t8    = b4v(k1)+b4v(k)


          du    = u(k1)-u(k)
          dv    = v(k1)-v(k)
          dei   = e(k1)-e(k)
          f2(k) = b1*du+b5*dv
          f2v(k)= d1*duj+d5*dvj
          f3(k) = b5*du+b2*dv
          f3v(k)= d7*duj+d2*dvj
          f4(k) = (t4+t7)*du+(t5+t6)*dv+t8*dei
          f4v(k) = 0.5*d1*0.5*(ujp1**2-ujm1**2)
          f4v(k) = f4v(k)+0.5*d2*0.5*(vjp1**2-vjm1**2)
          f4v(k) = f4v(k)+d5*0.5*(ujp1+ujm1)*dvj
          f4v(k) = f4v(k)+d7*0.5*(vjp1+vjm1)*duj
          f4v(k) = f4v(k)+d4*dej
        enddo

        if (nhalo.eq.1) then
          ke = kend - 1

          k=kend
          s(j,k,2) = s(j,k,2)+(f2(k)-f2(k-1))+0.5*(f2v(k-2)
     1          -4.*f2v(k-1)+3*f2v(k))
          s(j,k,3) = s(j,k,3)+(f3(k)-f3(k-1))+0.5*(f3v(k-2)
     1          -4.*f3v(k-1)+3*f3v(k))
          s(j,k,4) = s(j,k,4)+(f4(k)-f4(k-1))+0.5*(f4v(k-2)
     1          -4.*f4v(k-1)+3.*f4v(k))
        elseif (nhalo.gt.1) then
          ke = kend 
        endif
c
C$AD II-LOOP
        do k = kbeg,ke
          s(j,k,2) = s(j,k,2)+(f2(k)-f2(k-1))+0.5*(f2v(k+1)-f2v(k-1))
          s(j,k,3) = s(j,k,3)+(f3(k)-f3(k-1))+0.5*(f3v(k+1)-f3v(k-1))
          s(j,k,4) = s(j,k,4)+(f4(k)-f4(k-1))+0.5*(f4v(k+1)-f4v(k-1))
        enddo

      enddo

c
C$AD II-LOOP
      do k = kbeg,kend
C$AD II-LOOP
        do j = 1,jd
          ra    = 1./q(j,k,1)
          u(j)  = q(j,k,2)*ra
          v(j)  = q(j,k,3)*ra
          e(j)  = q(j,k,4)*ra-.5*(u(j)**2+v(j)**2)
          tt    = ggm1*e(j)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,nq)
          rkm    = 1./q(j,k-1,nq)
          rkp    = 1./q(j,k+1,nq)

          a4v(j) = (xx(j,k)**2+xy(j,k)**2)*rj
          a1v(j) = (a4v(j)+xx(j,k)**2/3.*rj)*vnu*dre
          a2v(j) = (a4v(j)+xy(j,k)**2/3.*rj)*vnu*dre
          a5v(j) = (xx(j,k)*xy(j,k)/3.*rj)*vnu*dre
          a4v(j) = a4v(j)*gkpr*gkap*dre

          d1v(j) = (4./3.*xx(j,k-1)*yx(j,k-1)+xy(j,k-1)*yy(j,k-1))*rkm
          d1v(j) = d1v(j)+(4./3.*xx(j,k+1)*yx(j,k+1)+xy(j,k+1)*
     1    yy(j,k+1))*rkp
          d1v(j) = vnu*dre*d1v(j)

          d2v(j) = (4./3.*xy(j,k-1)*yy(j,k-1)+xx(j,k-1)*yx(j,k-1))*rkm
          d2v(j) = d2v(j)+(4./3.*xy(j,k+1)*yy(j,k+1)+xx(j,k+1)*
     1    yx(j,k+1))*rkp
          d2v(j) = vnu*dre*d2v(j)

          d4v(j) = (xy(j,k-1)*yy(j,k-1)+xx(j,k-1)*yx(j,k-1))*rkm
          d4v(j) = d4v(j)+(xy(j,k+1)*yy(j,k+1)+xx(j,k+1)*
     1    yx(j,k+1))*rkp
          d4v(j) = d4v(j)*gkpr*gkap*dre

          d5v(j) = (-2./3.*xy(j,k-1)*yx(j,k-1)+xx(j,k-1)*yy(j,k-1))*rkm
          d5v(j) = d5v(j)+(-2./3.*xy(j,k+1)*yx(j,k+1)+xx(j,k+1)*
     1    yy(j,k+1))*rkp
          d5v(j) = vnu*dre*d5v(j)

          d7v(j) = (-2./3.*xx(j,k-1)*yy(j,k-1)+xy(j,k-1)*yx(j,k-1))*rkm
          d7v(j) = d7v(j)+(-2./3.*xx(j,k+1)*yy(j,k+1)+xy(j,k+1)*
     1    yx(j,k+1))*rkp
          d7v(j) = vnu*dre*d7v(j)
        enddo

C$AD II-LOOP
        do j = 1,jd-1
          j1    = j+1

          a1    = a1v(j1)+a1v(j)
          a2    = a2v(j1)+a2v(j)
          a5    = a5v(j1)+a5v(j)
          a4    = a4v(j1)+a4v(j)
          d1    = d1v(j)
          d2    = d2v(j)
          d4    = d4v(j)
          d5    = d5v(j)
          d7    = d7v(j)

          ukp1  = q(j,k+1,2)/q(j,k+1,1)
          ukm1  = q(j,k-1,2)/q(j,k-1,1)
          duk   = 0.5*(ukp1-ukm1)

          vkp1  = q(j,k+1,3)/q(j,k+1,1)
          vkm1  = q(j,k-1,3)/q(j,k-1,1)
          dvk   = 0.5*(vkp1-vkm1)

          ekp1  = q(j,k+1,4)/q(j,k+1,1)-.5*(ukp1**2+vkp1**2)
          ekm1  = q(j,k-1,4)/q(j,k-1,1)-.5*(ukm1**2+vkm1**2)
          dek   = 0.5*(ekp1-ekm1)

          t4    = u(j1)*a1v(j1)+u(j)*a1v(j)
          t5    = v(j1)*a2v(j1)+v(j)*a2v(j)
          t6    = u(j1)*a5v(j1)+u(j)*a5v(j)
          t7    = v(j1)*a5v(j1)+v(j)*a5v(j)
          t8    = a4v(j1)+a4v(j)


          du    = u(j1)-u(j)
          dv    = v(j1)-v(j)
          dei   = e(j1)-e(j)
          f2(j) = a1*du+a5*dv
          f2v(j)= d1*duk+d7*dvk
          f3(j) = a5*du+a2*dv
          f3v(j)= d5*duk+d2*dvk
          f4(j) = (t4+t7)*du+(t5+t6)*dv+t8*dei
          f4v(j) = 0.5*d1*0.5*(ukp1**2-ukm1**2)
          f4v(j) = f4v(j)+0.5*d2*0.5*(vkp1**2-vkm1**2)
          f4v(j) = f4v(j)+d5*0.5*(vkp1+vkm1)*duk
          f4v(j) = f4v(j)+d7*0.5*(ukp1+ukm1)*dvk
          f4v(j) = f4v(j)+d4*dek
        enddo

        if (nhalo.eq.1) then
          je = jend - 1

          j=jend
          s(j,k,2) = s(j,k,2)+(f2(j)-f2(j-1))+0.5*(f2v(j-2)
     1          -4.*f2v(j-1)+3*f2v(j))
          s(j,k,3) = s(j,k,3)+(f3(j)-f3(j-1))+0.5*(f3v(j-2)
     1          -4.*f3v(j-1)+3*f3v(j))
          s(j,k,4) = s(j,k,4)+(f4(j)-f4(j-1))+0.5*(f4v(j-2)
     1          -4.*f4v(j-1)+3*f4v(j))

        elseif (nhalo.gt.1) then
          je = jend 
        endif
c
C$AD II-LOOP
        do j = jbeg,je
          s(j,k,2) = s(j,k,2)+(f2(j)-f2(j-1))+0.5*(f2v(j+1)-f2v(j-1))
          s(j,k,3) = s(j,k,3)+(f3(j)-f3(j-1))+0.5*(f3v(j+1)-f3v(j-1))
          s(j,k,4) = s(j,k,4)+(f4(j)-f4(j-1))+0.5*(f4v(j+1)-f4v(j-1))
        enddo

      enddo

      return
      end

c***********************************************************************
      subroutine vmutur(x,y,q,s,turmu,xx,xy,yx,yy,ug,vg,jd,kd)
c
c  turbulent eddy viscosity. model is algebraic model of baldwin
c  and lomax
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd)
      real x(jd,kd), y(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)

      ! local variables
      real :: fmax(mdim),smax(mdim),vmax(mdim),vmin(mdim),
     <          grdnrm(mdim),ra(mdim),snm(mdim),snp(mdim),
     <          vrtm(mdim),vrtp(mdim)

      real :: sn(jd,kd),vort(jd,kd)
      real :: qu(jd,kd),qv(jd,kd)

      integer kedge,kedge1,j1,j2,k,j
      real dre,duj,vjp1,vjm1,ejp1,dx2,dy2
      real c1,c2,c3,c4,c5,re,wmu,ux,vx,uy,vy,tx,ty
      real rhomuw,uu,dx,dy,fl,flp,flm,dfm,dfp,dsm,dsp
      real am,bm,si,fli,dvi2,t1,t2,fwake,t3,flkeb,slen,tmi
      real dv2,fkleb

c***********************************************************************
      data c1,c2,c3,c4,c5/0.4,26.0,0.01688,1.00,0.3/
c***********************************************************************
c
      re = rey
c
c..compute vorticity and normal distances from wall for entire flowfield
c
      kedge   = (3*kd)/4
      kedge1  = kedge + 1
      j1      = jtail1
      j2      = jtail2
      dx2     = 0.5
      dy2     = 0.5
c
c..zero out turmu
c
      do 12 k = 1,kd
      do 12 j = 1,jd
        turmu(j,k) = 0.0
   12 continue
c
c..set wall molecular viscosity
c
      wmu = 1.0
c
cgrs..for compatibility with turbulence model
c
        do 120 k = 1,kedge1+1
        do 120 j = j1-1,j2+1
          qu(j,k)  = q(j,k,2)  - q(j,k,1)*ug(j,k)
          qv(j,k)  = q(j,k,3)  - q(j,k,1)*vg(j,k)
 120    continue
c
        k   = 1
        do 1 j = j1,j2
c
          ux  = (qu(j+1,k)/q(j+1,k,1) - qu(j-1,k)/q(j-1,k,1)) * dx2
          vx  = (qv(j+1,k)/q(j+1,k,1) - qv(j-1,k)/q(j-1,k,1)) * dx2
          uy  = -(3.0*qu(j,k)/q(j,k,1) - 4.0*qu(j,k+1)/q(j,k+1,1) +
     &                qu(j,k+2)/q(j,k+2,1))*dy2
          vy  = -(3.0*qv(j,k)/q(j,k,1) - 4.0*qv(j,k+1)/q(j,k+1,1) +
     &                qv(j,k+2)/q(j,k+2,1))*dy2
c     
          tx  =  xy(j,k)*ux -xx(j,k)*vx +yy(j,k)*uy -yx(j,k)*vy
          ty  =  xx(j,k)*vx -xy(j,k)*ux +yx(j,k)*vy -yy(j,k)*uy
c
          vort(j,k) = sqrt(tx**2 +ty**2)
          sn(j,k)  = 0.0
    1   continue
        k   = 1
        do 2 j = j1,j2
          rhomuw  = q(j,k,1)*q(j,k,nq)/wmu
          ra(j)   = sqrt( re*rhomuw*vort(j,k) )/c2
c
cgrs..for moving blades
cgrs..   uu = sqrt(uwal**2 + vwal**2)/q(j,k,1)
c
          uu = sqrt(qu(j,k)**2 +qv(j,k)**2)/q(j,k,1)
c
          fmax(j)  = 1.e-3
          vmax(j)  = uu
          vmin(j)  = uu
          vrtp(j)  = 0.0
          vrtm(j)  = 0.0
          grdnrm(j) = sqrt(yx(j,k)**2 +yy(j,k)**2)
    2   continue
c
        do 3 k = 2,kedge1
        do 3 j = j1,j2
          ux  = (qu(j+1,k)/q(j+1,k,1) - qu(j-1,k)/q(j-1,k,1)) * dx2
          vx  = (qv(j+1,k)/q(j+1,k,1) - qv(j-1,k)/q(j-1,k,1)) * dx2
          uy  = -(3.0*qu(j,k)/q(j,k,1) - 4.0*qu(j,k+1)/q(j,k+1,1) +
     &                qu(j,k+2)/q(j,k+2,1))*dy2
          vy  = -(3.0*qv(j,k)/q(j,k,1) - 4.0*qv(j,k+1)/q(j,k+1,1) +
     &                qv(j,k+2)/q(j,k+2,1))*dy2
c     
          tx  =  xy(j,k)*ux -xx(j,k)*vx +yy(j,k)*uy -yx(j,k)*vy
          ty  =  xx(j,k)*vx -xy(j,k)*ux +yx(j,k)*vy -yy(j,k)*uy
c
          vort(j,k) = sqrt(tx**2 +ty**2)
c
          dx  = x(j,k) -x(j,1)
          dy  = y(j,k) -y(j,1)
c
          sn(j,k)  = (yx(j,1)*dx +yy(j,1)*dy)/grdnrm(j)
    3   continue
        k   = 2
        do 4 j = j1,j2
          snp(j)   = sn(j,k+1)
          smax(j)  = sn(j,k)
          snm(j)   = sn(j,k-1)
    4   continue
c
c..compute fmax, smax, vmax, vmin
c
        do 11 k = 2,kedge
        do 11 j = j1,j2
c
          fl       = sn(j,k)*vort(j,k)*(1.0 -exp(-ra(j)*sn(j,k)))
          fmax(j)  = amax1(fmax(j),fl)
ccray          smax(j)  = cvmgt(sn(j,k),smax(j),fmax(j).eq.fl)
ccray          snp(j)   = cvmgt(sn(j,k+1),snp(j),fmax(j).eq.fl)
ccray          snm(j)   = cvmgt(sn(j,k-1),snm(j),fmax(j).eq.fl)
ccray          vrtp(j)  = cvmgt(vort(j,k+1),vrtp(j),fmax(j).eq.fl)
ccray          vrtm(j)  = cvmgt(vort(j,k-1),vrtm(j),fmax(j).eq.fl)
          if( fmax(j).eq.fl ) then
            smax(j)  = sn(j,k)
            snp(j)   = sn(j,k+1)
            snm(j)   = sn(j,k-1)
            vrtp(j)  = vort(j,k+1)
            vrtm(j)  = vort(j,k-1)
          end if
c
          uu = sqrt(qu(j,k)**2 +qv(j,k)**2)/q(j,k,1)
c
cgrs..for moving blades
cgrs..    uu = sqrt(uwal**2 + vwal**2)/q(j,k,1)
c
          vmax(j)  = amax1(uu,vmax(j))
          vmin(j)  = amin1(uu,vmin(j))
   11   continue
c
c..interpolation to improve estimate of fmax and smax
c
        do 21 j = j1,j2
          flp = snp(j)*vrtp(j)*(1.0 -exp(-ra(j)*snp(j)))
          flm = snm(j)*vrtm(j)*(1.0 -exp(-ra(j)*snm(j)))
          flp = amax1(flp,1.e-3)
          flm = amax1(flm,1.e-3)
          dfm = fmax(j) -flm
          dfp = flp -fmax(j)
          dsm = smax(j) -snm(j)
          dsp = snp(j) -smax(j)
          am  = (dsp**2*dfm +dsm**2*dfp)/(dsp*dsm*(dsp+dsm))
          bm  = (dsm*dfp -dsp*dfm)/(dsp*dsm*(dsp+dsm))
ccray            si  = cvmgt(smax(j) -0.5*am/(bm+1.e-21),smax(j),bm.lt.0.0 )
ccray            fli = cvmgt(fmax(j) -0.25*am**2/(bm+1.e-21),fmax(j),bm.lt.0.0)
          if( bm.lt.0.0 ) then
            si  = smax(j) -0.5*am/(bm+1.e-21)
            fli = fmax(j) -0.25*am**2/(bm+1.e-21)
          else
            si  = smax(j)
            fli = fmax(j)
          end if
ccray          fmax(j)  = cvmgt(fmax(j),fli,fli.lt.fmax(j)
ccray     *                       .or. si.lt.snm(j) .or. si.gt.snp(j))
ccray          smax(j)  = cvmgt(smax(j),si,fli.lt.fmax(j)
ccray     *                      .or. si.lt.snm(j) .or. si.gt.snp(j))
          if( fli.lt.fmax(j).or. si.lt.snm(j) .or. si.gt.snp(j) ) then
            fmax(j)  = fmax(j)
            smax(j)  = smax(j)
          else
            fmax(j)  = fli
            smax(j)  = si
          end if
   21   continue
c
c..compute outer eddy viscosity
c
        do 31 k = 2,kedge
        do 31 j = j1,j2
          dv2   = (vmax(j) -vmin(j))**2
          t1    = smax(j)*fmax(j)
          t2    = c4*smax(j)*dv2/fmax(j)
          fwake = amin1(t1,t2)
          t3    = (c5*sn(j,k)/smax(j))
          fkleb = 1.0 +5.5*amin1(1.e5,t3)**6
          turmu(j,k) = c3*q(j,k,1)*q(j,k,nq)*fwake/fkleb
   31   continue
c
c..compute inner eddy viscosity and set final eddy viscosity
c
        do 41 k = 2,kedge
        do 41 j = j1,j2
          slen  = c1*sn(j,k)*(1.0 -exp(-ra(j)*sn(j,k)))
          tmi   = q(j,k,1)*q(j,k,nq)*slen**2*vort(j,k)
          turmu(j,k) = amin1(turmu(j,k),tmi)*re
   41   continue

      return
      end


c*************************************************************************
