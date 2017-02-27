!TS TEF Time metrics computation
!
! author: asitav
! date  : Jun 2, 2015
!*************************************************************************
      subroutine gettimemetrics_tef_TS(x,y,xv,yv,ug,vg,ugv,vgv,
     .                                 jd,kd,time)
C
C***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real time
      real x(jd,kd),y(jd,kd),ug(jd,kd),vg(jd,kd)
      real xv(jmax,kmax),yv(jmax,kmax),ugv(jmax,kmax),vgv(jmax,kmax)


c ..   local variables

       integer j,k,ihar,jtail1v,jtail2v,jlev,jv,kv,fint
       real cs,ss,xtmp,ytmp,xo,yo

C .....Asitav
       real::thetaf,psi_rot
       real::ttef11,ttef12,ttef13,ttef14,ttef31,ttef32,ttef33,ttef34
       real::f1,f2,f3,val,angle,angle0,dist,xc,yc
       real::y1,y2,delx1,delx2

       real:: xsurf(jmax),ysurf(jmax)

       !derivatives
       real:: dthetaf,dpsi_rot
       real:: dttef11,dttef12,dttef13,dttef14,dttef31,dttef32,dttef33,dttef34
       real:: dcs,dss

       real,parameter::rmax=2.5, rmin=0.2

       !some constants
       !----------------------------
       pi = 4.0*atan(1.0)
       jtail1v = jtail(1)
       jtail2v = jmax-jtail(1)+1
       jlev    = (jtail1v+jtail2v)/2

       !update for cell center co-ord
       jtail1 = jtail1 + nhalo
       jtail2 = jd-jtail1+1
       jle    = (jtail1+jtail2)/2
       !----------------------------

       !compute state
       !-------------------------------------
       psi_rot = rotf*time

       thetaf  = theta_f0
       do ihar=1,nharmflap
          thetaf=thetaf+ampFlap(ihar)
     &         *sin(ihar*psi_rot-phiFlap(ihar))
       enddo

!       write(6,'(A,2(F12.5,x))') "psi_rot,TEF angle ", 
!     .        psi_rot*180./pi,thetaf*180./pi
       !-------------------------------------

       !time derivative of thetaf
       !--------------------------
       dpsi_rot = rotf

       dthetaf  = 0.0
       do ihar=1,nharmflap
          dthetaf=dthetaf+ihar*dpsi_rot*ampFlap(ihar)
     &         *cos(ihar*psi_rot-phiFlap(ihar))
       enddo
       !--------------------------
       
       do j=1,jmax
        xsurf(j) = xv(j,1) 
        ysurf(j) = yv(j,1)
       end do

       xc = (pxc*(xv(jtail1v,1) - xv(jlev,1)) + xv(jlev,1))

C      !locate hinge
C      !------------
       hinge: do jv=jtail1v,jtail2v
        delx1 = xv(jv  ,1)-xc
        delx2 = xv(jv+1,1)-xc
        if( (delx1*delx2 < 0.0) .and. (jv < jlev) ) then
         y1 = 0.5*(yv(jv,1)+yv(jv+1,1))
        else if( (delx1*delx2 < 0.0) .and. (jv > jlev ) ) then
         y2 = 0.5*(yv(jv,1)+yv(jv+1,1))
         exit hinge
        else
         cycle hinge
        end if
       end do hinge
c
c check jaina
c
        y1=(y1+y2)*0.5
        y2=y1

        !------------
        f3 = 1.0 !along z-dir
        !------------

        !----------------------------------------------------------------------
        !now find time metics for vertex and cell centers
        
        !------------------------
        !vertex time metrics
        !------------------------
        do k = 1,kmax
         do j = 1,jmax
C........ Point of rotation
          if(j<jlev) then
           yc = y1
          else if(j>jlev) then
           yc = y2
          end if

          dist=sqrt((xv(j,k)-xsurf(j))**2+(yv(j,k)-ysurf(j))**2);
          if (j > jtail2v) then
           dist=sqrt((xv(j,k)-xsurf(jtail2v))**2+(yv(j,k)-ysurf(jtail2v))**2);
          else if (j < jtail1v) then
           dist=sqrt((xv(j,k)-xsurf(jtail1v))**2+(yv(j,k)-ysurf(jtail1v))**2);
          end if

          if (dist > rmax) then
           f1=0.0;
          else if (dist < rmin) then
           f1=1.0;
          else
           val=dist-rmin;
           val=val/(rmax-rmin)*pi;
           f1=(1+cos(val))*0.5;
          end if

          angle =atan2( (yv(j,k)-yc), (xv(j,k)-xc) );
          angle0=atan2( (yv(jtail1v,1)-yc), (xv(jtail1v,1)-xc) );
          angle = angle - angle0;

          if (abs(angle) > pi*0.5+dela) then
           f2=0.0
          elseif (abs(angle)<pi*0.5-dela)  then
           f2=1.0
          else
           val=abs(angle)-pi*0.5+dela
           val=pi*val/(2.0*dela)
           f2=(1+cos(val))*0.5
          end if

          cs=cos( (thetaf) *f3*f2*f1)
          ss=sin( (thetaf) *f3*f2*f1)

          !time derivatives
          !-----------------------------
          dcs=-(dthetaf*f3*f2*f1)*sin( (thetaf) *f3*f2*f1)
          dss= (dthetaf*f3*f2*f1)*cos( (thetaf) *f3*f2*f1)
          !-----------------------------

C  .......The t-matrix for te flap
          ttef11 = cs
          ttef12 = 0.0
          ttef13 = ss
          ttef14 = xc*(1. - cs) - yc*ss

          ttef31 = -ss
          ttef32 = 0.0
          ttef33 = cs
          ttef34 = xc*ss + yc*(1-cs)

          !time derivatives
          !---------------------------------
          dttef11 = dcs
          dttef12 = 0.0
          dttef13 = dss
          dttef14 = xc*(-dcs) - yc*dss

          dttef31 = -dss
          dttef32 = 0.0
          dttef33 = dcs
          dttef34 = xc*dss + yc*(-dcs)
          !---------------------------------
C  .......move grid
          xo=xv(j,k) 
          yo=yv(j,k) 
          ugv(j,k) = xo*dttef11+yo*dttef13 + dttef14
          vgv(j,k) = xo*dttef31+yo*dttef33 + dttef34

         end do
        end do

        !------------------------
        !cell center time metrics
        !------------------------
        do k = 1,kd
         do j = 1,jd
C........ Point of rotation
          if(j<jle) then
           yc = y1
          else if(j>jle) then
           yc = y2
          end if

          dist=sqrt((x(j,k)-x(j,1))**2+(y(j,k)-y(j,1))**2);
          if (j > jtail2) then
           dist=sqrt((x(j,k)-x(jtail2,1))**2+(y(j,k)-y(jtail2,1))**2);
          else if (j < jtail1) then
           dist=sqrt((x(j,k)-x(jtail1,1))**2+(y(j,k)-y(jtail1,1))**2);
          end if

          if (dist > rmax) then
           f1=0.0;
          else if (dist < rmin) then
           f1=1.0;
          else
           val=dist-rmin;
           val=val/(rmax-rmin)*pi;
           f1=(1+cos(val))*0.5;
          end if

          angle =atan2( (y(j,k)-yc), (x(j,k)-xc) );
          if(j<jle) then
            angle0=atan2( (y(jtail1,1)-yc), (x(jtail1,1)-xc) );
          else
            angle0=atan2( (y(jtail2,1)-yc), (x(jtail2,1)-xc) );
          end if
          angle = angle - angle0;

          if (abs(angle) > pi*0.5+dela) then
           f2=0.0
          elseif (abs(angle)<pi*0.5-dela)  then
           f2=1.0
          else
           val=abs(angle)-pi*0.5+dela
           val=pi*val/(2.0*dela)
           f2=(1+cos(val))*0.5
          end if

          cs=cos( (thetaf) *f3*f2*f1)
          ss=sin( (thetaf) *f3*f2*f1)

          !time derivatives
          !-----------------------------
          dcs=-(dthetaf*f3*f2*f1)*sin( (thetaf) *f3*f2*f1)
          dss= (dthetaf*f3*f2*f1)*cos( (thetaf) *f3*f2*f1)
          !-----------------------------

C  .......The t-matrix for te flap
          ttef11 = cs
          ttef12 = 0.0
          ttef13 = ss
          ttef14 = xc*(1. - cs) - yc*ss

          ttef31 = -ss
          ttef32 = 0.0
          ttef33 = cs
          ttef34 = xc*ss + yc*(1-cs)

          !time derivatives
          !---------------------------------
          dttef11 = dcs
          dttef12 = 0.0
          dttef13 = dss
          dttef14 = xc*(-dcs) - yc*dss

          dttef31 = -dss
          dttef32 = 0.0
          dttef33 = dcs
          dttef34 = xc*dss + yc*(-dcs)
          !---------------------------------
C  .......move grid
          xo=x(j,k) 
          yo=y(j,k) 
          ug(j,k) = xo*dttef11+yo*dttef13 + dttef14
          vg(j,k) = xo*dttef31+yo*dttef33 + dttef34

         end do
        end do

      end subroutine gettimemetrics_tef_TS
!*************************************************************************
