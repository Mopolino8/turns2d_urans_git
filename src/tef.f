C     ******************************************************************
C     latest updated Jun 1, 2015 by Asitav
C     
C     Rotates the TEF of a 3d-wing
C     Changes the base grid (xg,yg,zg)
C     ********************************
C     Stripped off of 3d rigid_flap code: 
C      /nfs/asitav/shreyas_CFDCSD/codes/overturns_slat/src/overturns/deform.f
C     ********************************

      subroutine move_tef(init,x,y,xv,yv,xold,yold,xole,yole,ug,vg,
     .                    jd,kd,im)
C     x->xg
C     y->yg ->1 station
C     z->zg
C
C***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init,im
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xold(jmax,kmax),yold(jmax,kmax),xole(jmax,kmax),yole(jmax,kmax)

c ..   local variables

       integer j,k,ihar,jtail1v,jtail2v,jlev,jv,kv
       real cs,ss,xtmp,ytmp,xo,yo,csold,ssold,xnew,ynew

C .....Asitav
       real theta_fnew
       real::f1,f2,f3,val,yplane,angle,angle0,theta,dist,xc,yc
       real ttef11,ttef12,ttef13,ttef14,ttef31,ttef32,ttef33,ttef34
       real::dtheta_f1,dtheta_f2,y1,y2,delx1,delx2,psi_rot
       real, allocatable :: xsurf(:),ysurf(:)

       real,parameter::rmax=2.5, rmin=0.2

       real, save :: theta_prev

       allocate(xsurf(jmax),ysurf(jmax))

       pi = 4.0*atan(1.0)
       if(init.eq.1) then
         jtail1 = jtail(im)+nhalo
         jtail2 = jd-jtail1+1
       end if
       jtail1v = jtail1-nhalo
       jtail2v = jtail2-nhalo
       !print*,'info: ',jtail1v,jtail2v
       jlev    = jle-nhalo

       theta_fnew = theta_f0
       psi_rot    = rf_tef*totime

       do ihar=1,nharmflap
          theta_fnew=theta_fnew+ampFlap(ihar)
     &         *sin(ihar*psi_rot-phiFlap(ihar))
!amm     &         *cos(ihar*psi_rot-phiFlap(ihar))
       enddo

       if (init.eq.1) then
          dtheta_f2=theta_finit !theta_fnew
          dtheta_f1=0.0D0
       else
          dtheta_f2=theta_fnew
          dtheta_f1=-theta_prev
       endif
       !write(6,'(A,2(F12.5,x))') "psi_rot,TEF angle ", psi_rot*180./pi,theta_fnew*180./pi

       do j=1,jmax
        xsurf(j) = xv(j,1) 
        ysurf(j) = yv(j,1)
       end do

       xc = (pxc*(xv(jtail1v,1) - xv(jlev,1)) + xv(jlev,1))

C       locate hinge
C       ------------
       hinge: do jv=jtail1v,jtail2v
        delx1 = xv(jv  ,1)-xc
        delx2 = xv(jv+1,1)-xc
        if( (delx1*delx2 < 0.0) .and. (jv < jlev) ) then
         y1 = 0.5*(yv(jv,1)+yv(jv+1,1))
        else if( (delx1*delx2 < 0.0) .and. (jv > jlev ) ) then
         y2 = 0.5*(y(jv,1)+y(jv+1,1))
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

        !vectex flap
        !---------------------------------------------------------------
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
          else if (j < jtail1) then
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

          cs=cos( (dtheta_f1+dtheta_f2) *f3*f2*f1)
          ss=sin( (dtheta_f1+dtheta_f2) *f3*f2*f1)

C  .......The t-matrix for te flap
          ttef11 = cs
          ttef12 = 0.0
          ttef13 = ss
          ttef14 = xc*(1. - cs) - yc*ss

          ttef31 = -ss
          ttef32 = 0.0
          ttef33 = cs
          ttef34 = xc*ss + yc*(1-cs)
C  .......move grid
          !amm xo=xv(j,k) 
          !amm yo=yv(j,k) 
          !amm xv(j,k) = xo*ttef11+yo*ttef13 + ttef14
          !amm yv(j,k) = xo*ttef31+yo*ttef33 + ttef34

          !..store the old mesh for time metrics
          xole(j,k) = xold(j,k)
          yole(j,k) = yold(j,k)

          xold(j,k) = xv(j,k)
          yold(j,k) = yv(j,k)

          !..update mesh
          xv(j,k) = xold(j,k)*ttef11+yold(j,k)*ttef13 + ttef14
          yv(j,k) = xold(j,k)*ttef31+yold(j,k)*ttef33 + ttef34
         end do
        end do
        !---------------- end vertex flap ------------------------------

        !cell center flap
        !---------------------------------------------------------------
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

          cs=cos( (dtheta_f1+dtheta_f2) *f3*f2*f1)
          ss=sin( (dtheta_f1+dtheta_f2) *f3*f2*f1)

          csold=cos(-1.0*(dtheta_f1+dtheta_f2) *f3*f2*f1)
          ssold=sin(-1.0*(dtheta_f1+dtheta_f2) *f3*f2*f1)

C  .......The t-matrix for te flap
          ttef11 = cs
          ttef12 = 0.0
          ttef13 = ss
          ttef14 = xc*(1. - cs) - yc*ss

          ttef31 = -ss
          ttef32 = 0.0
          ttef33 = cs
          ttef34 = xc*ss + yc*(1-cs)
C  .......move grid
          xnew = x(j,k)*ttef11+y(j,k)*ttef13 + ttef14
          ynew = x(j,k)*ttef31+y(j,k)*ttef33 + ttef34

          !save grid velocity
          if(ntac.eq.1) then
            ug(j,k)=(xnew-x(j,k))/dt
            vg(j,k)=(ynew-y(j,k))/dt
          else
            !The t-matrix for te flap old time step
            ttef11 = csold
            ttef12 = 0.0
            ttef13 = ssold
            ttef14 = xc*(1. - csold) - yc*ssold

            ttef31 = -ssold
            ttef32 = 0.0
            ttef33 = csold
            ttef34 = xc*ssold + yc*(1-csold)

            xo = x(j,k)*ttef11+y(j,k)*ttef13 + ttef14
            yo = x(j,k)*ttef31+y(j,k)*ttef33 + ttef34

            ug(j,k)=(1.5*xnew-2*x(j,k)+0.5*xo)/dt
            vg(j,k)=(1.5*ynew-2*y(j,k)+0.5*yo)/dt
          end if
          
          x(j,k) = xnew
          y(j,k) = ynew

         end do
        end do
        !---------------------------------------------------------------

       theta_prev=theta_fnew
       !do jv=jtail1v,jtail2v
       ! write(1000,*)xv(jv,1),yv(jv,1),jv
       !end do

       deallocate(xsurf,ysurf)
      end subroutine move_tef
