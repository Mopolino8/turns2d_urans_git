C     ******************************************************************
C     latest updated Jun 1, 2015 by Asitav
C     
C     Rotates the TEF of a 3d-wing
C     Changes the base grid (xg,yg,zg)
C     ********************************
C     Stripped off of 3d rigid_flap code: 
C      /nfs/asitav/shreyas_CFDCSD/codes/overturns_slat/src/overturns/deform.f
C     ********************************

      subroutine move_tef_ts(xv,yv,xv1,yv1,time)
C     x->xg
C     y->yg ->1 station
C     z->zg
C
C***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      real time
      real xv(jmax,kmax),yv(jmax,kmax),xv1(jmax,kmax),yv1(jmax,kmax)

c ..   local variables

       integer j,k,ihar,jtail1v,jtail2v,jlev,jv,kv,fint
       real cs,ss,xtmp,ytmp,xo,yo

C .....Asitav
       real theta_fnew,theta_fold
       real::f1,f2,f3,val,yplane,angle,angle0,theta,dist,xc,yc
       real ttef11,ttef12,ttef13,ttef14,ttef31,ttef32,ttef33,ttef34
       real::dtheta_f1,dtheta_f2,y1,y2,delx1,delx2,psi_rot
       real:: xsurf(jmax),ysurf(jmax)

       real,parameter::rmax=2.5, rmin=0.2

       !amm real, save :: theta_prev !not required in TS

       pi = 4.0*atan(1.0)
       jtail1 = jtail(1)
       jtail2 = jmax-jtail1+1
       jtail1v = jtail1!-nhalo
       jtail2v = jtail2!-nhalo
       jlev    = jle!-nhalo

       theta_fnew = theta_f0
       !amm psi_rot    = rf_tef*time
       psi_rot    = rotf*time

       do ihar=1,nharmflap
          theta_fnew=theta_fnew+ampFlap(ihar)
     &         *sin(ihar*psi_rot-phiFlap(ihar))
!amm     &         *cos(ihar*psi_rot-phiFlap(ihar))
       enddo

       !amm if (init) then
       !amm    dtheta_f2=theta_finit !theta_fnew
       !amm    dtheta_f1=0.0D0
       !amm else
       !amm    dtheta_f2=theta_fnew
       !amm    dtheta_f1=-theta_prev
       !amm endif
       dtheta_f2=theta_fnew !theta_fnew for TS
       dtheta_f1=0.0D0 !assume started from zero flap 
       write(6,'(A,2(F12.5,x))') "psi_rot,TEF angle ", 
     .        psi_rot*180./pi,theta_fnew*180./pi

       do j=1,jmax
        xsurf(j) = xv1(j,1) 
        ysurf(j) = yv1(j,1)
       end do

       xc = (pxc*(xv1(jtail1v,1) - xv1(jlev,1)) + xv1(jlev,1))

C       locate hinge
C       ------------
       hinge: do jv=jtail1v,jtail2v
        delx1 = xv1(jv  ,1)-xc
        delx2 = xv1(jv+1,1)-xc
        if( (delx1*delx2 < 0.0) .and. (jv < jlev) ) then
         y1 = 0.5*(yv1(jv,1)+yv1(jv+1,1))
        else if( (delx1*delx2 < 0.0) .and. (jv > jlev ) ) then
         y2 = 0.5*(yv1(jv,1)+yv1(jv+1,1))
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

        do k = 1,kmax
         do j = 1,jmax
C........ Point of rotation
          if(j<jlev) then
           yc = y1
          else if(j>jlev) then
           yc = y2
          end if

          dist=sqrt((xv1(j,k)-xsurf(j))**2+(yv1(j,k)-ysurf(j))**2);
          if (j > jtail2) then
           dist=sqrt((xv1(j,k)-xsurf(jtail2))**2+(yv1(j,k)-ysurf(jtail2))**2);
          else if (j < jtail1) then
           dist=sqrt((xv1(j,k)-xsurf(jtail1))**2+(yv1(j,k)-ysurf(jtail1))**2);
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

          angle =atan2( (yv1(j,k)-yc), (xv1(j,k)-xc) );
          angle0=atan2( (yv1(jtail1v,1)-yc), (xv1(jtail1v,1)-xc) );
          !amm if(j<jle) then
          !amm  angle0=atan2( (yv1(jtail1v,1)-yc), (xv1(jtail1v,1)-xc) );
          !amm else
          !amm  angle0=atan2( (yv1(jtail2v,1)-yc), (xv1(jtail2v,1)-xc) );
          !amm end if
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
          xo=xv1(j,k) 
          yo=yv1(j,k) 
          xv(j,k) = xo*ttef11+yo*ttef13 + ttef14
          yv(j,k) = xo*ttef31+yo*ttef33 + ttef34

         end do
        end do

       !!amm theta_prev=theta_fnew
       !fint = int(time*rotf*nspec/2.0/pi)
       !do jv=jtail1v,jtail2v
       !  write(1000+fint,*)xv(jv,1),yv(jv,1)
       !end do

      end subroutine move_tef_ts
!*************************************************************************

