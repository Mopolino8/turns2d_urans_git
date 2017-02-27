c***********************************************************************
      subroutine move(init,x,y,xx,xy,ug,yx,yy,vg,jd,kd)
c
c  this stores the old surface metrics as well as rotates the grid and 
c  updates the metrics.  
c
c  if iunst=1
c      rf = reduced frequency
c      maxang = amplitude of pitching oscillation
c         pos/neg about origin
c  if iunst=2
c      rf = time for pitching ramp
c      maxang = change in pitching angle
c  if iunst=3
c      read in angle for each time step
c  if iunst=4
c      rf = reduced frequency
c      maxang = amplitude of pitching oscillation
c         pos only about orgin
c  if iunst=5
c      interpolation between grids, unsteady but not simply transformation
c      allowed here
c      all metrics get recalculated after adapt subroutine call
c      just keep this so iunst is not zero, steady state
c  if iunst=6
c      step change in grid velocity (modify, grid vel metric, keep
c                                    grid postion metric const)
c      phasestart=total time for pitching motion	
c

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd)

      ! local variables
      real,allocatable :: xtmpm(:,:),ytmpm(:,:)
      real thnm1,vfact,angnew,angold,dummy,t
      real scalnew,scalold,sn,cn,sno,cno,xold,yold,xnm1,ynm1
      integer j,k,init
      
      allocate(xtmpm(jd,kd),ytmpm(jd,kd))

      thnm1 = 0.0

c for grid interpolation
      if(iunst.eq.5) then
        thetao = 0.0
        thetan = 0.0
      endif

c for steady problems
      if(iunst.eq.0) then
        thetao = 0.0
        thetan = 0.0
      endif

c for ramps
      if(iunst.eq.3) then
        if(init.eq.1) then
          thnm1 = 0.0
          thetao = 0.0
          read(17,*) thetan     
          vfact  = 0.0
        elseif(init.eq.2) then
          thnm1 = thetao
          thetao = thetan 
          read(17,*) thetan
          vfact = 1.0
        else
          thnm1 = thetao
          thetao = thetan
          read(17,*) thetan
          vfact = 1.0
        endif
      endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c this is oscillations about origin both positive and negative
      if (iunst.eq.1) then
        angnew = totime*rf
        if(init.eq.1) then
          angold = 0.
          vfact=0.
        else
          angold = (totime-dt)*rf
          vfact=1.
        endif
        thetan = sin(angnew)*(pi/180.)*angmax
        thetao = sin(angold)*(pi/180.)*angmax
      endif
c
c this is oscillations purely positive from origin (le of airfoil, not 1/4c)
      if (iunst.eq.4) then
        angnew = totime*rf
        if(init.eq.1) then
          angold = 0.
          vfact=0.
        else
          angold = (totime-dt)*rf
          vfact=1.
        endif
        thetan = sin(angnew-pi/2.)*(pi/180.)*angmax/2.
     c                                      +angmax/2.*(pi/180.)
        thetao = sin(angold-pi/2.)*(pi/180.)*angmax/2.
     c                                      +angmax/2.*(pi/180.)
      endif
      if (iunst.eq.6) then
        vfact=1.
        thetao = 0.0
        thetan = 0.0
	if (istep.eq.0) then
          read(200) dummy,dummy
          read(200)((xtmpm(j,k),j=1,jd),k=1,kd),
     <             ((ytmpm(j,k),j=1,jd),k=1,kd)
          close(200)
	endif
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (iunst.eq.2) then
        t = totime/rf
        if(t.gt.1.) t=1.
        if(t.lt.0.) t=0.
        scalnew = (10.-15.*t+6.*t*t)*t**3
        if(init.eq.1) then
          scalold = 0.
          vfact=0.
        else
          t = (totime-dt)/rf
          if(t.gt.1.) t=1.
          if(t.lt.0.) t=0.
          scalold = (10.-15.*t+6.*t*t)*t**3
          vfact=1.
        endif
        thetan = scalnew*(pi/180.)*angmax
        thetao = scalold*(pi/180.)*angmax
      endif
******************************************************************
******************************************************************
******************************************************************
c now calculating metrics of transformation for grid
c
      dang = thetan*180./pi
c
      sn = sin(thetan-thetao)
      cn = cos(thetan-thetao)
      sno = sin(thnm1-thetao)
      cno = cos(thnm1-thetao)
c
      do 81 k = 1,kd
      do 81 j = 1,jd
c..grid
        xold  = x(j,k)
        yold  = y(j,k)
        if (iunst.ne.4) then    !for rotation about origin
          x(j,k)  = +sn*yold+cn*xold
          y(j,k)  = -sn*xold+cn*yold
        else                    !for rotation about xac
          x(j,k)  = (+sn*(yold-yac)+cn*(xold-xac))+xac
          y(j,k)  = (-sn*(xold-xac)+cn*(yold-yac))+yac
        endif
        xnm1    = +sno*yold+cno*xold
        ynm1    = -sno*xold+cno*yold
c..grid velocities
        ug(j,k) = (x(j,k)-xold)/dt*vfact
        vg(j,k) = (y(j,k)-yold)/dt*vfact
        if(iunst.eq.3 .and. ntac.ge.2 .and. istep.gt.1) then
          ug(j,k) = (1.5*x(j,k)-2.*xold+0.5*xnm1)/dt*vfact
          vg(j,k) = (1.5*y(j,k)-2.*yold+0.5*ynm1)/dt*vfact
        endif
	if (iunst.eq.6) then		!for instataneous change in pitch rate
c	  ug(j,k) = (xtmpm(j,k)-xold)/phasestart*vfact
c         vg(j,k) = (ytmpm(j,k)-yold)/phasestart*vfact
c         ug(j,k) = (xold-xtmpm(j,k))/phasestart*vfact
c         vg(j,k) = (yold-ytmpm(j,k))/phasestart*vfact
c         ug(j,k) = (xtmpm(j,k)-xold)/dt*vfact
c         vg(j,k) = (ytmpm(j,k)-yold)/dt*vfact
	if ((k.lt.60).and.((j.gt.58).and.(j.lt.160))) then
          ug(j,k) = -(xtmpm(j,k)-xold)/.125*vfact*fsmach
          vg(j,k) = -(ytmpm(j,k)-yold)/.125*vfact*fsmach
	else
	  ug(j,k) = 0.
	  vg(j,k) = 0.
	endif
	

	endif
c..xi metrics       
        xold  = xx(j,k)
        yold  = xy(j,k)
        xx(j,k)  = +sn*yold+cn*xold
        xy(j,k)  = -sn*xold+cn*yold
c..eta metrics   
        xold  = yx(j,k)
        yold  = yy(j,k)
        yx(j,k)  = +sn*yold+cn*xold
        yy(j,k)  = -sn*xold+cn*yold
 81   continue
c
      return
      end


c***********************************************************************
      subroutine move_new(init,x,y,xv,yv,xold,yold,xole,yole,ug,vg,jd,kd)
c
c  this stores the old surface metrics as well as rotates the grid and 
c  updates the metrics.  
c
c  if iunst=1
c      interpolation between grids, unsteady but not simply transformation
c      allowed here
c      all metrics get recalculated after adapt subroutine call
c      just keep this so iunst is not zero, steady state

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xold(jmax,kmax),yold(jmax,kmax),xole(jmax,kmax),yole(jmax,kmax)
      ! local variables

      real dtheta,sn,cn,snold,cnold,xnew,ynew,ang,xo,yo
      integer j,k
c***********************************

!     jaina, changing this subroutine, it looked too rough
!     can do better for ug(j,k),vg(j,k) with analytical expressions

      pi = 4.0*atan(1.0)
c     
      dtheta=angmax*(sin(rotf*istep*dt)-sin(rotf*(istep-1)*dt))
      sn = sin(dtheta*pi/180.)
      cn = cos(dtheta*pi/180.)
      snold = sin(-1.*dtheta*pi/180.)
      cnold = cos(-1.*dtheta*pi/180.)
c
      do 81 k = 1,kd
      do 81 j = 1,jd
c..grid

        xnew=(x(j,k)-xac)*cn+(y(j,k)-yac)*sn+xac
        ynew=-(x(j,k)-xac)*sn+(y(j,k)-yac)*cn+yac

	if(ntac.eq.1) then
 	ug(j,k)=(xnew-x(j,k))/dt
   	vg(j,k)=(ynew-y(j,k))/dt
	else
        xo=(x(j,k)-xac)*cnold+(y(j,k)-yac)*snold+xac
        yo=-(x(j,k)-xac)*snold+(y(j,k)-yac)*cnold+yac
 	ug(j,k)=(1.5*xnew-2*x(j,k)+0.5*xo)/dt
 	vg(j,k)=(1.5*ynew-2*y(j,k)+0.5*yo)/dt
	endif

        x(j,k)=xnew
        y(j,k)=ynew

  81	continue

      do k = 1,kmax
      do j = 1,jmax

c..store the old mesh for time metrics

        xole(j,k) = xold(j,k)
        yole(j,k) = yold(j,k)

        xold(j,k) = xv(j,k)
        yold(j,k) = yv(j,k)

c..update mesh

        xv(j,k) = (xold(j,k)-xac)*cn+(yold(j,k)-yac)*sn+xac
        yv(j,k) = -(xold(j,k)-xac)*sn+(yold(j,k)-yac)*cn+yac

      enddo
      enddo

	ang=atan2(y(jle,1)-y(jtail1,1),x(jtail1,1)-x(jle,1))

c	print *,'ang=',ang*180/pi



      return
      end

c***********************************************************************
      subroutine move_cyclo(init,x,y,xv,yv,xold,yold,xole,yole,ug,vg,
     &                      jd,kd,im)
c
c     This is the new grid motion subroutine which calculates the blade
c     kinematics of cycloidal rotor.

c     There are two options for grid motion:
c       fourbar_toggle = 0: This is a simple sinusoidal pitching motion
c       fourbar_toggle = 1: This is the four-bar pitching motion for 
c                           Moble's cyclocopter

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init,im
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xold(jmax,kmax),yold(jmax,kmax),xole(jmax,kmax),yole(jmax,kmax)
      ! local variables

      real xp,yp,dtheta,sn,cn,snold,cnold,xnew,ynew,ang
      real l1,l2,l3,l4,psi_rot,psi_rotold
      real ph,ph_old,theta,theta0,phi,dpsi
      real sih,sih_old,a11,a11_old,b11,b11_old,c11,c11_old,cp,cp_old
      real thetav,thetav_old
      integer j,k
c***********************************************************************

      pi = 4.0*atan(1.0)

! psi is rotational angle as measured from -z axis
      if (nblade.ne.0) then
        psi_rot = rotf*totime + (im-1)*2*pi/nblade
        psi_rotold = rotf*(totime-dt) + (im-1)*2*pi/nblade
      else
        psi_rot = rotf*totime
        psi_rotold = rotf*(totime-dt)
      endif
      dpsi = psi_rot - psi_rotold

      if (init.eq.1) dpsi = 0
c.. rotate background grid in global sense about origin

      do 181 k = 1,kd
      do 181 j = 1,jd
        xnew = cos(dpsi)*(x(j,k)-xac) - sin(dpsi)*(y(j,k)-yac) + xac
        ynew = sin(dpsi)*(x(j,k)-xac) + cos(dpsi)*(y(j,k)-yac) + yac

        if (iunst.eq.2) then
          ug(j,k) = -rotf*(ynew-yac)
          vg(j,k) =  rotf*(xnew-xac)
        else
          ug(j,k) = (xnew-x(j,k))/dt
          vg(j,k) = (ynew-y(j,k))/dt
        endif

        x(j,k) = xnew
        y(j,k) = ynew

  181 continue

      do k = 1,kmax
      do j = 1,jmax

c..store the old mesh for time metrics

        xole(j,k) = xold(j,k)
        yole(j,k) = yold(j,k)

        xold(j,k) = xv(j,k)
        yold(j,k) = yv(j,k)

c..update mesh

        xv(j,k) = cos(dpsi)*(xold(j,k)-xac) - sin(dpsi)*(yold(j,k)-yac) + xac
        yv(j,k) = sin(dpsi)*(xold(j,k)-xac) + cos(dpsi)*(yold(j,k)-yac) + yac

      enddo
      enddo

      if (.not.bodyflag(im)) go to 82

      theta0 = ang0*pi/180.
      theta = angmax*pi/180.
      phi = phase*pi/180.

      if (fourbar_toggle.EQ.1) then

c..fourbar motion specific to Moble's cyclocopter

        l1 = 6.75

        if (abs(theta-0.*pi/180).lt.1e-6) then
          l2 = 0.0
        elseif (abs(theta-2.*pi/180).lt.1e-6) then
          l2 = 0.0244
        elseif (abs(theta-4.*pi/180).lt.1e-6) then
          l2 = 0.0487
        elseif (abs(theta-6.*pi/180).lt.1e-6) then
          l2 = 0.0729
        elseif (abs(theta-10.*pi/180).lt.1e-6) then
          l2 = 0.1212
        elseif (abs(theta-15.*pi/180).lt.1e-6) then
          l2 = 0.1806
        elseif (abs(theta-20.*pi/180).lt.1e-6) then
          l2 = 0.2387
        elseif (abs(theta-25.*pi/180).lt.1e-6) then
          l2 = 0.295
        elseif (abs(theta-30.*pi/180).lt.1e-6) then
          l2 = 0.349
        elseif (abs(theta-35.*pi/180).lt.1e-6) then
          l2 = 0.4001
        elseif (abs(theta-40.*pi/180).lt.1e-6) then
          l2 = 0.55
        else
          print*,'Cannot determine l2 value for cycloidal motion'
        endif

        l3 = 6.769
        l4 = 1.26

        sih = psi_rot + pi/2 + phi
        sih_old = psi_rotold + pi/2 + phi
        a11 = sin(sih)
        a11_old = sin(sih_old)
        b11 = cos(sih)+l1/l2
        b11_old = cos(sih_old)+l1/l2
        c11 = l1/l4*cos(sih)+(l1*l1+l2*l2+l4*l4-l3*l3)/(2.*l2*l4)
        c11_old = l1/l4*cos(sih_old)+(l1*l1+l2*l2+l4*l4-l3*l3)/(2.*l2*l4)
        cp = -(l1/l3*cos(sih)+(l1*l1+l2*l2+l3*l3-l4*l4)/(2.*l2*l3))
        cp_old = -(l1/l3*cos(sih_old)+(l1*l1+l2*l2+l3*l3-l4*l4)/(2.*l2*l3))

        thetav = 2.*atan((a11-sqrt(a11*a11+b11*b11-c11*c11))/(b11+c11))
        thetav_old = 2.*atan((a11_old-sqrt(a11_old*a11_old+b11_old*
     &                b11_old-c11_old*c11_old))/(b11_old+c11_old))

c.. calculate pitch angle at current and previous timestep
        ph = pi/2. + thetav
        ph_old = pi/2. + thetav_old

      else

c.. simple sinusoidal pitching code
      ph = theta0 + theta*cos(psi_rot + pi/2 + phi)
      ph_old = theta0 + theta*cos(psi_rotold + pi/2 + phi)

      endif

      if (init.eq.1) ph_old = 0.

      dtheta = ph - ph_old

      xp = rartio*cos(psi_rot-pi/2.) + xac
      yp = rartio*sin(psi_rot-pi/2.) + yac

      sn = sin(dtheta)
      cn = cos(dtheta)
      snold = sin(-dtheta)
      cnold = cos(-dtheta)
c
      do 81 k = 1,kd
      do 81 j = 1,jd
c..grid

        xnew=(x(j,k)-xp)*cn+(y(j,k)-yp)*sn+xp
        ynew=-(x(j,k)-xp)*sn+(y(j,k)-yp)*cn+yp

        ug(j,k) = ug(j,k) + (xnew-x(j,k))/dt
        vg(j,k) = vg(j,k) + (ynew-y(j,k))/dt

        x(j,k)=xnew
        y(j,k)=ynew

  81  continue

      do k = 1,kmax
      do j = 1,jmax
        xnew=(xv(j,k)-xp)*cn+(yv(j,k)-yp)*sn+xp
        ynew=-(xv(j,k)-xp)*sn+(yv(j,k)-yp)*cn+yp

        xv(j,k)=xnew
        yv(j,k)=ynew
      enddo
      enddo


  82  continue

      if (init.eq.1) then
        do j = 1,jd
          do k = 1,kd
            ug(j,k) = 0.
            vg(j,k) = 0.
          enddo
        enddo
      endif

      ang=atan2(y(jle,1)-y(jtail1,1),x(jtail1,1)-x(jle,1))

      return
      end

c***********************************************************************
      subroutine rotateq(q,qtn,qtnm1,jd,kd)
c
c***********************************************************************
      use params_global
      implicit none

      integer jd,kd
      real q(jd,kd,nq),qtn(jd,kd,nv),qtnm1(jd,kd,nv)

      integer j,k
      real qtmp2,qtmp3,srot,cs,ss

      srot = rotf*dt
      cs = cos(srot)
      ss = sin(srot)

      do j=1,jd
      do k=1,kd
        qtmp2 = q(j,k,2)*cs - q(j,k,3)*ss
        qtmp3 = q(j,k,3)*cs + q(j,k,2)*ss
        q(j,k,2) = qtmp2
        q(j,k,3) = qtmp3
        qtmp2 = qtn(j,k,2)*cs - qtn(j,k,3)*ss
        qtmp3 = qtn(j,k,3)*cs + qtn(j,k,2)*ss
        qtn(j,k,2) = qtmp2
        qtn(j,k,3) = qtmp3
        qtmp2 = qtnm1(j,k,2)*cs - qtnm1(j,k,3)*ss
        qtmp3 = qtnm1(j,k,3)*cs + qtnm1(j,k,2)*ss
        qtnm1(j,k,2) = qtmp2
        qtnm1(j,k,3) = qtmp3
      enddo
      enddo

      end subroutine rotateq

c***********************************************************************

c***********************************************************************
      subroutine pitch(init,x,y,xv,yv,xold,yold,xole,yole,ug,vg,jd,kd)
c
c  this stores the old surface metrics as well as rotates the grid and 
c  updates the metrics.  
c
c  if iunst=1
c      interpolation between grids, unsteady but not simply transformation
c      allowed here
c      all metrics get recalculated after adapt subroutine call
c      just keep this so iunst is not zero, steady state

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xold(jmax,kmax),yold(jmax,kmax),xole(jmax,kmax),yole(jmax,kmax)
      ! local variables

      real dtheta,sn,cn,snold,cnold,xnew,ynew,ang,xo,yo
      integer j,k
c***********************************

!     jaina, changing this subroutine, it looked too rough
!     can do better for ug(j,k),vg(j,k) with analytical expressions

      pi = 4.0*atan(1.0)
c     
      !amm dtheta=angmax*(sin(rotf*istep*dt)-sin(rotf*(istep-1)*dt))
      dtheta=-angmax*(cos(rotf*istep*dt)-cos(rotf*(istep-1)*dt))
      sn = sin(dtheta*pi/180.)
      cn = cos(dtheta*pi/180.)
      snold = sin(-1.*dtheta*pi/180.)
      cnold = cos(-1.*dtheta*pi/180.)
c
      do 81 k = 1,kd
      do 81 j = 1,jd
c..grid

        xnew=(x(j,k)-xac)*cn+(y(j,k)-yac)*sn+xac
        ynew=-(x(j,k)-xac)*sn+(y(j,k)-yac)*cn+yac

	if(ntac.eq.1) then
 	ug(j,k)=(xnew-x(j,k))/dt
   	vg(j,k)=(ynew-y(j,k))/dt
	else
        xo=(x(j,k)-xac)*cnold+(y(j,k)-yac)*snold+xac
        yo=-(x(j,k)-xac)*snold+(y(j,k)-yac)*cnold+yac
 	ug(j,k)=(1.5*xnew-2*x(j,k)+0.5*xo)/dt
 	vg(j,k)=(1.5*ynew-2*y(j,k)+0.5*yo)/dt
	endif

        x(j,k)=xnew
        y(j,k)=ynew

  81	continue

      do k = 1,kmax
      do j = 1,jmax

c..store the old mesh for time metrics

        xole(j,k) = xold(j,k)
        yole(j,k) = yold(j,k)

        xold(j,k) = xv(j,k)
        yold(j,k) = yv(j,k)

c..update mesh

        xv(j,k) = (xold(j,k)-xac)*cn+(yold(j,k)-yac)*sn+xac
        yv(j,k) = -(xold(j,k)-xac)*sn+(yold(j,k)-yac)*cn+yac

      enddo
      enddo

	ang=atan2(y(jle,1)-y(jtail1,1),x(jtail1,1)-x(jle,1))

        !print*,'ang=',ang*180/pi


      return
      end
