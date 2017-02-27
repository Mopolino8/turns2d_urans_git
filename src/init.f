c*********************************************************************
      subroutine initia(q,s,x,y,xv,yv,xx,xy,yx,yy,ug,vg,ugv,vgv,turmu,
     &     xold,yold,xole,yole,iblank,tspec,im,jd,kd)

c
c  initialize the variables to default values, read inputs, check them
c  and set up the working files.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,im
      
      real q(jd,kd,nspec,nq),s(jd,kd,nspec,nv)
      real x(jd,kd,nspec),y(jd,kd,nspec),xv(jmax,kmax,nspec),yv(jmax,kmax,nspec)
      real xx(jd,kd,nspec),xy(jd,kd,nspec),yx(jd,kd,nspec),yy(jd,kd,nspec)
      real ug(jd,kd,nspec),vg(jd,kd,nspec)
      real ugv(jmax,kmax,nspec),vgv(jmax,kmax,nspec)
      real turmu(jd,kd,nspec)
      integer iblank(jd,kd,nspec)
      real xold(jmax,kmax,nspec),yold(jmax,kmax,nspec)
      real xole(jmax,kmax,nspec),yole(jmax,kmax,nspec)
      real tspec(nspec)

      ! local variables

      integer i,j,k,jtmp,ktmp,nmesh,nsp,ib
      character*40, qfile, integer_string
      integer logq
      integer lq(nspec)

c**   first executable statement


      pi = 4.0*atan(1.)
      dang = 0.
      do 882 nsp=1,nspec
        do 882 j = 1, jd
          do 882 k=1,kd
            turmu(j,k,nsp) = 0.
            ug(j,k,nsp) = 0.
            vg(j,k,nsp) = 0.
            iblank(j,k,nsp)=1
882   continue
      
      do 884 nsp=1,nspec
        do 884 j=1,jd
          do 884 k=1,kd
            s(j,k,nsp,1) = 0.
            s(j,k,nsp,2) = 0.
            s(j,k,nsp,3) = 0.
            s(j,k,nsp,4) = 0.
884   continue

c..setup bodyflag

      bodyflag(im) = .FALSE.
      jtail(im) = 1 
      do ib=1,nbc_all(im)
        if (ibtyp_all(ib,im).eq.5) then
          bodyflag(im) = .TRUE.
          jtail(im) = jbcs_all(ib,im)
        endif
      enddo

c..setup q and x

      istep0  = 0
      do nsp=1,nspec
        !amm tspec(nsp) = 2*pi*(nsp-1)/nspec/abs(rf)
        if(rf.ne.0) then
          tspec(nsp) = 2*pi*(nsp-1)/nspec/abs(rf)
        elseif(rf_tef.ne.0 .and. motiontype.eq.3) then
          tspec(nsp) = 2*pi*(nsp-1)/nspec/abs(rf_tef)
        end if
      enddo

	print*,'NUM GRIDS=========== ',num_grids

!      if (iread.eq.0) then
         call grid(xv(:,:,1),yv(:,:,1),iblank(:,:,1),jd,kd,1)
!      else
!        call grid(xv(:,:,1),yv(:,:,1),iblank(:,:,1),jd,kd,4)
!      endif

!$OMP PARALLEL IF(NSPEC > 1)
!$OMP DO
      do nsp = 2,nspec
        call find_coord_TS(xv(:,:,nsp),yv(:,:,nsp),xv(:,:,1),yv(:,:,1),tspec(nsp),jd,kd)
      enddo
!$OMP END DO

!$OMP DO
!$OMP& PRIVATE(j,k,logq,qfile,integer_string)
      spectralloop: do nsp = 1,nspec

c..initialize mesh from earlier times 

        do j = 1,jmax
          do k = 1,kmax
            xold(j,k,nsp) = xv(j,k,nsp)
            yold(j,k,nsp) = yv(j,k,nsp)
            xole(j,k,nsp) = xold(j,k,nsp)
            yole(j,k,nsp) = yold(j,k,nsp)
          enddo
        enddo

        call find_cellcenter(x(:,:,nsp),y(:,:,nsp),
     &                       xv(:,:,nsp),yv(:,:,nsp),jd,kd,im)

        if(iread.eq.0.and.iunst.gt.0) then !called only in unsteady simulation
          if (motiontype.eq.1) then
c           call move_new(1,x,y,xv,yv,xold,yold,xole,yole,ug,vg,jd,kd)
          elseif (motiontype.eq.2) then
            call move_cyclo(1,x,y,xv,yv,xold,yold,xole,yole,ug,vg,jd,kd,im)
          elseif (motiontype.eq.3.and.im.eq.1) then
            call move_tef(1,x,y,xv,yv,xold,yold,xole,yole,ug,vg,jd,kd,im)
          elseif(im.eq.1) then
            stop 'Error: Unknown motion type'
          endif
        endif
 
c..compute the metrics and the jacobians using finite volume formula

        call metfv(q(:,:,nsp,:),x(:,:,nsp),y(:,:,nsp),
     &             xv(:,:,nsp),yv(:,:,nsp),
     &             xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),jd,kd,im)
      
        if (timespectral) then
          call find_timemetrics_TS(x(:,:,nsp),y(:,:,nsp),
     &                             xv(:,:,nsp),yv(:,:,nsp),tspec(nsp),
     &                             ug(:,:,nsp),vg(:,:,nsp),
     &                             ugv(:,:,nsp),vgv(:,:,nsp),jd,kd)
        endif
          
        call qzero( q(:,:,nsp,:),jd,kd)

        if( iread .gt. 0 ) then

          logq = 3
          if (timespectral) then
            write(integer_string,*) nsp
            qfile = 'fort_TS'//trim(adjustl(integer_string))//'.'
          else
            qfile = 'fort.'
          endif

          write(integer_string,*) logq
          qfile = trim(adjustl(qfile))//trim(adjustl(integer_string))

          lq(nsp) = 1000*nsp + logq
          open(unit=lq(nsp),file=qfile,status='unknown',
     &                                 form='unformatted')

          if (im.eq.1) then
            if(num_grids.gt.1) read(lq(nsp)) nmesh
            read(lq(nsp)) (jtmp,ktmp,i=1,num_grids)
          endif

          call restr2( q(:,:,nsp,:),jd,kd,lq(nsp) )
          if(iturb.eq.1.or.iturb.eq.2) call restr_turb(q(:,:,nsp,:),jd,kd,lq(nsp))

          if (im.eq. num_grids) read(lq(nsp)) istep0

        endif

c..divide q by the jacobian, q/q5 

        call qdivj( q(:,:,nsp,:),jd,kd )
      enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL
 
      return
      end

c***********************************************************************
      subroutine qdivj( q,jd,kd )
c
c  divide q by jacobian
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)

      ! local variables
      integer j,k

c***  first executable statement

      do 71 k = 1,kd
      do 71 j = 1,jd
        q(j,k,1) = q(j,k,1)/q(j,k,nq)
        q(j,k,2) = q(j,k,2)/q(j,k,nq)
        q(j,k,3) = q(j,k,3)/q(j,k,nq)
        q(j,k,4) = q(j,k,4)/q(j,k,nq)
   71 continue

      return
      end

c***********************************************************************
      subroutine qzero( q,jd,kd )
c
c  set initial values of the flow variables to their free-stream
c  values for an impulsive start. otherwise the initial values
c  are read from a restart file.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)

      ! local variables
      real ruinf,rvinf,rtkeinf,rtomegainf,tu,re_theta
      integer j,k

c**   first executable statement

      ruinf = rinf*uinf
      rvinf = rinf*vinf
      do 11 k=1,kd
      do 11 j=1,jd
        q(j,k,1) = rinf
        q(j,k,2) = ruinf
        q(j,k,3) = rvinf
        q(j,k,4) = einf
  11  continue

      if (iturb.eq.1) then
        do k=1,kd
        do j=1,jd
          q(j,k,5) = vnuinf
        enddo
        enddo
       
      elseif (iturb.eq.2) then
        rtkeinf = rinf*tkeinf
        rtomegainf = rinf*tomegainf
        do k=1,kd
        do j=1,jd
          q(j,k,5) = rtkeinf
          q(j,k,6) = rtomegainf
        enddo
        enddo

        if (itrans.eq.1) then
          itmcinf = 1. 

          tu = 100.*sqrt(2./3*tkeinf)/fsmach
          if(tu.le.1.3) then
            re_theta = (1173.51-589.428*tu+0.2196/(tu*tu))
          else
            re_theta = 331.5*(tu-0.5658)**(-0.671)
          endif
          retinf = max(re_theta,20.)/rey

          do k=1,kd
          do j=1,jd
            q(j,k,7) = rinf*itmcinf
            q(j,k,8) = rinf*retinf
          enddo
          enddo
        endif
      endif
c
      return
      end

c*************************************************************************
