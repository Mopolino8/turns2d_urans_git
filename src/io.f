c**********************************************************************
      subroutine astore(x,y,iblank,q,iturb,itrans,jgmx,kgmx,nmesh,
     &     ipointer,igrd,iqdim,jd,kd,istep0,ts,nspec)
c
c     make calls to store the solution and eddy viscosity 
c     in multizoned plot3d file format
c**********************************************************************      
      implicit none
c**********************************************************************      

      integer nmesh,igrd,iqdim,jd,kd,im
      real x(igrd),y(igrd),q(iqdim)
      integer iturb,itrans,nspec
      integer iblank(igrd)
      integer jgmx(nmesh),kgmx(nmesh)
      integer ipointer(nmesh,5)
      integer istep0
      logical ts
      character*40, gfile, qfile, turfile, integer_string
      
      ! local variables

      integer logg,logq,logtur
      integer lg(nspec),lq(nspec),ltur(nspec)
      integer i,nsp,ig,igq,igs,igv,igb

c***  first executable statement

      logg=9  ! grid     file => fort.9
      logq=8  ! solution file => fort.8
      logtur=10 ! turbulence file => SA.dat or komega.dat 

!      rewind(logg)
!      rewind(logq)

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(gfile,qfile,integer_string)
      spectralloop1: do nsp = 1,nspec

        if (ts) then
          write(integer_string,*) nsp
          gfile = 'fort_TS'//trim(adjustl(integer_string))//'.'
          qfile = 'fort_TS'//trim(adjustl(integer_string))//'.'
        else
          gfile = 'fort.'
          qfile = 'fort.'
        endif

        write(integer_string,*) logq
        qfile = trim(adjustl(qfile))//trim(adjustl(integer_string))

        write(integer_string,*) logg
        gfile = trim(adjustl(gfile))//trim(adjustl(integer_string))

        lg(nsp) = 1000*nsp + logg
        lq(nsp) = 1000*nsp + logq
        ltur(nsp) = 1000*nsp + logtur

        open(unit=lg(nsp),file=gfile,status='unknown',
     &                                 form='unformatted')
        open(unit=lq(nsp),file=qfile,status='unknown',
     &                                 form='unformatted')

        if(nmesh.gt.1) write(lg(nsp)) nmesh
        write(lg(nsp)) (jgmx(i),kgmx(i),i=1,nmesh)
        if(nmesh.gt.1) write(lq(nsp)) nmesh
        write(lq(nsp)) (jgmx(i),kgmx(i),i=1,nmesh)

      enddo spectralloop1
!$OMP END PARALLEL DO

      ! store solution

      do im=1,nmesh
         call set_pointers_globals(im,ipointer,ig,
     &        igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
         call store(q(igq),x(ig),y(ig),iblank(ig),jd,kd,lg,lq)
      enddo

      ! store turbulent quantities

      if(iturb.eq.1) then

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(turfile,integer_string)
      spectralloop2: do nsp = 1,nspec

        if (ts) then
          write(integer_string,*) nsp
          turfile = 'SA_TS'//trim(adjustl(integer_string))
        else
          turfile = 'SA'
        endif

        turfile = trim(adjustl(turfile))//'.dat'

        open(unit=ltur(nsp),file=turfile,status='unknown',
     c                                               form='formatted')
        if (itrans.eq.0) then
          write(ltur(nsp),*) "TITLE = ""Turbulence Model solution File"""
          write(ltur(nsp),*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""Vnu_t"""
        else
          write(ltur(nsp),*) "TITLE = ""Turbulence and Transition Model
     & solution File"""
          write(ltur(nsp),*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""Vnu_t""
     & ""Intermittancy"" ""Re_theta"""
        endif

      enddo spectralloop2
!$OMP END PARALLEL DO

        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,
     &         igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
          call storevnu(x(ig),y(ig),iblank(ig),q(igq),jd,kd,lq,ltur)
        enddo

      elseif(iturb.eq.2) then

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(turfile,integer_string)
      spectralloop3: do nsp = 1,nspec

        if (ts) then
          write(integer_string,*) nsp
          turfile = 'komega_TS'//trim(adjustl(integer_string))
        else
          turfile = 'komega'
        endif

        turfile = trim(adjustl(turfile))//'.dat'

        open(unit=ltur(nsp),file=turfile,status='unknown',
     c                                               form='formatted')
        if (itrans.eq.0) then
         write(ltur(nsp),*) "TITLE = ""Turbulence Model solution File"""
          write(ltur(nsp),*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""K""
     & ""Omega""" 
        else
          write(ltur(nsp),*) "TITLE = ""Turbulence and Transition Model
     & solution File"""
          write(ltur(nsp),*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""K""
     & ""Omega"" ""Intermittancy"" ""Re_theta"""
        endif

      enddo spectralloop3
!$OMP END PARALLEL DO

        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,
     &         igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
          call storekomega(x(ig),y(ig),iblank(ig),q(igq),jd,kd,lq,ltur)

        enddo
      endif

!$OMP PARALLEL DO IF(NSPEC > 1)
      spectralloop4: do nsp = 1,nspec

      write(lq(nsp)) istep0

      close(lg(nsp))
      close(lq(nsp))
      close(ltur(nsp))

      enddo spectralloop4
!$OMP END PARALLEL DO

      return
      end

c***********************************************************************
      subroutine storevnu(x,y,iblank,q,jd,kd,lq,ltur)

c      write the turbulent viscosity out
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,lq(nspec),ltur(nspec)
      real q(jd,kd,nspec,nq), x(jd,kd,nspec), y(jd,kd,nspec)
      integer iblank(jd,kd,nspec)

      ! local variables
      integer j,k,nsp,n
      integer logq,logtur

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(j,k,n,logq,logtur)
      spectralloop: do nsp = 1,nspec

      logq = lq(nsp) 
      logtur = ltur(nsp)

      !asitav !save in 'store', bug if run with overset 
      !am write(logq) (((q(j,k,nsp,n),j=1,jd),k=1,kd),n=5,nv)

      write(logtur,*) "ZONE"
      write(logtur,*) "I = ",jd,", J = ",kd
      write(logtur,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
      write(logtur,'(3F22.13)') ((x(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((y(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3I10)') ((iblank(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((q(j,k,nsp,5),j=1,jd),k=1,kd)
      if (itrans.eq.1) then
        write(logtur,'(3F22.13)') ((q(j,k,nsp,6),j=1,jd),k=1,kd)
        write(logtur,'(3F22.13)') ((q(j,k,nsp,7),j=1,jd),k=1,kd)
      endif

      enddo spectralloop
!$OMP END PARALLEL DO

      return
      end

c***********************************************************************
      subroutine storekomega(x,y,iblank,q,jd,kd,lq,ltur)

c      write the turbulent viscosity out
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,lq(nspec),ltur(nspec)
      real q(jd,kd,nspec,nq), x(jd,kd,nspec), y(jd,kd,nspec)
      integer iblank(jd,kd,nspec)

      ! local variables
      integer j,k,nsp,n
      integer logq,logtur

      real,allocatable :: rho(:,:)

!$OMP PARALLEL IF(NSPEC > 1)
!$OMP& PRIVATE(rho)
      allocate(rho(jd,kd))

!$OMP DO
!$OMP& PRIVATE(j,k,n,logq,logtur)
      spectralloop: do nsp = 1,nspec

      do j = 1,jd
      do k = 1,kd
        rho(j,k) = q(j,k,nsp,1)*q(j,k,nsp,nq)
      enddo
      enddo

      logq = lq(nsp) 
      logtur = ltur(nsp)

      write(logq) (((q(j,k,nsp,n),j=1,jd),k=1,kd),n=5,nv)

      write(logtur,*) "ZONE"
      write(logtur,*) "I = ",jd,", J = ",kd
      write(logtur,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
      write(logtur,'(3F22.13)') ((x(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((y(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3I10)') ((iblank(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((q(j,k,nsp,5)/rho(j,k),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((q(j,k,nsp,6)/rho(j,k),j=1,jd),k=1,kd)
      if (itrans.eq.1) then
        write(logtur,'(3F22.13)') ((q(j,k,nsp,7)/rho(j,k),j=1,jd),k=1,kd)
        write(logtur,'(3F22.13)') ((q(j,k,nsp,8)/rho(j,k),j=1,jd),k=1,kd)
      endif

      enddo spectralloop
!$OMP END DO

      deallocate(rho)
!$OMP END PARALLEL

      return
      end

c***********************************************************************
      subroutine store( q,x,y,iblank,jd,kd,lg,lq)
c
c  write the solution out
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,lg(nspec),lq(nspec)
      real q(jd,kd,nspec,nq), x(jd,kd,nspec), y(jd,kd,nspec)
      integer iblank(jd,kd,nspec)
      ! local variables
      integer j,k,nsp,n
      integer logg,logq
      real fstip,reypr

c***  first executable statement

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(j,k,n,logg,logq,fstip,reypr)
      spectralloop: do nsp = 1,nspec
      do 11 k = 1,kd
      do 11 j = 1,jd
        q(j,k,nsp,1) = q(j,k,nsp,1)*q(j,k,nsp,nq)
        q(j,k,nsp,2) = q(j,k,nsp,2)*q(j,k,nsp,nq)
        q(j,k,nsp,3) = q(j,k,nsp,3)*q(j,k,nsp,nq)
        q(j,k,nsp,4) = q(j,k,nsp,4)*q(j,k,nsp,nq)
   11 continue

      logg = lg(nsp)
      logq = lq(nsp) 

      fstip = fsmach
      reypr = rey * fsmach
      !am write(logq) fstip,alfa,reypr,totime
      write(logq) fstip,alfa,reypr,totime
      write(logq) (((q(j,k,nsp,n),j=1,jd),k=1,kd),n=1,4)

      !asitav !save here instead of storevnu, bug in overset otherwise
      if(iturb.eq.1 .or. iturb.eq.2) 
     >  write(logq) (((q(j,k,nsp,n),j=1,jd),k=1,kd),n=5,nv)

      if(num_grids.gt.1) then
        write(logg)((x(j,k,nsp),j=1,jd),k=1,kd),
     <             ((y(j,k,nsp),j=1,jd),k=1,kd),
     <             ((iblank(j,k,nsp),j=1,jd),k=1,kd)
      else
        write(logg)((x(j,k,nsp),j=1,jd),k=1,kd),
     <             ((y(j,k,nsp),j=1,jd),k=1,kd)
      endif
    

c..scale q back with jacobian

      call qdivj( q(:,:,nsp,:),jd,kd)
      
      enddo spectralloop
!$OMP END PARALLEL DO

      return
      end

c***********************************************************************
      subroutine grid(xv,yv,iblank,jd,kd,logg)
c
c  read grid from disk
c
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,logg
      real xv(jmax,kmax),yv(jmax,kmax)
      integer iblank(jd,kd)
      
      ! local variables
      integer j,k

c***  first executable statement

	if(num_grids.gt.1) then
      
      read(logg)((xv(j,k),j=1,jmax),k=1,kmax),
     <     ((yv(j,k),j=1,jmax),k=1,kmax)

	else
      read(logg)((xv(j,k),j=1,jmax),k=1,kmax),
     <     ((yv(j,k),j=1,jmax),k=1,kmax)

	do k=1,kd
		do j=1,jd
		iblank(j,k)=1
		enddo
	enddo


      endif
c
      return
      end

c***********************************************************************
      subroutine restr2( q,jd,kd,logq )
c
c  read the restart file from unit 3 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq
      real q(jd,kd,nq),reypr
      integer j,k
      ! local variables
      integer n
      real alf,fstip
c***  first executable statement

      read(logq) fstip,alf,reypr,totime
      read(logq) (((q(j,k,n),j=1,jd),k=1,kd),n=1,4)

      return
      end


c***********************************************************************
      subroutine restr_turb( q,jd,kd,logq )

c     read in turbulence variables from restart file
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq
      real q(jd,kd,nq)

      ! local variables
      integer j,k,n

      read(logq) (((q(j,k,n),j=1,jd),k=1,kd),n=5,nv)

      return
      end

c***********************************************************************
      subroutine movie(q,x,y,iblank,ug,vg,jd,kd,init,logg,logq)
c
c  this writes data out to file by planes
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init,logg,logq
      real x(jd,kd), y(jd,kd), ug(jd,kd), vg(jd,kd)
      integer iblank(jd,kd)
      real q(jd,kd,nq)
      
      ! local variables
      integer j,k
      real ppp(jd,kd)
c***********************************************************************

      if(init.eq.1) then
        write(logg) jd,kd,int((nsteps - istep0-isin)/nmovie)
        write(logq) jd,kd,int((nsteps - istep0-isin)/nmovie)
        write(logq) fsmach,alfa,rey,totime
      else
        do 1 k = 1, kd
        do 1 j = 1, jd
          ppp(j,k) = gm1*( q(j,k,4)-0.5*(q(j,k,2)**2+q(j,k,3)**2)
     <                                    /q(j,k,1))*q(j,k,nq)
    1   continue


	
        write(logg) ((x(j,k),j=1,jd),k=1,kd),
     <            ((y(j,k),j=1,jd),k=1,kd),
     <            ((0.0000001*float(istep0/nmovie),j=1,jd),k=1,kd),
     <	          ((iblank(j,k),j=1,jd),k=1,kd)	


        write(logq) ((q(j,k,1)*q(j,k,nq),j=1,jd),k=1,kd),
     <   ((q(j,k,2)/q(j,k,1),j=1,jd),k=1,kd),
     <   ((q(j,k,3)/q(j,k,1),j=1,jd),k=1,kd),
     <   ((ppp(j,k)/(q(j,k,1)*q(j,k,nq))**1.4,j=1,jd),k=1,kd),
     <   (((ppp(j,k)-1./1.4)/(0.5*fsmach*fsmach),j=1,jd),k=1,kd)
      endif

      return
      end


      subroutine write_reystress(q,turmu,x,y,xx,xy,yx,yy,jd,kd,im)
c*************************************************************
      use params_global
c*************************************************************
      implicit none
c*************************************************************
      integer jd,kd,im
      real q(jd,kd,nspec,nq)
      real turmu(jd,kd,nspec)
      real x(jd,kd,nspec),y(jd,kd,nspec)
      real xx(jd,kd,nspec),xy(jd,kd,nspec),yx(jd,kd,nspec)
      real yy(jd,kd,nspec)
!     local variables
      integer js,je,ks,ke
      real, allocatable :: tauxx(:,:),tauxy(:,:),tauyy(:,:)
      real, allocatable :: u(:,:),v(:,:)
      real, allocatable :: ux(:,:),uy(:,:),vx(:,:),vy(:,:)
      real, allocatable :: vmul(:,:)
      integer j,k,jm1,km1,jp,kp,k3,k2,nsp,logrs1,logrs2
      real usi,vsi,ueta,veta,sfdiv,rei
      character*60 :: int_str,fileplt,fileg

      allocate(tauxx(jd,kd),tauxy(jd,kd),tauyy(jd,kd))
      allocate(u(jd,kd),v(jd,kd))
      allocate(ux(jd,kd),uy(jd,kd),vx(jd,kd),vy(jd,kd))
      allocate(vmul(jd,kd))

      do nsp=1,nspec

         call lamvis(q(:,:,nsp,:),vmul,jd,kd)
         do j=1, jd
            do k=1,kd
               u(j,k) = q(j,k,nsp,2)/q(j,k,nsp,1)
               v(j,k) = q(j,k,nsp,3)/q(j,k,nsp,1)
            end do
         end do
         tauxy(:,:) = 0.0
         ks = nhalo+1
         ke = kd - nhalo
         
         js = nhalo+1
         je = jd - nhalo
         
         rei = 1./rey
         do j = js,je 
            jp = j + 1
            jm1 = j - 1
            do k = ks,ke
               k2 = k + 1
               k3 = k2 + 1
               kp = k + 1
               km1 = k - 1
               usi  = 0.5*(u(jp,k)-u(jm1,k))
               vsi  = 0.5*(v(jp,k)-v(jm1,k))
               
               ueta = 0.5*(u(j,kp)-u(j,km1))
               veta = 0.5*(v(j,kp)-v(j,km1))
               
               ux(j,k) = xx(j,k,nsp)*usi + yx(j,k,nsp)*ueta
               uy(j,k) = xy(j,k,nsp)*usi + yy(j,k,nsp)*ueta
               
               vx(j,k) = xx(j,k,nsp)*vsi + yx(j,k,nsp)*veta
               vy(j,k) = xy(j,k,nsp)*vsi + yy(j,k,nsp)*veta
               
               sfdiv = 2./3*(ux(j,k) + vy(j,k))
               tauxy(j,k) = -turmu(j,k,nsp) * ( uy(j,k) + vx(j,k) )* rei/fsmach**2
               
            enddo
         enddo
         
         if(nspec.le.1) then
           fileplt='reystress.plt'
           fileg='reystress.g'
           logrs1 = 97000+nspec
           logrs2 = 98000+nspec
         else
           write(int_str,'(i7)')nsp
           fileplt='reystress_TS'//trim(adjustl(int_str))//'.plt'
           fileg='reystress_TS'//trim(adjustl(int_str))//'.g'
           logrs1 = 97000+nsp
           logrs2 = 98000+nsp
         end if
         
         if(im.eq.1) then
           open(logrs1, file = fileplt, form = 'formatted')
           write(logrs1,*) "TITLE = ""Turbulence Model solution File"""
           write(logrs1,*) 
     >       "VARIABLES = ""X"" ""Y"" ""u"" ""v"" ""reystress"""
         endif !im
         write(logrs1,*) "ZONE"
         write(logrs1,*) "I = ",jd,", J = ",kd
         write(logrs1,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
         write(logrs1,'(F22.13)') ((x(j,k,nsp),j=1,jd),k=1,kd)
         write(logrs1,'(F22.13)') ((y(j,k,nsp),j=1,jd),k=1,kd)
         write(logrs1,*) ((u(j,k),j=1,jd),k=1,kd)
         write(logrs1,*) ((v(j,k),j=1,jd),k=1,kd)
         write(logrs1,*) ((tauxy(j,k),j=1,jd),k=1,kd)
         if(im.eq.num_grids) close(logrs1)
         
         if(im.eq.1) open(logrs2, file = fileg, form = 'unformatted')
         write(logrs2) jd, kd
         write(logrs2) ((x(j,k,nsp),j=1,jd),k=1,kd)
         write(logrs2) ((y(j,k,nsp),j=1,jd),k=1,kd)
         write(logrs2) ((u(j,k),j=1,jd),k=1,kd)
         write(logrs2) ((v(j,k),j=1,jd),k=1,kd)
         write(logrs2) ((tauxy(j,k),j=1,jd),k=1,kd)
         if(im.eq.num_grids) close(logrs2)

      end do
      
      deallocate(tauxx,tauxy,tauyy,u,v,ux,uy,vx,vy,vmul)
      return
      end

c************************************************************************

