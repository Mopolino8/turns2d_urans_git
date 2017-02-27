c**********************************************************************
      subroutine astore(x,y,iblank,q,sb,iturb,itrans,jgmx,kgmx,nmesh,
     &     ipointer,igrd,iqdim,isdim,jd,kd,istep0,ts,nspec)
c
c     make calls to store the solution and eddy viscosity 
c     in multizoned plot3d file format
c**********************************************************************      
      implicit none
c**********************************************************************      

      integer nmesh,igrd,iqdim,isdim,jd,kd,im
      real x(igrd),y(igrd),q(iqdim),sb(isdim)
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

      logg=19  ! grid     file => fort.9
      logq=18  ! solution file => fort.8
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

!        open(unit=lg(nsp),file=gfile,status='unknown',
!     &                                 form='unformatted')
        open(unit=lq(nsp),file=qfile,status='unknown',
     &                                 form='unformatted')

!        if(nmesh.gt.1) write(lg(nsp)) nmesh
!        write(lg(nsp)) (jgmx(i),kgmx(i),i=1,nmesh)
        if(nmesh.gt.1) write(lq(nsp)) nmesh
        write(lq(nsp)) (jgmx(i),kgmx(i),i=1,nmesh)

      enddo spectralloop1
!$OMP END PARALLEL DO

      ! store solution

      do im=1,nmesh
         call set_pointers_globals(im,ipointer,ig,
     &        igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
         call store(q(igq),sb(igs),x(ig),y(ig),iblank(ig),jd,kd,lg,lq)
      enddo

      ! store turbulent quantities

      if(iturb.eq.1) then

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(turfile,integer_string)
      spectralloop2: do nsp = 1,nspec

        if (ts) then
          write(integer_string,*) nsp
          turfile = 'SAadj_TS'//trim(adjustl(integer_string))
        else
          turfile = 'SAadj'
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
          call storevnu(x(ig),y(ig),iblank(ig),q(igq),sb(igs),jd,kd,
     &                  lq,ltur)
        enddo

      elseif(iturb.eq.2) then

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(turfile,integer_string)
      spectralloop3: do nsp = 1,nspec

        if (ts) then
          write(integer_string,*) nsp
          turfile = 'komegaAdj_TS'//trim(adjustl(integer_string))
        else
          turfile = 'komegaAdj'
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
          call storekomega(x(ig),y(ig),iblank(ig),q(igq),sb(igs),jd,kd,lq,ltur)

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
      subroutine storevnu(x,y,iblank,q,sb,jd,kd,lq,ltur)

c      write the turbulent viscosity out
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,lq(nspec),ltur(nspec)
      real q(jd,kd,nspec,nq), sb(jd,kd,nspec,nv)
      real x(jd,kd,nspec), y(jd,kd,nspec)
      integer iblank(jd,kd,nspec)

      ! local variables
      integer j,k,nsp,n
      integer logq,logtur

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(j,k,n,logq,logtur)
      spectralloop: do nsp = 1,nspec

      logq = lq(nsp) 
      logtur = ltur(nsp)

      !am !done in store subroutine
      !am write(logq) (((sb(j,k,nsp,n),j=1,jd),k=1,kd),n=5,nv)
      !am !print*,'nsp,sum(sb) ',nsp,sum(sb(:,:,nsp,5))

      write(logtur,*) "ZONE"
      write(logtur,*) "I = ",jd,", J = ",kd
      write(logtur,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
      write(logtur,'(3F22.13)') ((x(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((y(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3I10)') ((iblank(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3E22.13)') ((sb(j,k,nsp,5),j=1,jd),k=1,kd)
      if (itrans.eq.1) then
        write(logtur,'(3F22.13)') ((sb(j,k,nsp,6),j=1,jd),k=1,kd)
        write(logtur,'(3F22.13)') ((sb(j,k,nsp,7),j=1,jd),k=1,kd)
      endif

      enddo spectralloop
!$OMP END PARALLEL DO

      return
      end

c***********************************************************************
      subroutine storekomega(x,y,iblank,q,sb,jd,kd,lq,ltur)

c      write the turbulent viscosity out
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,lq(nspec),ltur(nspec)
      real q(jd,kd,nspec,nq), sb(jd,kd,nspec,nv)
      real x(jd,kd,nspec), y(jd,kd,nspec)
      integer iblank(jd,kd,nspec)

      ! local variables
      integer j,k,nsp,n
      integer logq,logtur

!$OMP PARALLEL IF(NSPEC > 1)
!$OMP DO
!$OMP& PRIVATE(j,k,n,logq,logtur)
      spectralloop: do nsp = 1,nspec

      logq = lq(nsp) 
      logtur = ltur(nsp)

      write(logq) (((sb(j,k,nsp,n),j=1,jd),k=1,kd),n=5,nv)

      write(logtur,*) "ZONE"
      write(logtur,*) "I = ",jd,", J = ",kd
      write(logtur,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
      write(logtur,'(3F22.13)') ((x(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((y(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3I10)') ((iblank(j,k,nsp),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((sb(j,k,nsp,5),j=1,jd),k=1,kd)
      write(logtur,'(3F22.13)') ((sb(j,k,nsp,6),j=1,jd),k=1,kd)
      if (itrans.eq.1) then
        write(logtur,'(3F22.13)') ((sb(j,k,nsp,7),j=1,jd),k=1,kd)
        write(logtur,'(3F22.13)') ((sb(j,k,nsp,8),j=1,jd),k=1,kd)
      endif

      enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL

      return
      end

c***********************************************************************
      subroutine debugstore(q,qb,jd,kd,logq)
c
c  write the solution out
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,logq
      real q(jd,kd,nq), qb(jd,kd,nv)
      ! local variables
      integer j,k,n
      real fstip,reypr

c***  first executable statement

      fstip = fsmach
      reypr = rey * fsmach
      write(logq) jd,kd
      write(logq) fstip,alfa,reypr,totime
      write(logq) (((qb(j,k,n),j=1,jd),k=1,kd),n=1,4)

      close(logq)

      return
      end

c***********************************************************************
      subroutine store( q,sb,x,y,iblank,jd,kd,lg,lq)
c
c  write the solution out
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,lg(nspec),lq(nspec)
      real q(jd,kd,nspec,nq), sb(jd,kd,nspec,nv)
      real x(jd,kd,nspec), y(jd,kd,nspec)
      integer iblank(jd,kd,nspec)
      ! local variables
      integer j,k,nsp,n
      integer logg,logq
      real fstip,reypr

c***  first executable statement

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(j,k,n,logg,logq,fstip,reypr)
      spectralloop: do nsp = 1,nspec
!      do k = 1,kd
!      do j = 1,jd
!        sb(j,k,1) = sb(j,k,1)/q(j,k,nq)
!        sb(j,k,2) = sb(j,k,2)/q(j,k,nq)
!        sb(j,k,3) = sb(j,k,3)/q(j,k,nq)
!        sb(j,k,4) = sb(j,k,4)/q(j,k,nq)
!      enddo
!      enddo

      logg = lg(nsp)
      logq = lq(nsp) 

      fstip = fsmach
      reypr = rey * fsmach
      write(logq) fstip,alfa,reypr,totime
      write(logq) (((sb(j,k,nsp,n),j=1,jd),k=1,kd),n=1,4)
      print*,'nsp,sum(sb) ',nsp,sum(sb(:,:,nsp,5))
      
      !asitav !save here instead of storevnu, bug in overset otherwise
      if(iturb.eq.1 .or. iturb.eq.2) 
     >  write(logq) (((sb(j,k,nsp,n),j=1,jd),k=1,kd),n=5,nv)

!	if(num_grids.gt.1) then
!      write(logg)((x(j,k),j=1,jd),k=1,kd),
!     <           ((y(j,k),j=1,jd),k=1,kd),
!     <           ((iblank(j,k),j=1,jd),k=1,kd)
!	else
!      write(logg)((x(j,k),j=1,jd),k=1,kd),
!     <           ((y(j,k),j=1,jd),k=1,kd)
!	endif
    

c..scale q back with jacobian

!      do k = 1,kd
!      do j = 1,jd
!        sb(j,k,1) = sb(j,k,1)*q(j,k,nq)
!        sb(j,k,2) = sb(j,k,2)*q(j,k,nq)
!        sb(j,k,3) = sb(j,k,3)*q(j,k,nq)
!        sb(j,k,4) = sb(j,k,4)*q(j,k,nq)
!      enddo
!      enddo

      enddo spectralloop
!$OMP END PARALLEL DO

      return
      end

c***********************************************************************
      subroutine read_solution( q,jd,kd,logq )
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
      integer alf,fstip
      character*40, qfile, integer_string
c***  first executable statement

      read(logq) fstip,alf,reypr,totime
      read(logq) (((q(j,k,n),j=1,jd),k=1,kd),n=1,4)
      if (iturb.eq.1.or.iturb.eq.2) 
     &         read(logq) (((q(j,k,n),j=1,jd),k=1,kd),n=5,nv)

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
      subroutine restr2( sb,jd,kd,logq )
c
c  read the restart file from unit 3 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,nsp,logq
      real sb(jd,kd,nv),reypr
      integer j,k
      ! local variables
      integer n
      real alf,fstip
c***  first executable statement

      read(logq) fstip,alf,reypr,totime
      read(logq) (((sb(j,k,n),j=1,jd),k=1,kd),n=1,4)

      return
      end

c***********************************************************************
      subroutine restr_turb( sb,jd,kd,logq )

c     read in eddy viscosity values from restart file
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq,j,k,n
      real sb(jd,kd,nv)

      read(logq) (((sb(j,k,n),j=1,jd),k=1,kd),n=5,nv)

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

c************************************************************************
