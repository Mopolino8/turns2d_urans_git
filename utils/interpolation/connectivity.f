c********************************************************************
      module connectivity
        use grid
        implicit none

        integer :: idsize

        character*40,allocatable :: bc_file(:)
        integer,allocatable :: ndonor(:),nfringe(:),iisptr(:),iieptr(:)
        integer,allocatable :: idonor(:,:,:),imesh(:,:,:),ibc(:,:)
        real,allocatable :: frac(:,:,:)
        integer,allocatable :: iblank(:,:,:)

      contains

c********************************************************************
      subroutine connect
        use ihc
c********************************************************************

        if(.not.gridinit) call init_grid()
        call init_variables
        call connect2d(x,y,xc,yc,vol,iblank,idonor,frac,imesh,ibc, 
     &      iisptr,iieptr,ndonor,nfringe,jmax,kmax,nhalo,ngrids,bc_file)
        !call write_connect_data
      end subroutine connect

c********************************************************************
      subroutine init_variables
c..reading/setting the inputs
c********************************************************************
        
        integer j,k,n,jk

        allocate(bc_file(ngrids))

        bc_file(1) = 'interbc.inp'
        bc_file(2) = 'fineinterbc.inp'

        open(unit=1,file=bc_file(1))
        open(unit=2,file=bc_file(2))

        write(1,'(A)') '$IHCBCINP' 
        write(1,'(A)') '  NBC    =  1,'
        write(1,'(A)') '  IBTYP  =  1,'
        write(1,'(A)') '  IBDIR  =  1,'
        write(1,'(A)') '  JBCS   =  1,'
        write(1,'(A)') '  JBCE   = -1,'
        write(1,'(A)') '  KBCS   =  1,'
        write(1,'(A)') '  KBCE   = -1,'
        write(1,'(A)') '  IBPROC =  0,'
        write(1,'(A)') '$END'

        write(2,'(A)') '$IHCBCINP' 
        write(2,'(A)') '  NBC    =  1,'
        write(2,'(A)') '  IBTYP  =  3,'
        write(2,'(A)') '  IBDIR  =  1,'
        write(2,'(A)') '  JBCS   =  1,'
        write(2,'(A)') '  JBCE   = -1,'
        write(2,'(A)') '  KBCS   =  1,'
        write(2,'(A)') '  KBCE   = -1,'
        write(2,'(A)') '  IBPROC =  0,'
        write(2,'(A)') '$END'

        close(unit=1)
        close(unit=2)

        allocate(ndonor(ngrids),nfringe(ngrids))
        allocate(iisptr(ngrids),iieptr(ngrids))

        idsize = 0
        do n = 1,ngrids
          jk = jmax(n)*kmax(n)
          if(jk.gt.idsize) idsize = jk 
        enddo

        allocate(imesh(idsize,2,ngrids),idonor(idsize,2,ngrids))
        allocate(ibc(idsize,ngrids),frac(idsize,2,ngrids))

        allocate(iblank(jmc,kmc,ngrids))
        iblank = 1

      end subroutine init_variables

c********************************************************************
      subroutine write_connect_data
!...writing out the connectivity
c********************************************************************

        integer j,k,n,nf,nd,logi
        character*40,integer_string,tmp

        do n=1,ngrids
          write(integer_string,*) n-1
          integer_string = adjustl(integer_string)

          tmp = 'inter.'//trim(integer_string)
          open(unit=1,file=tmp,form='unformatted',status='unknown')
          
          logi = 1
          write(logi) nfringe(n),ndonor(n),iieptr(n),iisptr(n)
          write(logi) (idonor(nd,1,n),nd=1,ndonor(n)),
     &                (idonor(nd,2,n),nd=1,ndonor(n)),
     &                (frac(nd,1,n),nd=1,ndonor(n)),
     &                (frac(nd,2,n),nd=1,ndonor(n))
          write(logi) (imesh(nf,1,n),nf=1,nfringe(n)),
     &                (imesh(nf,2,n),nf=1,nfringe(n)),
     &                (ibc(nf,n),nf=1,nfringe(n))
          write(logi) ((iblank(j,k,n),j=1,jd(n)),k=1,kd(n))
          close(logi)
        enddo

      end subroutine write_connect_data

c********************************************************************
      end module connectivity
