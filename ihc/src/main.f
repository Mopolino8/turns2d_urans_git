c********************************************************************
      module global
        integer ngrids
        integer,allocatable :: jmax(:),kmax(:)
        integer,allocatable :: jd(:),kd(:)
        integer,allocatable :: jbeg(:),kbeg(:)
        integer,allocatable :: jend(:),kend(:)
        integer :: nhalo
        integer,allocatable :: iblank(:,:,:)
        real,allocatable :: x(:,:,:),y(:,:,:)
        real,allocatable :: xc(:,:,:),yc(:,:,:)
        real,allocatable :: vol(:,:,:)
        character*40,allocatable :: grid_file(:),bc_file(:) !relocate
        integer,allocatable :: ndonor(:),nfringe(:),iisptr(:),iieptr(:)
        integer,allocatable :: idonor(:,:,:),imesh(:,:,:),ibc(:,:)
        real,allocatable :: frac(:,:,:)
      end module global

c********************************************************************
      program main
        use global
        use ihc
        implicit none
c********************************************************************

        call init_variables
        call connect2d(x,y,xc,yc,vol,iblank,idonor,frac,imesh,ibc, 
     &      iisptr,iieptr,ndonor,nfringe,jmax,kmax,nhalo,ngrids,bc_file)
        !call write_connect_data
        call write_grid
      end program main

c********************************************************************
      subroutine init_variables
c..reading/setting the inputs
c
        use global
        implicit none

c********************************************************************
        
        integer j,k,n,jm,km,jmc,kmc,jk,idsize
        integer logg,logi
c        character*40,allocatable :: grid_file(:)

        print*,'enter the number of grid files'
        read(5,*) ngrids
        allocate(grid_file(ngrids),bc_file(ngrids))

        print*, 'enter the grid files'
        do n = 1,ngrids
         read(5,*) grid_file(n) 
         grid_file(n) = adjustl(grid_file(n))
        enddo

        print*, 'enter the bc files'
        do n = 1,ngrids
         read(5,*) bc_file(n) 
         bc_file(n) = adjustl(bc_file(n))
        enddo

c..processing the grid

c....determining the dimensions for grid

        allocate(jmax(ngrids),kmax(ngrids))
        allocate(jd(ngrids),kd(ngrids))
        allocate(jbeg(ngrids),jend(ngrids))
        allocate(kbeg(ngrids),kend(ngrids))

        nhalo = 2

        do n = 1,ngrids
          open(unit=1,file=grid_file(n),form='unformatted',
     >                                    status='unknown')
          logg=1
          read(logg) jmax(n),kmax(n)
          rewind(logg)
          close(logg)
        
          jd(n) = jmax(n) + 2*nhalo - 1
          kd(n) = kmax(n) + 2*nhalo - 1

          jbeg(n) = nhalo + 1
          kbeg(n) = nhalo + 1
          jend(n) = jd(n)-nhalo
          kend(n) = kd(n)-nhalo

        enddo

        jm = 0
        km = 0
        do n = 1,ngrids
          if(jmax(n).gt.jm) jm = jmax(n) 
          if(kmax(n).gt.km) km = kmax(n) 
        enddo
 
        jmc = 0
        kmc = 0
        do n = 1,ngrids
          if(jd(n).gt.jmc) jmc = jd(n) 
          if(kd(n).gt.kmc) kmc = kd(n) 
        enddo
 
c....reading the grid

        allocate(x(jm,km,ngrids),y(jm,km,ngrids))
        allocate(xc(jmc,kmc,ngrids),yc(jmc,kmc,ngrids))
        allocate(iblank(jmc,kmc,ngrids))

        x = 0.; y = 0.; xc = 0.; yc = 0.; iblank = 1
        do n = 1,ngrids
          open(unit=1,file=grid_file(n),form='unformatted',
     >                                    status='unknown')
          logg=1
          read(logg) jmax(n),kmax(n)
          read(logg) ((x(j,k,n),j=1,jmax(n)),k=1,kmax(n)),
     >               ((y(j,k,n),j=1,jmax(n)),k=1,kmax(n))
          close(logg)

        enddo
        call find_cellcenter()
        
c...compute cell volume

        allocate(vol(jmc,kmc,ngrids))
        call metfv()

c...initialize some of the connectivity variables

        allocate(ndonor(ngrids),nfringe(ngrids))
        allocate(iisptr(ngrids),iieptr(ngrids))

        idsize = 0
        do n = 1,ngrids
          jk = jmax(n)*kmax(n)
          if(jk.gt.idsize) idsize = jk 
        enddo

        allocate(imesh(idsize,2,ngrids),idonor(idsize,2,ngrids))
        allocate(ibc(idsize,ngrids),frac(idsize,2,ngrids))

      end subroutine init_variables

c*********************************************************************
      subroutine find_cellcenter()
c.. find cell center
c*********************************************************************
      use global
      implicit none
c*********************************************************************
      integer j,k,n,jv,kv,nh 

      do n = 1,ngrids
        do j = jbeg(n),jend(n)
          do k = kbeg(n),kend(n)
            jv = j - nhalo + 1
            kv = k - nhalo + 1
            xc(j,k,n) = 0.25*(x(jv,kv,n)   + x(jv-1,kv,n)+
     &                        x(jv,kv-1,n) + x(jv-1,kv-1,n))
            yc(j,k,n) = 0.25*(y(jv,kv,n)   + y(jv-1,kv,n)+
     &                        y(jv,kv-1,n) + y(jv-1,kv-1,n))
          enddo
        enddo

!...boundary values (need not be accurate - should not be used anywhere) 

        j = jbeg(n) - 1
        do k = kbeg(n),kend(n)
          jv = j - nhalo + 1
          kv = k - nhalo + 1
          xc(j,k,n) = xc(j+1,k,n) - 0.5*( x(jv+1,kv,n)   - x(jv,kv,n) + 
     &                                    x(jv+1,kv-1,n) - x(jv,kv-1,n))
          yc(j,k,n) = yc(j+1,k,n) - 0.5*( y(jv+1,kv,n)   - y(jv,kv,n) + 
     &                                    y(jv+1,kv-1,n) - y(jv,kv-1,n))
        enddo

        j = jend(n) + 1
        do k = kbeg(n),kend(n)
          jv = j - nhalo + 1
          kv = k - nhalo + 1
          xc(j,k,n) = xc(j-1,k,n)-0.5*(x(jv-2,kv,n)   - x(jv-1,kv,n) + 
     &                                 x(jv-2,kv-1,n) - x(jv-1,kv-1,n))
          yc(j,k,n) = yc(j-1,k,n)-0.5*(y(jv-2,kv,n)   - y(jv-1,kv,n) + 
     &                                 y(jv-2,kv-1,n) - y(jv-1,kv-1,n))
        enddo
 
        k = kbeg(n) - 1
        do j = jbeg(n),jend(n)
          jv = j - nhalo + 1
          kv = k - nhalo + 1
          xc(j,k,n) = xc(j,k+1,n)-0.5*(x(jv,kv+1,n) - x(jv,kv,n) + 
     &                                 x(jv-1,kv+1,n) - x(jv-1,kv,n))
          yc(j,k,n) = yc(j,k+1,n)-0.5*(y(jv,kv+1,n) - y(jv,kv,n) + 
     &                                 y(jv-1,kv+1,n) - y(jv-1,kv,n))
        enddo

        k = kend(n) + 1
        do j = jbeg(n),jend(n)
          jv = j - nhalo + 1
          kv = k - nhalo + 1
          xc(j,k,n) = xc(j,k-1,n)-0.5*(x(jv,kv-2,n)   - x(jv,kv-1,n) + 
     &                                 x(jv-1,kv-2,n) - x(jv-1,kv-1,n))
          yc(j,k,n) = yc(j,k-1,n)-0.5*(y(jv,kv-2,n)   - y(jv,kv-1,n) + 
     &                                 y(jv-1,kv-2,n) - y(jv-1,kv-1,n))
        enddo

        j = jbeg(n) - 1
        k = kbeg(n) - 1
        xc(j,k,n) = 2*xc(j+1,k,n) - xc(j+2,k,n)
        yc(j,k,n) = 2*yc(j+1,k,n) - yc(j+2,k,n)

        k = kend(n) + 1
        xc(j,k,n) = 2*xc(j+1,k,n) - xc(j+2,k,n)
        yc(j,k,n) = 2*yc(j+1,k,n) - yc(j+2,k,n)

        j = jend(n) + 1
        k = kbeg(n) - 1
        xc(j,k,n) = 2*xc(j,k+1,n) - xc(j,k+2,n)
        yc(j,k,n) = 2*yc(j,k+1,n) - yc(j,k+2,n)

        k = kend(n) + 1
        xc(j,k,n) = 2*xc(j,k-1,n) - xc(j,k-2,n)
        yc(j,k,n) = 2*yc(j,k-1,n) - yc(j,k-2,n)

        do nh = 2,nhalo
          j = jbeg(n) - nh
          do k = kbeg(n) - nh + 1,kend(n) + nh - 1
            xc(j,k,n) = 2*xc(j+1,k,n) - xc(j+2,k,n)
            yc(j,k,n) = 2*yc(j+1,k,n) - yc(j+2,k,n)
          enddo

          j = jend(n) + nh
          do k = kbeg(n) - nh + 1,kend(n) + nh - 1
            xc(j,k,n) = 2*xc(j-1,k,n) - xc(j-2,k,n)
            yc(j,k,n) = 2*yc(j-1,k,n) - yc(j-2,k,n)
          enddo
 
          k = kbeg(n) - nh
          do j = jbeg(n) - nh + 1,jend(n) + nh - 1
            xc(j,k,n) = 2*xc(j,k+1,n) - xc(j,k+2,n)
            yc(j,k,n) = 2*yc(j,k+1,n) - yc(j,k+2,n)
          enddo

          k = kend(n) + nh
          do j = jbeg(n) - nh + 1,jend(n) + nh - 1
            xc(j,k,n) = 2*xc(j,k-1,n) - xc(j,k-2,n)
            yc(j,k,n) = 2*yc(j,k-1,n) - yc(j,k-2,n)
          enddo

          j = jbeg(n) - nh
          k = kbeg(n) - nh
          xc(j,k,n) = 2*xc(j+1,k,n) - xc(j+2,k,n)
          yc(j,k,n) = 2*yc(j+1,k,n) - yc(j+2,k,n)

          k = kend(n) + nh
          xc(j,k,n) = 2*xc(j+1,k,n) - xc(j+2,k,n)
          yc(j,k,n) = 2*yc(j+1,k,n) - yc(j+2,k,n)

          j = jend(n) + nh
          k = kbeg(n) - nh
          xc(j,k,n) = 2*xc(j,k+1,n) - xc(j,k+2,n)
          yc(j,k,n) = 2*yc(j,k+1,n) - yc(j,k+2,n)

          k = kend(n) + nh
          xc(j,k,n) = 2*xc(j,k-1,n) - xc(j,k-2,n)
          yc(j,k,n) = 2*yc(j,k-1,n) - yc(j,k-2,n)

        enddo
      enddo

      end subroutine find_cellcenter

c***********************************************************************
      subroutine metfv()
c***********************************************************************
      use global
      implicit none
c***********************************************************************
      integer j,k,n,jc,kc,j1,k1,nh,nneg
      real dx1,dy1,dx2,dy2

      do n = 1,ngrids
c..find the volume near j,k

        do k=2,kmax(n)
          do j=2,jmax(n)
            jc = j + nhalo - 1 
            kc = k + nhalo - 1 
            dx1 = x(j,k,n)-x(j-1,k-1,n)
            dy1 = y(j,k,n)-y(j-1,k-1,n)
            dx2 = x(j-1,k,n)-x(j,k-1,n)
            dy2 = y(j-1,k,n)-y(j,k-1,n)
            vol(jc,kc,n) = 0.5*( dx1*dy2 -dx2*dy1 )
          enddo
        enddo

c..boundary values
c
        do nh = 1,nhalo
          j = jbeg(n) - nh
          j1 = jbeg(n) + nh - 1
          do k = kbeg(n) - nh + 1,kend(n) + nh - 1
            vol(j,k,n) = vol(j1,k,n)
          enddo

          j = jend(n) + nh
          j1 = jend(n) - nh + 1
          do k = kbeg(n) - nh + 1,kend(n) + nh - 1
            vol(j,k,n) = vol(j1,k,n)
          enddo
 
          k = kbeg(n) - nh
          k1 = kbeg(n) + nh - 1
          do j = jbeg(n) - nh,jend(n) + nh
            vol(j,k,n) = vol(j,k1,n)
          enddo

          k = kend(n) + nh
          k1 = kend(n) - nh + 1
          do j = jbeg(n) - nh,jend(n) + nh
            vol(j,k,n) = vol(j,k1,n)
          enddo
        enddo
c
c..check for negative jacobians
c
        nneg = 0
        do k = 1,kd(n)
          do  j = 1,jd(n)
            if( vol(j,k,n).le.0.0 ) then
              nneg = nneg+1
            end if
          enddo
        enddo
c
        if(nneg .ne. 0) then
          write(6,*) nneg, ' negative jacobians in block'
          do k = 1,kd(n)
            do j = 1,jd(n)
              if( vol(j,k,n).le.0.0 ) then
                write(6,603) vol(j,k,n), j, k
              end if
            enddo
          enddo
        endif

      enddo

  603 format( ' ',10x,'negative jacobian = ',1p,e10.3,1x,'at j,k =',
     &                   2i5,5x)

      end subroutine metfv

c********************************************************************
      subroutine write_connect_data
!...writing out the connectivity
        use global
        implicit none

c********************************************************************

        integer j,k,n,nf,nd,logi
        character*40,integer_string,tmp

        do n=1,ngrids
          write(integer_string,*) n-1
          integer_string = adjustl(integer_string)

          !tmp = trim(grid_file(n))//'_inter'!.'//trim(integer_string)
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
      subroutine write_grid
!...writing out the connectivity
        use global
        implicit none

c********************************************************************

        integer j,k,n,nf,nd,logi
        character*40,integer_string,tmp

        do n=1,ngrids
          write(integer_string,*) n-1
          integer_string = adjustl(integer_string)

          !tmp = trim(grid_file(n))//'_inter'!.'//trim(integer_string)
          tmp = 'g.'//trim(integer_string)//'_new'
          open(unit=1,file=tmp,form='unformatted',status='unknown')
          
          logi = 1
          write(logi) jend(n)-jbeg(n)+1,kend(n)-kbeg(n)+1
          write(logi) ((xc(j,k,n),j=jbeg(n),jend(n)),k=kbeg(n),kend(n)),
     &                ((yc(j,k,n),j=jbeg(n),jend(n)),k=kbeg(n),kend(n)),
     &             ((iblank(j,k,n),j=jbeg(n),jend(n)),k=kbeg(n),kend(n))
          close(logi)
        enddo

      end subroutine write_grid

c********************************************************************
