!********************************************************************
      module grid
        use input
        implicit none

        integer ngrids
        integer,allocatable :: jmax(:),kmax(:)
        integer,allocatable :: jd(:),kd(:)
        integer,allocatable :: jbeg(:),kbeg(:)
        integer,allocatable :: jend(:),kend(:)
        integer :: jm,km,jmc,kmc
        real,allocatable :: x(:,:,:),y(:,:,:)
        real,allocatable :: xc(:,:,:),yc(:,:,:)
        real,allocatable :: vol(:,:,:)
        character*40,allocatable :: grid_file(:)
        
        logical :: gridinit

        data ngrids /1/
        data gridinit /.false./
        
      contains
        
!********************************************************************
      subroutine init_grid()
!..reads grid and computes cell-center and volume        
!********************************************************************

        integer j,k,n
        integer logg

        gridinit = .true.

        if (interpolate_space) ngrids = 2

        if(allocated(grid_file)) deallocate(grid_file)
        allocate(grid_file(ngrids))

        grid_file(1) = 'fort.1'
        if (interpolate_space) grid_file(2) = 'finefort.1'

!....determining the dimensions for grid

        if(allocated(jmax)) deallocate(jmax)
        if(allocated(kmax)) deallocate(kmax)
        if(allocated(jd)) deallocate(jd)
        if(allocated(kd)) deallocate(kd)
        if(allocated(jbeg)) deallocate(jbeg)
        if(allocated(jend)) deallocate(jend)
        if(allocated(kbeg)) deallocate(kbeg)
        if(allocated(kend)) deallocate(kend)

        allocate(jmax(ngrids),kmax(ngrids))
        allocate(jd(ngrids),kd(ngrids))
        allocate(jbeg(ngrids),jend(ngrids))
        allocate(kbeg(ngrids),kend(ngrids))

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
 
!....reading the grid
        if (interpolate_space) then ! grid files not needed if no space interpolation

          if(allocated(x)) deallocate(x)
          if(allocated(y)) deallocate(y)
          if(allocated(xc)) deallocate(xc)
          if(allocated(yc)) deallocate(xc)

          allocate(x(jm,km,ngrids),y(jm,km,ngrids))
          allocate(xc(jmc,kmc,ngrids),yc(jmc,kmc,ngrids))

          x = 0.; y = 0.; xc = 0.; yc = 0.
          do n = 1,ngrids
            open(unit=1,file=grid_file(n),form='unformatted',
     >                                      status='unknown')
            logg=1
            read(logg) jmax(n),kmax(n)
            read(logg) ((x(j,k,n),j=1,jmax(n)),k=1,kmax(n)),
     >                 ((y(j,k,n),j=1,jmax(n)),k=1,kmax(n))
            close(logg)

          enddo
          call find_cellcenter()
        
!...compute cell volume

          if(allocated(vol)) deallocate(vol)

          allocate(vol(jmc,kmc,ngrids))
          call metfv()
        endif

      end subroutine init_grid

!*********************************************************************
      subroutine find_cellcenter()
!.. find cell center
!*********************************************************************
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

!***********************************************************************
      subroutine metfv()
!***********************************************************************
      integer j,k,n,jc,kc,j1,k1,nh,nneg
      real dx1,dy1,dx2,dy2

      do n = 1,ngrids
!..find the volume near j,k

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

!..boundary values
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
!..check for negative jacobians
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

!********************************************************************
      subroutine write_grid()
!********************************************************************
        integer j,k,n
        character*40,gfile

        gfile = 'finefort.9'
        gfile = trim(adjustl(gfile))

        open(unit=1,file=gfile,form='unformatted',
     >                                    status='unknown')
        n = 2
        write(1) jd(n),kd(n)
        write(1) ((xc(j,k,n),j=1,jd(n)),k=1,kd(n)),
     >           ((yc(j,k,n),j=1,jd(n)),k=1,kd(n))
        close(1)

      end subroutine write_grid

      end module grid

!********************************************************************
