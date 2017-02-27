!********************************************************************
module ihc
  implicit none
!********************************************************************

  integer,private :: ngrids

  type logic2D
     logical, allocatable :: arr(:,:)
  end type logic2D

  type integer2D
     integer, allocatable :: arr(:,:)
  end type integer2D

  type real2D
     real, allocatable :: arr(:,:)
  end type real2D

  type(integer2D),allocatable,private :: cell_bc(:),ichk(:)
  type(real2D),allocatable,private    :: vold(:),volr(:)
  type(logic2D),allocatable :: wallpoint(:),forcerecv(:),noblank(:)

  integer,allocatable,private :: jmax(:),kmax(:)
  integer,allocatable,private :: jmaxc(:),kmaxc(:)
  real,pointer,private    :: x(:,:,:),y(:,:,:)
  real,pointer,private    :: xc(:,:,:),yc(:,:,:)
  integer,pointer,private :: iblank(:,:,:)

  character*40,allocatable,private :: bc_file(:)

  integer,allocatable,private :: check_alloc(:)

  real,private :: overlap_factor

  private :: init_ihc,read_bc,boundconnect,recvconnect,holeconnect, &
             finddonor,find_fractions,bilinearinterp,matrixinv2x2

contains

!********************************************************************
  subroutine init_ihc(jmx,kmx,ngrds,bcfile)
!********************************************************************
    integer :: ngrds
    integer :: jmx(:),kmx(:)
    character*40 :: bcfile(:)

!...local variables
    integer :: n

    ngrids = ngrds

    if(allocated(cell_bc)) then
      do n = 1,size(cell_bc)
        if(allocated(cell_bc(n)%arr)) deallocate(cell_bc(n)%arr)
      enddo
      deallocate(cell_bc)
    endif

    if(allocated(ichk)) then 
      do n = 1,size(ichk)
        if(allocated(ichk(n)%arr)) deallocate(ichk(n)%arr)
      enddo
      deallocate(ichk)
    endif

    if(allocated(vold)) then 
      do n = 1,size(vold)
        if(allocated(vold(n)%arr)) deallocate(vold(n)%arr)
      enddo
      deallocate(vold)
    endif

    if(allocated(volr)) then 
      do n = 1,size(volr)
        if(allocated(volr(n)%arr)) deallocate(volr(n)%arr)
      enddo
      deallocate(volr)
    endif

    if(allocated(wallpoint)) then 
      do n = 1,size(wallpoint)
        if(allocated(wallpoint(n)%arr)) deallocate(wallpoint(n)%arr)
      enddo
      deallocate(wallpoint)
    endif

    if(allocated(forcerecv)) then 
      do n = 1,size(forcerecv)
        if(allocated(forcerecv(n)%arr)) deallocate(forcerecv(n)%arr)
      enddo
      deallocate(forcerecv)
    endif

    if(allocated(noblank)) then 
      do n = 1,size(noblank)
        if(allocated(noblank(n)%arr)) deallocate(noblank(n)%arr)
      enddo
      deallocate(noblank)
    endif

    if(allocated(check_alloc)) deallocate(check_alloc)
    if(allocated(jmax)) deallocate(jmax)
    if(allocated(kmax)) deallocate(kmax)
    if(allocated(jmaxc)) deallocate(jmaxc)
    if(allocated(kmaxc)) deallocate(kmaxc)
    if(allocated(bc_file)) deallocate(bc_file)

    nullify(x,y,xc,yc,iblank)

    allocate(jmax(ngrids),kmax(ngrids))
    allocate(jmaxc(ngrids),kmaxc(ngrids))
    allocate(bc_file(ngrids))

    jmax = jmx
    kmax = kmx
    jmaxc = jmax - 1
    kmaxc = kmax - 1
    bc_file = bcfile

    allocate(vold(ngrids),volr(ngrids))
    allocate(cell_bc(ngrids),ichk(ngrids))
    allocate(wallpoint(ngrids),forcerecv(ngrids),noblank(ngrids))

    do n = 1,ngrids
      allocate(vold(n)%arr(0:jmax(n),0:kmax(n)))
      allocate(volr(n)%arr(0:jmax(n),0:kmax(n)))
      allocate(cell_bc(n)%arr(0:jmax(n),0:kmax(n)))
      allocate(wallpoint(n)%arr(jmax(n),kmax(n)))
      allocate(forcerecv(n)%arr(0:jmax(n),0:kmax(n)))
      allocate(noblank(n)%arr(0:jmax(n),0:kmax(n)))
      allocate(ichk(n)%arr(0:jmax(n),0:kmax(n)))
    enddo

    overlap_factor = 0.95

    !allocate(x(jm,km,ngrids),y(jm,km,ngrids))
    !allocate(iblank(jm,km,ngrids))

    allocate(check_alloc(1))

    call read_bc
 
  end subroutine init_ihc

!********************************************************************
  subroutine read_bc
!..bc types
!...1 --> not a receiver 
!...2 --> iblank = 0  
!...3 --> is a receiver 
!...4 --> iblank.ne.0  
!...5 --> wall point 
!...6 --> stronger immunization, iblank = 1 (is a combination of bc's 1 and 4)
!...7 --> combination of 3 and 4 bc's
!********************************************************************
    integer :: j,k,n,ib
    integer :: js,je,ks,ke
    integer :: nbc,ibtyp(25),ibdir(25)
    integer :: jbcs(25),kbcs(25)
    integer :: jbce(25),kbce(25)
    integer :: ibproc(25)

    namelist/ihcbcinp/ nbc,ibtyp,ibdir,jbcs,jbce,kbcs,kbce,ibproc

    do n = 1,ngrids
      cell_bc(n)%arr = 0
      wallpoint(n)%arr = .false.
      forcerecv(n)%arr = .false.
      noblank(n)%arr = .false.
    enddo

    do n = 1,ngrids
      open(unit=21,file=bc_file(n),status='unknown')
      read(21,ihcbcinp)
      close(21)

      do ib = 1,nbc
        js = jbcs(ib)
        je = jbce(ib)
        ks = kbcs(ib)
        ke = kbce(ib)

        if (ibtyp(ib).eq.5) then
          if(js.lt.0) js = jmax(n) + js + 1
          if(ks.lt.0) ks = kmax(n) + ks + 1
          if(je.lt.0) je = jmax(n) + je + 1
          if(ke.lt.0) ke = kmax(n) + ke + 1
          do j = js,je
          do k = ks,ke
            wallpoint(n)%arr(j,k) = .true.
          enddo
          enddo
        else
          if(js.lt.0) js = jmaxc(n) + js + 1
          if(ks.lt.0) ks = kmaxc(n) + ks + 1
          if(je.lt.0) je = jmaxc(n) + je + 1
          if(ke.lt.0) ke = kmaxc(n) + ke + 1
          do j = js,je
          do k = ks,ke
            cell_bc(n)%arr(j,k) = ibtyp(ib)
          enddo
          enddo
        endif
      enddo 
    enddo

  end subroutine read_bc

!*********************************************************************
  subroutine do_connectihc(xy,xyc,vol,iblnk,jmx,kmx,nhalo,imesh,idonor,frac,&
              ibc,ndonor,nfringe,iisptr,iieptr,ngrds,init)

    integer,optional :: init
    integer :: ngrds
    integer :: jmx(:),kmx(:),nhalo
    integer :: ndonor(:),nfringe(:)
    real,target :: xy(:,:,:,:),xyc(:,:,:,:)
    real :: vol(:,:,:)
    integer,target :: iblnk(:,:,:)
    integer :: imesh(:,:,:),idonor(:,:,:)
    integer :: iisptr(:),iieptr(:)
    real    :: frac(:,:,:) 
    integer :: ibc(:,:)
    
    integer :: n,jbeg,kbeg
    character*40 :: integer_string

    real,pointer :: x1(:,:,:),y1(:,:,:)
    real,pointer :: xc1(:,:,:),yc1(:,:,:)
    character*40,allocatable :: bcfile(:)

    x1 => xy(1,:,:,:)
    y1 => xy(2,:,:,:)
    xc1 => xyc(1,:,:,:)
    yc1 => xyc(2,:,:,:)

    allocate(bcfile(ngrds))

    do n=1,ngrds
      write(integer_string,*) n-1
      integer_string = adjustl(integer_string)
      bcfile(n) = 'ihcbc.'//trim(integer_string)
    enddo

    call connect2d(x1,y1,xc1,yc1,vol,iblnk,idonor,frac,imesh,ibc,&
              iisptr,iieptr,ndonor,nfringe,jmx,kmx,nhalo,ngrds,bcfile,init)

  end subroutine do_connectihc

!*********************************************************************
  subroutine connect2d(x1,y1,xc1,yc1,vol,iblnk,idonor,frac,imesh,ibc,&
              iisptr,iieptr,ndonor,nfringe,jmx,kmx,nhalo,ngrds,bcfile,init)
    use octree
!********************************************************************

    integer,optional :: init
    integer :: ngrds,nhalo
    integer :: jmx(:),kmx(:)
    integer :: ndonor(:),nfringe(:)
    integer :: iisptr(:),iieptr(:)
    real,target :: x1(:,:,:),y1(:,:,:)
    real,target :: xc1(-nhalo+1:,-nhalo+1:,:),yc1(-nhalo+1:,-nhalo+1:,:)
    real,target :: vol(-nhalo+1:,-nhalo+1:,:)
    integer,target :: iblnk(-nhalo+1:,-nhalo+1:,:)
    integer :: imesh(:,:,:),idonor(:,:,:)
    real    :: frac(:,:,:) 
    integer :: ibc(:,:)
    character*40 :: bcfile(:)

    integer jk
    integer j,k,j1,k1,m,n,nn,n1
    integer idsize,count,llmin,llmax
    integer,allocatable :: nholept(:),iholept(:,:,:)
    integer,allocatable :: ibctmp(:,:),idontmp(:,:)
    real,allocatable :: diholept(:,:,:)
    logical,allocatable :: ioverlap(:,:)

!...allocate variables under user request

    if(present(init)) then
      if (init.eq.1) call init_ihc(jmx,kmx,ngrds,bcfile)
    endif

!...allocate variables if not initiated yet

    if(.not.allocated(check_alloc)) call init_ihc(jmx,kmx,ngrds,bcfile)
     
!...set pointers

    x => x1   ! x,y,xc,yc set as pointers to avoid additional memory usage
    y => y1   ! have to be careful not to change x,y,xc,yc in any of the
    xc => xc1 ! subroutines associated with this module 
    yc => yc1

    iblank => iblnk ! iblank values do change in this module

!...some allocatation

    idsize = 0
    do n = 1,ngrids
      jk = jmax(n)*kmax(n)
      if(jk.gt.idsize) idsize = jk 
    enddo

    allocate(nholept(ngrids))
    allocate(iholept(idsize,5,ngrids))
    allocate(diholept(idsize,2,ngrids))
    allocate(ibctmp(idsize,ngrids),idontmp(idsize,ngrids))

    allocate(ioverlap(ngrids,ngrids))

!...assign volume of the cells

    do n = 1,ngrids
      do j = 0,jmax(n)
      do k = 0,kmax(n)
        vold(n)%arr(j,k) = vol(j,k,n) 
        volr(n)%arr(j,k) = vol(j,k,n) 
      enddo
      enddo
    enddo

!...reassign cell volume based on the bc

    do n = 1,ngrids
      do j = 1,jmaxc(n)
      do k = 1,kmaxc(n)
        if (cell_bc(n)%arr(j,k).eq.1) volr(n)%arr(j,k) = 0.
        if (cell_bc(n)%arr(j,k).eq.1) vold(n)%arr(j,k) = 0.
        if (cell_bc(n)%arr(j,k).eq.2) volr(n)%arr(j,k) = -1.e30
        if (cell_bc(n)%arr(j,k).eq.2) vold(n)%arr(j,k) = 1.e30
        if (cell_bc(n)%arr(j,k).eq.3) volr(n)%arr(j,k) = 1.e30
        if (cell_bc(n)%arr(j,k).eq.3) vold(n)%arr(j,k) = 1.e30
        if (cell_bc(n)%arr(j,k).eq.6) volr(n)%arr(j,k) = 0.
        if (cell_bc(n)%arr(j,k).eq.6) vold(n)%arr(j,k) = 0.
        if (cell_bc(n)%arr(j,k).eq.7) volr(n)%arr(j,k) = 1.e30
        if (cell_bc(n)%arr(j,k).eq.7) vold(n)%arr(j,k) = 1.e30

        if (cell_bc(n)%arr(j,k).eq.3.or.cell_bc(n)%arr(j,k).eq.7) then
          forcerecv(n)%arr(j,k) = .true.
        endif

        if (cell_bc(n)%arr(j,k).eq.4.or.cell_bc(n)%arr(j,k).eq.6.or. &
            cell_bc(n)%arr(j,k).eq.7) then
          noblank(n)%arr(j,k) = .true.
        endif

      enddo
      enddo
    enddo

!...generate octree data
    llmin = 2
    llmax = 5
    call make_octree(x,y,jmax,kmax,ngrids,llmin,llmax)

!...initialize iblank array

    iblank = 1
    do n = 1,ngrids
      do j = 1,jmaxc(n)
      do k = 1,kmaxc(n)
        if (cell_bc(n)%arr(j,k).eq.2) iblank(j,k,n) = 0
      enddo
      enddo
    enddo

    do n = 1,ngrids
      ichk(n)%arr = 0
    enddo

!...check for overlapping meshes 

    ioverlap = .false.

    do n = 1,ngrids
      call boundconnect(n,ioverlap)
    enddo

!...do connectivity for forced receivers first

    nholept = 0

    do n = 1,ngrids
      call recvconnect(n,ioverlap,nholept,iholept,diholept,idsize)
    enddo

!...do connectivity for all other points

    do n = 1,ngrids
      call holeconnect(n,ioverlap,nholept,iholept,diholept,idsize)
    enddo

!... get connectivity info....
    ndonor = 0
    nfringe = 0

    do n=1,ngrids
     do m=1,nholept(n)
       n1 = iholept(m,3,n)

       j = iholept(m,1,n)
       k = iholept(m,2,n)
       if (iblank(j,k,n).ne.0) then
         j1 = iholept(m,4,n)
         k1 = iholept(m,5,n)
         nfringe(n) = nfringe(n) + 1
         imesh(nfringe(n),1,n) = j + nhalo
         imesh(nfringe(n),2,n) = k + nhalo

         ndonor(n1) = ndonor(n1) + 1
         ibctmp(nfringe(n),n) = ndonor(n1)
         idontmp(nfringe(n),n) = n1

         idonor(ndonor(n1),1,n1) = j1 + nhalo
         idonor(ndonor(n1),2,n1) = k1 + nhalo
         frac(ndonor(n1),1,n1) = diholept(m,1,n)
         frac(ndonor(n1),2,n1) = diholept(m,2,n)
         iblank(j,k,n) = -1
       endif
     enddo
    enddo
    
    !global donor cell ptrs
    iisptr = 0
    iieptr = 0
    iisptr(1) = 1
    iieptr(1) = ndonor(1)
    do n=2,ngrids
      iisptr(n) = iisptr(n-1) + ndonor(n-1) 
      iieptr(n) = iisptr(n) + ndonor(n) -1 
    end do 

    do n=1,ngrids
     do j=1,nfringe(n)
       n1 = idontmp(j,n)
       count = 0
       if(n1.gt.1) then
        do nn=2,n1
         count = count+ndonor(nn-1)
        end do
       end if
       ibc(j,n) = count+ibctmp(j,n)
     end do
    end do

  end subroutine connect2d

!********************************************************************
  subroutine boundconnect(ng,ioverlap)
    use octree
!********************************************************************

    integer :: ng
    logical :: ioverlap(ngrids,ngrids)

    integer :: j,k,n,nface
    integer :: lvl,nlvl,lvlcheck,nout
    integer :: jl,kl,js,ks,jstep,kstep

    real    :: xbar(2),xp(2)
    real,allocatable :: xg(:,:,:)

    integer :: nstop,ntime
    integer,allocatable :: jst(:),kst(:)
    integer,allocatable :: jet(:),ket(:)
    integer,allocatable :: box(:)

    allocate(jst(2),kst(2))
    allocate(jet(2),ket(2))
    allocate(box(2))

    do n = 1,ngrids

    if (ioverlap(n,ng)) ioverlap(ng,n) = .true.

    if (n.ne.ng.and..not.ioverlap(ng,n)) then
      allocate(xg(jmax(n),kmax(n),2))

      do j = 1,jmax(n)
        do k = 1,kmax(n)
          xg(j,k,1) = x(j,k,n)
          xg(j,k,2) = y(j,k,n)
        enddo
      enddo

      jl = 2; kl = 2
      js = -1; ks = -1

      j = 1; k = 1
      jstep = 1; kstep = 1
      
      nface = 1
      do while (.true.) !loop for all boundary points 

        if (nface.eq.1) j = 1
        if (nface.eq.2) j = jmaxc(ng)
        if (nface.eq.3) k = 1
        if (nface.eq.4) k = kmaxc(ng)

        xp(1) = xc(j,k,ng)
        xp(2) = yc(j,k,ng)

        xbar(1) = (xp(1)-xs(1,n))/scale(1,n)
        xbar(2) = (xp(2)-xs(2,n))/scale(2,n)

        if (xbar(1).ge.0.and.xbar(2).ge.0.and. &
            xbar(1).lt.1.and.xbar(2).lt.1) then

          call boxindex(xbar,lvlmin,nout,2)
          if (empty_flag(nout,n)) go to 10

          js = jl
          ks = kl
          call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)

          if (js.le.0) then

            call boxindex(xbar,lvlmax,nout,2)
            lvlcheck = lvlpoint(nout,n)
            nlvl = 2
            if (lvlcheck.eq.lvlmin) nlvl = 1

            do lvl = 1,nlvl 
              call boxindex(xbar,lvlcheck-lvl+1,nout,2)
              box(lvl) = nout
              jst(lvl) = istart(nout,1,n)
              kst(lvl) = istart(nout,2,n)
              jet(lvl) = iend(nout,1,n)
              ket(lvl) = iend(nout,2,n)
            enddo

            nstop = 0
            ntime = 0

            do while(nstop.lt.2*nlvl.and.ntime.lt.1)

              do lvl = 1,nlvl 

                nout = box(lvl)

                if (js.le.0.and.kst(lvl).gt.0) then
                  js = abs(jst(lvl))
                  ks = abs(kst(lvl))
                  call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)
                  jst(lvl) = js
                  kst(lvl) = ks

                  if (js.le.0.and.ks.le.0) then
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

                if (js.le.0.and.ket(lvl).gt.0) then
                  js = abs(jet(lvl))
                  ks = abs(ket(lvl))
                  call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)
                  jet(lvl) = js
                  ket(lvl) = ks

                  if (js.le.0.and.ks.le.0) then
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

              enddo
              ntime = ntime + 1

            enddo
          endif

          if (js.gt.0) then
            jl = js
            kl = ks
            ioverlap(ng,n) = .true.
            ioverlap(n,ng) = .true.
            exit
          endif
        endif

 10     if (nface.eq.1.or.nface.eq.3) then
          k = k + kstep
          if (k.gt.kmaxc(ng).or.k.lt.1) then
            nface = nface + 1
            kstep = -1*kstep
            k = k + kstep
            cycle
          endif
        endif

        if (nface.eq.2.or.nface.eq.4) then
          j = j + jstep 
          if (j.gt.jmaxc(ng).or.j.lt.1) then
            nface = nface + 1
            jstep = -1*jstep
            j = j + jstep
            cycle
          endif
        endif

        if (nface.eq.5) exit
      enddo

      deallocate(xg)
    endif

    enddo

  end subroutine boundconnect

!********************************************************************
  subroutine recvconnect(ng,ioverlap,nholept,iholept,diholept,idsize)
    use octree
!********************************************************************

    integer :: ng,idsize
    integer :: nholept(ngrids)
    integer :: iholept(idsize,5,ngrids)
    real    :: diholept(idsize,2,ngrids)
    logical :: ioverlap(ngrids,ngrids)

    integer :: j,k,n,nc
    integer :: lvl,nlvl,lvlcheck,nout,ncount,nbodycross
    integer :: jl,kl,js,ks,jstep,kstep
    real    :: voldonor
    real    :: scalei1,scalei2
    integer :: nstop,ntime

    real    :: xbar(2),xp(2),frac1(2)
    integer,allocatable :: idonorptr(:,:)
    real,allocatable :: xg(:,:,:),xgc(:,:,:)
    real,allocatable :: volrecv(:,:)
 
    integer,allocatable :: jst(:),kst(:)
    integer,allocatable :: jet(:),ket(:)
    integer,allocatable :: box(:)

    allocate(jst(2),kst(2))
    allocate(jet(2),ket(2))
    allocate(box(2))

    allocate(volrecv(jmaxc(ng),kmaxc(ng)))
    allocate(idonorptr(jmaxc(ng),kmaxc(ng)))

    do j = 1,jmaxc(ng)
      do k = 1,kmaxc(ng)
        volrecv(j,k) = overlap_factor*volr(ng)%arr(j,k)
      enddo
    enddo

    idonorptr = 0

    ncount = nholept(ng)

    do n = 1,ngrids
    if (n.ne.ng.and.ioverlap(ng,n)) then
      allocate(xg(jmax(n),kmax(n),2))
      allocate(xgc(0:jmax(n),0:kmax(n),2))

      do j = 1,jmax(n)
        do k = 1,kmax(n)
          xg(j,k,1) = x(j,k,n)
          xg(j,k,2) = y(j,k,n)
        enddo
      enddo

      do j = 0,jmax(n)
        do k = 0,kmax(n)
          xgc(j,k,1) = xc(j,k,n)
          xgc(j,k,2) = yc(j,k,n)
        enddo
      enddo

      jl = 2; kl = 2
      js = -1; ks = -1

      j = 1; k = 1
      jstep = 1; kstep = 1
      
      scalei1 = 1./scale(1,n)
      scalei2 = 1./scale(2,n)

      do while (.true.) !loop for all points 

        nbodycross = 0

        if (.not.forcerecv(ng)%arr(j,k)) go to 10

        xp(1) = xc(j,k,ng)
        xp(2) = yc(j,k,ng)

        xbar(1) = (xp(1)-xs(1,n))*scalei1
        xbar(2) = (xp(2)-xs(2,n))*scalei2

        if (volrecv(j,k).ge.0..and. &
            xbar(1).ge.0.and.xbar(2).ge.0.and. &
            xbar(1).lt.1.and.xbar(2).lt.1) then

          call boxindex(xbar,lvlmin,nout,2)
          if (empty_flag(nout,n)) go to 10

          js = jl
          ks = kl
          call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)

          if (js.le.0) then
            if (ks.le.0.and.wallpoint(n)%arr(abs(js),abs(ks))) then
              nbodycross = nbodycross + 1
            endif

            call boxindex(xbar,lvlmax,nout,2)
            lvlcheck = lvlpoint(nout,n)
            nlvl = 2
            if (lvlcheck.eq.lvlmin) nlvl = 1

            do lvl = 1,nlvl
              call boxindex(xbar,lvlcheck-lvl+1,nout,2)
              box(lvl) = nout
              jst(lvl) = istart(nout,1,n)
              kst(lvl) = istart(nout,2,n)
              jet(lvl) = iend(nout,1,n)
              ket(lvl) = iend(nout,2,n)
            enddo

            nstop = 0
            ntime = 0

            do while(nstop.lt.2*nlvl.and.ntime.lt.5)

              do lvl = 1,nlvl

                nout = box(lvl)

                if (js.le.0.and.kst(lvl).gt.0) then
                  js = abs(jst(lvl))
                  ks = abs(kst(lvl))
                  call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)
                  jst(lvl) = js
                  kst(lvl) = ks

                  if (js.le.0.and.ks.le.0) then
                    if (wallpoint(n)%arr(abs(js),abs(ks))) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

                if (js.le.0.and.ket(lvl).gt.0) then
                  js = abs(jet(lvl))
                  ks = abs(ket(lvl))
                  call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)
                  jet(lvl) = js
                  ket(lvl) = ks

                  if (js.le.0.and.ks.le.0) then
                    if (wallpoint(n)%arr(abs(js),abs(ks))) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

              enddo
              ntime = ntime + 1

            enddo 
          endif

          if (js.gt.0) then

            voldonor = vold(n)%arr(js,ks)
            if (forcerecv(n)%arr(js-1,ks-1) .or. &
                forcerecv(n)%arr(js  ,ks-1) .or. &
                forcerecv(n)%arr(js+1,ks-1) .or. &
                forcerecv(n)%arr(js-1,ks  ) .or. &
                forcerecv(n)%arr(js  ,ks  ) .or. &
                forcerecv(n)%arr(js+1,ks  ) .or. &
                forcerecv(n)%arr(js-1,ks+1) .or. &
                forcerecv(n)%arr(js  ,ks+1) .or. &
                forcerecv(n)%arr(js+1,ks+1)) voldonor = 1.e30

            if (voldonor.lt.volrecv(j,k)) then
              volrecv(j,k) = voldonor
              if (idonorptr(j,k).eq.0) then
                ncount = ncount + 1
                nc = ncount
                nholept(ng) = ncount
                idonorptr(j,k) = ncount
              else
                nc = idonorptr(j,k)
              endif
              iholept(nc,1,ng) = j
              iholept(nc,2,ng) = k
              iholept(nc,3,ng) = n
              iholept(nc,4,ng) = js
              iholept(nc,5,ng) = ks
              vold(ng)%arr(j-1, k-1) = 1.e30
              vold(ng)%arr(j  , k-1) = 1.e30
              vold(ng)%arr(j+1, k-1) = 1.e30
              vold(ng)%arr(j-1, k  ) = 1.e30
              vold(ng)%arr(j  , k  ) = 1.e30
              vold(ng)%arr(j+1, k  ) = 1.e30
              vold(ng)%arr(j-1, k+1) = 1.e30
              vold(ng)%arr(j  , k+1) = 1.e30
              vold(ng)%arr(j+1, k+1) = 1.e30
              call find_fractions(xp,xgc,frac1,js,ks, &
                                          jmax(n),kmax(n))
              diholept(nc,1,ng) = frac1(1)
              diholept(nc,2,ng) = frac1(2)

            endif

            jl = js
            kl = ks
          endif

          if (js.le.0.and.nbodycross.gt.2*nlvl) then
            if (.not.noblank(ng)%arr(j,k)) then
              iblank(j,k,ng) = 0
              volrecv(j,k) = -1.e30
            endif
          endif

          ichk(ng)%arr(j,k) = 1
        endif

 10     k = k + kstep
        if (k.gt.kmaxc(ng).or.k.lt.1) then
          kstep = -1*kstep
          k = k + kstep
          j = j + jstep
        endif
        if (j.gt.jmaxc(ng)) exit
      enddo

      deallocate(xg,xgc)
    endif
    enddo

    do nc = 1,ncount
      n  = iholept(nc,3,ng)
      js = iholept(nc,4,ng)
      ks = iholept(nc,5,ng)
      volr(n)%arr(js-1, ks-1) = 0.
      volr(n)%arr(js  , ks-1) = 0.
      volr(n)%arr(js+1, ks-1) = 0.
      volr(n)%arr(js-1, ks  ) = 0.
      volr(n)%arr(js  , ks  ) = 0.
      volr(n)%arr(js+1, ks  ) = 0.
      volr(n)%arr(js-1, ks+1) = 0.
      volr(n)%arr(js  , ks+1) = 0.
      volr(n)%arr(js+1, ks+1) = 0.
    enddo

  end subroutine recvconnect

!********************************************************************
  subroutine holeconnect(ng,ioverlap,nholept,iholept,diholept,idsize)
    use octree
!********************************************************************

    integer :: ng,idsize
    integer :: nholept(ngrids)
    integer :: iholept(idsize,5,ngrids)
    real    :: diholept(idsize,2,ngrids)
    logical :: ioverlap(ngrids,ngrids)

    integer :: j,k,n,nc
    integer :: lvl,nlvl,lvlcheck,nout,ncount,nbodycross
    integer :: jl,kl,js,ks,jstep,kstep
    real    :: voldonor
    real    :: scalei1,scalei2

    real    :: xbar(2),xp(2),frac1(2)
    integer,allocatable :: idonorptr(:,:)
    real,allocatable :: xg(:,:,:),xgc(:,:,:)
    real,allocatable :: volrecv(:,:)

    integer :: nstop,ntime
    integer,allocatable :: jst(:),kst(:)
    integer,allocatable :: jet(:),ket(:)
    integer,allocatable :: box(:)

    allocate(jst(2),kst(2))
    allocate(jet(2),ket(2))
    allocate(box(2))

    allocate(volrecv(jmaxc(ng),kmaxc(ng)))
    allocate(idonorptr(jmax(ng),kmaxc(ng)))

    do j = 1,jmaxc(ng)
      do k = 1,kmaxc(ng)
        volrecv(j,k) = overlap_factor*volr(ng)%arr(j,k)
      enddo
    enddo

    idonorptr = 0

    ncount = nholept(ng)

    do n = 1,ngrids
    if (n.ne.ng.and.ioverlap(ng,n)) then
      allocate(xg(jmax(n),kmax(n),2))
      allocate(xgc(0:jmax(n),0:kmax(n),2))

      do j = 1,jmax(n)
        do k = 1,kmax(n)
          xg(j,k,1) = x(j,k,n)
          xg(j,k,2) = y(j,k,n)
        enddo
      enddo

      do j = 0,jmax(n)
        do k = 0,kmax(n)
          xgc(j,k,1) = xc(j,k,n)
          xgc(j,k,2) = yc(j,k,n)
        enddo
      enddo

      jl = 2; kl = 2
      js = -1; ks = -1

      j = 1; k = 1
      jstep = 1; kstep = 1
      
      scalei1 = 1./scale(1,n)
      scalei2 = 1./scale(2,n)

      do while (.true.) !loop for all points 

        nbodycross = 0

        if (volrecv(j,k).lt.0..or.ichk(ng)%arr(j,k).eq.1) go to 10

        xp(1) = xc(j,k,ng)
        xp(2) = yc(j,k,ng)

        xbar(1) = (xp(1)-xs(1,n))*scalei1
        xbar(2) = (xp(2)-xs(2,n))*scalei2

        if (xbar(1).ge.0.and.xbar(2).ge.0.and. &
            xbar(1).lt.1.and.xbar(2).lt.1) then

          call boxindex(xbar,lvlmin,nout,2)
          if (empty_flag(nout,n)) go to 10

          js = jl
          ks = kl
          call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)

          if (js.le.0) then

            if (ks.le.0.and.wallpoint(n)%arr(abs(js),abs(ks))) then
              nbodycross = nbodycross + 1
            endif

            call boxindex(xbar,lvlmax,nout,2)
            lvlcheck = lvlpoint(nout,n)
            nlvl = 2
            if (lvlcheck.eq.lvlmin) nlvl = 1

            do lvl = 1,nlvl 
              call boxindex(xbar,lvlcheck-lvl+1,nout,2)
              box(lvl) = nout
              jst(lvl) = istart(nout,1,n)
              kst(lvl) = istart(nout,2,n)
              jet(lvl) = iend(nout,1,n)
              ket(lvl) = iend(nout,2,n)
            enddo

            nstop = 0
            ntime = 0

            do while(nstop.lt.2*nlvl.and.ntime.lt.5)

              do lvl = 1,nlvl 

                nout = box(lvl)

                if (js.le.0.and.kst(lvl).gt.0) then
                  js = abs(jst(lvl))
                  ks = abs(kst(lvl))
                  call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)
                  jst(lvl) = js
                  kst(lvl) = ks

                  if (js.le.0.and.ks.le.0) then
                    if (wallpoint(n)%arr(abs(js),abs(ks))) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

                if (js.le.0.and.ket(lvl).gt.0) then
                  js = abs(jet(lvl))
                  ks = abs(ket(lvl))
                  call finddonor(xp,xg,js,ks,jmax(n),kmax(n),100)
                  jet(lvl) = js
                  ket(lvl) = ks

                  if (js.le.0.and.ks.le.0) then
                    if (wallpoint(n)%arr(abs(js),abs(ks))) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

              enddo
              ntime = ntime + 1

            enddo
          endif

          if (js.gt.0) then

            voldonor = vold(n)%arr(js,ks)
            if (forcerecv(n)%arr(js-1,ks-1) .or. &
                forcerecv(n)%arr(js  ,ks-1) .or. &
                forcerecv(n)%arr(js+1,ks-1) .or. &
                forcerecv(n)%arr(js-1,ks  ) .or. &
                forcerecv(n)%arr(js  ,ks  ) .or. &
                forcerecv(n)%arr(js+1,ks  ) .or. &
                forcerecv(n)%arr(js-1,ks+1) .or. &
                forcerecv(n)%arr(js  ,ks+1) .or. &
                forcerecv(n)%arr(js+1,ks+1)) voldonor = 1.e30

            if (voldonor.lt.volrecv(j,k)) then
              volrecv(j,k) = voldonor
              if (idonorptr(j,k).eq.0) then
                ncount = ncount + 1
                nc = ncount
                nholept(ng) = ncount
                idonorptr(j,k) = ncount
              else
                nc = idonorptr(j,k)
              endif
              iholept(nc,1,ng) = j
              iholept(nc,2,ng) = k
              iholept(nc,3,ng) = n
              iholept(nc,4,ng) = js
              iholept(nc,5,ng) = ks
              vold(ng)%arr(j-1, k-1) = 1.e30
              vold(ng)%arr(j  , k-1) = 1.e30
              vold(ng)%arr(j+1, k-1) = 1.e30
              vold(ng)%arr(j-1, k  ) = 1.e30
              vold(ng)%arr(j  , k  ) = 1.e30
              vold(ng)%arr(j+1, k  ) = 1.e30
              vold(ng)%arr(j-1, k+1) = 1.e30
              vold(ng)%arr(j  , k+1) = 1.e30
              vold(ng)%arr(j+1, k+1) = 1.e30
              call find_fractions(xp,xgc,frac1,js,ks, &
                                      jmax(n),kmax(n))
              diholept(nc,1,ng) = frac1(1)
              diholept(nc,2,ng) = frac1(2)
            endif

            jl = js
            kl = ks
          endif

          if (js.le.0.and.nbodycross.gt.2*nlvl) then
            if (.not.noblank(ng)%arr(j,k)) then
              iblank(j,k,ng) = 0
              volrecv(j,k) = -1.e30
            endif
          endif
        endif

  10    k = k + kstep
        if (k.gt.kmaxc(ng).or.k.lt.1) then
          kstep = -1*kstep
          k = k + kstep
          j = j + jstep
        endif
        if (j.gt.jmaxc(ng)) exit
      enddo

      deallocate(xg,xgc)
    endif
    enddo

    print*,'Number of fringe points in mesh',ng,'is', ncount

    do nc = 1,ncount
      n  = iholept(nc,3,ng)
      js = iholept(nc,4,ng)
      ks = iholept(nc,5,ng)
      volr(n)%arr(js-1, ks-1) = 0.
      volr(n)%arr(js  , ks-1) = 0.
      volr(n)%arr(js+1, ks-1) = 0.
      volr(n)%arr(js-1, ks  ) = 0.
      volr(n)%arr(js  , ks  ) = 0.
      volr(n)%arr(js+1, ks  ) = 0.
      volr(n)%arr(js-1, ks+1) = 0.
      volr(n)%arr(js  , ks+1) = 0.
      volr(n)%arr(js+1, ks+1) = 0.
    enddo

  end subroutine holeconnect

!********************************************************************
      subroutine finddonor(xp,x,js,ks,jmax,kmax,maxsearch)
!*********************************************************************
      integer :: js,ks,jmax,kmax,maxsearch
      real :: x(jmax,kmax,2)
      real :: xp(2)

      integer :: jj,kk,j,k,m,n,nsearch
      real    :: dum
      logical :: Inside,Outside,alter,notoutside
      integer :: movej,movek
      integer :: movejp,movekp
      integer,allocatable :: i_(:)
      real,allocatable :: xc(:,:)
      real,allocatable :: p(:),q(:),r(:),pqr(:)
      integer :: mo(4),ma(4),mb(4)

      allocate(i_(2))
      allocate(xc(4,2))
      allocate(p(2),q(2),r(2),pqr(4)) ! 2**Ndim

      DATA mo(1),mo(2),mo(3),mo(4) /1,1,2,3/
      DATA ma(1),ma(2),ma(3),ma(4) /3,2,4,4/
      DATA mb(1),mb(2),mb(3),mb(4) /2,3,1,1/

!..Get starting cell index (previous cell) j,k for the search
        j = min(js,jmax-1)
        k = min(ks,kmax-1)

        Inside = .FALSE.                    !assumed false to enter the loop
        Outside = .FALSE.                   !assumed false initially
        notoutside = .FALSE.                !assumed false for maxsearch
        nsearch = 0
        if (maxsearch.eq.0) maxsearch = 60

        movejp = 0; movekp = 0
        do while (.not.Inside.and..not.Outside)
          Inside = .TRUE.                   !assumed true intially
          jj = j
          kk = k

          do n = 1,2
            xc(1,n) = x(jj,  kk  ,n)
            xc(2,n) = x(jj+1,kk  ,n)
            xc(3,n) = x(jj,  kk+1,n)
            xc(4,n) = x(jj+1,kk+1,n)
          enddo

          movej  = 0; movek  = 0
!..Cross+dot product pqr=(OPxOQ).OR for in/outside cell test of point xp
          DO m=1,4                       !2*Ndim = number of cell faces
            DO n=1,2
              dum = xc(mo(m),n)
              p(n) = xc(ma(m),n) - dum
              q(n) = xc(mb(m),n) - dum
              r(n) = xp(n) - dum
            ENDDO
            pqr(m) = (p(1)*q(2)-p(2)*q(1))*(p(1)*r(2)-p(2)*r(1))
            IF(pqr(m) .lt. 0.) THEN  !If outside, get neighboring cell index
              Inside = .FALSE.
              IF(m .eq. 1) then
                j = j-1
                movej = movej-1
              endif
              IF(m .eq. 2) then
                k = k-1
                movek = movek-1
              endif
              IF(m .eq. 3) then
                j = j+1
                movej = movej+1
              endif
              IF(m .eq. 4) then
                k = k+1
                movek = movek+1
              endif
            ENDIF
          ENDDO

          alter = .false.
          if ((j.lt.1.or.j.ge.jmax).and.movek.ne.0) then
            j = j - movej
            movej = 0
            alter = .true.
          endif
          if ((k.lt.1.or.k.ge.kmax).and.movej.ne.0) then
            k = k - movek
            movek = 0
            alter = .true.
          endif

          if ((movejp+movej).eq.0.and.(movekp+movek).eq.0.and. &
              alter) then
            movej = 0
            movek = 0
          endif

          if (movej.eq.0.and.movek.eq.0.and..not.inside) &
            outside = .true.

          movejp = movej
          movekp = movek

          if (j.lt.1.or.k.lt.1.or.j.ge.jmax.or.k.ge.kmax) Outside = .TRUE.

          if(nsearch .ge. 3) then
            if(3*int((nsearch-1)/3) .eq. nsearch-1) then
              i_(1) = j
              i_(2) = k
            else
              if(i_(1).eq.j .and. i_(2).eq.k) &
                inside = .true.
            endif
          endif
          nsearch = nsearch + 1
          if (nsearch.gt.maxsearch) then
            outside = .true.
            notoutside = .true.
          endif

        enddo

        if (Inside.and..not.outside) then
          js = j
          ks = k
        else
          js = -abs(j-min(movej,0))
          ks = -abs(k-min(movek,0))
          if (notoutside) ks = -ks !make ks positive if stopped by maxsearch
        endif

      end subroutine finddonor

!********************************************************************
      subroutine find_fractions(xp,x,frac,jd,kd,jmax,kmax)
!..Interpolate coordinates frac in cell jd,kd by using Newton iteration
!********************************************************************
      integer :: jd,kd,jmax,kmax
      real :: x(0:jmax,0:kmax,2),frac(2)
      real :: xp(2)

      integer :: n,niter
      real    :: resid
      real,allocatable :: xv(:,:),B(:,:)

      real :: onefourth
      real djm1,dj0,djp1,dkm1,dk0,dkp1
      real ddjm1,ddj0,ddjp1,ddkm1,ddk0,ddkp1
      real w1,w2,w3,w4,w5,w6,w7,w8,w9
      real dw11,dw21,dw31,dw41,dw51,dw61,dw71,dw81,dw91
      real dw12,dw22,dw32,dw42,dw52,dw62,dw72,dw82,dw92

      allocate(xv(9,2),B(2,3))

      onefourth = 1./4

!..get coordinates of stencil points for interpolation

      do n = 1,2
        xv(1,n) = x(jd-1,kd-1,n)
        xv(2,n) = x(jd,  kd-1,n)
        xv(3,n) = x(jd+1,kd-1,n)
        xv(4,n) = x(jd-1,kd  ,n)
        xv(5,n) = x(jd,  kd  ,n)
        xv(6,n) = x(jd+1,kd  ,n)
        xv(7,n) = x(jd-1,kd+1,n)
        xv(8,n) = x(jd,  kd+1,n)
        xv(9,n) = x(jd+1,kd+1,n)
      enddo

!..initial frac guess (linear with no cross terms)

      do n=1,2
        b(n,1) = 0.5*(xv(6,n) - xv(4,n))
        b(n,2) = 0.5*(xv(8,n) - xv(2,n))
        b(n,3) = xp(n) - xv(5,n)
      enddo

      call matrixinv2x2(b)

      do n=1,2
        frac(n) = b(n,3)
        if((frac(n)-1)*(frac(n)+1) .gt. 0.) frac(n) = .5
      enddo

!..newton iteration with cross terms

      niter = 0
      resid = 1.e10

      do while (niter.lt.10.and.resid.gt.1e-15)

        djm1 = 1.+frac(1)
        dkm1 = 1.+frac(2)
        dj0  = frac(1)
        dk0  = frac(2)
        djp1 = 1.-frac(1)
        dkp1 = 1.-frac(2)

        ddjm1 =  1.
        ddkm1 =  1.
        ddj0  =  1.
        ddk0  =  1.
        ddjp1 = -1.
        ddkp1 = -1.
        
        w1   =    (dj0  * djp1 * dk0  * dkp1)*onefourth
        w2   = -2*(djm1 * djp1 * dk0  * dkp1)*onefourth
        w3   = -  (djm1 * dj0  * dk0  * dkp1)*onefourth
        w4   = -2*(dj0  * djp1 * dkm1 * dkp1)*onefourth
        w5   =  4*(djm1 * djp1 * dkm1 * dkp1)*onefourth
        w6   =  2*(djm1 * dj0  * dkm1 * dkp1)*onefourth
        w7   = -  (dj0  * djp1 * dkm1 * dk0 )*onefourth
        w8   =  2*(djm1 * djp1 * dkm1 * dk0 )*onefourth
        w9   =    (djm1 * dj0  * dkm1 * dk0 )*onefourth

        dw11 =    ((ddj0  * djp1  * dk0  * dkp1) + &
                   (dj0   * ddjp1 * dk0  * dkp1))*onefourth
        dw21 = -2*((ddjm1 * djp1  * dk0  * dkp1) + &
                   (djm1  * ddjp1 * dk0  * dkp1))*onefourth
        dw31 = -  ((ddjm1 * dj0   * dk0  * dkp1) + &
                   (djm1  * ddj0  * dk0  * dkp1))*onefourth
        dw41 = -2*((ddj0  * djp1  * dkm1 * dkp1) + &
                   (dj0   * ddjp1 * dkm1 * dkp1))*onefourth
        dw51 =  4*((ddjm1 * djp1  * dkm1 * dkp1) + &
                   (djm1  * ddjp1 * dkm1 * dkp1))*onefourth
        dw61 =  2*((ddjm1 * dj0   * dkm1 * dkp1) + &
                   (djm1  * ddj0  * dkm1 * dkp1))*onefourth
        dw71 = -  ((ddj0  * djp1  * dkm1 * dk0 ) + &
                   (dj0   * ddjp1 * dkm1 * dk0 ))*onefourth
        dw81 =  2*((ddjm1 * djp1  * dkm1 * dk0 ) + &
                   (djm1  * ddjp1 * dkm1 * dk0 ))*onefourth
        dw91 =    ((ddjm1 * dj0   * dkm1 * dk0 ) + &
                   (djm1  * ddj0  * dkm1 * dk0 ))*onefourth

        dw12 =    ((dj0  * djp1 * ddk0  * dkp1) + &
                   (dj0  * djp1 * dk0   * ddkp1))*onefourth
        dw22 = -2*((djm1 * djp1 * ddk0  * dkp1) + &
                   (djm1 * djp1 * dk0   * ddkp1))*onefourth
        dw32 = -  ((djm1 * dj0  * ddk0  * dkp1) + &
                   (djm1 * dj0  * dk0   * ddkp1))*onefourth
        dw42 = -2*((dj0  * djp1 * ddkm1 * dkp1) + &
                   (dj0  * djp1 * dkm1  * ddkp1))*onefourth
        dw52 =  4*((djm1 * djp1 * ddkm1 * dkp1) + &
                   (djm1 * djp1 * dkm1  * ddkp1))*onefourth
        dw62 =  2*((djm1 * dj0  * ddkm1 * dkp1) + &
                   (djm1 * dj0  * dkm1  * ddkp1))*onefourth
        dw72 = -  ((dj0  * djp1 * ddkm1 * dk0 ) + &
                   (dj0  * djp1 * dkm1  * ddk0 ))*onefourth
        dw82 =  2*((djm1 * djp1 * ddkm1 * dk0 ) + &
                   (djm1 * djp1 * dkm1  * ddk0 ))*onefourth
        dw92 =    ((djm1 * dj0  * ddkm1 * dk0 ) + &
                   (djm1 * dj0  * dkm1  * ddk0 ))*onefourth

        do n=1,2
          b(n,1) =    &
                    dw11*xv(1,n)  &
                  + dw21*xv(2,n)  &
                  + dw31*xv(3,n)  &
                  + dw41*xv(4,n)  &
                  + dw51*xv(5,n)  &
                  + dw61*xv(6,n)  &
                  + dw71*xv(7,n)  &
                  + dw81*xv(8,n)  &
                  + dw91*xv(9,n)

          b(n,2) =    &
                    dw12*xv(1,n)  &
                  + dw22*xv(2,n)  &
                  + dw32*xv(3,n)  &
                  + dw42*xv(4,n)  &
                  + dw52*xv(5,n)  &
                  + dw62*xv(6,n)  &
                  + dw72*xv(7,n)  &
                  + dw82*xv(8,n)  &
                  + dw92*xv(9,n)

          b(n,3) = xp(n) -  & 
                  ( w1*xv(1,n)   &
                  + w2*xv(2,n)   &
                  + w3*xv(3,n)   &
                  + w4*xv(4,n)   &
                  + w5*xv(5,n)   &
                  + w6*xv(6,n)   &
                  + w7*xv(7,n)   &
                  + w8*xv(8,n)   &
                  + w9*xv(9,n))

        enddo
        call matrixinv2x2(b)
        do n=1,2
          frac(n) = frac(n) + b(n,3)
          if((frac(n)-1.)*(frac(n)+1.) .gt. 0.) frac(n) = .5
        enddo

          resid = 0.
          do n=1,2
            resid = resid + (b(n,3))**2
          enddo

          niter = niter + 1
      enddo

!...fix fractions if they don't fall in [-1 to 1] range.

      if((frac(1)+1.0001)*(frac(1)-1.0001).gt.0. .or. &
         (frac(2)+1.0001)*(frac(2)-1.0001).gt.0. ) then
        if(frac(1) .lt. -1.0001) frac(1) = -0.99
        if(frac(2) .lt. -1.0001) frac(2) = -0.99
        if(frac(1) .gt. 1.0001) frac(1) = 0.99
        if(frac(2) .gt. 1.0001) frac(2) = 0.99
      end if

      end subroutine find_fractions

!********************************************************************
        subroutine bilinearinterp(ideriv,dpsi,n,f,fint)

!!  bilinear interpolation
!!  (returns fint)

!********************************************************************
        real dpsi(2),f(2,2,2),fint
        integer j,k,m, n,ideriv
        real g(2,2)

        do m=1,2
        if(m .ne. ideriv) then
          g(1,m) = 1 - dpsi(m)
          g(2,m) = dpsi(m)
        else
          g(1,m) = -1
          g(2,m) = 1
        endif
        enddo

        fint = 0.
        do k=1,2
        do j=1,2
          fint = fint + f(j,k,n) * g(j,1)*g(k,2)
        enddo
        enddo

        end subroutine bilinearinterp

!********************************************************************
        subroutine matrixinv2x2(a)
!
!  cramer's method to solve ax=b
!  (both utilizes & returns a)
!!
!! a(2,3) --- matrix to be inverted
!! 
!! a(:,3) contains the right hand side vector. the final solution is also stored
!! in a(:,3).
!!

!********************************************************************

        real a(2,3),x,y,det
!
        det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
        x = a(1,3)/det
        y = a(2,3)/det
        a(1,3) =   x*a(2,2) - y*a(1,2)
        a(2,3) = - x*a(2,1) + y*a(1,1)
!
        end subroutine matrixinv2x2

!********************************************************************

end module ihc

!********************************************************************
