!********************************************************************
module octree
  implicit none
!********************************************************************

  integer :: lvlmin,lvlmax
  integer,private :: nbox
  real,pointer :: xs(:,:),xe(:,:),scale(:,:)
  logical,pointer :: empty_flag(:,:)
  integer,pointer :: lvlpoint(:,:)
  integer,pointer :: istart(:,:,:),iend(:,:,:)

  integer,allocatable,private :: check_alloc(:)

  private :: init_octree,localboxindex,index_ltog,interleave,parent
contains

!********************************************************************
  subroutine init_octree(ngrids,llmin,llmax)
   
    integer :: ngrids,llmin,llmax,lvl
  
    lvlmin = llmin
    lvlmax = llmax

    nbox = 0
    do lvl = 0,lvlmax
      nbox = nbox + ishft(1,2*lvl)
    enddo
   
    nullify(xs,xe,scale,empty_flag,lvlpoint,istart,iend)
    if(allocated(check_alloc)) deallocate(check_alloc)

    allocate(xs(2,ngrids),xe(2,ngrids),scale(2,ngrids))
    allocate(empty_flag(nbox,ngrids),lvlpoint(nbox,ngrids))
    allocate(istart(nbox,2,ngrids),iend(nbox,2,ngrids))

    allocate(check_alloc(1))
  
  end subroutine init_octree

!********************************************************************
  subroutine make_octree(x,y,jmax,kmax,ngrids,llmin,llmax,init)

!********************************************************************

    integer,optional :: init
    integer :: ngrids,llmin,llmax
    integer :: jmax(:),kmax(:)
    real    :: x(:,:,:),y(:,:,:)

!.. local variables
    integer :: j,k,n
    integer :: lvl,nout,npout,nmout,nbox,npbox,nsbox,nebox
    integer :: maxbox,maxlev
    real    :: eps_x,eps_y
    real    :: xbar(2)

!...allocate variables under user request

    if(present(init)) then
      if (init.eq.1) call init_octree(ngrids,llmin,llmax)
    endif

!...allocate variables if not initiated yet

    if (.not.allocated(check_alloc)) call init_octree(ngrids,llmin,llmax)

    empty_flag = .true.
    istart = 1
    iend = 1

!... scaling

    do n = 1,ngrids
      xs(1,n) = minval(x(1:jmax(n),1:kmax(n),n))
      xe(1,n) = maxval(x(1:jmax(n),1:kmax(n),n))

      xs(2,n) = minval(y(1:jmax(n),1:kmax(n),n))
      xe(2,n) = maxval(y(1:jmax(n),1:kmax(n),n))

      eps_x=(xe(1,n)-xs(1,n))*0.0001
      eps_y=(xe(2,n)-xs(2,n))*0.0001
      
      xs(1,n) = xs(1,n) - eps_x
      xs(2,n) = xs(2,n) - eps_y
      xe(1,n) = xe(1,n) + eps_x
      xe(2,n) = xe(2,n) + eps_y

      scale(1,n) = xe(1,n) - xs(1,n)
      scale(2,n) = xe(2,n) - xs(2,n)

      maxlev = 10

      do j=1,jmax(n)
      do k=1,kmax(n)

        xbar(1) = (x(j,k,n)-xs(1,n))/scale(1,n)
        xbar(2) = (y(j,k,n)-xs(2,n))/scale(2,n)

        lvl = lvlmax
        call boxindex(xbar,lvl,nout,2)
        if (empty_flag(nout,n)) then
          istart(nout,1,n) = j 
          istart(nout,2,n) = k 
          call boxindex(xbar,maxlev,nsbox,2)
 
          iend(nout,1,n) = j 
          iend(nout,2,n) = k 
          nebox = nsbox
          empty_flag(nout,n) = .false.
        else
          call boxindex(xbar,maxlev,nmout,2)
          if (nmout.lt.nsbox) then
            istart(nout,1,n) = j 
            istart(nout,2,n) = k 
            nsbox = nmout
          elseif (nmout.gt.nebox) then
            iend(nout,1,n) = j 
            iend(nout,2,n) = k 
            nebox = nmout
          endif
        endif

      enddo
      enddo

      do lvl = lvlmax,1,-1
        maxbox = 2**(2*lvl)
        do nbox = 0,maxbox-1
          call parent(nbox,npbox,2)
          call index_ltog(nbox,nout,lvl,2)
          call index_ltog(npbox,npout,lvl-1,2)

          if (.not.empty_flag(nout,n)) then
            if (empty_flag(npout,n)) then
              istart(npout,1,n) = istart(nout,1,n) 
              istart(npout,2,n) = istart(nout,2,n) 
            endif
            empty_flag(npout,n) = .false.
            iend(npout,1,n) = iend(nout,1,n) 
            iend(npout,2,n) = iend(nout,2,n) 
          endif
        enddo
      enddo

      maxbox = 2**(2*lvlmax)

      do nbox = 0,maxbox-1
        call index_ltog(nbox,nout,lvlmax,2)
        lvlpoint(nout,n) = lvlmax

        if (empty_flag(nout,n)) then
          npbox = nbox

          do lvl = lvlmax-1,lvlmin,-1
            call parent(npbox,npbox,2)
            call index_ltog(npbox,npout,lvl,2)
            if (empty_flag(npout,n)) then
              cycle
            else
              lvlpoint(nout,n) = lvl
              exit
            endif
          enddo 
        endif

      enddo

    enddo

  end subroutine make_octree

!*********************************************************************
  subroutine boxindex(x,l,n,d)
! finding the global box index given a point
!*********************************************************************
    integer :: n,l,d
    real :: x(d)

    call localboxindex(x,l,n,d)
    call index_ltog(n,n,l,d)

    end subroutine boxindex

!*********************************************************************
  subroutine localboxindex(x,l,n,d)
! finding the box index given a point
! This gives per level box number
!*********************************************************************
    integer :: n,l,d
    real :: x(d)

!.. local variables
    integer :: dim
    integer :: ncoord(d)

    do dim = 1,d
      ncoord(dim) = (ishft(1,l)*x(dim))
    enddo
    call interleave(ncoord,n,d)

    end subroutine localboxindex

!*********************************************************************
  subroutine index_ltog(n,m,l,d)
! convert local box index to global 
!*********************************************************************
    integer :: n,m,l,d
    
    m = n + 1
    do l = 0,l-1
      m = m + ishft(1,d*l)
    enddo
    
  end subroutine index_ltog

!*********************************************************************
  subroutine interleave(n,nout,d)
! subroutine for interleaving
!*********************************************************************
    integer :: nout,d
    integer :: n(d)

!.. local variables
    integer :: ntemp(d)
    integer :: dim,count,bit

    ntemp = n
    nout = 0
    count = 0

    do while (maxval(ntemp).ne.0)
      do dim = d,1,-1
        bit = ntemp(dim) - 2*ishft(ntemp(dim),-1)
        nout = nout + ishft(bit,count)
        count = count+1
        ntemp(dim) = ishft(ntemp(dim),-1)
      enddo
    enddo

  end subroutine interleave

!********************************************************************
  subroutine parent(n,m,d)
!
! finding the parent
!*********************************************************************
   integer :: n,m,d

   m = ishft(n,-d)

   end subroutine parent

end module octree

!********************************************************************
