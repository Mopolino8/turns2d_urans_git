!********************************************************************
module solution
  use input
  use grid
  use interpolators
  implicit none

  integer :: nv
  integer,allocatable :: jdim(:),kdim(:)
  real,allocatable :: q(:,:,:,:,:,:)
  real :: fsmach,alfa,rey,totime
  integer :: istep
  character*40,allocatable :: soln_prefix(:)
  character*40 :: suffix

  data nv /4/
 
  contains

!********************************************************************
  subroutine init_solution(flowsolution)
!********************************************************************

    logical flowsolution

    integer j,k,jtmp,ktmp,n,ntime,nstoch,nst,nsp
    character*40, append1,append,sfile,integer_string

    if(allocated(jdim)) deallocate(jdim)
    if(allocated(kdim)) deallocate(kdim)
    if(allocated(soln_prefix)) deallocate(soln_prefix)
    allocate(jdim(2),kdim(2))
    allocate(soln_prefix(2))

    soln_prefix(1) = 'fort'
    soln_prefix(2) = 'finefort'
    
    if (flowsolution) then
      suffix='.8'
    else
      suffix='.18'
    endif

    ntime = max(ntime_in,ntime_out)
    nstoch = max(nstochpts_in,nstochpts_out)
 
    if(allocated(q)) deallocate(q)
    allocate(q(jmc,kmc,ntime,nstoch,nv,2))

    append1 = ''
    append = ''
    do nst = 1,nstochpts_in
      if (interpolate_stochastic) then
        write(integer_string,*) nst
        append1 = '_S'//trim(adjustl(integer_string))
      endif

      do nsp = 1,ntime_in
        if (interpolate_time) then
          write(integer_string,*) nsp
          append = trim(adjustl(append1))//'_TS'//trim(adjustl(integer_string))
        endif

        sfile = &
        trim(adjustl(soln_prefix(1)))//trim(adjustl(append))//trim(adjustl(suffix))

        open(unit=1,file=sfile,form='unformatted',status='unknown') 
        read(1) jtmp,ktmp
        jdim(1) = jtmp; kdim(1) = ktmp
        read(1) fsmach,alfa,rey,totime
        read(1) (((q(j,k,nsp,nst,n,1),j=1,jtmp),k=1,ktmp),n=1,nv)
        read(1) istep
        close(1)
      enddo
    enddo

    if (treset) then
      totime = 0.0
      istep = 0
    endif
    
    if (interpolate_space) then
      jdim(2) = jd(2)
      kdim(2) = kd(2)
    else
      jdim(2) = jdim(1)
      kdim(2) = kdim(1)
    endif

  end subroutine init_solution

!********************************************************************
  subroutine interpolate_in_space(flowsolution)
!********************************************************************

    integer :: nsp,nst
    logical :: flowsolution

    if (flowsolution) then
      do nsp = 1,ntime_in
        do nst = 1,nstochpts_in
          call spatial_interpolator_muscl(q(:,:,nsp,nst,:,:),nv)
        enddo
      enddo
    else
      do nsp = 1,ntime_in
        do nst = 1,nstochpts_in
          call spatial_interpolator_bicubic(q(:,:,nsp,nst,:,:),nv)
        enddo
      enddo
    endif

    q(:,:,:,:,:,1) = q(:,:,:,:,:,2)

  end subroutine interpolate_in_space

!********************************************************************
  subroutine interpolate_in_time()
!********************************************************************

    integer :: j,k,n,nst
    
    if (ntime_in.eq.ntime_out) return

    do n = 1,nv
      do nst = 1,nstochpts_in
        do k = 1,kdim(2)
          do j = 1,jdim(2)
            call fourier_interpolator(q(j,k,:,nst,n,1),q(j,k,:,nst,n,2), &
                                      ntime_in,ntime_out)
          enddo
        enddo
      enddo
    enddo

    q(:,:,:,:,:,1) = q(:,:,:,:,:,2)

  end subroutine interpolate_in_time

!********************************************************************
  subroutine interpolate_in_stochastic()
!********************************************************************
    integer :: j,k,n,nsp
    
    do n = 1,nv
      do nsp = 1,ntime_out
        do k = 1,kdim(2)
          do j = 1,jdim(2)
            call spline_interpolator(stochpts_in,q(j,k,nsp,:,n,1),  &
                                     stochpts_out,q(j,k,nsp,:,n,2), &
                                     nstochpts_in,nstochpts_out)
          enddo
        enddo
      enddo
    enddo

    q(:,:,:,:,:,1) = q(:,:,:,:,:,2)

  end subroutine interpolate_in_stochastic

!********************************************************************
  subroutine write_output()
!********************************************************************

   integer j,k,jtmp,ktmp,n,nst,nsp
   character*40, append1,append,sfile,integer_string

   append1 = ''
   append = ''
   do nst = 1,nstochpts_out
     if (interpolate_stochastic) then
       write(integer_string,*) nst
       append1 = '_S'//trim(adjustl(integer_string))
     endif

     do nsp = 1,ntime_out
       if (interpolate_time) then
         write(integer_string,*) nsp
         append = trim(adjustl(append1))//'_TS'//trim(adjustl(integer_string))
       endif

       sfile = &
       trim(adjustl(soln_prefix(2)))//trim(adjustl(append))//trim(adjustl(suffix))

       open(unit=1,file=sfile,form='unformatted',status='unknown') 
       jtmp = jdim(2); ktmp = kdim(2)
       write(1) jtmp,ktmp
       write(1) fsmach,alfa,rey,totime
!..Printing q(:,:,:,:,:,1), since this should be same as q(:,:,:,:,:,2) at this point
!..Also takes care of the situation when no interpolation is requested
       write(1) (((q(j,k,nsp,nst,n,1),j=1,jtmp),k=1,ktmp),n=1,nv)
       write(1) istep
       close(1)
     enddo
   enddo

  end subroutine write_output

end module solution

!********************************************************************
