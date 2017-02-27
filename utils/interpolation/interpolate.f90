!********************************************************************
program interpolate
!********************************************************************
   use input
   use grid
   use solution
   use connectivity
   implicit none

   logical :: flowsolution

   call read_inputs()
   call init_grid()    !if no spatial interpolation, only initializes dimension

   if (interpolate_space) then
     call connect()    ! call connectivity
     call write_grid() ! write fine cell-centered grid for viewing interpolated data
   endif

!..Read in flow solution and interpolate
   flowsolution = .true. 
   call init_solution(flowsolution)
   call interpolate_solution(flowsolution)
   call write_output()

   if (interpolate_adjoints) then
!..Read in adjoint solution and interpolate
     flowsolution = .false.
     call init_solution(flowsolution)
     call interpolate_solution(flowsolution)
     call write_output()
   endif

end program interpolate

!********************************************************************
subroutine interpolate_solution(flowsolution)
!********************************************************************
   use input
   use solution
   implicit none

   logical :: flowsolution

   if (interpolate_space) then
     call interpolate_in_space(flowsolution)
   endif

   if (interpolate_time) then
     call interpolate_in_time()
   endif

   if (interpolate_stochastic) then
     call interpolate_in_stochastic()
   endif

end subroutine interpolate_solution

!********************************************************************
