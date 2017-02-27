module input

  logical :: interpolate_space, interpolate_time, interpolate_stochastic 
  logical :: interpolate_adjoints, treset
  
  integer :: ntime_in, ntime_out
  
  integer :: nstochpts_in, nstochpts_out
  real    :: stochpts_in(25), stochpts_out(25)

  integer :: nhalo
  
  data interpolate_space, interpolate_time, interpolate_stochastic &
  /.false.,.false.,.false./
  data interpolate_adjoints, treset /.true., .false./
  data ntime_in, ntime_out /1,1/
  data nstochpts_in, nstochpts_out /1,1/
  data nhalo /2/
  
  namelist/inputs/ interpolate_space, interpolate_time, interpolate_stochastic, &
    interpolate_adjoints, treset, &
    ntime_in, ntime_out, nstochpts_in, nstochpts_out, stochpts_in, stochpts_out, &
    nhalo

contains

!********************************************************************
subroutine read_inputs
!********************************************************************

  read(5,inputs)

end subroutine read_inputs

!********************************************************************
end module input
