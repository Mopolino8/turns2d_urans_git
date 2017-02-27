
      subroutine objectivef_ts(jd,kd,q,opt_obj,im)
      use params_global
      implicit none
      
      integer :: jd,kd,im
      real :: q(jd,kd,nspec,nq)
      real :: opt_obj
      
      !local variables
      real :: coeff(jd), tmp_norm
      integer :: j, k, np,nobjf,nsp


      opt_obj = 0.0 !sums up all Cp target points for all time instances
      coeff  = 0.0
      !tmp_norm = 0.0
!$OMP PARALLEL REDUCTION(+:opt_obj) 
!$OMP& IF(NSPEC > 1)
!$OMP& PRIVATE(coeff,tmp_norm,np,nobjf)
!$OMP DO ORDERED
      spectralloop: do nsp = 1,nspec
      
        if (obj_ftype .eq. obj_ftype_cp) then
           call calc_cp(jd, kd, q(:,:,nsp,:), coeff)
        else
           print *, "The objective type is not defined"
           stop
        end if
        

        if (obj_pointwise) then
           !am opt_obj = coeff(obj_points(1)) -  obj_f(obj_points(1)) !!serial
           opt_obj = opt_obj + coeff(obj_points(1)) -  obj_f(nsp) ! 1 Cp target for each nsp
        else
           do np = 1, obj_npoints
              nobjf = jd*(nsp-1)+obj_points(np) !jd per each TS
              call calc_norm(coeff(obj_points(np)), obj_f(nobjf), obj_normtype, tmp_norm)
              opt_obj = opt_obj + tmp_norm
              !opt_obj = opt_obj + (coeff(obj_points(np)) - obj_f(nobjf))**2
           end do

           !possible bug with ovset+TS !! 
           !needs to go out since this is called only when im=1
           !----------------------------------------------------
           if(obj_if_reg .and. nsp.eq.nspec) then 
              do k=1,kd                           
                 do j=1,jd
                    opt_obj = opt_obj + obj_reg_fac*(obj_coeff_prod(j,k,im) - 1.0)**2
                 end do
              end do
           end if
           !----------------------------------------------------

        end if

      enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL
      print *, "Objective function value: ", opt_obj
      write(747, *) opt_obj
      end subroutine  objectivef_ts
