C   Adjoint version of objectivef_ts 
C   Author: Asitav (Sep 3, 2015)
      subroutine objectivef_ts_ad(jd, kd, q, qb, opt_obj, opt_objb)
      use params_global
      use params_sensitivity
      implicit none
      include 'DIFFSIZES.inc'
      integer jd, kd
      real q(jd, kd, nspec, nq)
      real qb(jd, kd, nspec, nq)
      real opt_obj
      real opt_objb

      !local variables
      real coeff(jd), tmp_norm
      real tmp_normb
      integer j, k, np,nobjf,nsp
      real coeffb(jd)
c
!$OMP  PARALLEL IF(NSPEC > 1)
!$OMP& PRIVATE(coeff,coeffb,tmp_norm,tmp_normb,np,nobjf)
!$OMP DO ORDERED
      spectralloop: do nsp = nspec,1,-1
C
        IF (obj_ftype .EQ. obj_ftype_cp) THEN
          CALL CALC_CP(jd, kd, q(:,:,nsp,:), coeff)
C
          IF (obj_pointwise) THEN
            coeffb = 0.0
            !am coeffb(obj_points(1)) = coeffb(obj_points(1)) + opt_objb
            coeffb(obj_points(1)) = opt_objb
          ELSE
            coeffb = 0.0
            DO np=obj_npoints,1,-1
              nobjf = jd*(nsp-1)+obj_points(np) !jd per each TS
              tmp_normb = opt_objb
              CALL CALC_NORM_BO(coeff(obj_points(np)), 
     +                          coeffb(obj_points(np)), obj_f(nobjf), 
     +                          obj_normtype, tmp_norm, tmp_normb)
            ENDDO
          END IF
          CALL CALC_CP_BO(jd, kd, q(:,:,nsp,:), qb(:,:,nsp,:), coeff, coeffb)
        ELSE
          STOP
        END IF

      enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL

      opt_objb = 0.0

      END subroutine objectivef_ts_ad
