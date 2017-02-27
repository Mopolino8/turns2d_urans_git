c**********************************************************************
      subroutine compute_forces_bq(jd, kd, x, y, xv, yv, q, qb, xx, xy, 
     +                             yx, yy, cfx_tot, cfx_totb, cfy_tot, 
     +                             cfy_totb, cm_tot, cm_totb, cl_tot, 
     +                             cl_totb, cd_tot, cd_totb, cpower_tot
     +                             , cpower_totb, im)
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer jd, kd, im
      real q(jd, kd, nspec, nq), x(jd, kd, nspec), y(jd, kd, nspec)
      real qb(jd, kd, nspec, nq)
      real xv(jmax, kmax, nspec), yv(jmax, kmax, nspec)
      real xx(jd, kd, nspec), xy(jd, kd, nspec)
      real yx(jd, kd, nspec), yy(jd, kd, nspec)
c local variables
      real cfx_tot, cfy_tot, cm_tot, cl_tot, cd_tot, cpower_tot
      real cfx_totb, cfy_totb, cm_totb, cl_totb, cd_totb, cpower_totb
c
      integer nsp
      real x0, y0, cfx, cfy, cm0, cpower
      real cfxb, cfyb, cm0b, clb, cdb, cmb, cpowerb
      real ca, sa
c
c***  first executable statement
      pi = 4.0*atan(1.0)
      ca = cos(pi*alfa/180.0)
      sa = sin(pi*alfa/180.0)
c
      x0 = 0.25
      y0 = 0.00
c
!$OMP PARALLEL IF(NSPEC > 1)
!$OMP& PRIVATE(cfxb,cfyb,cm0b,clb,cdb,cmb,cpowerb)
!$OMP DO
      do nsp=nspec,1,-1

        !am cpowerb = cpower_totb
        !am cmb = cm_totb
        !am cdb = cd_totb
        !am clb = cl_totb
        if(if_obj .and. if_objtot) then
          print '(A,I7,3(x,E20.12))','nsp, airloads: ',
     >          nsp,clobj(nsp),cdobj(nsp),cmobj(nsp)
          cpowerb = cpower_totb
          if(obj_ftype.eq.obj_ftype_cltot) then
            cmb = 2.0*cmobj(nsp)*cm_totb
            cdb = 2.0*cdobj(nsp)*cd_totb
            clb = 2.0*(clobj(nsp)-obj_f(nsp))*cl_totb
          elseif(obj_ftype.eq.obj_ftype_cdtot) then
            cmb = 2.0*cmobj(nsp)*cm_totb
            cdb = 2.0*(cdobj(nsp)-obj_f(nsp))*cd_totb
            clb = 2.0*clobj(nsp)*cl_totb
          elseif(obj_ftype.eq.obj_ftype_cmtot) then
            cmb = 2.0*(cmobj(nsp)-obj_f(nsp))*cm_totb
            cdb = 2.0*cdobj(nsp)*cd_totb
            clb = 2.0*clobj(nsp)*cl_totb
          else
            print '(A,I7,x,A)','OBJ TYPE: ',obj_ftype,' not supported '
            stop 
          end if
        else
          cpowerb = cpower_totb
          cmb = cm_totb
          cdb = cd_totb
          clb = cl_totb
        endif
        cfyb = x0*cmb + ca*clb + sa*cdb + cfy_totb
        cfxb = ca*cdb - sa*clb - y0*cmb + cfx_totb
        cm0b = cmb
        call force2d_bq(jd, kd, x(:, :, nsp), y(:, :, nsp), xv(:, :, nsp
     +                  ), yv(:, :, nsp), q(:, :, nsp, :), qb(:, :, nsp
     +                  , :), xx(:, :, nsp), xy(:, :, nsp), yx(:, :, nsp
     +                  ), yy(:, :, nsp), cfx, cfxb, cfy, cfyb, cm0, 
     +                  cm0b, cpower, cpowerb)
      enddo
!$OMP END DO
!$OMP END PARALLEL

      cfx_totb = 0.0
      cm_totb = 0.0
      cfy_totb = 0.0
      cd_totb = 0.0
      cpower_totb = 0.0
      cl_totb = 0.0

      end subroutine compute_forces_bq 

c**********************************************************************
