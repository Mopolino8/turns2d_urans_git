c********1*********2*********3*********4*********5*********6*********7**
      program umturns2d
c                                                                      c
c  the code was cleaned-up primarily by j. sitaraman to f90 type 
c  constructs and  options for oversetting and other good stuff to work
c
c  last modified 09/05/2005 by karthik.
c
c  j.d baeder did previously claim to have cleaned it up in 1991 :)    c
c  to structure, vectorize and streamline to make all options working. c
c                                                                      c
c  many new limiters etc. added in 1996, along with better metrics     c
c                                                                      c
c**********************************************************************c
c
c
c     note: b.c. elements must be configured for each new grid topology
c           (currently c-grid)
c     note: mdim must be .ge. jdim and kdim
c
c         tape1:    grid (input)
c         tape3:    restart file (input)
c         tape5:    input file
c         tape6:    output file of summary data
c         tape7:    output file of run norms
c         tape8:    restart file (output - suitable for plot3d)
c         tape9:    grid (output)
c         tape11:   output file of cl,cd,cmp
c         tape17:   file of alpha(t) (input)
c 
c*********************************************************************** 
c 
c   input variable to airfoil in namelist inputs:
c   
c     iread  = tells whether to read in a restart file from unit 3
c            = 1 => restart
c            = 0 => initial run with no restart
c     jmax   = points in wrap around direction
c     kmax   = points in normal direction
c     jtail1 = j location of beginning airfoil (at tail on underside)
c     half   = symmetrical airfoil?
c            = 0 no symmetry assumed, 
c                        jtail2 = jmax-jtail1+1 & jle =(jmax+1)/2
c            = 1 symmetry assumed, jtail2 = jmax-1 & jle = jmax-1
c   
c   
c     npnorm = output residual to fort.7 every npnorm iterations
c              output force to fort.11 every npnorm its.
c              output spanwise forces to fort.12 every npnorm its. (rotor only)
c     nrest  = write out restart to fort.8 every nrest iterations
c     nsteps = total number of steps at end of this run
c              (typically 200-2000 larger than last run)
c   
c   
c     fsmach = free stream mach number (0 for hover)
c     alfa   = angle of attack for freestream (collective set by grid)
c     rey    = reynolds number
c     invisc = euler or navier-stokes
c            = .false. => navier-stokes
c            = .true. => euler
c     lamin  = is the flow laminar
c            = .true. => laminar
c            = .false. => fully turbulent
c   
c   
c     iunst  = type of unsteady flow
c            = 0 => steady
c            = 1 => unsteady pitching oscillation
c            = 2 => unsteady pitching ramp
c                         (need to input grid motion)
c            = 3 => prescribed pitching
c     ntac   = temporal order of accuracy
c            = 1 => first order
c            = 2 => second order
c            = 3 => third order 
c     itnmax = number of newton iterations per time step
c     dt     = time step size, should be less than 0.10 for time accuracy
c              (typically 50.0 if steady, 0.05 for time-accurate)
c     timeac = how does dt vary in space
c            = 1. => constant time step everywhere
c            otherwise space-varying dt
c     totime = total time for unsteady calculation (reset from restart file)
c   
c   
c     ilhs   = left hand side scheme
c            = 1 => LU-SGS algorithm
c            = 2 => ARC-2D with second and fourth order dissipation
c            = 3 => ARC-2D with upwind
c    iprecon = use low Mach preconditioning or not
c            = .true. => use preconditioning
c            = .false. => no preconditioning
c     Mp    = low Mach preconditioning factor
c     epse   = dissipation for implicit side (usually 0.01)
c     irhsy  = spatial order of accuracy
c            = -1 => 1st order
c            = -2 => 2nd order
c            = -3 => 3rd order
c     ilim   = type of limiting
c            =  0 => no limiting at all (muscl scheme)
c            <  0 => no limiting in k-direction (muscl scheme)
c            >  0 => limiting in both directions
c abs(ilim)  =  1 => differentiable limiter for 3rd order (koren's)
c            =  2 => differentiable limiter for 2nd order (van albada)
c            =  3 => chakravarthy-osher minmod
c            =  4 => sonic-a scheme of hyunh et al.
c            =  5 => sonic extension of chakravarthy-osher minmod
c            =  7 => cubic interpolation with no limiting
c            =  8 => quartic interpolation with no limiting
c            =  9 => quadratic reconstruction with no limiting (pade' scheme)
c    
c
c     jint   = calculate solution on only every jint points in j-direction (1)
c     kint   = calculate solution on only every kint points in k-direction (1)
c    
c     rf     = reduced frequency or time for unsteady motion
c     angmax = maximum change in angle of attack
c 
c************end prologue ********************************************
      use params_global
      use params_adjoints
      use params_sensitivity
      use ihc
c*********************************************************************
      implicit none
c*********************************************************************
      ! allocatable arrays

      real, allocatable :: s(:),q(:),sb(:),qb(:),qb1(:),qtn(:),qtnm1(:),qnewt(:)
      real, allocatable :: x(:),y(:),xv(:),yv(:)
      real, allocatable :: tscale(:),bt(:)
      integer, allocatable :: iblank(:)
      real, allocatable :: xx(:),xy(:),yx(:),yy(:)
      real, allocatable :: ug(:),vg(:),ugv(:),vgv(:),turmu(:)
      real, allocatable :: xole(:),yole(:),xold(:),yold(:)
      integer, allocatable :: jgmx(:),kgmx(:),jgmxv(:),kgmxv(:)
      integer, allocatable :: ipointer(:,:)

      real,allocatable :: tspec(:),Ds(:,:)

      ! arrays for overset meshing

      integer              :: Nj,Nk,idsize,j,k
      integer,allocatable  :: ndonor(:,:),nfringe(:,:)
      integer,allocatable  :: iisptr(:,:),iieptr(:,:)
      integer, allocatable :: imesh(:,:,:,:),idonor(:,:,:,:)
      integer, allocatable :: ibc(:,:,:)
      real, allocatable    :: frac(:,:,:,:)
      real, allocatable    :: xgl(:,:,:),ygl(:,:,:)
      real, allocatable    :: xglv(:,:,:),yglv(:,:,:),volg(:,:,:)
      integer, allocatable :: ibgl(:,:,:)
      character*40,allocatable :: ihcbcfile(:)
      !am
      integer              :: maxjd,maxkd

      ! local scalar variables
      
      integer ig,igq,igs,igv,igb
      integer nmesh,jd,kd,im,nsp
      integer ii,mstop
      character*40 bcfile,fprefix,integer_string

      real resmax,resrho,rsum,resold
      real cfx_tot,cfy_tot,cl_tot,cd_tot,cm_tot,cpower_tot
      logical file_exists
      real opt_objb, opt_obj
c** first executable statement

      nmesh=6
      allocate(jgmx(nmesh),kgmx(nmesh)) ! allocate mesh pointers with 
      allocate(jgmxv(nmesh),kgmxv(nmesh)) ! dummy number of meshes to start
                                        ! dummy number of meshes to start

c..  initialize data and read inputs

      bcfile = "bc.inp"
      fprefix = "fortadsens"
      call read_inputs(jgmx,kgmx,nmesh,bcfile)

c..  compute time-spectral coefficients 

      allocate(tspec(nspec),Ds(nspec,nspec))
      call computeTScoefs(nspec,Ds)

c*************** memory allocation for multiple mesh pointers ***********

      allocate(ipointer(nmesh,5))
      call determine_size(jgmx,kgmx,nspec,nq,nv,ipointer,nmesh,
     &     igrd,igrdv,igrdb,iqdim,isdim,mdim,nadd)

c***************** memory allocation block for flow variables**********

      allocate(s(isdim),q(iqdim),sb(isdim),qb(iqdim),qb1(iqdim))
      allocate(qtn(isdim),qtnm1(isdim),qnewt(isdim))
      allocate(x(igrd),y(igrd),xv(igrdv),yv(igrdv))
      allocate(tscale(igrd),bt(igrd))
      allocate(iblank(igrd))
      allocate(xx(igrd),xy(igrd),yx(igrd),yy(igrd))
      allocate(ug(igrd),vg(igrd),turmu(igrd))
      allocate(ugv(igrdv),vgv(igrdv))
      allocate(xole(igrdv),yole(igrdv))
      allocate(xold(igrdv),yold(igrdv))

c***************** memory allocation block for connectivity variables**********

      if(num_grids.gt.1) then

        Nj=jgmx(1); Nk=kgmx(1); idsize=jgmx(1)*kgmx(1)
        do im=2,nmesh
          if(Nj.le.jgmx(im)) Nj=jgmx(im)
          if(Nk.le.kgmx(im)) Nk=kgmx(im)
          if(idsize.le.jgmx(im)*kgmx(im)) idsize=jgmx(im)*kgmx(im)
        end do

        allocate(xgl(Nj,Nk,nmesh),ygl(Nj,Nk,nmesh),volg(Nj,Nk,nmesh))
        allocate(ibgl(Nj,Nk,nmesh))

        allocate(imesh(idsize,2,nmesh,nspec),idonor(idsize,2,nmesh,nspec))
        allocate(frac(idsize,2,nmesh,nspec),ibc(idsize,nmesh,nspec))
        allocate(nfringe(nmesh,nspec),ndonor(nmesh,nspec))
        allocate(iisptr(nmesh,nspec),iieptr(nmesh,nspec))
        allocate(ihcbcfile(nmesh))
         
        Nj=jgmx(1)-nadd; Nk=kgmx(1)-nadd
        do im=2,nmesh
          if(Nj.le.jgmx(im)-nadd) Nj=jgmx(im)-nadd
          if(Nk.le.kgmx(im)-nadd) Nk=kgmx(im)-nadd
        end do

        allocate(xglv(Nj,Nk,nmesh),yglv(Nj,Nk,nmesh))

      endif
c********* end memory allocation block *******************************

      do im=1,nmesh
         call set_pointers_globals(im,ipointer,ig,
     &        igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
         call initia(q(igq),qb(igq),s(igs),sb(igs),
     &        x(ig),y(ig),xv(igv),yv(igv),xx(ig),xy(ig),yx(ig),
     &        yy(ig),ug(ig),vg(ig),ugv(igv),vgv(igv),turmu(ig),
     &        xold(igv),yold(igv),xole(igv),yole(igv),iblank(ig),tspec,im,jd,kd)

         !include overset info !asitav
         !-----------------------
         if (num_grids.gt.1) then
          jgmxv(im) = jmax
          kgmxv(im) = kmax

          write(integer_string,*) im-1
          integer_string = adjustl(integer_string)
          ihcbcfile(im) = 'ihcbc.'//trim(integer_string)
         endif
         !-----------------------
      enddo

      !perform overset connectivity !asitav
      !====================================
      if(num_grids.gt.1) then

!!!$OMP PARALLEL IF (NSPEC > 1)
!!!$OMP DO
!!!$OMP& PRIVATE(im,j,k,ig,igq,igs,igv,igb,jd,kd,jmax,kmax)
!!!$OMP& FIRSTPRIVATE(xgl,ygl,ibgl,volg,xglv,yglv)
        spectralloop1: do nsp = 1,nspec

c.....collect grids
        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,
     &         igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
          do k=1,kd
            do j=1,jd
              xgl(j,k,im)=x(ig-1 + jd*kd*(nsp-1) + jd*(k-1) + j)
              ygl(j,k,im)=y(ig-1 + jd*kd*(nsp-1) + jd*(k-1) + j)
              volg(j,k,im)=1./q(igq - 1 + jd*kd*nspec*(nq-1) + 
     &                                    jd*kd*(nsp-1) + jd*(k-1) + j)
            enddo
          enddo

          do k=1,kmax
            do j=1,jmax
              xglv(j,k,im)=xv(igv-1 + jmax*kmax*(nsp-1)+jmax*(k-1)+j)
              yglv(j,k,im)=yv(igv-1 + jmax*kmax*(nsp-1)+jmax*(k-1)+j)
            enddo
          enddo
        enddo

c......call connectivity now

        call connect2d(xglv,yglv,xgl,ygl,volg,ibgl,
     &        idonor(:,:,:,nsp),frac(:,:,:,nsp),imesh(:,:,:,nsp),
     &        ibc(:,:,nsp),iisptr(:,nsp),iieptr(:,nsp),ndonor(:,nsp),
     &        nfringe(:,nsp),jgmxv,kgmxv,nhalo,nmesh,ihcbcfile)
        
c......connectivity info to all meshes
        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,
     &         igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
          do k=1,kd
            do j=1,jd
              iblank(ig-1 + jd*kd*(nsp-1) + jd*(k-1) + j) = ibgl(j,k,im)
            enddo
          enddo
        enddo 

        enddo spectralloop1
!!!$OMP END DO
!!!$OMP END PARALLEL

        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,
     &         igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
          call update_halo_iblanks(iblank(ig),jd,kd)
        enddo 

      endif !num_grids.gt.1
      
      do im=1, nmesh
         call set_pointers_globals(im,ipointer,ig,
     &        igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
         !am if(if_obj .and. bodyflag(im)) then
         if(if_obj .and. im.eq.1) then
            call obj_setup(jd,kd)
         end if
      end do

      !allocate and read beta
      !----------------------
      if(if_obj) then
        maxjd = jgmx(1)
        maxkd = kgmx(1)

        do im=2,nmesh
          if(maxjd.le.jgmx(im)) maxjd = jgmx(im)
          if(maxkd.le.kgmx(im)) maxkd = kgmx(im)
        end do
        allocate(obj_coeff_prod(maxjd,maxkd,nmesh))

        !read beta
        open(2211,file='beta.opt',form='formatted')
        do im=1,nmesh
          read(2211,*)((obj_coeff_prod(j,k,im),j=1,jgmx(im)),k=1,kgmx(im)) 
        end do
      end if
      !----------------------



      !=== end connectivity ===============

c......set the adjoints objective function

      inquire(file='adjoints.inp',exist=file_exists)

      if (file_exists) then
        open(unit=10,file='adjoints.inp',status='unknown')
        read(10,adinputs)
        close(10)
        if (cfx_totb + cfy_totb + cl_totb + cd_totb + cm_totb + 
     &      cpower_totb .eq. 0.0)  then
          cl_totb = 1.0
        endif
      else
        cfx_totb = 0.0; cfy_totb = 0.0; cl_totb = 1.0; cd_totb = 0.0
        cm_totb = 0.0; cpower_totb = 0.0
      endif

      !initialize c?_totb from obj_ftype definition
      !--------------------------------------------
      if(if_obj .and. if_objtot) then
        cfx_totb = 0.0; cfy_totb = 0.0; cl_totb = 0.0; cd_totb = 0.0
        cm_totb = 0.0; cpower_totb = 0.0

        opt_objb = 1.0
        
        if(obj_ftype .eq. obj_ftype_cltot) then
          cl_totb=opt_objb
        elseif(obj_ftype .eq. obj_ftype_cdtot) then
          cd_totb=opt_objb
        elseif(obj_ftype .eq. obj_ftype_cmtot) then
          cm_totb=opt_objb
        else !default is cl
          cl_totb=opt_objb
        end if
      end if
      !--------------------------------------------

      fsmachb = 0.; fmtipb = 0.; alfab = 0.

      qb1 = 0.
      do im=1,nmesh
        reyb = 0.
        call set_pointers_globals(im,ipointer,ig,
     &       igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
        !am if (bodyflag(im)) then
        if (bodyflag(im) .and. im.eq.1) then !not for wt wall
 
          !get opt_obj from q state
          !===================================================
          call compute_forces(jd,kd,x(ig),y(ig),xv(igv),yv(igv),
     >     q(igq),xx(ig),xy(ig),yx(ig),yy(ig),
     >     cfx_tot,cfy_tot,cm_tot,cl_tot,cd_tot,cpower_tot,im,fprefix)

          !asitav write out opt objective value
          !------------------------------------
          if(if_obj .and. if_objtot) then
            if(obj_ftype.eq.obj_ftype_cltot) then
              opt_obj = cl_tot
            elseif(obj_ftype.eq.obj_ftype_cdtot) then
              opt_obj = cd_tot
            elseif(obj_ftype.eq.obj_ftype_cmtot) then
              opt_obj = cm_tot
            else
              print '(A,I7,x,A)','OBJ TYPE: ',obj_ftype,
     >                           ' not supported '
              stop
            end if
            write(747,*)opt_obj
          end if
          !------------------------------------
          !===================================================

          call compute_forces_bq(jd,kd,x(ig),y(ig),xv(igv),yv(igv),
     +                 q(igq),qb1(igq),xx(ig),xy(ig),yx(ig),yy(ig),
     +                 cfx_tot, cfx_totb, cfy_tot, 
     +                 cfy_totb, cm_tot, cm_totb, cl_tot, 
     +                 cl_totb, cd_tot, cd_totb, cpower_tot,
     +                 cpower_totb, im)
        endif
        if (if_obj .and. .not.if_objtot .and. im.eq.1) then !not for wt wall
           opt_objb = 1.0
           fsmachb = 0.0
           fmtipb = 0.0
           alfab = 0.0
           !this might cause an issue with overset mesh!!!!!! you can probably reset qb in the bodyflag if condition
           call reset_qb(jd, kd, qb1(igq))
!am           call objectivef_bo(jd, kd, x(ig), y(ig), xv(igv), yv(igv), 
!am     +          q(igq), qb1(igq), xx(ig), xy(ig), yx(ig), yy(ig), opt_obj, opt_objb)
           call objectivef_ts_ad(jd,kd,q(igq),qb1(igq),opt_obj,opt_objb)
        end if
        

        if (fmtip.eq.0) fsmachb = fsmachb - rey*reyb/fsmach
      enddo

c..now perform the flow solution for each iteration or time step 

      do im=1,nmesh
         call set_pointers_globals(im,ipointer,ig,
     &   igq,igs,igv,igb,jd,kd,jgmx,kgmx,nmesh)
      
         if (.not. iprecon) Mp = 1.0

         call time_step(q(igq),xx(ig),xy(ig),yx(ig),yy(ig),
     <          ug(ig),vg(ig),tscale(ig),bt(ig),iblank(ig),jd,kd)

         call step(q(igq),qb(igq),qb1(igq),qtn(igs),qtnm1(igs),qnewt(igs),! q (at time steps)
     &        s(igs),sb(igs),                                      ! s
     &        x(ig),y(ig),xv(igv),yv(igv),iblank(ig),              ! grid
     &        xx(ig),xy(ig),yx(ig),yy(ig),                         ! metrics
     &        ug(ig),vg(ig),ugv(igv),vgv(igv),                     ! grid velocities
     &        xold(igv),yold(igv),xole(igv),yole(igv),             ! fine mesh 
     &        turmu(ig),                                           ! eddy viscosity
     &        tscale(ig),bt(ig),Ds,im,jd,kd,resmax,resrho,rsum) ! timestep, dimensions

          write(6,*),"Sensitivities with respect to  .."
          write(6,*) "Freestream      -->",fsmachb
          write(6,*) "Angle of attack -->",alfab
          write(6,*)
         
      enddo
        
      end

c**********************************************************************
      subroutine compute_forces_bq(jd, kd, x, y, xv, yv, q, qb, xx, xy, 
     +                             yx, yy, cfx_tot, cfx_totb, cfy_tot, 
     +                             cfy_totb, cm_tot, cm_totb, cl_tot, 
     +                             cl_totb, cd_tot, cd_totb, cpower_tot
     +                             , cpower_totb, im)
c**********************************************************************
      use params_global
      use params_sensitivity
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
      real ca, sa, cab, sab
c
c***  first executable statement
      pi = 4.0*atan(1.0)
      ca = cos(pi*alfa/180.0)
      sa = sin(pi*alfa/180.0)
c
      x0 = 0.25
      y0 = 0.00
c
!$OMP PARALLEL REDUCTION(+:alfab)
!$OMP& IF(NSPEC > 1)
!$OMP& PRIVATE(cfx,cfy,cm0,cpower)
!$OMP& PRIVATE(cfxb,cfyb,cm0b,clb,cdb,cmb,cpowerb,cab,sab)
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

        call force2d(jd,kd,x(:,:,nsp),y(:,:,nsp),
     <       xv(:,:,nsp),yv(:,:,nsp),q(:,:,nsp,:),
     <       xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     <       cfx,cfy,cm0,cpower,nsp)

        cab = cfy*clb + cfx*cdb
        sab = cfy*cdb - cfx*clb
        alfab = alfab + pi*COS(pi*(alfa/180.0))*sab/180.0 - pi*SIN(pi*(
     +    alfa/180.0))*cab/180.0
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
      subroutine initia(q,qb,s,sb,x,y,xv,yv,xx,xy,yx,yy,ug,vg,ugv,vgv,turmu,
     &     xold,yold,xole,yole,iblank,tspec,im,jd,kd)

c
c  initialize the variables to default values, read inputs, check them
c  and set up the working files.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,im
      
      real q(jd,kd,nspec,nq),qb(jd,kd,nspec,nq)
      real s(jd,kd,nspec,nv),sb(jd,kd,nspec,nv)
      real x(jd,kd,nspec),y(jd,kd,nspec),xv(jmax,kmax,nspec),yv(jmax,kmax,nspec)
      real xx(jd,kd,nspec),xy(jd,kd,nspec),yx(jd,kd,nspec),yy(jd,kd,nspec)
      real ug(jd,kd,nspec),vg(jd,kd,nspec)
      real ugv(jmax,kmax,nspec),vgv(jmax,kmax,nspec)
      real turmu(jd,kd,nspec)
      integer iblank(jd,kd,nspec)
      real xold(jmax,kmax,nspec),yold(jmax,kmax,nspec)
      real xole(jmax,kmax,nspec),yole(jmax,kmax,nspec)
      real tspec(nspec)

      ! local variables

      integer i,j,k,jtmp,ktmp,nmesh,nsp,ib
      character*40, prefix,qfile,psifile,integer_string
      integer logq,logpsi
      integer lq(nspec),lpsi(nspec)

c**   first executable statement


      pi = 4.0*atan(1.)
      dang = 0.
      do 882 nsp=1,nspec
        do 882 j = 1, jd
          do 882 k=1,kd
            turmu(j,k,nsp) = 0.
            ug(j,k,nsp) = 0.
            vg(j,k,nsp) = 0.
            iblank(j,k,nsp)=1
882   continue
      
      do 884 nsp=1,nspec
        do 884 j=1,jd
          do 884 k=1,kd
            s(j,k,nsp,1) = 0.
            s(j,k,nsp,2) = 0.
            s(j,k,nsp,3) = 0.
            s(j,k,nsp,4) = 0.
884   continue

      sb = 0.
      qb = 0.

c..setup bodyflag

      bodyflag(im) = .FALSE.
      jtail(im) = 1 
      do ib=1,nbc_all(im)
        if (ibtyp_all(ib,im).eq.5) then
          bodyflag(im) = .TRUE.
          jtail(im) = jbcs_all(ib,im)
        endif
      enddo

c..setup q and x

      istep0  = 0
      do nsp=1,nspec
        !amm tspec(nsp) = 2*pi*(nsp-1)/nspec/abs(rf) !asitav
        if(rf.ne.0) then
          tspec(nsp) = 2*pi*(nsp-1)/nspec/abs(rf)
        elseif(rf_tef.ne.0 .and. motiontype.eq.3) then
          tspec(nsp) = 2*pi*(nsp-1)/nspec/abs(rf_tef)
        end if
      enddo

!      if (iread.eq.0) then
         call grid(xv(:,:,1),yv(:,:,1),iblank(:,:,1),jd,kd,1)
!      else
!        call grid(xv(:,:,1),yv(:,:,1),iblank(:,:,1),jd,kd,4)
!      endif

!$OMP PARALLEL IF(NSPEC > 1)
!$OMP DO
      do nsp = 2,nspec
        call find_coord_TS(xv(:,:,nsp),yv(:,:,nsp),xv(:,:,1),yv(:,:,1),tspec(nsp),jd,kd)
      enddo
!$OMP END DO

!$OMP DO
!$OMP& PRIVATE(j,k,logq,logpsi,prefix,qfile,psifile,integer_string)
      spectralloop: do nsp = 1,nspec

c..initialize mesh from earlier times 

        do j = 1,jmax
          do k = 1,kmax
            xold(j,k,nsp) = xv(j,k,nsp)
            yold(j,k,nsp) = yv(j,k,nsp)
            xole(j,k,nsp) = xold(j,k,nsp)
            yole(j,k,nsp) = yold(j,k,nsp)
          enddo
        enddo

        call find_cellcenter(x(:,:,nsp),y(:,:,nsp),
     &                       xv(:,:,nsp),yv(:,:,nsp),jd,kd,im)


c..compute the metrics and the jacobians using finite volume formula

        call metfv(q(:,:,nsp,:),x(:,:,nsp),y(:,:,nsp),
     &             xv(:,:,nsp),yv(:,:,nsp),
     &             xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),jd,kd,im)
      
        if (timespectral) then
          call find_timemetrics_TS(x(:,:,nsp),y(:,:,nsp),
     &                             xv(:,:,nsp),yv(:,:,nsp),tspec(nsp),
     &                             ug(:,:,nsp),vg(:,:,nsp),
     &                             ugv(:,:,nsp),vgv(:,:,nsp),jd,kd)
        endif
          
!      call qzero( q(:,:,nsp,:),jd,kd)

        logq = 8
        logpsi = 18
        if (timespectral) then
          write(integer_string,*) nsp
          prefix = 'fort_TS'//trim(adjustl(integer_string))//'.'
        else
          prefix = 'fort.'
        endif

        write(integer_string,*) logq
        qfile = trim(adjustl(prefix))//trim(adjustl(integer_string))
        write(integer_string,*) logpsi
        psifile = trim(adjustl(prefix))//trim(adjustl(integer_string))

        lq(nsp) = 1000*nsp + logq
        open(unit=lq(nsp),file=qfile,status='unknown',
     &                               form='unformatted')

        if (im.eq.1) then
          if(num_grids.gt.1) read(lq(nsp)) nmesh
          read(lq(nsp)) (jtmp,ktmp,i=1,num_grids)
        endif

        call read_solution( q(:,:,nsp,:),jd,kd,lq(nsp) )

        lpsi(nsp) = 1000*nsp + logpsi
        open(unit=lpsi(nsp),file=psifile,status='unknown',
     &                               form='unformatted')

        if (im.eq.1) then
          if(num_grids.gt.1) read(lpsi(nsp)) nmesh
          read(lpsi(nsp)) (jtmp,ktmp,i=1,num_grids)
        endif

        call restr2( sb(:,:,nsp,:),jd,kd,lpsi(nsp) )
        if(iturb.eq.1.or.iturb.eq.2) call restr_turb(sb(:,:,nsp,:),jd,kd,lpsi(nsp))

c..divide q by the jacobian, q/q5 

        call qdivj( q(:,:,nsp,:),jd,kd )
      enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL
 
      return
      end

c***********************************************************************
      subroutine step(q,qb,qb1,qtn,qtnm1,qnewt,s,sb,
     &     x,y,xv,yv,iblank,
     &     xx,xy,yx,yy,ug,vg,ugv,vgv,
     &     xold,yold,xole,yole,
     &     turmu,tscale,bt,Ds,
     &     im,jd,kd,resmax,resrho,rsum)

c  note that the boundaries are updated explicitly
c
c***********************************************************************
      use params_global
      use params_sensitivity
c***********************************************************************
      implicit none
c***********************************************************************
      integer im,jd,kd
      real resmax,resrho,rsum
      real q(jd,kd,nspec,nq), qb(jd,kd,nspec,nq), qb1(jd,kd,nspec,nq)
      real s(jd,kd,nspec,nv), sb(jd,kd,nspec,nv)
      real tscale(jd,kd,nspec), bt(jd,kd,nspec)
      real turmu(jd,kd,nspec)
      real x(jd,kd,nspec),y(jd,kd,nspec)
      real xv(jmax,kmax,nspec),yv(jmax,kmax,nspec)
      real xx(jd,kd,nspec),xy(jd,kd,nspec),yx(jd,kd,nspec)
      real yy(jd,kd,nspec),ug(jd,kd,nspec),vg(jd,kd,nspec)
      real ugv(jmax,kmax,nspec),vgv(jmax,kmax,nspec)
      real xold(jmax,kmax,nspec),yold(jmax,kmax,nspec)
      real xole(jmax,kmax,nspec),yole(jmax,kmax,nspec)
      real qtn(jd,kd,nspec,nv),qtnm1(jd,kd,nspec,nv)
      real qnewt(jd,kd,nspec,nv)
      integer iblank(jd,kd,nspec)
      real Ds(nspec,nspec)

      integer nsp
      real cs,ss,csb,ssb

      real,allocatable :: turmub(:,:,:),vmul(:,:,:),vmulb(:,:,:)

      allocate(turmub(jd,kd,nspec),vmul(jd,kd,nspec),vmulb(jd,kd,nspec))

c***  first executable statement

      turmub = 0.0
      vmulb = 0.0

c..zero s array 
      qb = 0.
      !amm qb(:,:,:,1:4) = 0. !bug?
      cs = COS(pi*alfa/180.0)
      ss = SIN(pi*alfa/180.0)

      spectralloop1: do nsp = 1,nspec

      uinfb = 0.; vinfb = 0.; einfb = 0.; reyb = 0.
c
c..compute the right hand side and store in s array
c
      call RHSUP_BQ(q(:,:,nsp,:), qb(:,:,nsp,:), 
     +            s(:,:,nsp,:), sb(:,:,nsp,:), 
     +            xx(:,:,nsp), xy(:,:,nsp), yx(:,:,nsp), yy(:,:,nsp), 
     +            x(:,:,nsp), y(:,:,nsp), xv(:,:,nsp), yv(:,:,nsp), 
     +            xold(:,:,nsp), yold(:,:,nsp), xole(:,:,nsp), yole(:,:,nsp), 
     +            iblank(:,:,nsp), ugv(:,:,nsp), vgv(:,:,nsp), 
     +            jd, kd,1,im,bt)
!      if (iunst.eq.2) call momsource(q(:,:,nsp,:),s(:,:,nsp,:),jd,kd)
!      if (axisymmetric) call axisource(q(:,:,nsp,:),s(:,:,nsp,:),x(:,:,nsp),y(:,:,nsp),jd,kd)
c
c..compute viscous fluxes
c
      if( .not. invisc ) then
        if (iturb.eq.1) then
          call lamvis(q(:,:,nsp,:),vmul(:,:,nsp),jd,kd)
          call c_turm(q(:,:,nsp,:),turmu(:,:,nsp),vmul(:,:,nsp),jd,kd)
        endif

        call VISRHS_BQ(turmu(:,:,nsp), turmub(:,:,nsp), q(:,:,nsp,:), qb(:,:,nsp,:),
     +              s(:,:,nsp,:), sb(:,:,nsp,:), 
     +              xx(:,:,nsp), xy(:,:,nsp), yx(:,:,nsp), yy(:,:,nsp),
     +              ug(:,:,nsp), vg(:,:,nsp), tscale(:,:,nsp), 
     +              iblank(:,:,nsp), jd, kd)

        if (iturb.eq.1) then
          call C_TURM_BQ(q(:,:,nsp,:), qb(:,:,nsp,:), turmu(:,:,nsp), 
     +                turmub(:,:,nsp), vmul(:,:,nsp), vmulb(:,:,nsp), jd, kd)
          call LAMVIS_BQ(q(:,:,nsp,:), qb(:,:,nsp,:), vmul(:,:,nsp), 
     +                vmulb(:,:,nsp), jd, kd)
          call vmu_sa_ad(q(:,:,nsp,:),qb(:,:,nsp,:),
     >                   s(:,:,nsp,:),sb(:,:,nsp,:),turmu(:,:,nsp),
     >                   x(:,:,nsp),y(:,:,nsp),xv(:,:,nsp),yv(:,:,nsp),
     >                   xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     >                   ug(:,:,nsp),vg(:,:,nsp),jd,kd,
     >                   tscale(:,:,nsp),iblank(:,:,nsp),im,0)
        endif
      endif

      if (timespectral) call timespectralrhs_bq(q,qb,s,sb,Ds,jd,kd,nsp)

      call BC_BQ(q(:,:,nsp,:), qb(:,:,nsp,:), x(:,:,nsp), y(:,:,nsp), 
     +        xx(:,:,nsp), xy(:,:,nsp), yx(:,:,nsp), yy(:,:,nsp), 
     +        ug(:,:,nsp), vg(:,:,nsp), im, jd, kd, bt(:,:,nsp))

      !am print*,'im, uinfb,vinfb, einfb:',im,uinfb,vinfb,einfb

      fsmachb = fsmachb - (cs*uinfb + fsmach*einfb + ss*vinfb)
      if (fmtip.eq.0 ) fsmachb = fsmachb + rey*reyb/fsmach
      ssb = fsmach*vinfb
      csb = fsmach*uinfb
      alfab = alfab - (pi*COS(pi*(alfa/180.0))*ssb/180.0 - pi*SIN(pi*(alfa/
     +  180.0))*csb/180.0)

      enddo spectralloop1

      return
      end

c***********************************************************************
      subroutine qdivj( q,jd,kd )
c
c  divide q by jacobian
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)

      ! local variables
      integer j,k

c***  first executable statement

      do 71 k = 1,kd
      do 71 j = 1,jd
        q(j,k,1) = q(j,k,1)/q(j,k,nq)
        q(j,k,2) = q(j,k,2)/q(j,k,nq)
        q(j,k,3) = q(j,k,3)/q(j,k,nq)
        q(j,k,4) = q(j,k,4)/q(j,k,nq)
   71 continue

      return
      end

c***********************************************************************
      subroutine qzero( q,jd,kd )
c
c  set initial values of the flow variables to their free-stream
c  values for an impulsive start. otherwise the initial values
c  are read from a restart file.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)

      ! local variables
      real ruinf,rvinf,rtkeinf,rtomegainf,tu,re_theta
      integer j,k

c**   first executable statement

      ruinf = rinf*uinf
      rvinf = rinf*vinf
      do 11 k=1,kd
      do 11 j=1,jd
        q(j,k,1) = rinf
        q(j,k,2) = ruinf
        q(j,k,3) = rvinf
        q(j,k,4) = einf
  11  continue

      if (iturb.eq.1) then
        do k=1,kd
        do j=1,jd
          q(j,k,5) = vnuinf
        enddo
        enddo
       
      elseif (iturb.eq.2) then
        rtkeinf = rinf*tkeinf
        rtomegainf = rinf*tomegainf
        do k=1,kd
        do j=1,jd
          q(j,k,5) = rtkeinf
          q(j,k,6) = rtomegainf
        enddo
        enddo

        if (itrans.eq.1) then
          itmcinf = 1. 

          tu = 100.*sqrt(2./3*tkeinf)/fsmach
          if(tu.le.1.3) then
            re_theta = (1173.51-589.428*tu+0.2196/(tu*tu))
          else
            re_theta = 331.5*(tu-0.5658)**(-0.671)
          endif
          retinf = max(re_theta,20.)/rey

          do k=1,kd
          do j=1,jd
            q(j,k,7) = rinf*itmcinf
            q(j,k,8) = rinf*retinf
          enddo
          enddo
        endif
      endif
c
      return
      end

c***********************************************************************
