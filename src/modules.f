c*******************************************************************
c..   global parameters module
c*******************************************************************
      module params_global

c..   pi

      real :: pi

c  parameters that define critical points for c-h mesh


      integer half,jle,jtail1,jtail2
      
c..parameters for space and time
c..dx1,dy1,dz1 always set to 1
c..h,hd calculated based on input in initia

      real :: dx1, dy1, dz1, h, hd !dels

c  parameters that define the flow
c..fsmach,rey set by input
c..gamma,gm1,ggm1,pr,tinf set to fixed values in initia

      real ::  gm1, gamma, ggm1, fsmach, fmtip, rartio, pr, rey, tinf


c  parameters that define the freestream values, set in initia

      real ::  einf, htinf, pinf, rinf, uinf, vinf, ainf, vnuinf, 
     c         tkeinf, tomegainf, itmcinf, retinf, tuinf

c..parameters for connectivity
      logical :: static_conn

c  parameters that define grid spatial dimensions and time steps
c..jmax,kmax,nq,nv,nsteps set by input
c..istep0 is initial time step read in with q-file
c..jm,km calculated based on input in initia

      integer :: jmax,kmax,jm,km,jbeg,kbeg,jend,kend
      integer :: nq,nv,nmv,nturb
      integer :: istep0,nsteps

c  parameters that define time spectral calculation
c..timespectral defines whether calculation is spectral in time
c..nspec defines number of time spectral instances
      logical :: timespectral
      integer :: nspec

c  parameters for type of limiting and pressure bc
c..ilim set by input
      integer ::   ilim
c..parameters for LHS
      integer ::   ilhs,idual,idual_turb,iconstant_cfl
      real    ::   dtpseudo(100),cfltur,dualtime
c..parameters for preconditioning
      logical ::   iprecon
      real    ::   Mp
c  parameters used for grid coarsening
c..jint,kint set by input
      integer :: jint,kint

c  parameters used for spatial order of accuracy
c..irhsy set by input

      integer :: irhsy

c  parameters for describing rotor conditions
c..totime set by input
c..rf set by input
c..angmax set by input
c..dang is current angle of attack difference

      real ::  totime,rf,ang0,angmax,phase,dang,rotf,fourbar_toggle
      real ::  xac,yac

c parameter for describing axisymmetric flow

      logical :: axisymmetric

c  parameters for restart, storing out restart and writing residual
c..iread,nrest and npnorm set by input

      integer:: iread, nrest, npnorm

c  parameters for damping left-hand-side
c..epse set by input

      real :: epse

c  parameters for time 
c..dt,timeac,iunst,ntac,itnmax set by input
c..istep is current time step in the run

      real    :: dt, timeac
      integer :: iunst, ntac, itnmax, istep, itn
      integer :: motiontype

c  parameters for viscous flow
c..invisc,lamin set by input
c..iturb set by input
c..rmue set in initia
c..nturiter set in vmu_sa

      real    :: rmue
      logical :: invisc, lamin
      integer :: iturb, itrans
      integer :: nturiter

C..Parameter for boundary conditions
      real    :: jtail(6)
      integer :: irot_dir_all(6)
      integer :: nbc_all(6),ibtyp_all(25,6),ibdir_all(25,6)
      integer :: jbcs_all(25,6),kbcs_all(25,6)
      integer :: jbce_all(25,6),kbce_all(25,6)
      logical :: bodyflag(6)
      real    :: xac_all(6),yac_all(6)

      real    :: thetan,thetao
      integer :: nmovie,isin
      integer :: num_grids,nblade
      
      real    :: theta_col

C..Asitav TEF variables
C..--------------------
      integer :: iteflap,nharmflap
      real    :: pxc,dela,rf_tef,theta_f0,theta_finit
      real    :: ampFlap(5),phiFlap(5)
C..--------------------

      ! dimension values
      integer :: jdim,kdim,mdim,isdim,iqdim,nhalo,nadd
      integer :: igrd,igrdv,igrdb

C.. Anand
      ! Optimization/inverse framework stuff
      logical :: if_obj
      ! obj type codes
      integer :: obj_ftype_cp, obj_ftype_cf, obj_ftype_u, obj_ftype_v, obj_ftype_nu
      data obj_ftype_cp,obj_ftype_cf, obj_ftype_u, obj_ftype_v,obj_ftype_nu /1990,1991,1992,1993,1994/
      integer :: obj_ftype_cltot,obj_ftype_cdtot,obj_ftype_cmtot
      real, allocatable :: clobj(:),cdobj(:),cmobj(:)
      logical :: if_objtot !when objective is integrated forces (Cltot,Cdtot,Cmtot)
      data obj_ftype_cltot,obj_ftype_cdtot,obj_ftype_cmtot /2990,2991,2992/

      data if_obj/.false./
      data if_objtot/.false./


      logical :: obj_pointwise
      integer :: obj_points(1000)
      integer :: obj_ftype, obj_npoints
      integer :: obj_normtype
      real, allocatable :: obj_f(:), obj_coeff_prod(:,:,:) !now defined in Turns
      logical :: obj_samod
      real, allocatable :: obj_cpx(:),obj_cpy(:)
      
      logical :: obj_if_reg
      real :: obj_reg_fac
      data obj_pointwise,obj_ftype,obj_npoints,obj_normtype/.false.,1990,-1,2/
      data obj_if_reg, obj_reg_fac/.false.,0.0/


C..--------


      
      data dx1,dy1 / 2*1.0 /
      data gamma,pr,rmue,tinf /1.4,0.72,1.0,400./
      data iread /0/
      data num_grids,nblade /1,1/
      data jmax,kmax,jtail1,nhalo,nadd,half /109,36,19,1,1,0/
      data nq,nv,nmv,nturb /5,4,4,0/
      data nsteps,nrest,npnorm /500,1000,25/
      data timespectral,nspec /.false.,1/
      data fsmach,fmtip,alfa,rartio,rey,invisc,lamin,iturb, itrans
     <           /0.16,0.0,0.0,1.,3900000.,.true.,.true.,1,0/
      data vnuinf,tuinf /0.1,0.2/
      data iunst,ntac,itnmax,dt,timeac,motiontype /0,1,1,.1,1.,1/
      data epse,irhsy,ilhs,ilim,totime,rf,ang0,angmax,phase /0.01,-3,1,0,0.,0.,0.,0.,0./
      data rotf,fourbar_toggle,xac,yac /0.,0,0.,0./
      data axisymmetric /.false./
      data jint,kint /1,1/
      data static_conn /.false./
      data nmovie,isin /50,0/
      data ideform /0/
      data idual,idual_turb,dualtime,cfltur,iconstant_cfl/0,0,1.0,1.0,0/
      data iprecon,Mp/.false.,1.0/ 
!asitav
      data iteflap /0/


      

      namelist/inputs/ iread,num_grids,nblade,
     & jmax,kmax,nhalo,jtail1,half,nsteps,nrest,npnorm,
     & timespectral,nspec,
     & fsmach,alfa,rey,invisc,lamin,iturb,itrans,iunst,motiontype,
     & fsmach,fmtip,rartio,alfa,rey,invisc,lamin,iturb,itrans,iunst,
     & vnuinf,tuinf,ntac,itnmax,dt,dtpseudo,cfltur,timeac,
     & epse,irhsy,ilim,ilhs,idual,idual_turb,iconstant_cfl,iprecon,Mp,
     & totime,rf,ang0,angmax,phase,jint,kint,static_conn,
     & fourbar_toggle,axisymmetric,
     & nmovie,isin,ideform,
     & iteflap,
     & if_obj,if_objtot




      end module params_global

c**************************************************************************
