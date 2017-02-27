c***********************************************************************
! rhs part
      subroutine step_ad_part1(q,qb,qb1,qtn,qtnm1,qnewt,s,sb,
     &     qbs, !asitav
     &     x,y,xv,yv,iblank,
     &     xx,xy,yx,yy,ug,vg,ugv,vgv,
     &     xold,yold,xole,yole,
     &     turmu,tscale,bt,Ds,
     &     im,jd,kd,resmax,resrho,rsum)

c  note that the boundaries are updated explicitly
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer im,jd,kd
      real resmax,resrho,rsum
      real q(jd,kd,nspec,nq), qb(jd,kd,nspec,nq), qb1(jd,kd,nspec,nq)
      real s(jd,kd,nspec,nv), sb(jd,kd,nspec,nv)
      real qbs(jd,kd,nspec,nv) !asitav !goes into part2
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

      integer k,j,nsp,n,jd1,kd1,nd1
      real smrs,smrs1,smrs2,smrs3,smrs4,volum
      !amm real,allocatable :: qbs(:,:,:,:)

c***  first executable statement
      !amm allocate(qbs(jd,kd,nspec,nv))

      qb(:,:,:,1:4) = 0. !in part 1, carries forth to part2
      qbs = 0.  !in part 1

!$OMP PARALLEL IF(NSPEC > 1)

      if (.not.lamin.and.iturb.eq.1) then
!$OMP DO ORDERED
        spectralloop1: do nsp = 1,nspec
          call vmu_sa_ad(q(:,:,nsp,:),qb(:,:,nsp,:),
     >                   s(:,:,nsp,:),sb(:,:,nsp,:),turmu(:,:,nsp),
     >                   x(:,:,nsp),y(:,:,nsp),xv(:,:,nsp),yv(:,:,nsp),
     >                   xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     >                   ug(:,:,nsp),vg(:,:,nsp),jd,kd,
     >                   tscale(:,:,nsp),iblank(:,:,nsp),im,1)
          qb(:,:,nsp,5) = 0. 
        enddo spectralloop1
!$OMP END DO
      endif

      call compute_residual_ad(q,qb,qb1,qtn,qtnm1,qnewt,qbs,s,sb,
     &     x,y,xv,yv,iblank,xx,xy,yx,yy,ug,vg,ugv,vgv,
     &     xold,yold,xole,yole,turmu,tscale,bt,Ds,im,jd,kd)


!break here 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! do_interpolations_bq outside
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!

!$OMP END PARALLEL

      return
      end subroutine step_ad_part1

c***********************************************************************
      subroutine compute_residual_ad(q,qb,qb1,qtn,qtnm1,qnewt,qbs,s,sb,
     &     x,y,xv,yv,iblank,
     &     xx,xy,yx,yy,ug,vg,ugv,vgv,
     &     xold,yold,xole,yole,
     &     turmu,tscale,bt,Ds,
     &     im,jd,kd)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer im,jd,kd
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
      real qnewt(jd,kd,nspec,nv),qbs(jd,kd,nspec,nv)
      integer iblank(jd,kd,nspec)
      real Ds(nspec,nspec)

      integer j,k,nsp,n
      real,allocatable :: turmub(:,:,:),vmul(:,:,:),vmulb(:,:,:)

      allocate(turmub(jd,kd,nspec),vmul(jd,kd,nspec),vmulb(jd,kd,nspec))

c***  first executable statement

      turmub = 0.
      vmulb = 0.0

!$OMP DO ORDERED
!$OMP& PRIVATE(j,k,n)
      spectralloop1: do nsp = 1,nspec
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
c..start newton iteration here at each time step for convergence
c
      if (itnmax.gt.1) then
        if (timespectral) then
          call newtonTS(sb(:,:,nsp,:),qnewt(:,:,nsp,:),qbs(:,:,nsp,:),
     &                  tscale(:,:,nsp),jd,kd)
!        else
!          call newton(q(:,:,nsp,:),qtn(:,:,nsp,:),
!     &            qtnm1(:,:,nsp,:),qnewt(:,:,nsp,:),s(:,:,nsp,:),jd,kd)
        endif
      endif
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
        endif
      endif

      if (timespectral) call timespectralrhs_bq(q,qb,s,sb,Ds,jd,kd,nsp)

      do j = 1,jd
      do k = 1,kd
        do n = 1,4
          qb(j,k,nsp,n) = qb(j,k,nsp,n) - qb1(j,k,nsp,n)
        enddo
      enddo
      enddo
c     
      call BC_BQ(q(:,:,nsp,:), qb(:,:,nsp,:), x(:,:,nsp), y(:,:,nsp), 
     +        xx(:,:,nsp), xy(:,:,nsp), yx(:,:,nsp), yy(:,:,nsp), 
     +        ug(:,:,nsp), vg(:,:,nsp), im, jd, kd, bt(:,:,nsp))

      enddo spectralloop1
!$OMP END DO
!      call debugstore(q,qb,jd,kd,100) 

      end subroutine compute_residual_ad

c***********************************************************************
