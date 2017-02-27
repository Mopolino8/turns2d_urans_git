c***********************************************************************
      subroutine step_ad(q,qb,qb1,qtn,qtnm1,qnewt,s,sb,
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
      real,allocatable :: qbs(:,:,:,:)

c***  first executable statement
      allocate(qbs(jd,kd,nspec,nv))

      qb(:,:,:,1:4) = 0.
      qbs = 0.

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

c
c..check on convergence
c..if timespectral, check only for first sub-iteration 
c
      tscheck: if (.not.timespectral.or.itn.eq.1) then

!$OMP DO ORDERED
!$OMP& PRIVATE(j,k,n,jd1,kd1,nd1)
!$OMP& PRIVATE(rsum,resrho,resmax,smrs,smrs1,smrs2,smrs3,smrs4,volum)
      spectralloop2: do nsp = 1,nspec

        rsum = 0.
        resrho  = 0.0
        resmax  = 0.0
        jd1=2
        kd1=2
        nd1=1

        do k = kbeg,kend
          do n = 1,nmv
            do j = jbeg,jend
              smrs = abs(qb(j,k,nsp,n))*max(iblank(j,k,nsp),0)
              if (smrs .gt. resmax) then
                jd1 = j
                kd1 = k
                nd1 = n
                resmax = smrs
              endif
            enddo
          enddo

          do j = jbeg,jend
            smrs1 = qb(j,k,nsp,1)*max(iblank(j,k,nsp),0)
            smrs2 = qb(j,k,nsp,2)*max(iblank(j,k,nsp),0)
            smrs3 = qb(j,k,nsp,3)*max(iblank(j,k,nsp),0)
            smrs4 = qb(j,k,nsp,4)*max(iblank(j,k,nsp),0)
            resrho  = resrho + qb(j,k,nsp,1)**2
            rsum = rsum + smrs1*smrs1 + smrs2*smrs2 + smrs3*smrs3
     <                  + smrs4*smrs4
          enddo
        enddo
c
        volum  = float(jd*kd)
        resrho = sqrt(resrho/volum)
        rsum   = sqrt(rsum/volum)
c
        if( mod(istep0,npnorm).eq.0 ) then
!$OMP ORDERED
           write(170+im,101) istep0,rsum,resrho,resmax,totime,theta_col
  101      format(i7,5(e18.10))
           write(6,102) jd1,kd1,nd1,resmax,resrho,rsum
  102      format('  j,k,n,rmax,l2rho,l2 =',3i4,3(x,e12.4))
!$OMP END ORDERED
        endif

c..stop the code if the l2norm is gone out of bounds
        if( rsum .gt.1.0e14 ) then
          write(6,602)  rsum
  602     format(' ',10x,'norm is out of bounds,'
     $           ,1x,'l2norm = ',e18.10,1x,'solution suspended' )
          stop 'norm'
        end if

      enddo spectralloop2
!$OMP END DO
 
      endif tscheck

!$OMP DO
!$OMP& PRIVATE(j,k,n)
      spectralloop3: do nsp = 1,nspec

!c..multiply by low Mach preconditioner matrix
!      if (iprecon)  call rhslom_ad(q(:,:,nsp,:),qb(:,:,nsp,:),jd,kd,
!     &                          jbeg,jend,kbeg,kend,bt(:,:,nsp))

c..add newton correction in time-spectral simulation
      do k = kbeg,kend
        do j = jbeg,jend
          do n = 1,nv
            qb(j,k,nsp,n) = qb(j,k,nsp,n) + qbs(j,k,nsp,n)
          enddo
        enddo
      enddo

c..multiply by time-step, loop through all points in case residuals
c..of halo points are not zero because of residual bc application
      do k = 1,kd
        do j = 1,jd
          qb(j,k,nsp,1) = qb(j,k,nsp,1)*tscale(j,k,nsp)
          qb(j,k,nsp,2) = qb(j,k,nsp,2)*tscale(j,k,nsp)
          qb(j,k,nsp,3) = qb(j,k,nsp,3)*tscale(j,k,nsp)
          qb(j,k,nsp,4) = qb(j,k,nsp,4)*tscale(j,k,nsp)
        enddo
      enddo
c
c..now do the implicit part
c  
      if(ilhs.eq.1) then

         if(.not.iprecon) then
           call ilu2d(q(:,:,nsp,:),qb(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                tscale(:,:,nsp))
         else
           call preilu2d(q(:,:,nsp,:),qb(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                tscale(:,:,nsp),bt(:,:,nsp))
         endif  

      elseif(ilhs.eq.2.or.ilhs.eq.3) then
         if(.not.iprecon) then 
            call arc2d(q(:,:,nsp,:),qb(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                iblank(:,:,nsp),tscale(:,:,nsp),bt(:,:,nsp))
         else
            call arc2d_precon(q(:,:,nsp,:),qb(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                iblank(:,:,nsp),tscale(:,:,nsp),bt(:,:,nsp))
         endif
      endif
 
      enddo spectralloop3
!$OMP END DO
c
c..update q with corrections
c
!$OMP DO
!$OMP& PRIVATE(j,k)
      spectralloop4: do nsp = 1,nspec

      do 31 k = kbeg,kend
      do 31 j = jbeg,jend
        sb(j,k,nsp,1) = sb(j,k,nsp,1) + qb(j,k,nsp,1)*max(iblank(j,k,nsp),0)
        sb(j,k,nsp,2) = sb(j,k,nsp,2) + qb(j,k,nsp,2)*max(iblank(j,k,nsp),0)
        sb(j,k,nsp,3) = sb(j,k,nsp,3) + qb(j,k,nsp,3)*max(iblank(j,k,nsp),0)
        sb(j,k,nsp,4) = sb(j,k,nsp,4) + qb(j,k,nsp,4)*max(iblank(j,k,nsp),0)
   31 continue

!      call debugstore(q,qb,jd,kd,101)

c..update b c
      call bc_ad(sb(:,:,nsp,:),im,jd,kd)

      enddo spectralloop4
!$OMP END DO
!$OMP END PARALLEL

      return
      end

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
