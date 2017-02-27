c***********************************************************************
! lhs part
!
      subroutine step_ad_part2(q,qb,qb1,qtn,qtnm1,qnewt,s,sb,
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
      real qbs(jd,kd,nspec,nv) !asitav !comes from part1
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

      !amm qb(:,:,:,1:4) = 0. !comes from part1
      !amm qbs = 0. !comes from part1

!$OMP PARALLEL IF(NSPEC > 1)

!amm       if (.not.lamin.and.iturb.eq.1) then
!amm !$OMP DO ORDERED
!amm         spectralloop1: do nsp = 1,nspec
!amm           call vmu_sa_ad(q(:,:,nsp,:),qb(:,:,nsp,:),
!amm      >                   s(:,:,nsp,:),sb(:,:,nsp,:),turmu(:,:,nsp),
!amm      >                   x(:,:,nsp),y(:,:,nsp),xv(:,:,nsp),yv(:,:,nsp),
!amm      >                   xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
!amm      >                   ug(:,:,nsp),vg(:,:,nsp),jd,kd,
!amm      >                   tscale(:,:,nsp),iblank(:,:,nsp),im,1)
!amm           qb(:,:,nsp,5) = 0. 
!amm         enddo spectralloop1
!amm !$OMP END DO
!amm       endif
!amm 
!amm       call compute_residual_ad(q,qb,qb1,qtn,qtnm1,qnewt,qbs,s,sb,
!amm      &     x,y,xv,yv,iblank,xx,xy,yx,yy,ug,vg,ugv,vgv,
!amm      &     xold,yold,xole,yole,turmu,tscale,bt,Ds,im,jd,kd)

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
      end subroutine step_ad_part2

c***********************************************************************
