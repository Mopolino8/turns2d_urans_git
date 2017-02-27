c***********************************************************************
      subroutine step(q,qtn,qtnm1,qnewt,s,
     &     x,y,xv,yv,iblank,
     &     xx,xy,yx,yy,ug,vg,ugv,vgv,
     &     xold,yold,xole,yole,
     &     turmu,tscale,bt,Ds,
     &     im,jd,kd)

c  note that the boundaries are updated explicitly
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer im,jd,kd
      real q(jd,kd,nspec,nq), s(jd,kd,nspec,nv)
      real tscale(jd,kd,nspec),bt(jd,kd,nspec)
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
      real resmax,resrho,rsum
      real smrs,smrs1,smrs2,smrs3,smrs4,volum
      real,allocatable :: ss(:,:,:,:)

c***  first executable statement
      allocate(ss(jd,kd,nspec,nv))

c..zero s array 
      s = 0.
      ss = 0.

!$OMP PARALLEL IF(NSPEC > 1)

      call compute_residual(q,qtn,qtnm1,qnewt,s,ss,
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
              smrs = abs(s(j,k,nsp,n))*max(iblank(j,k,nsp),0)
              if (smrs .gt. resmax) then
                jd1 = j
                kd1 = k
                nd1 = n
                resmax = smrs
              endif
            enddo
          enddo

          do j = jbeg,jend
            smrs1 = s(j,k,nsp,1)*max(iblank(j,k,nsp),0)
            smrs2 = s(j,k,nsp,2)*max(iblank(j,k,nsp),0)
            smrs3 = s(j,k,nsp,3)*max(iblank(j,k,nsp),0)
            smrs4 = s(j,k,nsp,4)*max(iblank(j,k,nsp),0)
            resrho  = resrho + s(j,k,nsp,1)**2
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
           write(70+im,101) istep0,rsum,resrho,resmax,totime,theta_col
  101      format(i7,5(e18.10))
           write(6,102) jd1,kd1,nd1,resmax,resrho,rsum
  102      format('  j,k,n,rmax,l2rho,l2 =',3i4,3(x,e12.4))
!$OMP END ORDERED
        endif

c..stop the code if the l2norm is gone out of bounds
        if( rsum .gt.1000.0 ) then
          write(6,602)  rsum
  602     format(' ',10x,'norm is out of bounds,'
     $           ,1x,'l2norm = ',f16.12,1x,'solution suspended' )
          stop 'norm'
        end if

      enddo spectralloop2
!$OMP END DO
 
      endif tscheck

!$OMP DO
!$OMP& PRIVATE(j,k,n)
      spectralloop3: do nsp = 1,nspec

c..multiply by low Mach preconditioner matrix
      if (iprecon)  call rhslom(q(:,:,nsp,:),s(:,:,nsp,:),jd,kd,
     &                          jbeg,jend,kbeg,kend,bt(:,:,nsp))

c..add newton correction in time-spectral simulation
      do k = kbeg,kend
        do j = jbeg,jend
          do n = 1,nv
           s(j,k,nsp,n) = s(j,k,nsp,n) + ss(j,k,nsp,n)
          enddo
        enddo
      enddo

c..update residual bc 
      call rbc(q(:,:,nsp,:),s(:,:,nsp,:),x(:,:,nsp),y(:,:,nsp),
     &         xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &         ug(:,:,nsp),vg(:,:,nsp),im,jd,kd,bt(:,:,nsp))

c..multiply by time-step, loop through all points in case residuals
c..of halo points are not zero because of residual bc application
      do k = 1,kd
        do j = 1,jd
          s(j,k,nsp,1) = s(j,k,nsp,1)*tscale(j,k,nsp)
          s(j,k,nsp,2) = s(j,k,nsp,2)*tscale(j,k,nsp)
          s(j,k,nsp,3) = s(j,k,nsp,3)*tscale(j,k,nsp)
          s(j,k,nsp,4) = s(j,k,nsp,4)*tscale(j,k,nsp)
        enddo
      enddo
c
c..now do the implicit part
c  
      if(ilhs.eq.1) then

         if(.not.iprecon) then
           call ilu2d(q(:,:,nsp,:),s(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                tscale(:,:,nsp))
         else
           call preilu2d(q(:,:,nsp,:),s(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                tscale(:,:,nsp),bt(:,:,nsp))
         endif  

      elseif(ilhs.eq.2.or.ilhs.eq.3) then
         if(.not.iprecon) then 
            call arc2d(q(:,:,nsp,:),s(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                iblank(:,:,nsp),tscale(:,:,nsp),bt(:,:,nsp))
         else
            call arc2d_precon(q(:,:,nsp,:),s(:,:,nsp,:),
     &                jd,kd,jbeg,jend,kbeg,kend,
     &                xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &                ug(:,:,nsp),vg(:,:,nsp),turmu(:,:,nsp),
     &                iblank(:,:,nsp),tscale(:,:,nsp),bt(:,:,nsp))
         endif
      endif
 
      enddo spectralloop3
!$OMP END DO

c....solve turbulence equations 
!$OMP DO ORDERED
      spectralloop4: do nsp = 1,nspec
        if (.not.lamin)
     >    call solve_turbulence_equations(turmu(:,:,nsp),q(:,:,nsp,:),
     >           s(:,:,nsp,:),x(:,:,nsp),y(:,:,nsp),xv(:,:,nsp),yv(:,:,nsp),
     >           xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     >           ug(:,:,nsp),vg(:,:,nsp),
     >           tscale(:,:,nsp),iblank(:,:,nsp),jd,kd,im,nsp)

      enddo spectralloop4
!$OMP END DO
c
c..update q with corrections
c
!$OMP DO
      spectralloop5: do nsp = 1,nspec

      do 31 k = kbeg,kend
      do 31 j = jbeg,jend
        q(j,k,nsp,1) = q(j,k,nsp,1) + s(j,k,nsp,1)*max(iblank(j,k,nsp),0)
        q(j,k,nsp,2) = q(j,k,nsp,2) + s(j,k,nsp,2)*max(iblank(j,k,nsp),0)
        q(j,k,nsp,3) = q(j,k,nsp,3) + s(j,k,nsp,3)*max(iblank(j,k,nsp),0)
        q(j,k,nsp,4) = q(j,k,nsp,4) + s(j,k,nsp,4)*max(iblank(j,k,nsp),0)
   31 continue

      enddo spectralloop5
!$OMP END DO
!$OMP END PARALLEL

      return
      end

c***********************************************************************
      subroutine compute_residual(q,qtn,qtnm1,qnewt,s,ss,
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
      real q(jd,kd,nspec,nq), s(jd,kd,nspec,nv)
      real ss(jd,kd,nspec,nv)
      real tscale(jd,kd,nspec),bt(jd,kd,nspec)
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

c***  first executable statement

!$OMP DO ORDERED
      spectralloop: do nsp = 1,nspec
c     
c..update b c
c
      call bc(q(:,:,nsp,:),x(:,:,nsp),y(:,:,nsp),
     &        xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &        ug(:,:,nsp),vg(:,:,nsp),im,jd,kd,bt(:,:,nsp))

c..compute the right hand side and store in s array
c
      call rhsup(q(:,:,nsp,:),s(:,:,nsp,:),
     &           xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     &           x(:,:,nsp),y(:,:,nsp),xv(:,:,nsp),yv(:,:,nsp),
     &           xold(:,:,nsp),yold(:,:,nsp),xole(:,:,nsp),yole(:,:,nsp),
     &           iblank(:,:,nsp),ugv(:,:,nsp),vgv(:,:,nsp),
     &           jd,kd,1,im,bt(:,:,nsp))

      if (iunst.eq.2) call momsource(q(:,:,nsp,:),s(:,:,nsp,:),jd,kd)
      if (axisymmetric) call axisource(q(:,:,nsp,:),s(:,:,nsp,:),x(:,:,nsp),y(:,:,nsp),jd,kd)
c
c..start newton iteration here at each time step for convergence
c
      if (itnmax.gt.1) then
        if (timespectral) then
          call newtonTS(q(:,:,nsp,:),qnewt(:,:,nsp,:),ss(:,:,nsp,:),
     &                  tscale(:,:,nsp),jd,kd)
        else
          call newton(q(:,:,nsp,:),qtn(:,:,nsp,:),
     &            qtnm1(:,:,nsp,:),qnewt(:,:,nsp,:),s(:,:,nsp,:),jd,kd)
        endif
      endif
c
c..compute viscous fluxes
      if( .not. invisc ) then

        call visrhs(turmu(:,:,nsp),q(:,:,nsp,:),s(:,:,nsp,:),
     >              xx(:,:,nsp),xy(:,:,nsp),yx(:,:,nsp),yy(:,:,nsp),
     >              ug(:,:,nsp),vg(:,:,nsp),tscale(:,:,nsp),
     >              iblank(:,:,nsp),jd,kd)
      endif

      if (timespectral) call timespectralrhs(q,s,Ds,jd,kd,nsp)

      enddo spectralloop
!$OMP END DO

      end subroutine compute_residual

c***********************************************************************
      subroutine solve_turbulence_equations(turmu,q,s,x,y,xv,yv,
     &             xx,xy,yx,yy,ug,vg,tscale,iblank,jd,kd,im,nsp)
c
c  solve turbulence equations 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,im,nsp
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd), tscale(jd,kd)
      real x(jd,kd),y(jd,kd),xv(jmax,kmax),yv(jmax,kmax),ug(jd,kd),vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      integer iblank(jd,kd)

      if (iturb.eq.0) then
        if (bodyflag(im)) call vmutur( x,y,q,s,turmu,xx,xy,yx,yy,ug,vg,jd,kd)
      elseif (iturb.eq.1) then
        call vmu_sa( q,s,turmu,x,y,xv,yv,xx,xy,yx,yy,
     <               ug,vg,jd,kd,tscale,iblank,im,nsp)
      elseif (iturb.eq.2) then
        call komegasst( q,s,turmu,x,y,xv,yv,xx,xy,yx,yy,
     <                  ug,vg,jd,kd,tscale,iblank,im)
      endif

      end subroutine solve_turbulence_equations      

c*************************************************************************
      subroutine monitor(mstop,q,jd,kd)
c
c  this subroutine checks for negative speed of sound
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nspec,nq)

      ! local variables
      integer mstop,nneg,k,j,nsp
      real asqmin,rho,qq,asq

c*** first executable statement

c
c..negative speed of sound check
c
      mstop = 0
      asqmin = 0.002
      nneg = 0
!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(j,k,rho,qq,asq)
      spectralloop: do nsp = 1,nspec
       do 910 k = 1, kd
        do 920 j = 1, jd
          rho = q(j,k,nsp,1)
          qq = q(j,k,nsp,2)**2 +q(j,k,nsp,3)**2
          asq = ggm1*(q(j,k,nsp,4) - .5*qq/rho)/rho
ccray          nneg = cvmgt(nneg+1,nneg,asq.le.asqmin)
c..asq+1 == asq checks for NaN's (Vinod)  
          if( asq.le.asqmin .or. asq+1.eq.asq ) then
            nneg = nneg+1
          end if
 920    continue
 910   continue
c
      if(nneg .ne. 0) then
        write(6,*) nneg, ' negative speed of sound'
        do 74 k = 1,kd
        do 74 j = 1,jd
          rho = q(j,k,nsp,1)
          qq = q(j,k,nsp,2)**2 +q(j,k,nsp,3)**2
          asq = ggm1*(q(j,k,nsp,4) - .5*qq/rho)/rho
          if( asq.le.asqmin .or. asq+1.eq.asq ) then
            mstop = 1
            write(6,601) j,k,nsp,asq
!            return
          end if
   74   continue
      end if
      enddo spectralloop
!$OMP END PARALLEL DO
c
 601  format(' j,k,asq from monitor ',3i5,f12.5)
c
      return
      end

c***********************************************************************
