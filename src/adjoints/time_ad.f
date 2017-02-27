c**********************************************************************
      subroutine time_step(q,xx,xy,yx,yy,ug,vg,tscale,bt,iblank,jd,kd)
c
c     time step calculation
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer jd,kd
      real q(jd,kd,nspec,nq)
      real xx(jd,kd,nspec),xy(jd,kd,nspec),yx(jd,kd,nspec),yy(jd,kd,nspec)
      real ug(jd,kd,nspec),vg(jd,kd,nspec)
      real tscale(jd,kd,nspec),bt(jd,kd,nspec)
      integer iblank(jd,kd,nspec)

      ! local variables
      integer j,k,nsp
      real cflj,cflk,eigj,eigk,hnew,jac_factor
      real u,v,e,uu,vv,dxdx,dydy,asq,lam_j,lam_k
      real cflmaxj,cflmaxk,oat

c..set time-accuracy
      oat = 1.0
      if (ntac.ge.2 .and. istep.gt.1) oat = 2./3.
      if (ntac.eq.3 .and. istep.gt.2) oat = 6./11.

c..if time accurate preconditioning use dual time stepping

      if (iprecon.and.(timeac.eq.1.and..not.timespectral)) then
!        iconstant_cfl = 0
        idual = 1
        idual_turb = 1
      endif

      if (invisc) then
        jac_factor = 0.01
      else
        jac_factor = 0.002
      endif

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(j,k,cflj,cflk,hnew,u,v,e,uu,vv,dxdx,dydy,asq,lam_j,lam_k)
      spectralloop: do nsp = 1,nspec
        if(iconstant_cfl.eq.0) then

          if(idual.eq.0) then
            do k=1,kd
               do j=1,jd
                  tscale(j,k,nsp) = oat*h*max(iblank(j,k,nsp),0)*
     <                            ( 1.0 + jac_factor*(1.-timeac)*sqrt(q(j,k,nsp,nq)))/
     <                            ( 1.0 + (1.-timeac)*sqrt(q(j,k,nsp,nq)))
                  bt(j,k,nsp) = 1.
              enddo
            enddo
          else
            do k=1,kd
               do j=1,jd
                  tscale(j,k,nsp) = max(iblank(j,k,nsp),0)*
     <                            ( 1.0 + jac_factor*sqrt(q(j,k,nsp,nq)))/
     <                            ( 1.0 + sqrt(q(j,k,nsp,nq)))
                  if (iprecon.and.1.eq.1) then
                    tscale(j,k,nsp)  = 100.*max(iblank(j,k,nsp),0)*
     <                            ( 1.0 + 0.005*sqrt(q(j,k,nsp,nq)))/
     <                            ( 1.0 + sqrt(q(j,k,nsp,nq)))
!                    if (tscale(j,k,nsp).gt.40) tscale(j,k,nsp) = 40.
                  endif
                  tscale(j,k,nsp)=tscale(j,k,nsp)*dualtime
                  bt(j,k,nsp) = 1./(1.+tscale(j,k,nsp)/h/oat)
                  tscale(j,k,nsp)=tscale(j,k,nsp)/(1.+tscale(j,k,nsp)/h/oat)
               enddo
            enddo
          endif

        else

          do k=1,kd
            do j=1,jd
              u=q(j,k,nsp,2)/q(j,k,nsp,1)
              v=q(j,k,nsp,3)/q(j,k,nsp,1)
              e=q(j,k,nsp,4)/q(j,k,nsp,1)
              uu=(u-ug(j,k,nsp))*xx(j,k,nsp)+(v-vg(j,k,nsp))*xy(j,k,nsp)
              vv=(u-ug(j,k,nsp))*yx(j,k,nsp)+(v-vg(j,k,nsp))*yy(j,k,nsp)
              dxdx=xx(j,k,nsp)**2 + xy(j,k,nsp)**2
              dydy=yx(j,k,nsp)**2 + yy(j,k,nsp)**2
              asq=ggm1*(e-0.5*(u**2+v**2))
              lam_j=abs(uu)+sqrt(asq*dxdx)
              lam_k=abs(vv)+sqrt(asq*dydy)
              if (idual.eq.0) then
                tscale(j,k,nsp) = h/(lam_j+lam_k)
                tscale(j,k,nsp)  = tscale(j,k,nsp)*max(iblank(j,k,nsp),0)
                bt(j,k,nsp) = 1.
              else
                hnew=dualtime/(lam_j+lam_k)
                tscale(j,k,nsp)=hnew*max(iblank(j,k,nsp),0)           
                tscale(j,k,nsp)=tscale(j,k,nsp)/(1.+tscale(j,k,nsp)/h/oat)
                bt(j,k,nsp) = 1./(1.+tscale(j,k,nsp)/h/oat)
              endif
            enddo
          enddo

        endif
      enddo spectralloop
!$OMP END PARALLEL DO

      return
      end

c***********************************************************************
      subroutine stqol(q,qtn,qtnm1,qnewt,jd,kd)
c     
c  store data at previous time steps
c
c*********************************************************************** 
      use params_global
c*********************************************************************** 
      implicit none
c*********************************************************************** 
      integer jd,kd
      real q(jd,kd,nspec,nv),qtn(jd,kd,nspec,nv),qtnm1(jd,kd,nspec,nv)
      real qnewt(jd,kd,nspec,nv)

      ! local variables
      integer j,k,nsp,n
c**   first executable statement

!$OMP PARALLEL DO IF(NSPEC > 1)
!$OMP& PRIVATE(j,k,n)
      spectralloop: do nsp = 1,nspec
        do k = 1, kd

          do j = 1, jd
            do n = 1,nv
              qnewt(j,k,nsp,n) = q(j,k,nsp,n)
            enddo
          enddo

          if ( (ntac.eq.2 .and. istep.gt.1) .or. 
     <         (ntac.eq.3 .and. istep.eq.2) ) then
            do j = 1, jd
              do n = 1,nv
                qnewt(j,k,nsp,n) = qnewt(j,k,nsp,n) +
     <                            (q(j,k,nsp,n)-qtn(j,k,nsp,n))/3.
              enddo
            enddo
          endif

          if (ntac.eq.3 .and. istep.gt.2) then
            do j = 1, jd 
              do n = 1,nv
                qnewt(j,k,nsp,n) = qnewt(j,k,nsp,n) + 
     <                         7.*(q(j,k,nsp,n)-qtn(j,k,nsp,n))/11.-
     <                         2.*(qtn(j,k,nsp,n)-qtnm1(j,k,nsp,n))/11.
              enddo
            enddo
          endif

          do j = 1, jd
            do n = 1,nv
              qtnm1(j,k,nsp,n) = qtn(j,k,nsp,n)
            enddo
          enddo

          do j = 1, jd
            do n = 1,nv
              qtn(j,k,nsp,n) = q(j,k,nsp,n)
            enddo
          enddo

        enddo
      enddo spectralloop
!$OMP END PARALLEL DO
c
      return
      end

c**********************************************************************
      subroutine newton(q,qtn,qtnm1,qnewt,s,jd,kd)
c
c  take into account time terms on rhs for newton iterations
c
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer jd,kd
      real q(jd,kd,nq),qtn(jd,kd,nv),qtnm1(jd,kd,nv)
      real qnewt(jd,kd,nv),s(jd,kd,nv)

      ! local variables
      integer j,k,n
      real oat,tac

c..set time-accuracy
      oat = 1.0
      if (ntac.ge.2 .and. istep.gt.1) oat = 2./3.
      if (ntac.eq.3 .and. istep.gt.2) oat = 6./11.

      tac = 1./(oat*h) 
      do k = kbeg, kend
      do j = jbeg, jend
        do n = 1,nv
          s(j,k,n) = s(j,k,n) - (q(j,k,n) - qnewt(j,k,n))*tac
        enddo
      enddo
      enddo
c     
      return
      end

c**********************************************************************
      subroutine newtonTS(q,qnewt,s,tscale,jd,kd)
c
c  take into account time terms on rhs for newton iterations
c
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer jd,kd
      real q(jd,kd,nv),qnewt(jd,kd,nv),s(jd,kd,nv)
      real tscale(jd,kd)

      ! local variables
      integer j,k,n
      real oat,tac

      do k = kbeg, kend
      do j = jbeg, jend
        if (tscale(j,k).ne.0) then
          do n = 1,nv
            s(j,k,n) = s(j,k,n) - (q(j,k,n) - qnewt(j,k,n))/tscale(j,k)
          enddo
        endif
      enddo
      enddo
c     
      return
      end

c***********************************************************************
