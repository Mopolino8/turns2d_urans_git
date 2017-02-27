c***********************************************************************
      subroutine timespectralrhs(q,s,Ds,jd,kd,nsp)
c..compute time-spectral coupling term
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,nsp
      real q(jd,kd,nspec,nq),s(jd,kd,nspec,nv),Ds(nspec,nspec)

      ! local variables
      integer :: j,k,n,nsp1

      do n = 1,nv
        do nsp1 = 1,nspec
          do k = kbeg,kend 
            do j = jbeg,jend
              s(j,k,nsp,n) = s(j,k,nsp,n) - rotf*Ds(nsp,nsp1)*q(j,k,nsp1,n)
            enddo
          enddo
        enddo
      enddo

      end subroutine timespectralrhs

c***********************************************************************
      subroutine momsource(q,s,jd,kd)
c
c  add source term to the rhs
c
c*************************************************************************

      use params_global
      implicit none

      integer jd,kd
      real q(jd,kd,nq),s(jd,kd,nv)

      integer j,k

C$AD II-LOOP
      do j=jbeg,jend
C$AD II-LOOP
      do k=kbeg,kend
        s(j,k,2) = s(j,k,2) + rotf*q(j,k,3)
        s(j,k,3) = s(j,k,3) - rotf*q(j,k,2)
      enddo
      enddo

      return
      end

c***********************************************************************
      subroutine axisource(q,s,x,y,jd,kd)
c
c  add axisymmetric source term to the rhs
c
c*************************************************************************

      use params_global
      implicit none

      integer jd,kd
      real q(jd,kd,nq),s(jd,kd,nv)
      real x(jd,kd),y(jd,kd)

      integer j,k
      real ylim,ppp,fac

      ylim = 0.1
C$AD II-LOOP
      do j=jbeg,jend
C$AD II-LOOP
      do k=kbeg,kend

        if ( y(j,k).gt.ylim) then
          fac = 1./y(j,k)
        else
          fac = 1./ylim
        endif

        s(j,k,1) = s(j,k,1) - fac*q(j,k,3)
        s(j,k,2) = s(j,k,2) - fac*q(j,k,2)*q(j,k,3)/q(j,k,1)
        s(j,k,3) = s(j,k,3) - fac*q(j,k,3)*q(j,k,3)/q(j,k,1)
        ppp = gm1*( q(j,k,4)-0.5*(q(j,k,2)**2+q(j,k,3)**2)/q(j,k,1) )
        s(j,k,4) = s(j,k,4) - fac*(q(j,k,4) + ppp)*q(j,k,3)/q(j,k,1)
      enddo
      enddo

      return
      end

c***********************************************************************
      subroutine rhslom(q,s,jd,kd,js,je,ks,ke,bt)
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c*********************************************************************c

         real q(jd,kd,nq),s(jd,kd,nv),bt(jd,kd)
         integer jd,kd,js,je,ks,ke

      ! local variables
         real :: tmp(4)

         integer j,k,n
         real u,v,ge,aSq,phiSq,bSq


c***  first executable statement

C$AD II-LOOP
         do j = js,je
C$AD II-LOOP
            do k = ks,ke
               u = q(j,k,2)/q(j,k,1)
               v = q(j,k,3)/q(j,k,1)
               phiSq = 0.5*( u*u + v*v )
               ge    = q(j,k,4)/q(j,k,1)*gamma-phiSq*gm1
               aSq   = (q(j,k,4)/q(j,k,1)-phiSq)*gm1*gamma
C$AD II-LOOP
               do n = 1,4 
                 tmp(n)=s(j,k,n)
               enddo
               s(j,k,1) = phiSq*tmp(1)-u*tmp(2)-v*tmp(3)+tmp(4)
               s(j,k,2) = u*s(j,k,1)
               s(j,k,3) = v*s(j,k,1)
               s(j,k,4) = ge*s(j,k,1)
      
               bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

C$AD II-LOOP
               do n = 1,4
                  s(j,k,n) = s(j,k,n)*(gm1*(bSq-1.)/aSq)
                  s(j,k,n) = s(j,k,n) + tmp(n)
               enddo

           enddo
        enddo

      return
      end

c************************************************************************

