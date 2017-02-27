c***********************************************************************
      subroutine rhslom_ad(q,s,jd,kd,js,je,ks,ke,bt)
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c*********************************************************************c

         real q(jd,kd,nq),s(jd,kd,nq),bt(jd,kd)
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
               s(j,k,1) = phiSq*tmp(1)+u*tmp(2)+v*tmp(3)+ge*tmp(4)
               s(j,k,2) = -u*s(j,k,1)
               s(j,k,3) = -v*s(j,k,1)
               s(j,k,4) = s(j,k,1)

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
      subroutine timespectralrhs_bq(q, qb, s, sb, Ds, jd, kd, nsp)
      use params_global
      implicit none
c
c..compute time-spectral coupling term
c
c***********************************************************************
c
      integer jd, kd, nsp
c local variables
      real q(jd, kd, nspec, nq), s(jd, kd, nspec, nv), Ds(nspec, nspec)
      real qb(jd, kd, nspec, nq), sb(jd, kd, nspec, nv)
c
      integer j, k, n, nsp1
      do n=nv,1,-1
        do nsp1=nspec,1,-1
          do k=kend,kbeg,-1
            do j=jend,jbeg,-1
              qb(j, k, nsp, n) = qb(j, k, nsp, n) - rotf*Ds(nsp1, nsp)
     +          *sb(j, k, nsp1, n)
            enddo
          enddo
        enddo
      enddo
      end

c*************************************************************************
