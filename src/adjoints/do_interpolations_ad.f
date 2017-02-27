c***********************************************************************
      subroutine do_interpolations_ad(qg,qgb,jmx,kmx,ibcg,imeshg,idonorg,
     c        fracg,nfringeg,ndonorg,iisptrg,iieptrg,idsize,qsize,nmesh)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer idsize,qsize,nmesh
      integer jmx(nmesh),kmx(nmesh)
      integer nfringeg(nmesh,nspec),ndonorg(nmesh,nspec)
      integer iisptrg(nmesh,nspec),iieptrg(nmesh,nspec)
      real qg(qsize)!qg is state 
      real qgb(qsize) !qgb is adjoint rhs
      integer imeshg(idsize,2,nmesh,nspec),idonorg(idsize,2,nmesh,nspec)
      integer ibcg(idsize,nmesh,nspec)
      real fracg(idsize,2,nmesh,nspec)

c..local variables

      integer bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im
      integer qptr,nfringe,ndonor,iisptr,iieptr,nsp,qgptr,js,je
      real djm1,dj0,djp1,dkm1,dk0,dkp1
      real w1,w2,w3,w4,w5,w6,w7,w8,w9,onefourth
      real,allocatable :: qbcb(:,:),qbc(:,:)
      real,allocatable :: qb(:,:,:),q(:,:,:) !q qtate

      onefourth = 1./4

!!!$OMP PARALLEL IF (NSPEC > 1)
!!
!!!$OMP DO
!!!$OMP& PRIVATE(bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im)
!!!$OMP& PRIVATE(qptr,nfringe,ndonor,iisptr,iieptr)
!!!$OMP& PRIVATE(djm1,dj0,djp1,dkm1,dk0,dkp1)
!!!$OMP& PRIVATE(w1,w2,w3,w4,w5,w6,w7,w8,w9)
!!!$OMP& PRIVATE(q,qbc)

      !first compute state qbc 

      spectralloop:do nsp = 1,nspec

        bcdim=iieptrg(nmesh,nsp)
        allocate(qbc(bcdim,nv))
        allocate(qbcb(bcdim,nq)) !not nv
        qbc = 0.0
        !------------------------------------------
        !begin initializing qbc (qbc analysis state)
        !------------------------------------------
        qptr = 1
        do im = 1,nmesh

          ndonor = ndonorg(im,nsp)
          iisptr = iisptrg(im,nsp)
          iieptr = iieptrg(im,nsp)
          jmax = jmx(im); kmax = kmx(im)
          
          allocate(q(jmax,kmax,nq))
         
c.....assign local qb from global values

          do n = 1,nq
            do k=1,kmax
              do j=1,jmax
                q(j,k,n) = qg(qptr-1 + jmax*kmax*nspec*(n-1) + 
     &                     jmax*kmax*(nsp-1) + jmax*(k-1) + j)
              enddo
            enddo
          enddo

          qptr = qptr + jmax*kmax*nspec*nq
   
          do id=1,ndonor

            ii=idonorg(id,1,im,nsp)
            jj=idonorg(id,2,im,nsp)

            iim=ii-1
            jjm=jj-1

            iip=ii+1
            jjp=jj+1

            djm1 = 1.+fracg(id,1,im,nsp)
            dkm1 = 1.+fracg(id,2,im,nsp)
            dj0  = fracg(id,1,im,nsp)
            dk0  = fracg(id,2,im,nsp)
            djp1 = 1.-fracg(id,1,im,nsp)
            dkp1 = 1.-fracg(id,2,im,nsp)

            w1   =    (dj0  * djp1 * dk0  * dkp1)*onefourth
            w2   = -2*(djm1 * djp1 * dk0  * dkp1)*onefourth
            w3   = -  (djm1 * dj0  * dk0  * dkp1)*onefourth
            w4   = -2*(dj0  * djp1 * dkm1 * dkp1)*onefourth
            w5   =  4*(djm1 * djp1 * dkm1 * dkp1)*onefourth
            w6   =  2*(djm1 * dj0  * dkm1 * dkp1)*onefourth
            w7   = -  (dj0  * djp1 * dkm1 * dk0 )*onefourth
            w8   =  2*(djm1 * djp1 * dkm1 * dk0 )*onefourth
            w9   =    (djm1 * dj0  * dkm1 * dk0 )*onefourth

c.....collect in global pointer qbc from pointer iisptr->iieptr

            do n=1,nmv
               qbc(iisptr-1+id,n)= 
     &                     w1*q(iim,jjm,n)*q(iim,jjm,nq) 
     &                   + w2*q(ii ,jjm,n)*q(ii ,jjm,nq)
     &                   + w3*q(iip,jjm,n)*q(iip,jjm,nq) 
     &                   + w4*q(iim,jj ,n)*q(iim,jj ,nq)
     &                   + w5*q(ii ,jj ,n)*q(ii ,jj ,nq)
     &                   + w6*q(iip,jj ,n)*q(iip,jj ,nq)
     &                   + w7*q(iim,jjp,n)*q(iim,jjp,nq)
     &                   + w8*q(ii ,jjp,n)*q(ii ,jjp,nq)
     &                   + w9*q(iip,jjp,n)*q(iip,jjp,nq)
            enddo

            do n=nmv+1,nv
               qbc(iisptr-1+id,n)= 
     &                     w1*q(iim,jjm,n)
     &                   + w2*q(ii ,jjm,n)
     &                   + w3*q(iip,jjm,n)
     &                   + w4*q(iim,jj ,n)
     &                   + w5*q(ii ,jj ,n)
     &                   + w6*q(iip,jj ,n)
     &                   + w7*q(iim,jjp,n)
     &                   + w8*q(ii ,jjp,n)
     &                   + w9*q(iip,jjp,n)
            enddo

          enddo

          deallocate(q)
        enddo
        !print*,'adj: sum(qbc) :',sum(qbc)
        !------------------------------------------
        !end initializing qbc (qbc analysis state)
        !------------------------------------------



        !begin adjoint step
        !--------------------------------------------------------

        qbcb = 0.0 !adjoint global qbcb

c...LOOP THROUGH ALL MESHES AND UPDATE VALUES TO QBC ARRAY 
c...from q_ad (here qb)


        !qbcb <---q_ad
        qptr = 1
        do im = 1,nmesh

          nfringe = nfringeg(im,nsp)
          jmax = jmx(im); kmax = kmx(im)

          allocate(qb(jmax,kmax,nq)) !not nv
          allocate(q(jmax,kmax,nq))
          qb = 0.d0
          q = 0.d0

c.....assign local q_ad (qb here) from global values qgb
c.....i.e. q_ad <-- qb1
          !state first
          do n=1,nq !check if nq or nv
            do k=1,kmax
              do j=1,jmax
                qgptr = qptr-1 + jmax*kmax*nspec*(n-1)
     &                         + jmax*kmax*(nsp-1) + jmax*(k-1) + j

                q(j,k,n) = qg(qgptr) !analysis state q
              enddo
            enddo
          enddo

          !adjoint 
          do n=1,nq !check if nq or nv
          !am? do n=1,nv !check if nq or nv
            do k=1,kmax
              do j=1,jmax
                qgptr = qptr-1 + jmax*kmax*nspec*(n-1)
     &                         + jmax*kmax*(nsp-1) + jmax*(k-1) + j

                qb(j,k,n) = qgb(qgptr)
                qgb(qgptr) = 0.0 !gets updated later from qb(:,:,:)
              enddo
            enddo
          enddo

c..... analysis: over write fringe points solution w/ donor global qbc solution
c..... adjoint : initialize to donor qbcb from fringe points q_ad (qb here)

          !qbcb<--q_ad
          do id=1,nfringe
            j = ibcg(id,im,nsp)
            ii = imeshg(id,1,im,nsp)
            jj = imeshg(id,2,im,nsp)

            !am? qb(ii,jj,nq) = 0.0 !asitav test
            do n = 1,nmv
               !amm q(ii,jj,n) = qbc(j,n)/q(ii,jj,nq)
               qbcb(j,n) = qbcb(j,n) + qb(ii,jj,n)/q(ii,jj,nq)
               qb(ii,jj,nq) = qb(ii,jj,nq)
     &                      - qbc(j,n)*qb(ii,jj,n)/q(ii,jj,nq)**2
               qb(ii,jj,n) = 0.0 !asitav test >>>>>>>>>>>>>>>>>>>>>>>
            enddo

            do n = nmv+1,nv !not nv?
               !amm q(ii,jj,n) = qbc(j,n) !analysis
               qbcb(j,n) = qbcb(j,n) + qb(ii,jj,n)
               qb(ii,jj,n) = 0.0 !asitav test >>>>>>>>>>>>>>>>>>>>>>>
            enddo

          enddo !id=1,nfringe
          !print*,'adj sum(qb): ',sum(qb)

          !update back qgb with updated qb values from above
          do n=1,nq !check if nq or nv
            do k=1,kmax
              do j=1,jmax
                qgptr = qptr-1 + jmax*kmax*nspec*(n-1)
     &                         + jmax*kmax*(nsp-1) + jmax*(k-1) + j
                qgb(qgptr) = qgb(qgptr) + qb(j,k,n)  !update
                !am? qgb(qgptr) = qb(j,k,n)  !assign?
              enddo
            enddo
          enddo
         
          qptr = qptr + jmax*kmax*nspec*nq !same dims as qb, not nv

          deallocate(qb)
          deallocate(q)
        enddo  !im=1,nmesh
        !print*,'adj sum(qgb,qbcb): ',sum(qgb),sum(qbcb)

c...analysis: LOOP THROUGH ALL THE MESHES AND COLLECT GLOBAL QBC
c...adjoint : LOOP THROUGH ALL THE MESHES AND COLLECT from GLOBAL QBC to local q_ad

        !(qgb here) <-- qbcb
        qptr = 1
        do im = 1,nmesh

          ndonor = ndonorg(im,nsp)
          iisptr = iisptrg(im,nsp)
          iieptr = iieptrg(im,nsp)
          jmax = jmx(im); kmax = kmx(im)
          
          allocate(qb(jmax,kmax,nq)) !not nv
          allocate(q(jmax,kmax,nq))
          qb = 0.d0
          q = 0.d0
         
c.....initialize local q_ad (qb here) from global values qgb?

!am?          do n = 1,nq !not nv?
!am?            do k=1,kmax
!am?              do j=1,jmax
!am?                qgptr = qptr-1 + jmax*kmax*nspec*(n-1)
!am?     &                         + jmax*kmax*(nsp-1) + jmax*(k-1) + j
!am?                qb(j,k,n) = qgb(qgptr)
!am?                !am? qgb(qgptr) = 0.0 !asitav gets updated later from qb
!am?              enddo
!am?            enddo
!am?          enddo

          !state
          do n = 1,nq !for state
            do k=1,kmax
              do j=1,jmax
                qgptr = qptr-1 + jmax*kmax*nspec*(n-1)
     &                         + jmax*kmax*(nsp-1) + jmax*(k-1) + j
                q(j,k,n) = qg(qgptr)
              enddo
            enddo
          enddo

          !update/accumulate q_ad (donors) from updated qbcb (or qb1)
          do id=1,ndonor

            ii=idonorg(id,1,im,nsp)
            jj=idonorg(id,2,im,nsp)

            iim=ii-1
            jjm=jj-1

            iip=ii+1
            jjp=jj+1

            djm1 = 1.+fracg(id,1,im,nsp)
            dkm1 = 1.+fracg(id,2,im,nsp)
            dj0  = fracg(id,1,im,nsp)
            dk0  = fracg(id,2,im,nsp)
            djp1 = 1.-fracg(id,1,im,nsp)
            dkp1 = 1.-fracg(id,2,im,nsp)

            w1   =    (dj0  * djp1 * dk0  * dkp1)*onefourth
            w2   = -2*(djm1 * djp1 * dk0  * dkp1)*onefourth
            w3   = -  (djm1 * dj0  * dk0  * dkp1)*onefourth
            w4   = -2*(dj0  * djp1 * dkm1 * dkp1)*onefourth
            w5   =  4*(djm1 * djp1 * dkm1 * dkp1)*onefourth
            w6   =  2*(djm1 * dj0  * dkm1 * dkp1)*onefourth
            w7   = -  (dj0  * djp1 * dkm1 * dk0 )*onefourth
            w8   =  2*(djm1 * djp1 * dkm1 * dk0 )*onefourth
            w9   =    (djm1 * dj0  * dkm1 * dk0 )*onefourth

c.....analysis: collect in global pointer qbcb from pointer iisptr->iieptr
c.....adjoint : collect from global pointer qbcb to q_ad (qb here) from pointer iisptr->iieptr

            do n=1,nmv
               qb(iim,jjm,n)  = qb(iim,jjm,n)  + w1*qbcb(iisptr-1+id,n)*q(iim,jjm,nq)
               qb(iim,jjm,nq) = qb(iim,jjm,nq) + w1*qbcb(iisptr-1+id,n)*q(iim,jjm,n)

               qb(ii ,jjm,n)  = qb(ii ,jjm,n)  + w2*qbcb(iisptr-1+id,n)*q(ii ,jjm,nq)
               qb(ii ,jjm,nq) = qb(ii ,jjm,nq) + w2*qbcb(iisptr-1+id,n)*q(ii ,jjm,n)

               qb(iip,jjm,n)  = qb(iip,jjm,n)  + w3*qbcb(iisptr-1+id,n)*q(iip,jjm,nq)
               qb(iip,jjm,nq) = qb(iip,jjm,nq) + w3*qbcb(iisptr-1+id,n)*q(iip,jjm,n)

               qb(iim,jj ,n)  = qb(iim,jj ,n)  + w4*qbcb(iisptr-1+id,n)*q(iim,jj ,nq)
               qb(iim,jj ,nq) = qb(iim,jj ,nq) + w4*qbcb(iisptr-1+id,n)*q(iim,jj ,n)

               qb(ii ,jj ,n)  = qb(ii ,jj ,n)  + w5*qbcb(iisptr-1+id,n)*q(ii ,jj ,nq)
               qb(ii ,jj ,nq) = qb(ii ,jj ,nq) + w5*qbcb(iisptr-1+id,n)*q(ii ,jj ,n)

               qb(iip,jj ,n)  = qb(iip,jj ,n)  + w6*qbcb(iisptr-1+id,n)*q(iip,jj ,nq)
               qb(iip,jj ,nq) = qb(iip,jj ,nq) + w6*qbcb(iisptr-1+id,n)*q(iip,jj ,n)

               qb(iim,jjp,n)  = qb(iim,jjp,n)  + w7*qbcb(iisptr-1+id,n)*q(iim,jjp,nq)
               qb(iim,jjp,nq) = qb(iim,jjp,nq) + w7*qbcb(iisptr-1+id,n)*q(iim,jjp,n)

               qb(ii ,jjp,n)  = qb(ii ,jjp,n)  + w8*qbcb(iisptr-1+id,n)*q(ii ,jjp,nq)
               qb(ii ,jjp,nq) = qb(ii ,jjp,nq) + w8*qbcb(iisptr-1+id,n)*q(ii ,jjp,n)

               qb(iip,jjp,n)  = qb(iip,jjp,n)  + w9*qbcb(iisptr-1+id,n)*q(iip,jjp,nq)
               qb(iip,jjp,nq) = qb(iip,jjp,nq) + w9*qbcb(iisptr-1+id,n)*q(iip,jjp,n)

            enddo

            do n=nmv+1,nv !same as analysis
               qb(iim,jjm,n) = qb(iim,jjm,n) + w1*qbcb(iisptr-1+id,n)
               qb(ii ,jjm,n) = qb(ii ,jjm,n) + w2*qbcb(iisptr-1+id,n)
               qb(iip,jjm,n) = qb(iip,jjm,n) + w3*qbcb(iisptr-1+id,n)
               qb(iim,jj ,n) = qb(iim,jj ,n) + w4*qbcb(iisptr-1+id,n)
               qb(ii ,jj ,n) = qb(ii ,jj ,n) + w5*qbcb(iisptr-1+id,n)
               qb(iip,jj ,n) = qb(iip,jj ,n) + w6*qbcb(iisptr-1+id,n)
               qb(iim,jjp,n) = qb(iim,jjp,n) + w7*qbcb(iisptr-1+id,n)
               qb(ii ,jjp,n) = qb(ii ,jjp,n) + w8*qbcb(iisptr-1+id,n)
               qb(iip,jjp,n) = qb(iip,jjp,n) + w9*qbcb(iisptr-1+id,n)
            enddo

          enddo !1,ndonor

c.....reassign qbcb or updated q_ad (qb here) to update global qb (containing all Ng meshes)
          do n=1,nq ! same as in analysis? not nv?
            do k=1,kmax
              do j=1,jmax
                qgptr = qptr-1 + jmax*kmax*nspec*(n-1)
     &                         + jmax*kmax*(nsp-1) + jmax*(k-1) + j
                qgb(qgptr) = qgb(qgptr) + qb(j,k,n) !update, not assign
                !am? qgb(qgptr) = qb(j,k,n) !assign
              enddo
            enddo
          enddo

          qptr = qptr + jmax*kmax*nspec*nq

          deallocate(qb)
          deallocate(q)
        enddo !im=1,nmesh
        !print*,'adj sum(qgb): ',sum(qgb)

        deallocate(qbcb)
        deallocate(qbc)

      enddo spectralloop
!!!$OMP END DO

!!!$OMP END PARALLEL

      end subroutine do_interpolations_ad
