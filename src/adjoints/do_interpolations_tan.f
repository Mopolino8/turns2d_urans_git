c***********************************************************************
      subroutine do_interpolations_tan(qg,qgb,jmx,kmx,ibcg,imeshg,idonorg,
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
      real qgb(qsize) !qgb is tangent rhs
      integer imeshg(idsize,2,nmesh,nspec),idonorg(idsize,2,nmesh,nspec)
      integer ibcg(idsize,nmesh,nspec)
      real fracg(idsize,2,nmesh,nspec)

c..local variables

      integer bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im
      integer qptr,nfringe,ndonor,iisptr,iieptr,nsp,qgptr,js,je
      real djm1,dj0,djp1,dkm1,dk0,dkp1
      real w1,w2,w3,w4,w5,w6,w7,w8,w9,onefourth
      real,allocatable :: qbcb(:,:),qbc(:,:) !qbc analysis state
      real,allocatable :: qb(:,:,:),q(:,:,:) !qb adj state, q analysis qtate

      onefourth = 1./4

!!!$OMP PARALLEL IF (NSPEC > 1)
!!
!!!$OMP DO
!!!$OMP& PRIVATE(bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im)
!!!$OMP& PRIVATE(qptr,nfringe,ndonor,iisptr,iieptr)
!!!$OMP& PRIVATE(djm1,dj0,djp1,dkm1,dk0,dkp1)
!!!$OMP& PRIVATE(w1,w2,w3,w4,w5,w6,w7,w8,w9)
!!!$OMP& PRIVATE(q,qbc)

      !first compute state qbc or qbc here

      spectralloop:do nsp = 1,nspec

        bcdim=iieptrg(nmesh,nsp)
        allocate(qbcb(bcdim,nv)) !not nv
        allocate(qbc(bcdim,nv))
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
         
c.....assign local q from global values

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

c.....collect in global pointer qbcb from pointer iisptr->iieptr

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
        !print*,'tan: sum(qbc) :',sum(qbc)
        !------------------------------------------
        !end initializing qbc (qbc analysis state)
        !------------------------------------------

        !begin tangent step
        !--------------------------------------------------------

        qbcb = 0.0 !tangent global qbcb !d(qbc)/dD

c...LOOP THROUGH ALL MESHES AND UPDATE VALUES TO QBC ARRAY 
c...from q_tan (here qb)

        !qbcb <---q_tan
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

c.....assign local q_tan (qb here) from global values qgb 
c.....i.e. q_ad <-- qgb
          do n=1,nq !check if nq or nv
            do k=1,kmax
              do j=1,jmax
                qgptr = qptr-1 + jmax*kmax*nspec*(n-1)
     &                         + jmax*kmax*(nsp-1) + jmax*(k-1) + j
                q(j,k,n) = qg(qgptr) !analysis state
                qb(j,k,n) = qgb(qgptr)
              enddo
            enddo
          enddo

          qptr = qptr + jmax*kmax*nspec*nq

          do id=1,ndonor

            ii=idonorg(id,1,im,nsp)
            jj=idonorg(id,2,im,nsp)

            !if(id.eq.1 .and. im.eq.2) then
            !  qb(ii,jj,1) = 1.0
            !end if

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

c.....collect in global pointer qbcb from pointer iisptr->iieptr

            !qb is d(q)/dD actually, q: analysis state
            do n=1,nmv
               qbcb(iisptr-1+id,n)= 
     &                     w1*qb(iim,jjm,n)*q(iim,jjm,nq) 
     &                   + w1*q(iim,jjm,n)*qb(iim,jjm,nq) 
     &                   + w2*qb(ii ,jjm,n)*q(ii ,jjm,nq)
     &                   + w2*q(ii ,jjm,n)*qb(ii ,jjm,nq)
     &                   + w3*qb(iip,jjm,n)*q(iip,jjm,nq) 
     &                   + w3*q(iip,jjm,n)*qb(iip,jjm,nq) 
     &                   + w4*qb(iim,jj ,n)*q(iim,jj ,nq)
     &                   + w4*q(iim,jj ,n)*qb(iim,jj ,nq)
     &                   + w5*qb(ii ,jj ,n)*q(ii ,jj ,nq)
     &                   + w5*q(ii ,jj ,n)*qb(ii ,jj ,nq)
     &                   + w6*qb(iip,jj ,n)*q(iip,jj ,nq)
     &                   + w6*q(iip,jj ,n)*qb(iip,jj ,nq)
     &                   + w7*qb(iim,jjp,n)*q(iim,jjp,nq)
     &                   + w7*q(iim,jjp,n)*qb(iim,jjp,nq)
     &                   + w8*qb(ii ,jjp,n)*q(ii ,jjp,nq)
     &                   + w8*q(ii ,jjp,n)*qb(ii ,jjp,nq)
     &                   + w9*qb(iip,jjp,n)*q(iip,jjp,nq)
     &                   + w9*q(iip,jjp,n)*qb(iip,jjp,nq)
            enddo

            do n=nmv+1,nv
               qbcb(iisptr-1+id,n)= 
     &                     w1*qb(iim,jjm,n)
     &                   + w2*qb(ii ,jjm,n)
     &                   + w3*qb(iip,jjm,n)
     &                   + w4*qb(iim,jj ,n)
     &                   + w5*qb(ii ,jj ,n)
     &                   + w6*qb(iip,jj ,n)
     &                   + w7*qb(iim,jjp,n)
     &                   + w8*qb(ii ,jjp,n)
     &                   + w9*qb(iip,jjp,n)
            enddo

          enddo !id=1,ndonor

          deallocate(qb)
          deallocate(q)
        enddo !im=1,nmesh

c...LOOP THROUGH ALL MESHES AND UPDATE VALUES FROM QBC ARRAY

!!test 
!!-----------------------
!        qbcb = 1.0
!        qgb = 1.0
!!-----------------------

        qptr = 1
        do im = 1,nmesh

          nfringe = nfringeg(im,nsp)
          jmax = jmx(im); kmax = kmx(im)

          allocate(qb(jmax,kmax,nq))
          allocate(q(jmax,kmax,nq))

c.....re-assign local qb from global values
          !qb is dq/dD actually
          do n=1,nq
            do k=1,kmax
              do j=1,jmax
                qb(j,k,n) = qgb(qptr-1 + jmax*kmax*nspec*(n-1) + 
     &                      jmax*kmax*(nsp-1) + jmax*(k-1) + j)
                q(j,k,n) = qg(qptr-1 + jmax*kmax*nspec*(n-1) + 
     &                      jmax*kmax*(nsp-1) + jmax*(k-1) + j)
              enddo
            enddo
          enddo

c.....over write fringe points solution w/ donor global qbc solution

          !qbcb is d(qbc)/dD actually
          do id=1,nfringe
            j = ibcg(id,im,nsp)
            ii = imeshg(id,1,im,nsp)
            jj = imeshg(id,2,im,nsp)

            do n = 1,nmv
               !amm q(ii,jj,n) = qbc(j,n)/q(ii,jj,nq)
               qb(ii,jj,n) = qbcb(j,n)/q(ii,jj,nq) 
     &                     - qbc(j,n)*qb(ii,jj,nq)/q(ii,jj,nq)**2
            enddo

            do n = nmv+1,nv
               qb(ii,jj,n) = qbcb(j,n)
            enddo

          enddo
       
c.....reassign qbcb to global qb (containing all Ng meshes)
          !am? do n=1,nv
          do n=1,nq
            do k=1,kmax
              do j=1,jmax
                qgb(qptr-1 + jmax*kmax*nspec*(n-1) + jmax*kmax*(nsp-1) + 
     &               jmax*(k-1) + j) = qb(j,k,n)
                !write(5000+im,*)qb(j,k,n)
              enddo
            enddo
          enddo

          qptr = qptr + jmax*kmax*nspec*nq

          deallocate(qb)
          deallocate(q)
        enddo 
        print*,'tan sum(qgb): ',sum(qgb)

        deallocate(qbcb)
        deallocate(qbc)

      enddo spectralloop
!!!$OMP END DO

!!!$OMP END PARALLEL

      end subroutine do_interpolations_tan
