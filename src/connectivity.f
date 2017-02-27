c***********************************************************************
      subroutine do_interpolations(qg,jmx,kmx,ibcg,imeshg,idonorg,
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
      real qg(qsize)
      integer imeshg(idsize,2,nmesh,nspec),idonorg(idsize,2,nmesh,nspec)
      integer ibcg(idsize,nmesh,nspec)
      real fracg(idsize,2,nmesh,nspec)

c..local variables

      integer bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im
      integer qptr,nfringe,ndonor,iisptr,iieptr,nsp
      real djm1,dj0,djp1,dkm1,dk0,dkp1
      real w1,w2,w3,w4,w5,w6,w7,w8,w9,onefourth
      real,allocatable :: qbc(:,:)
      real,allocatable :: q(:,:,:)

      onefourth = 1./4

!!!$OMP PARALLEL IF (NSPEC > 1)
!!
!!!$OMP DO
!!!$OMP& PRIVATE(bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im)
!!!$OMP& PRIVATE(qptr,nfringe,ndonor,iisptr,iieptr)
!!!$OMP& PRIVATE(djm1,dj0,djp1,dkm1,dk0,dkp1)
!!!$OMP& PRIVATE(w1,w2,w3,w4,w5,w6,w7,w8,w9)
!!!$OMP& PRIVATE(q,qbc)

      spectralloop:do nsp = 1,nspec

        bcdim=iieptrg(nmesh,nsp)
        allocate(qbc(bcdim,nv))

c...LOOP THROUGH ALL THE MESHES AND COLLECT GLOBAL QBC

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

c...LOOP THROUGH ALL MESHES AND UPDATE VALUES FROM QBC ARRAY

        qptr = 1
        do im = 1,nmesh

          nfringe = nfringeg(im,nsp)
          jmax = jmx(im); kmax = kmx(im)

          allocate(q(jmax,kmax,nq))

c.....re-assign local q from global values
          do n=1,nq
            do k=1,kmax
              do j=1,jmax
                q(j,k,n) = qg(qptr-1 + jmax*kmax*nspec*(n-1) + 
     &                      jmax*kmax*(nsp-1) + jmax*(k-1) + j)
              enddo
            enddo
          enddo

c.....over write fringe points solution w/ donor global qbc solution

          do id=1,nfringe
            j = ibcg(id,im,nsp)
            ii = imeshg(id,1,im,nsp)
            jj = imeshg(id,2,im,nsp)

            do n = 1,nmv
               q(ii,jj,n) = qbc(j,n)/q(ii,jj,nq)
            enddo

            do n = nmv+1,nv
               q(ii,jj,n) = qbc(j,n)
            enddo

          enddo
       
c.....reassign qbc to global q (containing all Ng meshes)
          do n=1,nv
            do k=1,kmax
              do j=1,jmax
                qg(qptr-1 + jmax*kmax*nspec*(n-1) + jmax*kmax*(nsp-1) + 
     &               jmax*(k-1) + j) = q(j,k,n)
              enddo
            enddo
          enddo

          qptr = qptr + jmax*kmax*nspec*nq

          deallocate(q)
        enddo 

        deallocate(qbc)

      enddo spectralloop
!!!$OMP END DO

!!!$OMP END PARALLEL

      return
      end

c*************************************************************************
      subroutine update_halo_iblanks(iblank,jd,kd)
c.. update halo iblank
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************
      integer jd,kd
      integer iblank(jd,kd,nspec)

! local variables
      integer j,k,nsp

!$OMP PARALLEL IF(NSPEC > 1)

!$OMP DO
!$OMP& PRIVATE(j,k)
      spectralloop: do nsp = 1,nspec
        do k = 1,kd
          do j = 1,jbeg - 1
            iblank(j,k,nsp) = iblank(jbeg,k,nsp)
          enddo
          
          do j = jend+1,jd
            iblank(j,k,nsp) = iblank(jend,k,nsp)
          enddo
        enddo

        do j = 1,jd
          do k = 1,kbeg - 1
            iblank(j,k,nsp) = iblank(j,kbeg,nsp)
          enddo
          
          do k = kend+1,kd
            iblank(j,k,nsp) = iblank(j,kend,nsp)
          enddo
        enddo

      enddo spectralloop
!$OMP END DO
!$OMP END PARALLEL

      end subroutine update_halo_iblanks

c*************************************************************************
