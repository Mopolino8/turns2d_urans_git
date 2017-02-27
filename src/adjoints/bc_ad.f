c***********************************************************************
      subroutine bc_ad(sb,im,jd,kd)
c
c  explicitly update the mesh boundaries
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer im,jd,kd
      real sb(jd,kd,nv)

c***  first executable statement
      integer j,k,js,je,ks,ke,idir,ib

      do ib=1,nbc_all(im)
        js = jbcs_all(ib,im)
        je = jbce_all(ib,im)
        ks = kbcs_all(ib,im)
        ke = kbce_all(ib,im)
        if(js.lt.0) js = jmax + js + 1
        if(ks.lt.0) ks = kmax + ks + 1
        if(je.lt.0) je = jmax + je + 1
        if(ke.lt.0) ke = kmax + ke + 1

        if (js.gt.1) js = js + nhalo
        if (ks.gt.1) ks = ks + nhalo
        if (je.lt.jmax) je = je + nhalo - 1
        if (ke.lt.kmax) ke = ke + nhalo - 1
        if (je.eq.jmax) je = jd
        if (ke.eq.kmax) ke = kd

        idir = ibdir_all(ib,im)

c.. inviscid wall bc (k = 1 only) 
        if (ibtyp_all(ib,im).eq.3  .or.
     &      ibtyp_all(ib,im).eq.4  .or.
     &      ibtyp_all(ib,im).eq.5  .or.
     &      ibtyp_all(ib,im).eq.10 .or.
     &      ibtyp_all(ib,im).eq.11 .or.
     &      ibtyp_all(ib,im).eq.46 .or.
     &      ibtyp_all(ib,im).eq.47) then

          call bcextp_ad(sb,jd,kd,js,je,ks,ke,idir)

c.. periodic bc
        elseif (ibtyp_all(ib,im).eq.22) then
          call bcperiodic_ad(sb,jd,kd,js,je,ks,ke,idir)

c.. wake bc
        elseif (ibtyp_all(ib,im).eq.51) then
          call bcwake_ad(sb,jd,kd,js,je,ks,ke,idir)

        endif
      enddo

      return
      end

c***********************************************************************
      subroutine bcextp_ad(sb,jd,kd,js,je,ks,ke,idir)
c
c..extrapolation bc
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real sb(jd,kd,nv)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1,j2,k2,jc,kc
      integer iadd,iadir
      real foso

      foso = 1.0
c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      if (iadir.eq.1) then
        do jc = 1,je-js+1
          if (idir.eq.1) then
            j = je - jc + 1
          elseif (idir.eq.-1) then
            j = js + jc - 1
          endif
          j1 = j + iadd
          j2 = j1 + iadd
          do k = ks,ke
            sb(j,k,1) = (1.+foso)*sb(j1,k,1) - foso*sb(j2,k,1)
            sb(j,k,2) = (1.+foso)*sb(j1,k,2) - foso*sb(j2,k,2)
            sb(j,k,3) = (1.+foso)*sb(j1,k,3) - foso*sb(j2,k,3)
            sb(j,k,4) = (1.+foso)*sb(j1,k,4) - foso*sb(j2,k,4)
          enddo
        enddo
      elseif (iadir.eq.2) then
        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
          elseif (idir.eq.-2) then
            k = ks + kc - 1
          endif
          k1 = k + iadd
          k2 = k1 + iadd
          do j = js,je
            sb(j,k,1) = (1.+foso)*sb(j,k1,1) - foso*sb(j,k2,1)
            sb(j,k,2) = (1.+foso)*sb(j,k1,2) - foso*sb(j,k2,2)
            sb(j,k,3) = (1.+foso)*sb(j,k1,3) - foso*sb(j,k2,3)
            sb(j,k,4) = (1.+foso)*sb(j,k1,4) - foso*sb(j,k2,4)
          enddo
        enddo
      endif
c
      return
      end

c***********************************************************************
      subroutine bcperiodic_ad(sb,jd,kd,js,je,ks,ke,idir)
c
c..periodic boundary for overlapping planes
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd
      real sb(jd,kd,nv)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
      
      iadd = sign(1,idir)

      pi = 4*atan(1.0)
      
      if(idir.eq.1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = j - js + 1
           jc = jd - 2*jj + jj1
 
           do k = ks,ke
             sb(j,k,1) = sb(jc,k,1)
             sb(j,k,2) = sb(jc,k,2)
             sb(j,k,3) = sb(jc,k,3)
             sb(j,k,4) = sb(jc,k,4)
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j + 1
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             sb(j,k,1) = sb(jc,k,1)
             sb(j,k,2) = sb(jc,k,2)
             sb(j,k,3) = sb(jc,k,3)
             sb(j,k,4) = sb(jc,k,4)
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks + 1
           kc = kd - 2*kk + kk1
 
           do j = js,je
             sb(j,k,1) = sb(j,kc,1)
             sb(j,k,2) = sb(j,kc,2)
             sb(j,k,3) = sb(j,kc,3)
             sb(j,k,4) = sb(j,kc,4)
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k + 1
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             sb(j,k,1) = sb(j,kc,1)
             sb(j,k,2) = sb(j,kc,2)
             sb(j,k,3) = sb(j,kc,3)
             sb(j,k,4) = sb(j,kc,4)
           enddo
        enddo

      endif

      return
      end

c***********************************************************************
      subroutine bcwake_ad(sb,jd,kd,js,je,ks,ke,idir)
c
c  for complete cgrid   
c  treatment for in the trailing wake at k = 1
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real sb(jd,kd,nv)

      ! local variables
      
      integer k,k1,j,jj,kc,iadd,iadir
      
c**   first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwake_ad'
      elseif(iadir.eq.2) then

        do kc = 1,ke-ks+1
          if (idir.eq.2) then
            k = ke - kc + 1
            k1 = ke + kc
          elseif (idir.eq.-2) then
            k = ks + kc - 1
            k1 = ks - kc
          endif

c..in the wake, average points above and below

          do j = js, je
            jj = jd - j + 1
            sb(j,k,1)  = sb(jj,k1,1)
            sb(j,k,2)  = sb(jj,k1,2)
            sb(j,k,3)  = sb(jj,k1,3)
            sb(j,k,4)  = sb(jj,k1,4)
            sb(jj,k,1) = sb(j,k1,1)
            sb(jj,k,2) = sb(j,k1,2)
            sb(jj,k,3) = sb(j,k1,3)
            sb(jj,k,4) = sb(j,k1,4)
          enddo
       enddo
      
      endif
c
      return
      end

c***********************************************************************
