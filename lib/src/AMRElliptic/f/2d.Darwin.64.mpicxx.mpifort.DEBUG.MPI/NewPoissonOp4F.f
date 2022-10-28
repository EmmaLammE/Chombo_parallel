      subroutine OPERATORLAP4(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlofphicomp
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dxinv, lap
      integer n,ncomp
      integer i,j
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          lap = ( 
     &   (16.0d0)*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + (16.0d0)*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + (16.0d0)*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + (16.0d0)*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -((30.0d0)*2)*phi(i,j,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
          lofphi(i,j,n) = alpha*phi(i,j,n)+beta*lap
      enddo
      enddo
      enddo
      return
      end
      subroutine OPERATORLAPRES4(
     &           r
     &           ,irlo0,irlo1
     &           ,irhi0,irhi1
     &           ,nrcomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrcomp
      integer irlo0,irlo1
      integer irhi0,irhi1
      REAL*8 r(
     &           irlo0:irhi0,
     &           irlo1:irhi1,
     &           0:nrcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dxinv, lap, lhs
      integer n,ncomp
      integer i,j
      ncomp = nphicomp
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          lap = ( 
     &   (16.0d0)*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + (16.0d0)*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + (16.0d0)*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + (16.0d0)*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -((30.0d0)*2)*phi(i,j,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
          lhs = alpha*phi(i,j,n) + beta*lap
          r(i,j,n) = rhs(i,j,n) - lhs
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRES4(
     &           res
     &           ,ireslo0,ireslo1
     &           ,ireshi0,ireshi1
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,alpha
     &           ,beta
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0,ireslo1
      integer ireshi0,ireshi1
      REAL*8 res(
     &           ireslo0:ireshi0,
     &           ireslo1:ireshi1,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      REAL*8 alpha
      REAL*8 beta
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 denom,dxinv,lofphi,lap
      integer n,ncomp
      integer i,j
      integer ii,jj
      ncomp = nphicomp
      dxinv = (1.0d0) / (dx*dx)
      denom = 2  *2
      do n = 0, ncomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          ii = i/2 
          jj = j/2 
          lap = ( 
     &   (16.0d0)*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + (16.0d0)*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + (16.0d0)*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + (16.0d0)*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -((30.0d0)*2)*phi(i,j,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
          lofphi = alpha*phi(i,j,n) + beta*lap
          res(ii,jj,n) = res(ii,jj,n)
     &                            + (rhs(i,j,n) - lofphi) / denom
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBLAPLACIAN4(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,tmp
     &           ,itmplo0,itmplo1
     &           ,itmphi0,itmphi1
     &           ,ntmpcomp
     &           ,redBlack
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      integer ntmpcomp
      integer itmplo0,itmplo1
      integer itmphi0,itmphi1
      REAL*8 tmp(
     &           itmplo0:itmphi0,
     &           itmplo1:itmphi1,
     &           0:ntmpcomp-1)
      integer redBlack
      REAL*8 dx2t, thD
      integer i,j
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = (12.0d0)*dx*dx
      thD  = (1.000d0 / 30.000d0)/2
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               tmp(i,j,n) = thD*( 
     &           (16.0d0)*phi(i+1,j,n) - phi(i+2,j,n)
     &         + (16.0d0)*phi(i-1,j,n) - phi(i-2,j,n)
     &         + (16.0d0)*phi(i,j+1,n) - phi(i,j+2,n)
     &         + (16.0d0)*phi(i,j-1,n) - phi(i,j-2,n)
     &         - dx2t*rhs(i,j,n) )
            enddo
          enddo
        else if (redBlack .eq. black) then
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,n) = thD*( 
     &           (16.0d0)*tmp(i+1,j,n) - tmp(i+2,j,n)
     &         + (16.0d0)*tmp(i-1,j,n) - tmp(i-2,j,n)
     &         + (16.0d0)*tmp(i,j+1,n) - tmp(i,j+2,n)
     &         + (16.0d0)*tmp(i,j-1,n) - tmp(i,j-2,n)
     &         - dx2t*rhs(i,j,n) )
            enddo
          enddo
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,n) = tmp(i,j,n)
            enddo
          enddo
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine SORLAPLACIAN4(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 dx2t, thD, tmp, omega
      integer i,j
      integer n,ncomp
      dx2t = (12.0d0)*dx*dx
      thD  = (1.000d0 / 30.000d0)/2
      omega = 0.47
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         tmp = thD*( 
     &        (16.0d0)*phi(i+1,j,n) - phi(i+2,j,n)
     &        + (16.0d0)*phi(i-1,j,n) - phi(i-2,j,n)
     &        + (16.0d0)*phi(i,j+1,n) - phi(i,j+2,n)
     &        + (16.0d0)*phi(i,j-1,n) - phi(i,j-2,n)
     &        - dx2t*rhs(i,j,n) )
         phi(i,j,n) = omega*tmp
     &        + ((1.0d0)-omega)*phi(i,j,n)
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBHELMHOLTZ4(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,tmp
     &           ,itmplo0,itmplo1
     &           ,itmphi0,itmphi1
     &           ,ntmpcomp
     &           ,alpha
     &           ,beta
     &           ,redBlack
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      integer ntmpcomp
      integer itmplo0,itmplo1
      integer itmphi0,itmphi1
      REAL*8 tmp(
     &           itmplo0:itmphi0,
     &           itmplo1:itmphi1,
     &           0:ntmpcomp-1)
      REAL*8 alpha
      REAL*8 beta
      integer redBlack
      REAL*8 dx2t, lambda, lap, dxinv, helm
      integer i,j
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = (12.0d0)*dx*dx
      dxinv = (1.0d0)/(dx*dx)
      lambda = (1.0d0)/(alpha - beta*(30.0d0)*2*(1.000d0 / 12.000d0)*dxi
     &nv)
      lambda = lambda*(0.60)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
          lap = ( 
     &   (16.0d0)*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + (16.0d0)*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + (16.0d0)*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + (16.0d0)*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -((30.0d0)*2)*phi(i,j,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
          helm = alpha*phi(i,j,n) + beta*lap
          tmp(i,j,n) = phi(i,j,n) +
     &      lambda*( rhs(i,j,n) - helm )
            enddo
          enddo
        else if (redBlack .eq. black) then
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               lap = ( 
     &           (16.0d0)*tmp(i+1,j,n) - tmp(i+2,j,n)
     &         + (16.0d0)*tmp(i-1,j,n) - tmp(i-2,j,n)
     &         + (16.0d0)*tmp(i,j+1,n) - tmp(i,j+2,n)
     &         + (16.0d0)*tmp(i,j-1,n) - tmp(i,j-2,n)
     &                     -((30.0d0)*2)*tmp(i,j,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
               helm = alpha*tmp(i,j,n) + beta*lap
               phi(i,j,n) = tmp(i,j,n) +
     &              lambda*( rhs(i,j,n) - helm )
            enddo
          enddo
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,n) = tmp(i,j,n)
            enddo
          enddo
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine NEWGETFLUX4(
     &           flux
     &           ,ifluxlo0,ifluxlo1
     &           ,ifluxhi0,ifluxhi1
     &           ,nfluxcomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,beta_dx
     &           ,a_idir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1
      integer ifluxhi0,ifluxhi1
      REAL*8 flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 beta_dx
      integer a_idir
      INTEGER ncomp,n
      integer ii, jj
      integer i , j 
      ncomp = nphicomp
      ii = CHF_ID(a_idir, 0)
      jj = CHF_ID(a_idir, 1)
      do n = 0, ncomp-1
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          flux(i,j,n) = beta_dx * (1.000d0 / 12.000d0) *
     &        ( (15.0d0)*phi(i,j,n)
     &           + phi(i-2*ii,j-2*jj,n)
     &           - phi(i+ii,j+jj,n)
     &           - (15.0d0)*phi(i-ii,j-jj,n) )
      enddo
      enddo
      enddo
      return
      end
      subroutine PROLONGLINEAR(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,ncoarsecomp
     &           ,ifineBoxlo0,ifineBoxlo1
     &           ,ifineBoxhi0,ifineBoxhi1
     &           ,icrseBoxlo0,icrseBoxlo1
     &           ,icrseBoxhi0,icrseBoxhi1
     &           ,r
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           0:ncoarsecomp-1)
      integer ifineBoxlo0,ifineBoxlo1
      integer ifineBoxhi0,ifineBoxhi1
      integer icrseBoxlo0,icrseBoxlo1
      integer icrseBoxhi0,icrseBoxhi1
      integer r
      INTEGER ncomp, n
      integer i ,j 
      integer ic,jc
      ncomp = nphicomp
      do n = 0, ncomp-1
      do j = ifineBoxlo1,ifineBoxhi1
      do i = ifineBoxlo0,ifineBoxhi0
           ic = i/r
           jc = j/r
           phi(i,j,n) =  phi(i,j,n) +
     &          coarse(ic,jc,n)
           if (ic.ne.icrseBoxhi0 .and.
     &         (ic*r.lt.i .or. ic.eq.icrseBoxlo0)) then
              phi(i,j,n) =  phi(i,j,n) +
     &             (coarse(ic+1,jc,n)
     &              - coarse(ic,jc,n))/r*(i+(0.500d0)-ic*r-(0.500d0)*r)
           else
              phi(i,j,n) =  phi(i,j,n) +
     &             (- coarse(ic-1,jc,n)
     &              + coarse(ic,jc,n))/r*(i+(0.500d0)-ic*r-(0.500d0)*r)
           endif
           if (jc.ne.icrseBoxhi1 .and.
     &         (jc*r.lt.j .or. jc.eq.icrseBoxlo1)) then
              phi(i,j,n) =  phi(i,j,n) +
     &             (coarse(ic,jc+1,n)
     &              - coarse(ic,jc,n))/r*(j+(0.500d0)-jc*r-(0.500d0)*r)
           else
              phi(i,j,n) =  phi(i,j,n) +
     &             (- coarse(ic,jc-1,n)
     &              + coarse(ic,jc,n))/r*(j+(0.500d0)-jc*r-(0.500d0)*r)
           endif
      enddo
      enddo
      enddo
      return
      end
