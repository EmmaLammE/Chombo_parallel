      subroutine MINFLATF(
     &           flattening
     &           ,iflatteninglo0,iflatteninglo1
     &           ,iflatteninghi0,iflatteninghi1
     &           ,zetaDir
     &           ,izetaDirlo0,izetaDirlo1
     &           ,izetaDirhi0,izetaDirhi1
     &           ,nzetaDircomp
     &           ,du
     &           ,idulo0,idulo1
     &           ,iduhi0,iduhi1
     &           ,nducomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iflatteninglo0,iflatteninglo1
      integer iflatteninghi0,iflatteninghi1
      REAL*8 flattening(
     &           iflatteninglo0:iflatteninghi0,
     &           iflatteninglo1:iflatteninghi1)
      integer nzetaDircomp
      integer izetaDirlo0,izetaDirlo1
      integer izetaDirhi0,izetaDirhi1
      REAL*8 zetaDir(
     &           izetaDirlo0:izetaDirhi0,
     &           izetaDirlo1:izetaDirhi1,
     &           0:nzetaDircomp-1)
      integer nducomp
      integer idulo0,idulo1
      integer iduhi0,iduhi1
      REAL*8 du(
     &           idulo0:iduhi0,
     &           idulo1:iduhi1,
     &           0:nducomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer iv
      REAL*8  sumdu,minflattot,minZetaDir
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         sumdu = 0
         do iv = 0,nducomp - 1
            sumdu = sumdu + du(i0,i1,iv)
         enddo
         if (sumdu .lt. 0) then
            minflattot = zetaDir(i0,i1,0)
            do iv = 1,nducomp - 1
               minZetaDir = zetaDir(i0,i1,iv)
               minflattot = min(minflattot,minZetaDir)
            enddo
            flattening(i0,i1) = minflattot
         else
            flattening(i0,i1) = (1.0d0)
         endif
      enddo
      enddo
      return
      end
      subroutine GETDPTWOF(
     &           delta2p
     &           ,idelta2plo0,idelta2plo1
     &           ,idelta2phi0,idelta2phi1
     &           ,delta1p
     &           ,idelta1plo0,idelta1plo1
     &           ,idelta1phi0,idelta1phi1
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idelta2plo0,idelta2plo1
      integer idelta2phi0,idelta2phi1
      REAL*8 delta2p(
     &           idelta2plo0:idelta2phi0,
     &           idelta2plo1:idelta2phi1)
      integer idelta1plo0,idelta1plo1
      integer idelta1phi0,idelta1phi1
      REAL*8 delta1p(
     &           idelta1plo0:idelta1phi0,
     &           idelta1plo1:idelta1phi1)
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
         delta2p(i0,i1) = delta1p(i0 +ioff0,i1 +ioff1)
     &     + delta1p(i0 -ioff0,i1 -ioff1)
      enddo
      enddo
      if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
            delta2p(i0,i1) = delta1p(i0 +ioff0,i1 +ioff1)
     &                             + delta1p(i0,i1)
      enddo
      enddo
         endif
      if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
            delta2p(i0,i1) = delta1p(i0,i1)
     &                             + delta1p(i0 -ioff0,i1 -ioff1)
      enddo
      enddo
      endif
      return
      end
      subroutine GETFLATF(
     &           zetaTwiddle
     &           ,izetaTwiddlelo0,izetaTwiddlelo1
     &           ,izetaTwiddlehi0,izetaTwiddlehi1
     &           ,delta1p
     &           ,idelta1plo0,idelta1plo1
     &           ,idelta1phi0,idelta1phi1
     &           ,delta2p
     &           ,idelta2plo0,idelta2plo1
     &           ,idelta2phi0,idelta2phi1
     &           ,smallp
     &           ,bulkMin
     &           ,ibulkMinlo0,ibulkMinlo1
     &           ,ibulkMinhi0,ibulkMinhi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer izetaTwiddlelo0,izetaTwiddlelo1
      integer izetaTwiddlehi0,izetaTwiddlehi1
      REAL*8 zetaTwiddle(
     &           izetaTwiddlelo0:izetaTwiddlehi0,
     &           izetaTwiddlelo1:izetaTwiddlehi1)
      integer idelta1plo0,idelta1plo1
      integer idelta1phi0,idelta1phi1
      REAL*8 delta1p(
     &           idelta1plo0:idelta1phi0,
     &           idelta1plo1:idelta1phi1)
      integer idelta2plo0,idelta2plo1
      integer idelta2phi0,idelta2phi1
      REAL*8 delta2p(
     &           idelta2plo0:idelta2phi0,
     &           idelta2plo1:idelta2phi1)
      REAL*8 smallp
      integer ibulkMinlo0,ibulkMinlo1
      integer ibulkMinhi0,ibulkMinhi1
      REAL*8 bulkMin(
     &           ibulkMinlo0:ibulkMinhi0,
     &           ibulkMinlo1:ibulkMinhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 d,r0,r1,ratio,strength
      data d  /0.33d0/
      data r0 /0.75d0/
      data r1 /0.85d0/
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         strength = abs(delta1p(i0,i1)/bulkMin(i0,i1))
         if (strength .ge. d) then
            ratio =     abs(delta1p(i0,i1)
     &            / max(abs(delta2p(i0,i1)),smallp))
            if (ratio .le. r0) then
               zetaTwiddle(i0,i1) = (1.0d0)
            else if (ratio .ge. r1) then
               zetaTwiddle(i0,i1) = 0
            else
               zetaTwiddle(i0,i1) = (1.0d0) - (ratio - r0)/(r1 - r0)
            endif
         else
            zetaTwiddle(i0,i1) = (1.0d0)
         endif
      enddo
      enddo
      return
      end
      subroutine GETGRADF(
     &           du
     &           ,idulo0,idulo1
     &           ,iduhi0,iduhi1
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idulo0,idulo1
      integer iduhi0,iduhi1
      REAL*8 du(
     &           idulo0:iduhi0,
     &           idulo1:iduhi1)
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL*8 u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1)
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
         du(i0,i1) = (0.500d0)*(u(i0 +ioff0,i1 +ioff1)
     &                           - u(i0 -ioff0,i1 -ioff1))
      enddo
      enddo
      if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
            du(i0,i1) = (u(i0 +ioff0,i1 +ioff1)
     &                         - u(i0,i1))
      enddo
      enddo
      endif
      if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
            du(i0,i1) = (u(i0,i1)
     &                         - u(i0 -ioff0,i1 -ioff1))
      enddo
      enddo
      endif
      return
      end
      subroutine MIN3PTSF(
     &           mindata
     &           ,imindatalo0,imindatalo1
     &           ,imindatahi0,imindatahi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer imindatalo0,imindatalo1
      integer imindatahi0,imindatahi1
      REAL*8 mindata(
     &           imindatalo0:imindatahi0,
     &           imindatalo1:imindatahi1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1)
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
         mindata(i0,i1) = min(data(i0,i1),
     &     data(i0 +ioff0,i1 +ioff1),
     &     data(i0 -ioff0,i1 -ioff1))
      enddo
      enddo
      if (hasLo .ne. 0) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
            mindata(i0,i1) = min(data(i0,i1),
     &        data(i0 +ioff0,i1 +ioff1))
      enddo
      enddo
      endif
      if (hasHi .ne. 0) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
            mindata(i0,i1) = min(data(i0,i1),
     &        data(i0 -ioff0,i1 -ioff1))
      enddo
      enddo
      endif
      return
      end
      subroutine SECONDSLOPEDIFFSF(
     &           deltaWC
     &           ,ideltaWClo0,ideltaWClo1
     &           ,ideltaWChi0,ideltaWChi1
     &           ,ndeltaWCcomp
     &           ,deltaWL
     &           ,ideltaWLlo0,ideltaWLlo1
     &           ,ideltaWLhi0,ideltaWLhi1
     &           ,ndeltaWLcomp
     &           ,deltaWR
     &           ,ideltaWRlo0,ideltaWRlo1
     &           ,ideltaWRhi0,ideltaWRhi1
     &           ,ndeltaWRcomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndeltaWCcomp
      integer ideltaWClo0,ideltaWClo1
      integer ideltaWChi0,ideltaWChi1
      REAL*8 deltaWC(
     &           ideltaWClo0:ideltaWChi0,
     &           ideltaWClo1:ideltaWChi1,
     &           0:ndeltaWCcomp-1)
      integer ndeltaWLcomp
      integer ideltaWLlo0,ideltaWLlo1
      integer ideltaWLhi0,ideltaWLhi1
      REAL*8 deltaWL(
     &           ideltaWLlo0:ideltaWLhi0,
     &           ideltaWLlo1:ideltaWLhi1,
     &           0:ndeltaWLcomp-1)
      integer ndeltaWRcomp
      integer ideltaWRlo0,ideltaWRlo1
      integer ideltaWRhi0,ideltaWRhi1
      REAL*8 deltaWR(
     &           ideltaWRlo0:ideltaWRhi0,
     &           ideltaWRlo1:ideltaWRhi1,
     &           0:ndeltaWRcomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer lvar
      integer ioff0,ioff1
      REAL*8 dWR,dWL
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do lvar = 0,numSlopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            dWR = W(i0 +ioff0,i1 +ioff1,lvar)
     &        - W(i0,i1,lvar)
            dWL = W(i0,i1,lvar)
     &           - W(i0 -ioff0,i1 -ioff1,lvar)
            deltaWR(i0,i1,lvar) = dWR
            deltaWL(i0,i1,lvar) = dWL
            deltaWC(i0,i1,lvar) = (0.500d0)*(dWR + dWL)
      enddo
      enddo
         if (hasLo .ne. 0) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               dWR = W(i0 +ioff0,i1 +ioff1,lvar)
     &           - W(i0,i1,lvar)
               deltaWC(i0,i1,lvar) = dWR
               deltaWL(i0,i1,lvar) = dWR
               deltaWR(i0,i1,lvar) = dWR
      enddo
      enddo
         endif
         if (hasHi .ne. 0) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               dWL = W(i0,i1,lvar)
     &           - W(i0 -ioff0,i1 -ioff1,lvar)
               deltaWC(i0,i1,lvar) = dWL
               deltaWL(i0,i1,lvar) = dWL
               deltaWR(i0,i1,lvar) = dWL
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine FOURTHSLOPEDIFFSF(
     &           delta4WC
     &           ,idelta4WClo0,idelta4WClo1
     &           ,idelta4WChi0,idelta4WChi1
     &           ,ndelta4WCcomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,delta2W
     &           ,idelta2Wlo0,idelta2Wlo1
     &           ,idelta2Whi0,idelta2Whi1
     &           ,ndelta2Wcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndelta4WCcomp
      integer idelta4WClo0,idelta4WClo1
      integer idelta4WChi0,idelta4WChi1
      REAL*8 delta4WC(
     &           idelta4WClo0:idelta4WChi0,
     &           idelta4WClo1:idelta4WChi1,
     &           0:ndelta4WCcomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer ndelta2Wcomp
      integer idelta2Wlo0,idelta2Wlo1
      integer idelta2Whi0,idelta2Whi1
      REAL*8 delta2W(
     &           idelta2Wlo0:idelta2Whi0,
     &           idelta2Wlo1:idelta2Whi1,
     &           0:ndelta2Wcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      integer lvar
      REAL*8 dWR,dWL
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do lvar = 0,numSlopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            dWR =       W(i0 +ioff0,i1 +ioff1,lvar)
     &          - delta2W(i0 +ioff0,i1 +ioff1,lvar)*(0.250d0)
            dWL =       W(i0 -ioff0,i1 -ioff1,lvar)
     &          + delta2W(i0 -ioff0,i1 -ioff1,lvar)*(0.250d0)
            delta4WC(i0,i1,lvar) = (2.000d0 / 3.000d0)*(dWR - dWL)
      enddo
      enddo
         if (hasLo .ne. 0) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               delta4WC(i0,i1,lvar) =
     &           delta2W(i0,i1,lvar)
      enddo
      enddo
         endif
         if (hasHi .ne. 0) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               delta4WC(i0,i1,lvar) =
     &           delta2W(i0,i1,lvar)
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine APPLYFLATF(
     &           dW
     &           ,idWlo0,idWlo1
     &           ,idWhi0,idWhi1
     &           ,ndWcomp
     &           ,flattening
     &           ,iflatteninglo0,iflatteninglo1
     &           ,iflatteninghi0,iflatteninghi1
     &           ,numSlopes
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndWcomp
      integer idWlo0,idWlo1
      integer idWhi0,idWhi1
      REAL*8 dW(
     &           idWlo0:idWhi0,
     &           idWlo1:idWhi1,
     &           0:ndWcomp-1)
      integer iflatteninglo0,iflatteninglo1
      integer iflatteninghi0,iflatteninghi1
      REAL*8 flattening(
     &           iflatteninglo0:iflatteninghi0,
     &           iflatteninglo1:iflatteninghi1)
      integer numSlopes
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer lvar
      do lvar = 0,numSlopes - 1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            dW(i0,i1,lvar) = flattening(i0,i1)
     &        * dW(i0,i1,lvar)
      enddo
      enddo
      enddo
      return
      end
      subroutine VANLEERLIMITERF(
     &           dW
     &           ,idWlo0,idWlo1
     &           ,idWhi0,idWhi1
     &           ,ndWcomp
     &           ,dWleft
     &           ,idWleftlo0,idWleftlo1
     &           ,idWlefthi0,idWlefthi1
     &           ,ndWleftcomp
     &           ,dWright
     &           ,idWrightlo0,idWrightlo1
     &           ,idWrighthi0,idWrighthi1
     &           ,ndWrightcomp
     &           ,numslopes
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndWcomp
      integer idWlo0,idWlo1
      integer idWhi0,idWhi1
      REAL*8 dW(
     &           idWlo0:idWhi0,
     &           idWlo1:idWhi1,
     &           0:ndWcomp-1)
      integer ndWleftcomp
      integer idWleftlo0,idWleftlo1
      integer idWlefthi0,idWlefthi1
      REAL*8 dWleft(
     &           idWleftlo0:idWlefthi0,
     &           idWleftlo1:idWlefthi1,
     &           0:ndWleftcomp-1)
      integer ndWrightcomp
      integer idWrightlo0,idWrightlo1
      integer idWrighthi0,idWrighthi1
      REAL*8 dWright(
     &           idWrightlo0:idWrighthi0,
     &           idWrightlo1:idWrighthi1,
     &           0:ndWrightcomp-1)
      integer numslopes
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer iv
      REAL*8 dWl,dWr,dWc,dWlim
      do iv = 0,numslopes - 1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            dWc = dW     (i0,i1,iv)
            dWl = dWleft (i0,i1,iv)
            dWr = dWright(i0,i1,iv)
            dWlim = min((2.0d0)*abs(dWl),(2.0d0)*abs(dWr))
            dWlim = min(dWlim,abs(dWc))
            if (dWl * dWr .lt. 0) then
               dWlim = 0
            else
               dWlim = dWlim*sign((1.0d0),dWl)
            endif
            dW(i0,i1,iv) = dWlim
      enddo
      enddo
      enddo
      return
      end
      subroutine EXTPRESERVINGVANLEERLIMITERF(
     &           dW
     &           ,idWlo0,idWlo1
     &           ,idWhi0,idWhi1
     &           ,ndWcomp
     &           ,dWleft
     &           ,idWleftlo0,idWleftlo1
     &           ,idWlefthi0,idWlefthi1
     &           ,ndWleftcomp
     &           ,dWright
     &           ,idWrightlo0,idWrightlo1
     &           ,idWrighthi0,idWrighthi1
     &           ,ndWrightcomp
     &           ,numSlopes
     &           ,idir
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndWcomp
      integer idWlo0,idWlo1
      integer idWhi0,idWhi1
      REAL*8 dW(
     &           idWlo0:idWhi0,
     &           idWlo1:idWhi1,
     &           0:ndWcomp-1)
      integer ndWleftcomp
      integer idWleftlo0,idWleftlo1
      integer idWlefthi0,idWlefthi1
      REAL*8 dWleft(
     &           idWleftlo0:idWlefthi0,
     &           idWleftlo1:idWlefthi1,
     &           0:ndWleftcomp-1)
      integer ndWrightcomp
      integer idWrightlo0,idWrightlo1
      integer idWrighthi0,idWrighthi1
      REAL*8 dWright(
     &           idWrightlo0:idWrighthi0,
     &           idWrightlo1:idWrighthi1,
     &           0:ndWrightcomp-1)
      integer numSlopes
      integer idir
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      integer iv
      REAL*8 dWL, dWLL, dWR, dWRR, dWC, dWlim
      REAL*8 dp1, dp2, dpmin
      REAL*8 dW2L, dW2C, dW2R, dW2lim
      REAL*8 sign2, sign1
      REAL*8 cvl, dWvl
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      cvl = (5.0d0) * (0.250d0)
      do iv = 0,numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            dWC  = dW     (i0,i1, iv)
            dWL  = dWleft (i0,i1, iv)
            dWLL = dWleft (i0 -ioff0,i1 -ioff1, iv)
            dWR  = dWright(i0,i1, iv)
            dWRR = dWright(i0 +ioff0,i1 +ioff1, iv)
            dp1 = dWL * dWR
            dp2 = dWLL * dWRR
            dpmin = min(dp1, dp2)
            if (dpmin .lt. 0) then
               dW2L = dWL - dWLL
               dW2C = (dWR - dWL) * (0.500d0)
               dW2R = dWRR - dWR
               sign2 = sign((1.0d0), dW2C)
               dW2lim = min(abs(dW2C),
     &              max(sign2*dW2L, (0.0d0)),
     &              max(sign2*dW2R, (0.0d0)))
               dWvl = cvl * (3.0d0) * (0.500d0) * dW2lim
               if (sign2 * dWC .lt. 0.0) then
                  dWlim = min(dWvl, (2.0d0) * abs(dWL))
               else
                  dWlim = min(dWvl, (2.0d0) * abs(dWR))
               endif
            else
               dWlim = (2.0d0) * min(abs(dWL), abs(dWR))
            endif
            sign1 = sign((1.0d0), dWC)
            dW(i0,i1, iv) = sign1 * min(abs(dWC), dWlim)
      enddo
      enddo
       enddo
       return
       end
      subroutine PLMNORMALPREDF(
     &           dWMinus
     &           ,idWMinuslo0,idWMinuslo1
     &           ,idWMinushi0,idWMinushi1
     &           ,ndWMinuscomp
     &           ,dWPlus
     &           ,idWPluslo0,idWPluslo1
     &           ,idWPlushi0,idWPlushi1
     &           ,ndWPluscomp
     &           ,dW
     &           ,idWlo0,idWlo1
     &           ,idWhi0,idWhi1
     &           ,ndWcomp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1
     &           ,ilambdahi0,ilambdahi1
     &           ,nlambdacomp
     &           ,dtbydx
     &           ,nSlope
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndWMinuscomp
      integer idWMinuslo0,idWMinuslo1
      integer idWMinushi0,idWMinushi1
      REAL*8 dWMinus(
     &           idWMinuslo0:idWMinushi0,
     &           idWMinuslo1:idWMinushi1,
     &           0:ndWMinuscomp-1)
      integer ndWPluscomp
      integer idWPluslo0,idWPluslo1
      integer idWPlushi0,idWPlushi1
      REAL*8 dWPlus(
     &           idWPluslo0:idWPlushi0,
     &           idWPluslo1:idWPlushi1,
     &           0:ndWPluscomp-1)
      integer ndWcomp
      integer idWlo0,idWlo1
      integer idWhi0,idWhi1
      REAL*8 dW(
     &           idWlo0:idWhi0,
     &           idWlo1:idWhi1,
     &           0:ndWcomp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL*8 lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           0:nlambdacomp-1)
      REAL*8 dtbydx
      integer nSlope
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer iv
      REAL*8  lmin, lmax, lambdaK
      do iv = 0,nSlope-1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            lmin = min(lambda(i0,i1,0       ),(0.0d0))
            lmax = max(lambda(i0,i1,nSlope-1),(0.0d0))
            lambdaK = lambda(i0,i1,iv)
            if (lambdaK .gt. 0) then
               dWMinus(i0,i1,iv) = dW(i0,i1,iv) *
     &              (-(0.500d0)) * ((1.0d0) + dtbydx * lmin   )
               dWPlus (i0,i1,iv) = dW(i0,i1,iv) *
     &              (0.500d0) * ((1.0d0) - dtbydx * lambdaK)
            elseif (lambdaK .lt. 0) then
               dWMinus(i0,i1,iv) = dW(i0,i1,iv) *
     &              (-(0.500d0)) * ((1.0d0) + dtbydx * lambdaK)
               dWPlus (i0,i1,iv) = dW(i0,i1,iv) *
     &              (0.500d0) * ((1.0d0) - dtbydx * lmax   )
            else
               dWMinus(i0,i1,iv) = dW(i0,i1,iv) *
     &              (-(0.500d0)) * ((1.0d0) + dtbydx * lmin   )
               dWPlus (i0,i1,iv) = dW(i0,i1,iv) *
     &              (0.500d0) * ((1.0d0) - dtbydx * lmax   )
            endif
      enddo
      enddo
      enddo
      return
      end
      subroutine DIVUEDGEF(
     &           divu
     &           ,idivulo0,idivulo1
     &           ,idivuhi0,idivuhi1
     &           ,uNorm
     &           ,iuNormlo0,iuNormlo1
     &           ,iuNormhi0,iuNormhi1
     &           ,duTan
     &           ,iduTanlo0,iduTanlo1
     &           ,iduTanhi0,iduTanhi1
     &           ,nduTancomp
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivulo0,idivulo1
      integer idivuhi0,idivuhi1
      REAL*8 divu(
     &           idivulo0:idivuhi0,
     &           idivulo1:idivuhi1)
      integer iuNormlo0,iuNormlo1
      integer iuNormhi0,iuNormhi1
      REAL*8 uNorm(
     &           iuNormlo0:iuNormhi0,
     &           iuNormlo1:iuNormhi1)
      integer nduTancomp
      integer iduTanlo0,iduTanlo1
      integer iduTanhi0,iduTanhi1
      REAL*8 duTan(
     &           iduTanlo0:iduTanhi0,
     &           iduTanlo1:iduTanhi1,
     &           0:nduTancomp-1)
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      integer dimTan
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
         divu(i0,i1) = uNorm(i0 +f2cHi0,i1 +f2cHi1)
     &                       - uNorm(i0 +f2cLo0,i1 +f2cLo1)
         do dimTan = 0, 2 -2
            divu(i0,i1) = divu(i0,i1)
     &           + (0.500d0)*(duTan(i0 +f2cHi0,i1 +f2cHi1,dimTan)
     &                 + duTan(i0 +f2cLo0,i1 +f2cLo1,dimTan))
         enddo
      enddo
      enddo
      if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
            divu(i0,i1) = divu(i0 +ioff0,i1 +ioff1)
      enddo
      enddo
      endif
      if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
            divu(i0,i1) = divu(i0 -ioff0,i1 -ioff1)
      enddo
      enddo
      endif
      return
      end
      subroutine ARTVISCF(
     &           F
     &           ,iFlo0,iFlo1
     &           ,iFhi0,iFhi1
     &           ,nFcomp
     &           ,U
     &           ,iUlo0,iUlo1
     &           ,iUhi0,iUhi1
     &           ,nUcomp
     &           ,divu
     &           ,idivulo0,idivulo1
     &           ,idivuhi0,idivuhi1
     &           ,coeff
     &           ,idir
     &           ,ifboxlo0,ifboxlo1
     &           ,ifboxhi0,ifboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nFcomp
      integer iFlo0,iFlo1
      integer iFhi0,iFhi1
      REAL*8 F(
     &           iFlo0:iFhi0,
     &           iFlo1:iFhi1,
     &           0:nFcomp-1)
      integer nUcomp
      integer iUlo0,iUlo1
      integer iUhi0,iUhi1
      REAL*8 U(
     &           iUlo0:iUhi0,
     &           iUlo1:iUhi1,
     &           0:nUcomp-1)
      integer idivulo0,idivulo1
      integer idivuhi0,idivuhi1
      REAL*8 divu(
     &           idivulo0:idivuhi0,
     &           idivulo1:idivuhi1)
      REAL*8 coeff
      integer idir
      integer ifboxlo0,ifboxlo1
      integer ifboxhi0,ifboxhi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      integer iv
      REAL*8 dv, uLo, uHi
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do iv = 0, nUcomp - 1
      do i1 = ifboxlo1,ifboxhi1
      do i0 = ifboxlo0,ifboxhi0
            dv = divu(i0,i1)
            if (dv .lt. 0) then
               uLo = U(i0 +f2cLo0,i1 +f2cLo1, iv)
               uHi = U(i0 +f2cHi0,i1 +f2cHi1, iv)
               F(i0,i1, iv) = F(i0,i1, iv) +
     &              coeff * dv * (uHi - uLo)
            endif
      enddo
      enddo
      enddo
      return
      end
      subroutine PPMLIMITERF(
     &           dWMinus
     &           ,idWMinuslo0,idWMinuslo1
     &           ,idWMinushi0,idWMinushi1
     &           ,ndWMinuscomp
     &           ,dWPlus
     &           ,idWPluslo0,idWPluslo1
     &           ,idWPlushi0,idWPlushi1
     &           ,ndWPluscomp
     &           ,numSlopes
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndWMinuscomp
      integer idWMinuslo0,idWMinuslo1
      integer idWMinushi0,idWMinushi1
      REAL*8 dWMinus(
     &           idWMinuslo0:idWMinushi0,
     &           idWMinuslo1:idWMinushi1,
     &           0:ndWMinuscomp-1)
      integer ndWPluscomp
      integer idWPluslo0,idWPluslo1
      integer idWPlushi0,idWPlushi1
      REAL*8 dWPlus(
     &           idWPluslo0:idWPlushi0,
     &           idWPluslo1:idWPlushi1,
     &           0:ndWPluscomp-1)
      integer numSlopes
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer iv
      REAL*8 dWl,dWh,dWc,d2W,s
      do iv = 0, numSlopes - 1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            dWl = dWMinus(i0,i1, iv)
            dWh = dWPlus (i0,i1, iv)
            if (dWl * dWh .lt. 0) then
               dWc = (dWh + dWl) * (0.500d0)
               d2W =  dWh - dWl
               s = sign((1.0d0),dWc)
               if (dWc * d2W .gt. 0) then
                  dWPlus (i0,i1,iv) =
     &                 s * min(-(2.0d0) * s * dWl, s * dWh)
               else
                  dWMinus(i0,i1,iv) =
     &                 s * min(s * dWl, -(2.0d0) * s * dWh)
               endif
            else
               dWPlus (i0,i1, iv) = 0
               dWMinus(i0,i1, iv) = 0
            endif
      enddo
      enddo
      enddo
      return
      end
      subroutine COLELLASEKORALIMITERF(
     &           dWMinus
     &           ,idWMinuslo0,idWMinuslo1
     &           ,idWMinushi0,idWMinushi1
     &           ,ndWMinuscomp
     &           ,dWPlus
     &           ,idWPluslo0,idWPluslo1
     &           ,idWPlushi0,idWPlushi1
     &           ,ndWPluscomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,diffW
     &           ,idiffWlo0,idiffWlo1
     &           ,idiffWhi0,idiffWhi1
     &           ,ndiffWcomp
     &           ,d2W
     &           ,id2Wlo0,id2Wlo1
     &           ,id2Whi0,id2Whi1
     &           ,nd2Wcomp
     &           ,dW2fcf
     &           ,idW2fcflo0,idW2fcflo1
     &           ,idW2fcfhi0,idW2fcfhi1
     &           ,ndW2fcfcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           ,limitC
     &           ,eps
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndWMinuscomp
      integer idWMinuslo0,idWMinuslo1
      integer idWMinushi0,idWMinushi1
      REAL*8 dWMinus(
     &           idWMinuslo0:idWMinushi0,
     &           idWMinuslo1:idWMinushi1,
     &           0:ndWMinuscomp-1)
      integer ndWPluscomp
      integer idWPluslo0,idWPluslo1
      integer idWPlushi0,idWPlushi1
      REAL*8 dWPlus(
     &           idWPluslo0:idWPlushi0,
     &           idWPluslo1:idWPlushi1,
     &           0:ndWPluscomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer ndiffWcomp
      integer idiffWlo0,idiffWlo1
      integer idiffWhi0,idiffWhi1
      REAL*8 diffW(
     &           idiffWlo0:idiffWhi0,
     &           idiffWlo1:idiffWhi1,
     &           0:ndiffWcomp-1)
      integer nd2Wcomp
      integer id2Wlo0,id2Wlo1
      integer id2Whi0,id2Whi1
      REAL*8 d2W(
     &           id2Wlo0:id2Whi0,
     &           id2Wlo1:id2Whi1,
     &           0:nd2Wcomp-1)
      integer ndW2fcfcomp
      integer idW2fcflo0,idW2fcflo1
      integer idW2fcfhi0,idW2fcfhi1
      REAL*8 dW2fcf(
     &           idW2fcflo0:idW2fcfhi0,
     &           idW2fcflo1:idW2fcfhi1,
     &           0:ndW2fcfcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      REAL*8 limitC
      REAL*8 eps
      integer i0,i1
      integer ioff0,ioff1
      integer iv
      REAL*8 WC, dWM, dWP, WL, WR, dWL, dWR
      REAL*8 atfcf, d2Wlim, rho
      REAL*8 d2WL, d2WC, d2WR
      REAL*8 sd2WL, sd2WC, sd2WR, sd2fcf
      REAL*8 diffWL, diffWR, dWminA, dWminF
      logical extrem, bigM, bigP
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do iv = 0,numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            dWM = dWMinus(i0,i1, iv)
            dWP = dWPlus (i0,i1, iv)
            extrem = .false.
            if (dWM * dWP .ge. 0) then
               extrem = .true.
            else
               bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
               bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
               if (bigM .or. bigP) then
                  WL  = W(i0 -ioff0,i1 -ioff1, iv)
                  WC  = W(i0,i1, iv)
                  WR  = W(i0 +ioff0,i1 +ioff1, iv)
                  dWL = WC - WL
                  dWR = WR - WC
                  diffWL = diffW(i0 -ioff0,i1 -ioff1, iv)
                  diffWR = diffW(i0 +ioff0,i1 +ioff1, iv)
                  dWminA = min(abs(dWL), abs(dWR))
                  dWminF = min(abs(diffWL), abs(diffWR))
                  if ( ( (dWminA .ge. dWminF) .and.
     &                 (dWL * dWR .le. 0) ) .or.
     &                 ( (dWminF .gt. dWminA) .and.
     &                 (diffWL * diffWR .le. 0) ) ) then
                     extrem = .true.
                  else
                     if (bigM) then
                        dWMinus(i0,i1, iv) = -(2.0d0) * dWP
                     endif
                     if (bigP) then
                        dWPlus (i0,i1, iv) = -(2.0d0) * dWM
                     endif
                  endif
               endif
            endif
            if (extrem) then
               d2WL  = d2W(i0 -ioff0,i1 -ioff1, iv)
               d2WC  = d2W(i0,i1, iv)
               d2WR  = d2W(i0 +ioff0,i1 +ioff1, iv)
               if ((d2WC - d2WL) * (d2WR - d2WC) .lt. 0) then
                  atfcf = dW2fcf(i0,i1, iv)
                  if (abs(atfcf) .lt. eps) then
                     rho = 0
                  else
                     sd2WL = sign((1.0d0), d2WL)
                     sd2WC = sign((1.0d0), d2WC)
                     sd2WR = sign((1.0d0), d2WR)
                     sd2fcf = sign((1.0d0), atfcf)
                     if ((sd2WC .eq. sd2WL) .and.
     &                    (sd2WC .eq. sd2WR) .and.
     &                    (sd2WC .eq. sd2fcf)) then
                        d2Wlim = sd2WC *
     &                       min(abs(atfcf),
     &                       limitC * abs(d2WC),
     &                       limitC * abs(d2WL),
     &                       limitC * abs(d2WR))
                        rho = d2Wlim / atfcf
                     else
                        rho = 0
                     endif
                  endif
                  dWPlus (i0,i1, iv) = dWP * rho
                  dWMinus(i0,i1, iv) = dWM * rho
               endif
            endif
      enddo
      enddo
         if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               dWM = dWMinus(i0,i1, iv)
               dWP = dWPlus (i0,i1, iv)
               extrem = .false.
               if (dWM * dWP .ge. 0) then
                  extrem = .true.
               else
                  bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
                  bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
                  if (bigM) then
                     dWMinus(i0,i1, iv) = -(2.0d0) * dWP
                  endif
                  if (bigP) then
                     dWPlus (i0,i1, iv) = -(2.0d0) * dWM
                  endif
               endif
      enddo
      enddo
         endif
         if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               dWM = dWMinus(i0,i1, iv)
               dWP = dWPlus (i0,i1, iv)
               extrem = .false.
               if (dWM * dWP .ge. 0) then
                  extrem = .true.
               else
                  bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
                  bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
                  if (bigM) then
                     dWMinus(i0,i1, iv) = -(2.0d0) * dWP
                  endif
                  if (bigP) then
                     dWPlus (i0,i1, iv) = -(2.0d0) * dWM
                  endif
               endif
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine PPMFACEVALUESF(
     &           WFace
     &           ,iWFacelo0,iWFacelo1
     &           ,iWFacehi0,iWFacehi1
     &           ,nWFacecomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,dW
     &           ,idWlo0,idWlo1
     &           ,idWhi0,idWhi1
     &           ,ndWcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nWFacecomp
      integer iWFacelo0,iWFacelo1
      integer iWFacehi0,iWFacehi1
      REAL*8 WFace(
     &           iWFacelo0:iWFacehi0,
     &           iWFacelo1:iWFacehi1,
     &           0:nWFacecomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer ndWcomp
      integer idWlo0,idWlo1
      integer idWhi0,idWhi1
      REAL*8 dW(
     &           idWlo0:idWhi0,
     &           idWlo1:idWhi1,
     &           0:ndWcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      integer iv
      REAL*8 WLeft,WRight
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do iv = 0,numSlopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            WRight =  W(i0 +f2cHi0,i1 +f2cHi1, iv)
     &        - dW(i0 +f2cHi0,i1 +f2cHi1,iv) / (3.0d0)
            WLeft  =  W(i0 +f2cLo0,i1 +f2cLo1, iv)
     &           + dW(i0 +f2cLo0,i1 +f2cLo1, iv) / (3.0d0)
            WFace(i0,i1, iv) = (WLeft + WRight) * (0.500d0)
      enddo
      enddo
      enddo
      if (hasLo .eq. 1) then
         do iv = 0,numSlopes-1
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               WRight =  W(i0 +f2cHi0,i1 +f2cHi1,iv)
     &           - dW(i0 +f2cHi0,i1 +f2cHi1,iv) * (0.500d0)
               WFace(i0,i1, iv) = WRight
      enddo
      enddo
         enddo
      endif
      if (hasHi .eq. 1) then
         do iv = 0,numSlopes-1
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               WLeft =  W(i0 +f2cLo0,i1 +f2cLo1, iv)
     &           + dW(i0 +f2cLo0,i1 +f2cLo1, iv) * (0.500d0)
               WFace(i0,i1,iv) = WLeft
      enddo
      enddo
         enddo
      endif
      return
      end
      subroutine PPMNORMALPREDF(
     &           dWMinus
     &           ,idWMinuslo0,idWMinuslo1
     &           ,idWMinushi0,idWMinushi1
     &           ,ndWMinuscomp
     &           ,dWPlus
     &           ,idWPluslo0,idWPluslo1
     &           ,idWPlushi0,idWPlushi1
     &           ,ndWPluscomp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1
     &           ,ilambdahi0,ilambdahi1
     &           ,nlambdacomp
     &           ,dtbydx
     &           ,nSlope
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndWMinuscomp
      integer idWMinuslo0,idWMinuslo1
      integer idWMinushi0,idWMinushi1
      REAL*8 dWMinus(
     &           idWMinuslo0:idWMinushi0,
     &           idWMinuslo1:idWMinushi1,
     &           0:ndWMinuscomp-1)
      integer ndWPluscomp
      integer idWPluslo0,idWPluslo1
      integer idWPlushi0,idWPlushi1
      REAL*8 dWPlus(
     &           idWPluslo0:idWPlushi0,
     &           idWPluslo1:idWPlushi1,
     &           0:ndWPluscomp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL*8 lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           0:nlambdacomp-1)
      REAL*8 dtbydx
      integer nSlope
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer iv
      REAL*8 sigMinus,sigPlus,sigmin,sigmax,lambdaK,dWl,dWh
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         sigmin = lambda(i0,i1, 0)*dtbydx
         sigmin = min(sigmin,(0.0d0))
         sigmin = -sigmin
         sigmax = lambda(i0,i1, nSlope-1)*dtbydx
         sigmax = max(sigmax,(0.0d0))
         do iv = 0,nSlope - 1
            lambdaK = lambda(i0,i1, iv)
            if (lambdaK .gt. 0) then
               sigMinus = sigmin
               sigPlus  = lambdaK*dtbydx
            else
               sigMinus = -lambdaK*dtbydx
               sigPlus  = sigmax
            endif
            dWl = dWMinus(i0,i1,iv)
            dWh = dWPlus (i0,i1,iv)
            dWMinus(i0,i1, iv) =
     &           dWl + sigMinus  * ((dWh - dWl)
     &           - (dWh + dWl) * ((3.0d0) - (2.0d0)*sigMinus)) * (0.500d
     &0)
            dWPlus (i0,i1, iv) =
     &           dWh + sigPlus * ((dWl - dWh)
     &           - (dWh + dWl) * ((3.0d0) - (2.0d0)*sigPlus )) * (0.500d
     &0)
         enddo
      enddo
      enddo
      return
      end
      subroutine GETSECONDDIFF(
     &           d2W
     &           ,id2Wlo0,id2Wlo1
     &           ,id2Whi0,id2Whi1
     &           ,nd2Wcomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nd2Wcomp
      integer id2Wlo0,id2Wlo1
      integer id2Whi0,id2Whi1
      REAL*8 d2W(
     &           id2Wlo0:id2Whi0,
     &           id2Wlo1:id2Whi1,
     &           0:nd2Wcomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      integer iv
      REAL*8 WC, WL, WR, dWL, dWR
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do iv = 0, numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            WL  = W(i0 -ioff0,i1 -ioff1, iv)
            WC  = W(i0,i1, iv)
            WR  = W(i0 +ioff0,i1 +ioff1, iv)
            dWL = WC - WL
            dWR = WR - WC
            d2W(i0,i1, iv) = dWR - dWL
      enddo
      enddo
         if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               d2W(i0,i1, iv) =
     &           d2W(i0 +ioff0,i1 +ioff1, iv)
      enddo
      enddo
         endif
         if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               d2W(i0,i1, iv) =
     &           d2W(i0 -ioff0,i1 -ioff1, iv)
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine GETD2CELL(
     &           d2W
     &           ,id2Wlo0,id2Wlo1
     &           ,id2Whi0,id2Whi1
     &           ,nd2Wcomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nd2Wcomp
      integer id2Wlo0,id2Wlo1
      integer id2Whi0,id2Whi1
      REAL*8 d2W(
     &           id2Wlo0:id2Whi0,
     &           id2Wlo1:id2Whi1,
     &           0:nd2Wcomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer ioff0,ioff1
      integer lvar
      REAL*8 WL, WR, WC
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do lvar = 0, numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            WL = W(i0 -ioff0,i1 -ioff1, lvar)
            WC = W(i0,i1, lvar)
            WR = W(i0 +ioff0,i1 +ioff1, lvar)
            d2W(i0,i1, lvar) = (WR - WC) - (WC - WL)
      enddo
      enddo
         if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               d2W(i0,i1, lvar) =
     &           d2W(i0 +ioff0,i1 +ioff1, lvar)
      enddo
      enddo
         endif
         if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               d2W(i0,i1, lvar) =
     &           d2W(i0 -ioff0,i1 -ioff1, lvar)
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine GETD2LIMFACE(
     &           d2Wlim
     &           ,id2Wlimlo0,id2Wlimlo1
     &           ,id2Wlimhi0,id2Wlimhi1
     &           ,nd2Wlimcomp
     &           ,d2Wc
     &           ,id2Wclo0,id2Wclo1
     &           ,id2Wchi0,id2Wchi1
     &           ,nd2Wccomp
     &           ,d2Wcfc
     &           ,id2Wcfclo0,id2Wcfclo1
     &           ,id2Wcfchi0,id2Wcfchi1
     &           ,nd2Wcfccomp
     &           ,numSlopes
     &           ,idir
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           ,limitC
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nd2Wlimcomp
      integer id2Wlimlo0,id2Wlimlo1
      integer id2Wlimhi0,id2Wlimhi1
      REAL*8 d2Wlim(
     &           id2Wlimlo0:id2Wlimhi0,
     &           id2Wlimlo1:id2Wlimhi1,
     &           0:nd2Wlimcomp-1)
      integer nd2Wccomp
      integer id2Wclo0,id2Wclo1
      integer id2Wchi0,id2Wchi1
      REAL*8 d2Wc(
     &           id2Wclo0:id2Wchi0,
     &           id2Wclo1:id2Wchi1,
     &           0:nd2Wccomp-1)
      integer nd2Wcfccomp
      integer id2Wcfclo0,id2Wcfclo1
      integer id2Wcfchi0,id2Wcfchi1
      REAL*8 d2Wcfc(
     &           id2Wcfclo0:id2Wcfchi0,
     &           id2Wcfclo1:id2Wcfchi1,
     &           0:nd2Wcfccomp-1)
      integer numSlopes
      integer idir
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      REAL*8 limitC
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      integer lvar
      REAL*8 d2WL, d2WR, atcfc
      REAL*8 sd2WL, sd2WR, sd2cfc
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do lvar = 0, numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            d2WL  = d2Wc(i0 +f2cLo0,i1 +f2cLo1, lvar)
            d2WR  = d2Wc(i0 +f2cHi0,i1 +f2cHi1, lvar)
            atcfc = d2Wcfc(i0,i1, lvar)
            sd2WL = sign((1.0d0), d2WL)
            sd2WR = sign((1.0d0), d2WR)
            sd2cfc = sign((1.0d0), atcfc)
            if ((sd2WL .eq. sd2cfc) .and. (sd2WR .eq. sd2cfc)) then
               d2Wlim(i0,i1, lvar) = sd2cfc *
     &              min(limitC * min(abs(d2WL), abs(d2WR)),
     &              abs(atcfc))
            else
               d2Wlim(i0,i1, lvar) = 0
            endif
      enddo
      enddo
      enddo
      return
      end
      subroutine COLELLASEKORAFACELIMITER(
     &           Wface
     &           ,iWfacelo0,iWfacelo1
     &           ,iWfacehi0,iWfacehi1
     &           ,nWfacecomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,d2Wc
     &           ,id2Wclo0,id2Wclo1
     &           ,id2Wchi0,id2Wchi1
     &           ,nd2Wccomp
     &           ,numSlopes
     &           ,idir
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           ,limitC
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nWfacecomp
      integer iWfacelo0,iWfacelo1
      integer iWfacehi0,iWfacehi1
      REAL*8 Wface(
     &           iWfacelo0:iWfacehi0,
     &           iWfacelo1:iWfacehi1,
     &           0:nWfacecomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer nd2Wccomp
      integer id2Wclo0,id2Wclo1
      integer id2Wchi0,id2Wchi1
      REAL*8 d2Wc(
     &           id2Wclo0:id2Wchi0,
     &           id2Wclo1:id2Wchi1,
     &           0:nd2Wccomp-1)
      integer numSlopes
      integer idir
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      REAL*8 limitC
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      integer lvar
      REAL*8 Wf, WL, WR, d2Wlim, dp
      REAL*8 dWL, dWR, d2WL, d2WR, d2Wcfc
      REAL*8 sd2WL, sd2WR, sd2cfc
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do lvar = 0, numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            Wf = Wface(i0,i1, lvar)
            WL = W(i0 +f2cLo0,i1 +f2cLo1, lvar)
            WR = W(i0 +f2cHi0,i1 +f2cHi1, lvar)
            dWR = WR - Wf
            dWL = Wf - WL
            dp = dWR * dWL
            if (dp .lt. 0) then
               d2WL  = d2Wc(i0 +f2cLo0,i1 +f2cLo1, lvar)
               d2WR  = d2Wc(i0 +f2cHi0,i1 +f2cHi1, lvar)
               d2Wcfc = (3.0d0) * (dWR - dWL)
               if ((d2Wcfc - d2WL) * (d2WR - d2Wcfc) .lt. 0) then
                  sd2WL = sign((1.0d0), d2WL)
                  sd2WR = sign((1.0d0), d2WR)
                  sd2cfc = sign((1.0d0), d2Wcfc)
                  if ((sd2WL .eq. sd2cfc) .and.
     &                 (sd2WR .eq. sd2cfc)) then
                     d2Wlim = sd2cfc *
     &                    min(limitC * min(abs(d2WL), abs(d2WR)),
     &                    abs(d2Wcfc))
                  else
                     d2Wlim = 0
                  endif
                  Wface(i0,i1, lvar) =
     &                 (WL + WR) * (0.500d0) - d2Wlim / (6.0d0)
               endif
            endif
      enddo
      enddo
      enddo
      return
      end
      subroutine FOURTHINTERPFACES(
     &           Wface
     &           ,iWfacelo0,iWfacelo1
     &           ,iWfacehi0,iWfacehi1
     &           ,nWfacecomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,inextLoBoxlo0,inextLoBoxlo1
     &           ,inextLoBoxhi0,inextLoBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,inextHiBoxlo0,inextHiBoxlo1
     &           ,inextHiBoxhi0,inextHiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nWfacecomp
      integer iWfacelo0,iWfacelo1
      integer iWfacehi0,iWfacehi1
      REAL*8 Wface(
     &           iWfacelo0:iWfacehi0,
     &           iWfacelo1:iWfacehi1,
     &           0:nWfacecomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer inextLoBoxlo0,inextLoBoxlo1
      integer inextLoBoxhi0,inextLoBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer inextHiBoxlo0,inextHiBoxlo1
      integer inextHiBoxhi0,inextHiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2c1Lo0,f2c1Lo1
      integer f2c2Lo0,f2c2Lo1
      integer f2c3Lo0,f2c3Lo1
      integer f2cHi0,f2cHi1
      integer f2c1Hi0,f2c1Hi1
      integer f2c2Hi0,f2c2Hi1
      integer f2c3Hi0,f2c3Hi1
      integer lvar
      REAL*8 WL, WR, WLL, WRR, WLLL, WRRR, WLLLL, WRRRR
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2c1Lo0= -2*CHF_ID(0, idir)
      f2c1Lo1= -2*CHF_ID(1, idir)
      f2c2Lo0= -3*CHF_ID(0, idir)
      f2c2Lo1= -3*CHF_ID(1, idir)
      f2c3Lo0= -4*CHF_ID(0, idir)
      f2c3Lo1= -4*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      f2c1Hi0= 1*CHF_ID(0, idir)
      f2c1Hi1= 1*CHF_ID(1, idir)
      f2c2Hi0= 2*CHF_ID(0, idir)
      f2c2Hi1= 2*CHF_ID(1, idir)
      f2c3Hi0= 3*CHF_ID(0, idir)
      f2c3Hi1= 3*CHF_ID(1, idir)
      do lvar = 0, numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            WLL = W(i0 +f2c1Lo0,i1 +f2c1Lo1, lvar)
            WL  = W(i0 +f2cLo0,i1 +f2cLo1,  lvar)
            WR  = W(i0 +f2cHi0,i1 +f2cHi1,  lvar)
            WRR = W(i0 +f2c1Hi0,i1 +f2c1Hi1, lvar)
            Wface(i0,i1, lvar) =
     &           ((7.0d0) * (WL + WR) - (WLL + WRR)) / (12.0d0)
      enddo
      enddo
         if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               WR    = W(i0 +f2cHi0,i1 +f2cHi1, lvar)
               WRR   = W(i0 +f2c1Hi0,i1 +f2c1Hi1, lvar)
               WRRR  = W(i0 +f2c2Hi0,i1 +f2c2Hi1, lvar)
               WRRRR = W(i0 +f2c3Hi0,i1 +f2c3Hi1, lvar)
               Wface(i0,i1, lvar) =
     &              ((((20.0d0) + (5.0d0))*WR - ((20.0d0) + (3.0d0))*WRR
     &)
     &              + (((10.0d0) + (3.0d0))*WRRR - (3.0d0)*WRRRR)) / (12
     &.0d0)
      enddo
      enddo
      do i1 = inextLoBoxlo1,inextLoBoxhi1
      do i0 = inextLoBoxlo0,inextLoBoxhi0
               WL   = W(i0 +f2cLo0,i1 +f2cLo1,  lvar)
               WR   = W(i0 +f2cHi0,i1 +f2cHi1,  lvar)
               WRR  = W(i0 +f2c1Hi0,i1 +f2c1Hi1, lvar)
               WRRR = W(i0 +f2c2Hi0,i1 +f2c2Hi1, lvar)
               Wface(i0,i1, lvar) =
     &              ((3.0d0)*WL + ((10.0d0) + (3.0d0))*WR - (5.0d0)*WRR 
     &+ WRRR) /
     &              (12.0d0)
      enddo
      enddo
         endif
         if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               WL    = W(i0 +f2cLo0,i1 +f2cLo1,  lvar)
               WLL   = W(i0 +f2c1Lo0,i1 +f2c1Lo1, lvar)
               WLLL  = W(i0 +f2c2Lo0,i1 +f2c2Lo1, lvar)
               WLLLL = W(i0 +f2c3Lo0,i1 +f2c3Lo1, lvar)
               Wface(i0,i1, lvar) =
     &              ((((20.0d0) + (5.0d0))*WL - ((20.0d0) + (3.0d0))*WLL
     &)
     &              + (((10.0d0) + (3.0d0))*WLLL - (3.0d0)*WLLLL)) / (12
     &.0d0)
      enddo
      enddo
      do i1 = inextHiBoxlo1,inextHiBoxhi1
      do i0 = inextHiBoxlo0,inextHiBoxhi0
               WR   = W(i0 +f2cHi0,i1 +f2cHi1,  lvar)
               WL   = W(i0 +f2cLo0,i1 +f2cLo1,  lvar)
               WLL  = W(i0 +f2c1Lo0,i1 +f2c1Lo1, lvar)
               WLLL = W(i0 +f2c2Lo0,i1 +f2c2Lo1, lvar)
               Wface(i0,i1, lvar) =
     &              ((3.0d0)*WR + ((10.0d0) + (3.0d0))*WL - (5.0d0)*WLL 
     &+ WLLL) /
     &                (12.0d0)
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine SECONDINTERPFACES(
     &           Wface
     &           ,iWfacelo0,iWfacelo1
     &           ,iWfacehi0,iWfacehi1
     &           ,nWfacecomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nWfacecomp
      integer iWfacelo0,iWfacelo1
      integer iWfacehi0,iWfacehi1
      REAL*8 Wface(
     &           iWfacelo0:iWfacehi0,
     &           iWfacelo1:iWfacehi1,
     &           0:nWfacecomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2c1Lo0,f2c1Lo1
      integer f2c2Lo0,f2c2Lo1
      integer f2c3Lo0,f2c3Lo1
      integer f2cHi0,f2cHi1
      integer f2c1Hi0,f2c1Hi1
      integer f2c2Hi0,f2c2Hi1
      integer f2c3Hi0,f2c3Hi1
      integer lvar
      REAL*8 WL, WR, WLL, WRR, WLLL, WRRR, WLLLL, WRRRR
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2c1Lo0= -2*CHF_ID(0, idir)
      f2c1Lo1= -2*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      f2c1Hi0= 1*CHF_ID(0, idir)
      f2c1Hi1= 1*CHF_ID(1, idir)
      do lvar = 0, numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            WLL = W(i0 +f2c1Lo0,i1 +f2c1Lo1, lvar)
            WL  = W(i0 +f2cLo0,i1 +f2cLo1,  lvar)
            WR  = W(i0 +f2cHi0,i1 +f2cHi1,  lvar)
            WRR = W(i0 +f2c1Hi0,i1 +f2c1Hi1, lvar)
            Wface(i0,i1, lvar) = (WL + WR) / (2.0d0)
      enddo
      enddo
         if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               WR  = W(i0 +f2cHi0,i1 +f2cHi1, lvar)
               WRR = W(i0 +f2c1Hi0,i1 +f2c1Hi1, lvar)
               Wface(i0,i1, lvar) = ((3.0d0)*WR - WRR) / (2.0d0)
      enddo
      enddo
         endif
         if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               WL  = W(i0 +f2cLo0,i1 +f2cLo1,  lvar)
               WLL = W(i0 +f2c1Lo0,i1 +f2c1Lo1, lvar)
               Wface(i0,i1, lvar) = ((3.0d0)*WL - WLL) / (2.0d0)
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine PPMFOURTHFACE(
     &           Wface
     &           ,iWfacelo0,iWfacelo1
     &           ,iWfacehi0,iWfacehi1
     &           ,nWfacecomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,d2W
     &           ,id2Wlo0,id2Wlo1
     &           ,id2Whi0,id2Whi1
     &           ,nd2Wcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nWfacecomp
      integer iWfacelo0,iWfacelo1
      integer iWfacehi0,iWfacehi1
      REAL*8 Wface(
     &           iWfacelo0:iWfacehi0,
     &           iWfacelo1:iWfacehi1,
     &           0:nWfacecomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer nd2Wcomp
      integer id2Wlo0,id2Wlo1
      integer id2Whi0,id2Whi1
      REAL*8 d2W(
     &           id2Wlo0:id2Whi0,
     &           id2Wlo1:id2Whi1,
     &           0:nd2Wcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2c1Lo0,f2c1Lo1
      integer f2c2Lo0,f2c2Lo1
      integer f2c3Lo0,f2c3Lo1
      integer f2cHi0,f2cHi1
      integer f2c1Hi0,f2c1Hi1
      integer f2c2Hi0,f2c2Hi1
      integer f2c3Hi0,f2c3Hi1
      integer iv
      REAL*8 WL, WLL, WLLL, WLLLL
      REAL*8 WR, WRR, WRRR, WRRRR
      REAL*8 d2WL, d2WR
      REAL*8 c0, c1, c2, c3
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2c1Lo0= -2*CHF_ID(0, idir)
      f2c1Lo1= -2*CHF_ID(1, idir)
      f2c2Lo0= -3*CHF_ID(0, idir)
      f2c2Lo1= -3*CHF_ID(1, idir)
      f2c3Lo0= -4*CHF_ID(0, idir)
      f2c3Lo1= -4*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      f2c1Hi0= 1*CHF_ID(0, idir)
      f2c1Hi1= 1*CHF_ID(1, idir)
      f2c2Hi0= 2*CHF_ID(0, idir)
      f2c2Hi1= 2*CHF_ID(1, idir)
      f2c3Hi0= 3*CHF_ID(0, idir)
      f2c3Hi1= 3*CHF_ID(1, idir)
      c0 = (20.0d0) + (7.0d0)
      c1 = -((20.0d0) + (8.0d0))
      c2 = (10.0d0) + (7.0d0)
      c3 = -(4.0d0)
      do iv = 0, numslopes - 1
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
            WL =     W(i0 +f2cLo0,i1 +f2cLo1, iv)
            WR =     W(i0 +f2cHi0,i1 +f2cHi1, iv)
            d2WL = d2W(i0 +f2cLo0,i1 +f2cLo1, iv)
            d2WR = d2W(i0 +f2cHi0,i1 +f2cHi1, iv)
            Wface(i0,i1, iv) = (WL + WR) * (0.500d0)
     &           - (d2WL + d2WR) / (12.0d0)
      enddo
      enddo
         if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               d2WR = d2W(i0 +f2cHi0,i1 +f2cHi1, iv)
               WR   =   W(i0 +f2cHi0,i1 +f2cHi1, iv)
               WRR  =   W(i0 +f2c1Hi0,i1 +f2c1Hi1, iv)
               WRRR =   W(i0 +f2c2Hi0,i1 +f2c2Hi1, iv)
               WRRRR =  W(i0 +f2c3Hi0,i1 +f2c3Hi1, iv)
               Wface(i0,i1, iv) =
     &              (c0*WR + c1*WRR + c2*WRRR + c3*WRRRR - d2WR) /
     &              (12.0d0)
      enddo
      enddo
         endif
         if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               d2WL = d2W(i0 +f2cLo0,i1 +f2cLo1, iv)
               WL   =   W(i0 +f2cLo0,i1 +f2cLo1, iv)
               WLL  =   W(i0 +f2c1Lo0,i1 +f2c1Lo1, iv)
               WLLL =   W(i0 +f2c2Lo0,i1 +f2c2Lo1, iv)
               WLLLL =  W(i0 +f2c3Lo0,i1 +f2c3Lo1, iv)
               Wface(i0,i1, iv) =
     &              (c0*WL + c1*WLL + c2*WLLL + c3*WLLLL - d2WL) /
     &              (12.0d0)
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine HODIVCOEF(
     &           divVel
     &           ,idivVello0,idivVello1
     &           ,idivVelhi0,idivVelhi1
     &           ,csq
     &           ,icsqlo0,icsqlo1
     &           ,icsqhi0,icsqhi1
     &           ,dir
     &           ,M0sq
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivVello0,idivVello1
      integer idivVelhi0,idivVelhi1
      REAL*8 divVel(
     &           idivVello0:idivVelhi0,
     &           idivVello1:idivVelhi1)
      integer icsqlo0,icsqlo1
      integer icsqhi0,icsqhi1
      REAL*8 csq(
     &           icsqlo0:icsqhi0,
     &           icsqlo1:icsqhi1)
      integer dir
      REAL*8 M0sq
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      integer c2fLo0,c2fLo1
      integer c2fHi0,c2fHi1
      REAL*8 csqmin, dvold
      f2cLo0= -1*CHF_ID(0, dir)
      f2cLo1= -1*CHF_ID(1, dir)
      f2cHi0= 0*CHF_ID(0, dir)
      f2cHi1= 0*CHF_ID(1, dir)
      c2fLo0= 0*CHF_ID(0, dir)
      c2fLo1= 0*CHF_ID(1, dir)
      c2fHi0= 1*CHF_ID(0, dir)
      c2fHi1= 1*CHF_ID(1, dir)
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
         csqmin = min(csq(i0 +f2cLo0,i1 +f2cLo1),
     &                csq(i0 +f2cHi0,i1 +f2cHi1) )
         dvold = divVel(i0,i1)
         divVel(i0,i1) = dvold *
     &        min(dvold**2/csqmin/M0sq, (1.0d0))
      enddo
      enddo
      if (hasLo .ne. 0) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
            csqmin = csq(i0,i1)
            dvold = divVel(i0 +c2fLo0,i1 +c2fLo1)
            divVel(i0 +c2fLo0,i1 +c2fLo1) = dvold *
     &           min(dvold**2/csqmin/M0sq, (1.0d0))
      enddo
      enddo
      endif
      if (hasHi .ne. 0) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
            csqmin = csq(i0,i1)
            dvold = divVel(i0 +c2fHi0,i1 +c2fHi1)
            divVel(i0 +c2fHi0,i1 +c2fHi1) = dvold *
     &           min(dvold**2/csqmin/M0sq, (1.0d0))
      enddo
      enddo
      endif
      return
      end
      subroutine HIGHORDERDIVCO(
     &           divVel
     &           ,idivVello0,idivVello1
     &           ,idivVelhi0,idivVelhi1
     &           ,csq
     &           ,icsqlo0,icsqlo1
     &           ,icsqhi0,icsqhi1
     &           ,idir
     &           ,M0sq
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1
     &           ,icenterBoxhi0,icenterBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivVello0,idivVello1
      integer idivVelhi0,idivVelhi1
      REAL*8 divVel(
     &           idivVello0:idivVelhi0,
     &           idivVello1:idivVelhi1)
      integer icsqlo0,icsqlo1
      integer icsqhi0,icsqhi1
      REAL*8 csq(
     &           icsqlo0:icsqhi0,
     &           icsqlo1:icsqhi1)
      integer idir
      REAL*8 M0sq
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      REAL*8 csqmin, dvold, fac
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
         csqmin = min(csq(i0 +f2cLo0,i1 +f2cLo1),
     &                csq(i0 +f2cHi0,i1 +f2cHi1))
         dvold = divVel(i0,i1)
         fac = dvold**2/csqmin/M0sq
         if (fac .lt. (1.0d0)) then
            divVel(i0,i1) = dvold * fac
         endif
      enddo
      enddo
      if (hasLo .ne. 0) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
            csqmin = csq(i0 +f2cHi0,i1 +f2cHi1)
            dvold = divVel(i0,i1)
            fac = dvold**2/csqmin/M0sq
            if (fac .lt. (1.0d0)) then
               divVel(i0,i1) = dvold * fac
            endif
      enddo
      enddo
      endif
      if (hasHi .ne. 0) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
            csqmin = csq(i0 +f2cLo0,i1 +f2cLo1)
            dvold = divVel(i0,i1)
            fac = dvold**2/csqmin/M0sq
            if (fac .lt. (1.0d0)) then
               divVel(i0,i1) = dvold * fac
            endif
      enddo
      enddo
      endif
      return
      end
      subroutine CHECKCUBICLIMITERF(
     &           dWMinus
     &           ,idWMinuslo0,idWMinuslo1
     &           ,idWMinushi0,idWMinushi1
     &           ,ndWMinuscomp
     &           ,dWPlus
     &           ,idWPluslo0,idWPluslo1
     &           ,idWPlushi0,idWPlushi1
     &           ,ndWPluscomp
     &           ,W
     &           ,iWlo0,iWlo1
     &           ,iWhi0,iWhi1
     &           ,nWcomp
     &           ,d2W
     &           ,id2Wlo0,id2Wlo1
     &           ,id2Whi0,id2Whi1
     &           ,nd2Wcomp
     &           ,dW2fcf
     &           ,idW2fcflo0,idW2fcflo1
     &           ,idW2fcfhi0,idW2fcfhi1
     &           ,ndW2fcfcomp
     &           ,numSlopes
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1
     &           ,iloBoxhi0,iloBoxhi1
     &           ,inextLoBoxlo0,inextLoBoxlo1
     &           ,inextLoBoxhi0,inextLoBoxhi1
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1
     &           ,ihiBoxhi0,ihiBoxhi1
     &           ,inextHiBoxlo0,inextHiBoxlo1
     &           ,inextHiBoxhi0,inextHiBoxhi1
     &           ,hasHi
     &           ,iinnerCenterBoxlo0,iinnerCenterBoxlo1
     &           ,iinnerCenterBoxhi0,iinnerCenterBoxhi1
     &           ,limitC
     &           ,C3
     &           ,eps
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndWMinuscomp
      integer idWMinuslo0,idWMinuslo1
      integer idWMinushi0,idWMinushi1
      REAL*8 dWMinus(
     &           idWMinuslo0:idWMinushi0,
     &           idWMinuslo1:idWMinushi1,
     &           0:ndWMinuscomp-1)
      integer ndWPluscomp
      integer idWPluslo0,idWPluslo1
      integer idWPlushi0,idWPlushi1
      REAL*8 dWPlus(
     &           idWPluslo0:idWPlushi0,
     &           idWPluslo1:idWPlushi1,
     &           0:ndWPluscomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           0:nWcomp-1)
      integer nd2Wcomp
      integer id2Wlo0,id2Wlo1
      integer id2Whi0,id2Whi1
      REAL*8 d2W(
     &           id2Wlo0:id2Whi0,
     &           id2Wlo1:id2Whi1,
     &           0:nd2Wcomp-1)
      integer ndW2fcfcomp
      integer idW2fcflo0,idW2fcflo1
      integer idW2fcfhi0,idW2fcfhi1
      REAL*8 dW2fcf(
     &           idW2fcflo0:idW2fcfhi0,
     &           idW2fcflo1:idW2fcfhi1,
     &           0:ndW2fcfcomp-1)
      integer numSlopes
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer inextLoBoxlo0,inextLoBoxlo1
      integer inextLoBoxhi0,inextLoBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer inextHiBoxlo0,inextHiBoxlo1
      integer inextHiBoxhi0,inextHiBoxhi1
      integer hasHi
      integer iinnerCenterBoxlo0,iinnerCenterBoxlo1
      integer iinnerCenterBoxhi0,iinnerCenterBoxhi1
      REAL*8 limitC
      REAL*8 C3
      REAL*8 eps
      integer i0,i1
      integer ii0,ii1
      integer iv
      REAL*8 WC, dWM, dWP, WLL, WRR, dWL, dWR
      REAL*8 atfcf, d2Wlim, rho
      REAL*8 d2WLL, d2WL, d2WC, d2WR, d2WRR
      REAL*8 d3WLL, d3WL, d3WR, d3WRR, d3Wmax, d3Wmin
      REAL*8 sd2WL, sd2WC, sd2WR, sd2fcf
      REAL*8 dWavgM, dWavgP
      REAL*8 prodE1, prodE2, prodD3
      logical bigM, bigP
      ii0=CHF_ID(0, idir)
      ii1=CHF_ID(1, idir)
      do iv = 0,numslopes - 1
      do i1 = iinnerCenterBoxlo1,iinnerCenterBoxhi1
      do i0 = iinnerCenterBoxlo0,iinnerCenterBoxhi0
            dWM = dWMinus(i0,i1, iv)
            dWP = dWPlus (i0,i1, iv)
            bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
            bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
            WLL = W(i0 -2*ii0,i1 -2*ii1, iv)
            WC  = W(i0,i1, iv)
            WRR = W(i0 +2*ii0,i1 +2*ii1, iv)
            dWavgM = WC - WLL
            dWavgP = WRR - WC
            prodE1 = dWM * dWP
            prodE2 = dWavgM * dWavgP
            if ( (prodE1 .ge. 0) .or. (prodE2 .le. 0)) then
               d2WL  = d2W(i0 -ii0,i1 -ii1, iv)
               d2WC  = d2W(i0,i1, iv)
               d2WR  = d2W(i0 +ii0,i1 +ii1, iv)
               atfcf = dW2fcf(i0,i1, iv)
               rho = 0
               if (abs(atfcf) .ge. eps) then
                  sd2WL = sign((1.0d0), d2WL)
                  sd2WC = sign((1.0d0), d2WC)
                  sd2WR = sign((1.0d0), d2WR)
                  sd2fcf = sign((1.0d0), atfcf)
                  if ( (sd2WL .eq. sd2WC) .and.
     &                 (sd2WR .eq. sd2WC) .and.
     &                 (sd2fcf .eq. sd2WC) ) then
                     d2Wlim = sd2WC *
     &                    min(abs(atfcf),
     &                    limitC * abs(d2WC),
     &                    limitC * abs(d2WL),
     &                    limitC * abs(d2WR))
                     rho = d2Wlim / atfcf
                  endif
               endif
               if (rho .lt. ((1.0d0) - eps)) then
                  d2WLL = d2W(i0 -2*ii0,i1 -2*ii1, iv)
                  d2WRR = d2W(i0 +2*ii0,i1 +2*ii1, iv)
                  d3WLL = d2WL  - d2WLL
                  d3WL  = d2WC  - d2WL
                  d3WR  = d2WR  - d2WC
                  d3WRR = d2WRR - d2WR
                  d3Wmin = min(d3WLL, d3WL, d3WR, d3WRR)
                  d3Wmax = max(d3WLL, d3WL, d3WR, d3WRR)
                  prodD3 = C3 * max(abs(d3Wmax), abs(d3Wmin))
     &                 - abs(d3Wmax - d3Wmin)
                  if (prodD3 .le. 0) then
                     if (prodE1 .gt. 0) then
                        dWMinus(i0,i1, iv) = dWM * rho
                        dWPlus (i0,i1, iv) = dWP * rho
                     elseif ( bigM ) then
                        dWMinus(i0,i1, iv) = dWM * rho
     &                       - (2.0d0) * dWP * ((1.0d0) - rho)
                     elseif ( bigP ) then
                        dWPlus (i0,i1, iv) = dWP * rho
     &                       - (2.0d0) * dWM * ((1.0d0) - rho)
                     endif
                  endif
               endif
            else
               if ( bigM ) then
                  dWMinus(i0,i1, iv) = -(2.0d0) * dWP
               endif
               if ( bigP ) then
                  dWPlus (i0,i1, iv) = -(2.0d0) * dWM
               endif
            endif
      enddo
      enddo
         if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
               dWM = dWMinus(i0,i1, iv)
               dWP = dWPlus (i0,i1, iv)
               bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
               bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
               prodE1 = dWM * dWP
               if (prodE1 .ge. 0) then
                  d2WC  = d2W(i0,i1, iv)
                  d2WR  = d2W(i0 +ii0,i1 +ii1, iv)
                  atfcf = dW2fcf(i0,i1, iv)
                  rho = 0
                  if (abs(atfcf) .ge. eps) then
                     sd2WC = sign((1.0d0), d2WC)
                     sd2WR = sign((1.0d0), d2WR)
                     sd2fcf = sign((1.0d0), atfcf)
                     if ( (sd2WR .eq. sd2WC) .and.
     &                    (sd2fcf .eq. sd2WC) ) then
                        d2Wlim = sd2WC *
     &                       min(abs(atfcf),
     &                       limitC * abs(d2WC),
     &                       limitC * abs(d2WR))
                        rho = d2Wlim / atfcf
                     endif
                  endif
                  if (rho .lt. ((1.0d0) - eps)) then
                     d2WRR = d2W(i0 +2*ii0,i1 +2*ii1, iv)
                     d3WR  = d2WR  - d2WC
                     d3WRR = d2WRR - d2WR
                     d3Wmin = min(d3WR, d3WRR)
                     d3Wmax = max(d3WR, d3WRR)
                     prodD3 = C3 * max(abs(d3Wmax), abs(d3Wmin))
     &                    - abs(d3Wmax - d3Wmin)
                     if (prodD3 .le. 0) then
                        if (prodE1 .gt. 0) then
                           dWMinus(i0,i1, iv) = dWM * rho
                           dWPlus (i0,i1, iv) = dWP * rho
                        elseif ( bigM ) then
                           dWMinus(i0,i1, iv) = dWM * rho
     &                          - (2.0d0) * dWP * ((1.0d0) - rho)
                        elseif ( bigP ) then
                           dWPlus (i0,i1, iv) = dWP * rho
     &                          - (2.0d0) * dWM * ((1.0d0) - rho)
                        endif
                     endif
                  endif
               else
                  if ( bigM ) then
                     dWMinus(i0,i1, iv) = -(2.0d0) * dWP
                  endif
                  if ( bigP ) then
                     dWPlus (i0,i1, iv) = -(2.0d0) * dWM
                  endif
               endif
      enddo
      enddo
      do i1 = inextLoBoxlo1,inextLoBoxhi1
      do i0 = inextLoBoxlo0,inextLoBoxhi0
               dWM = dWMinus(i0,i1, iv)
               dWP = dWPlus (i0,i1, iv)
               bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
               bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
               prodE1 = dWM * dWP
               if (prodE1 .ge. 0) then
                  d2WL  = d2W(i0 -ii0,i1 -ii1, iv)
                  d2WC  = d2W(i0,i1, iv)
                  d2WR  = d2W(i0 +ii0,i1 +ii1, iv)
                  atfcf = dW2fcf(i0,i1, iv)
                  rho = 0
                  if (abs(atfcf) .ge. eps) then
                     sd2WL = sign((1.0d0), d2WL)
                     sd2WC = sign((1.0d0), d2WC)
                     sd2WR = sign((1.0d0), d2WR)
                     sd2fcf = sign((1.0d0), atfcf)
                     if ( (sd2WL .eq. sd2WC) .and.
     &                    (sd2WR .eq. sd2WC) .and.
     &                    (sd2fcf .eq. sd2WC) ) then
                        d2Wlim = sd2WC *
     &                       min(abs(atfcf),
     &                       limitC * abs(d2WC),
     &                       limitC * abs(d2WL),
     &                       limitC * abs(d2WR))
                        rho = d2Wlim / atfcf
                     endif
                  endif
                  if (rho .lt. ((1.0d0) - eps)) then
                     d2WRR = d2W(i0 +2*ii0,i1 +2*ii1, iv)
                     d3WL  = d2WC  - d2WL
                     d3WR  = d2WR  - d2WC
                     d3WRR = d2WRR - d2WR
                     d3Wmin = min(d3WL, d3WR, d3WRR)
                     d3Wmax = max(d3WL, d3WR, d3WRR)
                     prodD3 = C3 * max(abs(d3Wmax), abs(d3Wmin))
     &                    - abs(d3Wmax - d3Wmin)
                     if (prodD3 .le. 0) then
                        if (prodE1 .gt. 0) then
                           dWMinus(i0,i1, iv) = dWM * rho
                           dWPlus (i0,i1, iv) = dWP * rho
                        elseif ( bigM ) then
                           dWMinus(i0,i1, iv) = dWM * rho
     &                          - (2.0d0) * dWP * ((1.0d0) - rho)
                        elseif ( bigP ) then
                           dWPlus (i0,i1, iv) = dWP * rho
     &                          - (2.0d0) * dWM * ((1.0d0) - rho)
                        endif
                     endif
                  endif
               else
                  if ( bigM ) then
                     dWMinus(i0,i1, iv) = -(2.0d0) * dWP
                  endif
                  if ( bigP ) then
                     dWPlus (i0,i1, iv) = -(2.0d0) * dWM
                  endif
               endif
      enddo
      enddo
         endif
         if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
               dWM = dWMinus(i0,i1, iv)
               dWP = dWPlus (i0,i1, iv)
               bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
               bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
               prodE1 = dWM * dWP
               if (prodE1 .ge. 0) then
                  d2WL  = d2W(i0 -ii0,i1 -ii1, iv)
                  d2WC  = d2W(i0,i1, iv)
                  atfcf = dW2fcf(i0,i1, iv)
                  rho = 0
                  if (abs(atfcf) .ge. eps) then
                     sd2WL = sign((1.0d0), d2WL)
                     sd2WC = sign((1.0d0), d2WC)
                     sd2fcf = sign((1.0d0), atfcf)
                     if ( (sd2WL .eq. sd2WC) .and.
     &                    (sd2fcf .eq. sd2WC) ) then
                        d2Wlim = sd2WC *
     &                       min(abs(atfcf),
     &                       limitC * abs(d2WC),
     &                       limitC * abs(d2WL))
                        rho = d2Wlim / atfcf
                     endif
                  endif
                  if (rho .lt. ((1.0d0) - eps)) then
                     d2WLL = d2W(i0 -2*ii0,i1 -2*ii1, iv)
                     d3WLL = d2WL  - d2WLL
                     d3WL  = d2WC  - d2WL
                     d3Wmin = min(d3WLL, d3WL)
                     d3Wmax = max(d3WLL, d3WL)
                     prodD3 = C3 * max(abs(d3Wmax), abs(d3Wmin))
     &                    - abs(d3Wmax - d3Wmin)
                     if (prodD3 .le. 0) then
                        if (prodE1 .gt. 0) then
                           dWMinus(i0,i1, iv) = dWM * rho
                           dWPlus (i0,i1, iv) = dWP * rho
                        elseif ( bigM ) then
                           dWMinus(i0,i1, iv) = dWM * rho
     &                          - (2.0d0) * dWP * ((1.0d0) - rho)
                        elseif ( bigP ) then
                           dWPlus (i0,i1, iv) = dWP * rho
     &                          - (2.0d0) * dWM * ((1.0d0) - rho)
                        endif
                     endif
                  endif
               else
                  if ( bigM ) then
                     dWMinus(i0,i1, iv) = -(2.0d0) * dWP
                  endif
                  if ( bigP ) then
                     dWPlus (i0,i1, iv) = -(2.0d0) * dWM
                  endif
               endif
      enddo
      enddo
      do i1 = inextHiBoxlo1,inextHiBoxhi1
      do i0 = inextHiBoxlo0,inextHiBoxhi0
               dWM = dWMinus(i0,i1, iv)
               dWP = dWPlus (i0,i1, iv)
               bigM = ( abs(dWM) .gt. (2.0d0) * abs(dWP) )
               bigP = ( abs(dWP) .gt. (2.0d0) * abs(dWM) )
               prodE1 = dWM * dWP
               if (prodE1 .ge. 0) then
                  d2WL  = d2W(i0 -ii0,i1 -ii1, iv)
                  d2WC  = d2W(i0,i1, iv)
                  d2WR  = d2W(i0 +ii0,i1 +ii1, iv)
                  atfcf = dW2fcf(i0,i1, iv)
                  rho = 0
                  if (abs(atfcf) .ge. eps) then
                     sd2WL = sign((1.0d0), d2WL)
                     sd2WC = sign((1.0d0), d2WC)
                     sd2WR = sign((1.0d0), d2WR)
                     sd2fcf = sign((1.0d0), atfcf)
                     if ( (sd2WL .eq. sd2WC) .and.
     &                    (sd2WR .eq. sd2WC) .and.
     &                    (sd2fcf .eq. sd2WC) ) then
                        d2Wlim = sd2WC *
     &                       min(abs(atfcf),
     &                       limitC * abs(d2WC),
     &                       limitC * abs(d2WL),
     &                       limitC * abs(d2WR))
                        rho = d2Wlim / atfcf
                     endif
                  endif
                  if (rho .lt. ((1.0d0) - eps)) then
                     d2WLL = d2W(i0 -2*ii0,i1 -2*ii1, iv)
                     d3WLL = d2WL  - d2WLL
                     d3WL  = d2WC  - d2WL
                     d3WR  = d2WR  - d2WC
                     d3Wmin = min(d3WLL, d3WL, d3WR)
                     d3Wmax = max(d3WLL, d3WL, d3WR)
                     prodD3 = C3 * max(abs(d3Wmax), abs(d3Wmin))
     &                    - abs(d3Wmax - d3Wmin)
                     if (prodD3 .le. 0) then
                        if (prodE1 .gt. 0) then
                           dWMinus(i0,i1, iv) = dWM * rho
                           dWPlus (i0,i1, iv) = dWP * rho
                        elseif ( bigM ) then
                           dWMinus(i0,i1, iv) = dWM * rho
     &                          - (2.0d0) * dWP * ((1.0d0) - rho)
                        elseif ( bigP ) then
                           dWPlus (i0,i1, iv) = dWP * rho
     &                          - (2.0d0) * dWM * ((1.0d0) - rho)
                        endif
                     endif
                  endif
               else
                  if ( bigM ) then
                     dWMinus(i0,i1, iv) = -(2.0d0) * dWP
                  endif
                  if ( bigP ) then
                     dWPlus (i0,i1, iv) = -(2.0d0) * dWM
                  endif
               endif
      enddo
      enddo
         endif
      enddo
      return
      end
