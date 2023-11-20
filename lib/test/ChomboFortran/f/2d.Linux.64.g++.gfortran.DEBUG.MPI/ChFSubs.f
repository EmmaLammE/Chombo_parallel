      subroutine TESTVIRC_1(
     & vectorI
     & ,ivectorIhi0
     & ,VectorCI
     & ,iVectorCIhi0
     & ,vectorR
     & ,ivectorRhi0
     & ,VectorCR
     & ,iVectorCRhi0
     & ,vCmplx
     & ,ivCmplxhi0
     & ,Eps
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ivectorIhi0
      integer vectorI(
     & 0:ivectorIhi0)
      integer iVectorCIhi0
      integer VectorCI(
     & 0:iVectorCIhi0)
      integer ivectorRhi0
      REAL*8 vectorR(
     & 0:ivectorRhi0)
      integer iVectorCRhi0
      REAL*8 VectorCR(
     & 0:iVectorCRhi0)
      integer ivCmplxhi0
      COMPLEX*16 vCmplx(
     & 0:ivCmplxhi0)
      REAL*8 Eps
      integer status
      REAL*8 R10 ,R11 ,R12
      parameter( R10 = 1.1d0 ,R11 = 22.22d0 ,R12 = 333.333d0 )
      REAL*8 tmp,tmp2,diff
      COMPLEX*16 ctmp
      status = 0
      if ( ivectorIhi0 .NE. 0 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got integer arg 1 ubound = '
     & ,ivectorIhi0,' instead of 0'
      endif
      if ( vectorI(0) .NE. 12345678 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got integer arg 1 = '
     & ,vectorI(0) ,' instead of 12345678'
      endif
      vectorI(0) = 87654321
      if ( ivectorCIhi0 .NE. 1 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const integer arg 2 ubound 
     &= '
     & ,ivectorCIhi0,' instead of 1'
      endif
      if ( vectorCI(0) .NE. 23456789 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const integer arg 2 [0] = '
     & ,vectorCI(0) ,' instead of 23456789'
      endif
      if ( vectorCI(1) .NE. 98765432 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const integer arg 2 [1] = '
     & ,vectorCI(1) ,' instead of 98765432'
      endif
      if ( ivectorRhi0 .NE. 2 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got real arg 3 ubound = '
     & ,ivectorRhi0,' instead of 2'
      endif
      if ( vectorR(0) .NE. R10 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got real arg 3 [0] = '
     & ,vectorR(0) ,' instead of ',R10
      endif
      vectorR(0) = 1.0d0 / vectorR(0)
      if ( vectorR(1) .NE. R11 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got real arg 3 [1] = '
     & ,vectorR(1) ,' instead of ',R11
      endif
      vectorR(1) = 1.0d0 / vectorR(1)
      if ( vectorR(2) .NE. R12 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got real arg 3 [2] = '
     & ,vectorR(2) ,' instead of ',R12
      endif
      vectorR(2) = 1.0d0 / vectorR(2)
      if ( ivectorCRhi0 .NE. 3 )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const real arg 4 ubound = '
     & ,ivectorCRhi0,' instead of 3'
      endif
      tmp = 2**20
      tmp = 1.0d0 + ( 1.0d0 / tmp )
      if ( vectorCR(0) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const real arg 4 [0] = '
     & ,vectorCR(0) ,' instead of ', tmp
      endif
      tmp = 1.0d0 / tmp
      diff = vectorCR(1) - tmp
      if ( ABS( diff ) .GT. Eps )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const real arg 4 [1] = '
     & ,vectorCR(1) ,' instead of ', tmp
     & , ', diff = ' ,diff ,".GT." ,Eps
      endif
      tmp = 2**20
      tmp = tmp * tmp
      tmp = 1.0d0 + ( 1.0d0 / tmp )
      if ( vectorCR(2) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const real arg 4 [2] = '
     & ,vectorCR(2) ,' instead of ', tmp
      endif
      tmp = 1.0d0 / tmp
      if ( vectorCR(3) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testVIR_1 got const real arg 4 [3] = '
     & ,vectorCR(3) ,' instead of ', tmp
      endif
      if ( ivCmplxhi0 .NE. 1 )then
        status = status + 1
        write(6,*) '   error: testVIR_1 got const integer arg 5 ubound =
     & '
     & ,ivCmplxhi0,' instead of 1'
      endif
      ctmp = vCmplx(0)
      tmp = ctmp
      tmp2= IMAG(ctmp)
      if (vectorCR(0) .NE. tmp)then
          status = status+1
          write(6,*) '    error: testVIR_1 wrong real part of complex  '
     & ,tmp ,' instead of ', vectorCR(0)
      endif
      vCmplx(0) = CMPLX(tmp2,tmp/2)
      ctmp = vCmplx(1)
      tmp = ctmp
      tmp2= IMAG(ctmp)
      if (vectorCR(3) .NE. tmp2)then
          status = status+1
          write(6,*) '    error: testVIR_1 wrong img part of complex  '
     & ,tmp ,' instead of ', vectorCR(3)
      endif
      vCmplx(1) = CMPLX(tmp2,tmp/2)
      return
      end
      subroutine TESTIR1D_1(
     & arrayI
     & ,iarrayIhi0
     & ,arrayCI
     & ,iarrayCIhi0
     & ,arrayR
     & ,iarrayRhi0
     & ,arrayCR
     & ,iarrayCRhi0
     & ,Eps
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrayIhi0
      integer arrayI(
     & 0:iarrayIhi0)
      integer iarrayCIhi0
      integer arrayCI(
     & 0:iarrayCIhi0)
      integer iarrayRhi0
      REAL*8 arrayR(
     & 0:iarrayRhi0)
      integer iarrayCRhi0
      REAL*8 arrayCR(
     & 0:iarrayCRhi0)
      REAL*8 Eps
      integer status
      REAL*8 R10 ,R11 ,R12
      parameter( R10 = 1.1d0 ,R11 = 22.22d0 ,R12 = 333.333d0 )
      REAL*8 tmp ,diff
      status = 0
      if ( iarrayIhi0 .NE. 0 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got integer arg 1 ubound = '
     & ,iarrayIhi0,' instead of 0'
      endif
      if ( arrayI(0) .NE. 12345678 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got integer arg 1 = '
     & ,arrayI(0) ,' instead of 12345678'
      endif
      arrayI(0) = 87654321
      if ( iarrayCIhi0 .NE. 1 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const integer arg 2 ubound
     & = '
     & ,iarrayCIhi0,' instead of 1'
      endif
      if ( arrayCI(0) .NE. 23456789 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const integer arg 2 [0] = 
     &'
     & ,arrayCI(0) ,' instead of 23456789'
      endif
      if ( arrayCI(1) .NE. 98765432 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const integer arg 2 [1] = 
     &'
     & ,arrayCI(1) ,' instead of 98765432'
      endif
      if ( iarrayRhi0 .NE. 2 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got real arg 3 ubound = '
     & ,iarrayRhi0,' instead of 2'
      endif
      if ( arrayR(0) .NE. R10 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got real arg 3 [0] = '
     & ,arrayR(0) ,' instead of ',R10
      endif
      arrayR(0) = 1.0d0 / arrayR(0)
      if ( arrayR(1) .NE. R11 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got real arg 3 [1] = '
     & ,arrayR(1) ,' instead of ',R11
      endif
      arrayR(1) = 1.0d0 / arrayR(1)
      if ( arrayR(2) .NE. R12 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got real arg 3 [2] = '
     & ,arrayR(2) ,' instead of ',R12
      endif
      arrayR(2) = 1.0d0 / arrayR(2)
      if ( iarrayCRhi0 .NE. 3 )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const real arg 4 ubound = 
     &'
     & ,iarrayCRhi0,' instead of 3'
      endif
      tmp = 2**20
      tmp = 1.0d0 + ( 1.0d0 / tmp )
      if ( arrayCR(0) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const real arg 4 [0] = '
     & ,arrayCR(0) ,' instead of ', tmp
      endif
      tmp = 1.0d0 / tmp
      diff = arrayCR(1) - tmp
      if ( ABS( diff ) .GT. Eps )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const real arg 4 [1] = '
     & ,arrayCR(1) ,' instead of ', tmp
     & , ', diff = ' ,diff ,".GT." ,Eps
      endif
      tmp = 2**20
      tmp = tmp * tmp
      tmp = 1.0d0 + ( 1.0d0 / tmp )
      if ( arrayCR(2) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const real arg 4 [2] = '
     & ,arrayCR(2) ,' instead of ', tmp
      endif
      tmp = 1.0d0 / tmp
      if ( arrayCR(3) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testIR1D_1 got const real arg 4 [3] = '
     & ,arrayCR(3) ,' instead of ', tmp
      endif
      return
      end
      subroutine TESTIR1D_2(
     & arrayI
     & ,iarrayIhi0
     & ,arrayCI
     & ,iarrayCIhi0
     & ,arrayR
     & ,iarrayRhi0
     & ,arrayCR
     & ,iarrayCRhi0
     & ,Eps
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrayIhi0
      integer arrayI(
     & 0:iarrayIhi0)
      integer iarrayCIhi0
      integer arrayCI(
     & 0:iarrayCIhi0)
      integer iarrayRhi0
      REAL*8 arrayR(
     & 0:iarrayRhi0)
      integer iarrayCRhi0
      REAL*8 arrayCR(
     & 0:iarrayCRhi0)
      REAL*8 Eps
      integer status
      REAL*8 R10 ,R11 ,R12
      parameter( R10 = 1.1d0 ,R11 = 22.22d0 ,R12 = 333.333d0 )
      REAL*8 tmp ,diff
      status = 0
      if ( iarrayIhi0 .NE. 1 )then
        status = status + 1
        write(6,*) '    error: testCALL got integer arg 1 ubound = '
     & ,iarrayIhi0,' instead of 1'
      endif
      if ( arrayI(0) .NE. 12345678 )then
        status = status + 1
        write(6,*) '    error: testCALL got integer arg 1 = '
     & ,arrayI(0) ,' instead of 12345678'
      endif
      arrayI(0) = 87654321
      if ( iarrayCIhi0 .NE. 1 )then
        status = status + 1
        write(6,*) '    error: testCALL got const integer arg 2 ubound =
     & '
     & ,iarrayCIhi0,' instead of 1'
      endif
      if ( arrayCI(0) .NE. 23456789 )then
        status = status + 1
        write(6,*) '    error: testCALL got const integer arg 2 [0] = '
     & ,arrayCI(0) ,' instead of 23456789'
      endif
      if ( arrayCI(1) .NE. 98765432 )then
        status = status + 1
        write(6,*) '    error: testCALL got const integer arg 2 [1] = '
     & ,arrayCI(1) ,' instead of 98765432'
      endif
      if ( iarrayRhi0 .NE. 2 )then
        status = status + 1
        write(6,*) '    error: testCALL got real arg 3 ubound = '
     & ,iarrayRhi0,' instead of 2'
      endif
      if ( arrayR(0) .NE. R10 )then
        status = status + 1
        write(6,*) '    error: testCALL got real arg 3 [0] = '
     & ,arrayR(0) ,' instead of ',R10
      endif
      arrayR(0) = 1.0d0 / arrayR(0)
      if ( arrayR(1) .NE. R11 )then
        status = status + 1
        write(6,*) '    error: testCALL got real arg 3 [1] = '
     & ,arrayR(1) ,' instead of ',R11
      endif
      arrayR(1) = 1.0d0 / arrayR(1)
      if ( arrayR(2) .NE. R12 )then
        status = status + 1
        write(6,*) '    error: testCALL got real arg 3 [2] = '
     & ,arrayR(2) ,' instead of ',R12
      endif
      arrayR(2) = 1.0d0 / arrayR(2)
      if ( iarrayCRhi0 .NE. 3 )then
        status = status + 1
        write(6,*) '    error: testCALL got const real arg 4 ubound = '
     & ,iarrayCRhi0,' instead of 3'
      endif
      tmp = 2**20
      tmp = 1.0d0 + ( 1.0d0 / tmp )
      if ( arrayCR(0) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testCALL got const real arg 4 [0] = '
     & ,arrayCR(0) ,' instead of ', tmp
      endif
      tmp = 1.0d0 / tmp
      diff = arrayCR(1) - tmp
      if ( ABS( diff ) .GT. Eps )then
        status = status + 1
        write(6,*) '    error: testCALL got const real arg 4 [1] = '
     & ,arrayCR(1) ,' instead of ', tmp
     & , ', diff = ' ,diff ,".GT." ,Eps
      endif
      tmp = 2**20
      tmp = tmp * tmp
      tmp = 1.0d0 + ( 1.0d0 / tmp )
      if ( arrayCR(2) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testCALL got const real arg 4 [2] = '
     & ,arrayCR(2) ,' instead of ', tmp
      endif
      tmp = 1.0d0 / tmp
      if ( arrayCR(3) .NE. tmp )then
        status = status + 1
        write(6,*) '    error: testCALL got const real arg 4 [3] = '
     & ,arrayCR(3) ,' instead of ', tmp
      endif
      return
      end
      subroutine TESTDSEL_1(
     & SpaceDim
     & ,Verbose
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer SpaceDim
      integer Verbose
      integer status
      integer I1 ,I2 ,I3, I4, I5, I6
      parameter( I1 = 1 ,I2 = 2 ,I3 = 3, I4 = 4, I5 = 5, I6 = 6 )
      integer tmp ,tmp3(1:6)
      status = 0
      tmp = I2
      if ( tmp .NE. SpaceDim )then
                if ( Verbose .GT. 0 )then
          write(6,*) '    error: testDSEL_1 simple value selection for C
     &H_SPACEDIM = '
     & ,2, ' got ',tmp ,' instead of ' ,SpaceDim
        endif
      endif
      tmp3(1) = 0
      tmp3(2) = 0
      tmp3(3) = 0
                   tmp3(2) = 2
      if ( tmp3(2) .EQ. 0 )then
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testDSEL_1 statement selection for CH_S
     &PACEDIM = '
     & ,2, 'got (:)=',tmp3 ,'instead of (', 2 ,') = '
     & ,2 +1
        endif
      endif
      return
      end
      subroutine TESTMDO_1(
     & a
     & ,ialo0,ialo1
     & ,iahi0,iahi1
     & ,ibalo0,ibalo1
     & ,ibahi0,ibahi1
     & ,Vol
     & ,Verbose
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialo0,ialo1
      integer iahi0,iahi1
      integer a(
     & ialo0:iahi0,
     & ialo1:iahi1)
      integer ibalo0,ibalo1
      integer ibahi0,ibahi1
      integer Vol
      integer Verbose
      integer status
      integer i1,i2,stride,sumpts
      REAL*8 suma
      stride = -2
      status = 0
      sumpts = 0
      suma = 0.0
      do i2 = ibalo1,ibahi1
      do i1 = ibalo0,ibahi0
        suma = suma + a(i1,i2)
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDO_1: simple 3arg MULTIDO got '
     & ,'volume ',sumpts,' instead of ',Vol
        endif
      endif
      sumpts = 0
      do i2 = ibalo1,ibahi1
      do i1 = ibalo0,ibahi0
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDO_1: 3arg MULTIDO with explicit '
     & ,'stride 1 got volume ',sumpts,' instead of ',Vol
        endif
      endif
      sumpts = 0
      do i2 = ibalo1,ibahi1,2
      do i1 = ibalo0,ibahi0,2
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol/(2**2) )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDO_1: 3arg MULTIDO with stride 2 '
     & ,'got volume ',sumpts,' instead of '
     & ,Vol/(2**2)
        endif
      endif
      sumpts = 0
      do i2 = ibahi1,ibalo1,-1
      do i1 = ibahi0,ibalo0,-1
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDO_1: 3arg MULTIDO with stride -1'
     & ,' got volume ',sumpts,' instead of ',Vol
        endif
      endif
      sumpts = 0
      do i2 = ibalo1,ibahi1
      do i1 = ibalo0,ibahi0
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDO_1: 3arg MULTIDO with variable '
     & ,'stride got volume ',sumpts,' instead of ',Vol
        endif
      endif
      return
      end
      subroutine TESTAMDO_1(
     & a
     & ,ialo0,ialo1
     & ,iahi0,iahi1
     & ,ibalo0,ibalo1
     & ,ibahi0,ibahi1
     & ,Vol
     & ,ioffsetboxlo0,ioffsetboxlo1
     & ,ioffsetboxhi0,ioffsetboxhi1
     & ,Verbose
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialo0,ialo1
      integer iahi0,iahi1
      integer a(
     & ialo0:iahi0,
     & ialo1:iahi1)
      integer ibalo0,ibalo1
      integer ibahi0,ibahi1
      integer Vol
      integer ioffsetboxlo0,ioffsetboxlo1
      integer ioffsetboxhi0,ioffsetboxhi1
      integer Verbose
      integer status
      integer i0,i1,stride,sumpts
      integer ii0,ii1
      REAL*8 suma
      stride = -2
      status = 0
      ii0=1
                ii1=1
      sumpts = 0
      suma = 0.0
      do i1 = ibalo1,ibahi1
      do i0 = ibalo0,ibahi0
        suma = suma + a(i0,i1)
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDOA_1: simple 3arg AUTOMULTIDO got
     & '
     & ,'volume ',sumpts,' instead of ',Vol
        endif
      endif
      sumpts = 0
      suma = 0.0
      do i1 = ioffsetboxlo1,ioffsetboxhi1
      do i0 = ioffsetboxlo0,ioffsetboxhi0
        suma = suma + a(i0-ii0,i1-ii1)
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDOA_1: offset 3arg AUTOMULTIDO '
     & ,'with OFFSETIX got volume ',sumpts,' instead of ',Vol
        endif
      endif
      sumpts = 0
      do i1 = ibalo1,ibahi1
      do i0 = ibalo0,ibahi0
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDOA_1: 3arg AUTOMULTIDO with expli
     &cit '
     & ,'stride 1 got volume ',sumpts,' instead of ',Vol
        endif
      endif
      sumpts = 0
      do i1 = ibalo1,ibahi1,2
      do i0 = ibalo0,ibahi0,2
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol/(2**2) )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDOA_1: 3arg AUTOMULTIDO with strid
     &e 2 '
     & ,'got volume ',sumpts,' instead of '
     & ,Vol/(2**2)
        endif
      endif
      sumpts = 0
      do i1 = ibahi1,ibalo1,-1
      do i0 = ibahi0,ibalo0,-1
        sumpts = sumpts + 1
      enddo
      enddo
      if ( sumpts .NE. Vol )then
        status = status + 1
        if ( Verbose .GT. 0 )then
          write(6,*) '    error: testMDOA_1: 3arg AUTOMULTIDO with strid
     &e -1'
     & ,' got volume ',sumpts,' instead of ',Vol
        endif
      endif
      return
      end
      subroutine TEST_CHFID(
     & status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer status
      integer dir1, dir2, testval
      do dir1=0, 2 -1
         do dir2=0, 2 -1
            testval = CHF_ID(dir1, dir2)
            if (dir1.eq.dir2) then
               if (testval.eq.1) then
                  status = status + 0
               else
                  status = status + 1
               endif
            else
               if (testval.eq.0) then
                  status = status + 0
               else
                  status = status + 1
               endif
            endif
         enddo
      enddo
      return
      end
      subroutine TESTCALL(
     & arrayI
     & ,iarrayIhi0
     & ,arrayCI
     & ,iarrayCIhi0
     & ,arrayR
     & ,iarrayRhi0
     & ,arrayCR
     & ,iarrayCRhi0
     & ,Eps
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrayIhi0
      integer arrayI(
     & 0:iarrayIhi0)
      integer iarrayCIhi0
      integer arrayCI(
     & 0:iarrayCIhi0)
      integer iarrayRhi0
      REAL*8 arrayR(
     & 0:iarrayRhi0)
      integer iarrayCRhi0
      REAL*8 arrayCR(
     & 0:iarrayCRhi0)
      REAL*8 Eps
      integer status
      integer int1,status2
      integer int2(1,1)
      call TESTIR1D_2(
     & arrayI
     & ,iarrayIhi0
     & ,arrayCI
     & ,iarrayCIhi0
     & ,arrayR
     & ,iarrayRhi0
     & ,arrayCR
     & ,iarrayCRhi0
     & ,Eps
     & ,status
     & )
      int1 = 23
      call TESTCALL_2( int1 ,status2 )
      status = status + status2
      if ( int1 .NE. 32 )then
        write(6,*) 'error: TESTCALL: plain Fortran call to TESTCALL_2 '
     & ,'returned with arg 1 = ',int1,' instead of 32.'
        status = status + 1
      endif
      int2(1,1) = 23
      call TESTCALL_2( int2(1,1) ,status2 )
      status = status + status2
      if ( int2(1,1) .NE. 32 )then
        write(6,*) 'error: TESTCALL: plain Fortran call to '
     & ,'TESTCALL_2() with CHF IX returned with arg 1 = '
     & ,int2(1,1),' instead of 32.'
        status = status + 1
      endif
      arrayI(iarrayIhi0) = 23
      call TESTCALL_2( arrayI(iarrayIhi0) ,status2 )
      if ( int2(1,1) .NE. 32 )then
        write(6,*) 'error: TESTCALL: plain Fortran call to '
     & ,'TESTCALL_2() with CHF UBOUND returned with arg 1 = '
     & ,int2(1,1),' instead of 32.'
        status = status + 1
      endif
      if ( Eps .LT. 0 ) call MAYDAY_ERROR(
     & )
      return
      end
      subroutine TEST_PARSE1(
     & status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer status
      status = 0
      return
      end
      integer function LENC( str )
      character str*(*)
      integer l
      do l = LEN(str),1,-1
        if ( str(l:l) .NE. ' ' )goto 10
      enddo
      l = 0
   10 LENC = l
      return
      end
      subroutine TEST_SHIFT_FIA(
     & crFIA
     & ,icrFIAlo0,icrFIAlo1
     & ,icrFIAhi0,icrFIAhi1
     & ,ncrFIAcomp
     & ,ifnBoxlo0,ifnBoxlo1
     & ,ifnBoxhi0,ifnBoxhi1
     & ,fnFIA
     & ,ifnFIAlo0,ifnFIAlo1
     & ,ifnFIAhi0,ifnFIAhi1
     & ,nfnFIAcomp
     & ,refRatio
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncrFIAcomp
      integer icrFIAlo0,icrFIAlo1
      integer icrFIAhi0,icrFIAhi1
      integer crFIA(
     & icrFIAlo0:icrFIAhi0,
     & icrFIAlo1:icrFIAhi1,
     & 0:ncrFIAcomp-1)
      integer ifnBoxlo0,ifnBoxlo1
      integer ifnBoxhi0,ifnBoxhi1
      integer nfnFIAcomp
      integer ifnFIAlo0,ifnFIAlo1
      integer ifnFIAhi0,ifnFIAhi1
      integer fnFIA(
     & ifnFIAlo0:ifnFIAhi0,
     & ifnFIAlo1:ifnFIAhi1,
     & 0:nfnFIAcomp-1)
      integer refRatio
      integer iFn0,iFn1
      integer iCr0,iCr1
      do iFn1 = ifnBoxlo1,ifnBoxhi1
      do iFn0 = ifnBoxlo0,ifnBoxhi0
            iCr0 = iFn0/refRatio
            iCr1 = iFn1/refRatio
         crFIA(iCr0,iCr1, 0) = crFIA(iCr0,iCr1, 0) +
     & fnFIA(iFn0,iFn1, 0)
      enddo
      enddo
      return
      end
      subroutine TEST_SHIFT_FIA1(
     & crFIA
     & ,icrFIAlo0,icrFIAlo1
     & ,icrFIAhi0,icrFIAhi1
     & ,ifnBoxlo0,ifnBoxlo1
     & ,ifnBoxhi0,ifnBoxhi1
     & ,fnFIA
     & ,ifnFIAlo0,ifnFIAlo1
     & ,ifnFIAhi0,ifnFIAhi1
     & ,refRatio
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icrFIAlo0,icrFIAlo1
      integer icrFIAhi0,icrFIAhi1
      integer crFIA(
     & icrFIAlo0:icrFIAhi0,
     & icrFIAlo1:icrFIAhi1)
      integer ifnBoxlo0,ifnBoxlo1
      integer ifnBoxhi0,ifnBoxhi1
      integer ifnFIAlo0,ifnFIAlo1
      integer ifnFIAhi0,ifnFIAhi1
      integer fnFIA(
     & ifnFIAlo0:ifnFIAhi0,
     & ifnFIAlo1:ifnFIAhi1)
      integer refRatio
      integer iFn0,iFn1
      integer iCr0,iCr1
      do iFn1 = ifnBoxlo1,ifnBoxhi1
      do iFn0 = ifnBoxlo0,ifnBoxhi0
            iCr0 = iFn0/refRatio
            iCr1 = iFn1/refRatio
         crFIA(iCr0,iCr1) = crFIA(iCr0,iCr1) +
     & fnFIA(iFn0,iFn1)
      enddo
      enddo
      return
      end
      subroutine TEST_SHIFT_FRA(
     & crFRA
     & ,icrFRAlo0,icrFRAlo1
     & ,icrFRAhi0,icrFRAhi1
     & ,ncrFRAcomp
     & ,ifnBoxlo0,ifnBoxlo1
     & ,ifnBoxhi0,ifnBoxhi1
     & ,fnFRA
     & ,ifnFRAlo0,ifnFRAlo1
     & ,ifnFRAhi0,ifnFRAhi1
     & ,nfnFRAcomp
     & ,refRatio
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncrFRAcomp
      integer icrFRAlo0,icrFRAlo1
      integer icrFRAhi0,icrFRAhi1
      REAL*8 crFRA(
     & icrFRAlo0:icrFRAhi0,
     & icrFRAlo1:icrFRAhi1,
     & 0:ncrFRAcomp-1)
      integer ifnBoxlo0,ifnBoxlo1
      integer ifnBoxhi0,ifnBoxhi1
      integer nfnFRAcomp
      integer ifnFRAlo0,ifnFRAlo1
      integer ifnFRAhi0,ifnFRAhi1
      REAL*8 fnFRA(
     & ifnFRAlo0:ifnFRAhi0,
     & ifnFRAlo1:ifnFRAhi1,
     & 0:nfnFRAcomp-1)
      integer refRatio
      integer iFn0,iFn1
      integer iCr0,iCr1
      do iFn1 = ifnBoxlo1,ifnBoxhi1
      do iFn0 = ifnBoxlo0,ifnBoxhi0
            iCr0 = iFn0/refRatio
            iCr1 = iFn1/refRatio
         crFRA(iCr0,iCr1, 0) = crFRA(iCr0,iCr1, 0) +
     & fnFRA(iFn0,iFn1, 0)
      enddo
      enddo
      return
      end
      subroutine TEST_SHIFT_FRA1(
     & crFRA
     & ,icrFRAlo0,icrFRAlo1
     & ,icrFRAhi0,icrFRAhi1
     & ,ifnBoxlo0,ifnBoxlo1
     & ,ifnBoxhi0,ifnBoxhi1
     & ,fnFRA
     & ,ifnFRAlo0,ifnFRAlo1
     & ,ifnFRAhi0,ifnFRAhi1
     & ,refRatio
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icrFRAlo0,icrFRAlo1
      integer icrFRAhi0,icrFRAhi1
      REAL*8 crFRA(
     & icrFRAlo0:icrFRAhi0,
     & icrFRAlo1:icrFRAhi1)
      integer ifnBoxlo0,ifnBoxlo1
      integer ifnBoxhi0,ifnBoxhi1
      integer ifnFRAlo0,ifnFRAlo1
      integer ifnFRAhi0,ifnFRAhi1
      REAL*8 fnFRA(
     & ifnFRAlo0:ifnFRAhi0,
     & ifnFRAlo1:ifnFRAhi1)
      integer refRatio
      integer iFn0,iFn1
      integer iCr0,iCr1
      do iFn1 = ifnBoxlo1,ifnBoxhi1
      do iFn0 = ifnBoxlo0,ifnBoxhi0
            iCr0 = iFn0/refRatio
            iCr1 = iFn1/refRatio
         crFRA(iCr0,iCr1) = crFRA(iCr0,iCr1) +
     & fnFRA(iFn0,iFn1)
      enddo
      enddo
      return
      end
