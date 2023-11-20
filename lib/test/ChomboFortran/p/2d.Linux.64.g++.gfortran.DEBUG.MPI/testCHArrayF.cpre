      subroutine RANK1CHARRAY(
     & arr
     & ,iarrlo0
     & ,iarrhi0
     & ,arrc
     & ,iarrclo0
     & ,iarrchi0
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrhi0
      integer arr(
     & iarrlo0:iarrhi0)
      integer iarrclo0
      integer iarrchi0
      integer arrc(
     & iarrclo0:iarrchi0)
      integer i0
      integer inc
      inc = 1
      do i0 = iarrclo0, iarrchi0
         arr(i0) = arrc(i0) + inc
         inc = inc + 1
      enddo
      return
      end
      subroutine RANK2CHARRAY(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrhi0
     & ,iarrhi1
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrchi0
     & ,iarrchi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrhi0
      integer iarrhi1
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1)
      integer iarrclo0
      integer iarrclo1
      integer iarrchi0
      integer iarrchi1
      integer arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1)
      integer i0, i1
      integer inc
      inc = 1
      do i1 = iarrclo1, iarrchi1
         do i0 = iarrclo0, iarrchi0
            arr(i0, i1) = arrc(i0, i1) + inc
            inc = inc + 1
         enddo
      enddo
      return
      end
      subroutine RANK3CHARRAY(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrlo2
     & ,iarrhi0
     & ,iarrhi1
     & ,iarrhi2
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrclo2
     & ,iarrchi0
     & ,iarrchi1
     & ,iarrchi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrlo2
      integer iarrhi0
      integer iarrhi1
      integer iarrhi2
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2)
      integer iarrclo0
      integer iarrclo1
      integer iarrclo2
      integer iarrchi0
      integer iarrchi1
      integer iarrchi2
      integer arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1,
     & iarrclo2:iarrchi2)
      integer i0, i1, i2
      integer inc
      inc = 1
      do i2 = iarrclo2, iarrchi2
         do i1 = iarrclo1, iarrchi1
            do i0 = iarrclo0, iarrchi0
               arr(i0, i1, i2) = arrc(i0, i1, i2) + inc
               inc = inc + 1
            enddo
         enddo
      enddo
      return
      end
      subroutine RANK4CHARRAY(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrlo2
     & ,iarrlo3
     & ,iarrhi0
     & ,iarrhi1
     & ,iarrhi2
     & ,iarrhi3
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrclo2
     & ,iarrclo3
     & ,iarrchi0
     & ,iarrchi1
     & ,iarrchi2
     & ,iarrchi3
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrlo2
      integer iarrlo3
      integer iarrhi0
      integer iarrhi1
      integer iarrhi2
      integer iarrhi3
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2,
     & iarrlo3:iarrhi3)
      integer iarrclo0
      integer iarrclo1
      integer iarrclo2
      integer iarrclo3
      integer iarrchi0
      integer iarrchi1
      integer iarrchi2
      integer iarrchi3
      integer arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1,
     & iarrclo2:iarrchi2,
     & iarrclo3:iarrchi3)
      integer i0, i1, i2, i3
      integer inc
      inc = 1
      do i3 = iarrclo3, iarrchi3
         do i2 = iarrclo2, iarrchi2
            do i1 = iarrclo1, iarrchi1
               do i0 = iarrclo0, iarrchi0
                  arr(i0, i1, i2, i3) = arrc(i0, i1, i2, i3) + inc
                  inc = inc + 1
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine RANK5CHARRAY(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrlo2
     & ,iarrlo3
     & ,iarrlo4
     & ,iarrhi0
     & ,iarrhi1
     & ,iarrhi2
     & ,iarrhi3
     & ,iarrhi4
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrclo2
     & ,iarrclo3
     & ,iarrclo4
     & ,iarrchi0
     & ,iarrchi1
     & ,iarrchi2
     & ,iarrchi3
     & ,iarrchi4
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrlo2
      integer iarrlo3
      integer iarrlo4
      integer iarrhi0
      integer iarrhi1
      integer iarrhi2
      integer iarrhi3
      integer iarrhi4
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2,
     & iarrlo3:iarrhi3,
     & iarrlo4:iarrhi4)
      integer iarrclo0
      integer iarrclo1
      integer iarrclo2
      integer iarrclo3
      integer iarrclo4
      integer iarrchi0
      integer iarrchi1
      integer iarrchi2
      integer iarrchi3
      integer iarrchi4
      integer arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1,
     & iarrclo2:iarrchi2,
     & iarrclo3:iarrchi3,
     & iarrclo4:iarrchi4)
      integer i0, i1, i2, i3, i4
      integer inc
      inc = 1
      do i4 = iarrclo4, iarrchi4
         do i3 = iarrclo3, iarrchi3
            do i2 = iarrclo2, iarrchi2
               do i1 = iarrclo1, iarrchi1
                  do i0 = iarrclo0, iarrchi0
                     arr(i0, i1, i2, i3, i4) =
     & arrc(i0, i1, i2, i3, i4) + inc
                     inc = inc + 1
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine RANK6CHARRAY(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrlo2
     & ,iarrlo3
     & ,iarrlo4
     & ,iarrlo5
     & ,iarrhi0
     & ,iarrhi1
     & ,iarrhi2
     & ,iarrhi3
     & ,iarrhi4
     & ,iarrhi5
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrclo2
     & ,iarrclo3
     & ,iarrclo4
     & ,iarrclo5
     & ,iarrchi0
     & ,iarrchi1
     & ,iarrchi2
     & ,iarrchi3
     & ,iarrchi4
     & ,iarrchi5
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrlo2
      integer iarrlo3
      integer iarrlo4
      integer iarrlo5
      integer iarrhi0
      integer iarrhi1
      integer iarrhi2
      integer iarrhi3
      integer iarrhi4
      integer iarrhi5
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2,
     & iarrlo3:iarrhi3,
     & iarrlo4:iarrhi4,
     & iarrlo5:iarrhi5)
      integer iarrclo0
      integer iarrclo1
      integer iarrclo2
      integer iarrclo3
      integer iarrclo4
      integer iarrclo5
      integer iarrchi0
      integer iarrchi1
      integer iarrchi2
      integer iarrchi3
      integer iarrchi4
      integer iarrchi5
      integer arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1,
     & iarrclo2:iarrchi2,
     & iarrclo3:iarrchi3,
     & iarrclo4:iarrchi4,
     & iarrclo5:iarrchi5)
      integer i0, i1, i2, i3, i4, i5
      integer inc
      inc = 1
      do i5 = iarrclo5, iarrchi5
         do i4 = iarrclo4, iarrchi4
            do i3 = iarrclo3, iarrchi3
               do i2 = iarrclo2, iarrchi2
                  do i1 = iarrclo1, iarrchi1
                     do i0 = iarrclo0, iarrchi0
                        arr(i0, i1, i2, i3, i4, i5) =
     & arrc(i0, i1, i2, i3, i4, i5) + inc
                        inc = inc + 1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine RANK7CHARRAY(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrlo2
     & ,iarrlo3
     & ,iarrlo4
     & ,iarrlo5
     & ,iarrlo6
     & ,iarrhi0
     & ,iarrhi1
     & ,iarrhi2
     & ,iarrhi3
     & ,iarrhi4
     & ,iarrhi5
     & ,iarrhi6
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrclo2
     & ,iarrclo3
     & ,iarrclo4
     & ,iarrclo5
     & ,iarrclo6
     & ,iarrchi0
     & ,iarrchi1
     & ,iarrchi2
     & ,iarrchi3
     & ,iarrchi4
     & ,iarrchi5
     & ,iarrchi6
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrlo2
      integer iarrlo3
      integer iarrlo4
      integer iarrlo5
      integer iarrlo6
      integer iarrhi0
      integer iarrhi1
      integer iarrhi2
      integer iarrhi3
      integer iarrhi4
      integer iarrhi5
      integer iarrhi6
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2,
     & iarrlo3:iarrhi3,
     & iarrlo4:iarrhi4,
     & iarrlo5:iarrhi5,
     & iarrlo6:iarrhi6)
      integer iarrclo0
      integer iarrclo1
      integer iarrclo2
      integer iarrclo3
      integer iarrclo4
      integer iarrclo5
      integer iarrclo6
      integer iarrchi0
      integer iarrchi1
      integer iarrchi2
      integer iarrchi3
      integer iarrchi4
      integer iarrchi5
      integer iarrchi6
      integer arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1,
     & iarrclo2:iarrchi2,
     & iarrclo3:iarrchi3,
     & iarrclo4:iarrchi4,
     & iarrclo5:iarrchi5,
     & iarrclo6:iarrchi6)
      integer i0, i1, i2, i3, i4, i5, i6
      integer inc
      inc = 1
      do i6 = iarrclo6, iarrchi6
         do i5 = iarrclo5, iarrchi5
            do i4 = iarrclo4, iarrchi4
               do i3 = iarrclo3, iarrchi3
                  do i2 = iarrclo2, iarrchi2
                     do i1 = iarrclo1, iarrchi1
                        do i0 = iarrclo0, iarrchi0
                           arr(i0, i1, i2, i3, i4, i5, i6) =
     & arrc(i0, i1, i2, i3, i4, i5, i6) + inc
                           inc = inc + 1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine RANKCHARRAYSPACEDIMPLUS1CONT(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrlo2
     & ,iarrhi0
     & ,iarrhi1
     & ,iarrhi2
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrclo2
     & ,iarrchi0
     & ,iarrchi1
     & ,iarrchi2
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrlo2
      integer iarrhi0
      integer iarrhi1
      integer iarrhi2
      REAL*8 arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2)
      integer iarrclo0
      integer iarrclo1
      integer iarrclo2
      integer iarrchi0
      integer iarrchi1
      integer iarrchi2
      REAL*8 arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1,
     & iarrclo2:iarrchi2)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer c0
      REAL*8 inc;
      inc = (0.500d0)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         do c0 = iarrclo0, iarrchi0
            arr(c0, i0,i1) = arrc(c0, i0,i1) + inc
            inc = inc + (0.500d0)
         enddo
      enddo
      enddo
      return
      end
      subroutine RANKCHARRAYSPACEDIMPLUS1DIST(
     & arr
     & ,iarrlo0
     & ,iarrlo1
     & ,iarrlo2
     & ,iarrhi0
     & ,iarrhi1
     & ,iarrhi2
     & ,arrc
     & ,iarrclo0
     & ,iarrclo1
     & ,iarrclo2
     & ,iarrchi0
     & ,iarrchi1
     & ,iarrchi2
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dimc0
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iarrlo0
      integer iarrlo1
      integer iarrlo2
      integer iarrhi0
      integer iarrhi1
      integer iarrhi2
      REAL*8 arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2)
      integer iarrclo0
      integer iarrclo1
      integer iarrclo2
      integer iarrchi0
      integer iarrchi1
      integer iarrchi2
      REAL*8 arrc(
     & iarrclo0:iarrchi0,
     & iarrclo1:iarrchi1,
     & iarrclo2:iarrchi2)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dimc0
      integer i0,i1
      integer c0
      REAL*8 inc;
      inc = (0.500d0)
      do c0 = 0, dimc0-1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            arr(i0,i1, c0) = arrc(i0,i1, c0) + inc
            inc = inc + (0.500d0)
      enddo
      enddo
      enddo
      return
      end
