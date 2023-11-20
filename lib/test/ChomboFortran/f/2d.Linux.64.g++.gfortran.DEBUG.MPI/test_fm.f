      subroutine TEST_FM1(
     & a_i
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer a_i
      if( a_i .EQ. 1 ) goto 1
      if( a_i .EQ. 2 ) goto 2
      if( a_i .EQ. 3 ) goto 10
      if( a_i .EQ. 4 ) goto 11
      if( a_i .EQ. 5 ) goto 12345
      a_i = -1
      return
    1 continue
      a_i = 101
      return
    2 continue
      a_i = 102
      return
   10 continue
      a_i = 103
      return
   11 a_i = 104
      return
12345 a_i = 105
      return
      end
