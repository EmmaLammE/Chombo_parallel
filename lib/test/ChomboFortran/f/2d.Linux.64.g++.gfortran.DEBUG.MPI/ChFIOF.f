      subroutine CHF_FIO(
     & Int1
     & ,Real1
     & ,proc
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer Int1
      REAL*8 Real1
      integer proc
      integer status
      integer unit
      integer int2
      REAL*8 real2
      integer curdigit
      integer curnum
      CHARACTER*11 digit
      CHARACTER*100 procstr
      integer num
      integer index
      CHARACTER*100 filename
      digit = '0123456789'
      num = 0
      curnum = proc
  100 curdigit = mod(curnum,10) + 1
      num = num + 1
      curnum = curnum / 10
      if (curnum.gt.0 .and. num.lt.100) goto 100
      index = num
      curnum = proc
  200 curdigit = mod(curnum,10) + 1
      procstr(index:index) = digit(curdigit:curdigit)
      index = index - 1
      curnum = curnum / 10
      if (curnum.gt.0) goto 200
      filename = 'CHFIOtest'
      unit = 19 + proc
      open( unit ,file=filename ,FORM='FORMATTED'
     & ,iostat=status )
      if( status .NE. 0 )then
        write(*,*) 'error: CHF_FIO: open failed, status = ',status
        return
      endif
      write( unit,*,iostat=status ) Int1 ,' is the integer.'
      write( unit,*,iostat=status ) Real1 ,' is the real.'
      close( unit ,status='KEEP' ,iostat=status )
      open( unit ,file=filename ,status='OLD' ,FORM='FORMATTED'
     & ,iostat=status )
      if( status .NE. 0 )then
        write(*,*) 'error: CHF_FIO: open/old failed, status = ',status
        return
      endif
      read( unit,*,iostat=status) int2
      read( unit,*,iostat=status) real2
      if( int2 .NE. Int1 ) then
         status = 3
         return
      endif
      if( real2 .NE. Real1 ) then
         status = 4
         return
      endif
      close( unit ,status='DELETE' ,iostat=status )
      status = 0
      return
      status = 5
      return
      end
      subroutine CHF_STDIO(
     & Int1
     & ,Real1
     & ,status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer Int1
      REAL*8 Real1
      integer status
      character indent*4
      indent = '    '
      print*,indent,'CHF_STDIO: test 1, integer, real values = ' ,Int1 ,
     &Real1
      status = 0
      return
      end
      subroutine CHF_MAYDAY(
     & status
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer status
      call MAYDAY_ERROR()
      call MAYDAY_ABORT()
      status = -123456
      return
      end
