#ifndef _CHFIOF_F_H_
#define _CHFIOF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCHF_FIO 
#define GUARDCHF_FIO 
// Prototype for Fortran procedure CHF_FIO ...
//
void FORTRAN_NAME( CHF_FIO ,chf_fio )(
      CHFp_CONST_INT(Int1)
      ,CHFp_CONST_REAL(Real1)
      ,CHFp_CONST_INT(proc)
      ,CHFp_INT(status) );

#define FORT_CHF_FIO FORTRAN_NAME( inlineCHF_FIO, inlineCHF_FIO)
#define FORTNT_CHF_FIO FORTRAN_NAME( CHF_FIO, chf_fio)

inline void FORTRAN_NAME(inlineCHF_FIO, inlineCHF_FIO)(
      CHFp_CONST_INT(Int1)
      ,CHFp_CONST_REAL(Real1)
      ,CHFp_CONST_INT(proc)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_CHF_FIO");
 FORTRAN_NAME( CHF_FIO ,chf_fio )(
      CHFt_CONST_INT(Int1)
      ,CHFt_CONST_REAL(Real1)
      ,CHFt_CONST_INT(proc)
      ,CHFt_INT(status) );
}
#endif  // GUARDCHF_FIO 

#ifndef GUARDCHF_STDIO 
#define GUARDCHF_STDIO 
// Prototype for Fortran procedure CHF_STDIO ...
//
void FORTRAN_NAME( CHF_STDIO ,chf_stdio )(
      CHFp_CONST_INT(Int1)
      ,CHFp_CONST_REAL(Real1)
      ,CHFp_INT(status) );

#define FORT_CHF_STDIO FORTRAN_NAME( inlineCHF_STDIO, inlineCHF_STDIO)
#define FORTNT_CHF_STDIO FORTRAN_NAME( CHF_STDIO, chf_stdio)

inline void FORTRAN_NAME(inlineCHF_STDIO, inlineCHF_STDIO)(
      CHFp_CONST_INT(Int1)
      ,CHFp_CONST_REAL(Real1)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_CHF_STDIO");
 FORTRAN_NAME( CHF_STDIO ,chf_stdio )(
      CHFt_CONST_INT(Int1)
      ,CHFt_CONST_REAL(Real1)
      ,CHFt_INT(status) );
}
#endif  // GUARDCHF_STDIO 

#ifndef GUARDCHF_MAYDAY 
#define GUARDCHF_MAYDAY 
// Prototype for Fortran procedure CHF_MAYDAY ...
//
void FORTRAN_NAME( CHF_MAYDAY ,chf_mayday )(
      CHFp_INT(status) );

#define FORT_CHF_MAYDAY FORTRAN_NAME( inlineCHF_MAYDAY, inlineCHF_MAYDAY)
#define FORTNT_CHF_MAYDAY FORTRAN_NAME( CHF_MAYDAY, chf_mayday)

inline void FORTRAN_NAME(inlineCHF_MAYDAY, inlineCHF_MAYDAY)(
      CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_CHF_MAYDAY");
 FORTRAN_NAME( CHF_MAYDAY ,chf_mayday )(
      CHFt_INT(status) );
}
#endif  // GUARDCHF_MAYDAY 

}

#endif
