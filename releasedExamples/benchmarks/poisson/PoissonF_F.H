#ifndef _POISSONF_F_H_
#define _POISSONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSIMPLESTLAPLACIAN 
#define GUARDSIMPLESTLAPLACIAN 
// Prototype for Fortran procedure simplestlaplacian ...
//
void FORTRAN_NAME( SIMPLESTLAPLACIAN ,simplestlaplacian )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) );

#define FORT_SIMPLESTLAPLACIAN FORTRAN_NAME( inlineSIMPLESTLAPLACIAN, inlineSIMPLESTLAPLACIAN)
#define FORTNT_SIMPLESTLAPLACIAN FORTRAN_NAME( SIMPLESTLAPLACIAN, simplestlaplacian)

inline void FORTRAN_NAME(inlineSIMPLESTLAPLACIAN, inlineSIMPLESTLAPLACIAN)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SIMPLESTLAPLACIAN");
 FORTRAN_NAME( SIMPLESTLAPLACIAN ,simplestlaplacian )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSIMPLESTLAPLACIAN 

}

#endif
