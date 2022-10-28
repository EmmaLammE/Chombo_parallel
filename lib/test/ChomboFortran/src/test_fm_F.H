#ifndef _TEST_FM_F_H_
#define _TEST_FM_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDTEST_FM1 
#define GUARDTEST_FM1 
// Prototype for Fortran procedure test_FM1 ...
//
void FORTRAN_NAME( TEST_FM1 ,test_fm1 )(
      CHFp_INT(a_i) );

#define FORT_TEST_FM1 FORTRAN_NAME( inlineTEST_FM1, inlineTEST_FM1)
#define FORTNT_TEST_FM1 FORTRAN_NAME( TEST_FM1, test_fm1)

inline void FORTRAN_NAME(inlineTEST_FM1, inlineTEST_FM1)(
      CHFp_INT(a_i) )
{
 CH_TIMELEAF("FORT_TEST_FM1");
 FORTRAN_NAME( TEST_FM1 ,test_fm1 )(
      CHFt_INT(a_i) );
}
#endif  // GUARDTEST_FM1 

}

#endif
