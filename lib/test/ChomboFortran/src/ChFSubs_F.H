#ifndef _CHFSUBS_F_H_
#define _CHFSUBS_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDTESTVIRC_1 
#define GUARDTESTVIRC_1 
// Prototype for Fortran procedure testVIRC_1 ...
//
void FORTRAN_NAME( TESTVIRC_1 ,testvirc_1 )(
      CHFp_VI(vectorI)
      ,CHFp_CONST_VI(VectorCI)
      ,CHFp_VR(vectorR)
      ,CHFp_CONST_VR(VectorCR)
      ,CHFp_VC(vCmplx)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) );

#define FORT_TESTVIRC_1 FORTRAN_NAME( inlineTESTVIRC_1, inlineTESTVIRC_1)
#define FORTNT_TESTVIRC_1 FORTRAN_NAME( TESTVIRC_1, testvirc_1)

inline void FORTRAN_NAME(inlineTESTVIRC_1, inlineTESTVIRC_1)(
      CHFp_VI(vectorI)
      ,CHFp_CONST_VI(VectorCI)
      ,CHFp_VR(vectorR)
      ,CHFp_CONST_VR(VectorCR)
      ,CHFp_VC(vCmplx)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TESTVIRC_1");
 FORTRAN_NAME( TESTVIRC_1 ,testvirc_1 )(
      CHFt_VI(vectorI)
      ,CHFt_CONST_VI(VectorCI)
      ,CHFt_VR(vectorR)
      ,CHFt_CONST_VR(VectorCR)
      ,CHFt_VC(vCmplx)
      ,CHFt_CONST_REAL(Eps)
      ,CHFt_INT(status) );
}
#endif  // GUARDTESTVIRC_1 

#ifndef GUARDTESTIR1D_1 
#define GUARDTESTIR1D_1 
// Prototype for Fortran procedure TESTIR1D_1 ...
//
void FORTRAN_NAME( TESTIR1D_1 ,testir1d_1 )(
      CHFp_I1D(arrayI,iarrayIhi0)
      ,CHFp_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFp_R1D(arrayR,iarrayRhi0)
      ,CHFp_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) );

#define FORT_TESTIR1D_1 FORTRAN_NAME( inlineTESTIR1D_1, inlineTESTIR1D_1)
#define FORTNT_TESTIR1D_1 FORTRAN_NAME( TESTIR1D_1, testir1d_1)

inline void FORTRAN_NAME(inlineTESTIR1D_1, inlineTESTIR1D_1)(
      CHFp_I1D(arrayI,iarrayIhi0)
      ,CHFp_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFp_R1D(arrayR,iarrayRhi0)
      ,CHFp_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TESTIR1D_1");
 FORTRAN_NAME( TESTIR1D_1 ,testir1d_1 )(
      CHFt_I1D(arrayI,iarrayIhi0)
      ,CHFt_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFt_R1D(arrayR,iarrayRhi0)
      ,CHFt_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFt_CONST_REAL(Eps)
      ,CHFt_INT(status) );
}
#endif  // GUARDTESTIR1D_1 

#ifndef GUARDTESTIR1D_2 
#define GUARDTESTIR1D_2 
// Prototype for Fortran procedure TESTIR1D_2 ...
//
void FORTRAN_NAME( TESTIR1D_2 ,testir1d_2 )(
      CHFp_I1D(arrayI,iarrayIhi0)
      ,CHFp_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFp_R1D(arrayR,iarrayRhi0)
      ,CHFp_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) );

#define FORT_TESTIR1D_2 FORTRAN_NAME( inlineTESTIR1D_2, inlineTESTIR1D_2)
#define FORTNT_TESTIR1D_2 FORTRAN_NAME( TESTIR1D_2, testir1d_2)

inline void FORTRAN_NAME(inlineTESTIR1D_2, inlineTESTIR1D_2)(
      CHFp_I1D(arrayI,iarrayIhi0)
      ,CHFp_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFp_R1D(arrayR,iarrayRhi0)
      ,CHFp_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TESTIR1D_2");
 FORTRAN_NAME( TESTIR1D_2 ,testir1d_2 )(
      CHFt_I1D(arrayI,iarrayIhi0)
      ,CHFt_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFt_R1D(arrayR,iarrayRhi0)
      ,CHFt_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFt_CONST_REAL(Eps)
      ,CHFt_INT(status) );
}
#endif  // GUARDTESTIR1D_2 

#ifndef GUARDTESTDSEL_1 
#define GUARDTESTDSEL_1 
// Prototype for Fortran procedure TESTDSEL_1 ...
//
void FORTRAN_NAME( TESTDSEL_1 ,testdsel_1 )(
      CHFp_CONST_INT(SpaceDim)
      ,CHFp_CONST_INT(Verbose)
      ,CHFp_INT(status) );

#define FORT_TESTDSEL_1 FORTRAN_NAME( inlineTESTDSEL_1, inlineTESTDSEL_1)
#define FORTNT_TESTDSEL_1 FORTRAN_NAME( TESTDSEL_1, testdsel_1)

inline void FORTRAN_NAME(inlineTESTDSEL_1, inlineTESTDSEL_1)(
      CHFp_CONST_INT(SpaceDim)
      ,CHFp_CONST_INT(Verbose)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TESTDSEL_1");
 FORTRAN_NAME( TESTDSEL_1 ,testdsel_1 )(
      CHFt_CONST_INT(SpaceDim)
      ,CHFt_CONST_INT(Verbose)
      ,CHFt_INT(status) );
}
#endif  // GUARDTESTDSEL_1 

#ifndef GUARDTESTMDO_1 
#define GUARDTESTMDO_1 
// Prototype for Fortran procedure TESTMDO_1 ...
//
void FORTRAN_NAME( TESTMDO_1 ,testmdo_1 )(
      CHFp_FIA1(a)
      ,CHFp_BOX(ba)
      ,CHFp_CONST_INT(Vol)
      ,CHFp_CONST_INT(Verbose)
      ,CHFp_INT(status) );

#define FORT_TESTMDO_1 FORTRAN_NAME( inlineTESTMDO_1, inlineTESTMDO_1)
#define FORTNT_TESTMDO_1 FORTRAN_NAME( TESTMDO_1, testmdo_1)

inline void FORTRAN_NAME(inlineTESTMDO_1, inlineTESTMDO_1)(
      CHFp_FIA1(a)
      ,CHFp_BOX(ba)
      ,CHFp_CONST_INT(Vol)
      ,CHFp_CONST_INT(Verbose)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TESTMDO_1");
 FORTRAN_NAME( TESTMDO_1 ,testmdo_1 )(
      CHFt_FIA1(a)
      ,CHFt_BOX(ba)
      ,CHFt_CONST_INT(Vol)
      ,CHFt_CONST_INT(Verbose)
      ,CHFt_INT(status) );
}
#endif  // GUARDTESTMDO_1 

#ifndef GUARDTESTAMDO_1 
#define GUARDTESTAMDO_1 
// Prototype for Fortran procedure TESTAMDO_1 ...
//
void FORTRAN_NAME( TESTAMDO_1 ,testamdo_1 )(
      CHFp_FIA1(a)
      ,CHFp_BOX(ba)
      ,CHFp_CONST_INT(Vol)
      ,CHFp_BOX(offsetbox)
      ,CHFp_CONST_INT(Verbose)
      ,CHFp_INT(status) );

#define FORT_TESTAMDO_1 FORTRAN_NAME( inlineTESTAMDO_1, inlineTESTAMDO_1)
#define FORTNT_TESTAMDO_1 FORTRAN_NAME( TESTAMDO_1, testamdo_1)

inline void FORTRAN_NAME(inlineTESTAMDO_1, inlineTESTAMDO_1)(
      CHFp_FIA1(a)
      ,CHFp_BOX(ba)
      ,CHFp_CONST_INT(Vol)
      ,CHFp_BOX(offsetbox)
      ,CHFp_CONST_INT(Verbose)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TESTAMDO_1");
 FORTRAN_NAME( TESTAMDO_1 ,testamdo_1 )(
      CHFt_FIA1(a)
      ,CHFt_BOX(ba)
      ,CHFt_CONST_INT(Vol)
      ,CHFt_BOX(offsetbox)
      ,CHFt_CONST_INT(Verbose)
      ,CHFt_INT(status) );
}
#endif  // GUARDTESTAMDO_1 

#ifndef GUARDTEST_CHFID 
#define GUARDTEST_CHFID 
// Prototype for Fortran procedure TEST_CHFID ...
//
void FORTRAN_NAME( TEST_CHFID ,test_chfid )(
      CHFp_INT(status) );

#define FORT_TEST_CHFID FORTRAN_NAME( inlineTEST_CHFID, inlineTEST_CHFID)
#define FORTNT_TEST_CHFID FORTRAN_NAME( TEST_CHFID, test_chfid)

inline void FORTRAN_NAME(inlineTEST_CHFID, inlineTEST_CHFID)(
      CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TEST_CHFID");
 FORTRAN_NAME( TEST_CHFID ,test_chfid )(
      CHFt_INT(status) );
}
#endif  // GUARDTEST_CHFID 

#ifndef GUARDTESTCALL 
#define GUARDTESTCALL 
// Prototype for Fortran procedure TESTCALL ...
//
void FORTRAN_NAME( TESTCALL ,testcall )(
      CHFp_I1D(arrayI,iarrayIhi0)
      ,CHFp_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFp_R1D(arrayR,iarrayRhi0)
      ,CHFp_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) );

#define FORT_TESTCALL FORTRAN_NAME( inlineTESTCALL, inlineTESTCALL)
#define FORTNT_TESTCALL FORTRAN_NAME( TESTCALL, testcall)

inline void FORTRAN_NAME(inlineTESTCALL, inlineTESTCALL)(
      CHFp_I1D(arrayI,iarrayIhi0)
      ,CHFp_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFp_R1D(arrayR,iarrayRhi0)
      ,CHFp_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFp_CONST_REAL(Eps)
      ,CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TESTCALL");
 FORTRAN_NAME( TESTCALL ,testcall )(
      CHFt_I1D(arrayI,iarrayIhi0)
      ,CHFt_CONST_I1D(arrayCI,iarrayCIhi0)
      ,CHFt_R1D(arrayR,iarrayRhi0)
      ,CHFt_CONST_R1D(arrayCR,iarrayCRhi0)
      ,CHFt_CONST_REAL(Eps)
      ,CHFt_INT(status) );
}
#endif  // GUARDTESTCALL 

#ifndef GUARDTEST_PARSE1 
#define GUARDTEST_PARSE1 
// Prototype for Fortran procedure TEST_PARSE1 ...
//
void FORTRAN_NAME( TEST_PARSE1 ,test_parse1 )(
      CHFp_INT(status) );

#define FORT_TEST_PARSE1 FORTRAN_NAME( inlineTEST_PARSE1, inlineTEST_PARSE1)
#define FORTNT_TEST_PARSE1 FORTRAN_NAME( TEST_PARSE1, test_parse1)

inline void FORTRAN_NAME(inlineTEST_PARSE1, inlineTEST_PARSE1)(
      CHFp_INT(status) )
{
 CH_TIMELEAF("FORT_TEST_PARSE1");
 FORTRAN_NAME( TEST_PARSE1 ,test_parse1 )(
      CHFt_INT(status) );
}
#endif  // GUARDTEST_PARSE1 

#ifndef GUARDTEST_SHIFT_FIA 
#define GUARDTEST_SHIFT_FIA 
// Prototype for Fortran procedure TEST_SHIFT_FIA ...
//
void FORTRAN_NAME( TEST_SHIFT_FIA ,test_shift_fia )(
      CHFp_FIA(crFIA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FIA(fnFIA)
      ,CHFp_CONST_INT(refRatio) );

#define FORT_TEST_SHIFT_FIA FORTRAN_NAME( inlineTEST_SHIFT_FIA, inlineTEST_SHIFT_FIA)
#define FORTNT_TEST_SHIFT_FIA FORTRAN_NAME( TEST_SHIFT_FIA, test_shift_fia)

inline void FORTRAN_NAME(inlineTEST_SHIFT_FIA, inlineTEST_SHIFT_FIA)(
      CHFp_FIA(crFIA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FIA(fnFIA)
      ,CHFp_CONST_INT(refRatio) )
{
 CH_TIMELEAF("FORT_TEST_SHIFT_FIA");
 FORTRAN_NAME( TEST_SHIFT_FIA ,test_shift_fia )(
      CHFt_FIA(crFIA)
      ,CHFt_BOX(fnBox)
      ,CHFt_CONST_FIA(fnFIA)
      ,CHFt_CONST_INT(refRatio) );
}
#endif  // GUARDTEST_SHIFT_FIA 

#ifndef GUARDTEST_SHIFT_FIA1 
#define GUARDTEST_SHIFT_FIA1 
// Prototype for Fortran procedure TEST_SHIFT_FIA1 ...
//
void FORTRAN_NAME( TEST_SHIFT_FIA1 ,test_shift_fia1 )(
      CHFp_FIA1(crFIA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FIA1(fnFIA)
      ,CHFp_CONST_INT(refRatio) );

#define FORT_TEST_SHIFT_FIA1 FORTRAN_NAME( inlineTEST_SHIFT_FIA1, inlineTEST_SHIFT_FIA1)
#define FORTNT_TEST_SHIFT_FIA1 FORTRAN_NAME( TEST_SHIFT_FIA1, test_shift_fia1)

inline void FORTRAN_NAME(inlineTEST_SHIFT_FIA1, inlineTEST_SHIFT_FIA1)(
      CHFp_FIA1(crFIA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FIA1(fnFIA)
      ,CHFp_CONST_INT(refRatio) )
{
 CH_TIMELEAF("FORT_TEST_SHIFT_FIA1");
 FORTRAN_NAME( TEST_SHIFT_FIA1 ,test_shift_fia1 )(
      CHFt_FIA1(crFIA)
      ,CHFt_BOX(fnBox)
      ,CHFt_CONST_FIA1(fnFIA)
      ,CHFt_CONST_INT(refRatio) );
}
#endif  // GUARDTEST_SHIFT_FIA1 

#ifndef GUARDTEST_SHIFT_FRA 
#define GUARDTEST_SHIFT_FRA 
// Prototype for Fortran procedure TEST_SHIFT_FRA ...
//
void FORTRAN_NAME( TEST_SHIFT_FRA ,test_shift_fra )(
      CHFp_FRA(crFRA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FRA(fnFRA)
      ,CHFp_CONST_INT(refRatio) );

#define FORT_TEST_SHIFT_FRA FORTRAN_NAME( inlineTEST_SHIFT_FRA, inlineTEST_SHIFT_FRA)
#define FORTNT_TEST_SHIFT_FRA FORTRAN_NAME( TEST_SHIFT_FRA, test_shift_fra)

inline void FORTRAN_NAME(inlineTEST_SHIFT_FRA, inlineTEST_SHIFT_FRA)(
      CHFp_FRA(crFRA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FRA(fnFRA)
      ,CHFp_CONST_INT(refRatio) )
{
 CH_TIMELEAF("FORT_TEST_SHIFT_FRA");
 FORTRAN_NAME( TEST_SHIFT_FRA ,test_shift_fra )(
      CHFt_FRA(crFRA)
      ,CHFt_BOX(fnBox)
      ,CHFt_CONST_FRA(fnFRA)
      ,CHFt_CONST_INT(refRatio) );
}
#endif  // GUARDTEST_SHIFT_FRA 

#ifndef GUARDTEST_SHIFT_FRA1 
#define GUARDTEST_SHIFT_FRA1 
// Prototype for Fortran procedure TEST_SHIFT_FRA1 ...
//
void FORTRAN_NAME( TEST_SHIFT_FRA1 ,test_shift_fra1 )(
      CHFp_FRA1(crFRA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FRA1(fnFRA)
      ,CHFp_CONST_INT(refRatio) );

#define FORT_TEST_SHIFT_FRA1 FORTRAN_NAME( inlineTEST_SHIFT_FRA1, inlineTEST_SHIFT_FRA1)
#define FORTNT_TEST_SHIFT_FRA1 FORTRAN_NAME( TEST_SHIFT_FRA1, test_shift_fra1)

inline void FORTRAN_NAME(inlineTEST_SHIFT_FRA1, inlineTEST_SHIFT_FRA1)(
      CHFp_FRA1(crFRA)
      ,CHFp_BOX(fnBox)
      ,CHFp_CONST_FRA1(fnFRA)
      ,CHFp_CONST_INT(refRatio) )
{
 CH_TIMELEAF("FORT_TEST_SHIFT_FRA1");
 FORTRAN_NAME( TEST_SHIFT_FRA1 ,test_shift_fra1 )(
      CHFt_FRA1(crFRA)
      ,CHFt_BOX(fnBox)
      ,CHFt_CONST_FRA1(fnFRA)
      ,CHFt_CONST_INT(refRatio) );
}
#endif  // GUARDTEST_SHIFT_FRA1 

}

#endif
