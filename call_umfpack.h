#ifndef CALL_UMFPACK_H_INCLUDED
#define CALL_UMFPACK_H_INCLUDED
void
callMWUMFPACK(int n, mwSignedIndex* ia, mwSignedIndex* ja, 
              double *sa, double *b, double *x);
#endif
