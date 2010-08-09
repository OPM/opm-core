#include <mex.h>
#include <umfpack.h>
#include "call_umfpack.h"

void sortRows(int n, mwSignedIndex* ia, mwSignedIndex* ja, 
              double *sa)
{
   int i,j,k,row;

   for(row=0; row<n; ++row)
   {
      for(k=ia[row]; k<ia[row+1]; ++k)
      {
         for(i=k+1; i<ia[row+1]; ++i)
         {
            if(ja[k] > ja[i])
            {
               mwSignedIndex itmp = ja[i];
               ja[i] = ja[k];
               ja[k] = itmp;
               
               double dtmp = sa[i];
               sa[i] = sa[k]; 
               sa[k] = dtmp;               
            }
         }
      }
   }
}


/*---------------------------------------------------------------------------*/
void
callMWUMFPACK(int n, mwSignedIndex* ia, mwSignedIndex* ja, 
             double *sa, double *b, double *x)
/*---------------------------------------------------------------------------*/
{
   sortRows(n, ia, ja, sa);
   /* double *null = (double *) NULL ; */
   void *Symbolic, *Numeric ;
#if 0
   long   N = 5 ;
   long   IA [ ] = {0, 2, 5, 9, 10, 12} ;
   long   JA [ ] = { 0,  1,  0,   2,  4,  1,  2,  3,   4,  2,  1,  4} ;
   double SA [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
   double B [ ] = {8., 45., -3., 3., 19.} ;
   double X [5] ;
#endif
   double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];
   
   umfpack_dl_defaults(Control);
#if 0
   umfpack_dl_symbolic (N, N, IA, JA, SA, &Symbolic, Control, Info) ;
   umfpack_dl_numeric  (IA, JA, SA, Symbolic, &Numeric, Control, Info) ;
   umfpack_dl_free_symbolic (&Symbolic);   
   umfpack_dl_solve (UMFPACK_A, IA,JA,SA, &X[0], &B[0], Numeric, Control, Info) ;
   umfpack_dl_free_numeric (&Numeric);
#else
#if 1
   int i;
   for(i=0; i<n+1; ++i)
   {
      mexPrintf("ia[%d] = %d\n", i, ia[i]);
   }
   for(i=0; i<ia[n]; ++i)
   {
      mexPrintf("ja[%d] = %d\n", i, ja[i]);
   }
   for(i=0; i<ia[n]; ++i)
   {
      mexPrintf("sa[%d] = %f\n", i, sa[i]);
   }
   for(i=0; i<n; ++i)
   {
      mexPrintf("b[%d] = %f\n", i, b[i]);
   }
#endif
   umfpack_dl_symbolic (n, n, ia, ja, sa, &Symbolic, Control, Info) ;
   umfpack_dl_numeric  (ia, ja, sa, Symbolic, &Numeric, Control, Info) ;
   umfpack_dl_free_symbolic (&Symbolic);   
   umfpack_dl_solve (UMFPACK_A, ia, ja, sa, &x[0], &b[0], Numeric, Control, Info) ;
   umfpack_dl_free_numeric (&Numeric);
#endif
}
