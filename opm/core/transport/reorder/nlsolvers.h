/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef NLSOLVERS_H
#define NLSOLVERS_H

#ifdef __cplusplus
extern "C" {
#endif

struct NonlinearSolverCtrl
{
   enum     Method {RIDDERS, REGULAFALSI, BISECTION} method;
   double  nltolerance;
   int     maxiterations;
   double  initialguess;
   int     iterations;    /* set by solver */
   double  residual;      /* set by solver */
};

double find_zero  (double (*G)(double, void*), void *data, struct NonlinearSolverCtrl *ctrl);
double bisection   (double (*)(double, void*), void*, struct NonlinearSolverCtrl *ctrl);
double ridders     (double (*)(double, void*), void*, struct NonlinearSolverCtrl *ctrl);
double regulafalsi (double (*)(double, void*), void*, struct NonlinearSolverCtrl *ctrl);

#ifdef __cplusplus
}
#endif

#endif /* NLSOLVERS_H */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
