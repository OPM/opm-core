/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef NLSOLVERS_H
#define NLSOLVERS_H

double bisection   (double (*)(double, void*), void*, double, int, int*);
double ridder      (double (*)(double, void*), void*, double, int, int*);
double regulafalsi (double (*)(double, void*), void*, double, int, int*);

#endif /* NLSOLVERS_H */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
