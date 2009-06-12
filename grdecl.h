#ifndef GRDECL_H
#define GRDECL_H

struct Grdecl{
  int           n;
  int           dims[3];
  const double *coord;
  const double *zcorn;
  const int    *actnum;
};


#endif
