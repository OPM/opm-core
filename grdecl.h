#ifndef GRDECL_H
#define GRDECL_H

struct Grdecl{
  int     n;
  int     dims[3];
  double *coord;
  double *zcorn;
  int    *actnum;
};


#endif
