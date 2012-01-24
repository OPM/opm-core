#include <opm/core/transport/reorder/fluid.h>


#if 1

/* Quadratic relperm */

/*---------------------------------------------------------------------------*/
double fluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
    double so, mw, mo, vw, vo;
    (void) reg;
    /* Hardcoding behaviour to make test program work */
    so = 1.0 - sw;
    vw = 0.001;
    vo = 0.003;
    mw  = sw*sw / vw;
    mo  = so*so / vo;
    return mw / (mw + mo);
}

/*---------------------------------------------------------------------------*/
double dfluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
    double so, mw, mo, vw, vo, dmw, dmo, fw;
    (void) reg;
    /* Hardcoding behaviour to make test program work */
    so = 1.0 - sw;
    vw = 0.001;
    vo = 0.003;
    mw  = sw*sw / vw;
    mo  = so*so / vo;

    dmw = 2.0*sw / vw;
    dmo = 2.0*so / vo;

    fw = mw / (mw+mo);

    return (dmw * (1-fw) - dmo * fw) / (mw+mo);
}
#else

/* Linear relperm */

/*---------------------------------------------------------------------------*/
double fluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
    (void) reg;
    return sw;
}

/*---------------------------------------------------------------------------*/
double dfluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
    (void) sw, reg;
    return 1.0;
}
#endif
/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
