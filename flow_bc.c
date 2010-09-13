#include <stdlib.h>

#include "flow_bc.h"


/* ---------------------------------------------------------------------- */
flowbc_t *
allocate_flowbc(size_t nf)
/* ---------------------------------------------------------------------- */
{
    size_t    i;
    flowbc_t *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->type  = malloc(nf * sizeof *new->type);
        new->bcval = malloc(nf * sizeof *new->bcval);

        if ((new->type == NULL) || (new->bcval == NULL)) {
            deallocate_flowbc(new);
            new = NULL;
        } else {
            for (i = 0; i < nf; i++) {
                new->type [i] = UNSET;
                new->bcval[i] = 0.0;
            }
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
void
deallocate_flowbc(flowbc_t *fbc)
/* ---------------------------------------------------------------------- */
{
    if (fbc != NULL) {
        free(fbc->bcval);
        free(fbc->type );
    }

    free(fbc);
}
