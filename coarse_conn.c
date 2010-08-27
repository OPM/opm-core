#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define GOLDEN_RAT (0.6180339887498949) /* (sqrt(5) - 1) / 2 */
#define IS_POW2(x) (((x) & ((x) - 1)) == 0)

struct hash_set {
    size_t  m;                  /* Table/set capacity (1<<p for some p) */
    int    *s;                  /* Set representation */
};

struct block_neighbour {
    int              b;         /* Neighbouring block */
    struct hash_set *fconns;    /* Constituent connections */
};

struct block_neighbours {
    int                    nneigh;  /* Number of neighbours. */
    int                    cpty;    /* Neighbour capacity. */
    struct block_neighbour **neigh; /* Actual neighbours (sorted on neigh[i]->b) */
};

/* ---------------------------------------------------------------------- */
static size_t
hash_set_size(size_t m)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    if (m == 0) {
        return 1;
    }

    if (IS_POW2(m)) {
        return m;
    }

    /* General case.  Use next power of two. */
    /* Algorithm due to
     *
     *    Warren Jr., Henry S. (2002). Hacker's Delight.
     *    Addison Wesley. pp. 48. ISBN 978-0201914658
     *
     * by way of Wikipedia. */
    m -= 1;
    for (i = 1; i < CHAR_BIT * sizeof m; i <<= 1) {
        m = m | (m >> i);
    }

    return m + 1;
}


/* ---------------------------------------------------------------------- */
static size_t
hash_set_idx(int k, size_t m)
/* ---------------------------------------------------------------------- */
{
    return floor(m * fmod(k * GOLDEN_RAT, 1.0));
}


/* ---------------------------------------------------------------------- */
static size_t
hash_set_insert_core(int k, size_t m, int *s)
/* ---------------------------------------------------------------------- */
{
    size_t h1, h2, i, j;

    assert ((0 < m) && (m < (size_t)(-1)));
    assert (IS_POW2(m));

    j = h1 = hash_set_idx(k, m);
    assert (h1 < m);

    if (s[j] == -1) { s[j] = k; }
    if (s[j] ==  k) { return j; }

    /* Double hash probing.  h2 relatively prime to 'm' */
    h2 = 2 * hash_set_idx(k, (m << 1) - 1) - 1;

    for (i = 1; (s[j] != -1) && (s[j] != k) && (i < m); i++) {
        j += h2;
        j &= m - 1;             /* Modulo m since IS_POW2(m). */
    }

    if (i < m) {
        s[j] = k;               /* Possibly no-op. */
    } else {
        j = m + 1;              /* Invalid.  Caveat emptor. */
    }

    return j;
}


/* ---------------------------------------------------------------------- */
static int
hash_set_expand(size_t m, struct hash_set *t)
/* ---------------------------------------------------------------------- */
{
    int ret, *s, *p;
    size_t i;

    assert (m > t->m);

    s = malloc(m * sizeof *s);
    if (s != NULL) {
        memset(s, -1, m * sizeof *s);

        for (i = 0; i < t->m; i++) {
            ret = hash_set_insert_core(t->s[i], m, s);
            assert (ret == t->s[i]);
        }

        p    = t->s;
        t->s = s;
        t->m = m;

        free(p);

        ret = m;
    } else {
        ret = -1;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
hash_set_deallocate(struct hash_set *t)
/* ---------------------------------------------------------------------- */
{
    if (t != NULL) {
        free(t->s);
    }

    free(t);
}


/* ---------------------------------------------------------------------- */
static struct hash_set *
hash_set_allocate(int m)
/* ---------------------------------------------------------------------- */
{
    size_t            sz;
    struct hash_set *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        sz = hash_set_size(m);
        new->s = malloc(sz * sizeof *new->s);

        if (new->s == NULL) {
            hash_set_deallocate(new);

            new = NULL;
        } else {
            memset(new->s, -1, sz * sizeof *new->s);
            new->m = sz;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static int
hash_set_insert(int k, struct hash_set *t)
/* ---------------------------------------------------------------------- */
{
    int ret;
    size_t i;

    assert (k >= 0);
    assert (t != NULL);
    assert (IS_POW2(t->m));

    i = hash_set_insert_core(k, t->m, t->s);
    if (i == t->m + 1) {
        /* Table full.  Preferable an infrequent occurrence.  Expand
         * table and re-insert key (if possible). */
        ret = hash_set_expand(t->m << 1, t);

        if (ret > 0) {
            i = hash_set_insert_core(k, t->m, t->s);
            assert (i < t->m);

            ret = k;
        }
    } else {
        ret = k;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
block_neighbour_deallocate(struct block_neighbour *bn)
/* ---------------------------------------------------------------------- */
{
    if (bn != NULL) {
        hash_set_deallocate(bn->fconns);
    }

    free(bn);
}


/* ---------------------------------------------------------------------- */
static struct block_neighbour *
block_neighbour_allocate(int nconn)
/* ---------------------------------------------------------------------- */
{
    struct block_neighbour *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        if (nconn > 0) {
            new->fconns = hash_set_allocate(nconn);

            if (new->fconns != NULL) {
                new->b  = -1;
            } else {
                block_neighbour_deallocate(new);
                new     = NULL;
            }
        } else {
            new->b      = -1;
            new->fconns = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static int
block_neighbour_insert_fconn(int fconn, struct block_neighbour *bn)
/* ---------------------------------------------------------------------- */
{
    int ret;

    assert (bn != NULL);

    ret = 0;

    if (bn->fconns != NULL) {
        ret = hash_set_insert(fconn, bn->fconns);
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
block_neighbours_deallocate(struct block_neighbours *bns)
/* ---------------------------------------------------------------------- */
{
    int i;

    if (bns != NULL) {
        if (bns->neigh != NULL) {
            for (i = bns->nneigh - 1; i >= 0; i--) {
                block_neighbour_deallocate(bns->neigh[i]);
            }
        }

        free(bns->neigh);
    }

    free(bns);
}


/* ---------------------------------------------------------------------- */
static struct block_neighbours *
block_neighbours_allocate(int nneigh)
/* ---------------------------------------------------------------------- */
{
    struct block_neighbours *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        if (nneigh > 0) {
            new->neigh = malloc(nneigh * sizeof *new->neigh);

            if (new->neigh != NULL) {
                new->nneigh = 0;
                new->cpty   = nneigh;
            } else {
                block_neighbours_deallocate(new);
                new = NULL;
            }
        } else {
            new->nneigh = 0;
            new->cpty   = 0;
            new->neigh  = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static int
block_neighbours_expand(int nneigh, struct block_neighbours *bns)
/* ---------------------------------------------------------------------- */
{
    int                      ret;
    struct block_neighbour **neigh;

    assert (bns != NULL);

    neigh = realloc(bns->neigh, nneigh * sizeof *neigh);

    if (neigh != NULL) {
        bns->neigh = neigh;
        bns->cpty  = nneigh;
        ret        = nneigh;
    } else {
        ret = -1;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static int
block_neighbours_insert_neighbour(int b, int fconn,
                                  struct block_neighbours *bns)
/* ---------------------------------------------------------------------- */
{
    assert (bns != NULL);

    
}
