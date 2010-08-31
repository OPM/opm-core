#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "coarse_conn.h"

/* ======================================================================
 * Macros
 * ====================================================================== */
#define GOLDEN_RAT (0.6180339887498949)              /* (sqrt(5) - 1) / 2 */
#define IS_POW2(x) (((x) & ((x) - 1)) == 0)
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))
#define MIN(a,b)   (-MAX(-(a), -(b)))

/* ======================================================================
 * Data structures
 * ====================================================================== */

/* Poor-man's unordered set (ind. key insert/all key extract only). */
struct hash_set {
    size_t  m;                  /* Table/set capacity (1<<p for some p) */
    int    *s;                  /* Set representation */
};


/* Individual block connection. */
struct block_neighbour {
    int              b;         /* Neighbouring block */
    struct hash_set *fconns;    /* Constituent connections */
};


/* Adjacency list of single block (directed graph) */
struct block_neighbours {
    int                    nneigh;  /* Number of neighbours. */
    int                    cpty;    /* Neighbour capacity. */
    struct block_neighbour **neigh; /* Actual neighbours (sorted on neigh[i]->b) */
};


/* ======================================================================
 * Operations
 * ====================================================================== */


/* Define a hash array size (1<<p) capable of holding a set of size 'm' */
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


/* Hash element 'k' into table of size 'm' (multiplication method) */
/* ---------------------------------------------------------------------- */
static size_t
hash_set_idx(int k, size_t m)
/* ---------------------------------------------------------------------- */
{
    double x = fmod(k * GOLDEN_RAT, 1.0);
    double y = floor(m * x);
    return y;
}


/* Insert element 'k' into set 's' of size 'm'
 * (open addressing, double probing). */
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
    h2 = 2 * hash_set_idx(k, MAX(m >> 1, 1)) + 1;

    for (i = 1; (s[j] != -1) && (s[j] != k) && (i < m); i++) {
        j += h2;
        j &= m - 1;             /* Modulo m since IS_POW2(m). */
    }

    if ((s[j] == -1) || (s[j] == k)) {
        s[j] = k;               /* Possibly no-op. */
    } else {
        j = m + 1;              /* Invalid.  Caveat emptor. */
    }

    return j;
}


/* Increase size of hash set 't' to hold 'm' elements whilst copying
 * existing elements.  This is typically fairly expensive. */
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
            assert (ret < m);
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


/* Release dynamic memory resources for hash set 't'. */
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


/* Construct an emtpy hash set capable of holding 'm' elements */
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


/* Insert element 'k' into hash set 't'. */
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
static size_t
hash_set_count_elms(const struct hash_set *set)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;

    n = 0;
    for (i = 0; i < set->m; i++) {
        n += set->s[i] != -1;
    }

    return n;
}


/* Relase dynamic memory resources for single block neighbour 'bn'. */
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


/* Construct empty block neighbour connection capable of holding
 * 'nconn' fine-scale connections (e.g., fine-scale interfaces).
 * The fine-scale table is not allocated unless nconn > 0. */
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
                new->b  = INT_MIN;
            } else {
                block_neighbour_deallocate(new);
                new     = NULL;
            }
        } else {
            new->b      = INT_MIN;
            new->fconns = NULL;
        }
    }

    return new;
}


/* Insert fine-scale connection 'fconn' into block neighbour
 * connection 'bn', but only if the bn->fconns table has been allocated.  */
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


/* Relase dynamic memory resources for single-block adjacency list 'bns'. */
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


/* Allocate a single-block adjacency list capable of holding 'nneigh'
 * connections. */
/* ---------------------------------------------------------------------- */
static struct block_neighbours *
block_neighbours_allocate(int nneigh)
/* ---------------------------------------------------------------------- */
{
    int                      i;
    struct block_neighbours *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        if (nneigh > 0) {
            new->neigh = malloc(nneigh * sizeof *new->neigh);

            if (new->neigh != NULL) {
                for (i = 0; i < nneigh; i++) { new->neigh[i] = NULL; }
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


/* Increase size of single-block adjacency list 'bns' to hold 'nneigh'
 * coarse-scale connections. */
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

        for (ret = bns->nneigh; ret < bns->cpty; ret++) {
            bns->neigh[ret] = NULL;
        }
    } else {
        ret = -1;
    }

    return ret;
}


/* Insert fine-scale connection 'fconn' into single-block adjacency
 * list 'bns' in slot corresponding to connection 'b'.
 *
 * New coarse-scale connections are assumed to hold 'expct_nconn'
 * fine-scale connections.*/
/* ---------------------------------------------------------------------- */
static int
block_neighbours_insert_neighbour(int b, int fconn, int expct_nconn,
                                  struct block_neighbours *bns)
/* ---------------------------------------------------------------------- */
{
    int i, j, p, t, nmove, ret;

    assert (bns != NULL);

    ret = 1;
    if ((bns->neigh == NULL) || (bns->cpty == 0)) {
        ret = block_neighbours_expand(1, bns);
    }

    if (ret == 1) {
        /* bns->neigh points to table containing at least one slot. */
        i = 0;
        j = bns->nneigh;

        while (i < j) {
            p = (i + j) / 2;

            assert (bns->neigh[p] != NULL);

            t = bns->neigh[p]->b;

            if      (t < b) { i = p + 1; }
            else if (t > b) { j = p + 0; }
            else            { i = j = p; }
        }

        if ((i < bns->nneigh) &&
            (bns->neigh[i] != NULL) && (bns->neigh[i]->b == b)) {
            ret = block_neighbour_insert_fconn(fconn, bns->neigh[i]);
        } else {
            if (bns->nneigh == bns->cpty) {
                assert (bns->cpty >= 1);
                ret = block_neighbours_expand(2 * bns->cpty, bns);
            }

            if (ret >= 0) {
                if (i < bns->nneigh) {
                    nmove = bns->nneigh - i;

                    memmove(bns->neigh + i + 1, bns->neigh + i + 0,
                            nmove * sizeof *bns->neigh);
                }

                bns->neigh[i] = block_neighbour_allocate(expct_nconn);

                if (bns->neigh[i] != NULL) {
                    ret = block_neighbour_insert_fconn(fconn, bns->neigh[i]);

                    bns->neigh[i]->b = b;
                    bns->nneigh += 1;
                } else {
                    ret = -1;
                }
            }
        }
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static int
count_blocks(int nc, const int *p)
/* ---------------------------------------------------------------------- */
{
    int i, max_blk;

    max_blk = -1;
    for (i = 0; i < nc; i++) {
        max_blk = MAX(max_blk, p[i]);
    }

    return max_blk + 1;
}


/* ---------------------------------------------------------------------- */
static int
derive_block_faces(int nfinef, int nblk, int expct_nconn,
                   const int *p, const int *neighbours,
                   struct block_neighbours **bns)
/* ---------------------------------------------------------------------- */
{
    int f, c1, b1, c2, b2, b_in, b_out;
    int ret;

    ret = 0;
    for (f = 0; (f < nfinef) && (0 <= ret); f++) {
        c1 = neighbours[2*f + 0];   b1 = (c1 >= 0) ? p[c1] : -1;
        c2 = neighbours[2*f + 1];   b2 = (c2 >= 0) ? p[c2] : -1;

        assert ((b1 >= 0) || (b2 >= 0));

        if ((b1 >= 0) && (b2 >= 0)) {
            b_in  = MIN(b1, b2);
            b_out = MAX(b1, b2);
        } else if (b1 >= 0) { /* (b2 == -1) */
            b_in  = b1;
            b_out = b2;
        } else {/*(b2 >= 0) *//* (b1 == -1) */
            b_in  = b2;
            b_out = b1;
        }

        if (b_in != b_out) {
            /* Block boundary */
            if (bns[b_in] == NULL) {
                bns[b_in] = block_neighbours_allocate(1);
            }

            if (bns[b_in] != NULL) {
                ret = block_neighbours_insert_neighbour(b_out, f,
                                                        expct_nconn,
                                                        bns[b_in]);
            } else {
                ret = -1;
            }
        }
    }

    if (ret >= 0) {
        ret = 0;

        for (b1 = 0; b1 < nblk; b1++) {
            if (bns[b1] != NULL) {
                ret += bns[b1]->nneigh;
            }
        }
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static size_t
coarse_topology_build_coarsef(int nblk, struct block_neighbours **bns,
                              int *neighbours, int *blkfacepos)
/* ---------------------------------------------------------------------- */
{
    int    b, n, coarse_f;
    size_t nsubf;

    coarse_f = 0;
    nsubf    = 0;

    for (b = 0; b < nblk; b++) {
        if (bns[b] != NULL) {
            for (n = 0; n < bns[b]->nneigh; n++) {
                neighbours[2*coarse_f + 0] = b;
                neighbours[2*coarse_f + 1] = bns[b]->neigh[n]->b;

                coarse_f      += 1;
                blkfacepos[b] += 1;

                if (bns[b]->neigh[n]->b >= 0) {
                    blkfacepos[bns[b]->neigh[n]->b] += 1;
                }

                if (bns[b]->neigh[n]->fconns != NULL) {
                    nsubf += hash_set_count_elms(bns[b]->neigh[n]->fconns);
                }
            }
        }
    }

    /* Derive end pointers */
    for (b = 1; b < nblk; b++) {
        blkfacepos[b] += blkfacepos[b - 1];
    }
    blkfacepos[nblk] = blkfacepos[nblk - 1];

    return nsubf;
}


/* ---------------------------------------------------------------------- */
static void
reverse_bins(int nbin, const int *pbin, int *elements)
/* ---------------------------------------------------------------------- */
{
    int b, i, j, tmp;

    for (b = 0; b < nbin; b++) {
        i = pbin[b + 0] + 0;
        j = pbin[b + 1] - 1;

        while (i < j) {
            /* Swap reverse (lower <-> upper) */
            tmp         = elements[i];
            elements[i] = elements[j];
            elements[j] = tmp;

            i += 1;             /* Increase lower bound */
            j -= 1;             /* Decrease upper bound */
        }
    }
}


/* ---------------------------------------------------------------------- */
static int
coarse_topology_build_final(int ncoarse_f, int nblk,
                            const int *neighbours,
                            int *blkfacepos, int *blkfaces,
                            struct block_neighbours **bns,
                            int *subfacepos, int *subfaces)
/* ---------------------------------------------------------------------- */
{
    int              coarse_f, b1, b2, n, subpos, subface_valid;
    size_t           i;
    struct hash_set *set;

    assert ((subfacepos == NULL) == (subfaces == NULL));

    for (coarse_f = 0; coarse_f < ncoarse_f; coarse_f++) {
        b1 = neighbours[2*coarse_f + 0];
        b2 = neighbours[2*coarse_f + 1];

        assert (b1 != b2);

        if (b1 >= 0) { blkfaces[-- blkfacepos[b1]] = coarse_f; }
        if (b2 >= 0) { blkfaces[-- blkfacepos[b2]] = coarse_f; }
    }
    assert (blkfacepos[0] == 0); /* Basic consistency */

    reverse_bins(nblk, blkfacepos, blkfaces);

    if (subfacepos != NULL) {
        coarse_f = 0;
        subpos   = 0;

        subface_valid = 1;

        for (b1 = 0; (b1 < nblk) && subface_valid; b1++) {
            for (n = 0; n < bns[b1]->nneigh; n++) {
                set = bns[b1]->neigh[n]->fconns;
                subface_valid = set != NULL;

                if (subface_valid) {
                    for (i = 0; i < set->m; i++) {
                        if (set->s[i] != -1) {
                            subfaces[subpos ++] = set->s[i];
                        }
                    }
                } else {
                    break;
                }

                subfacepos[++ coarse_f] = subpos;
            }
        }
    }

    return (subfacepos == NULL) || subface_valid;
}


/* ---------------------------------------------------------------------- */
static struct coarse_topology *
coarse_topology_build(int ncoarse_f, int nblk,
                      struct block_neighbours **bns)
/* ---------------------------------------------------------------------- */
{
    int                     subface_valid;
    size_t                  nsubf;
    struct coarse_topology *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->neighbours = malloc(2 * ncoarse_f * sizeof *new->neighbours);
        new->blkfacepos = calloc(nblk + 1      , sizeof *new->blkfacepos);

        new->blkfaces   = NULL;
        new->subfacepos = NULL;
        new->subfaces   = NULL;

        if ((new->neighbours == NULL) ||
            (new->blkfacepos == NULL)) {
            coarse_topology_destroy(new);
            new = NULL;
        } else {
            memset(new->neighbours, INT_MIN,
                   2 * ncoarse_f * sizeof *new->neighbours);

            nsubf = coarse_topology_build_coarsef(nblk, bns,
                                                  new->neighbours,
                                                  new->blkfacepos);

            if (nsubf > 0) {
                new->subfacepos = calloc(ncoarse_f + 1, sizeof *new->subfacepos);
                new->subfaces   = malloc(nsubf        * sizeof *new->subfaces);

                if ((new->subfacepos == NULL) || (new->subfaces == NULL)) {
                    free(new->subfaces);   new->subfaces   = NULL;
                    free(new->subfacepos); new->subfacepos = NULL;
                }
            }

            new->blkfaces = malloc(new->blkfacepos[nblk] * sizeof *new->blkfaces);

            if (new->blkfaces == NULL) {
                coarse_topology_destroy(new);
                new = NULL;
            } else {
                subface_valid = coarse_topology_build_final(ncoarse_f, nblk,
                                                            new->neighbours,
                                                            new->blkfacepos,
                                                            new->blkfaces,
                                                            bns,
                                                            new->subfacepos,
                                                            new->subfaces);

                if (!subface_valid) {
                    free(new->subfaces);   new->subfaces   = NULL;
                    free(new->subfacepos); new->subfacepos = NULL;
                } else {
                    new->nblocks = nblk;
                    new->nfaces  = ncoarse_f;
                }
            }
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
struct coarse_topology *
coarse_topology_create(int nc, int nf, int expct_nconn,
                       const int *p, const int *neighbours)
/* ---------------------------------------------------------------------- */
{
    int b, nblocks, ncoarse_f;

    struct block_neighbours **bns;
    struct coarse_topology   *topo;

    nblocks = count_blocks(nc, p);

    bns = malloc(nblocks * sizeof *bns);
    if (bns != NULL) {
        for (b = 0; b < nblocks; b++) {
            bns[b] = NULL;
        }

        ncoarse_f = derive_block_faces(nf, nblocks, expct_nconn,
                                       p, neighbours, bns);

        topo = coarse_topology_build(ncoarse_f, nblocks, bns);

        for (b = 0; b < nblocks; b++) {
            block_neighbours_deallocate(bns[b]);
        }

        free(bns);
    } else {
        topo = NULL;
    }

    return topo;
}


/* ---------------------------------------------------------------------- */
void
coarse_topology_destroy(struct coarse_topology *t)
/* ---------------------------------------------------------------------- */
{
    if (t != NULL) {
        free(t->subfaces);
        free(t->subfacepos);

        free(t->blkfaces);
        free(t->blkfacepos);

        free(t->neighbours);
    }

    free(t);
}
