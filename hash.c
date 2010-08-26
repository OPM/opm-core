#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "hash.h"

#define GOLDEN_RAT (0.6180339887498949) /* (sqrt(5) - 1) / 2 */
#define IS_POW2(x) (((x) & ((x) - 1)) == 0)


/* ---------------------------------------------------------------------- */
static size_t
hash_table_size(size_t m)
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
hash_table_idx(int k, size_t m)
/* ---------------------------------------------------------------------- */
{
    return floor(m * fmod(k * GOLDEN_RAT, 1.0));
}


/* ---------------------------------------------------------------------- */
static size_t
hash_table_insert_core(int k, size_t m, int *a)
/* ---------------------------------------------------------------------- */
{
    size_t h1, h2, i, j;

    assert ((0 < m) && (m < (size_t)(-1)));
    assert (IS_POW2(m));

    j = h1 = hash_table_idx(k, m);
    assert (h1 < m);

    if (a[j] == -1) { a[j] = k; }
    if (a[j] ==  k) { return j; }

    /* Double hash probing.  h2 relatively prime to 'm' */
    h2 = 2 * hash_table_idx(k, (m << 1) - 1) - 1;

    for (i = 1; (a[j] != -1) && (a[j] != k) && (i < m); i++) {
        j += h2;
        j &= m - 1;             /* Modulo m since IS_POW2(m). */
    }

    if (i < m) {
        a[j] = k;               /* Possibly no-op. */
    } else {
        j = m + 1;              /* Invalid.  Caveat emptor. */
    }

    return j;
}


/* ---------------------------------------------------------------------- */
static int
hash_table_expand(size_t m, struct hash_table *t)
/* ---------------------------------------------------------------------- */
{
    int ret, *a, *p;
    size_t i;

    assert (m > t->m);

    a = malloc(m * sizeof *a);
    if (a != NULL) {
        memset(a, -1, m * sizeof *a);

        for (i = 0; i < t->m; i++) {
            ret = hash_table_insert_core(t->a[i], m, a);
            assert (ret == t->a[i]);
        }

        p    = t->a;
        t->a = a;
        t->m = m;

        free(p);

        ret = m;
    } else {
        ret = -1;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
struct hash_table *
hash_table_allocate(int m)
/* ---------------------------------------------------------------------- */
{
    size_t             sz;
    struct hash_table *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        sz = hash_table_size(m);
        new->a = malloc(sz * sizeof *new->a);

        if (new->a == NULL) {
            hash_table_deallocate(new);

            new = NULL;
        } else {
            memset(new->a, -1, sz * sizeof *new->a);
            new->m = sz;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
void
hash_table_deallocate(struct hash_table *t)
/* ---------------------------------------------------------------------- */
{
    if (t != NULL) {
        free(t->a);
    }

    free(t);
}


/* ---------------------------------------------------------------------- */
int
hash_table_insert(int k, struct hash_table *t)
/* ---------------------------------------------------------------------- */
{
    int ret;
    size_t i;

    assert (k >= 0);
    assert (t != NULL);
    assert (IS_POW2(t->m));

    i = hash_table_insert_core(k, t->m, t->a);
    if (i == t->m + 1) {
        /* Table full.  Preferable an infrequent occurrence.  Expand
         * table and re-insert key (if possible). */
        ret = hash_table_expand(t->m << 1, t);

        if (ret > 0) {
            i = hash_table_insert_core(k, t->m, t->a);
            assert (i < t->m);

            ret = k;
        }
    } else {
        ret = k;
    }

    return ret;
}
