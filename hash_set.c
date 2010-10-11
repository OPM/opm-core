#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "hash_set.h"

/* ======================================================================
 * Macros
 * ====================================================================== */
#define GOLDEN_RAT (0.6180339887498949)              /* (sqrt(5) - 1) / 2 */
#define IS_POW2(x) (((x) & ((x) - 1)) == 0)
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))


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
            assert ((size_t) ret < m);
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
void
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
struct hash_set *
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
int
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
size_t
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
