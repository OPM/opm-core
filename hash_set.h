#ifndef HASH_SET_H_INCLUDED
#define HASH_SET_H_INCLUDED

#include <stddef.h>


/* ---------------------------------------------------------------------- */
/* Poor-man's unordered set (ind. key insert/all key extract only).       */
/* ---------------------------------------------------------------------- */
struct hash_set {
    size_t  m;                  /* Table/set capacity (1<<p for some p) */
    int    *s;                  /* Set representation */
};


/* ---------------------------------------------------------------------- */
/* Release dynamic memory resources for hash set 's'.                     */
/* ---------------------------------------------------------------------- */
void
hash_set_deallocate(struct hash_set *s);


/* ---------------------------------------------------------------------- */
/* Construct an emtpy hash set capable of holding 'm' elements            */
/* ---------------------------------------------------------------------- */
struct hash_set *
hash_set_allocate(int m);


/* ---------------------------------------------------------------------- */
/* Insert element 'k' into hash set 's'.                                  */
/* ---------------------------------------------------------------------- */
int
hash_set_insert(int k, struct hash_set *s);


/* ---------------------------------------------------------------------- */
/* Count number of valid keys in a hash set.                              */
/* ---------------------------------------------------------------------- */
size_t
hash_set_count_elms(const struct hash_set *set);


#endif  /* HASH_SET_H_INCLUDED */
