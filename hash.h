#ifndef HASH_H_INCLUDED
#define HASH_H_INCLUDED

#include <stddef.h>

struct hash_table {
    size_t  m;
    int    *a;
};

struct hash_table *
hash_table_allocate(int m);

void
hash_table_deallocate(struct hash_table *t);

int
hash_table_insert(int k, struct hash_table *t);

#endif  /* HASH_H_INCLUDED */
