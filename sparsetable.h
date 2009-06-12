#ifndef SPARSETABLE_H
#define SPARSETABLE_H




typedef struct{
  int m;        /* number of rows          */
  int *ptr;     /* row pointer of size m+1 */
  int position; /* first position in ptr that is not filled. */

  int n;        /* size of data            */
  void *data;   /* sparse table data       */
  

} sparse_table_t;

void             free_sparse_table     (sparse_table_t *tab);
sparse_table_t   *malloc_sparse_table  (int m, int n, int datasz);
sparse_table_t   *realloc_sparse_table (sparse_table_t *tab, int m, int n, int datasz);

#endif
