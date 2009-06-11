#ifndef FACETOPOLOGY_H
#define FACETOPOLOGY_H


void findconnections(int n, int *pts[4],
		     int *ptnumber,
		     int *intersectionlist,
		     int *neighbors,
                     int *work,
                     sparse_table_t *ftab);

#endif
