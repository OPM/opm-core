#ifndef COARSE_CONN_H_INCLUDED
#define COARSE_CONN_H_INCLUDED

struct coarse_topology {
    int nblocks;
    int nfaces;

    int *neighbours;

    int *blkfacepos;
    int *blkfaces;

    int *subfacepos;
    int *subfaces;
};


struct coarse_topology *
coarse_topology_create(int nc, int nf, int expct_nconn,
                       const int *p, const int *neighbours);


void
coarse_topology_destroy(struct coarse_topology *t);

#endif  /* COARSE_CONN_H_INCLUDED */
