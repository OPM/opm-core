#ifndef UNIQUEPOINTS_H
#define UNIQUEPOINTS_H


void finduniquepoints(const struct grdecl *g,  /* input */
		      int                 *p,  /* for each z0 in zcorn, z0 = z[p0] */
		      sparse_table_t      *z); /* list of uniq zcorn valules for each pillar*/

#endif
