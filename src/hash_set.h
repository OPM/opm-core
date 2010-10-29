/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_HASH_SET_HEADER_INCLUDED
#define OPM_HASH_SET_HEADER_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif  /* OPM_HASH_SET_HEADER_INCLUDED */
