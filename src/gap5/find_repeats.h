#ifndef _FIND_REPEATS_H_
#define _FIND_REPEATS_H_
#include "io_utils.h"

int
find_repeats(GapIO *io,
             int idir,
             int minmat,
	     f_int mask,
	     float percd,
	     int num_contigs,
	     contig_list_t *clist,
	     char *outfile);

int csmatch_load_repeats(GapIO *io, FILE *fp, int match_type);

#endif /* _FIND_REPEATS_H_ */
