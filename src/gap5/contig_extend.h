#ifndef _CONTIG_EXTEND_H_
#define _CONTIG_EXTEND_H_

#include <tg_gio.h>

int contig_extend(GapIO *io, tg_rec *contigs, int ncontigs, int min_depth,
		  int match_score, int mismatch_score);

int contig_trim(GapIO *io, tg_rec *contigs, int ncontigs, int min_depth);


/* Joint mode */
int contig_trim_and_extend(GapIO *io, tg_rec *contig, int ncontigs,
			   int do_trim, int do_extend,
			   int trim_depth, int ext_depth,
			   int ext_match_score, int ext_mismatch_score);

#endif /* _CONTIG_EXTEND_H_ */

