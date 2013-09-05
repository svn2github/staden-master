#ifndef AUTO_BREAK_H
#define AUTO_BREAK_H

#include "dstring.h"

dstring_t *auto_break_contigs(GapIO *io, int argc, contig_list_t *argv,
			      //double filter_score, int by_consensus);
			      int min_mqual,
			      int min_score, int good_score, int bad_score,
			      int lagre_score, int spanning_score,
			      int singleton_score);

#endif /* AUTO_BREAK_H */
