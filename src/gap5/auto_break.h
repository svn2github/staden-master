#ifndef AUTO_BREAK_H
#define AUTO_BREAK_H

#include "dstring.h"

dstring_t *auto_break_contigs(GapIO *io, int argc, contig_list_t *argv,
			      double filter_score, int filter_consensus,
			      int min_mqual, int min_score, int unique_mqual,
			      int good_score, int good_unique_score,
			      int bad_score, int bad_unique_score,
			      int large_score, int large_unique_score,
			      int spanning_score,
			      int singleton_score);

#endif /* AUTO_BREAK_H */
