#ifndef _SHUFFLE_PADS_H
#define _SHUFFLE_PADS_H

#include <tg_gio.h>
#include "io_utils.h"
#include "io_lib/hash_table.h"

int shuffle_contigs_io(GapIO *io, int ncontigs, contig_list_t *contigs,
		       int band, int soft_clips, int flush);

int remove_pad_columns(GapIO *io, int rargc, contig_list_t *rargv,
		       int percent_pad, int quiet);

/* Original left/right soft-clips */
typedef struct {
    tg_rec rec;
    int left;
    int right;
} soft_clips;

HashTable *coherent_soft_clips(GapIO *io, tg_rec crec, int start, int end,
			       int tag_only, int min_depth, int min_tag_length,
			       Array *tag_arr);

#endif /* _SHUFFLE_PADS_H */
