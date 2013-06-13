#ifndef _EDITOR_JOIN_H_
#define _EDITOR_JOIN_H_

#include "tg_gio.h"

typedef struct {
    int shift;
    int off0, off1;
    int pos0, pos1;
    int len0, len1;
    int *depad_to_pad0; // current pointer
    int *depad_to_pad1; // current pointer
    int *dp0, *dp1;     // initial malloced value
    int *res, *S;       // Myers/Miller alignment buffer and current loc in it.
    int match_len;
    int match_count;
} alignment_t;

int join_contigs(GapIO *io, tg_rec clrec, tg_rec crrec, int offset);

void alignment_free(alignment_t *a);

alignment_t *align_contigs(GapIO *io0, tg_rec c0, int pos0, int len0,
			   GapIO *io1, tg_rec c1, int pos1, int len1,  
			   int fixed_left,
			   int fixed_right);

int align_apply_edits(GapIO *io0, tg_rec c0,
		      GapIO *io1, tg_rec c1,
		      alignment_t *a);

#endif /* _EDITOR_JOIN_H_ */
