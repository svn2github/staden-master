/* Copyright Genome Research Limited (GRL). All rights reserved */

/*
 * TODO:
 * 
 * ----------------------------------------------------------------------
 * Allow sequences to move. Often we have alignments ending or starting like:
 *
 *  ACGGG
 *  AC*GGGTA
 *  AC*GGGTA
 *  ACGGGGTA
 *  AC*GGGTA
 *
 * The first sequence is reinforcing there being 4 Gs, but it actually only
 * has 3. The problem is that it cannot insert the pad as that changes the
 * sequence length.
 *
 * Solution. Let X be a specific base call (one of A, C, G, T, but always the
 * same member of that set).  X(n) is a run of 1 or more X.
 * Find cases where we have sequence*X(n) or X(n)*sequence.
 * Check the malign vector at the * to see if it also contains X. If so
 * trim * and X(n).
 *
 * ----------------------------------------------------------------------
 * Investigate 454 rate of miscall vs indel. Seems maybe we need to mirror
 * this and get the pad penalty much lower than a mismatch.
 *
 * ----------------------------------------------------------------------
 * Investigate the issue of reassigning confidence values during runs of
 * bases for 454 data. AGGGT may have confidence X 40 30 10 X if in the +ve
 * direction but X 10 30 40 X if in the -ve direction. After pad shuffling
 * we need to have the pads aligned against the low quality bases and not
 * the high quality ones. This means several things:
 *
 * 1. Reording the confidence of base-calls in a run
 * 2. Making sure the pads always end up at the same end (needs another
 *    algorithm after this one to do that).
 * 3. The pad confidence value cannot now just be the average of the two
 *    surrounding bases. Maybe the preceeding base confidence works.
 *
 * ----------------------------------------------------------------------
 * Remove the O(N^2) complexity code and make this O(N). The most obvious
 * case is inserting and deleting into the consensus. Currently this does
 * large scale memmoves over the entire contig, but in theory we can do
 * little more than local updates if we have the following:
 *
 * Consensus base structure:
 *     next/prev points
 *     counts[6]
 *     scores[6]
 *     orig_position
 *
 * Sequence fragment structure:
 *     Consensus base pointer (for left-most end)
 *     distance from last (relative offset rather than absolute)
 *     length
 *     sequence
 *
 * Then consensus pad insertion/deletion is just a matter of updating links.
 * Q: How do we handle removal of a base to which a sequence fragment is
 * pointing? I guess we need a list of fragments pointed to by the consensus
 * base (which is like the up/down pointers in ReAligner). Keeping this up to
 * date is a bit tricky.
 *
 * ----------------------------------------------------------------------
 */


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "tg_gio.h"
#include "align.h"
#include "dna_utils.h"
#include "align_lib.h"
#include "text_output.h"
#include "shuffle_pads.h"
#include "consensus.h"
#include "tg_contig.h"
#include "break_contig.h" /* contig_visible_start(), contig_visible_end() */
#include "io_lib/hash_table.h"

typedef struct {
    int pos;
    int size; /* +ve or -ve for ins/del */
} con_indel_t;

void print_malign(MALIGN *malign);
void print_moverlap(MALIGN *malign, MOVERLAP *o, int offset);

static int common_word_L(int *counts, char *seq, int len);
static int common_word_R(int *counts, char *seq, int len);
int *find_adapter(GapIO *io, int ncontigs, contig_list_t *contigs);
int rewrite_soft_clips(GapIO *io, tg_rec crec, int start, int end,
		       HashTable *h_clips);

/* Returns depadded 'pos' in s */
static int depad_clip(seq_t *s, int pos) {
    int i, p;
    for (i = p = 0; i < ABS(s->len) && i < pos; i++) {
	if (s->seq[i] != '*')
	    p++;
    }
    return p;
}

/* Repads 'pos' in s */
static int repad_clip(seq_t *s, int pos) {
    int i, p;
    for (i = p = 0; i < ABS(s->len) && p < pos; i++)
	if (s->seq[i] != '*')
	    p++;
    return i;
}


/*
 * Insert 'size' pads into a contig at position 'pos'.
 */
void malign_padcon(MALIGN *malign, int pos, int size, Array indels) {
    CONTIGL *cl = malign->contigl;
    con_indel_t *id;

    id = ARRP(con_indel_t, indels, ArrayMax(indels));
    id->pos = pos;
    id->size = size;

    for (; cl; cl = cl->next) {
	/* We do one of three things: nothing, insert, or shift */
	/* Nothing: */
	if (cl->mseg->offset+cl->mseg->length-1 < pos)
	    continue;

	/* Shift right: */
	if (cl->mseg->offset >= pos) {
	    cl->mseg->offset += size;
	    continue;
	}

	/* Insert */
	cl->mseg->length += size;
	cl->mseg->seq = (char *)realloc(cl->mseg->seq, cl->mseg->length+1);
	memmove(&cl->mseg->seq[pos - cl->mseg->offset + size],
		&cl->mseg->seq[pos - cl->mseg->offset],
		cl->mseg->length-size - (pos - cl->mseg->offset));
	memset(&cl->mseg->seq[pos - cl->mseg->offset], '*', size);
	cl->mseg->seq[cl->mseg->length] = 0;
    }

    malign_insert_scores(malign, pos, size);
}

/*
 * Returns the number of consensus pads added or -1 for error.
 * "*edited_p" is set to 0 or 1 to indicate if the sequence was edited.
 */
int edit_mseqs(MALIGN *malign, CONTIGL *cl, MOVERLAP *o, int cons_pos,
	       Array indels, int *edited_p) {
    int i, npads, poso;
    char *cp, *old_cp, *old_seq;
    int edited = 0;

    /* Cons vector */
    npads = 0;
    for (poso = i = 0; i < o->s1_len; i++) {
	if (o->S1[i] < 0) {
	    /*printf("S1:Ins %d pads at pos %d+%d=%d\n",
	      -o->S1[i], poso, cons_pos, poso+cons_pos);*/
	    malign_padcon(malign, poso+cons_pos+npads, -o->S1[i], indels);
	    npads += -o->S1[i];
	} else {
	    poso += o->S1[i];
	}
    }

    /* sequence */
    /* Trim leading pads */
    cp = o->seq2_out;
    while(*cp == '.') {
	cp++;
	cl->mseg->offset++;
    }

    //xfree(cl->mseg->seq);
    old_cp = old_seq = cl->mseg->seq;
    cl->mseg->seq = strdup(cp);
    for (cp = cl->mseg->seq; *cp; cp++) {
	if (*cp == '.')
	    *cp = '*';
	if (*old_cp) {
	    if (!edited && *old_cp != *cp)
		edited = 1;
	    old_cp++;
	}
    }
    free(old_seq);

    /* Back off trailing pads */
    while (cp > cl->mseg->seq && *(cp-1) == '*')
	cp--;

    cl->mseg->length = cp-cl->mseg->seq;

    /*
    printf("cl->mseg->seq=%.*s (len %d)\n",
	   cl->mseg->length, cl->mseg->seq, cl->mseg->length);
    */
    if (edited_p)
	*edited_p = edited;

    return npads;
}

static int CONTIGL_sort_func(const void *v1, const void *v2) {
    const CONTIGL *cl1 = *(const CONTIGL **)v1;
    const CONTIGL *cl2 = *(const CONTIGL **)v2;

    if (cl1->mseg->offset == cl2->mseg->offset) return cl1 > cl2 ? 1 : -1;

    return cl1->mseg->offset - cl2->mseg->offset;
}

/*
 * Realigning the sequences may change their start positions and hence break
 * the sorted-on-position property.
 * We make sure this is maintained here.
 */
static void resort_contigl(MALIGN *malign) {
    CONTIGL *cl, **sorted;
    int i, nele, noop = 1, last_offset = INT_MIN;

    /*
     * This list is almost sorted already, but in excessive depth areas we
     * need a decent sort algorithm. So convert to an array and sort.
     */
    for (nele = 0, cl = malign->contigl; cl; cl = cl->next, nele++) {
	if (cl->mseg->offset < last_offset)
	    noop = 0;
	last_offset = cl->mseg->offset;
    }

    if (noop)
	return;
    
    sorted = malloc(nele * sizeof(*sorted));
    if (!sorted)
	return;

    for (nele = 0, cl = malign->contigl; cl; cl = cl->next)
	sorted[nele++] = cl;

    qsort(sorted, nele, sizeof(*sorted), CONTIGL_sort_func);

    malign->contigl = sorted[0];
    for (i=0; i < nele-1; i++) {
	sorted[i]->next = sorted[i+1];
    }
    sorted[i]->next = NULL;

    free(sorted);

    last_offset = INT_MIN; noop = 1;
    for (nele = 0, cl = malign->contigl; cl; cl = cl->next, nele++) {
	if (cl->mseg->offset < last_offset)
	    noop = 0;
	last_offset = cl->mseg->offset;
    }

    return;
}

typedef struct cl_list {
    CONTIGL *cl;
    int offset;
    struct cl_list *next;
} cl_list;

/*
 * If running over a region then we may have this:
 *
 * A ------------- |                     |
 * B      ---------|------------         |
 * C               | --------------------|---
 * D               |            ---------|------------
 * E               |                     |   ------------------
 *
 * Read A and E aren't overlapping the region, therefore don't
 * get entered into MALIGN.
 * A pad at the start of B or end of D would be considered as 100%
 * pad column causing the complete removal, but read A and E may
 * not have a pad in that position.
 *
 * Solution: Remove pads only when between start..end and not
 * outside that range. Delete other pads later on.
 */
static void remove_pads(GapIO *io, MALIGN *malign, contig_t *c,
			int start, int end) {
    int i, removed = 0;
    CONTIGL *cl = malign->contigl;
    cl_list *head = NULL, *c2, *last, *next;
    int npads, depth;

    for (i = 0; i < malign->length; i++) {
	/* Add new seqs to the depth array as we meet them */
	while (cl && cl->mseg->offset == i) {
	    c2 = (cl_list *)xmalloc(sizeof(cl_list));
	    c2->next = head;
	    c2->offset = 0;
	    c2->cl = cl;
	    head = c2;
	    cl = cl->next;
	}

	/* Remove any sequences we've now passed, also counting pads */
	npads = 0;
	depth = 0;
	last = NULL;
	for (c2 = head; c2; c2 = next) {
	    next = c2->next;
	    if (c2->offset == c2->cl->mseg->length) {
		if (last)
		    last->next = c2->next;
		else
		    head = c2->next;
		xfree(c2);
		continue;
	    }
	    last = c2;
	    if (c2->cl->mseg->seq[c2->offset++] == '*')
		npads++;
	    depth++;
	}

	if (npads != depth || depth == 0)
	    continue;


	/* We have a column of pads, so remove it */
	if (i+1-removed >= start && i+1-removed <= end) {
	    //printf("Remove pad at %d\n", i+1-removed);
	    contig_delete_pad(io, &c, i+1-removed);

	    removed++;
	}
    }

    malign_recalc_scores(malign, 0, malign->length-1);
}

/*
 * Iterates through all sequences in a contig realigning them against the
 * consensus vector.
 *
 * It then adds the newly aligned sequence back into the consensus, editing the
 * sequence and tag positions/lengths too.
 * To do this we may need to shuffle the start position of sequences
 * downstream, and hence also move consensus tags.
 */
MALIGN *realign_seqs(int contig, MALIGN *malign, int band, Array indels) {
    CONTIGL *lastl = NULL, *contigl;
    int r;
    int old_start, old_end, new_start, new_end;
    int rstart, rend, rnum = 0, edited;
    MALIGN new_reg;
    int total_npads = 0;

    new_reg.nregion = 0;
    new_reg.region = NULL;

    //printf("=== Relign_seqs over %d regions\n", malign->nregion);

    rstart = malign->nregion ? malign->region[0].start : INT_MIN;
    rend   = malign->nregion ? malign->region[0].end   : INT_MAX;
    //printf("Checking reg %d: %d..%d\n", rnum, rstart, rend);

    /* FIXME
     * Keep track of n-cons-pads when checking regions.
     */

    /* Loop through all sequences in the contig */
    contigl = malign->contigl;
    while (contigl) {
	int len;
	MOVERLAP *o;
	ALIGN_PARAMS *p;
	int cons_pos;
	int npads;
#if 1
	if (contigl->mseg->offset > rend) {
	    if (++rnum >= malign->nregion) {
		//printf("Last region ended at %d\n", rend);
		break;
	    }
	    rstart = malign->region[rnum].start + total_npads;
	    rend   = malign->region[rnum].end + total_npads;
	    //printf("Checking reg %d: %d..%d\n", rnum, rstart, rend);
	}

	if (contigl->mseg->offset + contigl->mseg->length-1 < rstart) {
	    lastl = contigl;
	    contigl = contigl->next;
	    continue;
	}
#endif

	/* Obtain a depadded copy of this mseg */
	len = contigl->mseg->length;


	/* Remove sequence from malign */
	malign_remove_contigl(malign, lastl, contigl);


	/* Align sequence to malign */
	p = create_align_params();
	set_align_params (p,
			  band,
			  8, /*gap_open*/
			  8, /*gap_extend*/
			  /* EDGE_GAPS_COUNT, */
			  EDGE_GAPS_ZEROX | BEST_EDGE_TRACE,
			  RETURN_EDIT_BUFFERS | RETURN_SEQ |
			  RETURN_NEW_PADS,
			  0,  /*seq1_start*/
			  0,  /*seq2_start*/
			  0,  /*old pad sym*/
			  0,  /*new pad sym*/
			  0   /*set_job*/);

	o = create_moverlap();
	init_moverlap(o, malign, contigl->mseg->seq, malign->length, len);

	cons_pos = contigl->mseg->offset;
	o->malign_len = malign->length - cons_pos;

	/* 3 bases overhang to the right */
	if (o->malign_len > contigl->mseg->length+band/2+1)
	    o->malign_len = contigl->mseg->length+band/2+1;

	/* And 3 to the left */
	if (cons_pos > band/2+1) {
	    cons_pos -= band/2+1;
	    o->malign_len += band/2+1;
	    contigl->mseg->offset -= band/2+1;
	} else {
	    o->malign_len += cons_pos;
	    contigl->mseg->offset -= cons_pos;
	    cons_pos = 0;
	}

	{
	    char *old_cons   = malign->consensus;
	    int **old_scores = malign->scores;
	    int **old_counts = malign->counts;

	    malign->consensus += cons_pos;
	    malign->counts    += cons_pos;
	    malign->scores    += cons_pos;

	    /* fixed_malign(o, p); */
	    r = realigner_malign(o, p); /* o->score = alignment score */
	    //printf("Score = %f\n", o->score);
	    
	    /*
	    if (!r)
		print_moverlap(malign, o, cons_pos);
	    else
		puts("FAILED");
	    */

	    malign->consensus = old_cons;
	    malign->counts    = old_counts;
	    malign->scores    = old_scores;
	}

	/* Edit the sequence with the alignment */
	old_start = contigl->mseg->offset;
	old_end   = contigl->mseg->offset + contigl->mseg->length-1;
	edited = 0;
	if (r == 0 && o->S1)
	    npads = edit_mseqs(malign, contigl, o, cons_pos, indels, &edited);
	else
	    npads = 0;
	new_start = contigl->mseg->offset;
	new_end   = contigl->mseg->offset + contigl->mseg->length-1;

	/* Keep track of region adjustments as we edit the consensus */
	total_npads += npads;
	if (rend != INT_MAX)
	    rend += npads;

	/* Put sequence back */
	malign_add_contigl(malign, lastl, contigl);

	/*
	 * Check if malign->mseg has changed between removal and addition.
	 * Also count diffs here?
	 *
	 * If it's changed, call malign_add_region on a new region list
	 * so we can reduce the work load in the next pass.
	 *
	 * However, still need to keep track of diffs on sequences we're
	 * skipping?
	 */
	if (npads || edited) {
	    malign_add_region(&new_reg,
			      MIN(old_start, new_start),
			      MAX(old_end, new_end));
	}
	// TODO


	/* Update the malign structure */
	if (npads > 0) {
	    malign_recalc_scores(malign,
				 MIN(old_start, new_start),
				 MAX(old_end, new_end));
	}
	    
	/* TODO:
	 *
	 * X Realloc malign->consensus / malign->score
	 * X Move malign->consensus from here to end right by npads.
	 * X Move malign->score      " ...
	 * X Update malign->length
	 * X Recompute consensus and score over the length of this reading.
	 *
	 * If contigl was doubly linked (sorted on left and right ends
	 * separately) then we could chain left/right to only update
	 * those readings which overlap this region. For now we can
	 * just chain from left each time.  Not optimal (O(N^2) for
	 * full realignment method then) but workable perhaps.
	 *
	 * See get_malign_counts, scale_malign_scores and get_malign_consensus
	 */


	/*
	 * Check if the short-cut method gives the same result as rebuilding
	 * from scratch.
	 */
#if 0
	{
	    int i, j;
	    MALIGN *copy;
	    copy = contigl_to_malign(malign->contigl, -4, -4);

	    for (i = 0; i < copy->length; i++) {
		for (j = 0; j < copy->charset_size+2; j++) {
		    if (copy->scores[i][j] != malign->scores[i][j]) {
			printf("[%d][%d] = %d (should be %d)\n",
			       i, j,
			       malign->scores[i][j],
			       copy->scores[i][j]);
		    }
		}
	    }
	    copy->contigl = NULL;
	    destroy_malign(copy, 0);
	}
#endif

	destroy_moverlap(o);
	destroy_alignment_params(p); 

	lastl = contigl;
	contigl = contigl->next;
    }

    /* Swap regions over */
    if (0) {
	int i;
	printf("\nCur region = %d elements\n", malign->nregion);
	for (i = 0; i < malign->nregion; i++) {
	    printf("\t%d\t%d\n",
		   malign->region[i].start,
		   malign->region[i].end);
	}
	printf("\nNew region = %d elements\n", new_reg.nregion);
	for (i = 0; i < new_reg.nregion; i++) {
	    printf("\t%d\t%d\n",
		   new_reg.region[i].start,
		   new_reg.region[i].end);
	}
    }

    if (malign->region)
	free(malign->region);
    malign->region = new_reg.region;
    malign->nregion = new_reg.nregion;

    resort_contigl(malign);

    return malign;
}

/**
 * Builds and returns MALIGN from a Gap5 IO handle for the contig 'cnum'.
 */
MALIGN *build_malign(GapIO *io, tg_rec cnum, int start, int end) {
    CONTIGL *contig, *first_contig = NULL, *last_contig = NULL;
    int i, j;
    contig_iterator *citer;
    rangec_t *r;

    /* Expand start and end to the range covered by seqs overlapping
     * start .. end
     */

    {
	seq_t *s;
	citer = contig_iter_new(io, cnum, 0,
				CITER_FIRST | CITER_ICLIPPEDSTART,
				start, start);
	r = contig_iter_next(io, citer);
	if (r) {
	    s = cache_search(io, GT_Seq, r->rec);

	    start = ((s->len < 0) ^ r->comp)
		? r->end - s->right - 2
		: r->start + s->left - 2;
	}

	contig_iter_del(citer);
    }

    {
	seq_t *s;
	citer = contig_iter_new(io, cnum, 0,
				CITER_LAST | CITER_ICLIPPEDEND,
				end, end);
	r = contig_iter_next(io, citer);
	if (r) {
	    s = cache_search(io, GT_Seq, r->rec);

	    end = ((s->len < 0) ^ r->comp)
		? r->end - s->left + 2
		: r->start + s->right + 2;
	}

	contig_iter_del(citer);
    }
    
    //printf("Generating data for %d..%d\n", start, end);

    /* Generate contigl linked list */
    //citer = contig_iter_new(io, cnum, 1, CITER_FIRST, CITER_CSTART, CITER_CEND);
    citer = contig_iter_new(io, cnum, 0, CITER_FIRST, start, end);
    
    while ((r = contig_iter_next(io, citer))) {
	seq_t *s, *sorig;
	char *seq;
	int len;
	int left, l_shift, r_shift;

	assert((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ);

	contig = create_contig_link();
	contig->id = r->rec;
	contig->mseg = create_mseg();

	sorig = s = cache_search(io, GT_Seq, r->rec);
	/* Check for out-of-bounds clip points.  It shouldn't happen, but
	   gap5 databases have been seen with this problem, and we
	   don't want to crash if there are any. */
	if (s->left < 1)            s->left = 1;
	if (s->right > ABS(s->len)) s->right = ABS(s->len);

	/* Fix reads of zero length */
	if (s->right < s->left) {
	    sorig = s = cache_rw(io, s);
	    s->right = s->left;
	    if (s->right > ABS(s->len))
		s->left = s->right = ABS(s->len);
	}

	/* Copy from s->left to s->right into mseg */
	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	len = s->right - s->left + 1;
	if (NULL == (seq = malloc(len+1)))
	    return NULL;

	for (j = 0, i = s->left-1; i < s->right; i++, j++) {
	    /* Protect against the sequence containing "."; our pad sym */
	    if (s->seq[i] == '.')
		seq[j] = 'N';
	    else
		seq[j] = s->seq[i];
	}
	seq[j] = 0;

	init_mseg(contig->mseg, seq, len, r->start-1 + s->left-1);
	contig->mseg->comp  = (s != sorig);

	if (last_contig) {
	    last_contig->next = contig;
	} else {
	    first_contig = contig;
	}
	last_contig = contig;

	if (s != sorig)
	    free(s);
    }
    contig_iter_del(citer);

    /* for 454 data -6 to -10 seem to work fine */
    return contigl_to_malign(first_contig, -7, -7);
}

#define LLEN 80
struct clist {
    char *seq;
    int len;
    char line[LLEN];
};

void print_malign(MALIGN *malign) {
    int i, j;
    struct clist *depth = NULL;
    int ndepth = 0;
    CONTIGL *cl = malign->contigl;

    puts("MALIGN OUTPUT");
    for (i = 0; i < malign->length; i++) {
	/* Maintain a list of CONTIGLs covering this point */

	/* ... adding new items to the list */
	while (cl && cl->mseg->offset <= i) {
	    ndepth++;
	    /* runaway loops completely kills deskpros */
	    if (ndepth > 100000)
		abort();
	    depth = (struct clist *)realloc(depth, ndepth * sizeof(*depth));
	    depth[ndepth-1].seq = cl->mseg->seq;
	    *depth[ndepth-1].seq = tolower(*depth[ndepth-1].seq);
	    depth[ndepth-1].seq[cl->mseg->length-1] =
		tolower(depth[ndepth-1].seq[cl->mseg->length-1]);
	    depth[ndepth-1].len = cl->mseg->length;
	    memset(depth[ndepth-1].line, ' ', LLEN);
	    cl = cl->next;
	}

	for (j = 0; j < ndepth; j++) {
	    depth[j].line[i%LLEN] = (depth[j].seq) ? *depth[j].seq++ : ' ';
	    if (depth[j].len > 0 && --depth[j].len == 0) {
		depth[j].seq = NULL;
	    }
	}

	/* Print line, and remove items from depth as and when needed */
	if (i%LLEN == LLEN-1) {
	    for (j = LLEN * (int)(i/LLEN); j < i; j+=10)
		printf("%10d", j+10);
	    printf("\n");
	    for (j = 0; j < ndepth; j++) {
		if (!depth[j].seq) {
		    memmove(&depth[j], &depth[j+1],
			    (ndepth-(j+1)) * sizeof(depth[j]));
		    ndepth--;
		    j--;
		}
	    }
	    printf("\n");
	}
    }

    /* Print remainder of lines */
    if ((i-1)%LLEN != LLEN-1) {
	for (j = LLEN * (int)(i/LLEN); j < i; j+=10)
	    printf("%10d", j+10);
	printf("\n");
	for (j = 0; j < ndepth; j++) {
	    printf("%.*s\n", i - LLEN * (int)(i/LLEN), depth[j].line);
	}
	printf("\n");
    }

    free(depth);
}

void print_moverlap(MALIGN *malign, MOVERLAP *o, int offset) {
    int i, j;
    struct clist *depth = NULL;
    int ndepth = 0;
    CONTIGL *cl = malign->contigl;
    int s1op = 0, s2op = 0;
    int *S1 = o->S1;
    int *S2 = o->S2;
    char *seq = o->seq2;
    int cins = 0;

    for (i = offset; i < malign->length+offset; i++) {
	/* Maintain a list of CONTIGLs covering this point */

	/* ... adding new items to the list */
	for (; cl && cl->mseg->offset+cins <= i; cl = cl->next) {
	    if (cl->mseg->offset+cins + cl->mseg->length-1 < i)
		continue;
	    ndepth++;
	    /* runaway loops completely kills deskpros */
	    if (ndepth > 1000)
		abort();
	    depth = (struct clist *)realloc(depth, ndepth * sizeof(*depth));
	    depth[ndepth-1].seq = cl->mseg->seq + i-(cl->mseg->offset+cins);
	    depth[ndepth-1].len = cl->mseg->length - (i-(cl->mseg->offset+cins));
	    memset(depth[ndepth-1].line, ' ', LLEN);
	}

	if (!s1op) {
	    s1op = *S1++;
	    if (S1-o->S1 > o->s1_len)
		break;
	}
	if (!s2op) {
	    s2op = *S2++;
	    if (S2-o->S2 > o->s2_len)
		break;
	}

	printf("%4d: ", i);

	if (s1op < 0) {
	    /* Ins to consensus */
	    s1op++;
	    printf("%c\n", *seq++);
	    cins++;
	    continue;
	} else if (s2op > 0) {
	    /* Match/mismatch */
	    printf("%c ", *seq++);
	    s2op--;
	} else if (s2op < 0) {
	    /* Ins to sequence */
	    printf("  ");
	    s2op++;
	}

	s1op--;
	for (j = 0; j < ndepth; j++) {
	    printf("%c", *depth[j].seq++);
	    if (--depth[j].len == 0) {
		depth[j].seq = NULL;
		memmove(&depth[j], &depth[j+1],
			(ndepth-(j+1)) * sizeof(depth[j]));
		ndepth--;
		j--;
	    }
	}
	printf("\n");
    }

    free(depth);
}

#include <ctype.h>
int64_t malign_diffs(MALIGN *malign, int64_t *tot) {
    CONTIGL *cl;
    int64_t diff_count = 0, tot_count = 0;

    for (cl = malign->contigl; cl; cl = cl->next) {
	int i;

	/*
	for (i = 0; i < cl->mseg->length; i++, end_gaps++) {
	    if (cl->mseg->seq[i] != '*')
		break;
	}
	for (i = cl->mseg->length-1; i >= 0; i--, end_gaps++) {
	    if (cl->mseg->seq[i] != '*')
		break;
	}
	*/

#if 0
	for (i = 0; i < cl->mseg->length; i++) {
	    char c = toupper(malign->consensus[i+cl->mseg->offset]);
	    char s = toupper(cl->mseg->seq[i]);
	    if (c == '-')
		c = '*';

	    /*printf("%c", c==s ? '.' : s);*/
	    if (s != c)
		diff_count++;
	    tot_count++;
	}
#else
	/* See set_malign_lookup() */
	static int l[256] = {
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, /*   0-15 */
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, /*  16 */
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 4, 5, 5, /*  32 */
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, /*  48 */
	    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5, 5, /*  64 */
	    5, 5, 5, 5, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, /*  80 */
	    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5, 5, /*  96 */
	    5, 5, 5, 5, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, /* 112-127 */
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, /* 128 */
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5  /* 255 */
	};
	for (i = 0; i < cl->mseg->length; i++) {
	    unsigned char s = l[(uint8_t) cl->mseg->seq[i]];

	    /*printf("%c", c==s ? '.' : s);*/
	    diff_count += malign->scores[i+cl->mseg->offset][s];
	}
	tot_count  += 128 * cl->mseg->length;
#endif
    }

    if (tot)
	*tot = tot_count;
    return diff_count;
}

#if 0
static void update_consensus_tags(GapIO *io, int cnum, MALIGN *malign) {
    int i, last = 0;
    for (i = 0; i < malign->length; i++) {
	int p = malign->orig_pos[i];
	if (p == 0) {
	    /* Insertion */
	    shift_contig_tags(io, cnum, i+1, +1);
	} else {
	    if (p-last != 1) {
		/* Deletion */
		shift_contig_tags(io, cnum, i+1, 1-(p-last));
	    }
	    last = p;
	}
    }
}
#endif

/*
 * Moves tags on sequence 'srec' in contig 'crec' right by one if they
 * start at or beyond pos, and extends by one if they overlap pos.
 *
 * For optimisation purposes, we already know the sequence spans positions
 * start..end, so we use this for filtering our tag search.
 * However there is much else that can be optimised. We repeatedly query
 * contigs, bins, etc and almost certainly perform many contig iters over
 * the same region. (Functionality first, efficiency later.)
 */
static void tag_shift_for_insert(GapIO *io, tg_rec crec, tg_rec srec,
				 int start, int len, int pos, tg_rec brec,
				 int dist) {
    contig_iterator *ci;
    rangec_t *r;
    contig_t *c = cache_search(io, GT_Contig, crec);;
    int end = start + len-1;

    //printf("> tag in seq %"PRIrec" at %d+%d\n", srec, start, pos);

    cache_incr(io, c);

    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST | CITER_ISTART,
				 start+pos, end, GRANGE_FLAG_ISANNO);
    if (!ci) {
	cache_decr(io, c);
	return;
    }

    while ((r = contig_iter_next(io, ci))) {
	range_t r2, *r_out;
	anno_ele_t *a;
	bin_index_t *bin;

	if (r->pair_rec != srec)
	    continue;

	bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	r2.start    = (r->start >= start+pos) ? r->start+dist : r->start;
	r2.end      = r->end+dist;
	r2.mqual    = r->mqual;
	r2.rec      = r->rec;
	r2.pair_rec = r->pair_rec;
	r2.flags    = r->flags;
	bin = bin_add_to_range(io, &c, brec, &r2, &r_out, NULL, 0);

	a = cache_search(io, GT_AnnoEle, r->rec);
	if (a->bin != bin->rec /*||
	    a->bin_idx != r_out - ArrayBase(range_t, bin->rng)*/) {
	    /* Annotation moved bins */
	    a = cache_rw(io, a);
	    a->bin = bin->rec;
	    //a->bin_idx = r_out - ArrayBase(range_t, bin->rng);
	}
    }

    contig_iter_del(ci);
    cache_decr(io, c);
}

static void tag_shift_for_delete(GapIO *io, tg_rec crec, tg_rec srec,
				 int start, int len, int pos, tg_rec brec,
				 int dist) {
    contig_iterator *ci;
    rangec_t *r;
    contig_t *c = cache_search(io, GT_Contig, crec);;
    int end = start + len-1;

    //printf("< tag in seq %"PRIrec" at %d\n", srec, pos);

    cache_incr(io, c);

    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST | CITER_ISTART,
				 start+pos, end, GRANGE_FLAG_ISANNO);
    if (!ci) {
	cache_decr(io, c);
	return;
    }

    while ((r = contig_iter_next(io, ci))) {
	range_t r2, *r_out;
	anno_ele_t *a;
	bin_index_t *bin;

	if (r->pair_rec != srec)
	    continue;

	bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	r2.start    = (r->start > start+pos) ? r->start-dist : r->start;
	r2.end      = r->end-dist;
	r2.mqual    = r->mqual;
	r2.rec      = r->rec;
	r2.pair_rec = r->pair_rec;
	r2.flags    = r->flags;

	if (r2.end < r2.start) {
	    /* Tag entirely removed now, it must have been on a pad */
	    a = cache_search(io, GT_AnnoEle, r->rec);
	    a = cache_rw(io, a);
	    cache_deallocate(io, a);
	    continue;
	}
	bin = bin_add_to_range(io, &c, brec, &r2, &r_out, NULL, 0);

	a = cache_search(io, GT_AnnoEle, r->rec);
	if (a->bin != bin->rec /*||
	    a->idx != r_out - ArrayBase(range_t, bin->rng)*/) {
	    /* Annotation moved bins */
	    a = cache_rw(io, a);
	    a->bin = bin->rec;
	    //a->bin_idx = r_out - ArrayBase(range_t, bin->rng);
	}
    }

    cache_decr(io, c);
    contig_iter_del(ci);
}

/*
 * Takes a multiple alignment and updates the on-disk data structures to
 * match. This needs to correct confidence values, original positions and
 * tags too.
 */
void update_io(GapIO *io, tg_rec cnum, MALIGN *malign, Array indels) {
    CONTIGL *cl;
    tg_rec rnum;
    range_t r, *r_out;
    bin_index_t *bin;
    contig_t *c = cache_search(io, GT_Contig, cnum);
    size_t i, nindel;

    cache_incr(io, c);

    /*
     * To minimise number of data modifications we use a three step approach.
     *
     * Step 1: insert columns of pads, shifting reads as appropriate.
     * Step 2: edit sequence alignments as required, possibly involving
     *         moving sequences and/or adding and removing pads.
     * Step 3: remove columns of entire pads.
     *
     * This means that when we introduce a column of pads we don't have
     * to make edits to every single read position down stream, and can
     * instead make use of the optimised recursive bin functions to do this
     * for us.
     */

    /* Step 1: make indels */
    nindel = ArrayMax(indels);
    for (i = 0; i < nindel; i++) {
	con_indel_t *id = arrp(con_indel_t, indels, i);
	int j;

	if (id->size > 0) {
	    contig_insert_bases(io, &c, id->pos+1, '*', -1, id->size);
	} else {
	    for (j = 0; j < -id->size; j++) {
		contig_delete_pad(io, &c, id->pos+1);
	    }
	}
    }

    /* Step 2: edit alignments */
    for (cl = malign->contigl; cl; cl = cl->next) {
	seq_t *s, *sorig;
	int len, update_range = 0;
	int shift, orig_start;

	rnum = cl->id;
	
	sorig = cache_search(io, GT_Seq, rnum);
	cache_incr(io, sorig);
	s = dup_seq(sorig);
	if (cl->mseg->comp)
	    complement_seq_t(s);

	len = s->right - s->left + 1;

	sequence_get_position(io, s->rec, NULL, &orig_start, NULL, NULL);

	/* Check if sequence has changed. If so assign a new one */
	if (cl->mseg->length != len ||
	    memcmp(s->seq + s->left-1, cl->mseg->seq, cl->mseg->length) != 0) {
	    int newlen = s->left-1 + ABS(s->len) - s->right + cl->mseg->length;
	    int i, j, np;
	    char   *newseq  = malloc(newlen+1);
	    int8_t *newconf = malloc(newlen+1);

	    /* Build new seq/conf arrays */
	    memcpy(newseq,  s->seq,  s->left-1);
	    memcpy(newconf, s->conf, s->left-1);

	    memcpy(&newseq[s->left-1], cl->mseg->seq, cl->mseg->length);

	    /*
	     * Step through both old and new sequences working out how
	     * they differ. This will (*should*) be entire pad movements.
	     * i = index to old seq
	     * j = index to new seq
	     * np = number of pads added minus removed from old seq.
	     */
	    np = 0;
	    for (i =j =s->left-1;
		 i < ABS(s->len) && j < s->left-1 + cl->mseg->length;
		 ) {
		/* Bases match */
		if (toupper(newseq[j]) == toupper(s->seq[i]) ||
		    (s->seq[i] == '.' && newseq[j] == 'N')) {
		    if (isupper(s->seq[i]))
			newseq[j] = toupper(newseq[j]);
		    else
			newseq[j] = tolower(newseq[j]);
		    newconf[j] = s->conf[i];
		    i++, j++;
		    continue;
		}

		/* Pad removed */
		if (s->seq[i] == '*') {
		    i++;
		    tag_shift_for_delete(io, cnum, rnum, cl->mseg->offset,
					 s->right - s->left + 1,
					 i+np-- - (s->left-1),
					 s->bin, 1);
		    /*
		    if (io_length(io, rnum) < 0) {
			tag_shift_for_delete(io, rnum, r.length - i + 1);
		    } else {
			tag_shift_for_delete(io, rnum, i+np--);
		    }
		    */
		    continue;
		}

		/* Pad created */
		if (newseq[j] == '*') {
		    int k;
		    int ql = 0, qr = 0;
		    for (k = i-1; k >= 0; k--) {
			if (s->seq[k] != '*') {
			    ql = s->conf[k];
			    break;
			}
		    }
		    for (k = i+1; k < s->right; k++) {
			if (s->seq[k] != '*') {
			    qr = s->conf[k];
			    break;
			}
		    }
		    newconf[j] = MIN(ql, qr); /* min conf of neighbours */
		    j++;
		    tag_shift_for_insert(io, cnum, rnum, cl->mseg->offset,
					 cl->mseg->length,
					 i+ ++np - (s->left-1),
					 s->bin, 1);
		    /*
		    if (io_length(io, rnum) < 0) {
			tag_shift_for_insert(io, rnum, r.length - i + 1);
		    } else {
			tag_shift_for_insert(io, rnum, i+ ++np);
		    }
		    */
		    continue;
		}

		fprintf(stderr, "Alignment introduced non-pad character");
		abort();
	    }

	    /* Pads previously at the end of the reading & now removed */
	    while (i < s->right) {
		if (s->seq[i] == '*') {
		    i++;
		    tag_shift_for_delete(io, cnum, rnum, cl->mseg->offset,
					 s->right - s->left + 1,
					 i+np-- - (s->left-1),
					 s->bin, 1);
		    /*
		    if (io_length(io, rnum) < 0) {
			tag_shift_for_delete(io, rnum, r.length - i + 1);
		    } else {
			tag_shift_for_delete(io, rnum, i+np--);
		    }
		    */
		} else {
		    /* Error: clipped data that wasn't a pad */
		    abort();
		}
	    }

	    /* Should only be pads remaining in newseq, if anything */
	    s->right = j;
	    for (; j < s->left-1 + cl->mseg->length; j++) {
		if (newseq[j] != '*') {
		    fprintf(stderr, "Alignment introduced non-pad character");
		    abort();
		}
		newconf[j] = 0;
	    }

	    /* Append on the right hand cutoff data */
	    for (; i < ABS(s->len); i++, j++) {
		newseq[j]  = s->seq[i];
		newconf[j] = s->conf[i];
	    }
	    if (j != newlen) {
		abort();
	    }

	    /* Write it back out */
	    /* Copy newseq/newconf into seq_t */

	    s->seq = newseq;
	    s->conf = newconf;
	    update_range = 0;
	    if (ABS(s->len) != j) {
		/* Length change implies updating the range array too */
		s->len = s->len >= 0 ? j : -j;
		update_range = 1;
	    }

	    if (cl->mseg->comp)
		complement_seq_t(s);

	    /* The memcpy trashes the block pointer, so special care needed */
	    {
		sorig = cache_rw(io, sorig);
		void *blk = sorig->block;
		memcpy(sorig, s, sizeof(seq_t)); 
		sorig->block = blk;
	    }

	    if (update_range)
		sorig = cache_item_resize(sorig, sizeof(*sorig) +
					  sequence_extra_len(sorig));

	    sequence_reset_ptr(sorig);

	    if (s->name)
		memcpy(sorig->name,       s->name,       s->name_len+1);
	    if (s->trace_name)
		memcpy(sorig->trace_name, s->trace_name, s->trace_name_len+1);
	    if (s->alignment)
		memcpy(sorig->alignment,  s->alignment,  s->alignment_len+1);
	    memcpy(sorig->seq,  s->seq,  ABS(s->len));
	    memcpy(sorig->conf, s->conf, ABS(s->len));

	    xfree(newconf);
	    xfree(newseq);
	}

	{
	    int st, en, or;
	    sequence_get_position(io, s->rec, NULL, &st, &en, &or);
	    if (or ^ (sorig->len < 0)) {
		shift = ABS(sorig->len) - sorig->right;
	    } else {
		shift = sorig->left-1;
	    }
	    st += shift;
	    if (st != cl->mseg->offset+1) {
		update_range = 1;
	    }
	}

	free(s);

	if (update_range) {
	    int bin_changed = 0;
	    int dist;

	    /* Get old range and pair data */
	    s = sorig;
	    bin = cache_search(io, GT_Bin, s->bin);
	    r = *arrp(range_t, bin->rng, s->bin_index);
	    assert(r.rec == s->rec);

	    dist =  cl->mseg->offset + 1 - shift - orig_start;
	    if (dist > 0) {
		tag_shift_for_insert(io, cnum, rnum,
				     orig_start, ABS(s->len)+dist,
				     0, s->bin, dist);
	    } else if (dist < 0) {
		tag_shift_for_delete(io, cnum, rnum,
				     orig_start+dist, ABS(s->len)-dist,
				     0, s->bin, -dist);
	    }

	    /* Update range, tedious and slow way */
	    bin_remove_item(io, &c, GT_Seq, s->rec);

	    r.start = cl->mseg->offset + 1 - shift;
	    r.end   = r.start + ABS(s->len) - 1;
	    bin = bin_add_range(io, &c, &r, &r_out, NULL, 0);

	    /* Check if the new bin has a different complemented status too */
	    if (s->bin != bin->rec) {
		int old_comp = bin_get_orient(io, s->bin);
		int new_comp = bin_get_orient(io, bin->rec);

		if (new_comp != old_comp) {
		    //int tmp;
		    s = cache_rw(io, s);
		    s->len *= -1;
		    s->flags ^= SEQ_COMPLEMENTED;
		    //tmp = s->left;
		    //s->left  = ABS(s->len) - (s->right-1);
		    //s->right = ABS(s->len) - (tmp-1);
		}

		bin_changed = 1;
	    }
	
	    /* Update seq bin & bin_index fields */
	    s = cache_rw(io, s);
	    s->bin = bin->rec;
	    s->bin_index = r_out - ArrayBase(range_t, bin->rng);

	    if (bin_changed) {
		if (-1 == sequence_fix_anno_bins(io, &s)) {
		    verror(ERR_WARN, "update_io",
			   "sequence_fix_anno_bins() failure");
		}
	    }
	}

	cache_decr(io, sorig);
    }

    /* Step 3 (remove pad columns) done in calling function. */

    cache_decr(io, c);
}

#if 0
static int isort(const void *vp1, const void *vp2) {
    return *(const int *)vp2 - *(const int *)vp1;
}

/*
 * Specifically for 454 data this reassigns confidence values to bases in
 * a run of the same base type.
 * It also reassigns confidence values of pads to be the minimum confidence
 * of the surrounding base call.
 */
void reassign_confidence_values(GapIO *io, int cnum) {
    GContigs c;
    GReadings r;
    int rnum;
    int scores[1000]; /* FIXME: check if we overflow! */

    contig_read(io, cnum, c);
    for (rnum = c.left; rnum; rnum = r.right) {
	char last = 0;
	char *seq;
	int1 *conf;
	int i, j, k;
	int cl, cr;

	gel_read(io, rnum, r);
	seq = TextAllocRead(io, r.sequence);
	conf = DataAllocRead(io, r.confidence, 1);

	/* Rearrange confidence in runs of bases */
	for (i = 0; i < r.length; i++) {
	    /* Find first non-pad, at 'i' */
	    while (i < r.length && seq[i] == '*')
		i++;
	    k = 0;
	    scores[k++] = conf[i];
	    last = seq[i];

	    /* Count how many there are. First diff base at 'j' */
	    j = i+1;
	    while (j < r.length && (seq[j] == '*' || seq[j] == last)) {
		if (seq[j] != '*')
		    scores[k++] = conf[j];
		j++;
	    }
		   
	    if (k != 1) {
		/* We have a run of k items (from >='i' and <'j') */
		qsort(scores, k, sizeof(*scores), isort);
		
		/* Reassign */
		j = i; k = 0;
		while (j < r.length && (seq[j] == '*' || seq[j] == last)) {
		    if (seq[j] != '*')
			conf[j] = scores[k++];
		    j++;
		}
	    }

	    i = j-1;
	}

	/* Reassign confidences to pads */
	cl = 0;
	for (i = 0; i < r.length; i++) {
	    if (seq[i] == '*') {
		for (j = i+1; j < r.length && seq[j] == '*'; j++)
		    ;
		cr = j < r.length ? conf[j] : 0;
		/* conf[i] = MIN(cl, cr); */
		conf[i] = (cl+cr)/2;
	    } else {
		cl = conf[i];
	    }
	}

	DataWrite(io, r.confidence, conf, r.length, 1);
	xfree(seq);
	xfree(conf);
    }
}
#endif

/*
 * The start..end range constitutes a previous range of Ns in the
 * consensus, probably due to artificially scaffolding contigs by
 * butting them end to end with a few Ns inbetween.
 *
 * The realigner may correctly identify the overlap, or it may just
 * give us garbage. We check here whether the new overlap (assuming
 * there is one) is valid and if not we slip the reads relative to one
 * another to reintroduce the gap.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int validate_N(GapIO *io, HashTable *h_clips, tg_rec trec,
		      tg_rec crec, int start, int end) {
    contig_iterator *citer;
    rangec_t *r;
    consensus_t *cons;
    int failed_reads = 0, all_reads = 0;
    int x_start = start - 50;
    int x_end = end + 50;
    int mid = (start+end)/2;
    int left_most = start;
    int right_most = end;
    int i, snp, tot;

    if (!(cons = calloc(x_end-x_start+1, sizeof(*cons))))
	return -1;
    calculate_consensus(io, crec, x_start, x_end, cons);

    citer = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);
    while ((r = contig_iter_next(io, citer))) {
	int i, diff_l, diff_r, i_start;
	seq_t *s, *sorig;
	int ignore = 0;

	// Limit this check to only reads that we have extended
	if (!HashTableSearch(h_clips, (char *)&r->rec, sizeof(tg_rec)))
	    continue;

	s = sorig = cache_search(io, GT_Seq, r->rec);

	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	//printf("#%"PRIrec"\n", r->rec);

	// Left side
	i_start = MAX(s->left-1, x_start - r->start);
	for (diff_l = 0, i = i_start; i < s->right; i++) {
	    if (r->start + i > mid)
		break;

	    //printf("  Chk %c %c\n", s->seq[i], cons[r->start+i - x_start]);
	    if (s->seq[i] != "ACGT*N"[cons[r->start+i - x_start].call])
		diff_l++;
	}
	diff_l = (diff_l > .3*(i-i_start)); // boolean, either good or bad
	//printf("  L %d bases, %d diff\n", i-i_start, diff_l);
	
	if (r->start+s->right - mid > mid - (r->start+s->left-1))
	    if (left_most > r->start + s->left-1)
		left_most = r->start + s->left-1;

	if (i-i_start <= 0)
	    ignore = 1;

	// Right side
	for (diff_r = 0; i < s->right; i++) {
	    if (r->start + i > x_end)
		break;
	    //printf("  Chk %c %c\n", s->seq[i], cons[r->start+i - x_start]);
	    if (s->seq[i] != "ACGT*N"[cons[r->start+i - x_start].call])
		diff_r++;
	}
	diff_r = (diff_r > .3*(r->start+i - mid));
	//printf("  R %d bases, %d diff\n", r->start+i - mid, diff_r);

	if (r->start+s->right - mid < mid - (r->start+s->left-1))
	    if (right_most < r->start + s->right)
		right_most = r->start + s->right;

	if (r->start+i - mid <= 0)
	    ignore = 1;

	if (ignore)
	    continue;

	if (diff_l || diff_r)
	    failed_reads++;
	all_reads++;

	if (s != sorig)
	    free(s);
    }
    contig_iter_del(citer);

    tot = snp = 0;
    for (i = MAX(left_most,mid-20); i <= MIN(right_most,mid+20); i++) {
	if (i < x_start)
	    continue;
	if (i > x_end)
	    break;
	if (cons[i-x_start].scores[6]>0)
	    snp++;
	tot++;
    }


    printf("For %d..%d => %d of %d poor, snp=%d of %d\n", left_most, right_most, failed_reads, all_reads, snp, tot);

    /* If more than 1/4tr are poor then we shift, otherwise it's ok */
    if (failed_reads*3 <= all_reads && snp*4 < tot) {
	contig_t *c = cache_search(io, GT_Contig, crec);
	bin_remove_item(io, &c, GT_AnnoEle, trec);
	vmessage("Consensus N in =%"PRIrec" at %d..%d resolved/ignored\n",
		 crec, start, end);
    } else {
	char *n_cons;

	printf("Leftmost unclipped = %d, rightmost = %d\n",
	       left_most, right_most);

	// Unclip and see if it reintroduces the N
	citer = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);
	while ((r = contig_iter_next(io, citer))) {
	    seq_t *s;
	    HashItem *hi;
	    soft_clips *sc;

	    hi = HashTableSearch(h_clips, (char *)&r->rec, sizeof(tg_rec));
	    if (!hi) continue;
	    sc = (soft_clips *)hi->data.p;

	    s = cache_search(io, GT_Seq, r->rec);
	    s = cache_rw(io, s);
	    s->left  = repad_clip(s, sc->left);
	    s->right = repad_clip(s, sc->right);
	}
	contig_iter_del(citer);

	n_cons = malloc(x_end - x_start + 1);
	calculate_consensus_simple(io, crec, x_start, x_end, n_cons, NULL);
	printf("%d..%d %.*s\n", x_start, x_end, x_end - x_start+1, n_cons);
	for (i = 0; i <= x_end - x_start + 1; i++) {
	    if (n_cons[i] == 'N')
		break;
	}
	free(n_cons);

	if (i <= x_end - x_start + 1) {
	    // has N
	    contig_t *c = cache_search(io, GT_Contig, crec);
	    bin_remove_item(io, &c, GT_AnnoEle, trec);

	    vmessage("Consensus N in =%"PRIrec" at %d..%d "
		     "resolved by unclipping.\n",
		     crec, start, end);
	} else {
	    vmessage("Consensus N in =%"PRIrec" at %d..%d retained\n",
		     crec, start, end);
	}
    }


    free(cons);

    return 0;
}

/*
 * Checks soft-clipping by counting SNPs over a region and comparing
 * to the original unclipped version (TODO). If it is significantly
 * more then we conclude the soft-clips are masking a misassembly so
 * we keep the tag. Otherwise we remove the tag and consider the
 * problem as resolved.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int validate_clip(GapIO *io, HashTable *h_clips, tg_rec trec,
			 tg_rec crec, int start, int end) {
    consensus_t *cons;
    int i, snp = 0, expected_snp = 0;
    anno_ele_t *e = cache_search(io, GT_AnnoEle, trec);
    contig_iterator *citer;
    rangec_t *r;
    int mismatch = 0, tot = 0;

    /* Compute pair-wise consensus and look for high SNP rates */
    if (!(cons = calloc(end-start+1, sizeof(*cons))))
	return -1;

    if (e->comment)
	sscanf(e->comment, "SNPs=%d\n", &expected_snp);

    calculate_consensus(io, crec, start, end, cons);
    for (i = start; i <= end; i++) {
	if (cons[i-start].scores[6]>0)
	    snp++;
    }

    if ((snp-expected_snp) >= 0.3 * (end - start + 1)) {
	vmessage("Validation of coherent soft-clip =%"PRIrec
		 " at %d..%d failed\n", crec, start, end);
	free(cons);
	return 0;
    }

    // Passed easy validation, now look harder incase it's deep
    // and SNPs wouldn't be called.
    citer = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);
    while ((r = contig_iter_next(io, citer))) {
	HashItem *hi;
	soft_clips *c;
	int left, right;
	seq_t *s;

	if (!(hi = HashTableSearch(h_clips, (char *)&r->rec, sizeof(r->rec))))
	    continue;

	c = (soft_clips *)hi->data.p;
	s = cache_search(io, GT_Seq, r->rec);
	    
	left  = repad_clip(s, c->left);
	right = repad_clip(s, c->right);

	// Right end
	if (right != s->right) {
	    if ((s->len<0) ^ r->comp) {
		for (i = right; i < s->right; i++, tot++) {
		    if (r->end -i -1 >= end)
			continue;
		    if (r->end -i -1 < start)
			break;
		    if (toupper(complement_base(s->seq[i])) !=
			"ACGTN"[cons[r->end - i - start].call])
			mismatch++;
		}
	    } else {
		for (i = right; i < s->right; i++, tot++) {
		    if (r->start + i < start)
			continue;
		    if (r->start + i >= end)
			break;
		    if (toupper(s->seq[i]) !=
			"ACGTN"[cons[r->start + i - start].call])
			mismatch++;
		}
	    }
	}


	// Left end
	if (left != s->left) {
	    if ((s->len<0) ^ r->comp) {
		for (i = s->left-1; i <= left; i++, tot++) {
		    if (r->end -i -1 >= end)
			continue;
		    if (r->end -i -1 < start)
			break;
		    if (toupper(complement_base(s->seq[i])) !=
			"ACGTN"[cons[r->end - i - start].call])
			mismatch++;
		}
	    } else {
		for (i = s->left-1; i <= left; i++, tot++) {
		    if (r->start + i < start)
			continue;
		    if (r->start + i >= end)
			break;
		    if (toupper(s->seq[i]) !=
			"ACGTN"[cons[r->start + i - start].call])
			mismatch++;
		}
	    }
	}
    }

    contig_iter_del(citer);
    
    if (3*mismatch < tot) {
	contig_t *c = cache_search(io, GT_Contig, crec);
	
	bin_remove_item(io, &c, GT_AnnoEle, trec);
	vmessage("Validation of coherent soft-clip =%"PRIrec
		 " at %d..%d passed\n", crec, start, end);
    } else {
	vmessage("Validation of coherent soft-clip =%"PRIrec
		 " at %d..%d failed\n", crec, start, end);
    }

    free(cons);
    return 0;
}

/*
 * Validates tagged soft-clip regions to verify that the extension looks
 * valid. If it has a large degree of discrepancy and relatively few
 * spanning reads, then we assume it was two contigs in a scaffold
 * butted up end to end, so we slip one past the other to ensure we
 * still have no overlap.
 */
static void validate_clip_regions(GapIO *io, HashTable *h_clips,
				  Array tag_arr) {
    int i;

    for (i = 0; i < ArrayMax(tag_arr); i++) {
	tg_rec trec = arr(tg_rec, tag_arr, i), crec;
	int start, end;
	anno_ele_t *e = cache_search(io, GT_AnnoEle, trec);
	
	if (!e)
	    continue;

	anno_get_position(io, trec, &crec, &start, &end, NULL);
	if (e->tag_type == str2type("NCLP"))
	    validate_N(io, h_clips, trec, crec, start, end);
	else
	    validate_clip(io, h_clips, trec, crec, start, end);
    }
}


#define CHUNK_SIZE 32768

int shuffle_contigs_io(GapIO *io, int ncontigs, contig_list_t *contigs,
		       int band, int soft_clips, int flush) {
    int i; //, start;
    Array indels;
    int *counts = NULL;
    int c_shift = 0;
    
    set_malign_lookup(5);
    /* set_alignment_matrix("/tmp/nuc_matrix", "ACGTURYMWSKDHVB-*"); */

    indels = ArrayCreate(sizeof(con_indel_t), 0);

    if (soft_clips)
	counts = find_adapter(io, ncontigs, contigs);

    for (i = 0; i < ncontigs; i++) {
	tg_rec cnum = contigs[i].contig;
	int j;
	int64_t old_score, new_score, tot_score, orig_score;
	MALIGN *malign;
	contig_t *c;
	HashTable *h_clips = NULL;
	Array tag_arr;
	int sub_start, sub_end;

	if (DB_VERS(io) >= 5) {
	    c = cache_search(io, GT_Contig, cnum);
	    if (!c)
		continue;
	    if (c->nseqs / (c->end - c->start + 1) >= 100) {
		verror(ERR_WARN, "shuffle_contigs_io",
		       "Skipping contig %s due to excessive depth\n",
		       get_contig_name(io, cnum));
		continue;
	    }
	}

	/*
	 * The shuffle pads code (malign) comes from gap4 and has lots of
	 * assumptions that the contig goes from base 1 to base N.
	 * Fixing these assumptions is a lot of work, so for now we will take
	 * the cheat route of moving the contig to ensure the assumption
	 * is valid.
	 */
	c = cache_search(io, GT_Contig, cnum);
	c_shift = 1-c->start;

	if (c_shift != 0) {
	    if (move_contig(io, cnum, c_shift) != 0)
		return -1;
	}

	/*
	 * Iterate over 33k chunks from end going backwards. Decrement
	 * our sub-range by 32k each time. This reduces the maximum
	 * memory capacity and also prevents our multi-pass method
	 * from purging the in-memory cache on long contigs.
	 */
	sub_end = contigs[i].end;
	sub_start = sub_end - (CHUNK_SIZE+200);

	do {
	    contig_list_t cl;

	    if (sub_start < contigs[i].start)
		sub_start = contigs[i].start;

	    cl.contig = contigs[i].contig;
	    cl.start  = sub_start + c_shift;
	    cl.end    = sub_end   + c_shift;

	    vmessage("Shuffling pads for contig %s %d..%d\n",
		     get_contig_name(io, cnum),
		     sub_start + c_shift, sub_end + c_shift);

	    if (soft_clips)
		h_clips = coherent_soft_clips(io,
					      cl.contig,
					      cl.start,
					      cl.end,
					      counts,
					      0, 3, 5,
					      &tag_arr);

	    //printf("Shuffle #%"PRIrec" from %d..%d, shift %d\n",
	    //       contigs[i].contig, contigs[i].start, contigs[i].end, c_shift);

	    malign = build_malign(io,
				  cl.contig,
				  cl.start,
				  cl.end);
	    resort_contigl(malign);

	    malign_add_region(malign,
			      cl.start,
			      cl.end);

	    ArrayMax(indels) = 0;
	    orig_score = new_score = malign_diffs(malign, &tot_score);
	    vmessage("Initial score %.2f%% mismatches (%"PRId64
		     " mismatches)\n",
		     (100.0 * orig_score)/tot_score, orig_score/128);
	    if (flush)
		UpdateTextOutput();
	    //print_malign(malign);
	    do {
		old_score = new_score;
		malign = realign_seqs(cnum, malign, band, indels);
		//print_malign(malign);
		new_score = malign_diffs(malign, &tot_score);
		vmessage("  Consensus difference score: %"PRId64"\n",
			 new_score);
		if (flush)
		    UpdateTextOutput();
	    } while (new_score < old_score);

	    if (new_score < orig_score) {
		//print_malign(malign);
		update_io(io, cnum, malign, indels);

		/*
		 * It's possible the contig ends could move if a sequence that
		 * was previously the end of a contig has been moved such that
		 * it's no longer the contig end. This can lead to tags off the
		 * end of the contig, so trim them (reusing break_contig
		 * code).
		 */
		contig_visible_start(io, cnum, CITER_CSTART);
		contig_visible_end(io, cnum, CITER_CEND);
	    } else {
		vmessage("Could not reduce number of consensus "
			 "differences.\n");
	    }

	    destroy_malign(malign, 1);

	    vmessage("Final score %.2f%% mismatches\n",
		     (100.0 * new_score)/tot_score);

	    /*
	     * Sequences like
	     *   AGCT**GATGC
	     *             TGGATCGA
	     * can end up causing holes. We break the contig in this case to
	     * avoid minor database inconsistencies.
	     */
	    // remove_contig_holes(io, cnum);

	    /* reassign_confidence_values(io, cnum); */
	    //}

	    if (h_clips) {
		validate_clip_regions(io, h_clips, tag_arr);
		ArrayDestroy(tag_arr);
	    }

	    if (h_clips) {
		rewrite_soft_clips(io,
				   cl.contig,
				   cl.start,
				   cl.end,
				   h_clips);
		HashTableDestroy(h_clips, 1);
	    }

	    /* Remove pad columns */
	    if (soft_clips || new_score < orig_score) {
		remove_pad_columns(io, 1, &cl, 100, 1);
	    }

	    if (flush)
		cache_flush(io);
	    
	    sub_start -= CHUNK_SIZE;
	    sub_end   -= CHUNK_SIZE;
	} while (sub_end > contigs[i].start);

	/* Shift contig back */
	if (c_shift != 0) {
	    if (move_contig(io, contigs[i].contig, -c_shift) != 0)
		return -1;
	}
    }

    ArrayDestroy(indels);

    if (counts)
	free(counts);

    return 0;
}

/*
 * ----------------------------------------------------------------------
 * Remove Pad Columns. Sometimes we don't want to realign data, we just
 * want to remove (aligned) columns of pads.
 * ----------------------------------------------------------------------
 */
int remove_pad_columns(GapIO *io, int ncontigs, contig_list_t *contigs,
		       int percent_pad, int quiet) {
    int i;
    consensus_t *cons = NULL;
    size_t max_alloc = 0;

    for (i = 0; i < ncontigs; i++) {
	tg_rec cnum = contigs[i].contig;
	size_t len, j;
	int ndel = 0;
	contig_t *c;

	if (!quiet) {
	    vmessage("Processing contig %d of %d (#%"PRIrec")\n",
		     i+1, ncontigs, cnum);
	    UpdateTextOutput();
	}

	c = cache_search(io, GT_Contig, cnum);
	if (!c)
	    return -1;

	cache_incr(io, c);
	
	len = contigs[i].end - contigs[i].start + 1;
	if (max_alloc < len) {
	    max_alloc = len;
	    cons = realloc(cons, max_alloc * sizeof(*cons));
	}
	
	if (0 != calculate_consensus(io, cnum,
				     contigs[i].start, contigs[i].end,
				     cons)) {
	    free(cons);
	    cache_decr(io, c);
	    return -1;
	}

	for (j = 0; j < len; j++) {
	    if (cons[j].call != 4)
		continue;

	    if (100 * cons[j].counts[4] / cons[j].depth < percent_pad)
		continue;

	    if (!quiet)
		vmessage("  Removing column %d %d%% pad (%d of %d), conf. %f)\n",
			 (int)j+contigs[i].start,
			 100 * cons[j].counts[4] / cons[j].depth,
			 cons[j].counts[4], cons[j].depth,
			 cons[j].scores[cons[j].call]);

	    contig_delete_base(io, &c, contigs[i].start + j - ndel);
	    ndel++;
	}

	cache_decr(io, c);
    }

    if (cons)
	free(cons);

    return 0;
}


/*
 * ----------------------------------------------------------------------
 * Unclip matching data.
 *
 * This algorithm hunts down the softclipped data and builds a
 * histogram of values per consensus column. Any regions of high depth
 * and high coherency are deemed to be worthy of unclipping and
 * realigning.
 *
 * Almost always this ambiguity comes from misassemblies or collapsed
 * repeats, or at the very least it is valuable information we should
 * know about and tag.
 * ----------------------------------------------------------------------
 */

/*
 * Tags a region of consensus. 
 *
 * Returns the tag record number on success;
 *         -1 on failure
 */
tg_rec tag_softclip(GapIO *io, tg_rec crec, int start, int end,
		    int snp, double avg_depth, int (*depth)[7], int dir) {
    int j, d = 0;
    tg_rec r;
    char *comment = malloc(end-start+1 + 100), *cp;
    int type;

    if (!comment)
	return -1;

    cp = comment;
    if (depth) {
	cp += sprintf(comment, "SNPs=%d\nAvg. depth=%5.1f\n"
		      "Soft-clip consensus=", snp, avg_depth);
	for (j = start; j <= end; j++) {
	    *cp++ = (*depth++)[6];
	}
	*cp++ = 0;
	type = str2type("CLIP");
    } else {
	sprintf(comment, "Consensus N");
	type = str2type("NCLP");
    }

    r = anno_ele_add(io, GT_Contig, crec, 0, type, comment, start, end, dir);

    free(comment);

    return r;
}

/*
 * Returns a hash table of soft_clips structures, indexed on read names. 
 * To iterate use HashTableIterCreate.
 *
 * Also, if non-NULL, fills out clips array holding the tag recs.
 * These can be used to identify regions of interest for further
 * study.  Tags are added for both the coherent soft clips themselves
 * and also any Ns in consensus caused by contig gaps.  These are
 * important as we wish to preserve them unless the realignment is
 * good.
 *
 * Returns NULL on failure.
 *         Hash of soft_clips* on success; caller to free().
 */
HashTable *coherent_soft_clips(GapIO *io, tg_rec crec, int start, int end,
			       int *counts, int tag_only,
			       int min_depth, int min_tag_length,
			       Array *tag_arr) {
    seq_t *s;
    contig_iterator *citer;
    rangec_t *r;
    int (*Ldepth)[7]; // ACGTN* total
    int (*Rdepth)[7]; // ACGTN* total
    int i, j, changed, nc = 0, tag_alloc = 0;
    tg_rec *ctags = NULL;
    HashTable *h;
    int pass = 0;
    consensus_t *cons;

    static int L[256] = {
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //00
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //10
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 5, 4, 4, 4, 4, 4, //20
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //30
	4, 0, 4, 1, 4, 4, 4, 2,   4, 4, 4, 4, 4, 4, 4, 4, //40
	4, 4, 4, 4, 3, 3, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //50
	4, 0, 4, 1, 4, 4, 4, 2,   4, 4, 4, 4, 4, 4, 4, 4, //60
	4, 4, 4, 4, 3, 3, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //70
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //80
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //90
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //a0
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //b0
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //c0
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //d0
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, //e0
	4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4};//f0
    
    if (!(Ldepth = malloc((end - start + 1) * 7*sizeof(int))))
	return NULL;
    if (!(Rdepth = malloc((end - start + 1) * 7*sizeof(int)))) {
	free(Ldepth);
	return NULL;
    }

    h = HashTableCreate(128, HASH_DYNAMIC_SIZE);
    if (!h) {
	free(Ldepth);
	free(Rdepth);
	return NULL;
    }

    if (tag_arr) {
	if (!(*tag_arr = ArrayCreate(sizeof(tg_rec), 0)))
	    return NULL;
    }

    // Add N tags, to validate later
    if (!(cons = calloc(end-start+1, sizeof(*cons))))
	return NULL;
    calculate_consensus(io, crec, start, end, cons);
    for (i=start; i<end; i++) {
	if (cons[i-start].call == 5) {
	    tg_rec rec;
	    int j=i;
	    while(i < end && cons[i-start].call == 5)
		i++;
	    rec = tag_softclip(io, crec, j, --i, 0, 0, 0, '+');
	    if (tag_arr) ArrayPush(*tag_arr, tg_rec, rec);
	}
    }
    
 second_pass:
    memset(Ldepth, 0, (end - start + 1) * 7*sizeof(int));
    memset(Rdepth, 0, (end - start + 1) * 7*sizeof(int));
    changed = 0;

    /* Gather cutoff depth analysis */
    citer = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);

    while ((r = contig_iter_next(io, citer))) {
	seq_t *s, *sorig;

	s = sorig = cache_search(io, GT_Seq, r->rec);

	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	if (!common_word_L(counts, &s->seq[MAX(s->left-13, 0)],
			   MIN(12, s->left-1))) {
	    for (i = 0; i < s->left-1; i++) {
		if (r->start + i >= start &&
		    r->start + i <= end) {
		    if (s->conf[i] < 20)
			continue;
		    Ldepth[r->start+i - start][L[(unsigned char) s->seq[i]]]++;
		    Ldepth[r->start+i - start][6]++;
		}
	    }
	}

	if (!common_word_R(counts, &s->seq[s->right], ABS(s->len)-s->right)) {
	    for (i = s->right; i < ABS(s->len); i++) {
		if (r->start + i >= start &&
		    r->start + i <= end) {
		    if (s->conf[i] < 20)
			continue;
		    Rdepth[r->start+i - start][L[(unsigned char) s->seq[i]]]++;
		    Rdepth[r->start+i - start][6]++;
		}
	    }
	}

	if (s != sorig)
	    free(s);
    }

    contig_iter_del(citer);


    /* Compute cutoff consensus */
    for (j = 0; j < 2; j++) {
	int (*depth)[7] = j ? Rdepth : Ldepth;
	int tag_start, tag_depth = 0;

	for (i = start; i <= end; i++) {
	    int b, c = 0, m = 0, M = 0;
	    if (depth[i-start][6] < min_depth) {
		if (tag_depth && i-1 - tag_start + 1 >= min_tag_length && 
		    tag_depth/(i-1 - tag_start + 1.0) >= min_depth) {
		    tg_rec rec;
		    int snp = 0, x;
		    vmessage("Coherent %s softclip, length %5d depth %5.1f, "
			     "from %d to %d\n",
			     j ? "right" : "left", i-1 - tag_start + 1,
			     tag_depth/(i-1 - tag_start + 1.0),
			     tag_start, i-1);
		    for (x = start; x <= i-1; x++)
			if (cons[x-start].scores[6]>0)
			    snp++;
		    rec = tag_softclip(io, crec, tag_start, i-1, snp,
				       tag_depth/(i-1 - tag_start + 1.0),
				       &depth[tag_start-start], "-+"[j]);
		    if (tag_arr) ArrayPush(*tag_arr, tg_rec, rec);
		}

		depth[i-start][6] = 0;
		tag_depth = 0;
		continue;
	    }

	    if (m < depth[i-start][0]) M=m, m = depth[i-start][0],c=0; // A
	    else if (M < depth[i-start][0]) M = depth[i-start][0];

	    if (m < depth[i-start][1]) M=m, m = depth[i-start][1],c=1; // C
	    else if (M < depth[i-start][1]) M = depth[i-start][1];

	    if (m < depth[i-start][2]) M=m, m = depth[i-start][2],c=2; // G
	    else if (M < depth[i-start][2]) M = depth[i-start][2];

	    if (m < depth[i-start][3]) M=m, m = depth[i-start][3],c=3; // T
	    else if (M < depth[i-start][3]) M = depth[i-start][3];

	    if (m < depth[i-start][5]) M=m, m = depth[i-start][5],c=5; // *
	    else if (M < depth[i-start][5]) M = depth[i-start][5];

	    // At least 60% for 1 base
	    b = m*100 >= depth[i-start][6]*60 ? "ACGTN*"[c] : 'N';

	    // Or at least 90% for the top two base types.
	    if (b == 'N')
		b = (m+M)*100 >= depth[i-start][6]*90 ? "ACGTN*"[c] : 'N';

	    //printf("%6d: %2d %2d %2d %2d %2d %2d / %2d => %c\n",
	    //	   i,
	    //	   depth[i-start][0], depth[i-start][1], depth[i-start][2],
	    //	   depth[i-start][3], depth[i-start][4], depth[i-start][5],
	    //	   depth[i-start][6], b);

	    depth[i-start][6] = b;

	    if (b == 'N') {
		//printf("tag_depth %d,  start %d, len %d, avg_depth %f\n",
		//       tag_depth, tag_start, i-1 - tag_start+1,
		//       tag_depth/(i-1 - tag_start + 1.0));
		if (tag_depth && i-1 - tag_start + 1 >= min_tag_length && 
		    tag_depth/(i-1 - tag_start + 1.0) >= min_depth) {
		    tg_rec rec;
		    int snp = 0, x;
		    vmessage("Coherent %s softclip, length %5d depth %5.1f, "
			     "from %d to %d\n",
			     j ? "right" : "left", i-1 - tag_start + 1,
			     tag_depth/(i-1 - tag_start + 1.0),
			     tag_start, i-1);
		    for (x = start; x <= i-1; x++)
			if (cons[x-start].scores[6]>0)
			    snp++;
		    rec = tag_softclip(io, crec, tag_start, i-1, snp,
				       tag_depth/(i-1 - tag_start + 1.0),
				       &depth[tag_start-start], "-+"[j]);
		    if (tag_arr) ArrayPush(*tag_arr, tg_rec, rec);
		    tag_depth = 0;
		} else if (tag_depth) {
		    tag_depth = 0;
		}
	    } else {
		if (tag_depth == 0)
		    tag_start = i;
		tag_depth += depth[i-start][c];
	    }
	}

	if (tag_depth && i-1 - tag_start + 1 >= min_tag_length &&
	    tag_depth/(i-1 - tag_start + 1.0) >= min_depth) {
	    tg_rec rec;
	    int snp = 0, x;
	    vmessage("Coherent %s softclip, length %5d depth %5.1f, "
		     "from %d to %d\n",
		     j ? "right" : "left", i-1 - tag_start + 1,
		     tag_depth/(i-1 - tag_start + 1.0),
		     tag_start, i-1);
	    for (x = start; x <= i-1; x++)
		if (cons[x-start].scores[6]>0)
		    snp++;
	    rec = tag_softclip(io, crec, tag_start, i-1, snp,
			       tag_depth/(i-1 - tag_start + 1.0),
			       &depth[tag_start-start], "-+"[j]);
	    if (tag_arr) ArrayPush(*tag_arr, tg_rec, rec);
	}
    }


    if (tag_only) {
	free(Ldepth);
	free(Rdepth);
	return NULL;
    }


    /* Extend cutoffs where matching depth */
    citer = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);

    while ((r = contig_iter_next(io, citer))) {
	seq_t *s, *sorig;
	int new_l, new_r;
	int i_max, score, score_max;

	s = sorig = cache_search(io, GT_Seq, r->rec);
	new_l = s->left;
	new_r = s->right;

	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	// Left clip
	score = score_max = 0;
	for (i_max = i = s->left-2; i >= 0; i--) {
	    if (!(r->start + i >= start &&
		  r->start + i <= end))
		break;

	    if (s->seq[i] != Ldepth[r->start+i - start][6]) {
		if ((score-=3) < -6)
		    break;
	    } else {
		if (score_max < ++score)
		    score_max = score, i_max = i;
	    }
	}
	i = i_max;
	if (i < s->left-2) {
	    if (s == sorig)
		new_l = i+1;
	    else
		new_r = ABS(s->len) - i + 1;
	    
	    //printf("%"PRIrec"<%.*s\n", s->rec, s->left-2 -i, &s->seq[i+1]);
	}

	// Right clip
	score = score_max = 0;
	for (i_max = i = s->right; i < ABS(s->len); i++) {
	    if (!(r->start + i >= start &&
		  r->start + i <= end))
		break;

	    if (s->seq[i] != Rdepth[r->start+i - start][6]) {
		if ((score-=3) < -6)
		    break;
	    } else {
		if (score_max < ++score)
		    score_max = score, i_max = i+1;
	    }
	}
	i = i_max;
	if (i > s->right) {
	    if (s == sorig)
		new_r = i;
	    else
		new_l = ABS(s->len) - i + 1;
	    //printf("%"PRIrec">%.*s\n", s->rec, i - s->right, &s->seq[s->right]);
	}

	if (s != sorig)
	    free(s);

	/*
	 * This will produce inconsistencies like:
	 *   Seq 180892: left/right clips outside of sequence bounds.
	 *
	 * We will patch up the data later to fix these.
	 */
	if (new_r != sorig->right ||
	    new_l != sorig->left) {
	    HashData hd;
	    int new_rec;
	    soft_clips *c = malloc(sizeof(*c));

	    c->rec   = sorig->rec;
	    c->left  = depad_clip(sorig, sorig->left);
	    c->right = depad_clip(sorig, sorig->right);

	    hd.p = c;
	    HashTableAdd(h, (char *)&c->rec, sizeof(c->rec), hd, &new_rec);
	    if (!new_rec)
		free(c);

	    changed=1;

	    s = cache_rw(io, sorig);
	    s->right = new_r;
	    s->left  = new_l;
	}
    }

    contig_iter_del(citer);

    /*
     * We may have neighbouring blocks of coherent soft-clips due to SNPs.
     * This is easiest resolved with multiple passes.
     */
    if (++pass < 3 && changed)
	goto second_pass;

    free(Ldepth);
    free(Rdepth);
    free(cons);

    return h;
}


/*
 * Scans through a contig paying particular attention to the known
 * soft-clips.  We can ether extend a sequence if the soft clipped
 * data matches the consensus, or if it was previously extended by
 * coherent_soft_clip then we can increase soft-clipping back to the
 * former value if it disagrees with the new consensus.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int rewrite_soft_clips(GapIO *io, tg_rec crec, int start, int end,
		       HashTable *h_clips) {
    contig_iterator *citer;
    rangec_t *r;
    int i;
    char *cons;
    int updated = 0;

    vmessage("Extend soft-clips for contig =%"PRIrec" at %d..%d\n",
	     crec, start, end);

    /* Compute consensus to align against */
    if (!(cons = malloc(end - start + 1))) {
	return -1;
    }
    calculate_consensus_simple(io, crec, start, end, cons, NULL);

    citer = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);
    while ((r = contig_iter_next(io, citer))) {
	seq_t *sorig;
	HashItem *hi;
	seq_t *s = cache_search(io, GT_Seq, r->rec);
	int score = 0, best_score = 0, best_i = 0;
	int orig_right = s->right;
	int orig_left  = s->left;

	/*
	 * If this is a read we previously unclipped, then back up to
	 * that point and verify it from there.
	 *
	 * NOTE: We should perhaps remember the unpadded clip position
	 * and recompute the equivalent original s->right again as the
	 * read has been realigned since then with possibly fewer or
	 * more pads added. I think for now this is "good enough" though.
	 */
	if ((hi = HashTableSearch(h_clips, (char *)&r->rec, sizeof(r->rec)))) {
	    soft_clips *c = (soft_clips *)hi->data.p;
	    int p_l, p_r;
	    s = cache_rw(io, s);
	    
	    if (s->left < (p_l = repad_clip(s, c->left)))
		s->left = p_l;
	    if (s->right > (p_r = repad_clip(s, c->right)))
		s->right = p_r;
	    updated = 1;
	}

	// Right end
	if ((s->len<0) ^ r->comp) {
	    for (i = s->right; i < ABS(s->len); i++) {
		if (r->end -i -1 < start)
		    break;
		if (toupper(complement_base(s->seq[i])) ==
		    cons[r->end - i - start]) {
		    if (best_score < ++score) {
			best_score = score;
			best_i = i+1;
		    }
		} else {
		    if ((score -= 3) <= -6)
			break;
		}
		//printf("-%4d %7d: %c %c\n",
		//       i, r->end - i,
		//       complement_base(s->seq[i]),
		//       cons[r->end - i - start]);
	    }
	} else {
	    for (i = s->right; i < ABS(s->len); i++) {
		if (r->start + i > end)
		    break;
		if (toupper(s->seq[i]) == cons[r->start + i - start]) {
		    if (best_score < ++score) {
			best_score = score;
			best_i = i+1;
		    }
		} else {
		    if ((score -= 3) <= -6)
			break;
		}
		//printf("+%4d %7d: %c %c\n", i, r->start + i,
		//       s->seq[i], cons[r->start + i - start]);
	    }
	}

	if (best_score > 0) {
	    //vmessage("#%"PRIrec": Extend 5' end by %d\n",
	    //	     s->rec, best_i - s->right);
	    s = cache_rw(io, s);
	    while (best_i > 1 && s->seq[best_i-1] == '*')
		best_i--;
	    s->right = best_i;
	    updated = 1;
	}
	//if (s->right < orig_right) {
	//    printf("#%"PRIrec": Trim 5' end by %d\n",
	//	   s->rec, orig_right - s->right);
	//}


	// Left end
	best_i = best_score = score = 0;
	if ((s->len<0) ^ r->comp) {
	    for (i = s->left-2; i >= 0; i--) {
		if (r->end -i -1 > end)
		    break;
		if (toupper(complement_base(s->seq[i])) ==
		    cons[r->end - i - start]) {
		    if (best_score < ++score) {
			best_score = score;
			best_i = i+1;
		    }
		} else {
		    if ((score -= 3) <= -6)
			break;
		}
		//printf("-%4d %7d: %c %c %d\n",
		//       i, r->end - i,
		//       complement_base(s->seq[i]),
		//       cons[r->end - i - start], score);
	    }
	} else {
	    for (i = s->left-2; i >= 0; i--) {
		if (r->start + i < start)
		    break;
		if (toupper(s->seq[i]) == cons[r->start + i - start]) {
		    if (best_score < ++score) {
			best_score = score;
			best_i = i+1;
		    }
		} else {
		    if ((score -= 3) <= -6)
			break;
		}
		//printf("+%4d %7d: %c %c %d\n", i, r->start + i,
		//       s->seq[i], cons[r->start + i - start], score);
	    }
	}

	if (best_score > 0) {
	    //vmessage("#%"PRIrec": Extend 3' end by %d\n",
	    //	     s->rec, s->left - best_i);
	    s = cache_rw(io, s);
	    while (best_i < s->right && s->seq[best_i-1] == '*')
		best_i++;
	    s->left = best_i;
	    updated = 1;
	}
	//if (s->left > orig_left) {
	//    printf("#%"PRIrec": Trim 3' end by %d\n",
	//	   s->rec, s->left - orig_left);
	//}

	// Update the range? Not needed as start/end haven't changed.
	// However the consensus valid range flag may be incorrect.
    }
    contig_iter_del(citer);

    free(cons);

    return 0;
}

/*
 * Does a word usage scan on soft-clips to try and identify likely
 * adapter sequences.
 *
 * We compare common words in cutoffs vs common words in used portions
 * and identify the discrepancies.
 */
#define ADAPTER_WORD 12
#define ADAPTER_SIZE (1<<(2*ADAPTER_WORD))
#define ADAPTER_MASK (ADAPTER_SIZE-1)
static int L[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*   0-15 */
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  16 */
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  32 */
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  48 */
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, /*  64 */
    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  80 */
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, /*  96 */
    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 112-127 */
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 128 */
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  /* 255 */
};

static unsigned int hash_word(char *seq) {
    unsigned int w = 0, i;
    for (i = 0; i < ADAPTER_WORD; i++) {
	w <<= 2;
	w |= L[(unsigned char)seq[i]];
    }
    return w;
}

static char *W(unsigned int w) {
    static char buf[ADAPTER_WORD+1];
    int i;
    for (i = 0; i < ADAPTER_WORD; i++)
	buf[i] = "ACGT"[(w>>(2*(ADAPTER_WORD-i-1)))&3];
    buf[ADAPTER_WORD]=0;
    return buf;
}

static int common_word_L(int *counts, char *seq, int len) {
    unsigned int w = 0, i;
    if (len < 4 || !counts)
	return 0;

    for (i = 0; i < ADAPTER_WORD && i < len; i++) {
	w <<= 2;
	w |= L[(unsigned char)seq[i]];
    }
    return counts[w];
}

static int common_word_R(int *counts, char *seq, int len) {
    unsigned int w = 0, i;
    if (len < 4 || !counts)
	return 0;

    for (i = 0; i < ADAPTER_WORD && i < len; i++) {
	w <<= 2;
	w |= L[(unsigned char)seq[i]];
    }
    return counts[w << (2*(ADAPTER_WORD - MIN(ADAPTER_WORD,len)))];
}

/*
 * Returns a malloced array of ADAPTER_WORD bases long holding 1 for an
 * unusually common word and 0 for a normal/expected word usage.
 *
 * Returns NULL on failure.
 */
int *find_adapter(GapIO *io, int ncontigs, contig_list_t *contigs) {
    int i, j;
    int *counts_clip, *counts_used;
    uint64_t clip_tot = 0;
    uint64_t used_tot = 0;
    uint64_t t1, t2;

    counts_clip = calloc(ADAPTER_SIZE, sizeof(int));
    counts_used = calloc(ADAPTER_SIZE, sizeof(int));

    for (i = 0; i < ncontigs; i++) {
	tg_rec crec = contigs[i].contig;
	int start = contigs[i].start;
	int end = contigs[i].end;
	contig_iterator *citer;
	rangec_t *r;
	unsigned int w;

	citer = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);
	while ((r = contig_iter_next(io, citer))) {
	    seq_t *s, *sorig;

	    sorig = s = cache_search(io, GT_Seq, r->rec);

	    if (s->left < ADAPTER_WORD+1 &&
		ABS(s->len) - s->right < ADAPTER_WORD)
		continue;

	    if ((s->len < 0) ^ r->comp) {
		s = dup_seq(s);
		complement_seq_t(s);
	    }

	    // FIXME: needs to work on depadded sequence.
	    if (s->left > ADAPTER_WORD) {
		//printf("#%"PRIrec" L %.*s\n",
		//       s->rec, ADAPTER_WORD, &s->seq[s->left-1-ADAPTER_WORD]);
		counts_clip[hash_word(&s->seq[s->left-1-ADAPTER_WORD])]++;
		clip_tot++;
	    }

	    if (ABS(s->len) - s->right >= ADAPTER_WORD) {
		//printf("#%"PRIrec" R %.*s\n",
		//       s->rec, ADAPTER_WORD, &s->seq[s->right]);
		counts_clip[hash_word(&s->seq[s->right])]++;
		clip_tot++;
	    }

	    // First and last word of used portion.
	    if (s->right - s->left > ADAPTER_WORD) {
		w = hash_word(&s->seq[s->left-1]);
		//printf("#%"PRIrec" M %s\n", s->rec, W(w));
		counts_used[w]++;
		w = hash_word(&s->seq[s->right-ADAPTER_WORD]);
		//printf("#%"PRIrec" M %s\n", s->rec, W(w));
		counts_used[w]++;
		used_tot+=2;
	    }

//	    w = hash_word(&s->seq[s->left-1]);
//	    j = s->left-1 + ADAPTER_WORD;
//	    do {
//		counts_used[w]++;
//		//printf("#%"PRIrec" M %s\n", s->rec, W(w));
//		w <<= 2;
//		w |= L[s->seq[j]];
//		w &= (1<<(2*ADAPTER_WORD))-1;
//		used_tot++;
//	    } while (++j <= s->right);

	    if (s != sorig)
		free(s);
	}
	contig_iter_del(citer);
    }

    // Filter to common words in clips only
    t1 = clip_tot * 0.01;  // 1%
    t2 = used_tot * 0.005; // 0.5%
    if (clip_tot > 1000 && used_tot > 1000) {
	for (i = 0; i < ADAPTER_SIZE; i++) {
	    //if (counts_clip[i] > t1)
	    //	printf("Clip: %s %5.1f\n",
	    //	       W(i), 100.0 * counts_clip[i]/clip_tot);
	    //if (counts_used[i] > t2)
	    //	printf("Used: %s %5.1f\n",
	    //	       W(i), 100.0 * counts_used[i]/used_tot);

	    if (counts_clip[i]>t1 && counts_used[i]<t2) {
		counts_clip[i] = 1;
		vmessage("Discarding word %s as likely adpater (%5.1f%%)\n",
			 W(i), 100.0 * counts_clip[i]/clip_tot);
	    } else {
		counts_clip[i] = 0;
	    }
	}
    } else {
	memset(counts_clip, 0, ADAPTER_SIZE * sizeof(int));
    }

    // Expand clipped words to partial matches, down to 4 bp.
    for (i = 0; i < ADAPTER_SIZE; i++) {
	int j;
	if (!counts_used[i])
	    continue;
	
	for (j = 4; j < ADAPTER_WORD; j++) {
	    counts_used[i & ((1<<(2*j))-1)] = 1;
	    counts_used[(i << (2*(ADAPTER_WORD-j))) & ADAPTER_MASK] = 1;
	}
    }

    free(counts_used);

    return counts_clip;
}
