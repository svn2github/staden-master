/*
 * TODO:
 * - proper deallocation system (contigs, bins, seqs, annos)
 * - tags on consensus need to move if we're breaking a contig, dup otherwise?
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <io_lib/hash_table.h>

#include "dis_readings.h"
#include "misc.h"
#include "break_contig.h"
#include "bitmap.h"
#include "consensus.h"
#include "text_output.h"

/*-----------------------------------------------------------------------------
 * Databse checking code.
 */

/* Recursively check linkage on bins and nseqs relationships.
 * brec  = bin record
 * prec  = parent record, to validate against.
 * ptype = type of parent record.
 *
 * Returns nseqs on success
 *         -1 on failure.
 */
//static FILE *errfp = stderr;
static FILE *errfp = NULL;;

static int check_contig_bins_r(GapIO *io, tg_rec brec, int ptype, tg_rec prec){
    bin_index_t *b;
    tg_rec copy_c0, copy_c1;
    int copy_nseq, i, ns;

    /* Check prec/ptype */
    b = cache_search(io, GT_Bin, brec);

    if (b->parent != prec || b->parent_type != ptype) {
	fprintf(errfp, "ERROR: bin parent record/type mismatch for bin %"
		PRIrec" : parent = %"PRIrec"/%"PRIrec" type = %d/%d\n",
		brec, b->parent, prec, b->parent_type, ptype);
	abort();
	return -1;
    }

    /* Count number of sequences in this contig */
    for (ns = i = 0; b->rng && i < ArrayMax(b->rng); i++) {
	range_t *r = arrp(range_t, b->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
	    ns++;
    }

    /* Take copies of the data from b so we don't have to lock it and
     * can allow it to get purged from the hache.
     */
    copy_c0   = b->child[0];
    copy_c1   = b->child[1];
    copy_nseq = b->nseqs;

    /* Recurse and check nseqs */
    if (copy_c0) {
	if ((i = check_contig_bins_r(io, copy_c0, GT_Bin, brec)) == -1)
	    return -1;
	ns += i;
    }
    if (copy_c1) {
	if ((i = check_contig_bins_r(io, copy_c1, GT_Bin, brec)) == -1)
	    return -1;
	ns += i;
    }

    if (ns != copy_nseq) {
	fprintf(errfp, "ERROR: nseq mismatch for bin %"PRIrec" : %d/%d\n",
		brec, ns, copy_nseq);
	abort();
	return -1;
    }

    return ns;
}

int check_contig_bins(GapIO *io) {
    int i, ret = 0;

    errfp = stdout;

    printf("check_contig_bins start, ncontigs=%d\n", io->db->Ncontigs);
    if (io->db->Ncontigs <= 340)
	return 0;

    for (i = 0; i < io->db->Ncontigs; i++) {
	tg_rec crec = arr(tg_rec, io->contig_order, i);
	contig_t *c = cache_search(io, GT_Contig, crec);
	//printf("Check contig %d root %d\n", crec, c->bin);
	if (c->bin) {
	    if (check_contig_bins_r(io, c->bin, GT_Contig, crec) == -1)
		ret = -1;
	}
    }

    printf("check_contig_bins end => %d\n", ret);

    return ret;
}

/* Singular of above */
int check_contig_bin(GapIO *io, tg_rec crec) {
    contig_t *c = cache_search(io, GT_Contig, crec);
    errfp = stdout;

    printf("Check contig %"PRIrec" root %"PRIrec"\n", crec, c->bin);
    if (c->bin) {
	if (check_contig_bins_r(io, c->bin, GT_Contig, crec) == -1)
	    return -1;
    }

    return 0;
}



/*-----------------------------------------------------------------------------
 * Internal functions & data types
 */

typedef struct {
    tg_rec contig;
    tg_rec bin;
    int start;
    int end;
    int clipped_start;
    int clipped_end;
    int comp;
    tg_rec rec;       /* sequence record */
    range_t rng;      /* rng in original bin */
    rangec_t *anno;   /* Annotation ranges */
    int n_anno;
} r_pos_t;

/*
 * Source and destination maps for contigs. This allows us to work out
 * which contigs are best to move consensus annotations too based on size
 * of overlap. Alternatively if we choose to implement duplication of
 * consensus annotations later then we'll also need this information.
 */
typedef struct {
    tg_rec src;
    int src_start; /* clipped */
    int src_end;
    tg_rec dest;
    int dest_start; /* unclipped, but we know clipped starts at bp 1 */
    int dest_end;
} contig_map;


/*
 * Part 1/2 of the disassemble_readings implementation (see below).
 *
 * 1. Produce a table of which contig number and region each reading
 *    record is within. Also track the anno_eles for this seq too.
 *
 * 2. Remove these readings from their contigs.
 *
 * 3. If "remove" is true, ensure if it's a read-pair that the other
 *    end gets the pairing information and flags cleared.
 *
 * Writes back to pos.
 *
 * The 'remove' argument is a boolean indicating whether the sequence
 * is about to be removed from the database. (Ie move==0 in
 * disassemble_reads.) We use this to determine whether to update the
 * read-pairing statuses.
 *
 * FIXME: we need to do something with annotations too. Unlink them, but also
 * keep track of their original locations relative to the sequence location
 * so we can add them back if required (move==1 or move==2).
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int unlink_read(GapIO *io, tg_rec rec, r_pos_t *pos, int remove) {
    contig_t *c;
    bin_index_t *bin;
    seq_t *seq;

    //printf("%soving record #%d\n", remove ? "Rem" : "M", rec);

    /*
     * Check if brec hasn't changed?
     * If not, then we can cache the bin location and speed up the
     * bin_get_item_position. However to do this we must be sure that
     * bin locations aren't moving, so we restrict this to specific
     * algorithms.
     */

    /* Get location */
    if (bin_get_item_position(io, GT_Seq, rec,
			      &pos->contig,
			      &pos->start,
			      &pos->end,
			      &pos->comp,
			      &pos->bin,
			      &pos->rng,
			      (void **)&seq)) {
	return -1;
    }
    /* seq already has cache_incr on it */

    /* Catch 0 length reads as these can cause bugs to occur */
    if (seq->right < seq->left) {
	seq = cache_rw(io, seq);

	if (seq->left > 1)
	    seq->left--;
	else
	    seq->right++;

	if ((seq->len < 0) ^ pos->comp) {
	    pos->clipped_start =
		pos->start + ABS(seq->len) - (seq->right-1) - 1;
	    pos->clipped_end   =
		pos->start + ABS(seq->len) - (seq->left-1) - 1;
	} else {
	    pos->clipped_start =
		pos->start + seq->left-1;
	    pos->clipped_end   =
		pos->start + seq->right-1;
	}
    }

    if ((seq->len < 0) ^ pos->comp) {
	pos->clipped_start = pos->start + ABS(seq->len) - (seq->right-1) - 1;
	pos->clipped_end   = pos->start + ABS(seq->len) - (seq->left-1) - 1;
    } else {
	pos->clipped_start = pos->start + seq->left-1;
	pos->clipped_end   = pos->start + seq->right-1;
    }

    /*
     * We compute the list of annotations for this read later on, when we
     * know which other reads are covering the same contig region.
     */
    pos->n_anno = 0;
    pos->anno = NULL;

    /*
     * Remove from bin range array. Delay consistency checking until
     * we've removed everything.
     */
    if (!(bin = cache_search(io, GT_Bin, pos->bin))) {
	cache_decr(io, seq);
	return -1;
    }
    cache_incr(io, bin);

    if (!(c = cache_search(io, GT_Contig, pos->contig))) {
	cache_decr(io, seq);
	cache_decr(io, bin);
	return -1;
    }
    cache_incr(io, c);

    if (fast_remove_item_from_bin(io, &c, &bin, GT_Seq, rec, seq->bin_index)) {
	cache_decr(io, seq);
	cache_decr(io, bin);
	cache_decr(io, c);
	return -1;
    }

    cache_decr(io, seq);
    cache_decr(io, bin);
    cache_decr(io, c);

    /*
     * For read-pairs, unlink with rest of template.
     * If we're removing huge volumes of data (eg a repeat) with the other
     * ends being scattered to the four winds, this will be very slow.
     *
     * Maybe we are better off not removing the link and just handling
     * the subsequent errors we will get later on during normal gap5
     * usage.
     */
    if (remove && pos->rng.pair_rec &&
	(seq = cache_search(io, GT_Seq, pos->rng.pair_rec))) {
	range_t *r;

	if (seq->parent_rec == pos->rng.pair_rec) {
	    seq = cache_rw(io, seq);
	    seq->parent_type = 0;
	    seq->parent_rec = 0;
	}

	/* Pair is held in bin range too */
	bin = cache_search(io, GT_Bin, seq->bin);
	r = arrp(range_t, bin->rng, seq->bin_index);
	assert(r->rec == seq->rec);
	bin = cache_rw(io, bin);
	bin->flags |= BIN_RANGE_UPDATED;

	/* Fix other end's range_t */
	r->pair_rec = 0;
	r->flags &= ~GRANGE_FLAG_TYPE_MASK;
	r->flags |=  GRANGE_FLAG_TYPE_SINGLE;
	r->pair_timestamp = 0;
    }

    return 0;
}

void bin_destroy_recurse(GapIO *io, tg_rec rec) {
    bin_index_t *bin = cache_search(io, GT_Bin,rec);

    cache_incr(io, bin);
    if (bin->child[0]) bin_destroy_recurse(io, bin->child[0]);
    if (bin->child[1]) bin_destroy_recurse(io, bin->child[1]);
    cache_decr(io, bin);

    cache_rec_deallocate(io, GT_Bin, rec);
}

/*
 * Checks if a bin is truely empty. It's not just sufficient to check
 * nseqs==0 as we could have tags or refpos markers too.
 * This differs to bin_empty() in that we are working out whether there is
 * no data in this bin and the child bins. Whereas the former is simply
 * interested in the bin->rng contents, so data within that specific bin.
 *
 * Returns 1 if empty.
 *         0 if not.
 */
static int bin_plus_children_empty(bin_index_t *bin) {
    int i;

    if (bin->nseqs || bin->nrefpos || bin->nanno)
	return 0;

    if (!bin->rng)
	return 1;

    for (i = 0; i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);

	if (r->flags & (GRANGE_FLAG_ISCONS | GRANGE_FLAG_ISREF))
	    continue;

	if (!(r->flags & GRANGE_FLAG_UNUSED))
	    return 0;
    }

    return 1;
}

/*
 * Looks for contig gaps between start..end in contig and if it finds them,
 * breaking the contig in two.
 */
int remove_contig_holes(GapIO *io, tg_rec contig, int start, int end,
			int empty_contigs_only) {
    contig_t *c;
    bin_index_t *bin;
    contig_iterator *iter;
    rangec_t *r;
    int last;
    int contig_start, contig_end;

    /* Destroy contigs if they're now entirely empty */
    c = cache_search(io, GT_Contig, contig);
    cache_incr(io, c);

    bin = cache_search(io, GT_Bin, c->bin);
    if (bin_plus_children_empty(bin)) {
	puts("Removing empty contig");

	if (c->bin)
	    bin_destroy_recurse(io, c->bin);
	c->timestamp = io_timestamp_incr(io);
	cache_decr(io, c);
	contig_destroy(io, contig);
	return 0;
    }

    /* Used fast_remove_item_from_bin => invalidate the read pairs posn */
    c->timestamp = io_timestamp_incr(io);

    /* Invalidate any cached consensus copies */
    if (bin_invalidate_consensus(io, contig, start, end) != 0) {
	cache_decr(io, c);
	return -1;
    }

    /* Hole at left end */
    if (c->start == start) {
	iter = contig_iter_new(io, contig, 1, CITER_FIRST, start, end);
	if (iter) {
	    r = contig_iter_next(io, iter);
	    if (r) {
		c = cache_rw(io, c);
		start = c->start = r->start;
	    }
	    contig_iter_del(iter);
	}
    }

    /* Hole at right end */
    if (c->end == end) {
	iter = contig_iter_new(io, contig, 1, CITER_LAST | CITER_IEND,
			       start, end);
	if (iter) {
	    r = contig_iter_prev(io, iter);
	    if (r) {
		c = cache_rw(io, c);
		end = c->end = r->end;
	    }
	    contig_iter_del(iter);
	}
    }

    if (empty_contigs_only) {
	cache_decr(io, c);
	return 0;
    }


    /* Make sure start/end are within clipped contig coordinates */
    consensus_valid_range(io, contig, &contig_start, &contig_end);
    if (start < contig_start)
	start = contig_start;
    if (end > contig_end)
	end = contig_end;


    /*
     * Look for holes in the middle, using ICLIPPEDEND mode so the data
     * is sorted by clipped sequence coords rather than just the r->end
     * rangec_t elements we use with most iterators.
     */
    iter = contig_iter_new(io, contig, 0, CITER_LAST | CITER_ICLIPPEDEND,
			   start, end);
    if (!iter) {
	cache_decr(io, c);
	return 0;
    }

    last = end;
    while (iter && (r = contig_iter_prev(io, iter))) {
	seq_t *s = cache_search(io, GT_Seq, r->rec);
	int cstart, cend;

	if (!s) {
	    cache_decr(io, c);
	    return -1;
	}

	if ((s->len < 0) ^ r->comp) {
	    cstart = r->start + ABS(s->len) - (s->right-1) - 1;
	    cend   = r->start + ABS(s->len) - (s->left-1) - 1;
	} else {
	    cstart = r->start + s->left-1;
	    cend   = r->start + s->right-1;
	}

	//printf("Seq %d, from %d..%d clipped %d..%d\n",
	//       r->rec, r->start, r->end, cstart, cend);

	if (cend < last) {
	    int r;

	    vmessage("GAP from %d..%d; breaking.\n", cend, last);
	    if (!empty_contigs_only)
		r = break_contig(io, contig, last, 0);
	    else
		r = 0;

	    /* Who knows what impact break_contig has - restart to be safe */
	    contig_iter_del(iter);
	    if (r == -1) {
		cache_decr(io, c);
		return -1;
	    }

	    iter = contig_iter_new(io, contig, 0,
				   CITER_LAST | CITER_ICLIPPEDEND,
	    			   start, last);
	}
	if (last > cstart)
	    last = cstart;
    }
    if (iter)
	contig_iter_del(iter);

    cache_decr(io, c);

    return 0;
}

static GapIO *xio = NULL;


/* qsort callback */
static int pos_sort(const void *vp1, const void *vp2) {
    const r_pos_t *p1 = (const r_pos_t *)vp1;
    const r_pos_t *p2 = (const r_pos_t *)vp2;

#if 0
    /* For stable sorting regardless of rec deallocation, use this: */
    if (p1->contig != p2->contig) {
	contig_t *c1 = cache_search(xio, GT_Contig, p1->contig);
	contig_t *c2 = cache_search(xio, GT_Contig, p2->contig);
	return strcmp(c1->name, c2->name);
    }
#endif

    if (p1->contig != p2->contig)
	return p1->contig - p2->contig;

    return p1->start - p2->start;
}

static int pos_sort_end(const void *vp1, const void *vp2) {
    const r_pos_t *p1 = (const r_pos_t *)vp1;
    const r_pos_t *p2 = (const r_pos_t *)vp2;

#if 0
    if (p1->contig != p2->contig) {
	contig_t *c1 = cache_search(xio, GT_Contig, p1->contig);
	contig_t *c2 = cache_search(xio, GT_Contig, p2->contig);
	return strcmp(c1->name, c2->name);
    }
#endif

    if (p1->contig != p2->contig)
	return p1->contig - p2->contig;

    return p1->end - p2->end;
}


/*
 * Part 4 of the disassemble_readings implementation.
 * 
 * 4. Check for holes within the "source" contigs. We use the
 *    initial table for this, along with regions to look around.
 *
 * If empty_only is true we only check for entirely empty contigs (and
 * deallocate them), otherwise we also break contig where holes exist.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int fix_holes(GapIO *io, r_pos_t *pos, int npos,
		     int remove_holes) {
    int i, start, end;
    tg_rec contig;
    int cstart, cend;
    contig_t *c;

    if (npos == 0)
	return 0;

    /*
     * pos[] is sorted by start coordinate. We need to loop backwards
     * through it though as break_contig will be producing a new right
     * hand contig, so moving backwards means we're always marching
     * down the un-changing (in contig number terms) left-hand contig.
     *
     * Due to sorted start coord but ragged end coord, we sort by
     * end instead so the overlapping segment detector works.
     */
    qsort(pos, npos, sizeof(*pos), pos_sort_end);

    /*
     * Step through pos finding overlapping reads so we can get entire
     * spans where we've removed data. We use this to limit our hole
     * fixing to just that region, reducing the search time on huge
     * contigs.
     */
    contig = pos[npos-1].contig;
    start  = pos[npos-1].start;
    end    = pos[npos-1].end;
    c = cache_search(io, GT_Contig, contig);
    cend = MAX(c->end, end);
    contig_visible_end(io, contig, cend);
    
    for (i = npos-1; i >= 0; i--) {
	if (pos[i].contig != contig) {
	    // Also trims end-tags
	    c = cache_search(io, GT_Contig, contig);
	    cstart = MIN(c->start, start);
	    contig_visible_start(io, contig, cstart); 
	}
	if (pos[i].contig != contig ||
	    pos[i].end < start) {
	    remove_contig_holes(io, contig, start, end, !remove_holes);
	    start  = pos[i].start;
	    end    = pos[i].end;

	    if (pos[i].contig != contig) {
		contig = pos[i].contig;
		c = cache_search(io, GT_Contig, contig);
		cend = MAX(c->end, end);
		contig_visible_end(io, contig, cend);
	    }
	} else {
	    if (start > pos[i].start)
		start = pos[i].start;
	}
    }

    c = cache_search(io, GT_Contig, contig);
    cstart = MIN(c->start, start);
    contig_visible_start(io, contig, cstart); // Also trims end-tags
    remove_contig_holes(io, contig, start, end, !remove_holes);

    return 0;
}


/*
 * Removes a sequence from the database, removing the bin range entry
 * and taking out of the BTree name index too if appropriate.
 */
static int seq_deallocate(GapIO *io, r_pos_t *pos) {
    int i;
    seq_t *s;
    bin_index_t *b;
    range_t *r;
    contig_t *c;
    tg_rec root;

    /* Remove from sequence name btree index */
    if (!(s = cache_search(io, GT_Seq, pos->rec)))
	return -1;
    cache_incr(io, s);
    
    /* Remove from relevant seq_block array */
    cache_item_remove(io, GT_Seq, pos->rec);

    /* Remove name,rec pair from b+tree */
    root = io->iface->seq.index_del(io->dbh, s->name, s->rec);
    if (root != -1 && root != io->db->seq_name_index) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = root;
    }

    /* Deallocate seq struct itself */
    if (!(b = cache_search(io, GT_Bin, s->bin))) {
	cache_decr(io, s);
	return -1;
    }

    /* Remove from range array */
    r = arrp(range_t, b->rng, s->bin_index);
    if (!(r->flags & GRANGE_FLAG_UNUSED)) {
	assert(r->rec == s->rec);

	b = cache_rw(io, b);
	b->flags |= BIN_RANGE_UPDATED;
	r->flags |= GRANGE_FLAG_UNUSED;

	if (b->start_used == r->start || b->end_used == r->end)
	    bin_set_used_range(io, b);
    }

    cache_decr(io, s);

    //Already achieved via cache_item_remove
    //cache_rec_deallocate(io, GT_Seq, pos->rec);

    /* Remove annotations too */
    c = cache_search(io, GT_Contig, pos->contig);
    cache_incr(io, c);
    for (i = 0; i < pos->n_anno; i++) {
	bin_remove_item(io, &c, GT_AnnoEle, pos->anno[i].rec);
	cache_item_remove(io, GT_AnnoEle, pos->anno[i].rec);
	//cache_rec_deallocate(io, GT_AnnoEle, pos->anno[i].rec);
    }
    cache_decr(io, c);

    return 0;
}

/*
 * Form a new contig from the reads in pos[].
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int create_contig_from(GapIO *io, r_pos_t *pos, int npos,
			      Array cmap) {
    int i;
    contig_t *c_old, *c_new;
    char name[8192];
    static int last_count=1;
    static tg_rec last_contig=-1;
    int offset;
    bin_index_t *bin;
    int old_comp;
    int dest_start = INT_MAX, dest_end = INT_MIN;
    int src_start = INT_MAX, src_end = INT_MIN;
    contig_map *map;

    if (npos <= 0)
	return -1;

    vmessage("\n=== new contig ===\n");
    for (i = 0; i < npos; i++) {
	vmessage("%d\tCtg %"PRIrec"\t%d..%d\tseq %"PRIrec"\n",
		 i, pos[i].contig, pos[i].start, pos[i].end, pos[i].rec);
    }


    /* Pick a new contig name */
    c_old = cache_search(io, GT_Contig, pos[0].contig);
    cache_incr(io, c_old);
    if (last_contig != pos[0].contig) {
	last_contig = pos[0].contig;
	last_count = 1;
    }
    do {
	sprintf(name, "%s%%%d", contig_get_name(&c_old), last_count++);
    } while (contig_index_query(io, name) != -1);

    bin = cache_search(io, GT_Bin, c_old->bin);
    old_comp = (bin->flags & BIN_COMPLEMENTED) ? 1 : 0;


    /* Create the new contig */
    c_new = contig_new(io, name);
    cache_incr(io, c_new);



    /* Add in the sequences */
    offset = INT_MAX;
    for (i = 0; i < npos; i++) {
	if (offset > pos[i].clipped_start)
	    offset = pos[i].clipped_start;
    }
    offset--;
    for (i = 0; i < npos; i++) {
	range_t r, *r_out;
	seq_t *s;
	int j;

	r = pos[i].rng;
	r.start = pos[i].start - offset;
	r.end   = pos[i].end - offset;

	if (src_start > pos[i].clipped_start)
	    src_start = pos[i].clipped_start;
	if (src_end < pos[i].clipped_end)
	    src_end = pos[i].clipped_end;

	if (dest_start > r.start)
	    dest_start = r.start;
	if (dest_end < r.end)
	    dest_end = r.end;

	r.y     = 0;
	if (pos[i].comp)
	    r.flags ^= GRANGE_FLAG_COMP1;

	bin = bin_add_range(io, &c_new, &r, &r_out, NULL, 0);

	s = cache_search(io, GT_Seq, pos[i].rec);
	s = cache_rw(io, s);
	s->bin = bin->rec;
	s->bin_index = r_out - ArrayBase(range_t, bin->rng);
	if (pos[i].comp) {
	    s->len = -s->len;
	    s->flags ^= SEQ_COMPLEMENTED;
	}

	/* Similarly move the annotations too; much the same as seqs */
	for (j = 0; j < pos[i].n_anno; j++) {
	    anno_ele_t *a;

	    //printf("Seq %d; pos %d,   anno %d; pos %d-%d\n",
	    //	   pos[i].rec,
	    //	   pos[i].start,
	    //	   pos[i].anno[j].rec,
	    //	   pos[i].anno[j].start,
	    //	   pos[i].anno[j].end);

	    bin_remove_item(io, &c_old, GT_AnnoEle, pos[i].anno[j].rec);
	    r.start    = pos[i].anno[j].start - offset;
	    r.end      = pos[i].anno[j].end - offset;
	    r.rec      = pos[i].anno[j].rec;
	    r.mqual    = pos[i].anno[j].mqual;
	    r.pair_rec = pos[i].anno[j].pair_rec;
	    r.flags    = pos[i].anno[j].flags;
	    
	    bin = bin_add_to_range(io, &c_new, s->bin, &r, &r_out, NULL, 0);
	    a = cache_search(io, GT_AnnoEle, pos[i].anno[j].rec);
	    a = cache_rw(io, a);
	    a->bin = bin->rec;
	    //a->bin_idx = r_out - ArrayBase(range_t, bin->rng);
	}
    }

    ArrayRef(cmap, ArrayMax(cmap));
    map = arrp(contig_map, cmap, ArrayMax(cmap)-1);
    map->src        = c_old->rec;
    map->src_start  = src_start;
    map->src_end    = src_end;

    map->dest       = c_new->rec;
    map->dest_start = dest_start;
    map->dest_end   = dest_end;

    cache_decr(io, c_old);
    cache_decr(io, c_new);

    return 0;
}

/*
 * Moves a bunch of reads to a new contig. We determine the set of reads
 * based on overlaps. Calls create_contig_from to do the grunt work.
 */
static int move_reads(GapIO *io, r_pos_t *pos, int npos, Array cmap) {
    int i, start, end;
    tg_rec contig;
    int i_start, err = 0;

    if (!npos)
	return 0;

    i_start = 0;
    contig  = pos[0].contig;
    start   = pos[0].start;
    end     = pos[0].end;
    for (i = 1; i < npos; i++) {
	if (pos[i].contig != contig ||
	    pos[i].start > end) {
	    if (create_contig_from(io, &pos[i_start], i - i_start, cmap))
		err = 1;
	    i_start = i;
	    contig  = pos[i].contig;
	    start   = pos[i].start;
	    end     = pos[i].end;
	} else {
	    if (end < pos[i].end)
		end = pos[i].end;
	}
    }
    if (create_contig_from(io, &pos[i_start], i - i_start, cmap))
	err = 1;

    return err;
}


static Bitmap contig_hole_bitmap(GapIO *io, tg_rec contig,
				 int start, int end) {
    contig_iterator *iter;
    rangec_t *r;
    Bitmap hole = BitmapCreate(end - start + 1);
    int last;

    if (NULL == hole) return NULL;

    iter = contig_iter_new(io, contig, 0,
			   CITER_LAST | CITER_ICLIPPEDEND,
			   start, end);
    if (!iter) goto fail;

    last = end+1;
    while (NULL != (r = contig_iter_prev(io, iter))) {
	seq_t *s = cache_search(io, GT_Seq, r->rec);
	int cstart, cend;
		
	if (!s) goto fail;

	if ((s->len < 0) ^ r->comp) {
	    cstart = r->start + ABS(s->len) - (s->right-1) - 1;
	    cend   = r->start + ABS(s->len) - (s->left-1) - 1;
	} else {
	    cstart = r->start + s->left-1;
	    cend   = r->start + s->right-1;
	}

	if (cend < last) {
	    int i;
	    for (i = cend+1; i < last; i++)
		if (i >= start && i <= end)
		    BIT_SET(hole, i-start);
	}
	if (last > cstart)
	    last = cstart;
    }

    if (start < last) {
	int i;
	for (i = start; i < last; i++)
	    if (i >= start && i <= end)
		BIT_SET(hole, i-start);
    }

    contig_iter_del(iter);

    return hole;
 fail:
    if (NULL != iter) contig_iter_del(iter);
    if (NULL != hole) BitmapDestroy(hole);
    return NULL;
}

/*
 * Tidies up consensus annotations by copying them to their most appropriate
 * new contig. Tags which entirely fit within a new contig get moved.
 * Others get copied, but we still need to trim annotations from the source
 * if they're overlapping or in some rare cases are completely in holes
 * (we had an internal hole before, but it's now on the contig end).
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int copy_contig_anno(GapIO *io, Array cmap) {
    int n = ArrayMax(cmap), i, j, k;
    contig_t *c, *new_c;
    rangec_t *rc;
    Bitmap hole;

    for (i = 0; i < n; i++) {
	contig_map *map = arrp(contig_map, cmap, i);
	int nr;

	c = new_c = NULL;
	rc = NULL;
	hole = NULL;

	/* Convert dest from unclipped to clipped */
	map->dest_start = 1;
	map->dest_end = map->src_end - map->src_start + 1;
	
	//printf("Mapped %"PRIrec" @ %d..%d -> %"PRIrec" @ %d..%d\n",
	//       map->src,  map->src_start,  map->src_end,
	//       map->dest, map->dest_start, map->dest_end);

	/* Find annotations spanning this source region */
	c = cache_search(io, GT_Contig, map->src);
	if (NULL == c) goto fail;
	cache_incr(io, c);
	rc = contig_anno_in_range(io, &c, map->src_start, map->src_end,
				  0, &nr);
	if (NULL == rc) goto fail;

	/* Trim to only consensus tags */
	for (k = j = 0; j < nr; j++) {
	    if (rc[j].flags & GRANGE_FLAG_TAG_SEQ)
		continue;

	    /*
	     * This looks odd, but it's here because a straight
	     * "rc[k++] = rc[j]" is not technically legal C and neither
	     * does Valgrind like it.
	     *
	     * The problem arises that structure assignment may be
	     * implemented using memcpy (c9x draft 6.2.6.1) and memcpy
	     * has undefined behaviour when given overlapping objects
	     * (section 7.21.2.1). Hence an asignment of rc[0] = rc[0]
	     * generates a memcpy with src==dest which is undefined.
	     * Grrr.
	     */
	    if (k != j)
		rc[k++] = rc[j];
	    else
		k++;
	}
	nr = k;

	if (!nr) {
	    free(rc);
	    cache_decr(io, c);
	    continue;
	}

	/* Produce a bitmap of hole or no-hole for source contig */
	hole = contig_hole_bitmap(io, c->rec, map->src_start, map->src_end);
	if (NULL == hole) goto fail;

	new_c = cache_search(io, GT_Contig, map->dest);
	if (NULL == new_c) goto fail;
	cache_incr(io, new_c);

	/* Duplicate them onto the destination contig */
	for (j = 0; j < nr; j++) {
	    range_t new_r;
	    rangec_t *r = &rc[j];
	    anno_ele_t *a;
	    bin_index_t *bin;
	    int in_hole = 1;

	    if (r->start < map->src_start ||
		r->end   > map->src_end) {
		in_hole = 0;
	    } else {
		for (k = r->start; k <= r->end; k++) {
		    if (BIT_CHK(hole, k - map->src_start) == 0) {
			in_hole = 0;
			break;
		    }
		}
	    }

	    //if (in_hole) {
	    //	printf("  Mov tag %"PRIrec" pos %d..%d\n",
	    //	       r->rec, r->start, r->end);
	    //} else {
	    //	printf("  Dup tag %"PRIrec" pos %d..%d\n",
	    //	       r->rec, r->start, r->end);
	    //}

	    a = cache_search(io, GT_AnnoEle, r->rec);
	    if (NULL == a) goto fail;

	    new_r.start    = r->start - map->src_start + map->dest_start;
	    new_r.end      = r->end   - map->src_start + map->dest_start;
	    new_r.flags    = r->flags;
	    new_r.mqual    = r->mqual;
	    new_r.pair_rec = 0;

	    if (in_hole) {
		new_r.rec  = a->rec;
		bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	    } else {
		new_r.rec  = anno_ele_new(io, 0, GT_Contig, 0, 0,
					  r->mqual, a->direction, a->comment);
		if (new_r.rec < 0) goto fail;
	    }

	    if (new_r.start < map->dest_start)
		new_r.start = map->dest_start;
	    if (new_r.end   > map->dest_end)
		new_r.end   = map->dest_end;

	    bin = bin_add_range(io, &new_c, &new_r, NULL, NULL, 0);
	    if (NULL == bin) goto fail;

	    if (NULL == (a = cache_search(io, GT_AnnoEle, new_r.rec)))
		goto fail;
	    if (NULL == (a = cache_rw(io, a))) goto fail;
	    a->bin = bin->rec;
	}

	cache_decr(io, c);
	cache_decr(io, new_c);
	free(rc);
	BitmapDestroy(hole);
    }

    return 0;
    
 fail:
    if (c)     cache_decr(io, c);
    if (new_c) cache_decr(io, new_c);
    if (rc)    free(rc);
    if (hole)  BitmapDestroy(hole);
    return -1;
}

/*
 * Ensure contig extents are consistent.
 *
 * Identify minimum and maximum coordinates per contig. If these
 * match the known start and end of that contig then we need to 
 * recompute the contig extents.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int fix_contigs(GapIO *io, r_pos_t *pos, int nreads) {
    int *ends, i;
    HacheTable *c_hash = HacheTableCreate(256, HASH_DYNAMIC_SIZE |
					  HASH_NONVOLATILE_KEYS);
    HacheIter *iter;
    HacheItem *hi;

    /* Accumulate start/end values per contig */
    for (i = 0 ; i < nreads; i++) {
	hi = HacheTableQuery(c_hash, (char *)&pos[i].contig,
			     sizeof(pos[i].contig));
	if (hi) {
	    ends = (int *)hi->data.p;
	    if (ends[0] > pos[i].start)
		ends[0] = pos[i].start;
	    if (ends[1] < pos[i].end)
		ends[1] = pos[i].end;
	} else {
	    ends = (int *)malloc(2 * sizeof(int));
	    ends[0] = pos[i].start;
	    ends[1] = pos[i].end;
	}
    }

    /* Iterate through contigs to check if extents need fixing */
    iter = HacheTableIterCreate();
    while (NULL != (hi = HacheTableIterNext(c_hash, iter))) {
	tg_rec crec = *(tg_rec *)hi->key;
	contig_t *c = cache_search(io, GT_Contig, crec);
	int new_start, *ns;
	int new_end, *ne;

	ends = (int *)hi->data.p;
	    
	ns = ends[0] <= c->start ? &new_start : NULL;
	ne = ends[1] >= c->end   ? &new_end   : NULL;
	if (ns || ne) {
	    if (-1 != consensus_unclipped_range(io, c->rec, ns, ne)) {
		c = cache_rw(io, c);
		if (ns) c->start = *ns;
		if (ne) c->end   = *ne;
		if (ns || ne)
		    c->timestamp = io_timestamp_incr(io);
	    }
	}
    }
    HacheTableIterDestroy(iter);

    HacheTableDestroy(c_hash, 0);

    return 0;
}

/*
 * Find annotations attached to the reads in rnums[].
 *
 * The easy way to do this is to find the positions for each read (already
 * supplied in pos[i]), then query all annotations over that region, and
 * finally then search those annotations for ones attached to the correct
 * record.
 *
 * This is slow when nreads is large and if our reads are all from the same
 * region we end up repeating the same searches many times and iterating
 * through the same list of annotations many times.
 *
 * So the new algorithm is:
 * 1) cluster reads into overlapping sets
 * 2) per set: identify all annotations in that set range
 * 3) loop through annotation in the set assigning to their appropriate
 *    reads.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int anno_range(GapIO *io, tg_rec contig, int start, int end,
		      r_pos_t *pos, int rfrom, int rto) {
    rangec_t *anno;
    int n_anno;
    contig_t *c;
    HacheTable *rmap;
    HacheData hd;
    HacheItem *hi;
    int i, j;

    if (NULL == (c = cache_search(io, GT_Contig, contig)))
	return -1;

    /* Build a hash per sequence holding the count of annotations */
    if (NULL == (rmap = HacheTableCreate(rto-rfrom, HASH_DYNAMIC_SIZE |
					 HASH_NONVOLATILE_KEYS)))
	return -1;
    for (i = rfrom; i <= rto; i++) {
	hd.i = i;
	HacheTableAdd(rmap, (char *)&pos[i].rec, sizeof(tg_rec), hd, NULL);
    }

    /* Find annotations in this region */
    anno = contig_anno_in_range(io, &c, start, end, 0, &n_anno);
    if (NULL == anno) goto fail;

    /* Copy to the pos[] array */
    for (i = 0; i < n_anno; i++) {
	rangec_t *a;
	if (!(hi = HacheTableQuery(rmap, (char *)&anno[i].pair_rec,
				   sizeof(anno[i].pair_rec))))
	    continue;
	
	/* Could be more efficient, but not likely the slow bit now */
	j = hi->data.i;
	pos[j].n_anno++;
	a = realloc(pos[j].anno, pos[j].n_anno * sizeof(rangec_t));
	if (NULL == a) goto fail;
	pos[j].anno = a;
	pos[j].anno[pos[j].n_anno-1] = anno[i];
    }

    free(anno);

    HacheTableDestroy(rmap, 0);
    return 0;

 fail:
    if (NULL != rmap) HacheTableDestroy(rmap, 0);
    if (NULL != anno) free(anno);
    return -1;
}

static int find_annos(GapIO *io, r_pos_t *pos, int nreads) {
    int i, j;
    tg_rec contig;
    int start, end;

    if (nreads == 0)
	return 0;

    contig = pos[0].contig;
    start  = pos[0].start;
    end    = pos[0].end;

    for (j = 0, i = 1; i < nreads; i++) {
	if (pos[i].start > end || pos[i].contig != contig) {
	    anno_range(io, contig, start, end, pos, j, i-1);
	    j = i;
	    contig = pos[i].contig;
	    start  = pos[i].start;
	    end    = pos[i].end;
	} else if (pos[i].end > end) {
	    end = pos[i].end;
	}
    }
    anno_range(io, contig, start, end, pos, j, i-1);

    return 0;
}

/*-----------------------------------------------------------------------------
 * External interfaces
 */


/**
 * Removes a set of readings from either the contig or the database.
 *
 * When removing from the database we need to delete everything related to
 * that reading (annotations, template if not used elsewhere, etc). When
 * moving to a new contig, we prefer to keep all readings clustered together.
 *
 * Ie if we remove A & B (overlapping) from one contig and C from another
 * then we create two new contigs containing A & B in one and C in the other.
 *
 * When creating new contigs, we have the option of copying over any
 * overlapping consensus tags to the new contigs. This choice only refers to
 * consensus tags; reading tags are always copied.
 *
 * move == 0   => remove
 * move == 1   => split to new single-read contigs
 * move == 2   => move to new still-joined contigs
 *
 * Returns 0 on success
 *        -1 on failure
 */
int disassemble_readings(GapIO *io, tg_rec *rnums, int nreads, int move,
                         int remove_holes, int duplicate_tags)
{
    int i,err = 0;
    r_pos_t *pos;
    HacheTable *dup_hash, *bin_hash;
    HacheIter *iter;
    HacheItem *hi;
    Array cmap;

    vfuncheader("Disassemble_readings");
    //    check_contig_bins(io);

    if (nreads <= 0)
	return 0;

    /*
     * The plan:
     *
     * 1. Produce a table of which contig number and region each reading
     *    record is within.
     *
     * 2. Remove these readings from their contigs.
     * 2b. Ensure if it's a read-pair that the other end gets the
     *     pairing information and flags cleared.
     *
     * 3a. move==0: deallocate the reading record.
     * 3b. move==1: produce new contigs/bins for each read.
     * 3c. move==2: produce 1 new contig for each overlapping set of
     *              disassembled reads.
     *
     * 4. Move/copy any consensus annotations.
     *
     * 5. Check for holes within the "source" contigs. We use the
     *    initial table for this, along with regions to look around.
     *
     */

    if (NULL == (pos = calloc(nreads, sizeof(*pos))))
	return -1;

    /*
     * Remove from contig bin.
     * Also handles accumulation of contig/positions and removal of
     * duplicate reads.
     */
    dup_hash = HacheTableCreate(8192, HASH_DYNAMIC_SIZE |
				      HASH_NONVOLATILE_KEYS);
    bin_hash = HacheTableCreate(256, HASH_DYNAMIC_SIZE |
				HASH_NONVOLATILE_KEYS);
    cmap = ArrayCreate(sizeof(contig_map), 0);

    //puts("Part 1 - unlink read"); system("date");

    for (i = 0; i < nreads; i++) {
	HacheData hd;
	int n;

	/* Part 1, 2 and 2b: remove reads */
	hd.i = 0;
	HacheTableAdd(dup_hash, (char *)&rnums[i], sizeof(rnums[i]), hd, &n);
	if (!n) {
	    gio_debug(io, 1, "Skipping duplicate entry %"PRIrec"\n", rnums[i]);
	    pos[i].contig = 0;
	    continue;
	}
	if (rnums[i] == -1) {
	    continue;
	}

	/* This doesn't tidy up the bin, so remember those to fix later */
	if (unlink_read(io, rnums[i], &pos[i], move == 0)) {
	    verror(ERR_WARN, "disassemble_readings",
		   "Failed to unlink seq #%"PRIrec, rnums[i]);
	    //return -1;
	    rnums[i] = -1;
	    pos[i].contig = 0;
	    continue;
	}
	HacheTableAdd(bin_hash, (char *)&pos[i].bin, sizeof(pos[i].bin),
		      hd, NULL);

	pos[i].rec = rnums[i];
    }

    /* Ensure bins are consistent */
    iter = HacheTableIterCreate();
    while (NULL != (hi = HacheTableIterNext(bin_hash, iter))) {
	tg_rec brec = *(tg_rec *)hi->key;
	bin_index_t *bin = cache_search(io, GT_Bin, brec);
	bin_set_used_range(io, bin);
    }
    HacheTableIterDestroy(iter);

    HacheTableDestroy(dup_hash, 0);
    HacheTableDestroy(bin_hash, 0);


    //puts("Part 2 - Sort by pos"); system("date");

    /* Sort position table and drop the duplicate entries */
    xio = io;
    qsort(pos, nreads, sizeof(*pos), pos_sort);
    for (i = 0; i < nreads; i++)
	if (pos[i].contig)
	    break;
    if (i != 0) {
	memmove(&pos[0], &pos[i], (nreads-i) * sizeof(pos[0]));
	nreads -= i;
    }


    //puts("Part 2.1 - Find annos"); system("date");

    /* Identify annotations attached to the reads */
    find_annos(io, pos, nreads);


    //puts("Part 2.2 - Fix contig extents"); system("date");

    /* Ensure contig extents are consistent */
    fix_contigs(io, pos, nreads);


    //puts("Part 3 - Del or move"); system("date");

    /* Part 3: */
    switch (move) {
    case 0:
	/* 3a. Deallocate the record */
	for (i = 0; i < nreads; i++) {
	    if (seq_deallocate(io, &pos[i]))
		err = 1;
	}
	break;

    case 1:
	/* 3b. produce new contigs/bins for each read. */
	for (i = 0; i < nreads; i++) {
	    if (create_contig_from(io, &pos[i], 1, cmap))
		err = 1;
	}
	break;

    case 2:
	/* 3c. produce 1 new contig for each overlapping set of
	       disassembled reads. */
	if (move_reads(io, pos, nreads, cmap))
	    err = 1;
	break;

    default:
	fprintf(stderr, "Unexpected 'move' value %d\n", move);
	err = 1;
	return -1;
    }

    cache_flush(io);

    //puts("Part 4 - Copy annos"); system("date");

    /* Part 4: Move/copy any consensus annotations. */
    copy_contig_anno(io, cmap);

    /* Part 5: fix holes in source contigs */
    if (fix_holes(io, pos, nreads, remove_holes))
	/* Too late to undo, so keep going! */
	err = 1;

    cache_flush(io);

    //puts("Part 5 - Deallocate"); system("date");

    // check_contig_bins(io);

    for (i = 0; i < nreads; i++)
	if (pos[i].anno)
	    free(pos[i].anno);
    free(pos);

    ArrayDestroy(cmap);

    return err ? -1 : 0;
}

typedef struct {
    tg_rec r,p; /* read and pair */
} rec_rec;
static int tg_rec_cmp(const void *v1, const void *v2) {
    return ((const rec_rec *)v1)->r - ((const rec_rec *)v2)->r;
}

/*
 * As per disassemble readings, but removes entire contigs.
 *
 * This is substantially faster as it doesn't need to track a lot of the
 * changes to contig dimensions.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int disassemble_contigs(GapIO *io, tg_rec *cnums, int ncontigs) {
    int i;
    int ret = 0;
    int npairs = 0;
    HashTable *pairs;
    HashIter *iter;
    HashItem *hi;
    rec_rec *seqs;

    pairs = HashTableCreate(8192, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);

    for (i = 0; i < ncontigs; i++) {
	contig_t *c;
	contig_iterator *iter;
	rangec_t *r;

	vmessage("Processing contig %d of %d\n", i+1, ncontigs);
	UpdateTextOutput();

	iter = contig_iter_new_by_type(io, cnums[i], 1, CITER_FIRST,
				       CITER_CSTART, CITER_CEND,
				       GRANGE_FLAG_ISANY);

	if (!iter) {
	    verror(ERR_WARN, "disassemble_contigs",
		   "Failed to load contig #%"PRIrec, cnums[i]);
	    ret = 1;
	    continue;
	}

	/* Destroy contents of the contig */
	while (NULL != (r = contig_iter_next(io, iter))) {
	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;
	    
	    switch (r->flags & GRANGE_FLAG_ISMASK) {
	    case GRANGE_FLAG_ISSEQ: {
		seq_t *s = cache_search(io, GT_Seq, r->rec);
		tg_rec root;

		if (!s) {
		    ret = 1;
		    continue;
		}

		/* Remove from B+Tree */
		root = io->iface->seq.index_del(io->dbh, s->name, s->rec);
		if (root != -1 && root != io->db->seq_name_index) {
		    io->db = cache_rw(io, io->db);
		    io->db->seq_name_index = root;
		}

		/* Identify pairs */
		if (r->pair_rec) {
		    hi = HashTableSearch(pairs, (char *)&r->rec,
					 sizeof(tg_rec));
		    if (hi) {
			/* Other end already removed */
			HashTableDel(pairs, hi, 0);
			npairs--;
		    } else {
			/* Mark for update of other end */
			HashData hd;
			hd.i = r->rec;
			HashTableAdd(pairs, (char *)&r->pair_rec,
				     sizeof(tg_rec), hd, NULL);
			npairs++;
		    }
		}

		/* Remove from seq_block */
		cache_item_remove(io, GT_Seq, r->rec);

		break;
	    }

	    case GRANGE_FLAG_ISANNO:
		cache_item_remove(io, GT_AnnoEle, r->rec);
		break;
	    }
	}

	contig_iter_del(iter);

	/* Destroy the contig bin structure itself */
	c = cache_search(io, GT_Contig, cnums[i]);

	if (c && c->bin)
	    bin_destroy_recurse(io, c->bin);
	contig_destroy(io, cnums[i]);
	
	/* Reduce btree memory usage by flushing after each contig */
	cache_flush(io);
    }
    vmessage("Flushing deletions\n");
    UpdateTextOutput();
    cache_flush(io);

    /*
     * Sort pairs by sequence record number. This ensures that if we have
     * a very large amount of sequences to update then we do so in a cache
     * sensitive manner.
     */
    if (NULL == (seqs = (rec_rec *)xmalloc(npairs * sizeof(*seqs))))
	return -1;

    iter = HashTableIterCreate();
    i = 0;
    while (NULL != (hi = HashTableIterNext(pairs, iter))) {
	seqs[i].r = *(tg_rec *)hi->key;
	seqs[i++].p = hi->data.i;
    }
    assert(i == npairs);
    HashTableIterDestroy(iter);
    HashTableDestroy(pairs, 0);

    /* Sort the list */
    qsort(seqs, npairs, sizeof(*seqs), tg_rec_cmp);

    vmessage("Unlinking from read-pairs\n");
    UpdateTextOutput();
    for (i = 0; i < npairs; i++) {
	seq_t *s = cache_search(io, GT_Seq, seqs[i].r);
	bin_index_t *bin;
	range_t *r;

	if (!s)
	    continue;

	if (i % 1000 == 0) {
	    vmessage("    %d of %d\n", i, npairs);
	    UpdateTextOutput();
	    if (i % 10000 == 0)
		cache_flush(io);
	}


	if (s->parent_rec == seqs[i].p) {
	    s = cache_rw(io, s);
	    s->parent_type = 0;
	    s->parent_rec = 0;
	}

	/* Also held in the bin range too */
	bin = cache_search(io, GT_Bin, s->bin);
	if (!bin || !bin->rng)
	    continue;

	r = arrp(range_t, bin->rng, s->bin_index);
	assert(r->rec == s->rec);

	bin = cache_rw(io, bin);
	bin->flags |= BIN_RANGE_UPDATED;

	r->pair_rec = 0;
	r->flags &= ~GRANGE_FLAG_TYPE_MASK;
	r->flags |=  GRANGE_FLAG_TYPE_SINGLE;
    }

    xfree(seqs);
    cache_flush(io);

    return ret;
}
