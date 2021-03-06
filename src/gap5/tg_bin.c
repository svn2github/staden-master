#include <assert.h>
#include <string.h>
#include <math.h>

#include "xalloc.h"
#include "tg_gio.h"
#include "tg_tracks.h"
#include "consensus.h"

#define get_bin(io, bnum) ((bin_index_t *)cache_search((io), GT_Bin, (bnum)))

/*
 * Allocates a new bin record.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
#if 0
int bin_new(GapIO *io, int pos, int sz, int parent, int parent_type) {
    int rec;
    bin_index_t bin;

    /* Initialise disk struct */
    bin.pos         = pos;
    bin.size        = sz;
    bin.start_used  = 0;
    bin.end_used    = 0;
    bin.parent      = parent;
    bin.parent_type = parent_type;
    bin.child[0]    = 0;
    bin.child[1]    = 0;
    bin.rng	    = NULL;
    bin.rng_rec     = 0;
    bin.flags       = BIN_BIN_UPDATED;
    bin.track       = NULL;
    bin.track_rec   = 0;
    bin.nseqs       = 0;
    bin.rng_free    = -1;
    bin.nrefpos     = 0;
    bin.nanno       = 0;

    if (-1 == (rec = io->iface->bin.create(io->dbh, &bin)))
	return -1;

    return rec;
}
#endif

tg_rec bin_new(GapIO *io, int pos, int sz, tg_rec parent, int parent_type) {
    tg_rec rec;
    bin_index_t *bin;

    if (-1 == (rec = io->iface->bin.create(io->dbh, NULL)))
	return -1;

    /* Initialise disk struct */
    bin = get_bin(io, rec);
    bin = cache_rw(io, bin);
    bin->pos         = pos;
    bin->size        = sz;
    bin->start_used  = 0;
    bin->end_used    = 0;
    bin->parent      = parent;
    bin->parent_type = parent_type;
    bin->child[0]    = 0;
    bin->child[1]    = 0;
    bin->rng	    = NULL;
    bin->rng_rec     = 0;
    bin->flags       = BIN_BIN_UPDATED;
    bin->track       = NULL;
    bin->track_rec   = 0;
    bin->nseqs       = 0;
    bin->rng_free    = -1;
    bin->nrefpos     = 0;
    bin->nanno       = 0;

    return rec;
}


/*
 * Doubles up the number of bins by adding a new root node and duplicating
 * the tree.
 *
 * It takes the old root_id as an argument and returns the new one.
 * Returns -1 on failure.
 */
static bin_index_t *contig_extend_bins_right(GapIO *io, contig_t **c) {
    tg_rec old_root_id = contig_get_bin(c);
    bin_index_t *oroot = get_bin(io, old_root_id), *nroot;
    tg_rec root_id;
    size_t sz = oroot->size;

    cache_incr(io, oroot);
    if (!(oroot = cache_rw(io, oroot))) {
	cache_decr(io, oroot);
	return NULL;
    }

    /* Create a new root */
    sz *= 2;
    if (sz > INT_MAX)
	sz = INT_MAX;
    root_id = bin_new(io, oroot->pos, sz, oroot->parent, oroot->parent_type);
    nroot = get_bin(io, root_id);
    cache_incr(io, nroot);
    if (!(nroot = cache_rw(io, nroot))) {
	cache_decr(io, oroot);
	cache_decr(io, nroot);
	return NULL;
    }

    nroot->child[0] = old_root_id;
    nroot->child[1] = 0;
    nroot->nseqs    = oroot->nseqs;
    nroot->nrefpos  = oroot->nrefpos;
    nroot->nanno    = oroot->nanno;

    nroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, nroot);

    /* Move old left bin */
    oroot->parent = root_id;
    oroot->parent_type = GT_Bin;
    oroot->pos = 0;

    oroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, oroot);

    contig_set_bin(io, c, root_id);

    return nroot;
}

static bin_index_t *contig_extend_bins_left(GapIO *io, contig_t **c) {
    tg_rec old_root_id = contig_get_bin(c);
    bin_index_t *oroot = get_bin(io, old_root_id), *nroot;
    tg_rec root_id;
    int sz = oroot->size;

    cache_incr(io, oroot);
    if (!(oroot = cache_rw(io, oroot))) {
	cache_decr(io, oroot);
	return NULL;
    }

    /* Create a new root */
    root_id = bin_new(io, oroot->pos-sz, sz*2, oroot->parent, oroot->parent_type);
    nroot = get_bin(io, root_id);
    cache_incr(io, nroot);
    if (!(nroot = cache_rw(io, nroot))) {
	cache_decr(io, oroot);
	cache_decr(io, nroot);
	return NULL;
    }

    nroot->child[0] = 0;
    nroot->child[1] = old_root_id;
    nroot->nseqs    = oroot->nseqs;
    nroot->nrefpos  = oroot->nrefpos;
    nroot->nanno    = oroot->nanno;

    nroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, nroot);

    /* Move old right bin */
    oroot->parent = root_id;
    oroot->parent_type = GT_Bin;
    oroot->pos = sz;

    oroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, oroot);

    contig_set_bin(io, c, root_id);

    return nroot;
}

/*
 * Allocates and returns next range from a range array
 * Returns -1 for error.
 */
static int next_range(GapIO *io, bin_index_t *bin) {
    Array ra = bin->rng;

    bin = cache_rw(io, bin);
    bin->flags |= BIN_BIN_UPDATED | BIN_RANGE_UPDATED;

    if (bin->rng_free == -1) {
	if (NULL == ArrayRef(ra, ArrayMax(ra)))
	    return -1;

	return ArrayMax(ra)-1;

    } else {
	int tmp;
	range_t *r = arrp(range_t, bin->rng, bin->rng_free);

	assert(bin->rng_free < ArrayMax(ra));
	assert(r->flags & GRANGE_FLAG_UNUSED);
	
	tmp = bin->rng_free;
	bin->rng_free = (int)r->rec;
	return tmp;
    }
}

/*
 * This finds a bin suitable for a given range. If such a bin doesn't
 * exist it can optionally create one if the 'extend' flag is true.
 *
 * When a bin is found the absolute offset of that bin is returned
 * in 'offset_r' (may be NULL).
 *
 * Returns the bin pointer on success
 *         NULL on failure
 */
#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

//#define CACHE_LAST_BIN
bin_index_t *bin_for_range(GapIO *io, contig_t **c,
			   int start, int end, int extend,
			   int *offset_r,  int *comp_r) {
    int offset;
    bin_index_t *bin = get_bin(io, contig_get_bin(c));
    int complement = 0;
    int f_a, f_b;

#ifdef CACHE_LAST_BIN
    static tg_rec last_c = 0;
    static GapIO *last_io = NULL;
    bin_index_t *last_bin = NULL;
    static tg_rec last_bin_rec = 0;
    static int last_start = 0, last_end = 0;
    static int last_offset, last_complement;
#endif

    if (NULL == bin) return NULL;
    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    offset = bin->pos;
    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    //cache_incr(io, bin);

    /*
     * If we're trying to insert a new node beyond the total bounds of the
     * root node then we need to extend the bin structure first to either
     * the left and/or the right.
     */
    while (end >= bin->pos + bin->size) {
	//cache_decr(io, bin);
	if (extend) {
	    if (NULL == (bin = contig_extend_bins_right(io, c))) return NULL;
	}
	else
	    return NULL;
	//cache_incr(io, bin);

	complement = 0;
	/*
	 * But root has switched from complemented to uncomplemented, so the
	 * range may need fixing too.
	 * See tg_tcl.c for an example way of detecting and fixing this.
	 */
    }

    while (start < bin->pos) {
	//cache_decr(io, bin);
	if (extend) {
	    if (NULL == (bin = contig_extend_bins_left(io, c))) return NULL;
	}
	else
	    return NULL;
	//cache_incr(io, bin);

	complement = 0;
    }

    /*
     * In theory we can jump straight to a candidate starting bin, possibly
     * even returning it right here if it's the min bin size, saving about
     * 10% of our CPU time in this function
     */
#ifdef CACHE_LAST_BIN
    if (last_bin_rec && last_c == (*c)->rec && last_io == io) {
	last_bin = cache_search(io, GT_Bin, last_bin_rec);
	if (start >= last_start && end <= last_end) {
	    if (last_bin && last_bin->size <= io->min_bin_size &&
		!last_bin->child[0] && !last_bin->child[1]) {
		/* leaf node, so we can return right now */
		if (offset_r)
		    *offset_r = last_offset;
		if (comp_r)
		    *comp_r = last_complement;
		return last_bin;
	    }

	    /* Maybe a smaller bin, but start the search from here on */
	    bin = last_bin;
	    offset = last_offset;
	    complement = last_complement;
	    //cache_incr(io, bin);

	    if (complement) {
		f_a = -1;
		f_b = offset + bin->size-1;
	    } else {
		f_a = +1;
		f_b = offset;
	    }

	    goto jump;
	}
    }
#endif

    /* Now recurse down the bin hierachy searching for the smallest bin */
    offset = bin->pos;
#ifdef CACHE_LAST_BIN
 jump:
#endif
    cache_incr(io, bin);

    for (;;) {
	int i;
	bin_index_t *ch;

	if (complement) {
	    f_a = -1;
	    f_b = offset + bin->size-1;
	} else {
	    f_a = +1;
	    f_b = offset;
	}

	/* Find which child bin is most suitable */
	for (i = 0; i < 2;) {
	    if (bin->child[i] <= 0) {
		i++;
		continue;
	    }

	    if (NULL == (ch = get_bin(io, bin->child[i]))) goto error;

	    //	    if (start >= offset + ch->pos &&
	    //		end   <= offset + ch->pos + ch->size-1) {
	    if (start >= NMIN(ch->pos, ch->pos + ch->size-1) &&
		end   <= NMAX(ch->pos, ch->pos + ch->size-1)) {
		cache_decr(io, bin);
		bin = ch;
		cache_incr(io, ch);
		offset = NMIN(ch->pos, ch->pos + ch->size-1);

		if (bin->flags & BIN_COMPLEMENTED) {
		    complement ^= 1;
		}
		if (complement) {
		    f_a = -1;
		    f_b = offset + bin->size-1;
		} else {
		    f_a = +1;
		    f_b = offset;
		}

		i = 0; /* restart loop */
	    } else {
		i++;
	    }
	}

	if (!extend) {
	    if (offset_r)
		*offset_r = offset;
	    if (comp_r)
		*comp_r = complement;
	    cache_decr(io, bin);
	    return bin;
	}

	/*
	 * We now have the smallest bin available holding this range, but 
	 * as we delay creation of the sub-bins until required then we
	 * should perhaps create a child bin.
	 */
	if (bin->size > io->min_bin_size &&
	    (!bin->child[0] || !bin->child[1])) {
	    int free_slt = (!bin->child[0] ? 1 : 0) | (!bin->child[1] ? 2 : 0);
	    int pos  = 0, sz = 0, slot = -1;
	    switch (free_slt) {
	    case 1:
	    case 2:
		if (NULL == (ch = get_bin(io, bin->child[2 - free_slt])))
		    goto error;

		pos = ch->pos == 0 ? ch->size : 0;
		sz  = ch->pos == 0 ? bin->size - pos : ch->pos;
		slot = free_slt - 1;
		break;
	    case 3:
		sz = bin->size / 2;
		if (start >= NMIN(0, sz - 1) && end <= NMAX(0, sz - 1)) {
		    pos = 0;
		    slot = complement ? 1 : 0;
		} else {
		    pos = sz;
		    sz = bin->size - pos;
		    slot = complement ? 0 : 1;
		}
		break;
	    }
	    if (slot >= 0 &&
		start >= NMIN(pos, pos+sz-1) &&
		end   <= NMAX(pos, pos+sz-1)) {
		/* It would fit - create it and continue recursing */
		bin_index_t *binw = cache_rw(io, bin);
		if (NULL == binw) goto error;
		
		binw->child[slot] = bin_new(io, pos, sz, binw->rec, GT_Bin);
		if (0 == binw->child[slot]) goto error;
		binw->flags |= BIN_BIN_UPDATED;
		cache_decr(io, binw);		
		if (NULL == (bin = get_bin(io, binw->child[slot]))) goto error;
		cache_incr(io, bin);
		offset = NMIN(pos, pos+sz-1);
		continue;
	    }
	}

	/* At the smallest bin already */
	break;
    }

    if (offset_r)
	*offset_r = offset;
    if (comp_r)
	*comp_r = complement;

#ifdef CACHE_LAST_BIN
    last_io = io;
    last_c = (*c)->rec;
    last_bin_rec = bin->rec;
    last_start = offset;
    last_end = offset + bin->size-1;
    last_offset = offset;
    last_complement = complement;
#endif

    cache_decr(io, bin);
    return bin;

 error:
    cache_decr(io, bin);
    return NULL;
}

/*
 * Adds 'n' to the nseq counter for a bin and all parent bins chaining up
 * to the root node. 'n' may be negative too.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_incr_nseq(GapIO *io, bin_index_t *bin, int n) {
    contig_t *c;

    while (bin) {
	if (!(bin = cache_rw(io, bin)))
	    return -1;
	bin->nseqs += n;
	bin->flags |= BIN_BIN_UPDATED;
	bin->flags &= ~BIN_CONS_VALID;

	if (bin->parent_type != GT_Bin)
	    break;

	assert(bin->rec != bin->parent);

	bin = get_bin(io, bin->parent);
    }

    assert(bin->parent_type == GT_Contig);
    c = cache_search(io, GT_Contig, bin->parent);
    c = cache_rw(io, c);
    c->nseqs += n;

    return 0;
}


/*
 * Adds 'n' to the nrefpos counter for a bin and all parent bins chaining up
 * to the root node. 'n' may be negative too.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_incr_nrefpos(GapIO *io, bin_index_t *bin, int n) {
    contig_t *c;

    while (bin) {
	if (!(bin = cache_rw(io, bin)))
	    return -1;
	bin->nrefpos += n;
	bin->flags |= BIN_BIN_UPDATED;
	bin->flags &= ~BIN_CONS_VALID;

	if (bin->parent_type != GT_Bin)
	    break;

	assert(bin->rec != bin->parent);

	bin = get_bin(io, bin->parent);
    }

    assert(bin->parent_type == GT_Contig);
    c = cache_search(io, GT_Contig, bin->parent);
    c = cache_rw(io, c);
    c->nrefpos += n;

    return 0;
}


/*
 * Adds 'n' to the nanno counter for a bin and all parent bins chaining up
 * to the root node. 'n' may be negative too.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_incr_nanno(GapIO *io, bin_index_t *bin, int n) {
    contig_t *c;

    while (bin) {
	if (!(bin = cache_rw(io, bin)))
	    return -1;
	bin->nanno += n;
	bin->flags |= BIN_BIN_UPDATED;
	bin->flags &= ~BIN_CONS_VALID;

	if (bin->parent_type != GT_Bin)
	    break;

	assert(bin->rec != bin->parent);

	bin = get_bin(io, bin->parent);
    }

    assert(bin->parent_type == GT_Contig);
    c = cache_search(io, GT_Contig, bin->parent);
    c = cache_rw(io, c);
    c->nanno += n;

    return 0;
}


/*
 * Returns whether a bin is empty (has no ranges). This isn't as trivial
 * as it sounds as we may have bin->rng allocated and containing data, but
 * with all data elements being UNUSED.
 *
 * Note we don't indicate whether the bin hierarchy is empty (and an empty
 * bin may have children with data). We're just returning data about this
 * specific bin range array.
 *
 * Returns 1 for empty
 *         0 if not
 */
int bin_empty(bin_index_t *bin) {
    int i;

    if (bin->rng == NULL || ArrayMax(bin->rng) == 0)
	return 1;

    // if (bin->start_used || bin->end_used)
    // 	return 0;

    /* So we have a non-blank range array, but maybe it's all unused? */
    for (i = 0; i < ArrayMax(bin->rng); i++) {
        range_t *r = arrp(range_t, bin->rng, i);
	if (!(r->flags & GRANGE_FLAG_UNUSED))
	    return 0;
    }

    return 1;
}


/*
 * Adds a range to the contig.
 *
 * Delay_nseq controls whether we update bin_incr_nseq() immediately or
 * whether to delay until later.
 * A value of 0 will update nseq on each call.
 * A value of 1 will update nseq whenever we happen to write to a different
 * bin than before.
 * A value of -1 will update any pending nseqs and do no other action (so
 * c, r, r_out, etc are all ignored and NULL is returned).
 *
 * Returns the bin we added the range to on success
 *         NULL on failure
 */

/*
 * As per bin_add_range() but add the range to a specific bin.
 * We use this for annotations to keep them in the same bin as their
 * associated sequence.
 *
 * We share the same last_bin and incr_*value variables as used by
 * bin_add_range() above.
 *
 * Returns the bin we added the range to on success
 *         NULL on failure
 */
bin_index_t *bin_add_to_range(GapIO *io, contig_t **c, tg_rec brec, range_t *r,
			      range_t **r_out, int *complemented,
			      int delay_nseq) {

    bin_index_t *bin;
    range_t *r2;
    int nr, offset, comp;

    /* Tidy-up operation when adding ranges in bulk */
    if (delay_nseq == -1) {
	if (io->last_bin
	    && (io->incr_svalue || io->incr_rvalue || io->incr_avalue)) {
	    if (c) *c = cache_rw(io, *c);
	    bin = cache_search(io, GT_Bin, io->last_bin);
	    bin_incr_nseq(io, bin, io->incr_svalue);
	    bin_incr_nrefpos(io, bin, io->incr_rvalue);
	    bin_incr_nanno(io, bin, io->incr_avalue);
	}
	io->last_bin = 0;
	io->incr_svalue = io->incr_rvalue = io->incr_avalue = 0;

	return NULL;
    }

    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	if (contig_get_start(c) == contig_get_end(c)
	    && contig_get_start(c) == 0) {
	    contig_set_start(io, c, r->start);
	    contig_set_end(io, c, r->end);
	    (*c)->clipped_timestamp = 0;
	}
	if (contig_get_start(c) > r->start) {
	    contig_set_start(io, c, r->start);
	    (*c)->clipped_timestamp = 0;
	}

	if (contig_get_end(c) < r->end) {
	    contig_set_end(io, c, r->end);
	    (*c)->clipped_timestamp = 0;
	}
	
	/* Check if the sequence may have changed visible start/end even if
	   it was within the boundaries of used start/end */
	if ((*c)->clipped_timestamp == (*c)->timestamp
	    && (r->start < (*c)->clipped_start || r->end > (*c)->clipped_end)) {
	    (*c)->clipped_timestamp = 0;
	}
    }

    if (brec) {
	tg_rec ctg;

	if (!(bin = cache_search(io, GT_Bin, brec)))
	    return NULL;
	cache_incr(io, bin); /* prevent bin_get_pos() purging it */

	bin_get_position(io, bin, &ctg, &offset, &comp);
	if (r->start < offset || r->end > offset + bin->size) {
	    fprintf(stderr, "Range will not fit inside requested bin\n");
	}
    } else {
	if (!(bin = bin_for_range(io, c, r->start, r->end, 1, &offset,
				  &comp))) //complemented)))
	    return NULL;
    }

    if (complemented)
	*complemented = comp;

    if (!(bin = cache_rw(io, bin)))
	return NULL;
    if (brec)
	cache_decr(io, bin);

    /* Adjust start/end used in bin */
    if (!bin_empty(bin)) {
	if (comp) {
	    if (bin->start_used > bin->size-1 - (r->end - offset))
		bin->start_used = bin->size-1 - (r->end - offset);
	    if (bin->end_used   < bin->size-1 - (r->start - offset))
		bin->end_used   = bin->size-1 - (r->start - offset);
	} else {
	    if (bin->start_used > r->start - offset)
		bin->start_used = r->start - offset;
	    if (bin->end_used   < r->end - offset)
		bin->end_used   = r->end - offset;
	}
    } else {
	/* Initial case */
	if (comp) {
	    bin->start_used = bin->size-1 - (r->end - offset);
	    bin->end_used   = bin->size-1 - (r->start - offset);
	} else {
	    bin->start_used = r->start - offset;
	    bin->end_used   = r->end - offset;
	}
    }
    
    /* Update Range array */
    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
    bin->flags &= ~BIN_CONS_VALID;
    if (!bin->rng)
	bin->rng = ArrayCreate(sizeof(range_t), 0);

    nr = next_range(io, bin);
    r2 = arrp(range_t, bin->rng, nr);

    *r2 = *r; /* struct copy */
    r2->start -= offset;
    r2->end -= offset;

    if (comp) {
	int tmp = r2->start;
	r2->start = bin->size-1 - r2->end;
	r2->end   = bin->size-1 - tmp;
    }

    if (r_out)
	*r_out = r2;

    /* Update nseq in bins, delaying this to avoid needless writes */
    if (delay_nseq == 1 && bin->rec != io->last_bin
	&& (io->incr_svalue || io->incr_rvalue || io->incr_avalue)) {
	bin_index_t *b2 = cache_search(io, GT_Bin, io->last_bin);
	if (c) *c = cache_rw(io, *c);
	bin_incr_nseq(io, b2, io->incr_svalue);
	bin_incr_nrefpos(io, b2, io->incr_rvalue);
	bin_incr_nanno(io, b2, io->incr_avalue);
	io->incr_svalue = io->incr_rvalue = io->incr_avalue = 0;
	io->last_bin = bin->rec;
    }

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	if (delay_nseq == 0) {
	    *c = cache_rw(io, *c);
	    bin_incr_nseq(io, bin, 1);
	} else {
	    io->incr_svalue++;
	    io->last_bin = bin->rec;
	}
    }

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
	if (delay_nseq == 0) {
	    *c = cache_rw(io, *c);
	    bin_incr_nrefpos(io, bin, 1);
	} else {
	    io->incr_rvalue++;
	    io->last_bin = bin->rec;
	}
    }

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	if (delay_nseq == 0) {
	    *c = cache_rw(io, *c);
	    bin_incr_nanno(io, bin, 1);
	} else {
	    io->incr_avalue++;
	    io->last_bin = bin->rec;
	}
    }

    return bin;
}

/*
 * Original version. As per bin_add_to_range but can add to any bin instead
 * of a specific one.
 */
bin_index_t *bin_add_range(GapIO *io, contig_t **c, range_t *r,
			   range_t **r_out, int *complemented,
			   int delay_nseq) {
    return bin_add_to_range(io, c, 0, r, r_out, complemented, delay_nseq);
}

/*
 * Finds the contig number and position of a record number.
 *
 * If non-NULL r_out is filled with the associated range_t struct.
 *
 * If non-NULL i_out is filled with a pointer to the object referred to
 * by type/rec. (This is just for minor optimisations.) If returned it
 * will have had cache_incr() run on it, so the caller should
 * use cache_decr() to permit deallocation.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_get_item_position(GapIO *io, int type, tg_rec rec,
			  tg_rec *contig, int *start, int *end, int *orient,
			  tg_rec *brec, range_t *r_out, void **i_out) {
    bin_index_t *bin;
    tg_rec bnum;
    int i, offset1 = 0, offset2 = 0, found = 0;
    int comp = 0;
    int idx = -1;
    int orig_start, orig_end;
    tg_rec orig_bin;

    if (type == GT_AnnoEle) {
	anno_ele_t *a = cache_search(io, GT_AnnoEle, rec);
	if (!a)
	    return -1;

	if (i_out) {
	    cache_incr(io, a);
	    *i_out = a;
	}
	bnum = a->bin;
    } else if (type == GT_Seq) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);
	if (!s)
	    return -1;

	if (i_out) {
	    cache_incr(io, s);
	    *i_out = s;
	}
	bnum = s->bin;
	idx = s->bin_index;
    } else {
	fprintf(stderr, "Unsupported record type %d in bin_get_item_position\n",
		type);
	return -1;
    }
    
    if (brec)
	*brec = bnum;

    if (!bnum) goto fail;

    /* Find the position of this anno within the bin */
    bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
    if (NULL == bin) goto fail;
    if (idx != -1 && bin->rng && idx < ArrayMax(bin->rng)) {
	/* Check it's valid */
	range_t *r = arrp(range_t, bin->rng, idx);
	if (r->rec == rec) {
	    found = 1;
	    offset1 = r->start;
	    offset2 = r->end;
	    if (r_out)
		*r_out = *r;
	} else {
	    idx = -1;
	}
    }

    if (idx == -1) {
	for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	    range_t *r = arrp(range_t, bin->rng, i);
	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if (r->rec == rec) {
		found = 1;
		offset1 = r->start;
		offset2 = r->end;

		if (r_out)
		    *r_out = *r;
		break;
	    }
	}
    }

    if (!found) goto fail;

    orig_start = offset1;
    orig_end   = offset2;
    orig_bin   = bin->rec;

#if 0
    {
	bin_check_cache(io, bin);
	if (contig) *contig = bin->cached_contig;
	if (start)  *start  = bin->cached_abspos + offset1;
	if (end)    *start  = bin->cached_abspos + offset2;
	if (orient) *orient = bin->cached_orient;
	/*
	printf("Orig range=%d..%d final=%d..%d in =%"PRIrec"\n",
	       orig_start, orig_end,
	       bin->cached_abspos + orig_start,
	       bin->cached_abspos + orig_end,
	       bin->cached_contig);
	*/
	return 0;
    }
#endif


    /* Find the position of this bin relative to the contig itself */
    for (;;) {
	if (bin->flags & BIN_COMPLEMENTED) {
	    offset1 = bin->size-1 - offset1;
	    offset2 = bin->size-1 - offset2;
	    comp ^= 1;
	}
	offset1 += bin->pos;
	offset2 += bin->pos;

	if (bin->parent_type != GT_Bin)
	    break;

	bnum = bin->parent;
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	if (NULL == bin) goto fail;
    }

    assert(bin->parent_type == GT_Contig);

    if (contig)
	*contig = bin->parent;
    if (start)
	*start = offset1 < offset2 ? offset1 : offset2;
    if (end)
	*end = offset1 > offset2 ? offset1 : offset2;
    if (orient)
	*orient = comp;

#if 0
    /*
    tg_rec ccc = bin->parent;
    bin = cache_search(io, GT_Bin, orig_bin);
    bin_check_cache(io, bin);
    printf("Orig range=%d..%d final=%d..%d in =%"PRIrec
	   " bin_abs=%d..%d in =%"PRIrec"\n",
	   orig_start, orig_end,
	   offset1, offset2, ccc,
	   bin->cached_abspos + orig_start,
	   bin->cached_abspos + orig_end,
	   bin->cached_contig);
    */
#endif

    return 0;
 fail:
    if (i_out) {
	cache_decr(io, *i_out);
	*i_out = NULL;
    }
    return -1;
}

/*
 * Computes the bin orientation with respect to the contig.
 * Returns 1 for complemented
 *         0 for uncomplemented.
 */
int bin_get_orient(GapIO *io, tg_rec rec) {
    bin_index_t *bin = NULL;
    int comp = 0;

    /* Bubble up bins until we hit the root */
    while (rec) {
	bin = (bin_index_t *)cache_search(io, GT_Bin, rec);
	if (bin->flags & BIN_COMPLEMENTED)
	    comp ^= 1;
	if (bin->parent_type != GT_Bin)
	    break;
	rec = bin->parent;
    }

    assert(bin && bin->parent_type == GT_Contig);
    return comp;
}

/*
 * Removes a record referenced from a known bin.
 * This function is exhaustive in its checking, but slow. So it is suitable
 * for removing small numbers of items and not for removing large numbers
 * in a loop.
 *
 * See fast_remove_item_from_bin() and fixup_bin() for a better way to
 * do bulk removal.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item_from_bin(GapIO *io, contig_t **c, bin_index_t **binp,
			     int type, tg_rec rec) {
    bin_index_t *bin;
    int i, start = INT_MAX, end = INT_MIN, bin_idx = -1;
    int new_contig_range = 0;
    int seq_start = INT_MAX, seq_end = INT_MIN;
    int item_start = INT_MAX, item_end = INT_MIN;

    if (!(bin = cache_rw(io, *binp)))
	return -1;
    *binp = bin;

    bin->flags &= ~BIN_CONS_VALID;
    bin->flags |= BIN_BIN_UPDATED;

    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (r->rec != rec) {
	    /* New start/end boundaries */
	    if (start > r->start)
		start = r->start;
	    if (end   < r->end)
		end   = r->end;

	    /* And start/end boundaries of ISSEQ items */
	    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
		if (seq_start > r->start)
		    seq_start = r->start;
		if (seq_end   < r->end)
		    seq_end   = r->end;
	    }

	    continue;
	} else {
	    item_start = r->start;
	    item_end   = r->end;
	}

	bin_idx = i;
    }

    /* Found it, and also have start/end + seq_start/seq_end */
    if (bin_idx != -1) {
	range_t *r = arrp(range_t, bin->rng, bin_idx);

	/* Fix bin extents if needed */
	if (bin->start_used != start || bin->end_used != end) {
	    if (start == INT_MAX) {
		bin->start_used = bin->end_used = 0;
	    } else {
		bin->start_used = start;
		bin->end_used   = end;
	    }

	    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
		new_contig_range = 1;
		/*
		if (seq_start > r->start)
		    seq_start = r->start;
		if (seq_end < r->end)
		    seq_end = r->end;
		*/
	    }
	}

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ &&
	    (r->start < seq_start || r->end > seq_end))
	    new_contig_range = 1;

	/* Remove from bin */
	r->flags |= GRANGE_FLAG_UNUSED;
	r->rec = (tg_rec)bin->rng_free;
	r->pair_timestamp = 0;
	bin->rng_free = bin_idx;
	bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    *c = cache_rw(io, *c);
	    bin_incr_nseq(io, bin, -1);

	    /* Also fix pair's timestamp for where we are */
	    if (r->pair_rec) {
		seq_t *s;
		bin_index_t *b;
		range_t *r2;
		
		s = cache_search(io, GT_Seq, r->pair_rec);
		b = cache_search(io, GT_Bin, s->bin);
		b = cache_rw(io, b);
		r2 = arrp(range_t, b->rng, s->bin_index);
		assert(r2->rec == s->rec);
			 
		r2->pair_timestamp = 0;
	    }

	    /* 
	     * Invalidate clipped start/end - FIXME: should check if this is
	     * really necessary.
	     */
	    (*c)->clipped_timestamp = 0;
	    
	}

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
	    *c = cache_rw(io, *c);
	    bin_incr_nrefpos(io, bin, -1);
	}

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    *c = cache_rw(io, *c);
	    bin_incr_nanno(io, bin, -1);
	}
    }

    /*
     * If we edited the bin start_used / end_used values, then possibly
     * the contig start/end may also have changed. We need to check.
     */
    if (new_contig_range) {
	int comp = 0;
	tg_rec bnum;

	for (;;) {
	    int tmp;
	    if (bin->flags & BIN_COMPLEMENTED) {
		if (seq_start != INT_MAX) {
		    tmp = bin->size-1 - seq_start;
		    seq_start   = bin->size-1 - seq_end;
		    seq_end = tmp;
		}

		tmp = bin->size-1 - item_start;
		item_start   = bin->size-1 - item_end;
		item_end = tmp;
		comp ^= 1;
	    }
	    if (seq_start != INT_MAX) {
		seq_start += bin->pos;
		seq_end   += bin->pos;
	    }

	    item_start += bin->pos;
	    item_end   += bin->pos;

	    if (bin->parent_type != GT_Bin)
		break;

	    bnum = bin->parent;
	    bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	}

	/*
	 * FIXME: this isn't entirely catching all situations right now
	 * so we always redo both clip points just incase. It's less
	 * efficient, but robust.
	 */
	if (seq_start == INT_MAX || seq_end == INT_MIN) {
	    /* Blank bin, possibly blank contig */
	    int st, en;
	    if (-1 != consensus_unclipped_range(io, (*c)->rec, &st, &en)) {
		(*c) = cache_rw(io, *c);
		(*c)->start = st;
		(*c)->end   = en;
	    }
	} else if (seq_start  <= (*c)->start || seq_end  >= (*c)->end ||
		   item_start <= (*c)->start || item_end >= (*c)->end) {
	    /*
	     * Seq is at very end of contig. So we need to do find the
	     * new start/end and correct it.
	     */
	    int new_start, *ns;
	    int new_end, *ne;

	    (*c) = cache_rw(io, *c);

	    ns = (seq_start <= (*c)->start || item_start <= (*c)->start)
		? &new_start : NULL;
	    ne = (seq_end   >= (*c)->end   || item_end   >= (*c)->end)
		? &new_end   : NULL;
	    //ns = &new_start;
	    //ne = &new_end;

	    if (-1 != consensus_unclipped_range(io, (*c)->rec, ns, ne)) {
		if (ns) (*c)->start = *ns;
		if (ne) (*c)->end   = *ne;

#if 0
		/* Technically this doesn't invalidate object
		 * locations, just contig sizes. However it may be
		 * good to reuse c->timestamp for cached_clipped_start
		 * and end?
		 */
		if (ns || ne)
		    (*c)->timestamp = io_timestamp_incr(io);
#endif
	    }
	}
    }

    return 0;
}

/*
 * Removes a record referenced from a bin. As above but we only know the
 * record and not which bin it's part of
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item(GapIO *io, contig_t **c, int type, tg_rec rec) {
    bin_index_t *bin;
    int start, end;
    tg_rec cnum, bnum;

    if (-1 == bin_get_item_position(io, type, rec, &cnum, &start, &end,
				    NULL, &bnum, NULL, NULL))
	return -1;


    //if (!(bin = bin_for_range(io, c, start, end, 0, NULL, NULL)))
    //    return -1;
    bin = cache_search(io, GT_Bin, bnum);

    return bin_remove_item_from_bin(io, c, &bin, type, rec);
}

/*
 * A faster alternative to bin_remove_item_from_bin(). Note this does not
 * ensure the contig is consistent afterwards or that the bin ranges are
 * valid. These checks need to be done afterwards by calling
 * bin_set_used_range() and possibly consensus_unclipped_range().
 *
 * If bin_idx is already known, pass it in. Otherwise pass in -1.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int fast_remove_item_from_bin(GapIO *io, contig_t **c, bin_index_t **binp,
			      int type, tg_rec rec, int bin_idx) {
    bin_index_t *bin;

    if (!(bin = cache_rw(io, *binp)))
	return -1;
    *binp = bin;

    bin->flags &= ~BIN_CONS_VALID;
    bin->flags |= BIN_BIN_UPDATED;

    if (bin_idx != -1) {
	/*  Check validity */
	if (bin->rng) {
	    range_t *r = arrp(range_t, bin->rng, bin_idx);
	    if (r->rec != rec) {
		bin_idx = -1;
	    }
	} else {
	    bin_idx = -1;
	}
    }

    if (bin_idx == -1) {
	int i;

	for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	    range_t *r = arrp(range_t, bin->rng, i);
	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if (r->rec != rec)
		continue;

	    bin_idx = i;
	}
    }

    /* Found it */
    if (bin_idx != -1) {
	range_t *r = arrp(range_t, bin->rng, bin_idx);

	/* Remove from bin */
	r->flags |= GRANGE_FLAG_UNUSED;
	r->rec = (tg_rec)bin->rng_free;
	bin->rng_free = bin_idx;
	bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    *c = cache_rw(io, *c);
	    bin_incr_nseq(io, bin, -1);
	}

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
	    *c = cache_rw(io, *c);
	    bin_incr_nrefpos(io, bin, -1);
	}

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    *c = cache_rw(io, *c);
	    bin_incr_nanno(io, bin, -1);
	}
    }

    return 0;
}


/*
 * Call after updating range array to ensure that the bin start_used and
 * end_used values are correct.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_set_used_range(GapIO *io, bin_index_t *bin) {
    int i;
    int start = INT_MAX, end = INT_MIN;

    for (i = 0; i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (start > r->start)
	    start = r->start;
	if (end   < r->end)
	    end   = r->end;
    }
	
    if (start != INT_MAX) {
	if (bin->start_used == start && bin->end_used == end)
	    return 0;

	if (NULL == (bin = cache_rw(io, bin)))
	    return -1;

	bin->start_used = start;
	bin->end_used = end;
    } else {
	if (bin->start_used == 0 && bin->end_used == 0)
	    return 0;

	if (NULL == (bin = cache_rw(io, bin)))
	    return -1;

	bin->start_used = 0;
	bin->end_used = 0;
    }

    return 0;
}

/*
 * Removes the refpos marker at a specific position in the contig.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_refpos(GapIO *io, tg_rec crec, int pos) {
    tg_rec bin_rec;
    int bin_idx;
    rangec_t rc;
    range_t *r;
    bin_index_t *bin;

    if (find_refpos_marker(io, crec, pos, &bin_rec, &bin_idx, &rc) != 0)
	return -1;


    assert((rc.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS);

    bin = cache_search(io, GT_Bin, bin_rec);
    bin = cache_rw(io, bin);
    r = arrp(range_t, bin->rng, bin_idx);

    r->flags |= GRANGE_FLAG_UNUSED;
    r->rec = bin->rng_free;
    bin->rng_free = bin_idx;
    bin_incr_nrefpos(io, bin, -1);
    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;

    if (bin->start_used == r->start || bin->end_used == r->end)
	return bin_set_used_range(io, bin);

    return 0;
}


/*
 * Finds the contig number and position of the start of a bin.
 * The right position is obviously this + bin->size (regardless of
 * whether it has been complemented).
 *
 * Returns 0 on success (plus *contig & *pos)
 *        -1 on failure
 */
int bin_get_position(GapIO *io, bin_index_t *bin, tg_rec *contig, int *pos,
		     int *comp_p) {
    tg_rec bnum;
    int offset1 = 0, offset2 = bin->size-1;
    int comp = 0;

    /* Find the position of this bin relative to the contig itself */
    for (;;) {
	if (bin->flags & BIN_COMPLEMENTED) {
	    offset1 = bin->size-1 - offset1;
	    offset2 = bin->size-1 - offset2;
	    comp ^= 1;
	}
	offset1 += bin->pos;
	offset2 += bin->pos;

	if (bin->parent_type != GT_Bin)
	    break;

	bnum = bin->parent;
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
    }
    
    assert(bin->parent_type == GT_Contig);
    *contig = bin->parent;

    *pos = MIN(offset1, offset2);
    if (comp_p)
	*comp_p = comp;

    return 0;
}

/*
 * Some tracks are "official" cached objects and are deallocated as part of
 * the tg_cache system. Others are temporary on-the-fly structs generated
 * for the purpose of temporary display. These we need to free.
 *
 * This function knows which is which and does the appropriate thing.
 */
void track_free(track_t *t) {
    if (t->flag & TRACK_FLAG_FREEME) {
	if (t->data)
	    ArrayDestroy(t->data);
	free(t);
    }
}

/*
 * Creates a fake track struct, to be freed with track_free
 *
 * Returns track_t on success
 *         NULL on failure
 */
track_t *track_create_fake(int type, int size) {
    track_t *t = (track_t *)calloc(1, sizeof(*t));
    if (!t)
	return 0;
    t->type = type;
    t->nitems = size;
    t->item_size = sizeof(int);
    t->data = ArrayCreate(sizeof(int), size);
    t->flag |= TRACK_FLAG_FREEME;

    return t;
}

/*
 * Creates a track of a given type for this bin.
 * Note, this does not actually add it to the bin (but probably should
 * otherwise it's nothing more than the non-bin track_create).
 *
 * Returns track_t pointer on success
 *         NULL on failure
 */
track_t *bin_create_track(GapIO *io, bin_index_t *bin, int type) {
    tg_rec rec;
    track_t *t;

    if (-1 == (rec = io->iface->track.create(io->dbh, NULL)))
	return NULL;

    t = (track_t *)cache_search(io, GT_Track, rec);
    track_set_type(io, &t, type);

    return t;
}

/*
 * Adds a track to this bin.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_add_track(GapIO *io, bin_index_t **bin, track_t *track) {
    bin_index_t *n;
    int i;
    bin_track_t *bt;

    if (!(n = cache_rw(io, *bin)))
	return -1;
    *bin = n;

    /* Create new bin-track, or error if already found */
    if (!n->track) {
	n->track = ArrayCreate(sizeof(bin_track_t), 0);
	n->flags |= BIN_TRACK_UPDATED;
    }

    for (i = 0; i < ArrayMax(n->track); i++) {
	bt = arrp(bin_track_t, n->track, i);
	if (bt->type == track->type)
	    return -1;
    }

    /* Add the track pointer */
    bt = (bin_track_t *)ArrayRef(n->track, ArrayMax(n->track));
    bt->type  = track->type;
    bt->flags = 1;
    bt->rec   = track->rec;
    bt->track = track;

    return 0;
}

/*
 * Finds the track of a given type for this bin.
 *
 * Returns track_t pointer on success (do not free)
 *         NULL on failure (eg bin too small)
 */
track_t *bin_get_track(GapIO *io, bin_index_t *bin, int type) {
    int i;

    /* If it exists and is up to date, return it */
    if (bin->track) {
	for (i = 0; i < ArrayMax(bin->track); i++) {
	    bin_track_t *bt = arrp(bin_track_t, bin->track, i);
	    if (bt->type == type) {
		if (bt->track)
		    return bt->track;
		
		return (track_t *)cache_search(io, GT_Track, bt->rec);
	    }
	}
    }

    return NULL;
}

/*
 * A bit like bin_get_track, but this is designed to auto-generate and
 * update the track as desired. The expectation is that this will always
 * succeed and anything else is a fatal error.
 */
track_t *bin_query_track(GapIO *io, bin_index_t *bin, int type) {
    track_t *bin_recalculate_track(GapIO *io, bin_index_t *bin, int type);
    int i;

    /* If it exists and is up to date, return it */
    if (bin->track) {
	for (i = 0; i < ArrayMax(bin->track); i++) {
	    bin_track_t *bt = arrp(bin_track_t, bin->track, i);
	    if (bt->type == type && (bt->flags & TRACK_FLAG_VALID))
		return (track_t *)cache_search(io, GT_Track, bt->rec);
	}
    }

    /* Otherwise generate and maybe cache */
    return bin_recalculate_track(io, bin, type);
}

/*
 * Invalidates a track.
 * Returns 0 on success
 *        -1 on failure.
 */
int bin_invalidate_track(GapIO *io, bin_index_t *bin, int type) {
    int i;

    if (!bin->track)
	return 0;

    for (i = 0; i < ArrayMax(bin->track); i++) {
	bin_track_t *bt = arrp(bin_track_t, bin->track, i);
	if (bt->type == type || type == TRACK_ALL) {
	    if (NULL == (bin = cache_rw(io, bin)))
		return -1;
	    printf("Update track for rec %"PRIrec"\n", bin->rec);
	    bin->flags |= BIN_TRACK_UPDATED;
	    bt = arrp(bin_track_t, bin->track, i);
	    bt->flags &= ~TRACK_FLAG_VALID;
	}
    }

    return 0;
}

/*
 * ---------------------------------------------------------------------------
 * Track handling
 */

track_t *bin_recalculate_track(GapIO *io, bin_index_t *bin, int type) {
    int pos;
    tg_rec cnum;
    track_t *track, *child;
    int nele;
    double bpv;
    contig_t *c;

    /*
     * So we have a bin of a given size in which we wish to have at least
     * RD_ELEMENTS of track samples, but it's out of date. We query the
     * contig for track data at a resolution of double what we need and
     * then average/downsample it to generate our new bin stats.
     */
    bpv = (double)bin->size / RD_ELEMENTS;
    if (bpv < 1) bpv = 1;
    nele = bin->size / bpv;
    if (nele & 1) nele++;
    bpv = (double)bin->size / nele;

    /*
     * Bottom layer, so no point querying a child - we just create it
     * ourselves now.
     */
    if (bpv <= 2) { /* FIXME: was 1, need 1 to stop aliasing? */
	track_t *fake;
	tg_rec rec;
	int *depth;
	fake = track_create_fake(type, bin->size);
	fake->flag = TRACK_FLAG_FREEME;

	/* FIXME: type specific code here - or the size at least (int) */
	switch (type) {
	case TRACK_READ_DEPTH:
	    depth = track_read_depth_r1(io, bin);
	    break;
	default:
	    fprintf(stderr, "Unknown track type %d\n", type);
	    return NULL;
	}
	memcpy(ArrayBase(int, fake->data), depth, bin->size * sizeof(int));
	free(depth);

	rec = io->iface->track.create(io->dbh, fake);
	track = (track_t *)cache_search(io, GT_Track, rec);

	printf("Initialising track %"PRIrec" flag %d in bin %"PRIrec
	       " at bpv of 1\n",
	       rec, track->flag, bin->rec);

	bin_add_track(io, &bin, track);
	track_free(fake);

	return track;
    }

    /* Else use child bin tracks, in ever-decreasing circles */
    if (-1 == bin_get_position(io, bin, &cnum, &pos, NULL))
	return NULL;
    c = (contig_t *)cache_search(io, GT_Contig, cnum);
    child = contig_get_track(io, &c, pos, pos + bin->size-1, type, bpv);
    if (NULL == child)
	return NULL;

    track = bin_get_track(io, bin, type);
    if (!track) {
	track = bin_create_track(io, bin, type);
	bin_add_track(io, &bin, track);
    }

    /* Copy child 'fake track' to our real track */
    track_set_data(io, &track, ArrayCreate(sizeof(int), nele));
    track_set_nitems(io, &track, nele);
    track_set_item_size(io, &track, sizeof(int));
    memcpy(ArrayBase(int, track->data), ArrayBase(int, child->data),
	   nele * sizeof(int));

    track_free(child);
    //cache_decr(io, track);

    return track;
}


/*
 * Clears the BIN_CONS_VALID flag on any bins that may contain
 * consensus data over the region contig:start..end.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_invalidate_consensus(GapIO *io, tg_rec contig, int start, int end) {
    int i, nr;
    rangec_t *r;
    contig_t *c;
    
    if (NULL == (c = (contig_t *)cache_search(io, GT_Contig, contig)))
	return -1;

    /* Also update timestamp as invalidating consensus may also imply
     * changing consensus length. We may be able to improve on this when
     * we know the invalidation is purely substitutions, but this is a
     * good fallback position.
     */
    c = cache_rw(io, c);
    c->timestamp = io_timestamp_incr(io);

    r = contig_bins_in_range(io, &c, start, end, 0, CONS_BIN_SIZE/2, &nr);

    for (i = 0; i < nr; i++) {
	bin_index_t *bin = (bin_index_t *)cache_search(io, GT_Bin, r[i].rec);
	if (!bin)
	    return -1;

	if (bin->flags & BIN_CONS_VALID) {
	    bin = cache_rw(io, bin);
	    bin->flags |=  BIN_BIN_UPDATED;
	    bin->flags &= ~BIN_CONS_VALID;
	}
    }

    if (r)
	free(r);

    return 0;
}

#if 0
/*
 * Recursively fixes cached bin values.
 */
static int bin_fix_cache(GapIO *io, bin_index_t *bin) {
    if (gio_base(io)->db->timestamp == bin->timestamp)
	return 0;

    cache_incr(io, bin);

    printf("Fixing bin #%"PRIrec" time=%d vs %d\n", bin->rec,
	   bin->timestamp, gio_base(io)->db->timestamp);

    if (bin->parent_type == GT_Bin) {
	bin_index_t *p = cache_search(io, GT_Bin, bin->parent);
	if (!p) {
	    cache_decr(io, bin);
	    abort();
	    return -1;
	}
	cache_incr(io, p);

	if (bin_fix_cache(io, p) != 0) {
	    cache_decr(io, bin);
	    cache_decr(io, p);
	    abort();
	    return -1;
	}

	/* compute vs parent */
	bin->cached_contig = p->cached_contig;
	if (bin->flags & BIN_COMPLEMENTED) {
	    bin->cached_orient = p->cached_orient ^ 1;
	    bin->cached_abspos = p->cached_abspos + bin->size-1 - bin->pos;
	} else {
	    bin->cached_orient = p->cached_orient;
	    bin->cached_abspos = p->cached_abspos + bin->pos;
	}

	cache_decr(io, p);
    } else {
	/* Root bin */
	bin->cached_contig = bin->parent;
	bin->cached_orient = (bin->flags & BIN_COMPLEMENTED) ? 1 : 0;
	bin->cached_abspos = bin->pos;
    }

    bin->timestamp = gio_base(io)->db->timestamp;
    cache_decr(io, bin);

    return 0;
}

/*
 * Ensures the cached_abspos and cached_contig fields are up to date.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_check_cache(GapIO *io, bin_index_t *bin) {
    int ctg_correct = 0;

    if (!bin)
	return -1;

    /*
     * Check bin vs database first as this needs no IO.
     * If that fails, check if contig is valid and if that timestamp is
     * correc too.
     */
    if (bin->cached_contig) {
	if (bin->timestamp == gio_base(io)->db->timestamp) {
	    ctg_correct = 1;
	} else if (cache_exists(io, GT_Contig, bin->cached_contig)){
	    contig_t *c = cache_search(io, GT_Contig, bin->cached_contig);
	    if (bin->timestamp >= c->timestamp)
		ctg_correct = 1;
	}
    }

    if (ctg_correct)
	return 0;

    /* The bin is out of date. So correct it */
    return bin_fix_cache(io, bin);
}
#endif
