/*
 * File: readpair.c:
 * Version: 2.0 (1.0 == FORTRAN version)
 *
 * Author: Staden Package group
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: "Find read pair" code.
 *
 * Created: 17 March 1994
 * Updated:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#include "xalloc.h"
#include "tg_gio.h"
#include "gap4_compat.h"
#include "cs-object.h"
#include "contig_selector.h"
#include "newgap_cmds.h"
#include "gap_globals.h"
#include "misc.h"
#include "text_output.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "readpair.h"
#include "editor_view.h"
#include "tk-io-reg.h"
#include "io_lib/hash_table.h"

/*
 * Match callback.
 * 'obj' is a match contained within the 'repeat' list.
 */
void *readpair_obj_func(int job, void *jdata, obj_read_pair *obj,
			mobj_read_pair *template) {
    static char buf[200];
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(template->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(template->io, cs_id);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	if (io_rdonly(template->io) && ((obj->c1 > 0 && obj->c2 < 0) ||
					(obj->c1 < 0 && obj->c2 > 0))) {
	    return "Information\0Hide\0IGNORE\0"
		"IGNORE\0SEPARATOR\0Remove\0";
	} else {
	    return "Information\0Hide\0Invoke join editor *\0"
		        "SEPARATOR\0Remove\0";
	}
	break;

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */
	    start_message();

	    vmessage("Read pair:\n");
	    vmessage("    From contig %s(=%"PRIrec") at %d reading %s(#%"PRIrec")\n",
		     get_contig_name(template->io, ABS(obj->c1)),
		     ABS(obj->c1), obj->pos1,
		     get_read_name(template->io, obj->read1), obj->read1);
	    vmessage("    With contig %s(=%"PRIrec") at %d reading %s(#%"PRIrec")\n",
		     get_contig_name(template->io, ABS(obj->c2)),
		     ABS(obj->c2), obj->pos2,
		     get_read_name(template->io, obj->read2), obj->read2);
	    {
		seq_t *s;;

		s = cache_search(template->io, GT_Seq, obj->read1);
		vmessage("    Direction of first read is %swards\n",
			 (s->flags & SEQ_END_MASK) == SEQ_END_FWD
			 ? "for" : "back");

		s = cache_search(template->io, GT_Seq, obj->read2);
		vmessage("    Direction of second read is %swards\n",
			 (s->flags & SEQ_END_MASK) == SEQ_END_FWD
			 ? "for" : "back");
	    }
	    vmessage("    Length %d\n\n", obj->length);
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, (obj_match *)obj,
		     (mobj_repeat *)template, csplot_hash);
	    break;

	case -2: /* default */
	case 2: /* Join editor */
	    {
	        tg_rec cnum[2], llino[2];
		int pos[2];
		seq_t *s;
		int comp, shortest;

		cnum[0] = ABS(obj->c1);
		cnum[1] = ABS(obj->c2);

		/* Complement a contig if needed */
		if ((obj->c1 > 0) != (obj->c2 > 0)) {
		    if (cnum[0] == cnum[1]) {
			verror(ERR_WARN, "join_editor",
			       "cannot display the same contig in two "
			       "different orientations");
			break;
		    }
		    
		    shortest = (io_clength(template->io, cnum[0])
				< io_clength(template->io, cnum[1])) ? 0 : 1;
		    if (-1 == complement_contig(template->io, cnum[shortest]))
			if (-1 == complement_contig(template->io,
						    cnum[1 - shortest]))
			    return NULL;
		}

		llino[0] = obj->read1;
		llino[1] = obj->read2;

		comp = sequence_get_orient(template->io, obj->read1);
		s = cache_search(template->io, GT_Seq, obj->read1);
		if (NULL == s) return NULL;
		pos[0] = comp ? ABS(s->len) - s->right : s->right - 1;

		comp = sequence_get_orient(template->io, obj->read2);
		s = cache_search(template->io, GT_Seq, obj->read2);
		if (NULL == s) return NULL;
		pos[1] = comp ? ABS(s->len) - s->right : s->right - 1;
		join_contig(template->io, cnum, llino, pos);

		break;
	    }

	case 3: /* Remove */
	    obj_remove(GetInterp(), cs->window, (obj_match *)obj,
		       (mobj_repeat *)template, csplot_hash);
	    break;

	}
	break;

    case OBJ_GET_BRIEF:
	snprintf(buf, sizeof(buf),
		 "Read pair: %c#%"PRIrec"@%d (mq %d) with %c#%"PRIrec"@%d (mq %d), len %d",
		 obj->c1 > 0 ? '+' : '-', obj->read1, obj->pos1, obj->mq1,
		 obj->c2 > 0 ? '+' : '-', obj->read2, obj->pos2, obj->mq2,
		 obj->length);
	return buf;
    }

    return NULL;
}

void readpair_callback(GapIO *io, tg_rec contig, void *fdata, reg_data *jdata) {
    mobj_read_pair *r = (mobj_read_pair *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Find read pairs");
	break;


    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, (mobj_repeat *)r,
			csplot_hash, cs->window);
	break;


    case REG_COMPLEMENT:

	csmatch_complement(io, contig, r, csplot_hash, cs->window);
	break;


    case REG_GET_OPS:

	if (r->all_hidden)
	    jdata->get_ops.ops = "PLACEHOLDER\0PLACEHOLDER\0Information\0"
		"PLACEHOLDER\0Hide all\0Reveal all\0"
		"Save matches\0SEPARATOR\0Remove\0";
	else
	    jdata->get_ops.ops = "Use for 'Next'\0Reset 'Next'\0Information\0"
		"Configure\0Hide all\0Reveal all\0"
		"Save matches\0SEPARATOR\0Remove\0";
	break;


    case REG_INVOKE_OP:

	switch (jdata->invoke_op.op) {
	case 0: /* Next */
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(r), NULL);
	    break;
	case 1: /* Reset Next */
	    csmatch_reset_next((mobj_repeat *)r);
	    break;
	case 2: /* Information */
	    csmatch_info((mobj_repeat *)r, "Read pair");
	    break;
	case 3: /* Configure */
	    csmatch_configure(io, cs->window, (mobj_repeat *)r);
	    break;
	case 4: /* Hide all */
	    csmatch_hide(GetInterp(), cs->window, (mobj_repeat *)r,
			 csplot_hash);
	    break;
	case 5: /* Reveal all */
	    csmatch_reveal(GetInterp(), cs->window, (mobj_repeat *)r,
			   csplot_hash);
	    break;
	case 6: { /* Save */
	    char *fn;
	    if (Tcl_VarEval(GetInterp(), "tk_getSaveFile ", "-parent ",
			    cs->window, NULL) != TCL_OK)
		break;
	    fn = Tcl_GetStringResult(GetInterp());
	    if (fn && *fn)
		csmatch_save((mobj_generic *)r, fn);
	    break;
	}
	case 7: /* Remove */
	    csmatch_remove(io, cs->window,
			   (mobj_repeat *)r,
			   csplot_hash);
	    break;
	}
	break;


    case REG_PARAMS:

	jdata->params.string = r->params;
	break;


    case REG_NUMBER_CHANGE:

	csmatch_renumber(io, contig, jdata->number.number,
			 (mobj_repeat *)r, csplot_hash, cs->window);
	break;


    case REG_ORDER:

	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs->window,
		       (mobj_repeat *)r,
		       csplot_hash);
	break;


    case REG_DELETE:

	csmatch_contig_delete(io, (mobj_repeat *)r, contig,
			      cs->window, csplot_hash);
	break;

    case REG_LENGTH:
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_GENERIC:
	switch (jdata->generic.task) {
	    int ret;

	case TASK_CS_PLOT:
	    PlotRepeats(io, (mobj_repeat *)r);
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(r), NULL);
	    break;

	case TASK_CS_SAVE:
	    ret = csmatch_save((mobj_generic *) r, (char *)jdata->generic.data);
	    vTcl_SetResult(GetInterp(), "%d", ret);
	    break;
	}
    }
}

/*
 * Plot the templates which span more than one contig on contig selector plot.
 * The direction of the line is governed by the "direction" of the template
 * A line going from top right to bottom left (/) indicates the second contig
 * needs completing
 * A line going from top left to bottom right (\) indicates the second contig
 * does not need completing
 */
int PlotTempMatches(GapIO *io, read_pair_t *rp) {
    int id;
    mobj_read_pair *template;
    obj_read_pair *matches;
    int n_matches = 0;
    int max_matches = 64; /* dynamically grows */
    char *val;

    if (!rp)
	return 0;

    if (NULL == (template = (mobj_read_pair *)xmalloc(sizeof(mobj_read_pair))))
	return -1;
    if (NULL == (matches = (obj_read_pair *)xmalloc(max_matches *
						    sizeof(obj_read_pair))))
	return -1;

    /* Create cs-object plot array */
    while (rp->rec[0]) {
	matches[n_matches].func =
	    (void *(*)(int, void *, struct obj_match_t *,
		       struct mobj_repeat_t *))readpair_obj_func;
	matches[n_matches].data = template;
	/* Fix dir on c1/c2 */
	matches[n_matches].c1 = rp->contig[0];
	matches[n_matches].c2 = rp->contig[1];
	matches[n_matches].pos1 = rp->start[0];
	matches[n_matches].pos2 = rp->start[1];
	matches[n_matches].end1 = rp->end[0];
	matches[n_matches].end2 = rp->end[1];
	matches[n_matches].length = (ABS(rp->end[0] - rp->start[0]) + 
				     ABS(rp->end[1] - rp->start[1])) / 2;
	matches[n_matches].read1 = rp->rec[0];
	matches[n_matches].read2 = rp->rec[1];
	matches[n_matches].flags = 0;
	matches[n_matches].mq1 = rp->mqual[0];
	matches[n_matches].mq2 = rp->mqual[1];
	if (++n_matches >= max_matches) {
	    obj_read_pair *tmp;
	    
	    max_matches *= 2;
	    tmp = (obj_read_pair *)
		xrealloc(matches,
			 max_matches * sizeof(obj_read_pair));
	    if (NULL == tmp) {
		xfree(template);
		xfree(matches);
		return -1;
	    } else {
		matches = tmp;
	    }
	}

	rp++;
    }


    /*
     * Free memory and return if we've got no matches.
     */
    if (0 == n_matches) {
	xfree(template);
	xfree(matches);

	return 0;
    }

    template->num_match = n_matches;
    template->match = (obj_match *)matches;
    template->io = io;
    strcpy(template->tagname, CPtr2Tcl(template));

    val = get_default_string(GetInterp(), gap5_defs, "READPAIR.COLOUR");
    strcpy(template->colour, val);

    template->linewidth = get_default_int(GetInterp(), gap5_defs,
					  "READPAIR.LINEWIDTH");

    template->params = (char *)xmalloc(10);
    if (template->params)
	sprintf(template->params, "none");
    template->all_hidden = 0;
    template->current = -1;
    template->reg_func = readpair_callback;
    template->match_type = REG_TYPE_READPAIR;

    /*
     * Register the repeat search with each of the contigs used.
     * Currently we assume that this is all.
     */
    id = register_id();
    contig_register(io, 0, readpair_callback, (void *)template, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
		    REG_NUMBER_CHANGE | REG_ORDER | REG_GENERIC,
		    REG_TYPE_READPAIR);
    update_results(io);

    return id;
}

HashTable *create_lib_hash(tg_rec *library, int nlibrary) {
    int i;
    HashTable *lib_hash = HashTableCreate(16, HASH_FUNC_HSIEH | 
					  HASH_DYNAMIC_SIZE |
					  HASH_POOL_ITEMS);
    if (NULL == lib_hash) return NULL;

    for (i = 0; i < nlibrary; i++) {
	    HashData hd;
	    hd.i = 1;
	    if (!HashTableAdd(lib_hash, (char *)&library[i], sizeof(tg_rec),
			      hd, 0)) {
		HashTableDestroy(lib_hash, 0);
		return NULL;
	    }
    }

    return lib_hash;
}

/*
 * Identifies spanning read pairs within end_size of the end of each contig.
 * Library/nlibrary is used to filter read-pairs as coming only from specific
 * libraries. It can be specified as NULL/0 if no filtering is required,
 * otherwise we will only accept pairs that are contained within the lib
 * list.
 *
 * Mode is 0 for all vs all.
 *         1 for ends vs all.
 *	   2 for ends vs ends.
 */

read_pair_t *spanning_pairs(GapIO *io, int num_contigs,
			    contig_list_t *contig_array,
			    enum readpair_mode mode,
			    int end_size, int min_mq, int min_freq,
			    HashTable *lib_hash) {
    int i, j;
    HashTable *h = NULL;
    HashIter *iter = NULL;
    HashItem *hi;
    read_pair_t *pairs = NULL;
    int npairs = 0, nalloc = 0;
    int no_large_contigs = 1;
    int slow_check_libs = 0;
    pool_alloc_t *rp_pool = NULL;
    HashTable *ctg_hash = NULL;
    HashTable *ctgs_in_list = NULL;
    tg_rec ctg_pair[2];

    h = HashTableCreate(1024, HASH_FUNC_HSIEH |
			      HASH_DYNAMIC_SIZE |
			      HASH_POOL_ITEMS);
    if (NULL == h) return NULL;
    rp_pool = pool_create(sizeof(rangec_t));
    if (NULL == rp_pool) goto fail;

    /*
     * For filtering by frequency. We construct a key consisting of
     * two contig recs and then increment that value during pair generation.
     * Any contig-pair with insufficient depth of read-pair is then
     * rejected;
     */
    if (min_freq > 0) {
	ctg_hash = HashTableCreate(1024, HASH_FUNC_HSIEH | 
				         HASH_DYNAMIC_SIZE |
				         HASH_POOL_ITEMS);
	if (NULL == ctg_hash) goto fail;
    }

    if (end_all == mode) {
	/*
	 * Make a hash table of the contigs in contig_array.
	 * This is needed to filter out pairs that span from
	 * contigs in the list to ones that are not.
	 */
	ctgs_in_list = HashTableCreate(num_contigs,
				       HASH_FUNC_HSIEH | 
				       HASH_DYNAMIC_SIZE |
				       HASH_POOL_ITEMS);
	if (NULL == ctgs_in_list) goto fail;
	for (i = 0; i < num_contigs; i++) {
	    HashData hd;
	    hd.i = 0;
	    if (NULL == HashTableAdd(ctgs_in_list,
				     (char *) &contig_array[i].contig,
				     sizeof(contig_array[i].contig),
				     hd, NULL)) {
		goto fail;
	    }
	}
    }

    for (i = 0; i < num_contigs; i++) {
	contig_iterator *ci;
	rangec_t *r;
	contig_t *c;
	int cstart, cend;
	tg_rec crec = contig_array[i].contig;
	int large_contig;
	
	c = cache_search(io, GT_Contig, crec);
	if (NULL == c) goto fail;
	cstart = c->start;
	cend   = c->end;

	/*
	 * For end vs end, can we optimise this into one linear scan through
	 * the contig instead of having an iterator at each end?
	 * Small contigs allow this (we just discard the data in the middle).
	 *
	 * For end vs all, we have to trade off scanning through data that
	 * isn't at the end vs finding the "all" part of the match?
	 * sequence_get_position() is very slow so we prefer to scan through
	 * bins far more here than brute force identification of the
	 * other end.
	 *
	 * For all vs all, we have to do all the work anyway so just
	 * face up to it.
	 */
	switch (mode) {
	default: /* To keep gcc happy */
	case all_all:
	    large_contig = 0;
	    break;

	case end_end:
	    large_contig = (cend-cstart)/4  >= end_size ? 1 : 0;
	    break;

	case end_all:
	    large_contig = (cend-cstart)/100 >= end_size ? 1 : 0;
	    break;
	}
	if (large_contig)
	    no_large_contigs = 0;

	/*
	 * Iterate through seqs in contig.
	 *
	 * Pair up reads, removing internal pairs as we go.
	 * We stop accumulating into our pair hash after end_size bases,
	 * but keep going to remove pairs until 2*end_size.
	 */
	ci = contig_iter_new(io, crec, 1, CITER_FIRST,
			     CITER_CSTART, CITER_CEND);
	if (NULL == ci) goto fail;

	while (NULL != (r = contig_iter_next(io, ci))) {
	    if (!r->pair_rec)
		continue; /* unpaired */

	    if (lib_hash) {
		if (r->library_rec) {
		    if (!HashTableSearch(lib_hash, (char *)&r->library_rec,
					 sizeof(r->library_rec))) {
			continue;
		    }
		} else {
		    slow_check_libs = 1;
		}
	    }

	    if (large_contig && r->start - cstart > 2*end_size)
		break;

	    if ((hi = HashTableSearch(h, (char *)&r->pair_rec,
				      sizeof(r->pair_rec)))) {
		rangec_t *r2 = hi->data.p;
		if (r2->orig_rec == crec) {
		    /* internal read pair */
		    pool_free(rp_pool, r2);
		    HashTableDel(h, hi, 0);
		    continue;
		}
	    }

	    if (mode == all_all || mode == end_all ||
		r->start - cstart < end_size || cend - r->start < end_size) {
		HashData hd;
		rangec_t *r2 = pool_alloc(rp_pool);
		if (NULL == r2) goto fail;
		*r2 = *r;
		r2->orig_rec = crec; /* convenient place to store contig */
		if (r->start - cstart < end_size ||
		    cend - r->start < end_size)
		    r2->pair_start = 1; /* near end */
		else
		    r2->pair_start = 0;
		hd.p = r2;

		if (!HashTableAdd(h, (char *)&r->rec, sizeof(r->rec),
				  hd, 0)) goto fail;
	    }
	}

	contig_iter_del(ci);

	/* Other end, if not yet done (for large contigs) */
	if (large_contig) {
	    ci = contig_iter_new(io, crec, 1, CITER_FIRST,
				 MAX(cstart, cend - 2*end_size), CITER_CEND);
	    if (NULL == ci) goto fail;

	    while (NULL != (r = contig_iter_next(io, ci))) {
		if (!r->pair_rec)
		    continue; /* unpaired */

		if (lib_hash) {
		    if (r->library_rec) {
			if (!HashTableSearch(lib_hash, (char *)&r->library_rec,
					     sizeof(r->library_rec))) {
			    continue;
			}
		    } else {
			slow_check_libs = 1;
		    }
		}

		if ((hi = HashTableSearch(h, (char *)&r->pair_rec,
					  sizeof(r->pair_rec)))){
		    rangec_t *r2 = hi->data.p;
		    if (r2->orig_rec == crec) {
			/* internal read pair */
			pool_free(rp_pool, r2);
			HashTableDel(h, hi, 0);
		    }
		}

		if (mode == end_all || cend - r->start < end_size) {
		    HashData hd;
		    rangec_t *r2 = pool_alloc(rp_pool);
		    if (NULL == r2) goto fail;
		    *r2 = *r;
		    r2->orig_rec = crec;
		    if (r->start - cstart < end_size ||
			cend - r->start < end_size)
			r2->pair_start = 1; /* near end */
		    else
			r2->pair_start = 0;
		    hd.p = r2;

		    if (!HashTableAdd(h, (char *)&r->rec, sizeof(r->rec),
				      hd, 0)) goto fail;
		}
	    }

	    contig_iter_del(ci);
	}
    }

    /*
     * Our hash table ideally now contains pairs with rec1 -> rec2 and 
     * rec2 -> rec1. This will be the case for mode == end_end or all_all.
     * Any single rec1 -> rec2 with no rec2 present in the hash will be due
     * to a read-pair outside of the region (end_end) or outside of the
     * set of contigs we chose to scan (all_all).
     *
     * For mode end_all we may have had some large contigs where we chose
     * not to scan through all reads, only effectively adding rec1 for the
     * ends of those contigs. In that case it is expected rec2 may be absent
     * and instead we need to run sequence_get_position instead. This is
     * slower, which is why we only resort to this route on extra-large
     * contigs.
     */
    iter = HashTableIterCreate();
    if (NULL == iter) goto fail;
    while ((hi = HashTableIterNext(h, iter))) {
	tg_rec rec1, rec2;
	rangec_t *r1 = (rangec_t *)hi->data.p, *r2, r2_tmp;
	int dir[2];

	rec1 = r1->rec;
	rec2 = r1->pair_rec;

	if (rec1 == 0)
	    continue; // already added other end */

	if (!(hi = HashTableSearch(h, (char *)&rec2, sizeof(rec1)))) {
	    tg_rec contig;
	    int start, end, orient;
	    range_t pair_range;

	    if (mode != end_all || no_large_contigs)
		continue; // unpaired

	    /* Maybe rec2 was in a large contig that we didn't scan through */
	    /* Should we use pair_* from r1 here ? */
	    sequence_get_position2(io, rec2, &contig, &start, &end, &orient,
				   NULL, &pair_range, NULL);
	    if (!HashTableSearch(ctgs_in_list,
				 (char *)&contig, sizeof(contig))) {
		/* Filter out unwanted contigs */
		continue;
	    }
	    r2 = &r2_tmp;
	    r2->rec = rec2;
	    r2->orig_rec = contig;
	    r2->start = start;
	    r2->end = end;
	    r2->comp = orient;
	    r2->pair_start = 0; 
	    r2->mqual = pair_range.mqual;
	    r2->flags = pair_range.flags;
	} else {
	    r2 = (rangec_t *)hi->data.p;
	}


	if (r1->orig_rec == r2->orig_rec)
	    /* Same contig */
	    continue;

	if (mode == end_all && !(r1->pair_start || r2->pair_start))
	    /* Spanning pair, but not with a match near an end */
	    continue;

	/* Mapped to a repeat? */
	if (r1->mqual < min_mq || r2->mqual < min_mq)
	    continue;

	/*
	 * Filter by library for older database versions. We do this late as
	 * it's a slow process due to needing to load the sequence struct
	 * itself rather than just iterate through the range arrays.
	 * If library_rec is set for some rangec_t's then they should
	 * have been filtered earlier.
	 */
	if (slow_check_libs && !r1->library_rec) {
	    seq_t *s = cache_search(io, GT_Seq, r1->rec);
	    tg_rec lib;

	    /* Assume r1 and r2 are both in the same library. Should be! */
	    if (s->parent_type != GT_Library)
		continue;

	    lib = s->parent_rec;
	    if (!HashTableSearch(lib_hash, (char *)&lib, sizeof(lib))) {
		continue;
	    }
	}

	/* Accepted - add it to pairs[] array */
	if (npairs >= nalloc) {
	    read_pair_t *new_pairs;
	    nalloc = nalloc ? nalloc*2 : 1024;
	    new_pairs = realloc(pairs, nalloc * sizeof(*pairs));
	    if (!new_pairs) goto fail;
	    pairs = new_pairs;
	}

	/*
	  r1 flags = 0x61 (1+20+40) => paired +!comp1 + comp2 + 1fwd + 2rev
	  r2 flags = 0x15 (1+4+10)  => paired + comp1 +!comp2 + 1rev + 2fwd
	 */

	dir[0] = (!!(r1->flags & GRANGE_FLAG_COMP1) ==
		  ((r1->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD))
	    ^ r1->comp ? -1 : 1;
	dir[1] = (!!(r2->flags & GRANGE_FLAG_COMP1) ==
		  ((r2->flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD))
	    ^ r2->comp ? -1 : 1;

	pairs[npairs].rec[0] = rec1;
	pairs[npairs].rec[1] = rec2;
	pairs[npairs].contig[0] = dir[0] * r1->orig_rec;
	pairs[npairs].contig[1] = dir[1] * r2->orig_rec;
	pairs[npairs].start[0] = r1->start;
	pairs[npairs].start[1] = r2->start;
	pairs[npairs].end[0] = r1->end;
	pairs[npairs].end[1] = r2->end;
	pairs[npairs].mqual[0] = r1->mqual;
	pairs[npairs].mqual[1] = r2->mqual;

	if (min_freq > 0) {
	    HashData hd;

	    ctg_pair[0] = pairs[npairs].contig[0];
	    ctg_pair[1] = pairs[npairs].contig[1];
	    hd.i = 0;
	    hi = HashTableAdd(ctg_hash, (char *)ctg_pair, sizeof(ctg_pair),
			      hd, NULL);
	    if (NULL == hi) goto fail;
	    hi->data.i++;
	}

	r2->rec = 0;

	npairs++;
    }

    /* Now filter by contig pair frequency */
    if (min_freq > 0) {
	for (i = j = 0; i < npairs; i++) {
	    HashItem *hi;

	    ctg_pair[0] = pairs[i].contig[0];
	    ctg_pair[1] = pairs[i].contig[1];
	    hi = HashTableSearch(ctg_hash, (char *)ctg_pair, sizeof(ctg_pair));
	    if (!hi || hi->data.i < min_freq)
		continue;

	    pairs[j++] = pairs[i];
	}
	npairs = j;
    }

    /* Indicate end of list */
    if (npairs >= nalloc) {
	read_pair_t *new_pairs;
	nalloc = nalloc ? nalloc*2 : 1024;
	new_pairs = realloc(pairs, nalloc * sizeof(*pairs));
	if (!new_pairs) goto fail;
	pairs = new_pairs;
    }
    pairs[npairs].rec[0] = 0;
	
    HashTableIterDestroy(iter);
    HashTableDestroy(h, 0);
    pool_destroy(rp_pool);

    if (ctg_hash)
	HashTableDestroy(ctg_hash, 0);
    if (ctgs_in_list)
	HashTableDestroy(ctgs_in_list, 0);

    vmessage("Found %d read-pairs\n", npairs);

    return pairs;

 fail:
    if (iter) HashTableIterDestroy(iter);
    if (h) HashTableDestroy(h, 0);
    if (rp_pool) pool_destroy(rp_pool);
    if (ctg_hash) HashTableDestroy(ctg_hash, 0);
    if (ctgs_in_list) HashTableDestroy(ctgs_in_list, 0);
    if (pairs) free(pairs);
    return NULL;
}


/*
 * External callable functions - the C interface
 * ---------------------------------------------
 */
int find_read_pairs(GapIO *io, int num_contigs, contig_list_t *contig_array,
		    enum readpair_mode mode, int end_size, int min_mq,
		    int min_freq, tg_rec *library, int nlibrary) {
    
    read_pair_t *tarr;
    HashTable *lib_hash = NULL;
    int id;

    if (library) {
	lib_hash = create_lib_hash(library, nlibrary);
	if (NULL == lib_hash) return -1;
    }

    if (NULL == (tarr = spanning_pairs(io, num_contigs, contig_array,
				       mode, end_size, min_mq, min_freq,
				       lib_hash))) {
	if (NULL != lib_hash) HashTableDestroy(lib_hash, 0);
	return -1;
    }

    /* Find only those templates spanning multiple contigs and plot them. */
    id = PlotTempMatches(io, tarr);

    /* Tidy up */
    if (NULL != lib_hash) HashTableDestroy(lib_hash, 0);
    free(tarr);

    return id;
}

/*
 * Loads a file (already opened in fp) of find read-pairs results and
 * returns the registered ID.
 *
 * Returns -1 on failure.
 */
int csmatch_load_read_pairs(GapIO *io, FILE *fp) {
    tg_rec c1, c2;
    int pos1, pos2, end1, end2, length, mq1, mq2, n;
    tg_rec read1, read2;
    int asize = 0;
    obj_read_pair *r;
    char *val;
    int id;

    mobj_read_pair *m = calloc(1, sizeof(*m));
    if (!m)
	return -1;

    strcpy(m->tagname, CPtr2Tcl(m));
    m->num_match = 0;
    m->match = NULL;
    m->io = io;
    m->all_hidden = 0;
    m->current = -1;
    val = get_default_string(GetInterp(), gap5_defs, "READPAIR.COLOUR");
    strcpy(m->colour, val);
    m->linewidth = get_default_int(GetInterp(), gap5_defs,
				   "READPAIR.LINEWIDTH");
    m->match_type = REG_TYPE_READPAIR;
    m->reg_func = readpair_callback;

    while (11 == (n = fscanf(fp, "%"PRIrec" %d %d %"PRIrec" %d %d %d "
			    "%"PRIrec" %"PRIrec" %d %d\n",
			    &c1, &pos1, &end1, &c2, &pos2, &end2, &length,
			    &read1, &read2, &mq1, &mq2))) {
	contig_t *c;

	if (m->num_match >= asize) {
	    asize = asize ? asize*2 : 16;
	    m->match = realloc(m->match,
			       asize * sizeof(*m->match));
	    if (!m->match)
		return -1;
	}

	if (!cache_exists(io, GT_Contig, ABS(c1)) ||
	    !(c = cache_search(io, GT_Contig, ABS(c1)))) {
	    verror(ERR_WARN, "csmatch_load_read_pairs",
		   "Contig =%"PRIrec" does not exist", ABS(c1));
	    continue;
	}

	if (pos1 < c->start)
	    pos1 = c->start;
	if (end1 > c->end)
	    end1 = c->end;

	if (!cache_exists(io, GT_Contig, ABS(c2)) ||
	    !(c = cache_search(io, GT_Contig, ABS(c2)))) {
	    verror(ERR_WARN, "csmatch_load_read_pairs",
		   "Contig =%"PRIrec" does not exist", ABS(c2));
	    continue;
	}

	if (pos2 < c->start)
	    pos2 = c->start;
	if (end2 > c->end)
	    end2 = c->end;

	r = (obj_read_pair *)&m->match[m->num_match++];
	r->func =  (void *(*)(int, void *, struct obj_match_t *,
			      struct mobj_repeat_t *))readpair_obj_func;
	r->data = m;

	r->c1    = c1;
	r->c2    = c2;
	r->pos1  = pos1;
	r->pos2  = pos2;
	r->end1  = end1;
	r->end2  = end2;
	r->read1 = read1;
	r->read2 = read2;
	r->mq1   = mq1;
	r->mq2   = mq2;
	r->flags = 0; // fixme
    }
    if (n != EOF)
	verror(ERR_WARN, "csmatch_load_read_pairs",
	       "File malformatted or truncated");

    if (m->num_match) {
	id = register_id();
	contig_register(io, 0, readpair_callback, (void *)m, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ORDER | REG_GENERIC,
			REG_TYPE_READPAIR);
	update_results(io);

	return id;
    }

 err:
    if (m) {
	if (m->match)
	    free(m->match);

	free(m);
    }

    return -1;
}
