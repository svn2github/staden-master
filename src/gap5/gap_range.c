/*
 *
 * gap_range.c - a tcl command to hold the x range for tracks.  For example
 *           template_display and depth_display
 *
 */



#include <math.h>
#include <string.h>

#include "hache_table.h"
#include "gap_range.h"
#include "gap_cli_arg.h"
#include "tg_sequence.h"

/*** Tk_ConfigSpec custom handler ***/ 

/* parse the range command to get at the data */
static int range_cmd_parse(ClientData clientData, Tcl_Interp *interp, Tk_Window tkwin,
	    	    	    char *value, char *widgRec, int offset) {
			    
    Tcl_CmdInfo info;
    
    if (!Tcl_GetCommandInfo(interp, value, &info)) return TCL_ERROR;
    
    *(gap_range_t **)(widgRec + offset) = (gap_range_t *)info.objClientData;

    return TCL_OK;
}
    

/* return value for the range command */
static char *range_cmd_print(ClientData clientData, Tk_Window tkwin,
	    	    	     char *widgRec, int offset, Tcl_FreeProc **freeProcPtr) {
    gap_range_t *gr = *(gap_range_t **)(widgRec + offset);
    
    if (gr == NULL) 
    	return "gap range not set";
    else
    	return "gap range set";
}

/* define the non-standard option types */ 
Tk_CustomOption range_option = {
    (Tk_OptionParseProc *)range_cmd_parse,
    range_cmd_print, (ClientData)NULL
};


/*** Tcl command functions ***/ 

static int grange_cmd(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[]) {
    // don't really want it to do anything
    return TCL_OK;
}

static void grange_cmd_delete(ClientData clientData) {
    gap_range_destroy((gap_range_t *)clientData);
    ckfree(clientData);
}
			

typedef struct {
    GapIO *io;
    tg_rec cnum;
} grange_arg;

/* forward declaration */
static void update_filter(gap_range_t *gr);

static int tcl_range(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[]) {
    gap_range_t *gr;
    char name[30];
    
    grange_arg args;
    cli_args a[] = {
	{"-io",     ARG_IO,  1, NULL, offsetof(grange_arg, io)},
	{"-cnum",   ARG_REC, 1, NULL, offsetof(grange_arg, cnum)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;
	
    if (NULL == (gr = (gap_range_t *)ckalloc(sizeof(gap_range_t)))) {
    	return TCL_ERROR;
    }
    
    gr->io = args.io;
    gr->crec = args.cnum;
    
    gr->r             = NULL;
    gr->tl            = NULL;
    gr->depth         = NULL;
    gr->wx0           = -1;
    gr->wx1           = -1;
    gr->nr            = 0;
    gr->template_mode = -1;
    gr->width         = -1;
    gr->ntl           = -1;
    set_filter(gr, -1, -1, -1, -1, -1, 0, 5);
    update_filter(gr);
    
    sprintf(name, "grange=%p", gr);
    
    if (NULL == Tcl_CreateObjCommand(interp, name, grange_cmd,
    					(ClientData)gr,
					(Tcl_CmdDeleteProc *)grange_cmd_delete)) {
	return TCL_ERROR;
    }
    
    Tcl_SetObjResult(interp, Tcl_NewStringObj(name, -1));
    
    return TCL_OK;
}
    
    
int GRange_Init(Tcl_Interp *interp) {
    if (NULL == Tcl_CreateObjCommand(interp, "g5::range", tcl_range,
					(ClientData)NULL,
					(Tcl_CmdDeleteProc *)NULL)) {
	return TCL_ERROR;
    }
    
    return TCL_OK;
}

/*** range functions ***/

void gap_range_destroy(gap_range_t *gr) {
    if (gr->r) {
    	free(gr->r);
    }
    
    if (gr->tl) {
    	free(gr->tl);
    }
    
    if (gr->depth) {
    	free(gr->depth);
    }
}

int gap_range_test(gap_range_t *gr) {
//    printf("GR - TEST\n");
    printf("r %p wx0 %f wx1 %f nr %d\n", gr->r, gr->wx0, gr->wx1, gr->nr);
    
    return 1;
}




static int is_filter_change(gap_range_t *gr) {
    int changed = 0;
    
    if (gr->old_filter.filter   != gr->new_filter.filter ||
    	gr->old_filter.min_qual != gr->new_filter.min_qual ||
	gr->old_filter.max_qual != gr->new_filter.max_qual ||
	gr->old_filter.c_mode   != gr->new_filter.c_mode ||
	gr->old_filter.libs_ctr != gr->new_filter.libs_ctr ||
	gr->old_filter.valid_sd != gr->new_filter.valid_sd) {
	    
    	changed = 1;
    }
    
    return changed;
}

void set_filter(gap_range_t *gr, int filter, int min, int max, int mode,
		int accuracy, int libs_ctr, double valid_sd) {

    gr->new_filter.filter   = filter;
    gr->new_filter.min_qual = min;
    gr->new_filter.max_qual = max;
    gr->new_filter.c_mode   = mode;
    gr->new_filter.accuracy = accuracy;
    gr->new_filter.libs_ctr = libs_ctr;
    gr->new_filter.valid_sd = valid_sd;
}


static void update_filter(gap_range_t *gr) {

    gr->old_filter.filter   = gr->new_filter.filter;
    gr->old_filter.min_qual = gr->new_filter.min_qual;
    gr->old_filter.max_qual = gr->new_filter.max_qual;
    gr->old_filter.c_mode   = gr->new_filter.c_mode;
    gr->old_filter.accuracy = gr->new_filter.accuracy;
    gr->old_filter.libs_ctr = gr->new_filter.libs_ctr;
    gr->old_filter.valid_sd = gr->new_filter.valid_sd;
}
    


/*
   recalculate the sequence range if neccessary
   if out of memory, returns 1 with a null gr->r pointer
   returns 1 on change, 0 otherwise
*/
int gap_range_recalculate(gap_range_t *gr, int width, double new_wx0,
			  double new_wx1, int new_mode, int force) {
    int changed = 0;
    contig_t *c;

    /* only change if needed or if force is true */
    if (force || gr->r == NULL
	|| new_wx0 != gr->wx0 || new_wx1 != gr->wx1
	|| new_mode != gr->template_mode || gr->width != width
	|| gr->new_filter.accuracy != gr->old_filter.accuracy) {

	c = cache_search(gr->io, GT_Contig, gr->crec);
	if (NULL == c) goto fail;
	cache_incr(gr->io, c);

    	if (gr->r) free(gr->r);
	gr->r = contig_seqs_in_range(gr->io, &c, new_wx0, new_wx1,
				     new_mode, &gr->nr);
	cache_decr(gr->io, c);
	
	if (gr->r) {
	    tline *tmp_line = NULL;
	    gap_depth_t *new_depth = NULL;
	    size_t nlines = gr->nr > 0 ? gr->nr : 1;

	    gr->wx0  = new_wx0;
	    gr->wx1  = new_wx1;
	    gr->template_mode = new_mode;
	    update_filter(gr);
	    
	    tmp_line = realloc(gr->tl, nlines * sizeof(tline));
	    if (NULL == tmp_line) goto fail;
	    gr->tl = tmp_line;
	    
	    if (gr->width != width) {
		int nlibs = gr->io->db->Nlibraries;
		size_t w = width > 0 ? width : 1;
		w *= nlibs + 1;
		new_depth = realloc(gr->depth, w * sizeof(gap_depth_t));
		if (NULL == new_depth) goto fail;
	    	gr->width = width;
		gr->depth = new_depth;
	    }
	}
	
	changed = 1;
    }
    
    return changed;
 fail:
    if (gr->r) free(gr->r);
    gr->r = NULL;
    gr->nr = 0;
    return 1;
}

int gap_range_x(gap_range_t *gr, double ax_conv, double bx_conv, 
    	    	int fwd_col, int rev_col, int single_col, int *span_col,
		int inconsistent_col, int force, int reads_only,
		HashTable *lib_recs) {
    int i, j;
    double max_height = 0;
    int lib_type, lib_idx;
    HashTable *h, *rec_h;
    contig_t *c;
    int valid_pair, nlibs;
    double valid_sd;

    if (!is_filter_change(gr) && !force)
	return gr->ntl;

    if (!(c = cache_search(gr->io, GT_Contig, gr->crec)))
	return 0;
    cache_incr(gr->io, c);

    update_filter(gr);
    gr->ntl = 0;
    nlibs = gr->io->db->Nlibraries;
    memset(gr->depth, 0, gr->width * sizeof(gap_depth_t) * (nlibs+1));
    valid_sd = gr->new_filter.valid_sd;

    h = HashTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    rec_h = HashTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);

    /* Build a mapping from library records to 1..N index */
    for (i = 0; i < gr->io->db->Nlibraries; i++) {
	HashData hd;
	hd.i = i+1;
	HashTableAdd(rec_h, (char *)arrp(tg_rec, gr->io->library, i),
		     sizeof(tg_rec), hd, NULL);
    }
    

    for (i = 0; i < gr->nr; i++) {
	int sta, end;
	int r_sta, r_end;
	int span   = 0;
	int single = 0;
	int col;
	double mq;
	rangec_t *r  = &gr->r[i];
	tline    *tl = &gr->tl[gr->ntl];
	library_t *lib;

	if (lib_recs) {
	    if (!HashTableSearch(lib_recs, (char *)r->library_rec,
				 sizeof(r->library_rec)))
		continue;
	}
	    
	sta = r->start;
	end = r->end;
	mq  = r->mqual;
	    
	/* depth plot setting (sequence) */
	if (!(mq < gr->new_filter.min_qual || mq > gr->new_filter.max_qual)) {
	    double inc;

	    r_sta = (sta - bx_conv) * ax_conv; //world to raster conversion
	    r_end = (end - bx_conv) * ax_conv;

	    inc = ((end - sta) * ax_conv) / (r_end - r_sta + 1);

	    if (r_sta < 0)  	   r_sta = 0;
	    if (r_end >= gr->width) r_end = gr->width - 1;

	    for (j = r_sta; j <= r_end; j++) {
		gr->depth[j].s += inc;
		if (max_height < gr->depth[j].s)
		    max_height = gr->depth[j].s;
	    }
	}
	    
	/* deal with filters */
	if (gr->new_filter.filter & (FILTER_PAIRED | FILTER_SPANNING) && !r->pair_rec) continue;
		
	if (gr->new_filter.filter & FILTER_SINGLE && r->pair_rec) continue;

	tl->t_strand = ((r->flags & GRANGE_FLAG_END_REV) != 0) == ((r->flags & GRANGE_FLAG_COMP1) != 0);
	tl->rec = r->rec;
	    
	if (!reads_only && r->pair_rec) {
	    /* only draw once */
	    if (r->flags & (1<<30)) {
		r->flags &= ~(1<<30);
		continue;
	    }
		
	    if (r->pair_ind >= 0 && r->pair_ind < gr->nr) {
		gr->r[r->pair_ind].flags |= (1<<30);
	    }
		
	    /* accurate drawing, this will slow things down */
	    if (r->pair_timestamp < gr->io->db->timestamp) {
		if (gr->new_filter.accuracy) {
		    /* Pair was off-screen, so get more accurate results */
		    sequence_get_range_pair_position(gr->io, r, gr->crec, 0);
		} else if (r->pair_ind == -1 && r->pair_contig == gr->crec) {
		    /* Same contig, but out of date. Check boundary */
		    //if (r->pair_start > c->end || r->pair_end < c->start)
		    r->pair_contig = 0;
		}
	    }

	    //span = (gr->crec != r->pair_contig) ? r->pair_contig : 0;
	    span = (gr->crec != r->pair_contig);
		
	    if (gr->new_filter.filter & FILTER_SPANNING && !span) continue;
		
	    if (gr->new_filter.filter & FILTER_NONSPANNING && span) continue;

	    if (!span) {
		sta = r->start < r->pair_start ? r->start : r->pair_start;
		end = r->end > r->pair_end ? r->end : r->pair_end;
		r->flags |= GRANGE_FLAG_CONTIG;
		    
		switch (gr->new_filter.c_mode) {
		case 0:
		case 3:
		    mq = (r->mqual + r->pair_mqual) / 2.0;
		    break;
		case 1:
		    mq = r->mqual < r->pair_mqual ? r->mqual : r->pair_mqual;
		    break;
		case 2:
		    mq = r->mqual < r->pair_mqual ? r->pair_mqual : r->mqual;
		    break;
		}
	    }
	} else {
	    mq = r->mqual;
	    if (!reads_only) single = 1;
	}
	    
	/* now for some colour settings */
	if (mq < 0) mq = 0;
	if (mq > 255/3) mq = 255/3;
	    
	tl->mq = mq;
	col = (int)(mq*3 / 8);
	    
	if (single) col = single_col;
	if (span && span_col) 	col = span_col[(r->pair_contig * 15551) % 10];

	if (r->library_rec) {
	    HashItem *hi;

	    if ((hi = HashTableSearch(h, (char *)&r->library_rec,
				       sizeof(tg_rec)))) {
		lib = (library_t *)hi->data.p;
	    } else {
		HashData hd;

		lib = cache_search(gr->io, GT_Library, r->library_rec);
		update_library_stats(gr->io, lib->rec, 100, NULL, NULL, NULL);
		hd.p = lib;
		cache_incr(gr->io, lib);

		HashTableAdd(h, (char *)&r->library_rec, sizeof(tg_rec),
			     hd, NULL);
	    }
	    lib_type = lib->lib_type;
	    if ((hi = HashTableSearch(rec_h, (char *)&r->library_rec,
				      sizeof(tg_rec)))) {
		lib_idx = hi->data.i;
	    } else {
		lib_idx = 0;
	    }
	} else {
	    lib = NULL;
	    lib_type = LIB_T_INWARD;
	    lib_idx = 0;
	}
	    
	/* Check consistency */
	switch (lib_type) {
	case LIB_T_INWARD:
	    if (r->pair_rec && r->flags & GRANGE_FLAG_CONTIG) {
		if (!((((r->start < r->pair_start)/*^r->comp*/) &&
		       ((r->flags & GRANGE_FLAG_COMP1)) == 0 &&
		       ((r->flags & GRANGE_FLAG_COMP2)) != 0) ||
		      (((r->pair_start < r->start)/*^r->comp*/) &&
		       ((r->flags & GRANGE_FLAG_COMP1)) != 0 &&
		       ((r->flags & GRANGE_FLAG_COMP2)) == 0)))
		    col = inconsistent_col;
	    }
	    break;

	case LIB_T_OUTWARD:
	    if (r->pair_rec && r->flags & GRANGE_FLAG_CONTIG) {
		if (!((((r->start > r->pair_start)/*^r->comp*/) &&
		       ((r->flags & GRANGE_FLAG_COMP1)) == 0 &&
		       ((r->flags & GRANGE_FLAG_COMP2)) != 0) ||
		      (((r->pair_start > r->start)/*^r->comp*/) &&
		       ((r->flags & GRANGE_FLAG_COMP1)) != 0 &&
		       ((r->flags & GRANGE_FLAG_COMP2)) == 0)))
		    col = inconsistent_col;
	    }
	    break;

	case LIB_T_SAME:
	    if (r->pair_rec && r->flags & GRANGE_FLAG_CONTIG) {
		if (!((((r->flags & GRANGE_FLAG_COMP1)) == 0) ==
		      (((r->flags & GRANGE_FLAG_COMP2)) == 0)))
		    col = inconsistent_col;
	    }
	    break;

	default:
	    fprintf(stderr, "Unexpected lib_type %d\n", lib_type);
	}

	valid_pair = (!single && !span && col != inconsistent_col) ? 1 : 0;
	if (valid_pair) {
	    int isize_min = lib
		? lib->insert_size[lib_type] - valid_sd*lib->sd[lib_type]
		: INT_MAX;
	    int isize_max = lib
		? lib->insert_size[lib_type] + valid_sd*lib->sd[lib_type]
		: INT_MIN;
	    if (end - sta > isize_max || end-sta < isize_min)
		valid_pair = 0;
	}
	    
	if (gr->new_filter.filter & FILTER_CONSISTENT && col == inconsistent_col) continue;
	    
	if (gr->new_filter.filter & FILTER_INCONSISTENT && col != inconsistent_col) continue;
	    
	/* Filter */
	if (tl->mq < gr->new_filter.min_qual || tl->mq > gr->new_filter.max_qual) continue;

	/* depth plot template pairs */
	if (valid_pair) {
	    double inc, inc2;

	    r_sta = (sta - bx_conv) * ax_conv; // world to raster conversion
	    r_end = (end - bx_conv) * ax_conv;
		
	    inc = ((end - sta) * ax_conv) / (r_end - r_sta + 1);

	    if (r_sta < 0) r_sta = 0;
	    if (r_end >= gr->width) r_end = gr->width - 1;

	    if (lib_idx) {
		inc2 = inc * (end-sta+1);
		int J = lib_idx*gr->width;
		for (j = r_sta, J += r_sta; j <= r_end; j++, J++) {
		    gr->depth[j].t += inc; // aggregate
		    gr->depth[J].s += inc; // per library
		    gr->depth[J].t += inc2;
		    if (max_height < gr->depth[j].t)
			max_height = gr->depth[j].t;
		}
	    } else {
		for (j = r_sta; j <= r_end; j++) {
		    gr->depth[j].t += inc;
		    if (max_height < gr->depth[j].t)
			max_height = gr->depth[j].t;
		}
	    }
	}
    
	/* Generate line data */
	tl->col[0] = tl->col[1] = tl->col[2] = col;
	tl->x[0]   = sta;
	/* range in col 0 */
	tl->x[1]   = sta;
	/* range in col 1 */
	tl->x[2]   = end; 
	/* range in col 3 */
	tl->x[3]   = end;
	    
	/* Readings too? */
	if (gr->new_filter.c_mode == 3) {
	    
	    col = (r->flags & GRANGE_FLAG_END_MASK)
		== GRANGE_FLAG_END_FWD ? fwd_col : rev_col;

	    if (r->start == tl->x[0]) {
		tl->x[1]   = r->end;
		tl->col[0] = col;
	    } else if (r->end == tl->x[3]) {
		tl->x[2]   = r->start;
		tl->col[2] = col;
	    } else if (!span && r->pair_rec &&
		       r->pair_start < r->start && r->pair_end > r->end) {
		/* Containment, so colour middle */
		tl->x[1] = r->start;
		tl->x[2] = r->end;
		tl->col[1] = col;
	    } else {
		verror(ERR_WARN, "gap_range_x",
		       "error, start/end do not match template pos (single)");
		verror(ERR_WARN, "gap_range_x", "start %d/%d end %d/%d",
		       r->start, tl->x[0], r->end, tl->x[3]); 
	    }
		
	    if (!span && r->pair_rec && (r->pair_start || r->pair_end)) {
		
		col = (r->flags & GRANGE_FLAG_PEND_MASK)
		    == GRANGE_FLAG_PEND_FWD ? fwd_col : rev_col;

		if (r->pair_start == tl->x[0]) {
		    tl->x[1]   = r->pair_end;
		    tl->col[0] = col;
		} else if (r->pair_end == tl->x[3]) {
		    tl->x[2]   = r->pair_start;
		    tl->col[2] = col;
		} else if (r->pair_start > r->start && r->pair_end < r->end) {
		    /* Containment, so colour middle */
		    tl->x[1] = r->pair_start;
		    tl->x[2] = r->pair_end;
		    tl->col[1] = col;
		} else {
		    verror(ERR_WARN, "gap_range_x",
			   "error, start/end do not match template pos (pair)");
		    verror(ERR_WARN, "gap_range_x",
			   "start %d/%d end %d/%d",
			   r->pair_start, tl->x[0], r->pair_end, tl->x[3]);	
		}
	    }
		
	}

	gr->ntl++;
    }

    cache_decr(gr->io, c);

    if (h)
	HashTableDestroy(h, 0);
    if (rec_h)
	HashTableDestroy(rec_h, 0);

    gr->max_height = ceil(max_height);
    gr->max_height++;

    return gr->ntl;
}



/*
    reset to starting values, free any memory
    leave the gap and contig values
*/
void gap_range_reset(gap_range_t *gr) {
    
    if (gr->r) {
    	free(gr->r);
    }

    gr->r             = NULL;
    gr->wx0           = -1;
    gr->wx1           = -1;
    gr->nr            = 0;
    gr->template_mode = -1;
}










































