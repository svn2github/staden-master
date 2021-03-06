#ifndef _CS_OBJECT_H
#define _CS_OBJECT_H

#include <inttypes.h>
#include <tcl.h>

#include "tg_gio.h"

/*
 * ============================================================================
 * Object callback jobs
 * ============================================================================
 */

/* Returns a list of operations to display */
#define OBJ_LIST_OPERATIONS	1

/* An operation has been selected - do it */
#define OBJ_INVOKE_OPERATION	2

/* A brief one line report about this object */
#define OBJ_GET_BRIEF		3

/*
 * ============================================================================
 * Object flags
 * ============================================================================
 */

#define OBJ_FLAG_HIDDEN		(1 << 0)
#define OBJ_FLAG_VISITED	(1 << 1)
#define OBJ_FLAG_JOINED         (1 << 2)

/*
 * ============================================================================
 * Forward declarations to keep compiling happy.
 * ============================================================================
 */
struct mobj_repeat_t;
struct mobj_fij_t;

/*
 * ============================================================================
 * Object structure definitions
 * ============================================================================
 */

/*
 * An individual match structure as used by find internal joins and
 * find repeats.
 */
typedef struct obj_match_t {
    /* The following are shared between obj_match, obj_read_pair and obj_fij */
    void *(*func)(int job,
		  void *jdata,
		  struct obj_match_t *obj,
		  struct mobj_repeat_t *data);
    struct mobj_repeat_t *data;

    int inum;
    tg_rec c1, c2;
    int pos1, pos2;
    int end1, end2;
    int length;
    int flags;

    /* Below here are local to this strucure only */
    /* To pad out to equivalent length to obj_read_pair - not needed normally
     * except when trying to use our 'generic' PlotRepeats function with a
     * read pair result.
     */
    tg_rec rpos; /* This has to be the >= size of largest obj type */
    tg_rec read;
    int score;
} obj_match, obj_checkass;

/*
 * An individual match returned by find read pairs. ASSUMPTION: The start is as
 * per obj_match (see the comments within that structure). Maybe these common
 * bits should be folded into part of obj_generic. Will we really have
 * data to plot in the future that doesn't have c1, c2, (etc) ?
 */
typedef struct obj_read_pair_t {
    /* The following are shared between obj_match, obj_read_pair and obj_fij */
    void *(*func)(int job,
		  void *jdata,
		  struct obj_match_t *obj,
		  struct mobj_repeat_t *data);
    struct mobj_repeat_t *data;

    int inum;
    tg_rec c1, c2;
    int pos1, pos2;
    int end1, end2;
    int length;
    int flags;

    /* Below here are local to this strucure only */
    tg_rec read1, read2;
    uint16_t mq1, mq2;
} obj_read_pair;

/*
 * An individual match returned by find internal joins.
 * ASSUMPTION: The start is as per obj_match (see the comments within that
 * structure). Maybe these common bits should be folded into part of
 * obj_generic. Will we really have data to plot in the future that doesn't
 * have c1, c2, (etc) ?
 */
typedef struct obj_fij_t {
    /* The following are shared between obj_match, obj_read_pair and obj_fij */
    void *(*func)(int job,
		  void *jdata,
		  struct obj_fij_t *obj,
		  struct mobj_fij_t *data);
    struct mobj_fij_t *data;

    int inum;
    tg_rec c1, c2;
    int pos1, pos2;
    int end1, end2;
    int length;
    int flags;

    /* Below here are local to this strucure only */
    int score;
    int percent; /* *10000 */
    uint64_t dummy2; /* This has to be the >= size of largest obj type */
    uint32_t dummy1; /* This has to be the >= size of largest obj type */
} obj_fij;

/*
 * The generic object - a union of specific objects. All specific objects
 * contain func() and data as the first two fields. To use these in
 * conjunction with calling a job we use (eg):
 *     obj->call.func(OBJ_LIST_OPERATIONS, obj, obj->call.data);
 */
typedef union {
    struct call_t {
	void *(*func)(int job, void *jdata, void *obj, void *data);
	void *data;
    } call;

    obj_match     match;
    obj_read_pair read_pair;
    obj_fij	  fij;
} obj_generic;

/*
 * ============================================================================
 * Meta object definitions - objects maintaining lists of objects.
 * ============================================================================
 */

#define COLOUR_LEN 30

/*
 * The find repeats function - maintains a list of matches
 */
typedef struct mobj_repeat_t {
    int num_match;
    obj_match *match;
    char tagname[20];
    int linewidth;
    char colour[COLOUR_LEN];
    char *params;
    int all_hidden;
    int current;

    GapIO *io;
    int match_type;
    void (*reg_func)(GapIO *io, tg_rec contig, void *fdata,
		     reg_data *jdata);
} mobj_repeat, mobj_read_pair, mobj_find_oligo;

typedef struct mobj_fij_t {
    int num_match;
    obj_fij *match;
    char tagname[20];
    int linewidth;
    char colour[COLOUR_LEN];
    char *params;
    int all_hidden;
    int current;

    GapIO *io;
    int match_type;
    void (*reg_func)(GapIO *io, tg_rec contig, void *fdata,
		     reg_data *jdata);
    float max_mismatch;
    int min_length;
} mobj_fij;

typedef struct mobj_checkass_t {
    int num_match;
    obj_match *match;
    char tagname[20];
    int linewidth;
    char colour[COLOUR_LEN];
    char *params;
    int all_hidden;
    int current;

    GapIO *io;

    int match_type;
    void (*reg_func)(GapIO *io, tg_rec contig, void *fdata,
		     reg_data *jdata);

    /* Private to mobj_checkass only */
    int cutoffs;
} mobj_checkass;

typedef union {
    mobj_repeat     repeat;
    mobj_read_pair  read_pair;
    mobj_find_oligo find_oligo;
    mobj_fij        fij;
    mobj_checkass   checkass;
} mobj_generic;

/*
 * ============================================================================
 * Hashing for translation of tk canvas IDs to generic objects
 * ============================================================================
 */
#define HASHMODULUS 256

typedef void *HItemType;

typedef struct hashitem *HTablePtr;
typedef struct hashitem {
    int key;
    HItemType info;
    HTablePtr next;
} HTableItem;

void InitialiseHash(HTablePtr *T);
void HashInsert(HTablePtr *T, int newkey, HItemType newinfo);
HItemType *HashSearch(HTablePtr *T, int search_key);
void csmatch_reset_hash(HTablePtr T[], mobj_repeat *r);

char *obj_get_ops(int inum);
void obj_invoke_op(int inum, int op);
void obj_invoke_next(mobj_repeat *mobj);
int obj_get_next(mobj_repeat *mobj);
char *obj_get_brief(int inum);

void obj_hide(Tcl_Interp *interp, char *csplot, obj_match *obj,
	      mobj_repeat *r, HTablePtr T[]);
void obj_reveal(Tcl_Interp *interp, char *csplot, obj_match *obj,
		mobj_repeat *r, HTablePtr T[]);
void obj_remove(Tcl_Interp *interp, char *cs_plot, obj_match *obj,
		mobj_repeat *r, HTablePtr T[]);

void csmatch_join_to(GapIO *io, tg_rec contig, reg_join *j, mobj_repeat *r,
		     HTablePtr *T, char *cs_plot);
void csmatch_complement(GapIO *io, tg_rec contig, mobj_repeat *r,
			HTablePtr *T, char *cs_plot);
void csmatch_configure(GapIO *io, char *cs_plot, mobj_repeat *r);

void csmatch_remove(GapIO *io, char *cs_plot, 
		    mobj_repeat *reg_dat,
		    HTablePtr *T);
void csmatch_info(mobj_repeat *r, char *name);
void csmatch_invoke_next(mobj_repeat *r);
int csmatch_get_next(mobj_repeat *r);
void csmatch_reveal(Tcl_Interp *interp, char *cs_plot, mobj_repeat *r,
		    HTablePtr T[]);
void csmatch_hide(Tcl_Interp *interp, char *cs_plot, mobj_repeat *r,
		  HTablePtr T[]);


void csmatch_renumber(GapIO *io, tg_rec old_contig, tg_rec new_contig,
		      mobj_repeat *r, HTablePtr T[], char *cs_plot);

void csmatch_contig_delete(GapIO *io, mobj_repeat *r, tg_rec contig,
			   char *cs_plot, HTablePtr T[]);

void csmatch_replot(GapIO *io, mobj_repeat *r, HTablePtr T[], char *cs_plot);
void csmatch_reset_next(mobj_repeat *r);

int csmatch_save(mobj_generic *m, char *fn);
int csmatch_load(GapIO *io, char *fn);

#endif
