#ifndef _TG_FUNC_H_
#define _TG_FUNC_H_

#include <tg_index.h>

#include "string_alloc.h"

typedef struct {
    char *tname;
    tg_rec rec;
    tg_rec bin;
    int idx;
    tg_rec crec;
    int pos;
    int orient;
    int flags;
    int mq;
    int len;
} pair_loc_t;

typedef struct {
    bttmp_t *file;
    pair_loc_t *pair;
    string_alloc_t *name_pool;
    int index;
    int pair_size;
    int mq;
    int len;
} pair_queue_t;

typedef struct {
    pair_queue_t *que;
    int que_size;
    int working_size;
    int write_size;
    int count;
    HacheTable *phache;
    bttmp_t *finish;
    tg_rec max_bin;
} tg_pair_t;
    

int bttmp_build_index(GapIO *io, bttmp_store_t *bs, long work_size, long group_size);
bttmp_store_t *bttmp_store_initialise(long write_sz);
void bttmp_store_delete(bttmp_store_t *bs);

tg_rec save_sequence(GapIO *io, seq_t *seq, bin_index_t *bin, range_t *r_out);

tg_rec save_range_sequence(GapIO *io, seq_t *seq, uint8_t mapping_qual,
			   tg_pair_t *pair, int is_pair, char *tname,
			   contig_t *c, tg_args *a, int flags, library_t *lib,
			   tg_rec *bin_rec);

void create_new_contig(GapIO *io, contig_t **c, char *cname, int merge);

/*
 * Turns a comma separated list of data types into a bit-field.
 */
int parse_data_type(char *type);

int tg_index_file_type (char *fn);

void unescape_line(char *txt);

tg_pair_t *create_pair(int queue);
void finish_pairs(GapIO *io, tg_pair_t *pair, int link_pairs);
void delete_pair(tg_pair_t *pair);
long tg_get_line(char **line, long *length, FILE *fp);

#endif


