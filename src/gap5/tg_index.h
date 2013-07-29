#ifndef _TG_INDEX_H_
#define _TG_INDEX_H_

#include <stdio.h>
#include "string_alloc.h"

typedef struct {
    char *name;
    FILE *fp;
} bttmp_t;

typedef struct {
    bttmp_t *file;
    string_alloc_t *data_pool;
    char **data;
    int index;
    int size;
} bttmp_queue_t;

typedef struct {
    bttmp_queue_t *que;
    int que_size;
    int working_size;
    int write_size;
    int count;
} bttmp_sort_t;

typedef struct {
    int append;
    int no_tree;
    int merge_contigs;
    int fmt;
    char *out_fn;
    int pair_reads;
    int min_bin_size;
    int fast_mode;
    bttmp_sort_t *tmp;
    int data_type;
    int comp_mode;
    int repad;
    int store_unmapped;
    int sam_aux;
    int pair_queue;
    int store_refpos;
    int remove_dups;
    int version;
    int link_pairs;
    char *tmp_dir;
} tg_args;

#define DATA_SEQ	1
#define DATA_QUAL	2
#define DATA_NAME	4
#define DATA_ANNO	8
#define DATA_ALL	15
#define DATA_BLANK	0x100 /* not even dummy seqs */

#endif /* _TG_INDEX_H_ */
