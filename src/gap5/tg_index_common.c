/*
 * tg_index_common.c - common functions for tg_index
 *
 * Andrew Whitwham, October 2009
 * Wellcome Trust Sanger Institute
 *
 */

#include <staden_config.h>

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <math.h>

#include "tg_gio.h"
#include "tg_index_common.h"
#include "zfio.h"
#include "break_contig.h"  /* For contig_visible_{start,end} */

/* --------------------------------------------------------------------------
 * Temporary file handling for storing name + record.
 * This is used in the name B+Tree generation code as it's more efficient
 * to delay generation of the B+Tree until after adding all the sequence
 * records.
 */
 
static char *get_tmp_directory(void) {
    
    char *dir;
    
    /* Find a place to put the tmp files    
       First *nix then Windows */
    if (NULL == (dir = getenv("TMPDIR"))) {
    	if (NULL == (dir = getenv("TMP_DIR"))) {
    	    dir = getenv("TEMP");
	}
    }
    
    /* if tmp_dir is null then we will just use the default locations */
    
    return dir;
}


static void c_replace(char *instr, char old, char new) {
    int i;
    size_t len = strlen(instr);
    
    for (i = 0; i < len; i++) {
    	if (instr[i] == old) {
	    instr[i] = new;
	}
    }
}
    

static char *create_tmp_file_name(char *file_name) {
    char *name;
    char *dir;
    
    dir = get_tmp_directory();
    
    if (dir) {
    	char *start;
	int s_len;
    
	c_replace(file_name, '\\', '/');

	if (NULL == (start = strrchr(file_name, '/'))) {
    	    start = file_name;
	} else {
    	    start++;
	}

	s_len = strlen(dir) + strlen(start) + 2; // 2 for terminator and seperator
	name = (char *)malloc(s_len);
	sprintf(name, "%s/%s", dir, start);
    } else {
    	name = (char *)malloc(strlen(file_name) + 1);
	strcpy(name, file_name);
    }
    
    return name;
}


bttmp_t *bttmp_file_open(void) {
    bttmp_t *tmp;
    int fd;
    char file_name[L_tmpnam];

    /*
     * This emits a warning from gcc:
     *     the use of `tmpnam' is dangerous, better use `mkstemp'
     *
     * The problem is, mkstemp isn't standard while tmpnam is. Instead
     * we use tmpnam in a safe manner via open with O_CREAT|O_EXCL and then
     * convert this to a FILE pointer. This is basically what tmpfile does, 
     * but in our case we need to know the file name too so we can sort
     * it later on.
     */

    if (NULL == tmpnam(file_name)) {
    	perror("tmpnam()");
	return NULL;
    }
    
    if (NULL == (tmp = malloc(sizeof(bttmp_t)))) {
    	fprintf(stderr, "Error: unable to allocate memory for tmp_file struct\n");
	return NULL;
    }
    
    tmp->name = create_tmp_file_name(file_name);
    
    if (tmp->name == NULL) {
    	fprintf(stderr, "Error: unable to create tmp file name.\n");
    	free(tmp);
	return NULL;
    }
    
    if (-1 == (fd = open(tmp->name, O_RDWR|O_CREAT|O_EXCL, 0666))) {
	perror(tmp->name);
	free(tmp->name);
	free(tmp);
	return NULL;
    }

    if (NULL == (tmp->fp = fdopen(fd, "wb+"))) {
	perror(tmp->name);
	free(tmp->name);
	free(tmp);
	return NULL;
    }

    return tmp;
}

void bttmp_file_close(bttmp_t *tmp) {
    if (tmp && tmp->name) {
    	if (tmp->fp) {
	    fclose(tmp->fp);
	}
	
	unlink(tmp->name);
	free(tmp->name);
	free(tmp);
    }
}

static void bttmp_initialise_data(bttmp_data_t *d, long sz) {
    d->data_pool = string_pool_create(1024 * 1024);
    d->data  = malloc(sz * sizeof(char *));
    d->index = 0;
}


static void bttmp_delete_data(bttmp_data_t *d) {
    string_pool_destroy(d->data_pool);
    free(d->data);
}    


/* Sort the temporary file, and rewind to start */
void bttmp_file_sort(bttmp_t *tmp) {
    char new_tmp[L_tmpnam];
    char buf[100+2*L_tmpnam];

    if (!tmpnam(new_tmp)) {
	verror(ERR_WARN, "bttmp_file_sort",
	       "Failed to find a temporary file name.\n");
	return;
    }
    sprintf(buf, "sort < %s > %s", tmp->name, new_tmp);
    fclose(tmp->fp);

    /* Use unix sort for now */
    printf("buf=%s\n", buf);
    if (-1 == system(buf))
	perror(buf);
    printf("done\n");

    // unlink(tmp->name);
    strcpy(tmp->name, new_tmp);
    tmp->fp = fopen(tmp->name, "rb+");
}

/*
 * Repeatedly fetch lines from the temp file.
 * NB: non-reentrant. Value is valid only until the next call to this
 * function.
 *
 * Return name on success and fills out rec.
 *       NULL on EOF (*rec==0) or failure (*rec==1)
 */
char *bttmp_file_get(bttmp_t *tmp, tg_rec *rec) {
    static char line[8192];
    int64_t recno;

    if (!tmp->fp) {
	*rec = 1;
	return NULL;
    }

    if (fscanf(tmp->fp, "%s %"PRId64"\n", line, &recno) == 2) {
	*rec = recno;
	return line;
    }

    *rec = feof(tmp->fp) ? 0 : 1;
	
    return NULL;
}


bttmp_queue_t *bttmp_add_queue(bttmp_sort_t *bs, bttmp_t *f) {
    bs->que[bs->index].file = f;
    bs->que[bs->index].data_pool = NULL;
    bs->que[bs->index].data  = malloc(bs->working_size * sizeof(char *));
    bs->que[bs->index].index = 0;
    bs->que[bs->index].size  = bs->working_size;
    
    bs->que[bs->index].file->fp = fopen(bs->que[bs->index].file->name, "r");

    return &bs->que[bs->index++];
}


static int cmpstringp(const void *p1, const void *p2) {
    return strcmp(* (char * const *) p1, * (char * const *) p2);
}


static void bttmp_save_file(bttmp_store_t *s, int close) {
    int i;
    
    qsort(s->data.data, s->data.index, sizeof(s->data.data[0]), cmpstringp);
    
    s->files[s->file_no] = bttmp_file_open();
    
    for (i = 0; i < s->data.index; i++) {
    	fprintf(s->files[s->file_no]->fp, "%s\n", s->data.data[i]);
    }
    
    string_pool_destroy(s->data.data_pool);
    free(s->data.data);
    if (close) {
	fclose(s->files[s->file_no]->fp);
    } else {
	fflush(s->files[s->file_no]->fp);
	rewind(s->files[s->file_no]->fp);
    }
}


static void bttmp_add_file(bttmp_store_t *s) {
    s->file_no++;
    
    if (s->file_no == s->file_grow) {
    	s->file_grow += s->file_grow;
	s->files = realloc(s->files, s->file_grow * sizeof(bttmp_t *));
    }
    
    bttmp_initialise_data(&s->data, s->write_size);
}


/*
 * Stores a name and record in a temporary file suitable for sorting and
 * adding to the name index at a later stage.
 */
static void bttmp_file_store(bttmp_store_t *tmp,  size_t name_len, char *name, tg_rec rec) {
    char entry[1024];
    
    sprintf(entry, "%.*s %"PRIrec"", (int)name_len, name, rec);
    
    tmp->data.data[tmp->data.index++] = string_dup(tmp->data.data_pool, entry);
    
    if (tmp->data.index == tmp->write_size) {
    	bttmp_save_file(tmp, 1);
	bttmp_add_file(tmp);
    }
}


// must find a better name for it than this
bttmp_store_t *bttmp_store_initialise(long write_sz) {
    bttmp_store_t *tmp = malloc(sizeof(bttmp_store_t));
    
    if (tmp == NULL) {
    	fprintf(stderr, "Error: unable to malloc bttmp_store_t\n");
	return NULL;
    }
    
    tmp->write_size = write_sz;
    tmp->file_no = 0;
    tmp->file_grow = 1000;

    bttmp_initialise_data(&tmp->data, write_sz);
    tmp->files = malloc(tmp->file_grow * sizeof(bttmp_t *));
     
    return tmp;
}


void bttmp_store_delete(bttmp_store_t *s) {
    free(s->files);
    free(s);
}   


bttmp_sort_t *bttmp_sort_initialise(long group_size, long work_size) {
    bttmp_sort_t *sort = malloc(sizeof(bttmp_sort_t));
    sort->que          = calloc(group_size, sizeof(bttmp_queue_t));
    sort->que_size     = group_size;
    sort->working_size = work_size;
    sort->index        = 0;
    
    return sort;
}


static int bttmp_load_data(bttmp_queue_t *bq) {
    int i;
    char *line_in = NULL;
    long line_size = 0;
    
    if (bq->data_pool) {
    	string_pool_destroy(bq->data_pool);
    }
    
    bq->data_pool = string_pool_create(1024);
    
    for (i = 0; i < bq->size; i++) {
    	int ret = tg_get_line(&line_in, &line_size, bq->file->fp);
	
	if (ret <= 0) break;
	
	bq->data[i] = string_dup(bq->data_pool, line_in);
    }
    
    bq->size = i;
    bq->index = 0;
    
    free(line_in);
    
    return bq->size;
}


static void bttmp_get_next_entry(bttmp_queue_t *que) {
    que->index++;
    
    if (que->index == que->size) {
    	bttmp_load_data(que);
    }
}


// a binary tree sort to speed up the merging
typedef struct sort_node_s {
    struct sort_node_s *up;
    struct sort_node_s *child_left;
    struct sort_node_s *child_right;
    bttmp_queue_t      *data;
} sort_node;


sort_node *new_sort_node(sort_node *up, sort_node *left, sort_node *right) {
    sort_node *leaf = malloc(sizeof(sort_node));
    leaf->up = up;
    leaf->child_left = left;
    leaf->child_right = right;
    leaf->data = NULL;
    
    return leaf;
}


sort_node *add_sort_leaf(sort_node *last_leaf, bttmp_queue_t *val) {
    int tier = 0;
    int not_found = 1;
    
    // new tree
    if (last_leaf == NULL) {
	last_leaf = new_sort_node(NULL, NULL, NULL);
	last_leaf->data = val;
	return last_leaf;
    }
    
    // find where the new leaf should go
    while (not_found) {
	if (tier && last_leaf->child_left == NULL) { // create new left child
	    sort_node *leaf = new_sort_node(last_leaf, NULL, NULL);
	    last_leaf->child_left = leaf;
	    tier--;
	    last_leaf = leaf;

	    if (tier == 0) not_found = 0;
	} else if (tier && last_leaf->child_right == NULL) { // create a right child
	    sort_node *leaf = new_sort_node(last_leaf, NULL, NULL);
	    last_leaf->child_right = leaf;
	    tier--;
	    last_leaf = leaf;

	    if (tier == 0) not_found = 0;
	} else { // go up
	    if (last_leaf->up) {
	    	last_leaf = last_leaf->up;
		tier++;
	    } else { // new top
	    	sort_node *leaf = new_sort_node(NULL, last_leaf, NULL);
	    	last_leaf->up = leaf;
		last_leaf = leaf;
		tier++;
	    }
	}
    }
    
    last_leaf->data = val;
    
    return last_leaf;
}


sort_node *sort_tree_head(sort_node *leaf) {
    
    while (leaf->up) {
    	leaf = leaf->up;
    }
    
    return leaf;
}


void populate_sort_tree(sort_node *node) {
    if (node->child_left) {
    	populate_sort_tree(node->child_left);
    }
    
    if (node->child_right) {
    	populate_sort_tree(node->child_right);
    }
    
    if (node->child_left && node->child_right) {
    	bttmp_queue_t *left  = node->child_left->data; 
    	bttmp_queue_t *right = node->child_right->data;
	
    	node->data = left;
	
	if (right->size && strcmp(left->data[left->index], right->data[right->index]) > 0) {
	   node->data = right;
	}
    }
    
    return;
}


sort_node *delete_sort_tree(sort_node *node) {
    sort_node *child = NULL;
    
    if (node->child_left) {
    	child = delete_sort_tree(node->child_left);
	if (child) free(child);
    }
    
    if (node->child_right) {
    	child = delete_sort_tree(node->child_right);
	if (child) free(child);
    }
    
    return node;
}


static bttmp_t *bttmp_merge_sort(bttmp_sort_t *sort) {
    int i;
    int next;
    bttmp_t *output = bttmp_file_open();
    sort_node *head, *last_leaf = NULL;
    bttmp_queue_t dummy;
    
    for (i = 0; i < sort->index; i++) {
    	bttmp_load_data(&sort->que[i]);
	last_leaf = add_sort_leaf(last_leaf, &sort->que[i]);
	sort->que[i].node = last_leaf;
    }
    
    next = pow(2, ceil(log(sort->index) / log(2)));
    
    if (sort->index < next) { // balance the tree
    	dummy.size = 0;
	
    	for (i = 0; i < (next - sort->index); i++) {
   	    last_leaf = add_sort_leaf(last_leaf, &dummy);
	}
    }
    
    head = sort_tree_head(last_leaf);
    populate_sort_tree(head);
    
    // do the actual sorting
    while (head->data->size) {
	sort_node *node = head->data->node;
    	sort_node *left;
	sort_node *right;
    	bttmp_queue_t *que = head->data;
	
	fprintf(output->fp, "%s", que->data[que->index]);
	bttmp_get_next_entry(que);
	
	// redo tree with new value
	while (node->up) {
	    node = node->up;
	    left  = node->child_left;
	    right = node->child_right;
	    
	    node->data = left->data;
	    
	    if ((!node->data->size) || (right->data->size &&  
	    	strcmp(node->data->data[node->data->index], right->data->data[right->data->index]) > 0)) {
		node->data = right->data;
	    }
	    
	}
    }

    delete_sort_tree(head);
    free(head);
    printf("Partial done %s\n", output->name);
    rewind(output->fp);
    
    return output;
}


static void bttmp_reset_sort(bttmp_sort_t *s) {
    int i;
    
    for (i = 0; i < s->index; i++) {
	bttmp_file_close(s->que[i].file);
	s->que[i].index = 0;
	s->que[i].size  = s->working_size;
	free(s->que[i].data);
	s->que[i].data = NULL;
    }
    
    s->index = 0;
}


static long bttmp_write_index(GapIO *io, FILE *fp) {
    char *line_in  = NULL;
    long line_size = 0;
    long num_found = 0;
    
    puts("Building index: one dot per 10k reads");

    while (0 < tg_get_line(&line_in, &line_size, fp)) {
	char name[1024];
	tg_rec rec;
	int64_t recno;

	sscanf(line_in, "%s %"PRId64"\n", name, &recno);
	rec = recno;

	sequence_index_update(io, name, strlen(name), rec);

	num_found++;

	if (!(num_found % 10000)) {
	    putchar('.'); fflush(stdout);
	    cache_flush(io);
	}
    }

    cache_flush(io);
    free(line_in);
    
    printf("\nIndexed %ld reads\n", num_found);
    
    return num_found;
}
    

static void bttmp_sort_delete(bttmp_sort_t *bs) {
    int i;
    
    for (i = 0; i < bs->que_size; i++) {
	if (bs->que[i].data_pool) string_pool_destroy(bs->que[i].data_pool);
	if (bs->que[i].data) free(bs->que[i].data);
    }
    
    if (bs->que) free(bs->que);
    
    free(bs);
}
   

int bttmp_build_index(GapIO *io, bttmp_store_t *bs, long work_size, long group_size) {
    int round = 0;
    bttmp_sort_t *sort = bttmp_sort_initialise(group_size, work_size);
    
    bttmp_save_file(bs, bs->file_no); // save the last unfinished file
    bs->file_no++;
    
    printf("Sorting read names...\n");
    
    while (bs->file_no > 1) {
    	int i;
	long new_file_ind = 0;
	long group = 0;
	bttmp_t **new_files = malloc(((bs->file_no / group_size) + 1) * sizeof(bttmp_t *));
	
	for (i = 0; i < bs->file_no; i++) {
	    bttmp_add_queue(sort, bs->files[i]);
	    group++;
	    
	    if (group == group_size) {
	    	new_files[new_file_ind++] = bttmp_merge_sort(sort);
		bttmp_reset_sort(sort);
		group = 0;
	    }
	}
	
	if (group) { // incomplete group to sort
	    new_files[new_file_ind++] = bttmp_merge_sort(sort);
    	    bttmp_reset_sort(sort);
    	}
	
	free(bs->files);
	bs->files = new_files;
	bs->file_no = new_file_ind;
	
	round++;
	printf("...sort round %d done\n", round);
    }
    
    printf("Sorting done.\n");
        
    // should be left with one file
    bttmp_write_index(io, bs->files[0]->fp);
    bttmp_file_close(bs->files[0]);
    bttmp_sort_delete(sort);
    return 0;
}
    

#if 0
/* debugging functions */
static void fprint_pair(FILE *tp, char *name, pair_loc_t *p) {
    fprintf(tp, "%s %"PRIrec" %"PRIrec" %d %"PRIrec" %d\n",
	    name, p->rec, p->bin, p->idx, p->crec, p->pos);
} 

static void print_pair(pair_loc_t *p) {
    fprintf(stderr, "rec:%"PRIrec" bin:%"PRIrec" idx:%d crec:%"PRIrec
	    " pos:%d\n",
	    p->rec, p->bin, p->idx, p->crec, p->pos);
} 

static void print_range(range_t *r) {
    fprintf(stderr, "start:%d end:%d rec:%"PRIrec" mqual:%d pair_rec:%"PRIrec
	    " flags:%d\n",
	    r->start, r->end, r->rec, r->mqual, r->pair_rec, r->flags);
}
#endif

/* --------------------------------------------------------------------------
 * Read-pair and sequence storing functions. Common to all file format
 * parsers.
 */

/* save sequence, returns recno */
tg_rec save_sequence(GapIO *io, seq_t *seq, bin_index_t *bin, range_t *r_out) {

    seq->bin = bin->rec;
    seq->bin_index = r_out - ArrayBase(range_t, bin->rng);
    
    return sequence_new_from(io, seq);
}


typedef struct {
    char *tname;
    char *data;
} pair_data;


static int cmp_pair(const void *p1, const void *p2) {
    pair_data *first = (pair_data *)p1; 
    pair_data *second = (pair_data *)p2;
    
    return strcmp(first->tname, second->tname);
}


static char *pair_to_pooled_string(string_alloc_t *s_pool, pair_loc_t *p) {
    char holder[255];
    
    sprintf(holder,
	    "%"PRIrec" %"PRIrec" %d %"PRIrec" %d %d %d %d %d",
    	    p->rec, p->bin, p->idx, p->crec, p->pos, p->orient, p->flags,
	    p->len, p->mq);
	    
    return string_dup(s_pool, holder);
}


static pair_queue_t *add_pair_queue(tg_pair_t *pair) {
    pair_queue_t *tmp;
    int cur = pair->que_size;
    
    pair->que_size++;
    
    tmp = (pair_queue_t *)realloc(pair->que, sizeof(pair_queue_t) * pair->que_size);
    
    if (tmp) {
    	pair->que = tmp;
    } else {
    	return NULL;
    }
    
    if (NULL == (pair->que[cur].file = bttmp_file_open())) {
    	fprintf(stderr, "Cannot open tmp file in add_pair_queue\n");
    	return NULL;
    }
    
    pair->que[cur].pair      = NULL;
    pair->que[cur].name_pool = NULL;
    pair->que[cur].index     = 0;
    pair->que[cur].pair_size = 0;
    
    fprintf(stderr, "New queue added, no. %d\n", pair->que_size);
    
    return &pair->que[cur];
}


static void save_pair_data(tg_pair_t *pair) {
    HacheIter *iter;
    HacheItem *hi;
    pair_data *save_pair;
    string_alloc_t *str_pool;
    int i = 0, i_max;
    pair_queue_t *que;
    
    if (NULL == (save_pair = (pair_data *)malloc(sizeof(pair_data) * pair->write_size))) {
    	fprintf(stderr, "Can't allocate memory in save_pair_data\n");
	return;
    }
    
    str_pool = string_pool_create(1024 * 1024);
    
    if (NULL == str_pool) {
    	fprintf(stderr, "Can't allocate string pool memory in save_pair_data\n");
	return;
    }
    
    iter = HacheTableIterCreate();

    while ((hi = HacheTableIterNext(pair->phache, iter))) {
	save_pair[i].tname = string_dup(str_pool, hi->key);
	save_pair[i].data  = pair_to_pooled_string(str_pool, (pair_loc_t *)hi->data.p);
    	i++;	
    }
    
    i_max = i;
    
    HacheTableIterDestroy(iter);
    
    qsort(save_pair, pair->count, sizeof(pair_data), cmp_pair);
    
    que = add_pair_queue(pair);
    
    if (NULL == que) {
    	fprintf(stderr, "Can't create new pair queue\n");
	return;
    }

    for (i = 0; i < i_max; i++) {
    	fprintf(que->file->fp, "%s %s\n", save_pair[i].tname, save_pair[i].data);
    }

    if (HacheTableEmpty(pair->phache, 1)) {
    	// TEST TEST TEST need to put in proper fail
    	fprintf(stderr, "save_pair_data failed on HacheTableEmpty\n");
    }
    
    string_pool_destroy(str_pool);
    free(save_pair);
    fflush(que->file->fp);
}
    	

static void find_pair(GapIO *io, tg_pair_t *pair, tg_rec recno, char *tname,
	       bin_index_t *bin, contig_t *c, seq_t *seq, tg_args *a,
	       range_t *r_out, library_t *lib) {		
    int new = 0;
    HacheData hd;
    pair_loc_t *pl  = NULL;
    HacheItem *hi   = NULL;
    
    /* Add data for this end */
    pl = (pair_loc_t *)malloc(sizeof(*pl));
    pl->rec    = recno;
    pl->bin    = bin->rec;
    pl->crec   = c->rec;
    pl->pos    = seq->len >= 0 ? seq->pos : seq->pos - seq->len - 1;
    pl->len    = ABS(seq->len);
    pl->idx    = seq->bin_index;
    pl->orient = seq->len < 0;
    pl->flags  = seq->flags;
    pl->mq     = seq->mapping_qual;
    hd.p = pl;
    
    hi = HacheTableAdd(pair->phache, tname, strlen(tname), hd, &new);
        
    if (new) pair->count++;
    
    /* Pair existed already */
    if (!new) {
	pair_loc_t *po = (pair_loc_t *)hi->data.p;
	int st, en;
	
	/* We found one so update r_out now, before flush */
	st = po->pos;
	en = po->pos + (po->orient ? - (po->len-1) : po->len-1);
	r_out->flags &= ~GRANGE_FLAG_TYPE_MASK;
	r_out->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	r_out->pair_rec = po->rec;
	r_out->pair_start = MIN(st,en);
	r_out->pair_end   = MAX(st,en);
	r_out->pair_mqual = po->mq;
	r_out->pair_contig = po->crec;
	r_out->pair_timestamp = io->db->timestamp;
	if ((po->flags & SEQ_END_MASK) == SEQ_END_REV)
	    r_out->flags |= GRANGE_FLAG_PEND_REV;
	if (po->flags & SEQ_COMPLEMENTED)
	    r_out->flags |= GRANGE_FLAG_COMP2;
	
	if (!a->fast_mode) {
	    int st = pl->pos;
	    int en = pl->pos + (pl->orient ? - (pl->len-1) : pl->len-1);
	    bin_index_t *bo;
	    range_t *ro;

	    // Make backwards link only if it's still due to be written out.
	    bo = (bin_index_t *)cache_search_no_load(io, GT_Bin, po->bin);
	    if (bo && cache_lock_mode(io, bo) == G_LOCK_RW) {
		//bo = cache_rw(io, bo);
		bo->flags |= BIN_RANGE_UPDATED;
		ro = arrp(range_t, bo->rng, po->idx);
		ro->flags &= ~GRANGE_FLAG_TYPE_MASK;
		ro->flags |=  GRANGE_FLAG_TYPE_PAIRED;
		ro->pair_rec = pl->rec;
	    } else {
		fprintf(pair->finish->fp,
			"%"PRIrec" %d %"PRIrec" %d %d %d %d %"PRIrec"\n",
			po->bin, po->idx, pl->rec, pl->flags,
			MIN(st, en), MAX(st, en),
			pl->mq, pl->crec);
	    }
	
	    if (po->bin > pair->max_bin) pair->max_bin = po->bin;
	}
	
	if (lib) {
	    /* Increment insert size in library */
	    if (po->crec == pl->crec) {
		int isize = pl->pos - po->pos;
		int ltype;

		/*
		 * We know that 'seq' is the right-most sequence so
		 * this position minus previous position gives us the
		 * insert size and comparing this vs previous orient
		 * we can work out the type of the library. Note although
		 * this read is further right than the previous one, when
		 * overlapping it's possible for the 5' of this minus the
		 * 5' of previous to appear as a negative size.
		 *
		 * Types of pair orientations:
		 *
		 *                           isize   ltype
		 * |------->     <-------|   +ve     IN
		 *
		 * <-------|                 -ve     IN
		 *    |------->
		 *
		 * <-------|     |------->   +ve     OUT
		 *
		 * |------->                 +ve     OUT
		 *    <-------|
		 *
		 * <-------|     <-------|   +ve     SAME
		 *
		 * |------->     |------->   +ve     SAME
		 *
		 * <-------|                 +ve     SAME
		 *    <-------|
		 *
		 * |------->                 +ve     SAME
		 *    |------->
		 */
		if (pl->orient == po->orient) {
		    ltype = LIB_T_SAME;
		} else {
		    if ((isize >= 0 && pl->orient == 1 /* <----| */) ||
			(isize <  0 && pl->orient == 0 /* |----> */))
			ltype = LIB_T_INWARD;
		    else
			ltype = LIB_T_OUTWARD;
		}
		
		lib = cache_rw(io, lib);

		//if (pl->mq > 10 && po->mq > 10) /* good only */
		accumulate_library(io, lib, ltype, ABS(isize));
	    }
	}	    

	/* And, making an assumption, remove from hache */
	HacheTableDel(pair->phache, hi, 1);
	pair->count--;
	free(pl);
    }
    
    if (pair->write_size && pair->count >= pair->write_size) {
    	// too many open pairs, store some away till later
    	fprintf(stderr, "Stored pairs %d\n", pair->count);    
	save_pair_data(pair);
	pair->count = 0;
    }
}


void create_new_contig(GapIO *io, contig_t **c, char *cname, int merge) {

    if (*c) {
	/* Check for inconsistent tags, overlapping the contig end */
	contig_visible_start(io, (*c)->rec, CITER_CSTART);
	contig_visible_end  (io, (*c)->rec, CITER_CEND);
	cache_decr(io, *c);
    }	    

    if (merge) {
	if (NULL == (*c = find_contig_by_name(io, cname)))  {
	    *c = contig_new(io, cname);
	}
    } else {
	char cname2[1024];
	int cname_count = 0;

	/* Not merging, so avoid duplicating contig name */
	snprintf(cname2, sizeof(cname2), "%.*s",
		 (int) sizeof(cname2) - 16, cname);
	while ((*c = find_contig_by_name(io, cname2))) {
	    snprintf(cname2, sizeof(cname2), "%.*s:%d",
		     (int) sizeof(cname2) - 16, cname, ++cname_count);
	}

	if (strcmp(cname, cname2)) {
	    vmessage("Contig name '%s' already existed; renaming to '%s'\n",
		     cname, cname2);
	}

	*c = contig_new(io, cname2);
    }
    
    cache_incr(io, *c);
}

static void sort_file (bttmp_t *old_files[], int div) {
    bttmp_t *new_files[10];
    char line[100];
    int i;

    memset(new_files, 0, sizeof(bttmp_t *) * 10);
    
    for (i = 0; i < 10; i++) {
	new_files[i] = bttmp_file_open(); 
    }
    
    for (i = 0; i < 10; i++) {
    	if (old_files[i]) {
    	    rewind(old_files[i]->fp);

	    while (fgets(line, 100, old_files[i]->fp)) {
    		int bin;
		int mod;

		sscanf(line, "%d", &bin);

		bin = bin / div;

		if (bin) {
	    	    mod = bin % 10;
		} else {
	    	    mod = 0;
		}

		fputs(line, new_files[mod]->fp);
	    }
	    
	    bttmp_file_close(old_files[i]);
	}
	
	old_files[i] = new_files[i];
    }
}


/*
 * Turns a comma separated list of data types into a bit-field.
 */
int parse_data_type(char *type) {
    char *cp;
    int data_type = 0;

    do {
	cp = strchr(type, ',');

	if (0 == strncmp(type, "seq", 3))
	    data_type |= DATA_SEQ;
	else if (0 == strncmp(type, "qual", 4))
	    data_type |= DATA_QUAL;
	else if (0 == strncmp(type, "name", 4))
	    data_type |= DATA_NAME;
	else if (0 == strncmp(type, "anno", 4))
	    data_type |= DATA_ANNO;
	else if (0 == strncmp(type, "all",  3))
	    data_type = DATA_ALL;
	else if (0 == strncmp(type, "none", 4))
	    data_type = 0;
	else if (0 == strncmp(type, "blank", 4))
	    data_type = DATA_BLANK;
	else
	    fprintf(stderr, "Ignoring unknown data_type '%.*s'\n",
		    (int)(cp ? cp-type : strlen(type)), type);

	type = cp ? cp+1 : NULL;
    } while (type);

    return data_type;
}

/* ------------------------------------------------------------------------ */
/* Auto file type detection */
int tg_index_file_type (char *fn) {
    char *suffix = strrchr(fn, '.');
    char data[11];
    zfp *fp;

    /* By standard suffix */
    if (suffix) {
	if (0 == strcmp(suffix, ".gz")) {
	    char *suffix2, tmp;
	    tmp = *suffix;
	    *suffix = 0;
	    suffix2 = strrchr(fn, '.');
	    *suffix = tmp;
	    if (suffix2)
		suffix = suffix2;
	}

	suffix++;
	if (0 == strcmp(suffix, "bam") ||
	    0 == strcmp(suffix, "BAM"))
	    return 'b';

	if (0 == strcmp(suffix, "sam") ||
	    0 == strcmp(suffix, "sam.gz") ||
	    0 == strcmp(suffix, "SAM"))
	    return 's';
	
	if (0 == strcmp(suffix, "cram") ||
	    0 == strcmp(suffix, "CRAM"))
	    return 'e';

	if (0 == strcmp(suffix, "ace") ||
	    0 == strcmp(suffix, "ace.gz") ||
	    0 == strcmp(suffix, "ACE"))
	    return 'A';
	
	if (0 == strcmp(suffix, "baf") ||
	    0 == strcmp(suffix, "baf.gz") ||
	    0 == strcmp(suffix, "BAF"))
	    return 'B';
	
	if (0 == strcmp(suffix, "map") ||
	    0 == strcmp(suffix, "MAP") ||
	    0 == strcmp(suffix, "maq"))
	    return 'm';

	if (0 == strcmp(suffix, "caf") ||
	    0 == strcmp(suffix, "CAF"))
	    return 'C';

	if (0 == strcmp(suffix, "fna") ||
	    0 == strcmp(suffix, "FNA") || 
	    0 == strcmp(suffix, "fasta") || 
	    0 == strcmp(suffix, "FASTA")) {
	    return 'F';
	}

	if (0 == strcmp(suffix, "fnq") ||
	    0 == strcmp(suffix, "FNQ") || 
	    0 == strcmp(suffix, "fastq") || 
	    0 == strcmp(suffix, "FASTQ")) {
	    return 'Q';
	}
    }

    /* By contents */
    if (NULL == (fp = zfopen(fn, "rb"))) {
	perror(fn);
	return '?';
    }

    if (NULL == zfgets(data, 10, fp)) {
	zfclose(fp);
	return '?';
    }
    zfclose(fp);

    if (0 == strncmp(data, "BAM\001", 4))
	return 'b'; /* bam */

    if (0 == strncmp(data, "CRAM", 4))
	return 'e'; /* cram */

    if (0 == strncmp(data, "AS ", 3))
	return 'A'; /* ace */

    /* Gets trickier to detect from here on */
    if (0 == strncmp(data, "CO=", 3))
	return 'B'; /* baf */

    if (0 == strncmp(data, "@HD\t", 3) ||
	0 == strncmp(data, "@SQ\t", 3) ||
	0 == strncmp(data, "@RG\t", 3) ||
	0 == strncmp(data, "@PG\t", 3))
	return 's'; /* sam */


    /* These are now pretty tenuous */
    if (*data == '>') {
	return 'F'; /* fasta */
    }

    if (*data == '@') {
	return 'Q'; /* fastq */
    }

    /*
     * And if still not found, well it maybe maq, but is just as likely
     * a differently formatted sam or baf. Give up at this point.
     */

    return '?';
}


/* ------------------------------------------------------------------------ */

/*
 * Relaces \n with newline and \\ with \.
 * Modifies the line in-situ as it can never grow.
 */
void unescape_line(char *txt) {
    char *cp;
    for (cp = txt; *txt; txt++) {
	if (*txt != '\\') {
	    *cp++ = *txt;
	} else {
	    if (*++txt == 'n')
		*cp++ = '\n';
	    else
		*cp++ = *txt;
	    if (!*txt)
		break;
	}
    }
    *cp++ = 0;
}


tg_pair_t *create_pair(int queue) {
    tg_pair_t *pair;
    
    if (NULL == (pair = (tg_pair_t *)malloc(sizeof(tg_pair_t)))) {
    	return NULL;
    }
    
    pair->que          = NULL;
    pair->que_size     = 0;
    pair->working_size = 1000;  
    pair->write_size   = queue; 
    pair->count        = 0;
    pair->phache       = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
    pair->phache->name = "pair";
    
    if (NULL == (pair->finish = bttmp_file_open())) {
    	free(pair);
	return NULL;
    }
    
    pair->max_bin = 0;
    
    return pair;
}

    
tg_rec save_range_sequence(GapIO *io, seq_t *seq, uint8_t mapping_qual,
			   tg_pair_t *pair, int is_pair, char *tname,
			   contig_t *c, tg_args *a, int flags, library_t *lib,
			   tg_rec *bin_rec) {
    range_t r, *r_out;
    tg_rec recno;
    bin_index_t *bin;
    static tg_rec fake_recno = 1;
    int comp;

    /* Update sequence library, aka read-group */
    if (lib && !seq->parent_type) {
	seq->parent_type = GT_Library;
	seq->parent_rec = lib->rec;
    }

    /* Create range */
    r.start = seq->pos;
    r.end   = seq->pos + ABS(seq->len)-1;
    r.rec   = 0;
    r.mqual = mapping_qual;
    r.pair_rec = 0;
    r.flags = flags;
    r.library_rec = lib ? lib->rec : 0;
    r.y     = 0;

    r.pair_contig    = 0;
    r.pair_timestamp = 0;
    r.pair_start     = 0;
    r.pair_end       = 0;
    r.pair_mqual     = 0;

    /* Add the range to a bin, and see which bin it was */
    bin = bin_add_range(io, &c, &r, &r_out, &comp, 1);
    if (bin_rec)
	*bin_rec = bin->rec;

    /* Save sequence */
    if (a->data_type == DATA_BLANK) {
	recno = fake_recno++;
    } else {
	if (comp) {
	    complement_seq_t(seq);
	    seq->len = -seq->len;
	}

	recno = save_sequence(io, seq, bin, r_out);
    }

    if (is_pair) {
	find_pair(io, pair, recno, tname, bin, c, seq, a, r_out, lib);
    }

    if (a->tmp && (a->data_type & DATA_NAME))
	bttmp_file_store(a->tmp, seq->name_len, seq->name, recno);
    
    if (seq->name)
	free(seq->name);

    /* Link bin back to sequence too before it gets flushed */
    r_out->rec = recno;
    
    return recno;
}
    
    

   
static int sort_pair_file(tg_pair_t *pair) {
    bttmp_t  *old_files[11];
    int div = 1;
    int max_div = 10; /* temp, needs to be variable */
    int i = 0;
    bttmp_t *final;
    tg_rec max_bin = pair->max_bin; 
    
    memset(old_files, 0, sizeof(bttmp_t *) * 11);

    old_files[0] = pair->finish;
    
    
    while ((max_bin % max_div) != max_bin) {
    	max_div *= 10;
    }
    
    while (div < max_div) {
    	sort_file(old_files, div);
	div = div * 10;
    }
    
    /* gather files together here */
    
    final = bttmp_file_open();
    
    while (old_files[i]) {
    	char line[100];
     	rewind(old_files[i]->fp);
	
	while (fgets(line, 100, old_files[i]->fp)) {
  	    fputs(line, final->fp);
	}
 	
     	bttmp_file_close(old_files[i++]);
    }
    
    pair->finish = final;
    
    return 1;
}

/*
 * If we have singletons still in our pair struct and we're using append mode
 * then maybe they are pairs of data that was already in the database.
 * Check for this case.
 */
static void merge_pairs(GapIO *io, tg_pair_t *pair) {
    HacheIter *iter;
    HacheItem *hi;
    
    iter = HacheTableIterCreate();

    while ((hi = HacheTableIterNext(pair->phache, iter))) {
	/* FIXME: sort these */
	tg_rec *rec;
	int i, nr;
	char name[8192];
	seq_t *s2;
	bin_index_t *b2 = NULL;
	range_t *r = NULL;
	pair_loc_t *p;
	int st, en;
	tg_rec bin_rec, contig;
	int bin_idx;

	memcpy(name, hi->key, hi->key_len);
	name[hi->key_len] = 0;

	/* Hunt for all occurences of this name and see if one is singleton */
	rec = sequence_index_query_all(io, name, 0, &nr);
	if (!rec)
	    continue;

	//printf("%.*s can be paired with", hi->key_len, hi->key);
	for (i = 0; i < nr; i++) {
	    //printf(" #%"PRIrec, rec[i]);
	    if (!(s2 = cache_search(io, GT_Seq, rec[i])))
		continue;
	    cache_incr(io, s2);
	    if (!(b2 = cache_search(io, GT_Bin, s2->bin))) {
		cache_decr(io, s2);
		continue;
	    }

	    bin_rec = s2->bin;
	    bin_idx = s2->bin_index;

	    r = arrp(range_t, b2->rng, s2->bin_index);
	    cache_decr(io, s2);

	    assert(r->rec == s2->rec);
	    if (r->pair_rec == 0)
		break;
	}
	//printf("\n");

	free(rec);

	if (!r || r->pair_rec)
	    continue;

	/* We found a singleton, so link it up */
	p = (pair_loc_t *)hi->data.p;
	st = p->pos;
	en = p->pos + (p->orient ? -(p->len-1) : p->len-1);

	/* Link other end to us; done at end with others */
	fprintf(pair->finish->fp, "%"PRIrec" %d %"PRIrec" %d %d %d %d %"PRIrec"\n",
		bin_rec, bin_idx, p->rec, p->flags,
		MIN(st, en), MAX(st, en), p->mq, p->crec);

	/* Us to other end is harder; we don't know the loc */
	cache_incr(io, b2);
	bin_get_item_position(io, GT_Seq, r->rec, &contig, &st, &en,
			      NULL, NULL, NULL, NULL);

	fprintf(pair->finish->fp, "%"PRIrec" %d %"PRIrec" %d %d %d %d %"PRIrec"\n",
		p->bin, p->idx, r->rec, r->flags, st, en, r->mqual, contig);

	cache_decr(io, b2);
    }

    HacheTableIterDestroy(iter);
    fflush(pair->finish->fp);
}


static void complete_pairs(GapIO *io, tg_pair_t *pair) {
    bin_index_t *bo = NULL;
    range_t *ro;
    tg_rec current_bin = -1;
    char line[1024];
    int rec_count = 0;
    int total_count = 0;
    
    rewind(pair->finish->fp);
    
    while (fgets(line, 1024, pair->finish->fp)) {
	int idx, flags;
	tg_rec bin, rec, pair_contig;
	int pair_start, pair_end, pair_mqual;
	
        sscanf(line, "%"PRIrec" %d %"PRIrec" %d %d %d %d %"PRIrec,
	       &bin, &idx, &rec, &flags, &pair_start, &pair_end,
	       &pair_mqual, &pair_contig);

	//printf("%s", line);

	if (bin != current_bin) {
	    if (rec_count > 50000) {
	    	total_count += rec_count;
	    	cache_flush(io);
		
		fprintf(stderr, "%d pairs finished so far\n", total_count);
	    
		rec_count = 0;
	    }

	    bo = (bin_index_t *)cache_search(io, GT_Bin, bin);
	    bo = cache_rw(io, bo);
	    bo->flags |= BIN_RANGE_UPDATED;
    	    current_bin = bin;
	}

	ro = arrp(range_t, bo->rng, idx);
	ro->flags &= ~GRANGE_FLAG_TYPE_MASK;
	ro->flags |=  GRANGE_FLAG_TYPE_PAIRED;
	ro->pair_rec = rec;
	ro->pair_contig = pair_contig;
	ro->pair_start  = pair_start;
	ro->pair_end    = pair_end;
	ro->pair_mqual  = pair_mqual;
	ro->pair_timestamp = io->db->timestamp;
	if ((flags & SEQ_END_MASK) == SEQ_END_REV)
	    ro->flags |= GRANGE_FLAG_PEND_REV;
	if (flags & SEQ_COMPLEMENTED)
	    ro->flags |= GRANGE_FLAG_COMP2;
	
	rec_count++;

    }
    
    total_count += rec_count;
    
    fprintf(stderr, "%d pairs finished in total.\n", total_count);
    
    cache_flush(io);
}


/* an alternative to the Gnu specific getline */
long tg_get_line(char **line, long *length, FILE *fp) {
    char *in_line;
    long len;
    long offset = 0;
    long in_chars;
    
    if (line == NULL || fp == NULL || length == NULL) {
    	return -1;
    }
    
    if (*line == NULL || *length <= 0) {
    	if ((*line = malloc(256 * sizeof(char))) == NULL) {
	    return -1;
	}
	
	*length = 256;
    }
    
    in_line = *line;
    len     = *length;
    
    while ((fgets((in_line + offset), (len - offset), fp))) {
    	char *tmp = NULL;
	
    	in_chars = strlen(in_line);
	
    	offset = in_chars;
	
	// see if we have our full line
	if (*(in_line + offset - 1) == '\n') {
	    break;
	}
    	
    	len += len;
	tmp = realloc(in_line, len * sizeof(char));
	    
    	if (tmp) {
	    in_line = tmp;
	} else {
	    fprintf(stderr, "Memory error in get_line\n");
	    return -1;
	}
    }
    
    *line = in_line;
    *length = len;
    
    return offset;
}


static int load_data(pair_queue_t *pq) {
    int i;
    char *line_in = NULL;
    long line_size = 0;
    
    if (pq->name_pool) {
    	string_pool_destroy(pq->name_pool);
    }
    
    pq->name_pool = string_pool_create(1024);
    
    for (i = 0; i < pq->pair_size; i++) {
    	char name[1024];
    	int found;
    	int ret = tg_get_line(&line_in, &line_size, pq->file->fp);
	
	if (ret <= 0) break;
	
	found = sscanf(line_in, "%s %"PRId64" %"PRId64" %d %"PRId64
		       " %d %d %d %d %d",
		       name, &pq->pair[i].rec, &pq->pair[i].bin,
		       &pq->pair[i].idx, &pq->pair[i].crec, &pq->pair[i].pos,
		       &pq->pair[i].orient, &pq->pair[i].flags,
		       &pq->pair[i].len, &pq->pair[i].mq);
	    
	if (found != 10) {
	    fprintf(stderr, "Error found in line: %s\n", line_in);
	    break;
	}
	
	pq->pair[i].tname = string_dup(pq->name_pool, name);
    }
    
    pq->pair_size = i;
    pq->index = 0;
    
    if (line_in) free(line_in);
    
    return pq->pair_size;
}


static int initialise_queues(tg_pair_t *pair) {
    int i;
    
    for (i = 0; i < pair->que_size; i++) {
    	rewind(pair->que[i].file->fp);
	
	if (NULL == (pair->que[i].pair = (pair_loc_t *)malloc(sizeof(pair_loc_t) * pair->working_size))) {
	   fprintf(stderr, "Out of memory allocating pairs in initialise_queues\n");
	   return -1;
	}
	
	pair->que[i].name_pool = NULL;
	pair->que[i].index = 0;
	pair->que[i].pair_size = pair->working_size;
	
	if (0 == load_data(&pair->que[i])) {
	    fprintf(stderr, "Initial data load failed on file %d\n", i);
	    return -1;
	}
    }
    
    return 0;
}


static void get_next(pair_queue_t *que) {
    que->index++;
    
    if (que->index == que->pair_size) {
    	load_data(que);
    }
}


static void save_match_pair(tg_pair_t *pair, pair_loc_t *p1, pair_loc_t *p2) {
    int st, en;

    st = p2->pos;
    en = p2->pos + (p2->orient ? - (p2->len-1) : p2->len-1);
    fprintf(pair->finish->fp, "%"PRIrec" %d %"PRIrec" %d %d %d %d %"PRIrec"\n",
	    p1->bin, p1->idx, p2->rec, p2->flags,
	    MIN(st, en), MAX(st, en), p2->mq, p2->crec);
	
    st = p1->pos;
    en = p1->pos + (p1->orient ? - (p1->len-1) : p1->len-1);
    fprintf(pair->finish->fp, "%"PRIrec" %d %"PRIrec" %d %d %d %d %"PRIrec"\n",
	    p2->bin, p2->idx, p1->rec, p1->flags,
	    MIN(st, en), MAX(st, en), p1->mq, p1->crec);
}
	


static int find_saved_pairs(GapIO *io, tg_pair_t *pair) {
    int num_found = 0;
    int i;
    int still_looking = 1;

    // scan the files until there is nothing left to see
    
    while (still_looking) {
    	char *compare = NULL;
	int file = 0;
	int done = 0, match = 0;
	
	for (i = 0; i < pair->que_size; i++) {
	    pair_queue_t *que = &pair->que[i];
	
	    if (que->pair_size) {
	    	done++;
	    
	    	if (compare == NULL) {
		    compare = que->pair[que->index].tname;
		    file = i;
		} else {
		    int cmp = strcmp(compare, que->pair[que->index].tname);
		    
		    if (cmp > 0) {
		    	compare = que->pair[que->index].tname;
		    	file = i;
		    } else if (cmp == 0) {
		    	match = i;
			break;
		    }
		}
	    }
	}
	
	if (done) {
	    if (match) {
		save_match_pair(pair,
				&pair->que[match].pair[pair->que[match].index],
				&pair->que[file].pair[pair->que[file].index]);
		get_next(&pair->que[match]);
		num_found++;
	    }

	    get_next(&pair->que[file]);
	} else {
	    still_looking = 0;
	}
    }
    
    return num_found;
}


void finish_pairs(GapIO *io, tg_pair_t *pair, int link_pairs) {
    int ret = 1;
    
    if (pair->que_size) {
    	// deal with left over pairs
    	save_pair_data(pair);

    	fprintf(stderr, "Resolving pair queues. Total is %d\n", pair->que_size);
    
	initialise_queues(pair);
	ret = find_saved_pairs(io, pair);
    
    	fprintf(stderr, "%d pairs found\n", ret);
    } else {
    	fprintf(stderr, "No pair queue found\n");
    }
    
    if (link_pairs)
	merge_pairs(io, pair);
    ret = sort_pair_file(pair);
    
    if (ret) {
    	fprintf(stderr, "Run complete pairs.\n");
    	complete_pairs(io, pair);
    } else {
    	fprintf(stderr, "sort_pair_file failed");
    }
    
    fprintf(stderr, "Pairs complete\n");
}


void delete_pair(tg_pair_t *pair) {
    int i;
    
    for (i = 0; i < pair->que_size; i++) {
    	if (pair->que[i].file) bttmp_file_close(pair->que[i].file);
	
	if (pair->que[i].pair) free(pair->que[i].pair);
	
	if (pair->que[i].name_pool) string_pool_destroy(pair->que[i].name_pool);
    }
    
    if (pair->que) free(pair->que);
    
    if (pair->phache) HacheTableDestroy(pair->phache, 1);
    
    if (pair->finish) bttmp_file_close(pair->finish);
    
    free(pair);
}
    


