#include "tg_gio.h"

int find_oligos(GapIO *io,
		int num_contigs,
		contig_list_t *contig_array,
		float mis_match,
		char *string,
		int consensus_only,
		int in_cutoff);

int
find_oligo_file(GapIO *io,
		int num_contigs,
		contig_list_t *contig_array,
		float mis_match,
		char *file,
		int consensus_only,
		int in_cutoff);

void find_oligo_callback(GapIO *io, tg_rec contig, void *fdata, reg_data *jdata);

void *find_oligo_obj_func1(int job,
			   void *jdata,
			   obj_match *obj,
			   mobj_find_oligo *find_oligo);
void *find_oligo_obj_func2(int job,
			   void *jdata,
			   obj_match *obj,
			   mobj_find_oligo *find_oligo);

