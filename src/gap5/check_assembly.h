#ifndef _CHECK_ASSEMBLY_H_
#define _CHECK_ASSEMBLY_H_

/*
 * Scans through one or more contigs checking each reading for correct
 * assembly.
 *
 * Returns -1 for failure, 0 for success.
 */
int check_assembly(GapIO *io, int num_contigs, contig_list_t *contigs, 
		   int winsize, float maxperc, int ignore_N);


void *checkass_obj_func(int job, void *jdata, obj_checkass *obj,
			mobj_checkass *ca);

void check_assembly_callback(GapIO *io, tg_rec contig, void *fdata,
			     reg_data *jdata);
#endif
