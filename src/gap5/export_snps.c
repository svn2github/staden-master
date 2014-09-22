/*
 * NOTE: This is experimental and is here simply for purposes of evaluating
 * how well Gap5's consensus algorithm and shuffle pads is working. It has hard
 * coded values and contig names (eg "chr20") and is not suitable yet as a
 * general purpose VCF generator. It also has bugs in the output when dealing
 * with compound SNPs/Indels.
 */

/* ------------------------------------------------------------------------
 * export_snps implementation.
 *
 * This computes the consensus in heterozygous mode and then generates a
 * VCF file of the heterozygous differences.  Homozygous differences require
 * an external reference to call, which isn't yet utilised.
 */

#include <staden_config.h>

#include <tcl.h>
#include <assert.h>
#include <ctype.h>

#include <tg_gio.h>
#include "misc.h"
#include "export_snps.h"
#include "gap_cli_arg.h"
#include "list_proc.h"
#include "xalloc.h"
#include "consensus.h"
#include "dstring.h"
#include <io_lib/cram.h>

static int export_snps(GapIO *io, int cc, contig_list_t *cv,
		       int depad, char *fn, char *fai);

typedef struct {
    GapIO *io;
    char *inlist;
    char *outfile;
    char *fai;
    int   unpadded;
} es_arg;

int tcl_export_snps(ClientData clientData, Tcl_Interp *interp,
		    int objc, Tcl_Obj *CONST objv[])
{
    int rargc, format_code, res;
    contig_list_t *rargv;
    es_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,     offsetof(es_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,     offsetof(es_arg, inlist)},
	{"-outfile",    ARG_STR, 1, "out.vcf",offsetof(es_arg, outfile)},
	{"-fai",        ARG_STR, 1, "",       offsetof(es_arg, fai)},
	{"-unpadded",   ARG_INT, 1, "1",      offsetof(es_arg, unpadded)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs_extended(args.io, args.inlist, &rargc, &rargv);

    res = export_snps(args.io, rargc, rargv, args.unpadded, args.outfile, args.fai);

    free(rargv);

    return res == 0 ? TCL_OK : -1;
}

static int export_snps_contig(GapIO *io, tg_rec contig, int start, int end,
			      int depad, refs_t *ref, FILE *fp) {
    int *ref_pos = NULL, *ref_id = NULL;
    char *cons = NULL;
    float *qual = NULL, last_qual = 30;
    int ret = -1, i, start_pos, last_pos, last_id = -1;
    char *ref_seq = NULL;

    if (!(ref_pos = malloc((end-start+1) * sizeof(int))))
	goto error;

    if (!(ref_id = malloc((end-start+1) * sizeof(int))))
	goto error;

    if (!(cons = malloc(end-start+1)))
	goto error;

    if (!(qual = malloc((end-start+1)*sizeof(*qual))))
	goto error;

    if (calculate_consensus_simple_het(io, contig, start, end, cons, qual))
	goto error;

    if (padded_to_reference_array(io, contig, start, end, ref_pos,
				  ref_id, &start_pos, NULL))
	goto error;
    last_pos = start_pos;

    for (i = 0; i <= end-start; i++) {
	float score;
	int j;

	if (ref_pos[i] == INT_MIN) {
	    char ins[1024];
	    int upper = isupper(cons[i]);
	    int len = 0;

	    char ref_base = ref_seq
		? ref_seq[last_pos - start_pos]
		: (i ? toupper(cons[i-1]) : 'N');
	    char con_base = i ? toupper(cons[i-1]) : 'N';

	    score = 0;
	    while (ref_pos[i] == INT_MIN && i <= end-start) {
		if (cons[i] != '*' && isupper(cons[i]) == upper) {
		    ins[len++] = toupper(cons[i]);
		    score += qual[i];
		}
		i++;
	    }
	    i--;
	    if (len) {
		if (upper)
		    //fprintf(fp, "%d\tINS %d %.*s\t%.0f\n", last_pos, len,
		    //	    len,ins, score/len);
		    fprintf(fp,
			    "chr20\t%d\t.\t%c\t%c%.*s\t%.0f\tPASS\t.\tGT\t1/1\n",
			    last_pos,
			    ref_base, con_base,
			    len,ins, score/len);
		else
		    //fprintf(fp, "%d\tins %d %.*s\t%.0f\n", last_pos, len,
		    //	    len,ins, score/len);
		    fprintf(fp,
			    "chr20\t%d\t.\t%c\t%c%.*s\t%.0f\tPASS\t.\tGT\t0/1\n",
			    last_pos,
			    ref_base, con_base,
			    len,ins, score/len);
	    }
	    continue;

	} else {
	    //if (ref_id[i] != last_id) {
	    if (last_id != 20) {
		if (ref && ref->fp) {
		    ref_entry *e;
		    //e = ref->ref_id[ref_id[i]];
		    e = ref->ref_id[22]; // HACK FOR TESTING
		    printf("Loading %d..%d\n", start_pos, (int)e->length);
		    ref_seq = load_ref_portion(ref->fp, e, start_pos, e->length);
		}

		//last_id = ref_id[i];
		last_id = 20;
		for (j = i; ref_pos[j] == INT_MIN; j--) ;
		last_pos = ref_pos[j]-1;
	    }

	    if (ref_pos[i] > last_pos+1 || cons[i] == '*') {
		char ref_base = ref_seq
		    ? ref_seq[last_pos - start_pos]
		    : (i ? toupper(cons[i-1]) : 'N');
		char con_base = i ? toupper(cons[i-1]) : 'N';

		while (cons[i] == '*' && i <= end-start)
		    i++;
		//fprintf(fp, "%d\tDEL %d %.*s\t%.0f\n", last_pos,
		//	ref_pos[i]-(last_pos+1),
		//	ref_pos[i]-(last_pos+1), ref_seq + last_pos+1-start_pos,
		//	(last_qual+qual[i])/2);
		fprintf(fp, "chr20\t%d\t.\t%c%.*s\t%c\t%.0f\tPASS\t.\tGT\t1/1\n",
			last_pos, ref_base, 
			ref_pos[i]-(last_pos+1), ref_seq + last_pos+1-start_pos,
			con_base,
			(last_qual+qual[i])/2);
		//if (i && cons[i-1] == '*') i--;
	    }
	}

	switch(cons[i]) {
	    #define _A 'A'
	    #define _C 'C'
	    #define _G 'G'
	    #define _T 'T'
	    static char B1[] = {
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 00
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 10
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 20
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 30
	     //    A  B  C   D  E  F  G   H  I  J  K   L  M  N  O
		0,_A,_C,_C, _A, 0, 0,_G, _A, 0, 0,_G,  0,_A,_A, 0, // 40
	     // P  Q  R  S   T  U  V  W   X  Y  Z
		0, 0,_A,_C, _T,_T,_A,_A,  0,_C, 0, 0,  0, 0, 0, 0, // 50
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 60
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 70
	    };
	    static char B2[] = {
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 00
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 10
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 20
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 30
	     //    A  B  C   D  E  F  G   H  I  J  K   L  M  N  O
		0,_A,_G,_C, _G, 0, 0,_G, _C, 0, 0,_T,  0,_C,_C, 0, // 40
	     // P  Q  R  S   T  U  V  W   X  Y  Z
		0, 0,_G,_G, _T,_T,_C,_T,  0,_T, 0, 0,  0, 0, 0, 0, // 50
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 60
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, // 70
	    };
	    char b1, b2, r;

	case 'A': case 'C': case 'G': case 'T':
	    if (ref_seq && ref_seq[ref_pos[i] - start_pos] != cons[i]) {
		//fprintf(fp, "%d\tHOM %c %c\t%.0f\n", ref_pos[i],
		//	ref_seq[ref_pos[i] - start_pos], cons[i], qual[i]);
		fprintf(fp, "chr20\t%d\t.\t%c\t%c\t%.0f\tPASS\t.\tGT\t1/1\n",
			ref_pos[i],
			ref_seq[ref_pos[i] - start_pos], toupper(cons[i]), qual[i]);
	    }
	    break;

	case 'M':
	case 'R':
	case 'W':
	case 'S':
	case 'Y':
	case 'K':
	    b1 = B1[cons[i]];
	    b2 = B2[cons[i]];
	    r = toupper(ref_seq[ref_pos[i] - start_pos]);
	    //fprintf(fp, "%d\tHET %c ", ref_pos[i], r);
	    fprintf(fp, "chr20\t%d\t.\t%c\t", ref_pos[i], r);
	    if (b1 == r)	
		//fprintf(fp, "%c\t%.0f\n", b2, qual[i]);
		fprintf(fp, "%c\t%.0f\tPASS\t.\tGT\t0/1\n", b2, qual[i]);
	    else if (b2 == r)
		//fprintf(fp, "%c\t%.0f\n", b1, qual[i]);
		fprintf(fp, "%c\t%.0f\tPASS\t.\tGT\t0/1\n", b1, qual[i]);
	    else
		//fprintf(fp, "%c,%c\t%.0f\n", b1,b2, qual[i]);
		fprintf(fp, "%c,%c\t%.0f\tPASS\t.\tGT\t1/2\n", b1,b2, qual[i]);
	    break;

	case 'N':
	    //fprintf(fp, "%d\tunknown-call\n", ref_pos[i]);
	    break;

	case 'a':
	case 'c':
	case 'g':
	case 't': {
	    int len = 0, j;
	    float score = 0;

	    char ref_base = ref_seq
		? ref_seq[last_pos - start_pos]
		: (i ? toupper(cons[i-1]) : 'N');
	    char con_base = i ? toupper(cons[i-1]) : 'N';

	    // heterozygous deletion, see how long
	    // Assuming it's a simple case, track how many [acgt] match the ref
	    // and produce xACGT->x 0/1 call.
	    while (islower(cons[i]) &&
		   ref_pos[i] != INT_MIN &&
		   toupper(cons[i]) == ref_seq[ref_pos[i] - start_pos] &&
		   i <= end-start)
		score+=qual[i], i++;

	    
	    for (j=i; ref_pos[j] == INT_MIN; j--)
		;
	    len = ref_pos[j] - last_pos + i-j;

	    //fprintf(fp, "%d\tdel %.*s %.0f\n", last_pos,
	    //	    len, &ref_seq[last_pos - start_pos],
	    //	    score/len);
	    fprintf(fp, "chr20\t%d\t.\t%.*s\t%c\t%.0f\tPASS\t.\tGT\t0/1\n",
		    last_pos, len, &ref_seq[last_pos - start_pos],
		    con_base, score/len);
	    
	    //if (len) i--;
	    break;
	}
	}

	for (j = i; ref_pos[j] == INT_MIN; j--) ;
	last_pos = ref_pos[j];
	last_qual = qual[j];

	// Homozygous deletion+change?
	// Ref = GAAT. Seqs = GCCT + GT
	// GAAT -> GCCT,GT  1/2

	// What about compound homozygous + heterozygous deletion?
	// Ie one allele losing 3bp and the other losing 6bp.
	//
	// GTCTAAA -> GTCT,G 1/2

	// Also compound deletion + insertion?
	// GTCT -> G     1/1
	// G    -> GACC  1/1
	// Or is it just 3 SNPs in a row?
    }

    ret = 0;
 error:
    if (ref_pos)
	free(ref_pos);

    if (ref_id)
	free(ref_id);

    if (cons)
	free(cons);

    if (qual)
	free(qual);

    return ret;
}

static int export_snps(GapIO *io, int cc, contig_list_t *cv,
		       int depad, char *fn, char *fai) {
    int i, r = 0;
    FILE *fp;
    refs_t *refs = NULL;

    fp = fopen(fn, "w");
    if (!fp)
	return -1;
    
    fprintf(fp, "##fileformat=VCFv4.0\n");
    fprintf(fp, "##source=VCFWriter\n");
    fprintf(fp, "##INFO=<ID=OP,Number=1,Type=Integer,Description=\"Original position before normalization\">\n");
    fprintf(fp, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tgap5\n");

    if (fai && *fai)
	refs = refs_load_fai(NULL, fai, 1);

    for (i = 0; i < cc; i++) {
	r |= export_snps_contig(io, cv[i].contig, cv[i].start, cv[i].end,
				depad, refs, fp);
    }

    fclose(fp);
    if (refs)
	refs_free(refs);

    return r;
}
