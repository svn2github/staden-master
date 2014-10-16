/*
 * ---------------------------------------------------------------------------
 * Identifies SNP sites and forms haplotypes from these by linking SNPs
 * together based on spanning readings. 
 *
 * Trim anything that maybe dodgy (quality differs to background?).
 *
 * Iteratively build up haplotypes. Haplotype is [ACGT*-]+.
 * Ie a sequence or unknown (-).
 * A sequence matches a known haplotype if it disagrees with - only.
 * Eg ---ACGTA--- and ----CGTAC-- align and will be merged.
 *
 * For efficiency leading and trailing "-" are run-length encoded.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "misc.h"
#include "find_haplotypes.h"
#include "consensus.h"
#include "align.h"
#include "array.h"

/*
 * A string of haplotypic bases. This excludes non-SNP sites
 */
typedef struct haplotype_str {
    struct haplotype_str *next;
    char *snps;          // string of [ACGT*-]
    int *count;          // depth of snps[]
    int nsnps;
    int nseq;
    // start by ignoring this and use noddy implementation
    int start;           // number of leading "-"
    int end;             // (end-start) = len
    Array recs;
} haplotype_str;

/*
 * Compares 'snps' against the known haplotypes hs[0..n_hs-1].
 *
 * If it matches it incorporates it, otherwise it adds a new
 * haplotype.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int haplotype_str_add(haplotype_str **hs, char *snps, int nsnp, tg_rec rec) {
    haplotype_str *tmp, *last = NULL, *last_last = NULL;
    int overlap = 0;
    int i;

    // FIXME: need something better than linked list. Maybe a way to skip 
    // ones that are out of range (ended already).
    // Also have start/end to avoid making them so long.

    for (tmp = *hs; tmp; last_last = last, last = tmp, tmp = tmp->next) {
	for (i = 0; i < nsnp; i++) {
	    if (tmp->snps[i] != '-' && snps[i] != '-') {
		overlap=1;
		if (tmp->snps[i] != snps[i]) {
		    overlap = 0;
		    break;
		}
	}
	}
	if (i != nsnp)
	    continue;

	// Otherwise it matches, so update it
	if (overlap) {
	    for (i = 0; i < nsnp; i++) {
		if (snps[i] != '-') {
		    tmp->snps[i] = snps[i];
		    tmp->count[i]++;
		}
	    }
	    tmp->nseq++;

	    // Maintain sorted list by nseq
	    if (last && tmp->nseq > last->nseq) {
		if (last_last) {
		    last_last->next = tmp;
		    last->next = tmp->next;
		    tmp->next = last;
		} else {
		    *hs = tmp;
		    last->next = tmp->next;
		    tmp->next = last;
		}
	    }
	    ArrayPush(tmp->recs, tg_rec, rec);
	    return 0;
	}
    }


    // Hasn't been merged, so start a new haplotype string
    if (last)
	tmp = last->next = calloc(1, sizeof(*last->next));
    else 
	tmp = *hs = calloc(1, sizeof(*last->next));

    if (!tmp)
	return -1;

    tmp->snps = (char *)malloc(nsnp);
    tmp->count = (int *)calloc(nsnp, sizeof(int));
    tmp->nsnps = nsnp;
    tmp->nseq = 1;
    for (i = 0; i < nsnp; i++) {
	if ((tmp->snps[i] = snps[i]) != '-')
	    tmp->count[i] = 1;
    }

    tmp->recs = ArrayCreate(sizeof(tg_rec), 1);
    ArrayPush(tmp->recs, tg_rec, rec);

    return 0;
}

void haplotype_str_dump(haplotype_str *hs) {
    while (hs) {
	int i, d;
	for (i = d = 0; i < hs->nsnps; i++)
	    d += hs->count[i];
	printf("%5d %.*s\n", hs->nseq, hs->nsnps, hs->snps);
//	printf("%5d ", hs->nseq);
//	for (i = 0; i < hs->nsnps; i++)
//	    putchar('!'+MIN(90,hs->count[i]));
//	putchar('\n');

//	for (i = 0; i < ArrayMax(hs->recs); i++) {
//	    printf("\t#%"PRIrec"\n", arr(tg_rec, hs->recs, i));
//	}
	hs = hs->next;
    }
}


/*
 * A simple linked list of haplotypic sites, forming a doubly linked list so we can easily remove from it.
 */
typedef struct haplotype_pos {
    int pos;             // pos
    int score;           // FIXME: define this.
    struct haplotype_pos *prev, *next;
} haplotype_pos;

int add_haplotype_pos(haplotype_pos **phead, haplotype_pos **ptail, int pos) {
    haplotype_pos *p = calloc(1, sizeof(*p));
    if (!p)
	return -1;

    p->pos = pos;

    if (*ptail) {
	(*ptail)->next = p;
	p->prev = *ptail;
	*ptail = p;
    } else {
	*phead = *ptail = p;
    }

    return 0;
}

void del_haplotype_pos(haplotype_pos **phead, haplotype_pos **ptail,
		       haplotype_pos *p) {
    if (p == *phead)
	*phead = p->next;
    else
	p->prev->next = p->next;

    if (p == *ptail)
	*ptail = p->prev;
    else
	p->next->prev = p->prev;

    free(p);
}


static int find_haplotypes_single(GapIO *io, tg_rec crec, int start, int end,
				  Array rec_list) {
    consensus_t *cons = NULL;
    int ret = -1, i;
    haplotype_pos *phead = NULL, *ptail = NULL;
    int pass;
    rangec_t *rng = NULL;
    int nr;
    contig_t *c;
    int nsnps = 0;
    haplotype_str *hs = NULL;
    char *hstr = NULL;

    // Accumulate a list of haplotypes
    if (!(cons = calloc(end-start+1, sizeof(*cons))))
	goto err;

    if (-1 == calculate_consensus(io, crec, start, end, cons))
	goto err;

    for (i = start; i <= end; i++) {
	if (cons[i-start].scores[6]>10) { // FIXME: Or > some specified minimum haplotype score.
	    printf("Pos %5d: het %c/%c  score %d\n",
		   i,
		   "ACGT*"[cons[i-start].het_call / 5],
		   "ACGT*"[cons[i-start].het_call % 5],
		   (int)cons[i-start].scores[6]);

	    add_haplotype_pos(&phead, &ptail, i);
	    nsnps++;
	}
    }

    hstr = malloc(nsnps);

    c = cache_search(io, GT_Contig, crec);
    if (!c)
	goto err;

    rng = contig_seqs_in_range(io, &c, start, end,
			       CSIR_SORT_BY_X | CSIR_SORT_BY_CLIPPED, &nr);
    if (!rng)
	goto err;

    // Accumulate haplotypes
    {
	rangec_t *r;
	haplotype_pos *p1, *p2;
	int i;
	int snp_no = 0;

	p1 = phead;
	for (i = 0; i < nr; i++) {
	    rangec_t *r = &rng[i];
	    int left, right;
	    seq_t *s;
	    char b;
	    int snp_no2;

	    // FIXME: optimise, no need to reset all the while.
	    // Fill out bits we need and remember start/end.
	    memset(hstr, '-', nsnps); 

	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    s = cache_search(io, GT_Seq, r->rec);
	    if (s->right < s->left)
		continue; // ERROR: no unclipped bases.

	    if ((s->len < 0) ^ r->comp) {
		left = r->start + ABS(s->len) - (s->right-1) - 1;
		right = r->start + ABS(s->len) - (s->left-1) - 1;
	    } else {
		left = r->start + s->left - 1;
		right = r->start + s->right - 1;
	    }

	    while (p1 && p1->pos < left) {
		p1 = p1->next;
		snp_no++;
	    }
	    if (!p1)
		break;

	    if (right < p1->pos)
		continue;

	    snp_no2 = snp_no;
	    for (p2 = p1; p2 && p2->pos <= right; p2 = p2->next) {
		if ((s->len < 0) ^ r->comp) {
		    b = complement_base(s->seq[ABS(s->len)-1 - (p2->pos - r->start)]);
		} else {
		    b = s->seq[p2->pos - r->start];
		}

		hstr[snp_no2++] = b;
	    }

	    //printf("#%"PRIrec": %.*s\n", r->rec, nsnps, hstr);
	    haplotype_str_add(&hs, hstr, nsnps, r->rec);
	}
    }

    haplotype_str_dump(hs);

    {
	haplotype_str *h;
	for (h = hs; hs; hs = hs->next)
	    ArrayPush(rec_list, Array, hs->recs);
    }

    ret = 0;
 err:

    {
	haplotype_pos *curr, *next = NULL;
	for (curr = phead; curr; curr = next) {
	    next = curr->next;
	    free(curr);
	}
    }

    if (cons)
	free(cons);
    if (rng)
	free(rng);
    if (hstr)
	free(hstr);

    return ret;
}

/*
 * Splits readings into haplotypic groups and also returns haplotype consensus?
 * Works via lists? Files?
 *
 * Returns an Array of Arrat of seq record numbers.
 * Returns NULL for failure.a
 */
Array find_haplotypes(GapIO *io, contig_list_t *contigs, int ncontigs) {
    int i;
    Array rec_list = ArrayCreate(sizeof(Array), 0);

    for (i = 0; i < ncontigs; i++) {
	Array recs;
	printf("find_haplotypes =%"PRIrec"\t%d..%d\n",
	       contigs[i].contig, contigs[i].start, contigs[i].end);
	if (-1 == find_haplotypes_single(io, contigs[i].contig,
					 contigs[i].start, contigs[i].end,
					 rec_list)) {
	    // FIXME: free
	    return NULL;
	}
    }

    return rec_list;
}
